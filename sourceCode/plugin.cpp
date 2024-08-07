#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <simPlusPlus/Plugin.h>
#include "plugin.h"
#include "stubs.h"
#include "config.h"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "external/tinygltf/tiny_gltf.h"

using simUID = int;

// it seems most GLTF viewers won't work with double-precision floating point buffers...
#ifdef GLTF_DOUBLE
using gltfFloat = double;
static const int TINYGLTF_COMPONENT_TYPE_gltfFloat = TINYGLTF_COMPONENT_TYPE_DOUBLE;
#else
using gltfFloat = float;
static const int TINYGLTF_COMPONENT_TYPE_gltfFloat = TINYGLTF_COMPONENT_TYPE_FLOAT;
#endif

struct simPose3D
{
    int handle;
    std::array<double, 3> position;
    std::array<double, 4> orientation;
    bool visible;
};

struct simAnimTrack
{
    simUID uid;
    int nodeIndex;
    std::map<size_t, simPose3D> track;
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
    if(!v.empty())
    {
        out << '[';
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
        out << "\b\b]";
    }
    return out;
}

class Plugin : public sim::Plugin
{
public:
    void onInit()
    {
        if(!registerScriptStuff())
            throw std::runtime_error("failed to register script stuff");

        char *ps = std::getenv("COPPELIASIM_GLTF_BUFFER_PREVIEW");
        if(ps) bufferPreviewSize = std::atoi(ps);

        setExtVersion("glTF support");
        setBuildDate(BUILD_DATE);

        initAnimationFrames();
    }

    void onSimulationAboutToStart()
    {
        initAnimationFrames();
    }

    void onInstancePass(const sim::InstancePassFlags &flags)
    {
        readAnimationFrame();
    }

    bool getGLTFPose(int handle, int relTo, tinygltf::Node &node)
    {
        auto t = sim::getObjectPosition(handle, relTo);
        auto r = sim::getObjectQuaternion(handle, relTo);
        node.translation = {t[0], t[1], t[2]};
        node.rotation = {r[0], r[1], r[2], r[3]};
        return true;
    }

    template<typename T>
    void minMax(const std::vector<T> &v, int offset, int step, double &min, double &max)
    {
        for(int i = offset; i < v.size(); i += step)
        {
            if(i == offset || v[i] < min) min = v[i];
            if(i == offset || v[i] > max) max = v[i];
        }
    }

    template<typename T>
    void minMaxVec(const std::vector<T> &v, int step, std::vector<double> &min, std::vector<double> &max)
    {
        min.resize(step);
        max.resize(step);
        for(int i = 0; i < step; i++)
            minMax(v, i, step, min[i], max[i]);
    }

    template<typename T>
    void releaseBuffer(const T *b)
    {
        sim::releaseBuffer(b);
    }

    int getVisibleLayers()
    {
        return sim::getInt32Param(sim_intparam_visible_layers);
    }

    std::string getObjectName(int handle)
    {
        return sim::getObjectAlias(handle, 4);
    }

    int getObjectLayers(int handle)
    {
        try
        {
            return sim::getObjectInt32Param(handle, sim_objintparam_visibility_layer);
        }
        catch(...)
        {
            return 0;
        }
    }

    bool is(int handle, int param)
    {
        try
        {
            return sim::getObjectInt32Param(handle, param) != 0;
        }
        catch(...)
        {
            return false;
        }
    }

    bool isCompound(int handle)
    {
        return is(handle, sim_shapeintparam_compound);
    }

    bool isWireframe(int handle)
    {
        return is(handle, sim_shapeintparam_wireframe);
    }

    bool isVisible(int handle)
    {
        int parentHandle = handle;
        while(parentHandle != -1)
        {
            if(sim::getModelProperty(parentHandle) & sim_modelproperty_not_visible)
                return false;
            parentHandle = sim::getObjectParent(parentHandle);
        }

        int visibleLayers = getVisibleLayers();
        int layers = getObjectLayers(handle);
        return visibleLayers & layers;
    }

    bool isShape(int handle)
    {
        int objType = sim::getObjectType(handle);
        return objType == sim_sceneobject_shape;
    }

    bool isCamera(int handle)
    {
        int objType = sim::getObjectType(handle);
        return objType == sim_sceneobject_camera;
    }

    bool isLight(int handle)
    {
        int objType = sim::getObjectType(handle);
        return objType == sim_sceneobject_light;
    }

    std::vector<int> ungroupShapeCopy(int handle)
    {
        auto handles = sim::copyPasteObjects({handle}, 0);
        return sim::ungroupShape(handles[0]);
    }

    void simPose3D_get(simPose3D *p, int handle, int relTo)
    {
        p->handle = handle;
        p->position = sim::getObjectPosition(handle, relTo);
        p->orientation = sim::getObjectQuaternion(handle, relTo);
        p->visible = !isShape(handle) || (isVisible(handle) && !isWireframe(handle));
    }

    std::string buf2str(const void *data, size_t size)
    {
        auto b = reinterpret_cast<const unsigned char*>(data);

        std::string dataPreview;

        if(bufferPreviewSize)
        {
            dataPreview = "{";
            for(int i = 0; i < size; i++)
            {
                dataPreview += i ? " " : "";
                if(i >= bufferPreviewSize)
                {
                    if(i) dataPreview += "...";
                    break;
                }
                dataPreview += boost::lexical_cast<std::string>(int(b[i]));
            }
            dataPreview += "} ";
        }
        dataPreview += (boost::format("(%d bytes)") % size).str();

        return dataPreview;
    }

    void clear(clear_in *in, clear_out *out)
    {
        textureMap.clear();

        model.accessors.clear();
        model.animations.clear();
        model.buffers.clear();
        model.bufferViews.clear();
        model.materials.clear();
        model.meshes.clear();
        model.nodes.clear();
        model.textures.clear();
        model.images.clear();
        model.skins.clear();
        model.samplers.clear();
        model.cameras.clear();
        model.scenes.clear();
        model.lights.clear();
        model.extensionsUsed.clear();
        model.extensionsRequired.clear();

        model.asset.version = "2.0";
        model.asset.generator = "CoppeliaSim glTF plugin";

        model.nodes.push_back({});
        model.nodes[0].name = "Root node";
        model.nodes[0].matrix = {
           -1, 0, 0, 0,
            0, 0, 1, 0,
            0, 1, 0, 0,
            0, 0, 0, 1
        };

        model.scenes.push_back({});
        model.scenes[0].name = "Default Scene";
        model.scenes[0].nodes = {0};

        model.defaultScene = 0;

        frames.clear();
        times.clear();
    }

    void clear()
    {
        clear_in args;
        clear_out ret;
        clear(&args, &ret);
    }

    void loadASCII(loadASCII_in *in, loadASCII_out *out)
    {
        out->result = gltf.LoadASCIIFromFile(&model, &out->errors, &out->warnings, in->filepath);
    }

    void loadBinary(loadBinary_in *in, loadBinary_out *out)
    {
        out->result = gltf.LoadBinaryFromFile(&model, &out->errors, &out->warnings, in->filepath);
    }

    void saveASCII(saveASCII_in *in, saveASCII_out *out)
    {
        out->result = gltf.WriteGltfSceneToFile(&model, in->filepath, true, true, true, false);
    }

    void saveBinary(saveBinary_in *in, saveBinary_out *out)
    {
        out->result = gltf.WriteGltfSceneToFile(&model, in->filepath, true, true, true, true);
    }

    void serialize(serialize_in *in, serialize_out *out)
    {
        std::stringstream ss;
        gltf.WriteGltfSceneToStream(&model, ss, true, false);
        out->json = ss.str();
    }

    int addBuffer(const void *buffer, int size, const std::string &name)
    {
        int i = model.buffers.size();
        model.buffers.push_back({});
        model.buffers[i].data.resize(size);
        model.buffers[i].name = name + " buffer";
        sim::addLog(sim_verbosity_debug, "addBuffer: added buffer '%s' %s", name, buf2str(buffer, size));
        std::memcpy(model.buffers[i].data.data(), buffer, size);
        return i;
    }

    int addBufferView(int buffer, int byteLength, int byteOffset, const std::string &name)
    {
        int i = model.bufferViews.size();
        model.bufferViews.push_back({});
        model.bufferViews[i].buffer = buffer;
        model.bufferViews[i].byteLength = byteLength;
        model.bufferViews[i].byteOffset = byteOffset;
        model.bufferViews[i].name = name + " buffer view";
        return i;
    }

    int addAccessor(int bufferView, int byteOffset, int componentType, int type, int count, const std::vector<double> &minValues, const std::vector<double> &maxValues, const std::string &name)
    {
        int i = model.accessors.size();
        model.accessors.push_back({});
        model.accessors[i].bufferView = bufferView;
        model.accessors[i].byteOffset = byteOffset;
        model.accessors[i].componentType = componentType;
        model.accessors[i].type = type;
        model.accessors[i].count = count;
        model.accessors[i].minValues = minValues;
        model.accessors[i].maxValues = maxValues;
        model.accessors[i].name = name + " accessor";
        return i;
    }

    static void stbi_write_func_vector(void* context, void* data, int size)
    {
        for (int i = 0; i < size; ++i)
        {
            ((std::vector<unsigned char>*)context)->push_back(*(reinterpret_cast<unsigned char *>(data) + i));
        }
    }

    std::vector<unsigned char> raw2bmp(const unsigned char *data, int res[2], int bytesPerPixel)
    {
        int width = res[0];
        int height = res[1];
        int paddingSize = (4 - (width * bytesPerPixel) % 4) % 4;
        int rowSize = bytesPerPixel * width + paddingSize;
        int fileSize = 14 + 40 + rowSize * height;
        std::vector<unsigned char> buf(fileSize, 0);
        buf[ 0] = (unsigned char)('B');
        buf[ 1] = (unsigned char)('M');
        buf[ 2] = (unsigned char)(fileSize);
        buf[ 3] = (unsigned char)(fileSize >> 8);
        buf[ 4] = (unsigned char)(fileSize >> 16);
        buf[ 5] = (unsigned char)(fileSize >> 24);
        buf[10] = (unsigned char)(14 + 40);
        int h = 14;
        buf[h +  0] = (unsigned char)(40);
        buf[h +  4] = (unsigned char)(width);
        buf[h +  5] = (unsigned char)(width >> 8);
        buf[h +  6] = (unsigned char)(width >> 16);
        buf[h +  7] = (unsigned char)(width >> 24);
        buf[h +  8] = (unsigned char)(height);
        buf[h +  9] = (unsigned char)(height >> 8);
        buf[h + 10] = (unsigned char)(height >> 16);
        buf[h + 11] = (unsigned char)(height >> 24);
        buf[h + 12] = (unsigned char)(1);
        buf[h + 14] = (unsigned char)(bytesPerPixel * 8);
        h += 40;
        unsigned char byteOrder[4] = {2, 1, 0, 3};
        for(int i = height - 1; i >= 0; i--)
        {
            for(int j = 0; j < width * bytesPerPixel; j += bytesPerPixel)
                for(int k = 0; k < 4; k++)
                    buf[h++] = data[i * rowSize + j + byteOrder[k]];
            for(int j = 0; j < paddingSize; j++)
                buf[h++] = 0;
        }
        assert(h == fileSize);
        return buf;
    }

    std::vector<unsigned char> raw2jpg(const unsigned char *data, int res[2], int bytesPerPixel, int quality)
    {
        int width = res[0];
        int height = res[1];
        std::vector<unsigned char> buffer;
        auto result = stbi_write_jpg_to_func(stbi_write_func_vector, &buffer, width, height, bytesPerPixel, data, quality);
        return buffer;
    }

    std::vector<unsigned char> raw2png(const unsigned char *data, int res[2], int bytesPerPixel, int stride_bytes)
    {
        int width = res[0];
        int height = res[1];
        std::vector<unsigned char> buffer;
        auto result = stbi_write_png_to_func(stbi_write_func_vector, &buffer, width, height, bytesPerPixel, data, stride_bytes);
        return buffer;
    }

    std::vector<unsigned char> convertRawImage(const unsigned char *imgdata, int res[2])
    {
        if(exportTextureFormat == simgltf_texture_format_bmp)
            return raw2bmp(reinterpret_cast<const unsigned char *>(imgdata), res, 4);
        else if(exportTextureFormat == simgltf_texture_format_png)
            return raw2png(reinterpret_cast<const unsigned char *>(imgdata), res, 4, res[0]*4);
        else if(exportTextureFormat == simgltf_texture_format_jpg)
            return raw2jpg(reinterpret_cast<const unsigned char *>(imgdata), res, 4, 100);
        throw std::runtime_error("unsupported texture format");
    }

    int addImage(int id, const void *imgdata, int res[2], const std::string &objname)
    {
        if(textureMap.find(id) != textureMap.end())
        {
            sim::addLog(sim_verbosity_debug, "addImage: texture with id %d already loaded at index %d", id, textureMap[id]);
            return textureMap[id];
        }
        sim::addLog(sim_verbosity_debug, "addImage: loading texture of object %s with id %d %s", objname, id, buf2str(imgdata, res[0] * res[1] * 4));

        auto buf = convertRawImage(reinterpret_cast<const unsigned char *>(imgdata), res);
        std::string name = (boost::format("texture image %d [%s] (%dx%d, BMP %d bytes)") % id % objname % res[0] % res[1] % buf.size()).str();

        int b = addBuffer(buf.data(), buf.size(), name);

        int v = addBufferView(b, buf.size(), 0, name);

        int i = model.images.size();
        model.images.push_back({});
        model.images[i].bufferView = v;
        model.images[i].mimeType = "image/bmp";
        model.images[i].name = name;
        textureMap[id] = i;
        sim::addLog(sim_verbosity_debug, "addImage: loaded texture of object %s with id %d at index %d %s", objname, id, i, buf2str(imgdata, res[0] * res[1] * 4));
        sim::addLog(sim_verbosity_debug, "addImage: model now has %d images", model.images.size());
        return i;
    }

    void expandVertices(double *vertices, int verticesSize, int *indices, int indicesSize, double *normals, float *texCoords, std::vector<gltfFloat> &vertices2, std::vector<int> &indices2, std::vector<gltfFloat> &normals2, std::vector<gltfFloat> &texCoords2)
    {
        vertices2.resize(3 * indicesSize);
        indices2.resize(indicesSize);
        normals2.resize(3 * indicesSize);
        if(texCoords)
            texCoords2.resize(2 * indicesSize);
        else
            texCoords2.clear();
        int iv = 0, ii = 0, in = 0, it = 0;
        for(int i = 0; i < indicesSize; i += 3)
        {
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < 3; k++)
                    vertices2[iv++] = vertices[3*indices[i+j]+k];
            for(int j = 0; j < 3; j++)
                indices2[ii++] = i+j;
            for(int j = 0; j < 9; j++)
                normals2[in++] = normals[3*i+j];
            if(texCoords)
                for(int j = 0; j < 6; j++)
                    texCoords2[it++] = texCoords[2*i+j];
        }
    }

    int addMesh(int handle, const std::string &name)
    {
        sim::addLog(sim_verbosity_debug, "addMesh: %s: adding mesh for shape handle %d", name, handle);

        struct SShapeVizInfo info;
        int result = sim::getShapeViz(handle | sim_handleflag_extended, 0, &info);
        bool hasTexture = result == 2;
        sim::addLog(sim_verbosity_debug, "addMesh: %s: has texture: %d (result %d)", name, hasTexture, result);

        std::vector<gltfFloat> vertices2;
        std::vector<int> indices2;
        std::vector<gltfFloat> normals2;
        std::vector<gltfFloat> texCoords2;
        expandVertices(info.vertices, info.verticesSize, info.indices, info.indicesSize, info.normals, hasTexture ? info.textureCoords : 0L, vertices2, indices2, normals2, texCoords2);
        releaseBuffer(info.vertices);
        releaseBuffer(info.indices);
        releaseBuffer(info.normals);

        std::vector<float> diffuse(&info.colors[0], &info.colors[0] + 3);
        std::vector<float> specular(&info.colors[3], &info.colors[3] + 3);
        std::vector<float> emission(&info.colors[6], &info.colors[6] + 3);

        int bv = addBuffer(vertices2.data(), sizeof(gltfFloat) * vertices2.size(), name + " vertex");
        int bi = addBuffer(indices2.data(), sizeof(int) * indices2.size(), name + " index");
        int bn = addBuffer(normals2.data(), sizeof(gltfFloat) * normals2.size(), name + " normal");

        int vv = addBufferView(bv, sizeof(gltfFloat) * vertices2.size(), 0, name + " vertex");
        int vi = addBufferView(bi, sizeof(int) * indices2.size(), 0, name + " index");
        int vn = addBufferView(bn, sizeof(gltfFloat) * normals2.size(), 0, name + " normal");

        std::vector<double> vmin, vmax, imin, imax, nmin, nmax, tmin, tmax;
        minMaxVec(vertices2, 3, vmin, vmax);
        int av = addAccessor(vv, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3, vertices2.size() / 3, vmin, vmax, name + " vertex");
        minMaxVec(indices2, 1, imin, imax);
        int ai = addAccessor(vi, 0, TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT, TINYGLTF_TYPE_SCALAR, indices2.size(), imin, imax, name + " index");
        minMaxVec(normals2, 3, nmin, nmax);
        int an = addAccessor(vn, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3, normals2.size() / 3, nmin, nmax, name + " normal");

        int mesh = model.meshes.size();
        model.meshes.push_back({});
        model.meshes[mesh].name = name + " mesh";
        model.meshes[mesh].primitives.push_back({});
        model.meshes[mesh].primitives[0].attributes["POSITION"] = av;
        model.meshes[mesh].primitives[0].attributes["NORMAL"] = an;
        model.meshes[mesh].primitives[0].indices = ai;
        model.meshes[mesh].primitives[0].mode = TINYGLTF_MODE_TRIANGLES;

        int mat = model.materials.size();
        model.meshes[mesh].primitives[0].material = mat;
        model.materials.push_back({});
        model.materials[mat].name = name + " material";
        model.materials[mat].emissiveFactor = {emission[0], emission[1], emission[2]};
        double alpha = 1.0;
        if (info.transparency > 0.001) // is actually opacity when different from 0.0
            alpha = info.transparency;
        model.materials[mat].pbrMetallicRoughness.baseColorFactor = {diffuse[0], diffuse[1], diffuse[2], alpha};
        model.materials[mat].pbrMetallicRoughness.metallicFactor = 0.1;
        model.materials[mat].pbrMetallicRoughness.roughnessFactor = 0.5;
        model.materials[mat].doubleSided = ((info.options & 1) == 0);
        /*
        if (alpha < 1.0)
        {
             model.materials[mat].alphaMode = "BLEND";
             model.materials[mat].alphaCutoff = alpha;
        }
        */

        if(hasTexture)
        {
            int bt = addBuffer(texCoords2.data(), sizeof(gltfFloat) * texCoords2.size(), name + " texture coord");
            int vt = addBufferView(bt, sizeof(gltfFloat) * texCoords2.size(), 0, name + " texture coord");
            minMaxVec(texCoords2, 2, tmin, tmax);
            int at = addAccessor(vt, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC2, texCoords2.size() / 2, tmin, tmax, name + " texture coord");
            model.meshes[mesh].primitives[0].attributes["TEXCOORD_0"] = at;
            bool repU = info.textureOptions & 1;
            bool repV = info.textureOptions & 2;
            bool interp = info.textureOptions & 4;
            int sampler = model.samplers.size();
            model.samplers.push_back({});
            model.samplers[sampler].magFilter = interp ? TINYGLTF_TEXTURE_FILTER_LINEAR : TINYGLTF_TEXTURE_FILTER_NEAREST;
            model.samplers[sampler].minFilter = TINYGLTF_TEXTURE_FILTER_LINEAR_MIPMAP_LINEAR;
            model.samplers[sampler].wrapS = repU ? TINYGLTF_TEXTURE_WRAP_REPEAT : TINYGLTF_TEXTURE_WRAP_CLAMP_TO_EDGE;
            model.samplers[sampler].wrapT = repV ? TINYGLTF_TEXTURE_WRAP_REPEAT : TINYGLTF_TEXTURE_WRAP_CLAMP_TO_EDGE;
            int tex = model.textures.size();
            model.textures.push_back({});
            model.textures[tex].name = name + " texture";
            model.textures[tex].sampler = sampler;
            model.textures[tex].source = addImage(info.textureId, info.texture, info.textureRes, name);
            model.materials[mat].pbrMetallicRoughness.baseColorTexture.texCoord = 0; // will use TEXCOORD_0
            model.materials[mat].pbrMetallicRoughness.baseColorTexture.index = tex;
            releaseBuffer(info.textureCoords);
            releaseBuffer(info.texture);
        }

        return mesh;
    }

    void exportShape(exportShape_in *in, exportShape_out *out)
    {
        int obj = in->shapeHandle;
        out->nodeIndex = model.nodes.size();
        model.nodes.push_back({});
        model.nodes[in->parentNodeIndex].children.push_back(out->nodeIndex);
        model.nodes[out->nodeIndex].name = getObjectName(obj);
        getGLTFPose(obj, in->parentHandle, model.nodes[out->nodeIndex]);

        if(isCompound(obj))
        {
            for(int subObj : ungroupShapeCopy(obj))
            {
                if(isVisible(subObj) && isShape(subObj) && !isWireframe(subObj))
                {
                    exportShape_in args;
                    args._ = in->_;
                    args.shapeHandle = subObj;
                    args.parentHandle = obj;
                    args.parentNodeIndex = out->nodeIndex;
                    exportShape_out ret;
                    exportShape(&args, &ret);
                }
                sim::removeObjects({subObj});
            }
            return;
        }

        model.nodes[out->nodeIndex].mesh = addMesh(obj, model.nodes[out->nodeIndex].name);
    }

    void exportObject(exportObject_in *in, exportObject_out *out)
    {
        int visibleLayers = getVisibleLayers();
        int obj = in->objectHandle;

        if(isShape(obj) && isVisible(obj) && !isWireframe(obj))
        {
            exportShape_in args;
            args._ = in->_;
            args.shapeHandle = obj;
            exportShape_out ret;
            exportShape(&args, &ret);
            out->nodeIndex = ret.nodeIndex;
        }
        else if(isCamera(obj))
        {
            int cameraIndex = model.cameras.size();
            model.cameras.push_back({});
            model.cameras[cameraIndex].name = getObjectName(obj);
#if 0
            if(sim::getObjectInt32Param(obj, sim_cameraintparam_perspective_operation))
            {
                model.cameras[cameraIndex].type = "perspective";
                model.cameras[cameraIndex].perspective.aspectRatio = 16/9.;
                model.cameras[cameraIndex].perspective.yfov = sim::getObjectFloatParam(obj, sim_camerafloatparam_perspective_angle);
                model.cameras[cameraIndex].perspective.znear = sim::getObjectFloatParam(obj, sim_camerafloatparam_near_clipping);
                model.cameras[cameraIndex].perspective.zfar = sim::getObjectFloatParam(obj, sim_camerafloatparam_far_clipping);
            }
            else
            {
                model.cameras[cameraIndex].type = "orthographic";
                model.cameras[cameraIndex].orthographic.xmag = sim::getObjectInt32Param(obj, sim_cameraintparam_resolution_x);
                model.cameras[cameraIndex].orthographic.ymag = sim::getObjectInt32Param(obj, sim_cameraintparam_resolution_y);
                model.cameras[cameraIndex].orthographic.znear = sim::getObjectFloatParam(obj, sim_camerafloatparam_near_clipping);
                model.cameras[cameraIndex].orthographic.zfar = sim::getObjectFloatParam(obj, sim_camerafloatparam_far_clipping);
            }
#else
            model.cameras[cameraIndex].type = "perspective";
            model.cameras[cameraIndex].perspective.aspectRatio = 16/9.;
            model.cameras[cameraIndex].perspective.yfov = sim::getObjectFloatParam(obj, sim_camerafloatparam_perspective_angle);
            model.cameras[cameraIndex].perspective.znear = 0.001;
            model.cameras[cameraIndex].perspective.zfar = 1000.0;
#endif
            out->nodeIndex = model.nodes.size();
            model.nodes.push_back({});
            model.nodes[0].children.push_back(out->nodeIndex);
            int cameraOrientationFixNode = model.nodes.size();
            model.nodes.push_back({});
            model.nodes[out->nodeIndex].name = getObjectName(obj) + " camera node";
            model.nodes[out->nodeIndex].children.push_back(cameraOrientationFixNode);
            model.nodes[cameraOrientationFixNode].matrix = {-1,  0,  0,  0,
                                                             0,  1,  0,  0,
                                                             0,  0, -1,  0,
                                                             0,  0,  0,  1};
            model.nodes[cameraOrientationFixNode].camera = cameraIndex;
            model.nodes[cameraOrientationFixNode].name = getObjectName(obj) + " camera node [rot fix]";
            getGLTFPose(obj, -1, model.nodes[out->nodeIndex]);
        }
        else if(isLight(obj))
        {
            int lightIndex = model.lights.size();
            model.lights.push_back({});
            model.lights[lightIndex].name = getObjectName(obj);
            std::array<double, 3> diffuse;
            bool lightOn = sim::getLightParameters(obj, diffuse) & 1;
            model.lights[lightIndex].color = {diffuse[0], diffuse[1], diffuse[2]};
            model.lights[lightIndex].intensity = 1.0; // FIXME: where to get this value from sim?
            model.lights[lightIndex].type = "point"; // FIXME: where to get this value from sim?
            model.lights[lightIndex].range = 1.0; // FIXME: where to get this value from sim?
            out->nodeIndex = lightIndex;
        }
    }

    std::vector<int> getAllObjects()
    {
        std::vector<int> v;
        for(int obj : sim::getObjectsInTree(sim_handle_scene, sim_handle_all))
            if((isShape(obj) && isVisible(obj) && !isWireframe(obj)) || isCamera(obj))
                 v.push_back(obj);
        return v;
    }

    void exportAllObjects(exportAllObjects_in *in, exportAllObjects_out *out)
    {
        exportObjects_in args;
        args._ = in->_;
        args.objectHandles = getAllObjects();
        if(args.objectHandles.empty()) return;
        exportObjects_out ret;
        exportObjects(&args, &ret);
    }

    void exportSelectedObjects(exportSelectedObjects_in *in, exportSelectedObjects_out *out)
    {
        exportObjects_in args;
        args._ = in->_;
        args.objectHandles = sim::getObjectSel();
        if(args.objectHandles.empty()) return;
        exportObjects_out ret;
        exportObjects(&args, &ret);
    }

    void exportObjects(exportObjects_in *in, exportObjects_out *out)
    {
        exportObject_in args;
        args._ = in->_;
        exportObject_out ret;
        for(int obj : in->objectHandles)
        {
            args.objectHandle = obj;
            exportObject(&args, &ret);
        }
    }

    void exportAnimation(exportAnimation_in *in, exportAnimation_out *out)
    {
        model.animations.resize(1);

        // create time buffer:
        int n = times.size();
        int bt = addBuffer(times.data(), sizeof(gltfFloat) * n, "time");
        int vt = addBufferView(bt, sizeof(gltfFloat) * n, 0, "time");
        std::vector<double> tmin, tmax;
        minMaxVec(times, 1, tmin, tmax);
        int at = addAccessor(vt, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_SCALAR, n, tmin, tmax, "time");

        for(auto &frame : frames)
        {
            simAnimTrack &track = frame.second;

            // for each object we need two channels: translation and rotation
            int ip = model.animations[0].channels.size();
            model.animations[0].channels.push_back({});
            model.animations[0].samplers.push_back({});
            int ir = model.animations[0].channels.size();
            model.animations[0].channels.push_back({});
            model.animations[0].samplers.push_back({});
            // XXX: animate visibility with scale channel
            int is = model.animations[0].channels.size();
            model.animations[0].channels.push_back({});
            model.animations[0].samplers.push_back({});

            // create translation and rotation buffers:
            std::string name = model.nodes[track.nodeIndex].name;
            std::vector<gltfFloat> p(n * 3), r(n * 4), s(n * 3);
            for(int i = 0; i < n; i++)
            {
                bool exists = track.track.find(i) != track.track.end();
                for(int j = 0; j < 3; j++)
                    p[3 * i + j] = exists ? track.track[i].position[j] : 0.0;
                for(int j = 0; j < 4; j++)
                    r[4 * i + j] = exists ? track.track[i].orientation[j] : 0.0;
                for(int j = 0; j < 3; j++)
                    s[3 * i + j] = exists && track.track[i].visible ? 1.0 : 0.0;
            }

            int bp = addBuffer(p.data(), sizeof(gltfFloat) * n * 3, name + " position");
            int vp = addBufferView(bp, sizeof(gltfFloat) * n * 3, 0, name + " position");
            std::vector<double> pmin, pmax;
            minMaxVec(p, 3, pmin, pmax);
            int ap = addAccessor(vp, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3, n, pmin, pmax, name + " position");

            int br = addBuffer(r.data(), sizeof(gltfFloat) * n * 4, name + " rotation");
            int vr = addBufferView(br, sizeof(gltfFloat) * n * 4, 0, name + " rotation");
            std::vector<double> rmin, rmax;
            minMaxVec(r, 4, rmin, rmax);
            int ar = addAccessor(vr, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC4, n, rmin, rmax, name + " rotation");

            int bs = addBuffer(s.data(), sizeof(gltfFloat) * n * 3, name + " scale");
            int vs = addBufferView(bs, sizeof(gltfFloat) * n * 3, 0, name + " scale");
            std::vector<double> smin, smax;
            minMaxVec(s, 3, smin, smax);
            int as = addAccessor(vs, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3, n, smin, smax, name + " scale");

            // create samplers & channels:
            model.animations[0].samplers[ip].interpolation = "STEP";
            model.animations[0].samplers[ip].input = at;
            model.animations[0].samplers[ip].output = ap;
            model.animations[0].channels[ip].sampler = ip;
            model.animations[0].channels[ip].target_node = track.nodeIndex;
            model.animations[0].channels[ip].target_path = "translation";
            model.animations[0].samplers[ir].interpolation = "STEP";
            model.animations[0].samplers[ir].input = at;
            model.animations[0].samplers[ir].output = ar;
            model.animations[0].channels[ir].sampler = ir;
            model.animations[0].channels[ir].target_node = track.nodeIndex;
            model.animations[0].channels[ir].target_path = "rotation";
            model.animations[0].samplers[is].interpolation = "STEP";
            model.animations[0].samplers[is].input = at;
            model.animations[0].samplers[is].output = as;
            model.animations[0].channels[is].sampler = is;
            model.animations[0].channels[is].target_node = track.nodeIndex;
            model.animations[0].channels[is].target_path = "scale";
        }
    }

    void animationFrameCount(animationFrameCount_in *in, animationFrameCount_out *out)
    {
        out->count = frames.size();
    }

    void recordAnimation(recordAnimation_in *in, recordAnimation_out *out)
    {
        recordAnimationFlag = in->enable;
    }

    void initAnimationFrames()
    {
        clear();
    }

    void readAnimationFrame()
    {
        if(!recordAnimationFlag || sim::getSimulationState() != sim_simulation_advancing_running)
            return;

        double time = sim::getSimulationTime();
        size_t timeIndex = times.size();
        times.push_back(time);

        std::vector<int> allObjects = getAllObjects();
        for(int handle : allObjects)
        {
            long long int uid = sim::getObjectUid(handle);
            auto it = frames.find(uid);
            if(it == frames.end())
            {
                exportObject_in args;
                exportObject_out ret;
                args.objectHandle = handle;
                exportObject(&args, &ret);
                frames[uid].nodeIndex = ret.nodeIndex;
            }

            simPose3D_get(&frames[uid].track[timeIndex], handle, -1);
        }
    }

    void setExportTextureFormat(setExportTextureFormat_in *in, setExportTextureFormat_out *out)
    {
        exportTextureFormat = static_cast<TextureFormat>(in->textureFormat);
    }

    void getExportTextureFormat(getExportTextureFormat_in *in, getExportTextureFormat_out *out)
    {
        out->textureFormat = exportTextureFormat;
        out->formatName = textureformat_string(exportTextureFormat);
    }

private:
    tinygltf::TinyGLTF gltf;
    tinygltf::Model model;

    std::map<int, int> textureMap;

    // for animation data:
    std::map<simUID, simAnimTrack> frames;
    std::vector<gltfFloat> times;
    bool recordAnimationFlag = false;

    int bufferPreviewSize = 0;

    TextureFormat exportTextureFormat = simgltf_texture_format_bmp;
};

SIM_PLUGIN(Plugin)
#include "stubsPlusPlus.cpp"
