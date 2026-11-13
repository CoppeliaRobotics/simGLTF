#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <simPlusPlus-2/Plugin.h>
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

// it seems most GLTF viewers won't work with double-precision floating point buffers...
#ifdef GLTF_DOUBLE
using gltfFloat = double;
static const int TINYGLTF_COMPONENT_TYPE_gltfFloat = TINYGLTF_COMPONENT_TYPE_DOUBLE;
#else
using gltfFloat = float;
static const int TINYGLTF_COMPONENT_TYPE_gltfFloat = TINYGLTF_COMPONENT_TYPE_FLOAT;
#endif

using _handle_t = int;

struct simPose3D
{
    _handle_t handle;
    std::array<double, 3> position;
    std::array<double, 4> orientation;
    bool visible;
};

struct simAnimTrack
{
    _handle_t handle;
    int nodeIndex;
    _handle_t parentHandle;
    std::map<size_t, simPose3D> track;
    std::map<size_t, simPose3D> intrinsicTrack; // only used for joints
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

    bool getGLTFPose(_handle_t handle, _handle_t relTo, tinygltf::Node &node)
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

    int getVisibilityLayers()
    {
        return sim::getIntProperty(sim_handle_scene, "visibilityLayers");
    }

    std::string getObjectName(_handle_t handle)
    {
        return sim::getObjectAlias(handle, 4);
    }

    int getObjectLayers(_handle_t handle)
    {
        try
        {
            return sim::getIntProperty(handle, "layer");
        }
        catch(...)
        {
            return 0;
        }
    }

    bool isCompound(_handle_t handle)
    {
        return sim::getBoolProperty(handle, "compound");
    }

    bool isWireframe(_handle_t handle)
    {
        /*
        try
        {
            return sim::getObjectInt32Param(handle, sim_shapeintparam_wireframe) != 0;
        }
        catch(...)
        {
            return false;
        }
        */
        return false;
    }

    bool isVisible(_handle_t handle)
    {
        _handle_t parentHandle = handle;
        while(parentHandle != -1)
        {
            if(sim::getBoolProperty(parentHandle, "model.notVisible"))
                return false;
            parentHandle = sim::getHandleProperty(parentHandle, "parent");
        }

        int visibleLayers = getVisibilityLayers();
        int layers = getObjectLayers(handle);
        return visibleLayers & layers;
    }

    bool isShape(_handle_t handle)
    {
        std::string objType = sim::getStringProperty(handle, "objectType");
        return objType == "shape";
    }

    bool isJoint(_handle_t handle)
    {
        std::string objType = sim::getStringProperty(handle, "objectType");
        return objType == "joint";
    }

    bool isCamera(_handle_t handle)
    {
        std::string objType = sim::getStringProperty(handle, "objectType");
        return objType == "camera";
    }

    bool isLight(_handle_t handle)
    {
        std::string objType = sim::getStringProperty(handle, "objectType");
        return objType == "light";
    }

    std::vector<_handle_t> ungroupShapeCopy(_handle_t handle)
    {
        auto handles = sim::copyPasteObjects({handle}, 0);
        return sim::ungroupShape(handles[0]);
    }

    void simPose3D_get(simPose3D *p, _handle_t handle, _handle_t relTo)
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
        handleToNode.clear();
        handleToMeshNode.clear();
        handleToChildParent.clear();

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
        for(int i = 0; i < size; ++i)
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

    int addImage(_handle_t id, const void *imgdata, int res[2], const std::string &objname)
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

    int addMesh(_handle_t handle, const std::string &name)
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
        if(info.transparency > 0.001) // is actually opacity when different from 0.0
            alpha = info.transparency;
        model.materials[mat].pbrMetallicRoughness.baseColorFactor = {diffuse[0], diffuse[1], diffuse[2], alpha};
        model.materials[mat].pbrMetallicRoughness.metallicFactor = 0.1;
        model.materials[mat].pbrMetallicRoughness.roughnessFactor = 0.5;
        model.materials[mat].doubleSided = ((info.options & 1) == 0);
        /*
        if(alpha < 1.0)
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
        // parentHandle / parentNodeIndex from the stub are ignored:
        // the real CoppeliaSim hierarchy is used.
        out->nodeIndex = ensureNode(in->shapeHandle);
    }

    void exportObject(exportObject_in *in, exportObject_out *out)
    {
        out->nodeIndex = ensureNode(in->objectHandle);
    }

    void exportObjects(exportObjects_in *in, exportObjects_out *out)
    {
        for(_handle_t obj : in->objectHandles)
            ensureNode(obj);
    }

    void exportAllObjects(exportAllObjects_in *in, exportAllObjects_out *out)
    {
        for(_handle_t obj : getAllObjects())
            ensureNode(obj);
    }

    void exportSelectedObjects(exportSelectedObjects_in *in, exportSelectedObjects_out *out)
    {
        for(_handle_t h : sim::getHandleArrayProperty(sim_handle_scene, "selection"))
            ensureNode(static_cast<int>(h));
    }

    std::vector<_handle_t> getAllObjects()
    {
        std::vector<_handle_t> v;
        for(_handle_t obj : sim::getHandleArrayProperty(sim_handle_scene, "objects"))
            if((isShape(obj) && isVisible(obj) && !isWireframe(obj)) || isCamera(obj))
                 v.push_back(obj);
        return v;
    }

    void exportAnimation(exportAnimation_in *in, exportAnimation_out *out)
    {
        if(times.empty()) return;

        model.animations.resize(1);

        int n = times.size();
        int bt = addBuffer(times.data(), sizeof(gltfFloat) * n, "time");
        int vt = addBufferView(bt, sizeof(gltfFloat) * n, 0, "time");
        std::vector<double> tmin, tmax;
        minMaxVec(times, 1, tmin, tmax);
        int at = addAccessor(vt, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_SCALAR,
                            n, tmin, tmax, "time");

        for(auto &frame : frames)
        {
            simAnimTrack &track = frame.second;
            _handle_t handle = frame.first;                       // NOTE: use the key
            std::string name = model.nodes[track.nodeIndex].name;

            // Scale (visibility) goes on the mesh child if this shape was split,
            // otherwise on the node itself (leaf shape / camera / joint / light).
            auto mit = handleToMeshNode.find(handle);
            int scaleTarget = (mit != handleToMeshNode.end()) ? mit->second : track.nodeIndex;

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

            // --- translation ---
            {
                int c = model.animations[0].channels.size();
                model.animations[0].channels.push_back({});
                model.animations[0].samplers.push_back({});
                int b = addBuffer(p.data(), sizeof(gltfFloat) * n * 3, name + " position");
                int v = addBufferView(b, sizeof(gltfFloat) * n * 3, 0, name + " position");
                std::vector<double> mn, mx; minMaxVec(p, 3, mn, mx);
                int a = addAccessor(v, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3,
                                    n, mn, mx, name + " position");
                model.animations[0].samplers[c].interpolation = "STEP";
                model.animations[0].samplers[c].input  = at;
                model.animations[0].samplers[c].output = a;
                model.animations[0].channels[c].sampler     = c;
                model.animations[0].channels[c].target_node = track.nodeIndex;
                model.animations[0].channels[c].target_path = "translation";
            }

            // --- rotation ---
            {
                int c = model.animations[0].channels.size();
                model.animations[0].channels.push_back({});
                model.animations[0].samplers.push_back({});
                int b = addBuffer(r.data(), sizeof(gltfFloat) * n * 4, name + " rotation");
                int v = addBufferView(b, sizeof(gltfFloat) * n * 4, 0, name + " rotation");
                std::vector<double> mn, mx; minMaxVec(r, 4, mn, mx);
                int a = addAccessor(v, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC4,
                                    n, mn, mx, name + " rotation");
                model.animations[0].samplers[c].interpolation = "STEP";
                model.animations[0].samplers[c].input  = at;
                model.animations[0].samplers[c].output = a;
                model.animations[0].channels[c].sampler     = c;
                model.animations[0].channels[c].target_node = track.nodeIndex;
                model.animations[0].channels[c].target_path = "rotation";
            }

            // --- scale (target = mesh child if split, otherwise the node itself) ---
            {
                int c = model.animations[0].channels.size();
                model.animations[0].channels.push_back({});
                model.animations[0].samplers.push_back({});
                int b = addBuffer(s.data(), sizeof(gltfFloat) * n * 3, name + " scale");
                int v = addBufferView(b, sizeof(gltfFloat) * n * 3, 0, name + " scale");
                std::vector<double> mn, mx; minMaxVec(s, 3, mn, mx);
                int a = addAccessor(v, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3,
                                    n, mn, mx, name + " scale");
                model.animations[0].samplers[c].interpolation = "STEP";
                model.animations[0].samplers[c].input  = at;
                model.animations[0].samplers[c].output = a;
                model.animations[0].channels[c].sampler     = c;
                model.animations[0].channels[c].target_node = scaleTarget;
                model.animations[0].channels[c].target_path = "scale";
            }

            // for joints, also animate the child frame with the intrinsic pose
            if(isJoint(handle))
            {
                auto cpit = handleToChildParent.find(handle);
                if(cpit != handleToChildParent.end())
                {
                    int childFrameNode = cpit->second;

                    std::vector<gltfFloat> cp(n * 3), cr(n * 4);
                    for(int i = 0; i < n; i++)
                    {
                        bool exists = track.intrinsicTrack.find(i) != track.intrinsicTrack.end();
                        for(int j = 0; j < 3; j++)
                            cp[3 * i + j] = exists ? track.intrinsicTrack[i].position[j] : 0.0;
                        for(int j = 0; j < 4; j++)
                            cr[4 * i + j] = exists ? track.intrinsicTrack[i].orientation[j] : 0.0;
                    }

                    // --- intrinsic translation ---
                    {
                        int c = model.animations[0].channels.size();
                        model.animations[0].channels.push_back({});
                        model.animations[0].samplers.push_back({});
                        int b = addBuffer(cp.data(), sizeof(gltfFloat) * n * 3, name + " joint translation");
                        int v = addBufferView(b, sizeof(gltfFloat) * n * 3, 0, name + " joint translation");
                        std::vector<double> mn, mx; minMaxVec(cp, 3, mn, mx);
                        int a = addAccessor(v, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3,
                                            n, mn, mx, name + " joint translation");
                        model.animations[0].samplers[c].interpolation = "STEP";
                        model.animations[0].samplers[c].input  = at;
                        model.animations[0].samplers[c].output = a;
                        model.animations[0].channels[c].sampler     = c;
                        model.animations[0].channels[c].target_node = childFrameNode;
                        model.animations[0].channels[c].target_path = "translation";
                    }

                    // --- intrinsic rotation ---
                    {
                        int c = model.animations[0].channels.size();
                        model.animations[0].channels.push_back({});
                        model.animations[0].samplers.push_back({});
                        int b = addBuffer(cr.data(), sizeof(gltfFloat) * n * 4, name + " joint rotation");
                        int v = addBufferView(b, sizeof(gltfFloat) * n * 4, 0, name + " joint rotation");
                        std::vector<double> mn, mx; minMaxVec(cr, 4, mn, mx);
                        int a = addAccessor(v, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC4,
                                            n, mn, mx, name + " joint rotation");
                        model.animations[0].samplers[c].interpolation = "STEP";
                        model.animations[0].samplers[c].input  = at;
                        model.animations[0].samplers[c].output = a;
                        model.animations[0].channels[c].sampler     = c;
                        model.animations[0].channels[c].target_node = childFrameNode;
                        model.animations[0].channels[c].target_path = "rotation";
                    }
                }
            }
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
        if(!recordAnimationFlag || sim::getIntProperty(sim_handle_scene, "simulationState") != sim_simulation_advancing_running)
            return;

        double time = sim::getFloatProperty(sim_handle_scene, "simulationTime");
        size_t timeIndex = times.size();
        times.push_back(time);

        // Make sure nodes exist for all currently-visible objects
        // (and their ancestors — ensureNode recurses).
        for(_handle_t handle : getAllObjects())
            ensureNode(handle);

        // Record local poses for *every* node we know about, so that
        // dummy nodes (joints, sensors, invisible parents) are animated too.
        // std::map iterators are stable; no insertion happens during the loop.
        for(auto &pair : handleToNode)
        {
            _handle_t handle = pair.first;
            int nodeIndex = pair.second;

            if(frames.find(handle) == frames.end())
            {
                frames[handle].nodeIndex    = nodeIndex;
                frames[handle].parentHandle = getParentHandle(handle);
            }

            try
            {
                simPose3D_get(&frames[handle].track[timeIndex], handle, frames[handle].parentHandle);

                // for joints, also record the intrinsic pose (drives the child frame)
                if(isJoint(handle))
                {
                    auto pose = sim::getPoseProperty(handle, "intrinsicPose");
                    simPose3D &p = frames[handle].intrinsicTrack[timeIndex];
                    p.handle      = handle;
                    p.position    = {pose[0], pose[1], pose[2]};
                    p.orientation = {pose[3], pose[4], pose[5], pose[6]};
                    p.visible     = true;
                }
            }
            catch(...)
            {
                // Handle disappeared (deleted mid-sim). Leave the frame
                // missing so that exportAnimation emits scale=(0,0,0).
            }
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
    std::map<_handle_t, int> handleToNode;
    std::map<_handle_t, int> handleToMeshNode;
    std::map<_handle_t, int> handleToChildParent; // where to attach children (joints)

    _handle_t getParentHandle(_handle_t h)
    {
        return sim::getHandleProperty(h, "parent");
    }

    // If `shapeHandle` is a shape that has just become a parent, split its node
    // into a hierarchy dummy + a child mesh node carrying the visual content.
    // The child node is what gets the scale animation; the parent node stays at
    // scale 1 so that hierarchy children (joints, links, …) remain visible.
    void maybeSplitShapeNode(_handle_t shapeHandle, int shapeNodeIndex)
    {
        if(!isShape(shapeHandle)) return;
        if(handleToMeshNode.find(shapeHandle) != handleToMeshNode.end()) return;

        int meshNode = model.nodes.size();
        model.nodes.push_back({});
        model.nodes[meshNode].name = model.nodes[shapeNodeIndex].name + " [mesh]";

        // 1. Move the node's own mesh, if any (non-compound, visible shape).
        if(model.nodes[shapeNodeIndex].mesh >= 0)
        {
            model.nodes[meshNode].mesh       = model.nodes[shapeNodeIndex].mesh;
            model.nodes[shapeNodeIndex].mesh = -1;
        }

        // 2. Move any mesh-carrying children (a visible compound's sub-meshes)
        //    under the new mesh node, so that scaling the mesh child hides them
        //    while the compound node itself (and its non-mesh children, if any)
        //    are unaffected.
        auto &parentChildren = model.nodes[shapeNodeIndex].children;
        for(auto it = parentChildren.begin(); it != parentChildren.end();)
        {
            int child = *it;
            if(model.nodes[child].mesh >= 0)
            {
                model.nodes[meshNode].children.push_back(child);
                it = parentChildren.erase(it);
            }
            else
            {
                ++it;
            }
        }

        model.nodes[shapeNodeIndex].children.push_back(meshNode);
        handleToMeshNode[shapeHandle] = meshNode;
    }

    // Ensure a glTF node exists for the given handle (and its ancestors),
    // and return its node index.
    int ensureNode(_handle_t handle)
    {
        auto it = handleToNode.find(handle);
        if(it != handleToNode.end())
            return it->second;

        _handle_t parentHandle = getParentHandle(handle);

        int parentNodeIndex = 0;
        if(parentHandle != -1)
        {
            int parentPrimaryNode = ensureNode(parentHandle);
            maybeSplitShapeNode(parentHandle, parentPrimaryNode);
            // Children of a joint attach to the joint's child frame, not the joint frame
            auto cpit = handleToChildParent.find(parentHandle);
            parentNodeIndex = (cpit != handleToChildParent.end())
                            ? cpit->second
                            : parentPrimaryNode;
        }

        int nodeIndex = model.nodes.size();
        model.nodes.push_back({});
        model.nodes[parentNodeIndex].children.push_back(nodeIndex);
        model.nodes[nodeIndex].name = getObjectName(handle);
        getGLTFPose(handle, parentHandle, model.nodes[nodeIndex]);

        handleToNode[handle] = nodeIndex;

        // joints get an extra child-frame node that will carry the intrinsic pose
        if(isJoint(handle))
        {
            int childFrameNode = model.nodes.size();
            model.nodes.push_back({});
            model.nodes[nodeIndex].children.push_back(childFrameNode);
            model.nodes[childFrameNode].name = model.nodes[nodeIndex].name + " [joint child frame]";
            // Initial transform is identity; the animation will drive it.
            handleToChildParent[handle] = childFrameNode;
        }

        if(isShape(handle) && isVisible(handle) && !isWireframe(handle))
            exportShapeContent(handle, nodeIndex);
        else if(isCamera(handle))
            exportCameraContent(handle, nodeIndex);
        else if(isLight(handle))
            exportLightContent(handle, nodeIndex);
        // else: dummy node

        return nodeIndex;
    }

    void exportShapeContent(_handle_t obj, int nodeIndex)
    {
        if(isCompound(obj))
        {
            for(_handle_t subObj : ungroupShapeCopy(obj))
            {
                if(isVisible(subObj) && isShape(subObj) && !isWireframe(subObj))
                {
                    int subNode = model.nodes.size();
                    model.nodes.push_back({});
                    model.nodes[nodeIndex].children.push_back(subNode);
                    model.nodes[subNode].name = getObjectName(subObj);
                    getGLTFPose(subObj, obj, model.nodes[subNode]);
                    model.nodes[subNode].mesh = addMesh(subObj, model.nodes[subNode].name);
                }
                sim::removeObjects({subObj});
            }
            return;
        }
        model.nodes[nodeIndex].mesh = addMesh(obj, model.nodes[nodeIndex].name);
    }

    void exportCameraContent(_handle_t obj, int nodeIndex)
    {
        int cameraIndex = model.cameras.size();
        model.cameras.push_back({});
        model.cameras[cameraIndex].name = getObjectName(obj);
        model.cameras[cameraIndex].type = "perspective";
        model.cameras[cameraIndex].perspective.aspectRatio = 16 / 9.;
        model.cameras[cameraIndex].perspective.yfov =
            sim::getFloatProperty(obj, "viewAngle");
        auto clippingPlanes = sim::getFloatArrayProperty(obj, "clippingPlanes");
        model.cameras[cameraIndex].perspective.znear = clippingPlanes[0];
        model.cameras[cameraIndex].perspective.zfar  = clippingPlanes[1];

        model.nodes[nodeIndex].name = getObjectName(obj) + " camera node";

        // camera-as-such must live on a child node with an orientation fix,
        // so that any children of the camera in CoppeliaSim remain siblings
        // of the fix node (i.e. untouched by the fix).
        int fixNode = model.nodes.size();
        model.nodes.push_back({});
        model.nodes[nodeIndex].children.push_back(fixNode);
        model.nodes[fixNode].name = getObjectName(obj) + " [camera rot fix]";
        model.nodes[fixNode].matrix = { -1, 0, 0, 0,
                                         0, 1, 0, 0,
                                         0, 0,-1, 0,
                                         0, 0, 0, 1 };
        model.nodes[fixNode].camera = cameraIndex;
    }

    void exportLightContent(_handle_t obj, int nodeIndex)
    {
        int lightIndex = model.lights.size();
        model.lights.push_back({});
        model.lights[lightIndex].name = getObjectName(obj);
        auto diffuse = sim::getColorProperty(obj, "color.diffuse");
        model.lights[lightIndex].color     = { diffuse[0], diffuse[1], diffuse[2] };
        model.lights[lightIndex].intensity = 1.0;   // FIXME
        model.lights[lightIndex].type      = "point"; // FIXME
        model.lights[lightIndex].range     = 1.0;   // FIXME
        // glTF 2.0 core has no lights; data is currently orphaned.
        // Node stays a dummy (pose only), so hierarchy is preserved.
        model.nodes[nodeIndex].name = getObjectName(obj) + " [light]";
    }

    tinygltf::TinyGLTF gltf;
    tinygltf::Model model;

    std::map<_handle_t, int> textureMap;

    // for animation data:
    std::map<_handle_t, simAnimTrack> frames;
    std::vector<gltfFloat> times;
    bool recordAnimationFlag = false;

    int bufferPreviewSize = 0;

    TextureFormat exportTextureFormat = simgltf_texture_format_bmp;
};

SIM_PLUGIN(Plugin)
#include "stubsPlusPlus.cpp"
