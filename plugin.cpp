#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "simPlusPlus/Plugin.h"
#include "simPlusPlus/Handle.h"
#include "plugin.h"
#include "stubs.h"
#include "config.h"

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "external/tinygltf/tiny_gltf.h"

tinygltf::TinyGLTF gltf;

using sim::Handle;

struct simPose3D
{
    simInt handle;
    simFloat position[3];
    simFloat orientation[4];

    void get(int handle, int relTo)
    {
        this->handle = handle;
        simGetObjectPosition(handle, relTo, &position[0]);
        simGetObjectQuaternion(handle, relTo, &orientation[0]);
    }
};

struct simAnimFrame
{
    simFloat time;
    std::vector<simPose3D> poses;

    void read(simFloat time, const std::vector<int> &handles)
    {
        this->time = time;
        poses.resize(handles.size());
        size_t i = 0;
        for(const auto &handle : handles)
            poses[i++].get(handle, -1);
    }
};

// for animation data:
std::vector<int> handles;
std::map<int, size_t> handleIndex;
std::vector<simAnimFrame> frames;
std::map<int, size_t> nodeIndex;

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

template<> std::string Handle<tinygltf::Model>::tag()
{
    return PLUGIN_NAME ".model";
}

tinygltf::Model * getModel(const std::string &handle)
{
    tinygltf::Model *model = Handle<tinygltf::Model>::obj(handle);
    if(!model) throw std::runtime_error("invalid glTF model handle");
    return model;
}

void reportError(const std::string &message)
{
    std::string m = "glTF: " + message;
    simAddStatusbarMessage(m.c_str());
    std::cerr << m << std::endl;
}

void create(SScriptCallBack *p, const char *cmd, create_in *in, create_out *out)
{
    tinygltf::Model *model = new tinygltf::Model;

    model->asset.version = "2.0";
    model->asset.generator = "CoppeliaSim glTF plugin";

    model->samplers.push_back({});
    model->samplers[0].magFilter = TINYGLTF_TEXTURE_FILTER_LINEAR;
    model->samplers[0].minFilter = TINYGLTF_TEXTURE_FILTER_LINEAR_MIPMAP_LINEAR;
    model->samplers[0].wrapS = TINYGLTF_TEXTURE_WRAP_MIRRORED_REPEAT;
    model->samplers[0].wrapT = TINYGLTF_TEXTURE_WRAP_MIRRORED_REPEAT;

    model->nodes.push_back({});
    model->nodes[0].name = "Root node";
    model->nodes[0].matrix = {
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, 1, 0, 0,
        0, 0, 0, 1
    };

    model->scenes.push_back({});
    model->scenes[0].name = "Default Scene";
    model->scenes[0].nodes = {0};

    model->defaultScene = 0;

    out->handle = Handle<tinygltf::Model>::str(model);
}

void destroy(SScriptCallBack *p, const char *cmd, destroy_in *in, destroy_out *out)
{
    auto model = getModel(in->handle);
    delete model;
}

void loadASCII(SScriptCallBack *p, const char *cmd, loadASCII_in *in, loadASCII_out *out)
{
    tinygltf::Model *model = new tinygltf::Model;
    if(gltf.LoadASCIIFromFile(model, &out->errors, &out->warnings, in->filepath))
        out->handle = Handle<tinygltf::Model>::str(model);
    else
        delete model;
}

void loadBinary(SScriptCallBack *p, const char *cmd, loadBinary_in *in, loadBinary_out *out)
{
    tinygltf::Model *model = new tinygltf::Model;
    if(gltf.LoadBinaryFromFile(model, &out->errors, &out->warnings, in->filepath))
        out->handle = Handle<tinygltf::Model>::str(model);
    else
        delete model;
}

void saveASCII(SScriptCallBack *p, const char *cmd, saveASCII_in *in, saveASCII_out *out)
{
    auto model = getModel(in->handle);
    out->success = gltf.WriteGltfSceneToFile(model, in->filepath, true, true, true, false);
}

void saveBinary(SScriptCallBack *p, const char *cmd, saveBinary_in *in, saveBinary_out *out)
{
    auto model = getModel(in->handle);
    out->success = gltf.WriteGltfSceneToFile(model, in->filepath, true, true, true, true);
}

void serialize(SScriptCallBack *p, const char *cmd, serialize_in *in, serialize_out *out)
{
    auto model = getModel(in->handle);
    std::stringstream ss;
    gltf.WriteGltfSceneToStream(model, ss, true, false);
    out->json = ss.str();
}

bool getGLTFMatrix(int handle, int relTo, std::vector<double> &m)
{
    simFloat x[12];
    if(simGetObjectMatrix(handle, relTo, &x[0]) == -1)
        return false;
    m = {
        x[0], x[4], x[8], 0,
        x[1], x[5], x[9], 0,
        x[2], x[6], x[10], 0,
        x[3], x[7], x[11], 1
    };
    return true;
}

template<typename T>
void minMax(const T *v, simInt size, simInt offset, simInt step, double &min, double &max)
{
    for(int i = offset; i < size; i += step)
    {
        if(i == offset || v[i] < min) min = v[i];
        if(i == offset || v[i] > max) max = v[i];
    }
}

template<typename T>
void minMaxVec(const T *v, simInt size, simInt step, std::vector<double> &min, std::vector<double> &max)
{
    min.resize(step);
    max.resize(step);
    for(int i = 0; i < step; i++)
        minMax(v, size, i, step, min[i], max[i]);
}

template<typename T>
void releaseBuffer(const T *b)
{
    simReleaseBuffer(reinterpret_cast<const char *>(b));
}

int addBuffer(tinygltf::Model *model, void *buffer, int size, const std::string &name)
{
    int i = model->buffers.size();
    model->buffers.push_back({});
    tinygltf::Buffer &b = model->buffers[i];
    b.data.resize(size);
    b.name = name + " buffer";
    std::memcpy(b.data.data(), buffer, size);
    return i;
}

int addBufferView(tinygltf::Model *model, int buffer, int byteLength, int byteOffset, const std::string &name)
{
    int i = model->bufferViews.size();
    model->bufferViews.push_back({});
    tinygltf::BufferView &v = model->bufferViews[i];
    v.buffer = buffer;
    v.byteLength = byteLength;
    v.byteOffset = byteOffset;
    v.name = name + " buffer view";
    return i;
}

int addAccessor(tinygltf::Model *model, int bufferView, int byteOffset, int componentType, int type, int count, const std::vector<double> &minValues, const std::vector<double> &maxValues, const std::string &name)
{
    int i = model->accessors.size();
    model->accessors.push_back({});
    tinygltf::Accessor &a = model->accessors[i];
    a.bufferView = bufferView;
    a.byteOffset = byteOffset;
    a.componentType = componentType;
    a.type = type;
    a.count = count;
    a.minValues = minValues;
    a.maxValues = maxValues;
    a.name = name + " accessor";
    return i;
}

int addMesh(tinygltf::Model *model, int handle, const std::string &name)
{
    simFloat *vertices;
    simInt verticesSize;
    simInt *indices;
    simInt indicesSize;
    simFloat *normals;
    if(simGetShapeMesh(handle, &vertices, &verticesSize, &indices, &indicesSize, &normals) == -1)
        return -1;

    int bv = addBuffer(model, vertices, sizeof(simFloat) * verticesSize, name + " vertex");
    int bi = addBuffer(model, indices, sizeof(simInt) * indicesSize, name + " index");
    int bn = addBuffer(model, normals, sizeof(simFloat) * indicesSize * 3, name + " normal");

    int vv = addBufferView(model, bv, sizeof(simFloat) * verticesSize, 0, name + " vertex");
    int vi = addBufferView(model, bi, sizeof(simInt) * indicesSize, 0, name + " index");
    int vn = addBufferView(model, bn, sizeof(simFloat) * indicesSize * 3, 0, name + " normal");

    std::vector<double> vmin, vmax, imin, imax, nmin, nmax;
    minMaxVec(vertices, verticesSize, 3, vmin, vmax);
    int av = addAccessor(model, vv, 0, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3, verticesSize / 3, vmin, vmax, name + " vertex");
    minMaxVec(indices, indicesSize, 1, imin, imax);
    int ai = addAccessor(model, vi, 0, TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT, TINYGLTF_TYPE_SCALAR, indicesSize, imin, imax, name + " index");
    minMaxVec(normals, indicesSize * 3, 3, nmin, nmax);
    int an = addAccessor(model, vn, 0, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3, indicesSize, nmin, nmax, name + " normal");

    int i = model->meshes.size();
    model->meshes.push_back({});
    tinygltf::Mesh &mesh = model->meshes[i];
    mesh.name = name + " mesh";
    mesh.primitives.push_back({});
    mesh.primitives[0].attributes["POSITION"] = av;
    mesh.primitives[0].attributes["NORMAL"] = an;
    mesh.primitives[0].indices = ai;
    mesh.primitives[0].mode = TINYGLTF_MODE_TRIANGLES;
    return i;
}

int addBasicMaterial(tinygltf::Model *model, const std::vector<simFloat> &col1, const std::vector<simFloat> &col2, const std::string &name)
{
    model->materials.resize(model->materials.size() + 1);
    tinygltf::Material &material = model->materials[model->materials.size() - 1];
    material.name = name + " material";
    material.pbrMetallicRoughness.baseColorFactor = {col1[0], col1[1], col1[2], 1.0};
    material.pbrMetallicRoughness.metallicFactor = 0.1;
    material.pbrMetallicRoughness.roughnessFactor = 0.5;
    return model->materials.size() - 1;
}

int addTexture(tinygltf::Model *model, void *imgdata)
{
    //local i=#model.images
    //table.insert(model.images,{uri='data:application/octet-stream;base64,'..b64enc(imgdata)})
    //return i
    return 0;
}

void getObjectSelection(std::vector<simInt> &v)
{
    int selectionSize = simGetObjectSelectionSize();
    v.resize(selectionSize);
    simGetObjectSelection(v.data());
}

simInt getVisibleLayers()
{
    simInt v = 0;
    if(simGetInt32Parameter(sim_intparam_visible_layers, &v) == -1)
        return 0;
    else
        return v;
}

std::string getObjectName(simInt handle)
{
    simChar *name = simGetObjectName(handle);
    std::string ret;
    if(name)
    {
        ret = name;
        releaseBuffer(name);
    }
    return ret;
}

simInt getObjectLayers(simInt handle)
{
    simInt v = 0;
    if(simGetObjectInt32Parameter(handle, sim_objintparam_visibility_layer, &v) == 1)
        return v;
    else
        return 0;
}

bool is(simInt handle, simInt param)
{
    simInt v = 0;
    if(simGetObjectInt32Parameter(handle, param, &v) == 1)
        return v != 0;
    else
        return false;
}

bool isCompound(simInt handle)
{
    return is(handle, sim_shapeintparam_compound);
}

bool isWireframe(simInt handle)
{
    return is(handle, sim_shapeintparam_wireframe);
}

std::vector<simInt> ungroupShape(simInt handle)
{
    std::vector<simInt> ret;
    simInt count;
    simInt *shapes = simUngroupShape(handle, &count);
    if(shapes)
    {
        ret.resize(count);
        for(int i = 0; i < count; i++)
            ret[i] = shapes[i];
        releaseBuffer(shapes);
    }
    return ret;
}

std::vector<simInt> ungroupShapeCopy(simInt handle)
{
    simInt handles[1] = {handle};
    simCopyPasteObjects(handles, 1, 0);
    return ungroupShape(handles[0]);
}

std::vector<simFloat> getShapeColor(simInt handle, simInt colorComponent)
{
    std::vector<simFloat> ret;
    ret.resize(3);
    simGetShapeColor(handle, 0, colorComponent, ret.data());
    return ret;
}

void exportShape(SScriptCallBack *p, const char *cmd, exportShape_in *in, exportShape_out *out)
{
    auto model = getModel(in->handle);
    simInt obj = in->shapeHandle;
    out->nodeIndex = model->nodes.size();
    model->nodes.push_back({});
    model->nodes[in->parentNodeIndex].children.push_back(out->nodeIndex);
    model->nodes[out->nodeIndex].name = getObjectName(obj);
    getGLTFMatrix(obj, in->parentHandle, model->nodes[out->nodeIndex].matrix);

    if(isCompound(obj))
    {
        for(simInt subObj : ungroupShapeCopy(obj))
        {
            exportShape_in args;
            args.handle = in->handle;
            args.shapeHandle = subObj;
            args.parentHandle = obj;
            args.parentNodeIndex = out->nodeIndex;
            exportShape_out ret;
            exportShape(p, &args, &ret);
            simRemoveObject(subObj);
        }
        return;
    }

    std::vector<simFloat> col1 = getShapeColor(obj, sim_colorcomponent_ambient_diffuse);
    std::vector<simFloat> col2 = getShapeColor(obj, sim_colorcomponent_specular);
    model->nodes[out->nodeIndex].mesh = addMesh(model, obj, model->nodes[out->nodeIndex].name);
    model->meshes[model->nodes[out->nodeIndex].mesh].primitives[0].material = addBasicMaterial(model, col1, col2, model->nodes[out->nodeIndex].name);
}

void exportObject(SScriptCallBack *p, const char *cmd, exportObject_in *in, exportObject_out *out)
{
    auto model = getModel(in->handle);
    simInt visibleLayers = getVisibleLayers();
    simInt obj = in->objectHandle;
    simInt layers = getObjectLayers(obj);
    bool visible = visibleLayers & layers;

    simInt objType = simGetObjectType(obj);
    if(objType == sim_object_shape_type && visible)
    {
        exportShape_in args;
        args.handle = in->handle;
        args.shapeHandle = obj;
        exportShape_out ret;
        exportShape(p, &args, &ret);
        out->nodeIndex = ret.nodeIndex;
    }
    if(objType == sim_object_camera_type)
    {
        int cameraIndex = model->cameras.size();
        model->cameras.push_back({});
        model->cameras.back().type = "perspective";
        model->cameras.back().perspective.aspectRatio = 16/9.;
        simFloat a;
        if(simGetObjectFloatParameter(obj, sim_camerafloatparam_perspective_angle, &a) == 1)
            model->cameras.back().perspective.yfov = a;
        model->cameras.back().perspective.znear = 0.001;
        model->cameras.back().perspective.zfar = 1000;
        out->nodeIndex = model->nodes.size();
        model->nodes.push_back({});
        model->nodes[out->nodeIndex].camera = cameraIndex;
        model->nodes[out->nodeIndex].name = getObjectName(obj);
        getGLTFMatrix(obj, -1, model->nodes[out->nodeIndex].matrix);
    }

    nodeIndex[obj] = out->nodeIndex;
}

void getAllObjects(std::vector<simInt> &v)
{
    simInt allObjectsCount;
    simInt *allObjectsBuf = simGetObjectsInTree(sim_handle_scene, sim_handle_all, 0, &allObjectsCount);
    if(allObjectsBuf)
    {
        for(int i = 0; i < allObjectsCount; i++)
        {
            simInt obj = allObjectsBuf[i];
            simInt objType = simGetObjectType(obj);
            simInt visibleLayers = getVisibleLayers();
            simInt layers = getObjectLayers(obj);
            bool visible = visibleLayers & layers;
            if((objType == sim_object_shape_type && visible && !isWireframe(obj))
                    || objType == sim_object_camera_type)
                v.push_back(obj);
        }
        releaseBuffer(allObjectsBuf);
    }
}

void exportAllObjects(SScriptCallBack *p, const char *cmd, exportAllObjects_in *in, exportAllObjects_out *out)
{
    auto model = getModel(in->handle);
    exportObjects_in args;
    args.handle = in->handle;
    getAllObjects(args.objectHandles);
    if(args.objectHandles.empty()) return;
    exportObjects_out ret;
    exportObjects(p, &args, &ret);
}

void exportSelectedObjects(SScriptCallBack *p, const char *cmd, exportSelectedObjects_in *in, exportSelectedObjects_out *out)
{
    auto model = getModel(in->handle);
    exportObjects_in args;
    args.handle = in->handle;
    getObjectSelection(args.objectHandles);
    if(args.objectHandles.empty()) return;
    exportObjects_out ret;
    exportObjects(p, &args, &ret);
}

void exportObjects(SScriptCallBack *p, const char *cmd, exportObjects_in *in, exportObjects_out *out)
{
    auto model = getModel(in->handle);
    exportObject_in args;
    exportObject_out ret;
    for(simInt obj : in->objectHandles)
    {
        args.handle = in->handle;
        args.objectHandle = obj;
        exportObject(p, &args, &ret);
    }
}

void exportAnimation(SScriptCallBack *p, const char *cmd, exportAnimation_in *in, exportAnimation_out *out)
{
    auto model = getModel(in->handle);

    model->animations.resize(1);

    // create time buffer:
    int n = frames.size();
    simFloat t[n];
    for(int i = 0; i < n; i++) t[i] = frames[i].time;
    int bt = addBuffer(model, t, sizeof(simFloat) * n, "time");
    int vt = addBufferView(model, bt, sizeof(simFloat) * n, 0, "time");
    std::vector<double> tmin, tmax;
    minMaxVec(t, n, 1, tmin, tmax);
    int at = addAccessor(model, vt, 0, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_SCALAR, n, tmin, tmax, "time");

    for(simInt handle : handles)
    {
        if(nodeIndex.find(handle) == nodeIndex.end())
        {
            reportError((boost::format("Object with handle %d (%s) has no corresponding node") % handle % getObjectName(handle)).str());
            continue;
        }

        // for each object we need two channels: translation and rotation
        int ip = model->animations[0].channels.size();
        model->animations[0].channels.push_back({});
        model->animations[0].samplers.push_back({});
        int ir = model->animations[0].channels.size();
        model->animations[0].channels.push_back({});
        model->animations[0].samplers.push_back({});

        // create translation and rotation buffers:
        std::string name = getObjectName(handle);
        simFloat p[n * 3], r[n * 4];
        int hi = handleIndex[handle];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < 3; j++)
                p[3 * i + j] = frames[i].poses[hi].position[j];
            for(int j = 0; j < 4; j++)
                r[4 * i + j] = frames[i].poses[hi].orientation[j];
        }

        int bp = addBuffer(model, p, sizeof(simFloat) * n * 3, name + " position");
        int vp = addBufferView(model, bp, sizeof(simFloat) * n * 3, 0, name + " position");
        std::vector<double> pmin, pmax;
        minMaxVec(p, n * 3, 3, pmin, pmax);
        int ap = addAccessor(model, vp, 0, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3, n, pmin, pmax, name + " position");

        int br = addBuffer(model, r, sizeof(simFloat) * n * 4, name + " rotation");
        int vr = addBufferView(model, br, sizeof(simFloat) * n * 4, 0, name + " rotation");
        std::vector<double> rmin, rmax;
        minMaxVec(r, n * 4, 4, rmin, rmax);
        int ar = addAccessor(model, vr, 0, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC4, n, rmin, rmax, name + " rotation");

        // create samplers & channels:
        model->animations[0].samplers[ip].interpolation = "STEP";
        model->animations[0].samplers[ip].input = at;
        model->animations[0].samplers[ip].output = ap;
        model->animations[0].channels[ip].sampler = ip;
        model->animations[0].channels[ip].target_node = nodeIndex[handle];
        model->animations[0].channels[ip].target_path = "translation";
        model->animations[0].samplers[ir].interpolation = "STEP";
        model->animations[0].samplers[ir].input = at;
        model->animations[0].samplers[ir].output = ar;
        model->animations[0].channels[ir].sampler = ir;
        model->animations[0].channels[ir].target_node = nodeIndex[handle];
        model->animations[0].channels[ir].target_path = "rotation";
    }
}

void initAnimationFrames()
{
    handles.clear();
    handleIndex.clear();
    frames.clear();
    getAllObjects(handles);
    for(int i = 0; i < handles.size(); i++)
        handleIndex[handles[i]] = i;
}

void readAnimationFrame()
{
    if(simGetSimulationState() != sim_simulation_advancing_running)
        return;

    simFloat t = simGetSimulationTime();
    if(!frames.empty() && t < frames.back().time + 0.001)
        return;

    frames.push_back({});
    frames.back().read(t, handles);
}

class Plugin : public sim::Plugin
{
public:
    void onStart()
    {
        if(!registerScriptStuff())
            throw std::runtime_error("failed to register script stuff");

        simSetModuleInfo(PLUGIN_NAME, 0, "glTF support", 0);
        simSetModuleInfo(PLUGIN_NAME, 1, BUILD_DATE, 0);
    }

    void onSimulationAboutToStart()
    {
        initAnimationFrames();
    }

    void onInstancePass(const sim::InstancePassFlags &flags)
    {
        readAnimationFrame();
    }
};

SIM_PLUGIN(PLUGIN_NAME, PLUGIN_VERSION, Plugin)
