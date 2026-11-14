#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cctype>
#include <set>
#include <limits>
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

// ---------------------------------------------------------------------------
// Markers
//
// CoppeliaSim markers are not regular meshes but primitive descriptions:
// points, lines, tubes, one-point primitives (spheres, squares, ...), axes
// and custom meshes. glTF has no marker primitive, so every marker is
// translated into the closest glTF representation. Per-item colours are
// written as COLOR_0 vertex colours, including their alpha channel, and the
// notShaded flag becomes KHR_materials_unlit.
// ---------------------------------------------------------------------------

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Vec3 operator+(const Vec3 &a, const Vec3 &b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3 &a, const Vec3 &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(const Vec3 &a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator*(double s, const Vec3 &a) { return a * s; }

inline double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline Vec3 cross(const Vec3 &a, const Vec3 &b)
{
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

inline double length(const Vec3 &v) { return std::sqrt(dot(v, v)); }

inline Vec3 normalize(const Vec3 &v, const Vec3 &fallback = {0.0, 0.0, 1.0})
{
    const double l = length(v);
    return (l > 1e-12) ? v * (1.0 / l) : fallback;
}

// (qx, qy, qz, qw), as returned by the packedQuaternions marker property.
struct MarkerItemTransform
{
    Vec3 position;
    std::array<double, 4> orientation = {0.0, 0.0, 0.0, 1.0}; // qx, qy, qz, qw
    Vec3 scale = {1.0, 1.0, 1.0};
};

struct MarkerData
{
    int type = -1;
    std::string typeString;
    bool overlay = false;
    bool notShaded = false;
    std::vector<float> points;        // packedPoints: 3 floats per point
    std::vector<float> quaternions;   // packedQuaternions: 4 floats per item
    std::vector<float> sizes;         // packedSizes: 3 floats per item
    std::vector<unsigned char> colors; // colors: 4 bytes (RGBA) per point
    // custom markers only:
    std::vector<double> vertices;     // 3 x N matrix, row-major
    std::vector<int> indices;
    std::vector<double> normals;      // 3 x N matrix, row-major
    int vertexCount = 0;
    uint64_t meshHash = 0;            // custom mesh identity
    std::vector<int64_t> ids;         // per-item ids, unique within the marker
};

struct MarkerMeshBuilder
{
    int mode = TINYGLTF_MODE_TRIANGLES;
    std::vector<gltfFloat> positions;
    std::vector<gltfFloat> normals;
    std::vector<gltfFloat> colors;

    void add(const Vec3 &p, const Vec3 &n, const std::array<double, 4> &c)
    {
        positions.push_back(static_cast<gltfFloat>(p.x));
        positions.push_back(static_cast<gltfFloat>(p.y));
        positions.push_back(static_cast<gltfFloat>(p.z));
        normals.push_back(static_cast<gltfFloat>(n.x));
        normals.push_back(static_cast<gltfFloat>(n.y));
        normals.push_back(static_cast<gltfFloat>(n.z));
        colors.push_back(static_cast<gltfFloat>(c[0]));
        colors.push_back(static_cast<gltfFloat>(c[1]));
        colors.push_back(static_cast<gltfFloat>(c[2]));
        colors.push_back(static_cast<gltfFloat>(c[3]));
    }

    // Adds a triangle with a winding that agrees with its vertex normals, so
    // that the marker materials can stay single-sided (markers cull back
    // faces) regardless of how the primitive was generated.
    void addTriangle(const Vec3 &a, const Vec3 &b, const Vec3 &c,
                     const Vec3 &na, const Vec3 &nb, const Vec3 &nc,
                     const std::array<double, 4> &color)
    {
        if(dot(cross(b - a, c - a), na + nb + nc) < 0.0)
        {
            add(a, na, color);
            add(c, nc, color);
            add(b, nb, color);
        }
        else
        {
            add(a, na, color);
            add(b, nb, color);
            add(c, nc, color);
        }
    }

    // Preserves the caller's winding. Used for custom meshes, whose indices
    // already encode the intended front faces.
    void addTriangleRaw(const Vec3 &a, const Vec3 &b, const Vec3 &c,
                        const Vec3 &na, const Vec3 &nb, const Vec3 &nc,
                        const std::array<double, 4> &color)
    {
        add(a, na, color);
        add(b, nb, color);
        add(c, nc, color);
    }
};

// Stable identity for one marker item, derived from its data (points, colours,
// sizes and quaternions), as suggested by the marker API's lack of item ids.
struct MarkerItemDesc
{
    int64_t id = 0;
    int itemIndex = -1;
};

struct MarkerItemState
{
    int nodeIndex = -1;
    std::map<size_t, bool> visible; // animation time index -> item present
    std::map<size_t, MarkerItemTransform> transforms;
};

struct MarkerMeshKey
{
    uint64_t color = 0;
    int flags = 0;

    bool operator<(const MarkerMeshKey &other) const
    {
        if(flags != other.flags)
            return flags < other.flags;
        return color < other.color;
    }
};

struct MarkerGeometry
{
    int positionAccessor = -1;
    int normalAccessor = -1;
    int vertexCount = 0;
    int mode = TINYGLTF_MODE_TRIANGLES;
};

struct MarkerState
{
    int nodeIndex = -1;
    std::map<int64_t, MarkerItemState> items;
    std::map<MarkerMeshKey, int> itemMeshes; // custom marker geometry
    MarkerGeometry customGeometry;
    uint64_t customGeometryHash = 0;
};

struct MarkerUidBuilder
{
    uint64_t value = 1469598103934665603ULL;

    void addBytes(const void *data, size_t size)
    {
        const unsigned char *bytes = static_cast<const unsigned char *>(data);
        for(size_t i = 0; i < size; i++)
        {
            value ^= bytes[i];
            value *= 1099511628211ULL;
        }
    }

    template<typename T>
    void add(const T &v)
    {
        addBytes(&v, sizeof(T));
    }
};

struct BufferSlice
{
    int buffer = 0;
    size_t offset = 0;
    size_t length = 0;
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

    std::string getObjectName(_handle_t handle)
    {
        return sim::getObjectAlias(handle, 4);
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
        int visible = 0;
        if(simGetBoolProperty(handle, "visible", &visible) == sim_propertyret_ok)
            return visible != 0;
        sim::clearLastError();
        return false;
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

    bool isMarker(_handle_t handle)
    {
        std::string objType = sim::getStringProperty(handle, "objectType");
        return objType == "marker";
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
        p->visible = !(isShape(handle) || isMarker(handle)) || (isVisible(handle) && !isWireframe(handle));
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
        markerStates.clear();
        markerMeshCaches.clear();
        markerGeometryCache.clear();
        markerMaterials.clear();
        bufferSlices.clear();
        sharedBufferIndex = -1;

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

    static bool removeJsonIntProperty(std::string &json, const char *name)
    {
        const std::string key = std::string("\"") + name + "\"";
        bool changed = false;
        size_t search = 0;
        while(true)
        {
            const size_t pos = json.find(key, search);
            if(pos == std::string::npos)
                break;

            const size_t colon = json.find(':', pos + key.size());
            if(colon == std::string::npos)
                break;
            const size_t valueStart = json.find_first_of("-0123456789", colon + 1);
            if(valueStart == std::string::npos)
                break;
            size_t valueEnd = json.find_first_not_of("0123456789-+", valueStart);
            if(valueEnd == std::string::npos)
                valueEnd = json.size();

            size_t removeStart = pos;
            size_t removeEnd = valueEnd;

            size_t before = pos;
            while(before > 0 && std::isspace(static_cast<unsigned char>(json[before - 1])))
                before--;
            if(before > 0 && json[before - 1] == ',')
            {
                removeStart = before - 1;
            }
            else
            {
                size_t after = valueEnd;
                while(after < json.size() && std::isspace(static_cast<unsigned char>(json[after])))
                    after++;
                if(after < json.size() && json[after] == ',')
                    removeEnd = after + 1;
            }

            json.erase(removeStart, removeEnd - removeStart);
            changed = true;
            search = removeStart;
        }
        return changed;
    }

    static uint32_t readU32(const std::string &data, size_t offset)
    {
        return static_cast<uint32_t>(static_cast<unsigned char>(data[offset + 0])) |
               (static_cast<uint32_t>(static_cast<unsigned char>(data[offset + 1])) << 8) |
               (static_cast<uint32_t>(static_cast<unsigned char>(data[offset + 2])) << 16) |
               (static_cast<uint32_t>(static_cast<unsigned char>(data[offset + 3])) << 24);
    }

    static void appendU32(std::string &data, uint32_t value)
    {
        data.push_back(static_cast<char>(value & 0xff));
        data.push_back(static_cast<char>((value >> 8) & 0xff));
        data.push_back(static_cast<char>((value >> 16) & 0xff));
        data.push_back(static_cast<char>((value >> 24) & 0xff));
    }

    static bool sanitizeGlb(std::string &data)
    {
        constexpr uint32_t gltfMagic = 0x46546C67; // "glTF"
        constexpr uint32_t jsonChunk = 0x4E4F534A; // "JSON"

        if(data.size() < 20 || readU32(data, 0) != gltfMagic || readU32(data, 4) != 2)
            return false;

        std::string result(12, '\0');
        bool changed = false;
        size_t offset = 12;
        while(offset + 8 <= data.size())
        {
            const uint32_t chunkLength = readU32(data, offset);
            const uint32_t chunkType = readU32(data, offset + 4);
            if(offset + 8 + chunkLength > data.size())
                return false;

            std::string payload = data.substr(offset + 8, chunkLength);
            if(chunkType == jsonChunk && removeJsonIntProperty(payload, "wrapR"))
                changed = true;

            const uint32_t paddedLength = (static_cast<uint32_t>(payload.size()) + 3u) & ~3u;
            appendU32(result, paddedLength);
            appendU32(result, chunkType);
            result.append(payload);
            result.append(paddedLength - payload.size(), chunkType == jsonChunk ? ' ' : '\0');
            offset += 8 + chunkLength;
        }

        if(!changed || offset != data.size())
            return false;

        result[0] = data[0];
        result[1] = data[1];
        result[2] = data[2];
        result[3] = data[3];
        result[4] = data[4];
        result[5] = data[5];
        result[6] = data[6];
        result[7] = data[7];
        const uint32_t totalLength = static_cast<uint32_t>(result.size());
        result[8] = static_cast<char>(totalLength & 0xff);
        result[9] = static_cast<char>((totalLength >> 8) & 0xff);
        result[10] = static_cast<char>((totalLength >> 16) & 0xff);
        result[11] = static_cast<char>((totalLength >> 24) & 0xff);

        data.swap(result);
        return true;
    }

    static void sanitizeGltfFile(const std::string &filepath, bool binary)
    {
        std::ifstream in(filepath, std::ios::binary);
        if(!in)
            return;
        std::stringstream ss;
        ss << in.rdbuf();
        std::string data = ss.str();
        in.close();

        const bool changed = binary ? sanitizeGlb(data) : removeJsonIntProperty(data, "wrapR");
        if(!changed)
            return;

        std::ofstream out(filepath, std::ios::binary | std::ios::trunc);
        out.write(data.data(), static_cast<std::streamsize>(data.size()));
    }

    void saveASCII(saveASCII_in *in, saveASCII_out *out)
    {
        out->result = gltf.WriteGltfSceneToFile(&model, in->filepath, true, true, true, false);
        if(out->result)
            sanitizeGltfFile(in->filepath, false);
    }

    void saveBinary(saveBinary_in *in, saveBinary_out *out)
    {
        out->result = gltf.WriteGltfSceneToFile(&model, in->filepath, true, true, true, true);
        if(out->result)
            sanitizeGltfFile(in->filepath, true);
    }

    void serialize(serialize_in *in, serialize_out *out)
    {
        std::stringstream ss;
        gltf.WriteGltfSceneToStream(&model, ss, true, false);
        std::string json = ss.str();
        removeJsonIntProperty(json, "wrapR");
        out->json = json;
    }

    int addBuffer(const void *buffer, int size, const std::string &name)
    {
        if(sharedBufferIndex < 0)
        {
            // Put all generated data into a single buffer: glTF embeds each
            // buffer separately, and one big buffer avoids enormous base64
            // overhead in marker-heavy files.
            model.buffers.push_back({});
            sharedBufferIndex = static_cast<int>(model.buffers.size()) - 1;
            model.buffers[sharedBufferIndex].name = "CoppeliaSim glTF buffer";
        }

        auto &data = model.buffers[sharedBufferIndex].data;
        // glTF requires accessor byte offsets to be a multiple of the
        // component type size (4 bytes for the types used here). Also keep
        // bufferView offsets 4-byte aligned.
        const size_t offset = (data.size() + 3u) & ~static_cast<size_t>(3u);
        if(offset > data.size())
            data.resize(offset, 0);
        data.resize(offset + static_cast<size_t>(size));
        sim::addLog(sim_verbosity_debug, "addBuffer: added buffer '%s' %s", name, buf2str(buffer, size));
        std::memcpy(data.data() + offset, buffer, size);
        bufferSlices.push_back({sharedBufferIndex, offset, static_cast<size_t>(size)});
        return static_cast<int>(bufferSlices.size()) - 1;
    }

    int addBufferView(int buffer, int byteLength, int byteOffset, const std::string &name)
    {
        const BufferSlice &slice = bufferSlices.at(static_cast<size_t>(buffer));
        int i = model.bufferViews.size();
        model.bufferViews.push_back({});
        model.bufferViews[i].buffer = slice.buffer;
        model.bufferViews[i].byteLength = byteLength;
        model.bufferViews[i].byteOffset = slice.offset + static_cast<size_t>(byteOffset);
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
        // glTF only allows JPEG and PNG images.
        if(exportTextureFormat == simgltf_texture_format_jpg)
            return raw2jpg(reinterpret_cast<const unsigned char *>(imgdata), res, 4, 100);
        return raw2png(reinterpret_cast<const unsigned char *>(imgdata), res, 4, res[0]*4);
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
        const bool jpeg = exportTextureFormat == simgltf_texture_format_jpg;
        const std::string formatName = jpeg ? "JPG" : "PNG";
        std::string name = (boost::format("texture image %d [%s] (%dx%d, %s %d bytes)") % id % objname % res[0] % res[1] % formatName % buf.size()).str();

        int b = addBuffer(buf.data(), buf.size(), name);

        int v = addBufferView(b, buf.size(), 0, name);

        int i = model.images.size();
        model.images.push_back({});
        model.images[i].bufferView = v;
        model.images[i].mimeType = jpeg ? "image/jpeg" : "image/png";
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
        // info.transparency is 0.0 for opaque, otherwise it is the
        // transparency amount. glTF baseColorFactor alpha is opacity.
        const double transparency = std::max(0.0, std::min(1.0, static_cast<double>(info.transparency)));
        const double alpha = 1.0 - transparency;
        model.materials[mat].pbrMetallicRoughness.baseColorFactor = {diffuse[0], diffuse[1], diffuse[2], alpha};
        model.materials[mat].pbrMetallicRoughness.metallicFactor = 0.1;
        model.materials[mat].pbrMetallicRoughness.roughnessFactor = 0.5;
        model.materials[mat].doubleSided = ((info.options & 1) == 0);
        if(transparency > 0.0)
            model.materials[mat].alphaMode = "BLEND";

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

    // -----------------------------------------------------------------------
    // Marker support
    //
    // Marker data is fetched with the new simLib-2 property API. A marker is
    // a scene object (it has a pose in the hierarchy) whose items are stored
    // in its local frame. Every item becomes a child node with its own mesh,
    // so that animation can show/hide/move items independently of each other.
    // -----------------------------------------------------------------------

    std::string readMarkerStringProperty(_handle_t handle, const char *name)
    {
        char *value = nullptr;
        std::string result;
        if(simGetStringProperty(handle, name, &value) == sim_propertyret_ok && value != nullptr)
            result = value;
        else
            sim::clearLastError();
        if(value != nullptr)
            simReleaseBuffer(value);
        return result;
    }

    bool readMarkerBoolProperty(_handle_t handle, const char *name)
    {
        int value = 0;
        if(simGetBoolProperty(handle, name, &value) == sim_propertyret_ok)
            return value != 0;
        sim::clearLastError();
        return false;
    }

    static std::vector<unsigned char> readMarkerBufferProperty(_handle_t handle, const char *name)
    {
        char *value = nullptr;
        int length = 0;
        std::vector<unsigned char> result;
        if(simGetBufferProperty(handle, name, &value, &length) == sim_propertyret_ok &&
           value != nullptr && length > 0)
        {
            result.assign(reinterpret_cast<unsigned char *>(value),
                          reinterpret_cast<unsigned char *>(value) + length);
        }
        else
            sim::clearLastError();
        if(value != nullptr)
            simReleaseBuffer(value);
        return result;
    }

    std::vector<int> readMarkerIntArrayProperty(_handle_t handle, const char *name)
    {
        int *value = nullptr;
        int length = 0;
        std::vector<int> result;
        if(simGetIntArrayProperty(handle, name, &value, &length) == sim_propertyret_ok &&
           value != nullptr && length > 0)
            result.assign(value, value + length);
        else
            sim::clearLastError();
        if(value != nullptr)
            simReleaseBuffer(value);
        return result;
    }

    static std::vector<int64_t> readMarkerLongArrayProperty(_handle_t handle, const char *name)
    {
        long long int *value = nullptr;
        int length = 0;
        std::vector<int64_t> result;
        if(simGetLongArrayProperty(handle, name, &value, &length) == sim_propertyret_ok &&
           value != nullptr && length > 0)
            result.assign(value, value + length);
        else
            sim::clearLastError();
        if(value != nullptr)
            simReleaseBuffer(value);
        return result;
    }

    static std::vector<double> readMarkerFloatArrayProperty(_handle_t handle, const char *name)
    {
        double *value = nullptr;
        int length = 0;
        std::vector<double> result;
        if(simGetFloatArrayProperty(handle, name, &value, &length) == sim_propertyret_ok &&
           value != nullptr && length > 0)
            result.assign(value, value + length);
        else
            sim::clearLastError();
        if(value != nullptr)
            simReleaseBuffer(value);
        return result;
    }

    static std::vector<double> readMarkerDoubleBufferProperty(_handle_t handle, const char *name)
    {
        const std::vector<unsigned char> bytes = readMarkerBufferProperty(handle, name);
        if(bytes.empty() || bytes.size() % sizeof(double) != 0)
            return {};
        std::vector<double> result(bytes.size() / sizeof(double));
        std::memcpy(result.data(), bytes.data(), bytes.size());
        return result;
    }

    // Accept flat interleaved data (x,y,z per point/normal) irrespective of
    // how the property is exposed. The matrix form used by current builds is
    // row-major 3 x N, i.e. all x values, then all y values, then all z values.
    static std::vector<double> markerMatrixToFlat(const std::vector<double> &matrix, int rows, int cols)
    {
        if(rows == 3 && cols > 0)
        {
            std::vector<double> flat(static_cast<size_t>(cols) * 3);
            for(int i = 0; i < cols; i++)
            {
                flat[3 * i + 0] = matrix[i];
                flat[3 * i + 1] = matrix[cols + i];
                flat[3 * i + 2] = matrix[2 * cols + i];
            }
            return flat;
        }
        if((rows > 0 && cols == 3) || (!matrix.empty() && matrix.size() % 3 == 0))
            return matrix;
        return {};
    }

    static std::vector<double> readMarkerPointData(_handle_t handle, const char *name)
    {
        // Prefer the matrix form: it carries the layout explicitly. Current
        // builds expose all x values, then all y values, then all z values.
        int rows = 0;
        int cols = 0;
        const std::vector<double> matrix = readMarkerMatrixProperty(handle, name, rows, cols);
        if(!matrix.empty())
            return markerMatrixToFlat(matrix, rows, cols);

        std::vector<double> flat = readMarkerFloatArrayProperty(handle, name);
        if(!flat.empty())
            return flat;
        return readMarkerDoubleBufferProperty(handle, name);
    }

    static std::vector<double> readMarkerMatrixProperty(_handle_t handle, const char *name, int &rows, int &cols)
    {
        double *value = nullptr;
        rows = 0;
        cols = 0;
        std::vector<double> result;
        if(simGetMatrixProperty(handle, name, &value, &rows, &cols) == sim_propertyret_ok &&
           value != nullptr && rows > 0 && cols > 0)
        {
            result.assign(value, value + static_cast<size_t>(rows) * static_cast<size_t>(cols));
        }
        else
            sim::clearLastError();
        if(value != nullptr)
            simReleaseBuffer(value);
        return result;
    }

    static std::vector<float> unpackFloats(const std::vector<unsigned char> &bytes)
    {
        std::vector<float> result(bytes.size() / sizeof(float));
        if(!result.empty())
            std::memcpy(result.data(), bytes.data(), result.size() * sizeof(float));
        return result;
    }

    static int markerTypeFromString(const std::string &type)
    {
        if(type == "points")     return sim_markertype_points;
        if(type == "lines")      return sim_markertype_lines;
        if(type == "triangles")  return sim_markertype_triangles;
        if(type == "spheres")    return sim_markertype_spheres;
        if(type == "squares")    return sim_markertype_squares;
        if(type == "discs")      return sim_markertype_discs;
        if(type == "cubes")      return sim_markertype_cubes;
        if(type == "cylinders")  return sim_markertype_cylinders;
        if(type == "custom")     return sim_markertype_custom;
        if(type == "tubes")      return sim_markertype_tubes;
        if(type == "axes")       return sim_markertype_axes;
        return -1;
    }

    MarkerData readMarkerData(_handle_t handle)
    {
        MarkerData marker;
        marker.typeString = readMarkerStringProperty(handle, "marker.type");
        marker.type = markerTypeFromString(marker.typeString);
        marker.overlay = readMarkerBoolProperty(handle, "overlay");
        marker.notShaded = readMarkerBoolProperty(handle, "notShaded");
        marker.points = unpackFloats(readMarkerBufferProperty(handle, "packedPoints"));
        marker.quaternions = unpackFloats(readMarkerBufferProperty(handle, "packedQuaternions"));
        marker.sizes = unpackFloats(readMarkerBufferProperty(handle, "packedSizes"));
        marker.colors = readMarkerBufferProperty(handle, "colors");
        marker.ids = readMarkerLongArrayProperty(handle, "ids");
        if(marker.type == sim_markertype_custom)
        {
            marker.vertices = readMarkerPointData(handle, "vertices");
            marker.vertexCount = static_cast<int>(marker.vertices.size() / 3);
            marker.indices = readMarkerIntArrayProperty(handle, "indices");
            marker.normals = readMarkerPointData(handle, "normals");

            MarkerUidBuilder meshHash;
            if(!marker.vertices.empty())
                meshHash.addBytes(marker.vertices.data(), marker.vertices.size() * sizeof(double));
            if(!marker.indices.empty())
                meshHash.addBytes(marker.indices.data(), marker.indices.size() * sizeof(int));
            if(!marker.normals.empty())
                meshHash.addBytes(marker.normals.data(), marker.normals.size() * sizeof(double));
            marker.meshHash = meshHash.value;
        }
        return marker;
    }

    static std::array<double, 4> markerColorAt(const MarkerData &marker, int pointIndex)
    {
        const std::array<double, 4> fallback = {1.0, 1.0, 1.0, 1.0};
        const size_t base = static_cast<size_t>(pointIndex) * 4;
        if(base + 4 <= marker.colors.size())
        {
            const unsigned char *c = marker.colors.data() + base;
            return {c[0] / 255.0, c[1] / 255.0, c[2] / 255.0, c[3] / 255.0};
        }
        return fallback;
    }

    static Vec3 markerPointAt(const MarkerData &marker, int index)
    {
        const size_t base = static_cast<size_t>(index) * 3;
        return {marker.points[base + 0], marker.points[base + 1], marker.points[base + 2]};
    }

    static MarkerItemTransform markerItemTransform(const MarkerData &marker, int itemIndex, double defaultSize)
    {
        MarkerItemTransform transform;
        if(static_cast<size_t>(itemIndex) * 3 + 3 <= marker.points.size())
            transform.position = markerPointAt(marker, itemIndex);

        const size_t qbase = static_cast<size_t>(itemIndex) * 4;
        if(qbase + 4 <= marker.quaternions.size())
        {
            transform.orientation = {marker.quaternions[qbase + 0],
                                     marker.quaternions[qbase + 1],
                                     marker.quaternions[qbase + 2],
                                     marker.quaternions[qbase + 3]};
        }

        const size_t sbase = static_cast<size_t>(itemIndex) * 3;
        if(sbase + 3 <= marker.sizes.size())
        {
            transform.scale = {marker.sizes[sbase + 0],
                               marker.sizes[sbase + 1],
                               marker.sizes[sbase + 2]};
        }
        else
        {
            transform.scale = {defaultSize, defaultSize, defaultSize};
        }
        return transform;
    }

    static tinygltf::Value markerExtras(const MarkerData &marker)
    {
        // glTF has no native marker concept and no way to disable depth testing
        // per material, so the marker properties are kept in extras for tools
        // that consume this file. notShaded is also expressed with
        // KHR_materials_unlit (see addMarkerMaterial).
        tinygltf::Value::Object extras;
        extras["markerType"] = tinygltf::Value(marker.typeString);
        extras["overlay"] = tinygltf::Value(marker.overlay);
        extras["notShaded"] = tinygltf::Value(marker.notShaded);
        return tinygltf::Value(std::move(extras));
    }

    static int markerPointsPerItem(int type)
    {
        if(type == sim_markertype_lines || type == sim_markertype_tubes)
            return 2;
        if(type == sim_markertype_triangles)
            return 3;
        return 1;
    }

    // Points and spheres share the unit-sphere geometry; lines and tubes share
    // the unit-tube geometry.
    static int markerGeometryType(int type)
    {
        if(type == sim_markertype_points)
            return sim_markertype_spheres;
        if(type == sim_markertype_lines)
            return sim_markertype_tubes;
        return type;
    }

    static double markerDefaultItemSize(int type)
    {
        // Custom items always carry explicit sizes; axes too in practice, but
        // fall back to the marker's default item size.
        return (type == sim_markertype_custom) ? 1.0 : 0.005;
    }

    static uint64_t markerColorKey(const std::array<double, 4> &color)
    {
        uint64_t key = 0;
        for(int i = 0; i < 4; i++)
            key |= static_cast<uint64_t>(std::lround(color[i] * 255.0)) << (8 * i);
        return key;
    }

    static std::vector<MarkerItemDesc> buildMarkerItemDescs(const MarkerData &marker)
    {
        std::vector<MarkerItemDesc> items;
        if(marker.type < 0 || marker.type == sim_markertype_triangles || marker.points.empty())
            return items;

        const int ppi = markerPointsPerItem(marker.type);
        if(ppi <= 0)
            return items;

        const int itemCount = static_cast<int>(marker.points.size() / 3) / ppi;
        items.reserve(itemCount);
        for(int i = 0; i < itemCount; i++)
        {
            const int64_t id = (i < static_cast<int>(marker.ids.size()))
                                   ? marker.ids[i]
                                   : std::numeric_limits<int64_t>::min() + i;
            items.push_back({id, i});
        }
        return items;
    }

    int getMarkerItemMaterial(bool notShaded, bool transparent)
    {
        const int key = (notShaded ? 1 : 0) | (transparent ? 2 : 0);
        auto it = markerMaterials.find(key);
        if(it != markerMaterials.end())
            return it->second;

        int mat = model.materials.size();
        model.materials.push_back({});
        model.materials[mat].name = "marker material";
        // Marker colours are per-vertex (COLOR_0); the material itself is white.
        model.materials[mat].pbrMetallicRoughness.baseColorFactor = {1.0, 1.0, 1.0, 1.0};
        model.materials[mat].pbrMetallicRoughness.metallicFactor = 0.0;
        model.materials[mat].pbrMetallicRoughness.roughnessFactor = 0.5;
        model.materials[mat].doubleSided = false;
        if(transparent)
        {
            // Vertex colours carry the marker alpha.
            model.materials[mat].alphaMode = "BLEND";
        }
        if(notShaded)
        {
            // "notShaded" means plain colour without any light influence.
            model.materials[mat].extensions["KHR_materials_unlit"] =
                tinygltf::Value(tinygltf::Value::Object{});
            if(std::find(model.extensionsUsed.begin(), model.extensionsUsed.end(),
                         "KHR_materials_unlit") == model.extensionsUsed.end())
                model.extensionsUsed.push_back("KHR_materials_unlit");
        }
        markerMaterials[key] = mat;
        return mat;
    }

    MarkerGeometry addMarkerGeometry(const MarkerMeshBuilder &data, const std::string &name)
    {
        MarkerGeometry geometry;
        geometry.mode = data.mode;
        geometry.vertexCount = static_cast<int>(data.positions.size() / 3);

        int bp = addBuffer(data.positions.data(), sizeof(gltfFloat) * data.positions.size(), name + " marker point");
        int bn = addBuffer(data.normals.data(), sizeof(gltfFloat) * data.normals.size(), name + " marker normal");
        int vp = addBufferView(bp, sizeof(gltfFloat) * data.positions.size(), 0, name + " marker point");
        int vn = addBufferView(bn, sizeof(gltfFloat) * data.normals.size(), 0, name + " marker normal");

        std::vector<double> pmin, pmax, nmin, nmax;
        minMaxVec(data.positions, 3, pmin, pmax);
        minMaxVec(data.normals, 3, nmin, nmax);

        geometry.positionAccessor = addAccessor(vp, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3,
                                                data.positions.size() / 3, pmin, pmax, name + " marker point");
        geometry.normalAccessor = addAccessor(vn, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC3,
                                              data.normals.size() / 3, nmin, nmax, name + " marker normal");
        return geometry;
    }

    int addMarkerColorAccessor(const std::vector<gltfFloat> &colors, const std::string &name)
    {
        int b = addBuffer(colors.data(), sizeof(gltfFloat) * colors.size(), name + " marker color");
        int v = addBufferView(b, sizeof(gltfFloat) * colors.size(), 0, name + " marker color");
        std::vector<double> cmin, cmax;
        minMaxVec(colors, 4, cmin, cmax);
        return addAccessor(v, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat, TINYGLTF_TYPE_VEC4,
                           colors.size() / 4, cmin, cmax, name + " marker color");
    }

    int addMarkerMesh(const std::string &name, const MarkerGeometry &geometry, int colorAccessor, int material)
    {
        int mesh = model.meshes.size();
        model.meshes.push_back({});
        model.meshes[mesh].name = name + " marker mesh";
        model.meshes[mesh].primitives.push_back({});
        model.meshes[mesh].primitives[0].attributes["POSITION"] = geometry.positionAccessor;
        model.meshes[mesh].primitives[0].attributes["NORMAL"] = geometry.normalAccessor;
        model.meshes[mesh].primitives[0].attributes["COLOR_0"] = colorAccessor;
        model.meshes[mesh].primitives[0].mode = geometry.mode;
        model.meshes[mesh].primitives[0].material = material;
        return mesh;
    }

    static void buildMarkerPointUnit(const std::array<double, 4> &color, MarkerMeshBuilder &mesh)
    {
        // Points are rendered as small spheres with radius 0.0025.
        buildMarkerOnePointUnit(sim_markertype_spheres, color, mesh);
    }

    static void buildMarkerLineUnit(const std::array<double, 4> &c0,
                                    const std::array<double, 4> &c1,
                                    MarkerMeshBuilder &mesh)
    {
        // Lines are rendered as tubes with radius 0.0025, with a colour
        // gradient from the first to the second endpoint.
        buildMarkerTubeUnit(c0, mesh);
        for(size_t i = 0; i < mesh.positions.size() / 3; i++)
        {
            if(mesh.positions[3 * i + 2] > 0.5f)
            {
                mesh.colors[4 * i + 0] = static_cast<gltfFloat>(c1[0]);
                mesh.colors[4 * i + 1] = static_cast<gltfFloat>(c1[1]);
                mesh.colors[4 * i + 2] = static_cast<gltfFloat>(c1[2]);
                mesh.colors[4 * i + 3] = static_cast<gltfFloat>(c1[3]);
            }
        }
    }

    static void buildMarkerTubeUnit(const std::array<double, 4> &color, MarkerMeshBuilder &mesh)
    {
        mesh.mode = TINYGLTF_MODE_TRIANGLES;
        const int sides = 8;
        const double pi = 3.14159265358979323846;
        const double radius = 0.5;

        for(int s = 0; s < sides; s++)
        {
            const double a0 = 2.0 * pi * s / sides;
            const double a1 = 2.0 * pi * (s + 1) / sides;
            const Vec3 n0{std::cos(a0), std::sin(a0), 0.0};
            const Vec3 n1{std::cos(a1), std::sin(a1), 0.0};
            const Vec3 r0 = n0 * radius;
            const Vec3 r1 = n1 * radius;

            mesh.addTriangle({r0.x, r0.y, 0.0}, {r1.x, r1.y, 0.0}, {r1.x, r1.y, 1.0},
                             n0, n1, n1, color);
            mesh.addTriangle({r0.x, r0.y, 0.0}, {r1.x, r1.y, 1.0}, {r0.x, r0.y, 1.0},
                             n0, n1, n0, color);
        }
    }

    static void buildMarkerAxesUnit(MarkerMeshBuilder &mesh)
    {
        mesh.mode = TINYGLTF_MODE_LINE;
        const std::array<double, 4> red = {1.0, 0.0, 0.0, 1.0};
        const std::array<double, 4> green = {0.0, 1.0, 0.0, 1.0};
        const std::array<double, 4> blue = {0.0, 0.0, 1.0, 1.0};
        const Vec3 origin{0.0, 0.0, 0.0};
        const Vec3 normal{0.0, 0.0, 1.0};

        mesh.add(origin, normal, red);
        mesh.add({1.0, 0.0, 0.0}, normal, red);
        mesh.add(origin, normal, green);
        mesh.add({0.0, 1.0, 0.0}, normal, green);
        mesh.add(origin, normal, blue);
        mesh.add({0.0, 0.0, 1.0}, normal, blue);
    }

    static void buildMarkerOnePointUnit(int type, const std::array<double, 4> &color,
                                        MarkerMeshBuilder &mesh)
    {
        mesh.mode = TINYGLTF_MODE_TRIANGLES;
        const int slices = 12;
        const int stacks = 8;
        const double pi = 3.14159265358979323846;
        const Vec3 localUp{0.0, 0.0, 1.0};
        const auto triangle = [&](const Vec3 &a, const Vec3 &na,
                                  const Vec3 &b, const Vec3 &nb,
                                  const Vec3 &c, const Vec3 &nc)
        {
            mesh.addTriangle(a, b, c, na, nb, nc, color);
        };

        if(type == sim_markertype_squares)
        {
            triangle({-0.5, -0.5, 0.0}, localUp, {0.5, -0.5, 0.0}, localUp, {0.5, 0.5, 0.0}, localUp);
            triangle({-0.5, -0.5, 0.0}, localUp, {0.5, 0.5, 0.0}, localUp, {-0.5, 0.5, 0.0}, localUp);
        }
        else if(type == sim_markertype_discs)
        {
            for(int k = 0; k < slices; k++)
            {
                const double a0 = 2.0 * pi * k / slices;
                const double a1 = 2.0 * pi * (k + 1) / slices;
                triangle({0.0, 0.0, 0.0}, localUp,
                         {0.5 * std::cos(a0), 0.5 * std::sin(a0), 0.0}, localUp,
                         {0.5 * std::cos(a1), 0.5 * std::sin(a1), 0.0}, localUp);
            }
        }
        else if(type == sim_markertype_spheres)
        {
            for(int k = 0; k < stacks; k++)
            {
                const double lat0 = pi * (-0.5 + static_cast<double>(k) / stacks);
                const double lat1 = pi * (-0.5 + static_cast<double>(k + 1) / stacks);
                const double z0 = 0.5 * std::sin(lat0);
                const double r0 = 0.5 * std::cos(lat0);
                const double z1 = 0.5 * std::sin(lat1);
                const double r1 = 0.5 * std::cos(lat1);
                for(int j = 0; j < slices; j++)
                {
                    const double lng0 = 2.0 * pi * j / slices;
                    const double lng1 = 2.0 * pi * (j + 1) / slices;
                    const Vec3 a{r0 * std::cos(lng0), r0 * std::sin(lng0), z0};
                    const Vec3 b{r0 * std::cos(lng1), r0 * std::sin(lng1), z0};
                    const Vec3 c{r1 * std::cos(lng0), r1 * std::sin(lng0), z1};
                    const Vec3 d{r1 * std::cos(lng1), r1 * std::sin(lng1), z1};
                    triangle(a, a * 2.0, c, c * 2.0, d, d * 2.0);
                    triangle(a, a * 2.0, d, d * 2.0, b, b * 2.0);
                }
            }
        }
        else if(type == sim_markertype_cubes)
        {
            const Vec3 corners[8] = {
                {-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {-0.5, 0.5, -0.5},
                {-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}};
            const int faces[6][4] = {{0, 1, 2, 3}, {5, 4, 7, 6}, {4, 0, 3, 7},
                                     {1, 5, 6, 2}, {4, 5, 1, 0}, {3, 2, 6, 7}};
            const Vec3 faceNormals[6] = {{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0},
                                         {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                         {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}};
            for(int f = 0; f < 6; f++)
            {
                const Vec3 &n = faceNormals[f];
                const Vec3 &a = corners[faces[f][0]];
                const Vec3 &b = corners[faces[f][1]];
                const Vec3 &c = corners[faces[f][2]];
                const Vec3 &d = corners[faces[f][3]];
                triangle(a, n, b, n, c, n);
                triangle(a, n, c, n, d, n);
            }
        }
        else if(type == sim_markertype_cylinders)
        {
            for(int k = 0; k < slices; k++)
            {
                const double a0 = 2.0 * pi * k / slices;
                const double a1 = 2.0 * pi * (k + 1) / slices;
                const Vec3 b0{0.5 * std::cos(a0), 0.5 * std::sin(a0), 0.0};
                const Vec3 b1{0.5 * std::cos(a1), 0.5 * std::sin(a1), 0.0};
                const Vec3 n0{std::cos(a0), std::sin(a0), 0.0};
                const Vec3 n1{std::cos(a1), std::sin(a1), 0.0};

                triangle({b0.x, b0.y, 0.5}, n0, {b1.x, b1.y, 0.5}, n1, {b1.x, b1.y, -0.5}, n1);
                triangle({b0.x, b0.y, 0.5}, n0, {b1.x, b1.y, -0.5}, n1, {b0.x, b0.y, -0.5}, n0);
                // Caps.
                triangle({0.0, 0.0, 0.5}, localUp, {b1.x, b1.y, 0.5}, localUp, {b0.x, b0.y, 0.5}, localUp);
                triangle({0.0, 0.0, -0.5}, localUp * -1.0, {b0.x, b0.y, -0.5}, localUp * -1.0,
                         {b1.x, b1.y, -0.5}, localUp * -1.0);
            }
        }
    }

    static void buildMarkerCustomUnit(const MarkerData &marker,
                                      const std::array<double, 4> &color,
                                      MarkerMeshBuilder &mesh)
    {
        mesh.mode = TINYGLTF_MODE_TRIANGLES;
        const int vertexCount = marker.vertexCount;
        if(vertexCount <= 0 || marker.vertices.size() < static_cast<size_t>(vertexCount) * 3)
            return;

        // Vertices/normals are flat interleaved lists: x,y,z per point.
        const auto vertexAt = [&](int index) -> Vec3
        {
            return {marker.vertices[3 * index + 0],
                    marker.vertices[3 * index + 1],
                    marker.vertices[3 * index + 2]};
        };

        std::vector<int> generatedIndices;
        const std::vector<int> *indices = &marker.indices;
        if(marker.indices.empty())
        {
            generatedIndices.resize(vertexCount);
            for(int i = 0; i < vertexCount; i++)
                generatedIndices[i] = i;
            indices = &generatedIndices;
        }

        const int indexCount = static_cast<int>(indices->size());
        const int normalCount = static_cast<int>(marker.normals.size() / 3);
        // Custom markers normally provide one normal per triangle corner
        // (3 * triangle count). Accept per-vertex normals too, for older data.
        const bool cornerNormals = normalCount == indexCount;
        const bool vertexNormals = !cornerNormals && normalCount >= vertexCount;
        const auto normalAt = [&](int index) -> Vec3
        {
            return {marker.normals[3 * index + 0],
                    marker.normals[3 * index + 1],
                    marker.normals[3 * index + 2]};
        };

        for(int k = 0; k + 2 < indexCount; k += 3)
        {
            Vec3 p[3];
            Vec3 n[3];
            bool valid = true;
            for(int j = 0; j < 3; j++)
            {
                const int index = (*indices)[k + j];
                if(index < 0 || index >= vertexCount)
                {
                    valid = false;
                    break;
                }
                p[j] = vertexAt(index);
                if(cornerNormals)
                    n[j] = normalAt(k + j);
                else if(vertexNormals)
                    n[j] = normalAt(index);
            }
            if(!valid)
                continue;

            if(!cornerNormals && !vertexNormals)
            {
                const Vec3 faceNormal = normalize(cross(p[1] - p[0], p[2] - p[0]));
                n[0] = n[1] = n[2] = faceNormal;
            }
            mesh.addTriangleRaw(p[0], p[1], p[2], n[0], n[1], n[2], color);
        }
    }

    static std::array<double, 4> markerRotationFromZ(const Vec3 &direction)
    {
        const Vec3 z{0.0, 0.0, 1.0};
        const double d = dot(z, direction);
        if(d > 1.0 - 1e-9)
            return {0.0, 0.0, 0.0, 1.0};
        if(d < -1.0 + 1e-9)
            return {1.0, 0.0, 0.0, 0.0};

        const Vec3 axis = cross(z, direction);
        const double w = 1.0 + d;
        const double len = std::sqrt(dot(axis, axis) + w * w);
        return {axis.x / len, axis.y / len, axis.z / len, w / len};
    }

    static MarkerItemTransform markerItemNodeTransform(const MarkerData &marker, int itemIndex)
    {
        MarkerItemTransform transform;

        if(marker.type == sim_markertype_points)
        {
            transform.position = markerPointAt(marker, itemIndex);
            // Render points as small spheres with radius 0.0025 (diameter 0.005).
            transform.scale = {0.005, 0.005, 0.005};
            return transform;
        }

        if(marker.type == sim_markertype_lines || marker.type == sim_markertype_tubes)
        {
            const int p0Index = itemIndex * 2;
            const int p1Index = p0Index + 1;
            const Vec3 p0 = markerPointAt(marker, p0Index);
            const Vec3 p1 = markerPointAt(marker, p1Index);
            const Vec3 axis = p1 - p0;
            const double len = length(axis);

            transform.position = p0;
            if(len > 1e-12)
                transform.orientation = markerRotationFromZ(axis * (1.0 / len));

            if(marker.type == sim_markertype_tubes)
            {
                // Only the first size value is used, as the tube diameter.
                const double diameter = marker.sizes.empty() ? 0.005 : marker.sizes[0];
                transform.scale = {diameter, diameter, len};
            }
            else
            {
                // Lines are rendered as tubes with radius 0.0025 (diameter 0.005).
                transform.scale = {0.005, 0.005, len};
            }
            return transform;
        }

        return markerItemTransform(marker, itemIndex, markerDefaultItemSize(marker.type));
    }

    int getMarkerItemMesh(_handle_t handle, const MarkerData &marker, const MarkerItemDesc &desc)
    {
        MarkerMeshBuilder mesh;
        MarkerMeshKey key;
        bool transparent = false;
        const int itemIndex = desc.itemIndex;

        switch(marker.type)
        {
        case sim_markertype_points:
        {
            const std::array<double, 4> color = markerColorAt(marker, itemIndex);
            key.color = markerColorKey(color);
            transparent = color[3] < 1.0;
            buildMarkerPointUnit(color, mesh);
            break;
        }
        case sim_markertype_lines:
        {
            const std::array<double, 4> c0 = markerColorAt(marker, itemIndex * 2);
            const std::array<double, 4> c1 = markerColorAt(marker, itemIndex * 2 + 1);
            key.color = (markerColorKey(c0) << 32) | markerColorKey(c1);
            transparent = c0[3] < 1.0 || c1[3] < 1.0;
            buildMarkerLineUnit(c0, c1, mesh);
            break;
        }
        case sim_markertype_tubes:
        {
            const std::array<double, 4> color = markerColorAt(marker, itemIndex * 2);
            key.color = markerColorKey(color);
            transparent = color[3] < 1.0;
            buildMarkerTubeUnit(color, mesh);
            break;
        }
        case sim_markertype_axes:
            buildMarkerAxesUnit(mesh);
            break;
        case sim_markertype_squares:
        case sim_markertype_discs:
        case sim_markertype_spheres:
        case sim_markertype_cubes:
        case sim_markertype_cylinders:
        {
            const std::array<double, 4> color = markerColorAt(marker, itemIndex);
            key.color = markerColorKey(color);
            transparent = color[3] < 1.0;
            buildMarkerOnePointUnit(marker.type, color, mesh);
            break;
        }
        case sim_markertype_custom:
        {
            const std::array<double, 4> color = markerColorAt(marker, itemIndex);
            key.color = markerColorKey(color);
            transparent = color[3] < 1.0;
            buildMarkerCustomUnit(marker, color, mesh);
            break;
        }
        default:
            return -1;
        }

        if(mesh.positions.empty() || mesh.colors.size() != mesh.positions.size() / 3 * 4)
            return -1;

        key.flags = (marker.notShaded ? 1 : 0) | (transparent ? 2 : 0);

        MarkerGeometry geometry;
        if(marker.type == sim_markertype_custom)
        {
            MarkerState &state = markerStates[handle];
            if(state.customGeometry.positionAccessor < 0 || state.customGeometryHash != marker.meshHash)
            {
                state.customGeometry = addMarkerGeometry(mesh, "marker item");
                state.customGeometryHash = marker.meshHash;
                state.itemMeshes.clear();
            }
            geometry = state.customGeometry;
        }
        else
        {
            MarkerGeometry &cached = markerGeometryCache[markerGeometryType(marker.type)];
            if(cached.positionAccessor < 0)
                cached = addMarkerGeometry(mesh, "marker item");
            geometry = cached;
        }

        auto &cache = (marker.type == sim_markertype_custom)
                          ? markerStates[handle].itemMeshes
                          : markerMeshCaches[marker.type];
        auto it = cache.find(key);
        if(it != cache.end())
            return it->second;

        const int colorAccessor = addMarkerColorAccessor(mesh.colors, "marker item");
        const int material = getMarkerItemMaterial(marker.notShaded, transparent);
        const int meshIndex = addMarkerMesh("marker item", geometry, colorAccessor, material);
        cache[key] = meshIndex;
        return meshIndex;
    }

    int createMarkerItemNode(_handle_t handle, const MarkerData &marker, const MarkerItemDesc &desc)
    {
        MarkerState &state = markerStates[handle];
        if(state.nodeIndex < 0 || desc.itemIndex < 0)
            return -1;

        const int meshIndex = getMarkerItemMesh(handle, marker, desc);
        if(meshIndex < 0)
            return -1;

        const int nodeIndex = model.nodes.size();
        model.nodes.push_back({});
        model.nodes[state.nodeIndex].children.push_back(nodeIndex);
        model.nodes[nodeIndex].name = "marker item " + std::to_string(desc.id);

        const MarkerItemTransform transform = markerItemNodeTransform(marker, desc.itemIndex);
        const auto nearZero = [](double v) { return std::fabs(v) < 1e-12; };
        const auto nearOne = [](double v) { return std::fabs(v - 1.0) < 1e-12; };

        // Leave default TRS components out of the JSON: marker-heavy files
        // contain many item nodes and every omitted property matters.
        if(!nearZero(transform.position.x) || !nearZero(transform.position.y) ||
           !nearZero(transform.position.z))
        {
            model.nodes[nodeIndex].translation = {transform.position.x, transform.position.y,
                                                  transform.position.z};
        }
        if(!nearZero(transform.orientation[0]) || !nearZero(transform.orientation[1]) ||
           !nearZero(transform.orientation[2]) || !nearOne(transform.orientation[3]))
        {
            model.nodes[nodeIndex].rotation = {transform.orientation[0], transform.orientation[1],
                                              transform.orientation[2], transform.orientation[3]};
        }
        if(!nearOne(transform.scale.x) || !nearOne(transform.scale.y) ||
           !nearOne(transform.scale.z))
        {
            model.nodes[nodeIndex].scale = {transform.scale.x, transform.scale.y, transform.scale.z};
        }
        model.nodes[nodeIndex].mesh = meshIndex;

        MarkerItemState &item = state.items[desc.id];
        item.nodeIndex = nodeIndex;
        return nodeIndex;
    }

    static void markMarkerItemsInvisible(MarkerState &state, size_t timeIndex)
    {
        for(auto &pair : state.items)
            pair.second.visible[timeIndex] = false;
    }

    void recordMarkerItems(_handle_t handle, size_t timeIndex)
    {
        auto stateIt = markerStates.find(handle);
        if(stateIt == markerStates.end())
            return;
        MarkerState &state = stateIt->second;

        const MarkerData marker = readMarkerData(handle);
        const std::vector<MarkerItemDesc> items = buildMarkerItemDescs(marker);
        std::set<int64_t> seen;
        for(const auto &desc : items)
        {
            auto it = state.items.find(desc.id);
            if(it == state.items.end())
            {
                if(createMarkerItemNode(handle, marker, desc) < 0)
                    continue;
                it = state.items.find(desc.id);
                if(it == state.items.end())
                    continue;
            }

            const int meshIndex = getMarkerItemMesh(handle, marker, desc);
            if(meshIndex >= 0 && it->second.nodeIndex >= 0)
                model.nodes[it->second.nodeIndex].mesh = meshIndex;

            it->second.visible[timeIndex] = true;
            it->second.transforms[timeIndex] = markerItemNodeTransform(marker, desc.itemIndex);
            seen.insert(desc.id);
        }

        for(auto &pair : state.items)
        {
            if(seen.find(pair.first) == seen.end())
                pair.second.visible[timeIndex] = false;
        }
    }

    void exportMarkerContent(_handle_t obj, int nodeIndex)
    {
        const MarkerData marker = readMarkerData(obj);
        MarkerState &state = markerStates[obj];
        state.nodeIndex = nodeIndex;
        model.nodes[nodeIndex].extras = markerExtras(marker);

        if(marker.points.empty())
            return;

        if(marker.type == sim_markertype_triangles)
        {
            // TODO: triangle markers are not handled yet.
            sim::addLog(sim_verbosity_debug,
                        "exportMarkerContent: %s: triangle markers are not exported yet",
                        getObjectName(obj).c_str());
            return;
        }
        if(marker.type < 0)
        {
            sim::addLog(sim_verbosity_debug,
                        "exportMarkerContent: %s: unsupported marker type '%s'",
                        getObjectName(obj).c_str(), marker.typeString.c_str());
            return;
        }

        const std::vector<MarkerItemDesc> items = buildMarkerItemDescs(marker);
        for(const auto &desc : items)
        {
            if(state.items.find(desc.id) == state.items.end())
                createMarkerItemNode(obj, marker, desc);
        }
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
            if((isShape(obj) && isVisible(obj) && !isWireframe(obj)) ||
               (isMarker(obj) && isVisible(obj)) ||
               isCamera(obj))
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
            if(objectExists(handle) && isJoint(handle))
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

        // Marker items are handled like individual shapes: each item node gets
        // translation/rotation/scale tracks when it moves, and its scale is
        // zero while it is absent. Items are created lazily, so their earlier
        // frames also read as invisible.
        const auto addItemChannel = [&](const std::string &name,
                                        const std::vector<gltfFloat> &values,
                                        int components, int targetNode,
                                        const std::string &path)
        {
            int c = model.animations[0].channels.size();
            model.animations[0].channels.push_back({});
            model.animations[0].samplers.push_back({});
            int b = addBuffer(values.data(), sizeof(gltfFloat) * values.size(), name);
            int v = addBufferView(b, sizeof(gltfFloat) * values.size(), 0, name);
            std::vector<double> mn, mx;
            minMaxVec(values, components, mn, mx);
            int a = addAccessor(v, 0, TINYGLTF_COMPONENT_TYPE_gltfFloat,
                                components == 4 ? TINYGLTF_TYPE_VEC4 : TINYGLTF_TYPE_VEC3,
                                n, mn, mx, name);
            model.animations[0].samplers[c].interpolation = "STEP";
            model.animations[0].samplers[c].input  = at;
            model.animations[0].samplers[c].output = a;
            model.animations[0].channels[c].sampler     = c;
            model.animations[0].channels[c].target_node = targetNode;
            model.animations[0].channels[c].target_path = path;
        };

        const auto sameVec = [](const Vec3 &a, const Vec3 &b)
        {
            return std::fabs(a.x - b.x) < 1e-9 && std::fabs(a.y - b.y) < 1e-9 &&
                   std::fabs(a.z - b.z) < 1e-9;
        };
        const auto sameQuat = [](const std::array<double, 4> &a, const std::array<double, 4> &b)
        {
            return std::fabs(a[0] - b[0]) < 1e-9 && std::fabs(a[1] - b[1]) < 1e-9 &&
                   std::fabs(a[2] - b[2]) < 1e-9 && std::fabs(a[3] - b[3]) < 1e-9;
        };

        for(auto &markerPair : markerStates)
        {
            MarkerState &state = markerPair.second;
            if(state.nodeIndex < 0)
                continue;
            const std::string markerName = model.nodes[state.nodeIndex].name;

            for(auto &itemPair : state.items)
            {
                MarkerItemState &item = itemPair.second;
                if(item.nodeIndex < 0 || item.visible.empty())
                    continue;

                MarkerItemTransform fallback;
                bool haveFallback = false;
                for(int i = 0; i < n && !haveFallback; i++)
                {
                    auto vit = item.visible.find(static_cast<size_t>(i));
                    if(vit == item.visible.end() || !vit->second)
                        continue;
                    auto tit = item.transforms.find(static_cast<size_t>(i));
                    if(tit != item.transforms.end())
                    {
                        fallback = tit->second;
                        haveFallback = true;
                    }
                }
                if(!haveFallback)
                    continue;

                std::vector<gltfFloat> p(n * 3), r(n * 4), s(n * 3);
                bool anyInvisible = false;
                bool translationChanged = false;
                bool rotationChanged = false;
                bool scaleChanged = false;
                MarkerItemTransform current = fallback;

                for(int i = 0; i < n; i++)
                {
                    auto vit = item.visible.find(static_cast<size_t>(i));
                    const bool visible = (vit != item.visible.end() && vit->second);
                    auto tit = item.transforms.find(static_cast<size_t>(i));
                    if(visible && tit != item.transforms.end())
                        current = tit->second;
                    if(!visible)
                        anyInvisible = true;

                    p[3 * i + 0] = static_cast<gltfFloat>(current.position.x);
                    p[3 * i + 1] = static_cast<gltfFloat>(current.position.y);
                    p[3 * i + 2] = static_cast<gltfFloat>(current.position.z);
                    for(int j = 0; j < 4; j++)
                        r[4 * i + j] = static_cast<gltfFloat>(current.orientation[j]);
                    if(visible)
                    {
                        s[3 * i + 0] = static_cast<gltfFloat>(current.scale.x);
                        s[3 * i + 1] = static_cast<gltfFloat>(current.scale.y);
                        s[3 * i + 2] = static_cast<gltfFloat>(current.scale.z);
                        if(!sameVec(current.position, fallback.position))
                            translationChanged = true;
                        if(!sameQuat(current.orientation, fallback.orientation))
                            rotationChanged = true;
                        if(!sameVec(current.scale, fallback.scale))
                            scaleChanged = true;
                    }
                    else
                    {
                        s[3 * i + 0] = static_cast<gltfFloat>(0.0);
                        s[3 * i + 1] = static_cast<gltfFloat>(0.0);
                        s[3 * i + 2] = static_cast<gltfFloat>(0.0);
                    }
                }

                if(!anyInvisible && !translationChanged && !rotationChanged && !scaleChanged)
                    continue;

                if(translationChanged)
                    addItemChannel(markerName + " marker item translation", p, 3,
                                   item.nodeIndex, "translation");
                if(rotationChanged)
                    addItemChannel(markerName + " marker item rotation", r, 4,
                                   item.nodeIndex, "rotation");
                if(scaleChanged || anyInvisible)
                    addItemChannel(markerName + " marker item scale", s, 3,
                                   item.nodeIndex, "scale");
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
                // The object may already have been destroyed again before its
                // first animation frame was recorded.
                frames[handle].parentHandle = objectExists(handle) ? getParentHandle(handle) : -1;
            }

            try
            {
                simPose3D_get(&frames[handle].track[timeIndex], handle, frames[handle].parentHandle);

                if(markerStates.find(handle) != markerStates.end())
                    recordMarkerItems(handle, timeIndex);

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
                auto markerIt = markerStates.find(handle);
                if(markerIt != markerStates.end())
                    markMarkerItemsInvisible(markerIt->second, timeIndex);
            }
        }
    }

    void setExportTextureFormat(setExportTextureFormat_in *in, setExportTextureFormat_out *out)
    {
        TextureFormat format = static_cast<TextureFormat>(in->textureFormat);
        if(format != simgltf_texture_format_png && format != simgltf_texture_format_jpg)
        {
            // glTF 2.0 only supports PNG and JPEG images.
            sim::addLog(sim_verbosity_warnings,
                        "setExportTextureFormat: format %d is not valid in glTF, using PNG instead",
                        in->textureFormat);
            format = simgltf_texture_format_png;
        }
        exportTextureFormat = format;
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
    std::map<_handle_t, MarkerState> markerStates; // marker items, keyed by marker handle
    std::map<int, std::map<MarkerMeshKey, int>> markerMeshCaches; // marker type -> item mesh
    std::map<int, MarkerGeometry> markerGeometryCache; // marker type -> shared unit geometry
    std::map<int, int> markerMaterials; // (notShaded | transparent << 1) -> material
    std::vector<BufferSlice> bufferSlices;
    int sharedBufferIndex = -1;

    _handle_t getParentHandle(_handle_t h)
    {
        return sim::getHandleProperty(h, "parent");
    }

    // Objects can be created and destroyed while the animation is being
    // recorded. Use the throwing wrapper, but clear the error it leaves behind
    // when the handle no longer resolves (e.g. a marker that disappeared
    // mid-simulation).
    static bool objectExists(_handle_t h)
    {
        try
        {
            sim::getObjectType(h);
            return true;
        }
        catch(...)
        {
            sim::clearLastError();
            return false;
        }
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
        {
            // The node may have been created while the marker was invisible
            // (e.g. as an ancestor of a visible object); add its items now.
            if(isMarker(handle) && isVisible(handle) &&
               markerStates.find(handle) == markerStates.end())
                exportMarkerContent(handle, it->second);
            return it->second;
        }

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
        else if(isMarker(handle) && isVisible(handle))
            exportMarkerContent(handle, nodeIndex);
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

    TextureFormat exportTextureFormat = simgltf_texture_format_png;
};

SIM_PLUGIN(Plugin)
#include "stubsPlusPlus.cpp"
