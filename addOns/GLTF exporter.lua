sim = require 'sim'

function sysCall_info()
    return {autoStart = false, menu = 'Exporters\nGLTF exporter...'}
end

function sysCall_init()
    simGLTF = require 'simGLTF'
    local scenePath = sim.getStringParam(sim.stringparam_scene_path)
    local sceneName = sim.getStringParam(sim.stringparam_scene_name):match("(.+)%..+")
    if sceneName == nil then sceneName = 'untitled' end
    local fileName = simUI.fileDialog(
                         simUI.filedialog_type.save, 'Export to glTF...', scenePath,
                         sceneName .. '.gltf', 'glTF file', 'gltf'
                     )
    if fileName == nil then return end
    local fmt, fmtName = simGLTF.getExportTextureFormat()
    sim.addLog(
        sim.verbosity_scriptinfos, 'Texture export format is set to "' .. fmtName ..
            '". You can change that with simGLTF.setExportTextureFormat(format).'
    )
    simGLTF.clear()
    simGLTF.exportAllObjects()
    simGLTF.saveASCII(fileName)
    sim.addLog(sim.verbosity_scriptinfos, 'Exported glTF content to ' .. fileName)
    return {cmd = 'cleanup'}
end
