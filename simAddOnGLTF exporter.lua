sim = require 'sim'

function sysCall_info()
    return {autoStart = false, menu = 'Exporters\nGLTF exporter...'}
end

function sysCall_init()
    simGLTF = require 'simGLTF'
    local scenePath = sim.getStringParameter(sim.stringparam_scene_path)
    local sceneName = sim.getStringParameter(sim.stringparam_scene_name):match("(.+)%..+")
    if sceneName == nil then sceneName = 'untitled' end
    local fileName = sim.fileDialog(
                         sim.filedlg_type_save, 'Export to glTF...', scenePath,
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
