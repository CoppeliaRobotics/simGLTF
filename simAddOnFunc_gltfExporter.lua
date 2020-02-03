local scenePath=sim.getStringParameter(sim.stringparam_scene_path)
local sceneName=sim.getStringParameter(sim.stringparam_scene_name):match("(.+)%..+")
if sceneName==nil then sceneName='untitled' end
local fileName=sim.fileDialog(sim.filedlg_type_save,'Export to glTF...',scenePath,sceneName..'.gltf','glTF file','gltf')
if fileName==nil then return end
local exportAnimation=simGLTF.animationFrameCount()>0 and sim.msgBox(sim.msgbox_type_question,sim.msgbox_buttons_yesno,'Export to glTF...','Export animation content from last simulation run?')==sim.msgbox_return_yes
if exportAnimation then
    -- cleared automatically on simulation start
    simGLTF.exportAnimation()
else
    simGLTF.clear()
    simGLTF.exportAllObjects()
end
simGLTF.saveASCII(fileName)
sim.addStatusbarMessage('Exported glTF content to '..fileName)
