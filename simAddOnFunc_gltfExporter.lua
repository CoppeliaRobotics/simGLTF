local scenePath=sim.getStringParameter(sim.stringparam_scene_path)
local sceneName=sim.getStringParameter(sim.stringparam_scene_name):match("(.+)%..+")
if sceneName==nil then sceneName='untitled' end
local fileName=sim.fileDialog(sim.filedlg_type_save,'Export to glTF...',scenePath,sceneName..'.gltf','glTF file','gltf')
if fileName==nil then return end
local includeAnimation=simGLTF.animationFrameCount()>0 and sim.msgBox(sim.msgbox_type_question,sim.msgbox_buttons_yesno,'Export to glTF...','Include also animation content (from last simulation run)?')==sim.msgbox_return_yes
m=simGLTF.create()
simGLTF.exportAllObjects(m)
if includeAnimation then
    simGLTF.exportAnimation(m)
end
simGLTF.saveASCII(m,fileName)
sim.addStatusbarMessage('Exported glTF content to '..fileName)
simGLTF.destroy(m)
