local appPath=sim.getStringParameter(sim.stringparam_application_path)
local fileName=sim.fileDialog(sim.filedlg_type_save,'Export to glTF...',appPath,'scene.gltf','glTF file','gltf')
local includeAnimation=sim.msgBox(sim.msgbox_type_question,sim.msgbox_buttons_yesno,'Export to glTF...','Include also animation content (from last simulation run)?')==sim.msgbox_return_yes
if fileName~='' then
    m=simGLTF.create()
    simGLTF.exportAllObjects(m)
    if includeAnimation then
        simGLTF.exportAnimation(m)
    end
    simGLTF.saveASCII(m,fileName)
    simGLTF.destroy(m)
end
