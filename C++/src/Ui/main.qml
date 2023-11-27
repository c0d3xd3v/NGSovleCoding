// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3
import QtQuick.Window 2.15
import QtQuick.Dialogs 1.3

import "qrc:"

ApplicationWindow
{
    id: window
    title: qsTr("VTK QtQuick App")
    width: 800
    height: 550
    visible: true

    readonly property bool inPortrait: window.width < 750

    function loadMesh(path)
    {
        var hashString = vtkqtcontroller.loadFile(path);
        var comp = Qt.createComponent("MeshPane.qml");
        sidePanel.loadedMeshesListModel.append({"hashString" : hashString})
    }

    FileDialog {
        id: fileDialog
        title: "Please choose a file"
        folder: shortcuts.home
        nameFilters: ["surface mesh files (*.stl *.obj)"]
        onAccepted: loadMesh(fileDialog.fileUrls)
    }
    onBeforeRendering: {
        var w = renderScene.renderSceneItem.width
        var h = renderScene.renderSceneItem.height
        vtkqtcontroller.resize(w, h);
    }

    SidePanel{id:sidePanel}
    RenderScene{id: renderScene}
}
