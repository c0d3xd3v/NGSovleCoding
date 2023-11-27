// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3

Drawer {
    id: drawer
    width: 250
    height: parent.height
    modal: inPortrait
    interactive: inPortrait
    position: inPortrait ? 0 : 1
    visible: !inPortrait
    clip: true

    property var loadedMeshesListModel: ListModel {}

    Pane {
        id: sidePane
        anchors.fill: parent
        padding: 0
        clip: true

        MainToolbar{id:toolbar}

        ListView {
            leftMargin: 3
            rightMargin: 3
            topMargin: 3
            bottomMargin: 3
            anchors.top: toolbar.bottom
            anchors.right: parent.right
            anchors.left: parent.left
            anchors.bottom: parent.bottom
            id: loadedMeshesList
            clip: true
            spacing: 3
            delegate: MeshPane {}
            model: loadedMeshesListModel
        }
    }

    function deleteMesh(hashString, pane)
    {
        console.log("delete mesh " + hashString)
        vtkqtcontroller.removeObject(hashString)
        loadedMeshesListModel.remove(pane)
    }
}
