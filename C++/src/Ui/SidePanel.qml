import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

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

    Frame {
        id: sidePane
        anchors.fill: parent
        padding: 0
        clip: true

        MainToolbar{id:toolbar}
        ColumnLayout {
            anchors.top: toolbar.bottom
            Layout.preferredHeight: drawer.height
            ListView {
                leftMargin: 3
                rightMargin: 3
                topMargin: 3
                bottomMargin: 3
                Layout.fillHeight: true
                Layout.preferredHeight: sidePane.height - toolbar.height
                Layout.preferredWidth: sidePane.width
                highlight: Rectangle {
                    opacity: 1.0;
                    color: "transparent";
                    border.color: "lightsteelblue";
                    border.width: 1.5;
                    radius: 0
                }
                id: loadedMeshesList
                clip: true
                boundsBehavior: Flickable.StopAtBounds
                spacing: 3
                delegate:  MeshPane2 {
                    Layout.preferredHeight: 210
                    width: if(parent) parent.width; else 1
                    onHasFocus: function (index) {
                        loadedMeshesList.currentIndex = index
                    }
                }
                model: loadedMeshesListModel
            }
        }
    }

    function deleteMesh(hashString, pane)
    {
        loadedMeshesListModel.remove(pane)
        //console.log(pane.mesh_controller)
        //vtkqtcontroller.removeObject(hashString, pane.mesh_controller)
    }
}
