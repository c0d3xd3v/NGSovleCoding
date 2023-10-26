// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3
import QtQuick.Window 2.15
import QtQuick.Dialogs 1.3

import VTK 9.2

import "qrc:"

// window containing the application
ApplicationWindow {
    // title of the application
    title: qsTr("VTK QtQuick App")
    width: 800
    height: 550
    visible: true
    id: window
    readonly property bool inPortrait: window.width < 750

    FileDialog {
        id: fileDialog
        title: "Please choose a file"
        folder: shortcuts.home
        nameFilters: ["surface mesh files (*.stl *.obj)"]
        onAccepted: {
            var hashString = vtkqtcontroller.loadFile(fileDialog.fileUrls);
            var comp = Qt.createComponent("MeshPane.qml");
            var obj = comp.createObject(meshListPane, {hashString: hashString});
            obj.deleteMesh.connect(deleteMesh)
        }
    }

    function deleteMesh(hashString, pane)
    {
        //console.log("delete mesh " + hashString)
        vtkqtcontroller.removeObject(hashString)
        pane.destroy()
    }

    onBeforeRendering: {
        vtkqtcontroller.resize(renderItem.width, renderItem.height);
    }

    Drawer {
        id: drawer
        width: 250
        height: parent.height
        padding: 0
        modal: inPortrait
        interactive: inPortrait
        position: inPortrait ? 0 : 1
        visible: !inPortrait
        clip: true

        Pane {
            id: sidePane
            anchors.fill: parent
            leftPadding: 0
            padding: 0
            clip: true
            spacing: 0

            MainToolbar{id:toolbar}

            Pane {
                padding: 3
                id: listPane
                anchors.top: toolbar.bottom
                anchors.right: parent.right
                anchors.left: parent.left
                anchors.bottom: parent.bottom
                width: drawer.width
                ColumnLayout {
                    id: meshListPane
                }
            }
        }
    }

    // Instantiate the vtk render window
    VTKRenderWindow {
      id: vtkwindow
      width: window.width
      height: window.height
      //anchors.fill: window
      anchors.leftMargin: !inPortrait ? drawer.width : undefined
    }

    // add one or more vtk render items
    VTKRenderItem {
        id: renderItem
      objectName: "ConeView"
      width: vtkwindow.width
      height: vtkwindow.height
      anchors.fill: vtkwindow
      anchors.leftMargin: !inPortrait ? drawer.width : undefined
      // Provide the handle to the render window
      renderWindow: vtkwindow
    }
}
