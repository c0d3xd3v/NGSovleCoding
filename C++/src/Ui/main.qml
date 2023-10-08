// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3
import QtQuick.Window 2.15
import QtQuick.Dialogs 1.3
import Qt.labs.settings 1.0

// import the VTK module
import VTK 9.2

// window containing the application
ApplicationWindow {
    // title of the application
    title: qsTr("VTK QtQuick App")
    width: 600
    height: 400
    visible: true
    id: window
    readonly property bool inPortrait: window.width < 750


    FileDialog {
        id: fileDialog
        title: "Please choose a file"
        folder: shortcuts.home
        nameFilters: ["surface mesh files (*.stl *.obj)"]
        onAccepted: {
            console.log("You chose: " + fileDialog.fileUrls)
            //Qt.quit()
        }
        onRejected: {
            console.log("Canceled")
            //Qt.quit()
        }
    }

    onBeforeRendering: {
        vtkqtcontroller.resize(renderItem.width, renderItem.height);
    }

    Drawer {
        id: drawer
        width: 250
        height: parent.height

        modal: inPortrait
        interactive: inPortrait
        position: inPortrait ? 0 : 1
        visible: !inPortrait
        clip: true

        Pane {
            anchors.fill: parent
            leftPadding: 0
            padding: 0
            clip: true
            spacing: 0
            ToolBar {
                id:toolbar
                anchors.left: parent.left
                anchors.top: parent.Top
                width: drawer.width
                RowLayout {
                    anchors.left: parent.left
                    anchors.top: parent.top
                    anchors.right: parent.right
                    spacing: 0
                    ToolButton {
                        flat: false
                        height:64
                        width:64
                        icon.source: "qrc:icons/openfile.svg"
                        onClicked: {
                            fileDialog.visible = true
                            fileDialog.sidebarVisible = false
                        }
                    }
                    Item {
                        Layout.fillWidth: true
                    }
                    ToolButton {
                        height:64
                        width:64
                        flat: false
                        text: qsTr("⋮")
                        //onClicked: menu.open()
                    }
                }
            }
            Pane {
                anchors.top: toolbar.bottom
                anchors.left: parent.left
                anchors.bottom: parent.bottom
                anchors.right: parent.right
                ColumnLayout {
                    anchors.fill: parent
                    Label {
                        text: "Title"
                        elide: Label.ElideRight
                        horizontalAlignment: Qt.AlignHCenter
                        verticalAlignment: Qt.AlignVCenter
                        //width: drawer.width
                        Layout.fillWidth: true
                    }
                    Rectangle {
                        id: fileDlg
                        //width: drawer.width
                        Layout.fillWidth: true
                        Layout.fillHeight: true
                    }
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
