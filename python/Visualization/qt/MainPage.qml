import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs

Item {
    id: item2
    width: 240
    height: 400

    FileDialog {
           id: fileDialog
           title: "Please choose a file"
   //        currentFolder: shortcuts.home
           nameFilters: ["surface mesh files (*.stl *.obj)"]
           onAccepted: function() {
               MainCtrl.loadSurfaceMesh(currentFile)
           }
    }

    ToolBar {
        id: toolBar
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top

        RowLayout {
            id: rowLayout
            anchors.fill: parent

            ToolButton {
                id: toolButton
                text: qsTr("open file")
                display: AbstractButton.TextBesideIcon
                onClicked: fileDialog.visible = true
            }

            RenderSettingsButton {
            }

            ExportSettingsButton {
            }

            Item {
                id: item1
                Layout.fillHeight: true
                Layout.fillWidth: true
            }


        }
    }

    ScrollView {
        id: mainPageScrollView
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: toolBar.bottom
        anchors.bottom: parent.bottom
        anchors.rightMargin: 5
        //anchors.rightMargin: 5
        anchors.leftMargin: 5
        anchors.bottomMargin: 0
        anchors.topMargin: 0


        ColumnLayout {
            id: columnLayout
            x: 0
            y: 0
            width: mainPageScrollView.width
            height: 388
            spacing: 6

            RowLayout {
                id: rowLayout1
                height: 100
                Layout.fillWidth: true
                spacing: 5

                Label {
                    id: label
                    text: qsTr("Label")
                }

                TextField {
                    id: textField
                    Layout.margins: 5
                    Layout.topMargin: 15
                    Layout.fillWidth: true
                    placeholderText: qsTr("Text Field")
                }

            }



            RowLayout {
                id: rowLayout2
                height: 100
                Layout.fillWidth: true

                Label {
                    id: label1
                    text: qsTr("Label")
                }

                TextField {
                    id: textField1
                    Layout.margins: 5
                    Layout.fillWidth: true
                    placeholderText: qsTr("Text Field")
                }

            }


            RowLayout {
                id: rowLayout3
                width: 100
                height: 100
                Layout.fillWidth: true

                Label {
                    id: label2
                    text: qsTr("Label")
                }

                TextField {
                    id: textField2
                    Layout.margins: 5
                    Layout.fillWidth: true
                    placeholderText: qsTr("Text Field")
                }

            }






            Button {
                id: button
                text: qsTr("Button")
                Layout.fillWidth: true
            }




            Item {
                id: item3
                width: 200
                height: 200
                Layout.fillHeight: true
                Layout.fillWidth: true
            }






        }
    }

}
