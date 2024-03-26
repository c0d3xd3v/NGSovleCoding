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
            height: 290
            spacing: 6

            ColumnLayout {
                id: rowLayout1
                Layout.fillHeight: false
                Layout.fillWidth: true
                spacing: 0

                Label {
                    id: label
                    text: qsTr("Label")
                    Layout.margins: 5
                }

                TextField {
                    id: textField
                    Layout.bottomMargin: 0
                    leftPadding: 8
                    padding: 8
                    Layout.margins: 5
                    Layout.topMargin: 0
                    Layout.fillWidth: true
                    placeholderText: qsTr("Text Field")
                }

            }



            ColumnLayout {
                id: rowLayout2
                Layout.fillWidth: true

                Label {
                    id: label1
                    text: qsTr("Label")
                    Layout.bottomMargin: 0
                    Layout.margins: 5
                }

                TextField {
                    id: textField1
                    Layout.bottomMargin: 0
                    Layout.topMargin: 0
                    Layout.margins: 5
                    Layout.fillWidth: true
                    placeholderText: qsTr("Text Field")
                }

            }


            ColumnLayout {
                id: rowLayout3
                Layout.fillWidth: true
                spacing: 5

                Label {
                    id: label2
                    text: qsTr("Label")
                    Layout.bottomMargin: 0
                    Layout.margins: 5
                }

                TextField {
                    id: textField2
                    Layout.topMargin: 0
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
                Layout.fillHeight: true
                Layout.fillWidth: true
            }






        }
    }

}
