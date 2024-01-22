import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3

ToolBar {
    id:toolbar
    anchors.left: parent.left
    anchors.top: parent.Top
    width: drawer.width
    RowLayout {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        spacing: 0
        ToolButton {
            flat: false
            icon.height: 32
            icon.width: 32
            icon.color: "transparent"
            icon.source: "qrc:icons/openfile.png"
            onClicked: {
                fileDialog.visible = true
                //fileDialog.sidebarVisible = false
            }
        }
        Item {
            Layout.fillWidth: true
        }
    }
}

