import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Item {
    id: pane
    width: 240
    height: 400

    ToolBar {
        id: toolBar
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top

        RowLayout {
            id: rowLayout
            anchors.fill: parent

            BackButton {
                id: toolButton
            }

            Item {
                id: item1
                Layout.fillHeight: true
                Layout.fillWidth: true
            }
        }
    }

    ScrollView {
        id: scrollView
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: toolBar.bottom
        anchors.bottom: parent.bottom
        anchors.topMargin: 0
    }
}
