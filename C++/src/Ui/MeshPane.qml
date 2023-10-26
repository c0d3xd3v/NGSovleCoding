import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3

Pane {
    id: _meshPane
    padding: 0
    width: parent.parent.width
    height: 200
    property var hashString: null

    signal deleteMesh(var hashString, var pane)

    Frame {
        padding: 0
        anchors.fill: parent
        Frame{
            id: toolBar
            anchors.top: parent.top
            anchors.left: parent.left
            anchors.right: parent.right
            padding: 0
            RowLayout {
                spacing: 0
                anchors.right: parent.right
                anchors.left: parent.left
                ToolButton {
                    flat: false
                    height:64
                    width:64
                    icon.source: "qrc:/icons/close.svg"
                    onClicked: {
                        deleteMesh(hashString, _meshPane)
                    }
                }
                Item {
                    Layout.fillWidth: true
                }
                ToolButton {
                    flat: true
                    height:64
                    width:64
                    rotation: -90
                    icon.source: "qrc:/icons/arrow_back.svg"
                    onClicked: {
                        rotation = rotation * -1;
                    }
                }
            }
        }
        GridLayout {
            anchors.top: toolBar.bottom
            ToolButton {
                height:64
                width:64
                flat: true
                rotation: -90
                icon.source: "qrc:/icons/arrow_back.svg"
                onClicked: {
                    rotation = rotation * -1;
                }
            }
            ToolButton {
                height:64
                width:64
                flat: true
                rotation: -90
                icon.source: "qrc:/icons/arrow_back.svg"
                onClicked: {
                    rotation = rotation * -1;
                }
            }
        }
        Item {
            Layout.fillWidth: true
            Layout.fillHeight: true
        }
    }
}
