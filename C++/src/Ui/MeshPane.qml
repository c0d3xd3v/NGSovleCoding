import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3

Frame {
    height: 300
    anchors.right: parent.right
    anchors.left: parent.left
    padding: 0

    property var _meshPane : this

    Frame{
        id: toolBar
        anchors.right: parent.right
        anchors.left: parent.left
        height: clsBtn.height
        padding: 0
        MouseArea {
            anchors.fill: parent
            onClicked: {
                console.log(_meshPane.height)
                if(_meshPane.height == 300)
                    _meshPane.height = clsBtn.height;
                else _meshPane.height = 300;
            }
        }
        RowLayout {
            spacing: 0
            Layout.fillWidth: true
            ToolButton {
                id: clsBtn
                flat: false
                height:32
                width:32
                icon.source: "qrc:/icons/close.svg"
                onClicked: {
                    deleteMesh(hashString, index)
                }
            }
            Text {
                id:txt
                width: loadedMeshesList.width
                text: hashString
                horizontalAlignment: Qt.AlignHCenter
            }
        }
    }
}
