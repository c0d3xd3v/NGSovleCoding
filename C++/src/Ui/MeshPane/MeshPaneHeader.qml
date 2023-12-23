import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15

RowLayout {
    Layout.fillWidth: true
    Layout.preferredHeight: clsBtn.height
    Layout.alignment: Qt.AlignLeft
    spacing: 1
    ToolButton {
        padding: 0
        id: clsBtn
        flat: false
        icon.source: "qrc:/icons/close.svg"
        onClicked: {
            //console.log(mesh_controller)
            mesh_controller.cleanupRendering()
            deleteMesh(hashString, index)
        }
        Component.onCompleted: {parent.height=height}
    }
    MouseArea{
        Layout.preferredHeight: clsBtn.height
        Layout.fillWidth: true
        onClicked: {
            if(!root_item.focus){
                root_item.focus = true
                console.log("clicked")
            }else if(root_item.focus && root_item.height === ref_height)
            {
                root_item.height = clsBtn.height + 10.*parent.spacing
                console.log(ref_height)
            } else {
                root_item.height = ref_height
            }
        }
        Label {
            anchors.fill: parent
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter
            text: hashString
        }
    }
}
