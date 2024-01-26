import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "qrc:/MeshPane/"

Frame {
    //background:Rectangle {opacity: 1.0}
    padding: 5
    signal hasFocus(index: int)
    onFocusChanged: if(focus) hasFocus(index)
    property var root_item: this
    property int ref_height: 0
    Component.onCompleted: ref_height = height
    ColumnLayout {
        Layout.margins: 3
        anchors.fill: parent
        spacing: 5
        clip: true
        MeshPaneHeader{}
        ToolSeparator {
            Layout.fillWidth: true
            orientation: Qt.Horizontal
        }
        //MeshPaneHorizontalSeparator{}
        MeshPanePageSwipe{id:view}
        ToolSeparator {
            Layout.fillWidth: true
            orientation: Qt.Horizontal
        }
        //MeshPaneHorizontalSeparator{}
        MeshPanePageSelect{}
    }
}
