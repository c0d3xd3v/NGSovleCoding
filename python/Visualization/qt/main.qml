import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import QmlVtk 1.0

Window {
    width: 600
    height: 600
    visible: true

    RowLayout {
        anchors.fill: parent
        VTKItem {
            objectName: "ConeView"
            Layout.fillHeight: true
            Layout.fillWidth: true
        }
        /*
        VTKItem {
            objectName: "ConeView2"
            Layout.fillHeight: true
            Layout.fillWidth: true
        }
        */
    }
}
