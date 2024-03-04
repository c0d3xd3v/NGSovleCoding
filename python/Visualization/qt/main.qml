import QtQuick.Window
import VTK 9.3

Window {
    width: 600
    height: 600
    VTKRenderWindow {
        id: vtkwindow
        anchors.fill: parent
        VTKRenderItem {
         id: renderItem
         objectName: "ConeView"
         anchors.fill: vtkwindow
         renderWindow: vtkwindow
        }
    }
}
