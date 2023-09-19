// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Window 2.15

// import the VTK module
import VTK 9.2

// window containing the application
ApplicationWindow {
    // title of the application
    title: qsTr("VTK QtQuick App")
    width: 400
    height: 400
    visible: true

    Drawer {
        id: drawer
        width: 250
        height: parent.height

        Label {
            text: "Content goes here!"
            anchors.centerIn: parent
        }
    }

    // Instantiate the vtk render window
    VTKRenderWindow {
      id: vtkwindow
      width: parent.width
      height: parent.height
    }

    // add one or more vtk render items
    VTKRenderItem {
      objectName: "ConeView"
      width: parent.width
      height: parent.height
      // Provide the handle to the render window
      renderWindow: vtkwindow
    }

}
