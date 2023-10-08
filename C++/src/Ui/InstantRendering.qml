// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.3
import QtQuick.Window 2.15

// import the VTK module
import VTK 9.2

// window containing the application
ApplicationWindow {
    // title of the application
    title: qsTr("VTK QtQuick App")
    width: 600
    height: 400
    visible: true
    id: window
    readonly property bool inPortrait: window.width < 750

    onBeforeRendering: {
        vtkqtcontroller.resize(renderItem.width, renderItem.height);
    }

    // Instantiate the vtk render window
    VTKRenderWindow {
      id: vtkwindow
      width: window.width
      height: window.height
      //anchors.fill: window
      anchors.leftMargin: !inPortrait ? drawer.width : undefined
    }

    // add one or more vtk render items
    VTKRenderItem {
        id: renderItem
      objectName: "ConeView"
      width: vtkwindow.width
      height: vtkwindow.height
      anchors.fill: vtkwindow
      anchors.leftMargin: !inPortrait ? drawer.width : undefined
      // Provide the handle to the render window
      renderWindow: vtkwindow
    }
}
