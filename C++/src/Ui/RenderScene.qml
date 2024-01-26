import QtQuick

import VTK 9.3

// Instantiate the vtk render window
VTKRenderWindow {
  id: vtkwindow
  width: window.width
  height: window.height
  //anchors.fill: window
  anchors.leftMargin: !inPortrait ? sidePanel.width : undefined
  // add one or more vtk render items
  property var renderSceneItem: renderItem
  VTKRenderItem {
    id: renderItem
    objectName: "ConeView"
    width: vtkwindow.width
    height: vtkwindow.height
    anchors.fill: vtkwindow
    anchors.leftMargin: !inPortrait ? sidePanel.width : undefined
    // Provide the handle to the render window
    renderWindow: vtkwindow
  }
}
