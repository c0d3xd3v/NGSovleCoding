import sys
from PySide6 import QtCore, QtWidgets
from PySide6.QtGui import QIcon, QAction
from PySide6.QtWidgets import QToolBar

from vtk import vtkCommand
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from Visualization.vtkhelper import *

class VtkFEMMeshWidget(QVTKRenderWindowInteractor):
    def __init__(self, parent=None):
        QVTKRenderWindowInteractor.__init__(self, parent)
        self.renderWindow = self.GetRenderWindow()

        self.renderer = vtk.vtkRenderer()

        self.renderWindow.AddRenderer(self.renderer)

        self.renderWindowInteractor = self.renderWindow.GetInteractor()
        self.renderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
        self.renderWindowInteractor.Initialize()
        self.renderWindowInteractor.Start()

        self.om = self.setupOrientationIndicator()
        self.om.SetInteractor(self.renderWindowInteractor)
        self.om.EnabledOn()
        self.om.InteractiveOff()

        self.renderer.AddObserver(vtkCommand.ModifiedEvent, self.on_resize)

        self.renderWindowInteractor.Render()
        self.renderer.ResetCamera()


    def setupOrientationIndicator(self):
        colors = vtk.vtkNamedColors()
        labels = 'xyz'
        colors.SetColor("ParaViewBkg", [82, 87, 110, 255])
        axes = make_cube_actor(labels, colors)
        om = vtk.vtkOrientationMarkerWidget()
        om.SetOrientationMarker(axes)
        om.SetViewport(0, 0, 10.0, 10.0)
        return om

    def on_resize(self, obj, event):
        size = self.renderer.GetSize()
        w = size[0]
        h = size[1]
        if w != 0. and h != 0.:
            px = w - 25*3
            py = h/2 - 250/2;
            self.om.SetViewport(0, 0, 100./w, 100./h)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    frame = QtWidgets.QFrame()

    vl = QtWidgets.QVBoxLayout()
    vtkWidget = VtkFEMMeshWidget()
    vl.addWidget(vtkWidget)
    frame.setLayout(vl)
    frame.show()
    sys.exit(app.exec())
