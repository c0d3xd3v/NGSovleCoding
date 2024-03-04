import sys
import PySide6
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
        vtksnowwhite = [0.0, 0.0, 0.0]  # black
        vtk.vtkNamedColors().GetColorRGB("snow", vtksnowwhite)
        vtkskyblue = [0.0, 0.0, 0.0]  # black
        vtk.vtkNamedColors().GetColorRGB("skyblue", vtkskyblue)
        self.renderer.SetBackground(vtksnowwhite)
        self.renderer.SetBackground2(vtkskyblue)
        self.renderer.GradientBackgroundOn()
        self.renderer.SetGradientMode(vtk.vtkViewport.GradientModes.VTK_GRADIENT_VERTICAL)

        self.renderWindow.AddRenderer(self.renderer)

        self.renderWindowInteractor = self.renderWindow.GetInteractor()
        self.renderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
        self.renderWindowInteractor.Initialize()
        self.renderWindowInteractor.Start()

        self.om = self.setupOrientationIndicator()
        self.om.SetInteractor(self.renderWindowInteractor)
        self.om.EnabledOn()
        self.om.InteractiveOff()

        self.renderWindowInteractor.Render()
        self.renderer.ResetCamera()

        self.timerRequestObjects = []
        self.timer_count = 0
        def updateTimer(a, b):
            self.timer_count += 1
            win = self.renderer.GetRenderWindow()
            if win != None:
                for actor in self.timerRequestObjects:
                    actor.setTimer(self.timer_count/100.0)
                win.Render()
            self.timer_count = self.timer_count % 100
        self.renderWindowInteractor.AddObserver("TimerEvent", updateTimer)

    def registerTimerRequestForActor(self, actor):
        self.timerRequestObjects.append(actor)

    def setupOrientationIndicator(self):
        colors = vtk.vtkNamedColors()
        labels = 'xyz'
        colors.SetColor("ParaViewBkg", [82, 87, 110, 255])
        axes = make_cube_actor(labels, colors)
        om = vtk.vtkOrientationMarkerWidget()
        om.SetOrientationMarker(axes)
        om.SetViewport(0, 0, 10.0, 10.0)
        return om

    def resizeEvent(self, ev):
        QVTKRenderWindowInteractor.resizeEvent(self, ev)
        size = self.renderer.GetSize()
        w = size[0]
        h = size[1]
        if w != 0. and h != 0.:
            self.om.SetViewport(0, 0, 100./w, 100./h)
            self.renderWindowInteractor.Render()

    def changeEvent(self, ev):
        QVTKRenderWindowInteractor.changeEvent(self, ev)

    def showEvent(self, ev):
        QVTKRenderWindowInteractor.showEvent(self, ev)
        self.renderWindowInteractor.Render()
        self.renderer.ResetCamera()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    frame = QtWidgets.QFrame()

    vl = QtWidgets.QVBoxLayout()
    vtkWidget = VtkFEMMeshWidget()
    vl.addWidget(vtkWidget)
    frame.setLayout(vl)
    frame.show()

    sys.exit(app.exec())
