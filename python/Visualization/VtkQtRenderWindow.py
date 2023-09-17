import sys
from PySide2 import QtCore, QtWidgets
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QAction, QToolBar
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from Visualization.VtkModeshapeActor import *
from Visualization.VtkModeshapeScene import *

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent = None):
        QtWidgets.QMainWindow.__init__(self, parent)
 
        self.modeActor = None
        self.pressureActor = None
        self.scene = None

        self.clipOn = False
        self.animateOn = False

        self.toolbar = QToolBar("Toolbar")
        self.addToolBar(self.toolbar)

        toogleClipBA = QAction(QIcon("Visualization/icons/pqClip.svg"), "&toogle clipping", self)
        toogleClipBA.setStatusTip("toogle clipping")
        toogleClipBA.triggered.connect(self.toogleClip)
        toogleClipBA.setCheckable(True)
        self.toolbar.addAction(toogleClipBA)

        self.toogleAnimateBA = QAction(QIcon("Visualization/icons/pqVcrPlay.svg"), "&toogle animate", self)
        self.toogleAnimateBA.setStatusTip("toogle animate")
        self.toogleAnimateBA.triggered.connect(self.toogleAnimate)
        self.toogleAnimateBA.setCheckable(True)
        self.toolbar.addAction(self.toogleAnimateBA)

        self.frame = QtWidgets.QFrame()
 
        self.vl = QtWidgets.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)

        self.renderWindow = self.vtkWidget.GetRenderWindow()

        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.vtkWidget)
        self.show()
        self.iren.Initialize()
        self.iren.CreateRepeatingTimer(1)

    def GetRenderWindow(self):
        return self.renderWindow

    def GetInteractor(self):
        return self.iren

    def setSceneAndActor(self, modeActor, pressureActor, scene):
        self.modeActor = modeActor
        self.pressureActor = pressureActor
        self.scene = scene

    def toogleAnimate(self):
        if(self.modeActor != None and self.scene != None):
            if(self.animateOn == True):
                self.animateOn = False
                self.toogleAnimateBA.setIcon(QIcon("Visualization/icons/pqVcrPlay.svg"))
                self.scene.disableAnimate()
            else:
                self.animateOn = True
                self.toogleAnimateBA.setIcon(QIcon("Visualization/icons/pqVcrPause.svg"))
                self.scene.enableAnimate()

    def toogleClip(self):
        if(self.modeActor != None and self.scene != None):
            if(self.clipOn):
                self.clipOn = False
                self.pressureActor.disableClipping()
                self.scene.disableClippingWidget()
            else:
                self.clipOn = True
                self.pressureActor.enableClipping()
                self.scene.enableClippingWidget()


 