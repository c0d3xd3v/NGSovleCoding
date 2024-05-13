import locale
import igl
import vtk
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt, Signal, QTimer, QEvent, QObject, Slot, QUrl
from PySide6.QtQml import qmlRegisterType, QQmlApplicationEngine

from Visualization.qt.VTKItem import VTKItem
from Visualization.VtkNGSolve import gfuActor2

class MainCtrl(QObject):
    meshLoaded = Signal()
    def __init__(self):
        super().__init__()
        self.__engine = None
        self.__vtkitem = None
        self.poly_mesh = None
        self.mapper =None
        self.actor = None
        self.function_names = None
        self.animationTimer = QTimer()
        self.animationTimer.timeout.connect(self.animationTimeout)

    def setupInternal(self, engine: QQmlApplicationEngine, item : VTKItem):
        self.__engine = engine
        self.__vtkitem = item
        ctxt = self.__engine.rootContext()

    def reisterToQML(self):
        ctxt = self.__engine.rootContext()
        ctxt.setContextProperty("MainCtrl", self)

    @Slot(str)
    def loadSurfaceMesh(self, path):
        url = QUrl(path).toLocalFile()
        old_locale = locale.getlocale(locale.LC_NUMERIC)
        locale.setlocale(locale.LC_NUMERIC, "C")
        sv, sf = igl.read_triangle_mesh(url)
        locale.setlocale(locale.LC_NUMERIC, old_locale)

        self.poly_mesh = iglToVtkPolydata(sf, sv)
        self.mapper = vtk.vtkOpenGLPolyDataMapper()
        self.mapper.SetInputData(self.poly_mesh)
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.__vtkitem.clearAllActors()
        self.__vtkitem.renderer.renderer.AddActor(self.actor)
        self.__vtkitem.renderer.renderer.ResetCamera()
        self.__vtkitem.update()
        self.meshLoaded.emit()

    def loadMeshActor(self, mesh, gfu, nconv):

        self.actor, self.poly_mesh  = gfuActor2(mesh, gfu, nconv)

        pointData = self.poly_mesh.GetPointData()
        n = pointData.GetNumberOfArrays()
        self.function_names = []
        for k in range(n):
            array_name = pointData.GetArrayName(k)
            self.function_names.append(array_name)

        #self.actor.select_function("eigenmode0")

        self.__vtkitem.renderer.renderer.AddActor(self.actor.actor)
        self.__vtkitem.renderer.renderer.ResetCamera()
        self.__vtkitem.update()

        self.meshLoaded.emit()

        #self.animationTimer.start(100)

    @Slot(result=list)
    def getFunctionNames(self):
        return self.function_names

    @Slot(str)
    def selectFunctionByName(self, fname):
        self.actor.select_function(fname)
        self.__vtkitem.update()

    @Slot(bool)
    def toogleWireframe(self, toogle):
        print("toogle")
        if toogle:
            self.actor.enableEdges()
        else:
            self.actor.disableEdges()
        self.__vtkitem.update()

    @Slot()
    def animationTimeout(self):
        self.actor.setTimer(1)
        self.__vtkitem.update()

def Draw2(mesh, gfu, nconv, periodic_timer=False):
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import Qt, QTimer
    from PySide6.QtQml import QQmlApplicationEngine
    from PySide6.QtQuickControls2 import QQuickStyle
    from Visualization.qt.VTKItem import VTKItem
    #import vtk
    #from Visualization.vtkhelper import *

    QQuickStyle.setStyle("Fusion");
    app = QApplication()
    engine = QQmlApplicationEngine()

    mainctrl = MainCtrl()
    ctxt = engine.rootContext()
    ctxt.setContextProperty("MainCtrl", mainctrl)

    engine.load('Visualization/qt/main.qml')
    toplevel = engine.rootObjects()[0]
    item = toplevel.findChild(VTKItem, "ConeView")

    def object_created():
        mainctrl.setupInternal(engine, item)
        mainctrl.loadMeshActor(mesh, gfu, nconv)

    QTimer.singleShot(1000, lambda : object_created())

    app.exec()
