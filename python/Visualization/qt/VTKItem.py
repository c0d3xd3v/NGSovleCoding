import sys

import vtk
import igl

from PySide6.QtCore import Qt, Signal, QTimer, QEvent, QPointF, QObject, Slot, QUrl
from PySide6.QtGui import QMatrix4x4, QCursor, QMouseEvent, QWheelEvent, QHoverEvent, QGuiApplication
from PySide6.QtQuick import QQuickFramebufferObject, QSGSimpleRectNode, QQuickItem
from PySide6.QtOpenGL import QOpenGLFramebufferObject, QOpenGLFramebufferObjectFormat

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from Visualization.vtkhelper import *
from Visualization.qt.VTKItemFramebufferRenderer import FbItemRenderer

import locale


def cloneMouseEvent(event: QMouseEvent):
    return QMouseEvent(
        event.type(),
        event.localPos(),
        event.windowPos(),
        event.screenPos(),
        event.button(),
        event.buttons(),
        event.modifiers(),
        event.source(),
    )

def cloneWheelEvent(event: QWheelEvent):
    return QWheelEvent(
        event.position(),
        event.globalPosition(),
        event.pixelDelta(),
        event.angleDelta(),
        event.buttons(),
        event.modifiers(),
        event.phase(),
        event.inverted(),
        event.source(),
    )

class VTKItem(QQuickFramebufferObject):
    rendererInitialized = Signal()
    def __init__(self):
        super().__init__()
        self.renderer = None
        self.interactor = None
        self.renderWindow = None

        self.lastMouseButtonEvent: QMouseEvent = None
        self.lastMouseMoveEvent: QMouseEvent = None
        self.lastWheelEvent: QWheelEvent = None

        self.colors = vtk.vtkNamedColors()
        labels = 'xyz'

        self.colors.SetColor("ParaViewBkg", [82, 87, 110, 255])
        self.axes = make_cube_actor(labels, self.colors)
        self.om = vtk.vtkOrientationMarkerWidget()
        self.om.SetOrientationMarker(self.axes)

        self.setAcceptHoverEvents(True)
        self.setAcceptedMouseButtons(Qt.AllButtons)

        self.mouseIn = False
        self.zoomIn = False
        self.zoomOut = False
        self.mouseMove = False
        self.mouseLeftPress = False
        self.mouseLeftRelease = False

        QTimer.singleShot(1, lambda : self.object_created())

    def object_created(self):
        if self.renderer != None:
            renderer = self.renderer.renderer
            self.interactor = self.renderer.rwi
            #self.interactor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

            vtksnowwhite = [0.0, 0.0, 0.0]  # black
            vtk.vtkNamedColors().GetColorRGB("snow", vtksnowwhite)
            vtkskyblue = [0.0, 0.0, 0.0]  # black
            vtk.vtkNamedColors().GetColorRGB("skyblue", vtkskyblue)
            renderer.SetBackground(vtksnowwhite)
            renderer.SetBackground2(vtkskyblue)
            renderer.GradientBackgroundOn()
            renderer.SetGradientMode(vtk.vtkViewport.GradientModes.VTK_GRADIENT_VERTICAL)

            self.om.SetInteractor(self.interactor)
            self.om.EnabledOn()
            self.om.InteractiveOff()
            self.om.SetViewport(0, 0, 10.0, 10.0)
            self.om.Modified()

            self.update()
            self.rendererInitialized.emit()
        else:
            QTimer.singleShot(5, lambda : self.object_created())

    def mousePressEvent(self, event):
        event.accept()
        if self.mouseIn:
            self.mouseLeftPress = True
            #self.interactor.mousePressEvent(event)
            #self.interactor.Modified()
        self.update()
        #print(event)

    def mouseReleaseEvent(self, event):
        event.ignore()
        self.mouseLeftRelease = True
        #self.interactor.mousePressEvent(event)
        self.interactor.Modified()
        self.update()
        #print(event)

    def mouseMoveEvent(self, event):
        event.ignore()
        if self.mouseIn:
            self.mouseMove = True
            self.interactor.mousePressEvent(event)
            #self.interactor.Modified()
        self.update()
        #print(event)

    def hoverLeaveEvent(self, event: QHoverEvent):
        event.accept()
        self.mouseIn = False
        #self.interactor.Render()
        self.interactor.mousePressEvent(event)
        self.update()
        #print(event)

    def hoverEnterEvent(self, event: QHoverEvent):
        event.accept()
        self.mouseIn = True
        self.interactor.mousePressEvent(event)
        self.update()
        #print(event)

    def hoverMoveEvent(self, event: QHoverEvent):
        event.accept()
        if self.mouseIn:
            self.mouseMove = True
            pos = event.position()
            self.interactor.mousePressEvent(event)
            #self.interactor.Modified()
        #print(event)

    def wheelEvent(self, event: QWheelEvent):
        event.accept()
        #if self.mouseIn:
        #    self.interactor.wheelEvent(event)
        #    self.interactor.Modified()
        delta = event.angleDelta().y()

        if delta > 0:
            self.zoomIn = True
            self.zoomOut = False
        elif delta < 0:
            self.zoomIn = False
            self.zoomOut = True
        self.update()
        #print(event)

    def clearAllActors(self):
        actors = self.renderer.renderer.GetActors()
        actors.InitTraversal()
        actor = actors.GetNextItem()
        while actor:
            self.renderer.renderer.RemoveActor(actor)
            actor = actors.GetNextItem()

    def updatePaintNode(self, node, inOutData):
        width = self.width()
        height = self.height()
        #print(width, height)
        if width != 0. and height != 0.:
            self.om.SetViewport(0, 0, 100./width, 100./height)

        return QQuickFramebufferObject.updatePaintNode(self, node, inOutData)

    def createRenderer(self):
        self.renderer = FbItemRenderer()
        return self.renderer


if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import Qt, QTimer, QCoreApplication
    from PySide6.QtQml import QQmlApplicationEngine, qmlRegisterType
    from PySide6.QtQuickControls2 import QQuickStyle

    import vtk
    from Visualization.vtkhelper import *

    qmlRegisterType(VTKItem, "QmlVtk", 1, 0, "VTKItem")
    QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
    QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
    QQuickStyle.setStyle("Imagine");
    #sys_argv = sys.argv
    #sys_argv += ["-style", "Windows"]
    app = QApplication()
    engine = QQmlApplicationEngine()
    engine.load('Visualization/qt/main.qml')
    toplevel = engine.rootObjects()[0]
    item = toplevel.findChild(VTKItem, "ConeView")

    def object_created():
        item = toplevel.findChild(VTKItem, "ConeView")
        renderer = item.renderer.renderer

        sphere_source = vtk.vtkSphereSource()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere_source.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        renderer.AddActor(actor)

        renderer.ResetCamera()
        item.update()
    QTimer.singleShot(100, lambda : object_created())

    class MainCtrl(QObject):
        def __init__(self, engine: QQmlApplicationEngine, item : VTKItem):
            super().__init__()
            self.__engine = engine
            self.__vtkitem = item
            ctxt = self.__engine.rootContext()
            self.poly_mesh = None
            ctxt.setContextProperty("MainCtrl", self)
            self.mapper =None
            self.actor = None

        @Slot(str)
        def loadSurfaceMesh(self, path):
            url = QUrl(path).toLocalFile()
            print(url)
            print(self.__vtkitem)

            # Save the current locale
            old_locale = locale.getlocale(locale.LC_NUMERIC)

            # Set the LC_NUMERIC category to "C"
            locale.setlocale(locale.LC_NUMERIC, "C")
            sv, sf = igl.read_triangle_mesh(url)
            print(len(sv), len(sf))
            # Restore the old locale
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
            #print(self.poly_mesh)


    mainctrl = MainCtrl(engine, item)

    app.exec()

