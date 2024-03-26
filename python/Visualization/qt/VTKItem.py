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

def create_mouse_event_from_hover_event(hover_event: QHoverEvent) -> QMouseEvent:
    # Extract relevant information from the hover event
    pos = hover_event.pos()
    global_pos = hover_event.scenePosition()  # Global position
    button_state = Qt.NoButton  # Since it's a hover event, no button is pressed
    modifiers = hover_event.modifiers()
    buttons = Qt.MouseButton()

    # Create a QMouseEvent with the extracted information
    mouse_event = QMouseEvent(
        QEvent.MouseMove,
        pos,
        global_pos,
        button_state,
        buttons,
        modifiers
    )

    return mouse_event

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

        QTimer.singleShot(10, lambda : self.object_created())

    def object_created(self):
        if self.renderer != None:
            renderer = self.renderer.renderer
            self.interactor = self.renderer.rwi

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
            QTimer.singleShot(15, lambda : self.object_created())

    def mousePressEvent(self, event: QMouseEvent):
        self.__processMouseButtonEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent):
        self.__processMouseButtonEvent(event)

    def __processMouseButtonEvent(self, event: QMouseEvent):
        self.lastMouseButtonEvent = cloneMouseEvent(event)
        self.lastMouseButtonEvent.ignore()
        event.accept()
        self.update()

    def mouseMoveEvent(self, event):
        self.lastMouseMoveEvent = cloneMouseEvent(event)
        self.lastMouseMoveEvent.ignore()
        event.accept()
        self.update()

    def wheelEvent(self, event: QWheelEvent):
        self.lastWheelEvent = cloneWheelEvent(event)
        self.lastWheelEvent.ignore()
        event.accept()
        self.update()

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
    QQuickStyle.setStyle("Material");
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
    QTimer.singleShot(10, lambda : object_created())

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

            # Save the current locale
            old_locale = locale.getlocale(locale.LC_NUMERIC)

            # Set the LC_NUMERIC category to "C"
            locale.setlocale(locale.LC_NUMERIC, "C")
            sv, sf = igl.read_triangle_mesh(url)
            print("sv, sf : ", len(sv), len(sf))
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

