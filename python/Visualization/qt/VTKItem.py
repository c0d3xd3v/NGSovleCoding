import sys

import vtk
import igl
import numpy as np

from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt, Property, Signal, QTimer, QEvent, QPointF, QObject, Slot, QUrl, QPoint, QCoreApplication, QStringListModel
from PySide6.QtGui import QMatrix4x4, QCursor, QMouseEvent, QWheelEvent, QHoverEvent, QGuiApplication, QSurfaceFormat
from PySide6.QtQuick import QQuickFramebufferObject, QSGSimpleRectNode, QQuickItem
from PySide6.QtOpenGL import QOpenGLFramebufferObject, QOpenGLFramebufferObjectFormat
from PySide6.QtQml import qmlRegisterType, QQmlApplicationEngine

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from Visualization.vtkhelper import *
from Visualization.qt.VTKItemFramebufferRenderer import FbItemRenderer
from Visualization.VtkNGSolve import gfuActor, gfuActor2

import locale

def cloneMouseEvent(event: QMouseEvent):
    return QMouseEvent(
        event.type(),
        event.position(),
        event.scenePosition(),
        event.globalPosition(),
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

def convertToMouseEvent(
    eventType: QEvent.Type,
    localPos: QPointF,
    button: Qt.MouseButton,
    buttons: Qt.MouseButtons,
    modifiers: Qt.KeyboardModifiers,
):
    return QMouseEvent(eventType, localPos, button, buttons, modifiers)

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

    @Slot(float, float, int, int, int)
    def onMousePressed(
        self, x: float, y: float, button: int, buttons: int, modifiers: int
    ):
        self.lastMouseButtonEvent = convertToMouseEvent(
            QEvent.MouseButtonPress,
            QPointF(x, y),
            Qt.MouseButton(button),
            Qt.MouseButtons(buttons),
            Qt.KeyboardModifiers(modifiers),
        )
        self.lastMouseButtonEvent.ignore()
        self.update()

    @Slot(float, float, int, int, int)
    def onMouseReleased(
        self, x: float, y: float, button: int, buttons: int, modifiers: int
    ):
        self.lastMouseButtonEvent = convertToMouseEvent(
            QEvent.MouseButtonRelease,
            QPointF(x, y),
            Qt.MouseButton(button),
            Qt.MouseButtons(buttons),
            Qt.KeyboardModifiers(modifiers),
        )
        self.lastMouseButtonEvent.ignore()
        self.update()

    @Slot(float, float, int, int, int)
    def onMouseMove(
        self, x: float, y: float, button: int, buttons: int, modifiers: int
    ):
        self.lastMouseMoveEvent = convertToMouseEvent(
            QEvent.MouseMove,
            QPointF(x, y),
            Qt.MouseButton(button),
            Qt.MouseButtons(buttons),
            Qt.KeyboardModifiers(modifiers),
        )
        self.lastMouseMoveEvent.ignore()
        self.update()

    @Slot(QPoint, int, int, int, QPoint, float, float)
    def onMouseWheel(
        self,
        angleDelta: QPoint,
        buttons: int,
        inverted: int,
        modifiers: int,
        pixelDelta: QPoint,
        x: float,
        y: float,
    ):
        self.lastWheelEvent = QWheelEvent(
            QPointF(x, y),
            QPointF(x, y),
            pixelDelta,
            angleDelta,
            Qt.MouseButtons(buttons),
            Qt.KeyboardModifiers(modifiers),
            Qt.NoScrollPhase,
            bool(inverted),
        )
        self.lastWheelEvent.ignore()
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


def defaultFormat(stereo_capable):
  """ Po prostu skopiowałem to z https://github.com/Kitware/VTK/blob/master/GUISupport/Qt/QVTKRenderWindowAdapter.cxx
     i działa poprawnie bufor głębokości
  """
  fmt = QSurfaceFormat()
  fmt.setRenderableType(QSurfaceFormat.OpenGL)
  fmt.setVersion(3, 2)
  fmt.setProfile(QSurfaceFormat.CoreProfile)
  fmt.setSwapBehavior(QSurfaceFormat.DoubleBuffer)
  fmt.setRedBufferSize(8)
  fmt.setGreenBufferSize(8)
  fmt.setBlueBufferSize(8)
  fmt.setDepthBufferSize(8)
  fmt.setAlphaBufferSize(8)
  fmt.setStencilBufferSize(0)
  fmt.setStereo(stereo_capable)
  fmt.setSamples(0)

  return fmt

QSurfaceFormat.setDefaultFormat(defaultFormat(False))
qmlRegisterType(VTKItem, "QmlVtk", 1, 0, "VTKItem")
QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
