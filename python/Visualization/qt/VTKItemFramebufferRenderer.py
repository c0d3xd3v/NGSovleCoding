from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt, Signal, QTimer, QEvent, QPointF
from PySide6.QtGui import QMatrix4x4, QCursor, QMouseEvent, QWheelEvent
from PySide6.QtQuick import QQuickFramebufferObject, QSGSimpleRectNode
from PySide6.QtOpenGL import QOpenGLFramebufferObject, QOpenGLFramebufferObjectFormat
import vtk

from vtk import vtkCommand
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from Visualization.vtkhelper import *


class FbItemRenderer(QQuickFramebufferObject.Renderer):
    def __init__(self):
            super().__init__()
            self.rw = None
            self.rwi = None
            self.renderer = None
            self._fbo = None

            self.rw = vtk.vtkGenericOpenGLRenderWindow()
            self.rwi = vtk.vtkGenericRenderWindowInteractor()

            self.style = None
            self.style = vtk.vtkInteractorStyleTrackballCamera()
            #self.style.SetMotionFactor(2.)
            self.rwi.SetInteractorStyle(self.style)
            #self.rwi.SetScale(0.0001)

            self.mouseIn = False
            self.vtkitem = None

            self.rwi.SetRenderWindow(self.rw)
            self.renderer = vtk.vtkRenderer()
            self.rw.AddRenderer(self.renderer)
            self.window = None



    def createFramebufferObject(self, size):
        fmt = QOpenGLFramebufferObjectFormat()
        fmt.setAttachment(QOpenGLFramebufferObject.Depth)
        fbo = QOpenGLFramebufferObject(size, fmt)
        self._fbo = fbo
        self.isinit = False
        return fbo

    def synchronize(self, item: QQuickFramebufferObject):
        size = item.size()
        self.rw.SetSize(int(size.width()), int(size.height()))
        self.rwi.SetSize(int(size.width()), int(size.height()))
        self.mouseIn = item.mouseIn
        self.vtkitem = item

    def render(self):
        # Called with the FBO bound and the viewport set.
        # Issue OpenGL commands here.
        if not self.isinit:
            self.rw.SetOwnContext(False)
            self.rw.OpenGLInitContext()
            self.rw.SetIsCurrent(True)
            self.rw.SetReadyForRendering(True)
            self.rw.FramebufferFlipYOn()
            self.rwi.Initialize()
            self.rwi.Start()
            self.isinit = True

        size = self.vtkitem.size()
        #self.rw.SetSize(int(size.width()), int(size.height()))
        #self.rwi.SetSize(int(1.0*20), int(size.height()/size.width()*20))

        self.rw.SetIsCurrent(True)
        self.rw.SetReadyForRendering(True)
        #self.style.SetMotionFactor(int(size.width()/size.height()*20))

        if self.vtkitem.lastMouseButtonEvent and not self.vtkitem.lastMouseButtonEvent.isAccepted():
            self.__processMouseButtonEvent(self.vtkitem.lastMouseButtonEvent)
            self.vtkitem.lastMouseButtonEvent.accept()

        if self.vtkitem.lastMouseMoveEvent and not self.vtkitem.lastMouseMoveEvent.isAccepted():
            self.__processMouseMoveEvent(self.vtkitem.lastMouseMoveEvent)
            self.vtkitem.lastMouseMoveEvent.accept()

        if self.vtkitem.lastWheelEvent and not self.vtkitem.lastWheelEvent.isAccepted():
            self.__processWheelEvent(self.vtkitem.lastWheelEvent)
            self.vtkitem.lastWheelEvent.accept()

        self.rw.Render()
        self.rwi.Render()

    def __processMouseButtonEvent(self, event: QMouseEvent):
        ctrl, shift = self.__getCtrlShift(event)
        repeat = 0
        if event.type() == QEvent.MouseButtonDblClick:
            repeat = 1

        self.__setEventInformation(event.position(), ctrl, shift, chr(0), repeat, None)
        if (
            event.type() == QEvent.MouseButtonPress
            or event.type() == QEvent.MouseButtonDblClick
        ):
            if event.button() == Qt.LeftButton:
                self.rwi.LeftButtonPressEvent()
            elif event.button() == Qt.RightButton:
                self.rwi.RightButtonPressEvent()
            elif event.button() == Qt.MiddleButton:
                self.rwi.MiddleButtonPressEvent()
        elif event.type() == QEvent.MouseButtonRelease:
            if event.button() == Qt.LeftButton:
                self.rwi.LeftButtonReleaseEvent()
            elif event.button() == Qt.RightButton:
                self.rwi.RightButtonReleaseEvent()
            elif event.button() == Qt.MiddleButton:
                self.rwi.MiddleButtonReleaseEvent()

    def __processMouseMoveEvent(self, event: QMouseEvent):
        ctrl, shift = self.__getCtrlShift(event)
        self.__setEventInformation(event.position(), ctrl, shift, chr(0), 0, None)
        mf = self.style.GetMotionFactor()
        #self.style.SetMotionFactor(3.*mf)
        self.rwi.MouseMoveEvent()
        #self.style.SetMotionFactor(mf)

    def __processWheelEvent(self, event: QWheelEvent):
        ctrl, shift = self.__getCtrlShift(event)
        self.__setEventInformation(event.position(), ctrl, shift, chr(0), 0, None)

        delta = event.angleDelta().y()
        mf = self.style.GetMotionFactor()
        self.style.SetMotionFactor(mf*0.25)
        if delta > 0:
            self.rwi.MouseWheelForwardEvent()
        elif delta < 0:
            self.rwi.MouseWheelBackwardEvent()
        self.style.SetMotionFactor(mf)

    def __setEventInformation(self, positionPoint: QPointF, ctrl, shift, key, repeat=0, keysum=None):
        scale = self.__getPixelRatio()
        #if self.__fbo.mirrorVertically():
        (w, h) = self.rw.GetSize()
        y = h - positionPoint.y()
        x = positionPoint.x()

        self.rwi.SetEventInformation(
            int(round(x * scale)),
            int(round(y * scale)),
            ctrl,
            shift,
            key,
            repeat,
            keysum,
        )

    def __getCtrlShift(self, event):
        ctrl = shift = False

        if hasattr(event, "modifiers"):
            if event.modifiers() & Qt.ShiftModifier:
                shift = True
            if event.modifiers() & Qt.ControlModifier:
                ctrl = True
        else:
            if self.__saveModifiers & Qt.ShiftModifier:
                shift = True
            if self.__saveModifiers & Qt.ControlModifier:
                ctrl = True

        return ctrl, shift

    def __getPixelRatio(self):
        # Source: https://stackoverflow.com/a/40053864/3388962
        pos = QCursor.pos()
        for screen in QApplication.screens():
            rect = screen.geometry()
            if rect.contains(pos):
                return screen.devicePixelRatio()
        # Should never happen, but try to find a good fallback.
        return QApplication.devicePixelRatio()
