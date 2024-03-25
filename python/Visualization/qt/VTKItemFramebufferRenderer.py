from PySide6.QtCore import Qt, Signal, QTimer, QEvent, QPointF
from PySide6.QtGui import QMatrix4x4, QCursor, QMouseEvent
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
            self.rwi = QVTKRenderWindowInteractor()
            self.style = None
            self.style = vtk.vtkInteractorStyleTrackballCamera()
            #self.style.SetMotionFactor(25.0)
            self.rwi.SetInteractorStyle(self.style)
            self.rwi.setMouseTracking(True)

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
        self.style.SetMotionFactor(int(size.width()/size.height()*20))
        if self.mouseIn:
            if self.vtkitem.mouseMove:
                self.rwi.MouseMoveEvent()
                self.vtkitem.mouseMove = False
            if self.vtkitem.mouseLeftPress:
                self.rwi.LeftButtonPressEvent()
                self.vtkitem.mouseLeftPress = False
            if self.vtkitem.zoomIn:
                self.rwi.MouseWheelForwardEvent()
                self.vtkitem.zoomIn = False
            elif self.vtkitem.zoomOut:
                self.rwi.MouseWheelBackwardEvent()
                self.vtkitem.zoomOut = False

        self.rw.Render()
        self.rwi.Render()
