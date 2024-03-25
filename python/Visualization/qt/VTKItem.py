from PySide6.QtCore import Qt, Signal, QTimer, QEvent, QPointF
from PySide6.QtGui import QMatrix4x4, QCursor, QMouseEvent, QWheelEvent
from PySide6.QtQuick import QQuickFramebufferObject, QSGSimpleRectNode, QQuickItem
from PySide6.QtOpenGL import QOpenGLFramebufferObject, QOpenGLFramebufferObjectFormat
import vtk

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from Visualization.vtkhelper import *
from Visualization.qt.VTKItemFramebufferRenderer import FbItemRenderer

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

        QTimer.singleShot(1, lambda : self.object_created())

    def object_created(self):
        if self.renderer != None:
            renderer = self.renderer.renderer
            self.interactor = self.renderer.rwi
            self.renderWindow = self.renderer.rw

            self.interactor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

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
            QTimer.singleShot(1, lambda : self.object_created())

    def mousePressEvent(self, event):
        event.accept()
        self.interactor.mousePressEvent(event)

        self.update()
        #print(event)

    def mouseReleaseEvent(self, event):
        event.accept()
        self.interactor.mousePressEvent(event)
        self.update()
        #print(event)

    def mouseMoveEvent(self, event):
        event.accept()
        self.interactor.mousePressEvent(event)
        self.update()
        #print(event)

    def wheelEvent(self, event: QWheelEvent):
        event.accept()
        self.interactor.wheelEvent(event)
        delta = event.angleDelta().y()
        if delta > 0:
            self.interactor.MouseWheelForwardEvent()
        elif delta < 0:
            self.interactor.MouseWheelBackwardEvent()
        self.update()
        #print(event)

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

    import vtk
    from Visualization.vtkhelper import *

    qmlRegisterType(VTKItem, "QmlVtk", 1, 0, "VTKItem")
    QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
    QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

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
    #item.rendererInitialized.connect(lambda : object_created(), Qt.QueuedConnection)

    app.exec()

