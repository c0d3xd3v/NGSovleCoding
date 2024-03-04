
def Draw(femActor, periodic_timer=False):
    import sys
    from PySide6 import QtCore, QtWidgets
    from PySide6.QtGui import QIcon, QAction
    from PySide6.QtWidgets import QToolBar
    from Visualization.qt.VtkFEMMeshWidget import VtkFEMMeshWidget
    app = QtWidgets.QApplication(sys.argv)
    frame = QtWidgets.QFrame()
    vl = QtWidgets.QVBoxLayout()
    vtkWidget = VtkFEMMeshWidget()
    vl.addWidget(vtkWidget)
    frame.setLayout(vl)
    frame.show()
    vtkWidget.renderer.AddActor(femActor.actor)
    vtkWidget.renderer.ResetCamera()
    timerId = vtkWidget.renderWindowInteractor.CreateRepeatingTimer(500)
    if periodic_timer: vtkWidget.registerTimerRequestForActor(femActor)
    vtkWidget.renderWindowInteractor.Start()
    app.exec()
    qApp.shutdown()

def Draw2(femActor, periodic_timer=False):
    from PySide6.QtWidgets import QApplication
    from PySide6.QtQuick import QQuickView, QQuickItem
    from PySide6.QtQml import QQmlApplicationEngine
    from vtkmodules.qt import QVTKRenderWindowInteractor

    app = QApplication()
    engine = QQmlApplicationEngine()
    ctx = engine.rootContext()
    engine.load('Visualization/main.qml')

    toplevel = engine.rootObjects()
    win = toplevel[0]
    win.show()
    item = win.findChild(QQuickItem, "ConeView")
    print(item)
    sys.exit(app.exec())
