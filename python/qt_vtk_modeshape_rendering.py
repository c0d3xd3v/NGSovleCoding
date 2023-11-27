import vtkmodules.all as vtk

from Visualization.VtkModeshapeActor import *
from Visualization.VtkPressureFieldActor import *
from Visualization.VtkModeshapeScene import *
from Visualization.VtkQtRenderWindow import *

#import wildmeshing as wm

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow()
    renderWindow = window.GetRenderWindow()
    interactor = window.GetInteractor()

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("mode0.vtu")
    reader.Update()
    modeshapeData = reader.GetOutput()

    reader2 = vtk.vtkXMLUnstructuredGridReader()
    reader2.SetFileName("pressure0.vtu")
    reader2.Update()
    pressureData = reader2.GetOutput()
    pd = pressureData.GetPointData()

    modeshapeActor = VtkModeShapeActor()
    modeshapeActor.setDataset(modeshapeData)

    pressurefieldActor = VtkPressureFieldActor()
    pressurefieldActor.setDataset(pressureData)

    modeshapeScene = VtkModeshapeScene()
    modeshapeScene.setupScene(modeshapeActor, pressurefieldActor)
    modeshapeScene.setClippedActor(pressurefieldActor)
    modeshapeScene.renderer.AddActor(pressurefieldActor.actor)

    renderWindow.AddRenderer(modeshapeScene.renderer)
    window.setSceneAndActor(modeshapeActor, pressurefieldActor, modeshapeScene)

    modeshapeScene.setRenderWindowInteractor(interactor)
    interactor.AddObserver("TimerEvent", modeshapeScene.updateTimer)

    modeshapeActor.disableClipping()
    pressurefieldActor.disableClipping()
    modeshapeScene.disableAnimate()
    modeshapeScene.disableClippingWidget()

    sys.exit(app.exec())
