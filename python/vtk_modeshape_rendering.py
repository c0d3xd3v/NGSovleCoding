import vtkmodules.all as vtk

from Visualization.VtkModeshapeActor import *
from Visualization.VtkModeshapeScene import *
from Visualization.VtkRenderWindow import *

if __name__ == "__main__":
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("mode0.vtu")
    reader.Update()
    dataset = reader.GetOutput()

    modeshapeActor = VtkModeShapeActor()
    modeshapeActor.setDataset(dataset)

    modeshapeScene = VtkModeshapeScene()
    modeshapeScene.setupScene(modeshapeActor)

    renderWindow = VtkRenderWindow()
    renderWindow.addRenderer(modeshapeScene.renderer)
    renderWindow.AddTimerObserver(modeshapeScene.updateTimer)

    modeshapeScene.setRenderWindowInteractor(renderWindow.GetInteractor())

    renderWindow.start()
