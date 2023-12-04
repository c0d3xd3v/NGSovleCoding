import vtkmodules.all as vtk

class MyInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self):
        self.VisibleFilter = None
        self.AddObserver('LeftButtonPressEvent', self.OnLeftButtonDown)

    def SetVisibleFilter(self, vis):
        self.VisibleFilter = vis

    def OnLeftButtonDown(self, obj, event):
        self.VisibleFilter.Update()
        print("There are currently:", self.VisibleFilter.GetOutput().GetNumberOfPoints(), "visible.")
        super().OnLeftButtonDown()

#vtk.vtkObjectFactory.RegisterFactory(MyInteractor)

def main():
    colors = vtk.vtkNamedColors()

    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetCenter(5.0, 0, 0)
    sphereSource.Update()

    pointSource = vtk.vtkPointSource()
    pointSource.SetRadius(2.0)
    pointSource.SetNumberOfPoints(200)
    pointSource.Update()

    sphereMapper = vtk.vtkPolyDataMapper()
    sphereMapper.SetInputConnection(sphereSource.GetOutputPort())

    sphereActor = vtk.vtkActor()
    sphereActor.SetMapper(sphereMapper)
    sphereActor.GetProperty().SetColor(colors.GetColor3d("MistyRose"))

    pointsMapper = vtk.vtkPolyDataMapper()
    pointsMapper.SetInputConnection(pointSource.GetOutputPort())

    pointsActor = vtk.vtkActor()
    pointsActor.SetMapper(pointsMapper)
    pointsActor.GetProperty().SetColor(colors.GetColor3d("Honeydew"))

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName("SelectVisiblePoints")

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(sphereActor)
    renderer.AddActor(pointsActor)
    renderer.SetBackground(colors.GetColor3d("ivory_black"))

    selectVisiblePoints = vtk.vtkSelectVisiblePoints()
    selectVisiblePoints.SetInputConnection(pointSource.GetOutputPort())
    selectVisiblePoints.SetRenderer(renderer)
    selectVisiblePoints.Update()

    style = MyInteractor()
    renderWindowInteractor.SetInteractorStyle(style)
    style.SetVisibleFilter(selectVisiblePoints)

    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == "__main__":
    main()
