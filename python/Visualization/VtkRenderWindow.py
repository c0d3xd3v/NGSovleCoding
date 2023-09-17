import vtkmodules.all as vtk

class VtkRenderWindow():
    def __init__(self):
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self.renderWindowInteractor.SetRenderWindow(self.renderWindow)
        self.renderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
        self.renderWindowInteractor.CreateRepeatingTimer(10)

    def addRenderer(self, renderer):
        self.renderWindow.AddRenderer(renderer)

    def start(self):
        self.renderWindow.Render()
        self.renderWindowInteractor.Start()

    def GetInteractor(self):
        return self.renderWindowInteractor

    def AddTimerObserver(self, observer):
        self.renderWindowInteractor.AddObserver('TimerEvent', observer)
