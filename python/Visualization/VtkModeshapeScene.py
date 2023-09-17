import vtkmodules.all as vtk

from Visualization.vtkhelper import *
from Visualization.VtkModeshapeActor import *

class VtkModeshapeScene():
    def __init__(self):
        self.animateOn = False
        self.scalarBar = vtk.vtkScalarBarActor()
        self.om = vtk.vtkOrientationMarkerWidget()
        self.planeWidget = vtk.vtkImplicitPlaneWidget()
        self.colors = vtk.vtkNamedColors()
        self.renderer = vtk.vtkRenderer()
        self.modeshapeActor = None
        self.timer_count = 0

    def setupScene(self, modeshapeActor : VtkModeShapeActor, pressurefield):
        self.modeshapeActor = modeshapeActor
        self.clippedActor = pressurefield
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(82./255., 87./255., 110./255.)
        self.scalarBar = self.setupScalarBar(modeshapeActor.mapper)
        self.om = self.setupOrientationIndicator()
        self.planeWidget = self.setupClippingIndicator(pressurefield.dataset)
        self.renderer.AddActor(modeshapeActor.actor)
        self.renderer.AddActor2D(self.scalarBar)
        self.renderer.AddObserver(vtk.vtkCommand.ModifiedEvent, self.on_resize)
        camera = self.renderer.GetActiveCamera()
        center = modeshapeActor.GetCenter()
        pos = tuple([center[0], center[1], center[2] -1000.])
        camera.SetPosition(pos)
        camera.SetFocalPoint(center)

    def setRenderWindowInteractor(self, renderWindowInteractor):
        self.om.SetInteractor(renderWindowInteractor)
        self.om.EnabledOn()
        self.om.InteractiveOff()
        self.planeWidget.SetInteractor(renderWindowInteractor)
        #self.planeWidget.On() 

    def enableClippingWidget(self):
        self.planeWidget.On() 

    def disableClippingWidget(self):
        self.planeWidget.Off() 

    def enableAnimate(self):
        self.animateOn = True

    def disableAnimate(self):
        self.animateOn = False

    def setupScalarBar(self, mapper):
        scalarBar = vtk.vtkScalarBarActor()
        scalarBar.SetLookupTable(mapper.GetLookupTable());
        scalarBar.SetTitle("");
        scalarBar.SetNumberOfLabels(4);
        scalarBar.UnconstrainedFontSizeOn()
        scalarBar.SetMaximumHeightInPixels(250)
        scalarBar.SetMaximumWidthInPixels(25)
        return scalarBar

    def setClippedActor(self, clippedActor):
        self.clippedActor = clippedActor
        self.setupClippingIndicator(clippedActor.dataset)

    def setupOrientationIndicator(self):
        colors = vtk.vtkNamedColors()
        labels = 'xyz'
        colors.SetColor("ParaViewBkg", [82, 87, 110, 255])
        axes = make_cube_actor(labels, colors)
        om = vtk.vtkOrientationMarkerWidget()
        om.SetOrientationMarker(axes)
        om.SetViewport(0, 0, 0.2, 0.2)
        return om

    def setupClippingIndicator(self, dataset):
        planeWidget = vtk.vtkImplicitPlaneWidget()
        planeWidget.SetPlaceFactor( 1 )
        planeWidget.OriginTranslationOff()
        planeWidget.OutlineTranslationOff()
        planeWidget.SetInputData(dataset)
        planeWidget.PlaceWidget()
        planeWidget.SetOrigin(dataset.GetCenter())
        planeWidget.GetNormalProperty().SetOpacity(0)
        planeWidget.SetNormal(0, 0, 1)
        planeWidget.TubingOff()
        planeWidget.GetPlaneProperty().SetOpacity(0.1)
        planeWidget.AddObserver("InteractionEvent", self.clipCallback)
        return planeWidget

    def clipCallback(self, obj, event):
        obj.GetPlane(self.clippedActor.clipPlane)
        self.clippedActor.clipPlane.Modified()
        self.clippedActor.clipper.Update()
        self.clippedActor.gridToPolyData.Modified()
        self.clippedActor.gridToPolyData.Update()

    def on_resize(self, obj, event):
        size = self.renderer.GetSize()
        w = size[0]
        h = size[1]
        if w != 0. and h != 0.:
            px = w - 25*3
            py = h/2 - 250/2;
            self.scalarBar.SetPosition(px/w, py/h)
            self.scalarBar.SetPosition2(250., 25.)
            self.om.SetViewport(0, 0, 100./w, 100./h)

    def updateTimer(self, a, b):
        if self.animateOn == True:
            self.timer_count += 1
            win = self.renderer.GetRenderWindow()
            if win != None:
                self.modeshapeActor.setTimer(self.timer_count/100.0)
                self.clippedActor.setTimer(self.timer_count/100.0)
                win.Render()
            self.timer_count = self.timer_count % 100

