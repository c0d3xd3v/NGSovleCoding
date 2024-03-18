import vtkmodules.all as vtk
import igl

import Visualization.vtkhelper as vtkhelper

selectedPoints = None
renderer = None

def create_rectangle(xmin, ymin, xmax, ymax):
    # Erstelle Punkte für das Rechteck
    points = vtk.vtkPoints()
    points.InsertNextPoint(xmin, ymin, 0.0)  # Punkt 0
    points.InsertNextPoint(xmin+xmax, ymin, 0.0)  # Punkt 1
    points.InsertNextPoint(xmin+xmax, ymin+ymax, 0.0)  # Punkt 2
    points.InsertNextPoint(xmin, ymin+ymax, 0.0)  # Punkt 3

    # Erstelle Zellen für das Rechteck
    rectangle = vtk.vtkCellArray()
    rectangle.InsertNextCell(4)  # Ein Quadrat hat 4 Punkte
    rectangle.InsertCellPoint(0)
    rectangle.InsertCellPoint(1)
    rectangle.InsertCellPoint(2)
    rectangle.InsertCellPoint(3)

    return points, rectangle

def createVTKPointSet(points_in, ids):
    points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()

    for i in range(ids.GetNumberOfTuples()):
        p = points_in[ids.GetValue(i)]
        #print(p)
        pid = [0]
        pid[0] = points.InsertNextPoint(p)
        vertices.InsertNextCell(1, pid)

    return points, vertices

class MyInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self):
        self.VisibleFilter = None
        self.selected_points_actor = None
        self.selected_points_polydata = None
        self.sv = None
        self.AddObserver('LeftButtonPressEvent', self.OnLeftButtonDown)

    def SetVisibleFilter(self, vis):
        self.VisibleFilter = vis

    def SetSelectedPointsActorData(self, selected_points_actor, selected_points_polydata, sv):
        self.selected_points_actor = selected_points_actor
        self.selected_points_polydata = selected_points_polydata
        self.sv = sv

    def OnLeftButtonDown(self, obj, event):
        self.selected_points_actor.VisibilityOff()
        self.VisibleFilter.Update()
        points = self.VisibleFilter.GetOutput()
        self.selected_points_actor.VisibilityOn()

        ids = points.GetPointData().GetArray("Point_IDs")

        renderer = self.VisibleFilter.GetRenderer()
        # Erstelle eine vtkCoordinate-Instanz für Fensterkoordinaten
        window_coord = vtk.vtkCoordinate()

        cam = renderer.GetActiveCamera()
        print(cam.GetPosition())
        #renderer.Render()
        #camera = renderer.GetActiveCamera()
        #camera
        display_points = vtk.vtkPoints()
        for i in range(ids.GetNumberOfTuples()):
            p = self.sv[ids.GetValue(i)]
            renderer.SetWorldPoint(p[0], p[1], p[2], 1.0)
            renderer.WorldToDisplay()
            dspp = renderer.GetDisplayPoint()
            display_points.InsertNextPoint(dspp[0], dspp[1], 0)
            #print(dspp)
            #print(p)
            #pid = [0]
            #window_coord.SetCoordinateSystemToWorld()
            #window_coord.SetValue(p[0], p[1], p[2])
            #window_coord.SetCoordinateSystemToDisplay()
            #pixel_coords = window_coord.GetComputedValue(renderer)
            #print(pixel_coords)

        vtkpts, vtkvts = createVTKPointSet(self.sv, ids)
        self.selected_points_polydata.SetPoints(vtkpts)
        self.selected_points_polydata.SetVerts(vtkvts)

        #for i in range(ids.GetNumberOfTuples()):
        #    print(ids.GetValue(i))
        print("There are currently:", points.GetNumberOfPoints(), "visible.")

        super().OnLeftButtonDown()

#vtk.vtkObjectFactory.RegisterFactory(MyInteractor)

def main():
    colors = vtk.vtkNamedColors()

    sv, sf = igl.read_triangle_mesh("../data/ridex-Body004.obj")
    triangle_polydata = vtkhelper.iglToVtkPolydata(sf, sv)

    pointid = vtk.vtkIdFilter()
    pointid.SetInputData(triangle_polydata)
    pointid.SetCellIds(False)
    pointid.SetPointIds(True)
    pointid.SetPointIdsArrayName("Point_IDs")
    pointid.Update()

    sphereMapper = vtk.vtkPolyDataMapper()
    sphereMapper.SetInputData(pointid.GetOutput())

    sphereActor = vtk.vtkActor()
    sphereActor.SetMapper(sphereMapper)
    sphereActor.GetProperty().BackfaceCullingOn()
    sphereActor.GetProperty().SetColor(colors.GetColor3d("MistyRose"))

    selectedpoints = vtk.vtkPolyData()
    # Visualize
    selectedpoints_mapper = vtk.vtkPolyDataMapper()
    selectedpoints_mapper.SetInputData(selectedpoints)

    selectedpoints_actor = vtk.vtkActor()
    selectedpoints_actor.SetMapper(selectedpoints_mapper)
    selectedpoints_actor.GetProperty().SetColor(colors.GetColor3d('Tomato'))
    selectedpoints_actor.GetProperty().SetPointSize(5)
    selectedpoints_actor.GetProperty().EdgeVisibilityOn()
    selectedpoints_actor.GetProperty().SetRenderPointsAsSpheres(True)
    #selectedpoints_actor.GetProperty().SetSphereRadius(0.1)
    # Erstelle ein Rechteck PolyData
    points, rectangle = create_rectangle(100, 100, 100, 100)

    rectangle_polydata = vtk.vtkPolyData()
    rectangle_polydata.SetPoints(points)
    rectangle_polydata.SetPolys(rectangle)

    outline_filter = vtk.vtkOutlineFilter()
    outline_filter.SetInputData(rectangle_polydata)
    outline_filter.Update()

    # Erstelle den Mapper und den Actor
    mapper = vtk.vtkPolyDataMapper2D()
    mapper.SetInputData(outline_filter.GetOutput())

    # Erstelle einen vtkActor2D anstelle von vtkActor
    actor = vtk.vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d("Tomato"))

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName("SelectVisiblePoints")

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(sphereActor)
    renderer.AddActor(selectedpoints_actor)
    renderer.AddActor(actor)

    selectVisiblePoints = vtk.vtkSelectVisiblePoints()
    selectVisiblePoints.SetInputData(pointid.GetOutput())
    selectVisiblePoints.SetRenderer(renderer)
    selectVisiblePoints.SelectInvisibleOff()
    #selectVisiblePoints.SetTolerance(1e-3)
    selectVisiblePoints.Update()

    style = MyInteractor()
    renderWindowInteractor.SetInteractorStyle(style)
    style.SetVisibleFilter(selectVisiblePoints)
    style.SetSelectedPointsActorData(selectedpoints_actor, selectedpoints, sv)

    def on_resize(obj, event):
        w = 200
        h = 100
        rw_size = renderWindow.GetSize()
        xmin = int(rw_size[0]/2 - w/2)
        xmax = int(w)
        ymin = int(rw_size[1]/2 - h/2)
        ymax = int(h)
        #print(f'{xmin},{xmax},{ymin},{ymax}')
        points, rectangle = create_rectangle(xmin, ymin, xmax, ymax)
        rectangle_polydata.SetPoints(points)
        rectangle_polydata.SetPolys(rectangle)
        outline_filter.Update()
        selectVisiblePoints.SelectionWindowOn()
        selectVisiblePoints.SetSelection(xmin, xmin+xmax, ymin, ymin+ymax)
        #print(rw_size)
    renderWindow.AddObserver("RenderEvent", on_resize)

    #near_clipping_distance = 0.1
    #far_clipping_distance = 100.0
    #renderer.GetActiveCamera().SetClippingRange(near_clipping_distance, far_clipping_distance)

    renderer.UseHiddenLineRemovalOn()
    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == "__main__":
    main()
