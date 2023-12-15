import vtk

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

    # Erstelle PolyData und füge Punkte und Zellen hinzu
    rectangle_polydata = vtk.vtkPolyData()
    rectangle_polydata.SetPoints(points)
    rectangle_polydata.SetPolys(rectangle)

    return rectangle_polydata

def main():
    colors = vtk.vtkNamedColors()

    # Erstelle ein Rechteck PolyData
    rectangle_polydata = create_rectangle(100, 100, 100, 100)

    # Erstelle den Mapper und den Actor
    mapper = vtk.vtkPolyDataMapper2D()
    mapper.SetInputData(rectangle_polydata)

    # Erstelle einen vtkActor2D anstelle von vtkActor
    actor = vtk.vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d("Tomato"))

    # Erstelle die Szene mit Renderer und RenderWindow
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    # Erstelle den RenderWindowInteractor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Füge den Actor zur Szene hinzu
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d("LightGray"))

    # Setze die Kameraansicht
    renderer.ResetCamera()
#    renderer.GetActiveCamera().Zoom(1.2)

    # Setze den Fenstertitel und starte die Interaktion
    render_window.SetWindowName("VTK Rectangle Example (2D Actor)")
    render_window.Render()
    render_window_interactor.Start()

if __name__ == "__main__":
    main()
