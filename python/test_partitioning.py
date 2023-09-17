import igl
import vtkmodules.all as vtk

from Visualization.vtkhelper import *
from meshing.partitioning import *


sv, sf = igl.read_triangle_mesh("../data/skull.msh__sf.obj")
membership = spectralClustering(sf, 5)
print(np.unique(membership))

renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

triangle_polydata = iglToVtkPolydata(sf, sv)
triangle_polydata = addScalarCellData(triangle_polydata, membership)

triangle_mapper = vtk.vtkPolyDataMapper()
triangle_mapper.SetInputData(triangle_polydata)

triangle_actor = vtk.vtkActor()
triangle_actor.SetMapper(triangle_mapper)

renderer.AddActor(triangle_actor)
renderer.ResetCamera()
render_window.Render()

render_window_interactor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

render_window_interactor.Start()
