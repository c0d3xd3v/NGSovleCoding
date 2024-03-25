import sys
import ngsolve

from numerics.preconditioning import fe_preconditioning
from numerics.eigensystemsolver import solveEigensystem
from elasticity.eigenfrequencies import build_elasticity_system_on_fes, steel

from Visualization.VtkNGSolve import gfuActor
from Visualization.qt.drawutils import Draw, Draw2

import vtk

path      = sys.argv[1]
mesh      = ngsolve.Mesh(path)
count     = 10

ngsolve.SetNumThreads(8)
with ngsolve.TaskManager():
    solid_fes       = ngsolve.VectorH1(mesh, order=1, complex=True)
    a, b            = build_elasticity_system_on_fes(steel, solid_fes)
    a, b, pre, kapa = fe_preconditioning(solid_fes, a, b, "local")
    gfu, lams       = solveEigensystem(solid_fes, a, b, count, "lobpcg", pre)

print(gfu)
#ngsolve.VTKOutput(ma=mesh, coefs=[gfu.MDComponent(6).real], names=["gfu"], filename='test_file', subdivision=0, legacy=True).Do()

actor, polyData  = gfuActor(mesh, gfu, count, 3)
actor.select_function("eigenmode6")

appendFilter = vtk.vtkAppendFilter()
appendFilter.AddInputData(polyData)
appendFilter.Update()

unstructuredGrid = vtk.vtkUnstructuredGrid()
unstructuredGrid.ShallowCopy(appendFilter.GetOutput())

# Write the unstructured grid.

writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileVersion(vtk.vtkUnstructuredGridWriter.VTK_LEGACY_READER_VERSION_4_2)
writer.SetFileName("UnstructuredGrid.vtk")
writer.SetInputData(unstructuredGrid)
writer.Write()

Draw(actor, periodic_timer=True)
