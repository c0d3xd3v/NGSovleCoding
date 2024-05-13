from mpi4py import MPI
import sys

import vtk
import netgen
import ngsolve
from ngsolve import MPI_Init

from pde.eigensystemsolver import solveEigenmodes
from pde.eigenfrequencies import steel

from Visualization.VtkNGSolve import gfuActor2
from Visualization.qt.drawutils import Draw2


count = 30

comm = MPI.COMM_WORLD
if comm.rank == 0:
    path = sys.argv[1]
    mesh = ngsolve.Mesh(path)
    ngmesh = mesh.ngmesh
    ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
mesh = ngsolve.Mesh(ngmesh)

ngsolve.SetNumThreads(8)
with ngsolve.TaskManager():
    gfu, lams = solveEigenmodes(mesh, steel, 1, "", 15, "slepc_gd")

Draw2(mesh, gfu, len(lams) - 1, periodic_timer=True)
'''
_, polyData = gfuActor2(mesh, gfu, len(lams) - 1)

appendFilter = vtk.vtkAppendFilter()
appendFilter.AddInputData(polyData)
appendFilter.Update()

unstructuredGrid = vtk.vtkUnstructuredGrid()
unstructuredGrid.ShallowCopy(appendFilter.GetOutput())

writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileVersion(vtk.vtkUnstructuredGridWriter.VTK_LEGACY_READER_VERSION_4_2)
writer.SetFileName("UnstructuredGrid.vtk")
writer.SetInputData(unstructuredGrid)
writer.Write()
'''