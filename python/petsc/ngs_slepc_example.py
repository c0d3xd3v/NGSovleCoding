from mpi4py import MPI
import ngsolve, sys, time
import netgen
import vtk

import petsc4py.PETSc as psc

from Visualization.VtkNGSolve import gfuActor2
from Visualization.qt.drawutils import Draw2

from pde.eigenfrequencies import build_elasticity_system_on_fes, steel
from pde.slepc_solvers import SLEPcLOBPCG, SLEPcGD, SLEPcKrylovSchur

comm = MPI.COMM_WORLD
if comm.rank == 0:
    path = sys.argv[1]
    mesh = ngsolve.Mesh(path)
    ngmesh = mesh.ngmesh
    ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
mesh = ngsolve.Mesh(ngmesh)

fes = ngsolve.VectorH1(mesh, order=1, complex=True)
u,v = fes.TnT()
a, b = build_elasticity_system_on_fes(steel, fes)
a.Assemble()
b.Assemble()

gfu, eps = SLEPcGD(a, b, fes)

nconv = eps.getConverged()
its = eps.getIterationNumber()
eps_type = eps.getType()
nev, ncv, mpd = eps.getDimensions()
tol, maxit = eps.getTolerances()

Print = psc.Sys.Print
Print()
Print("******************************")
Print("*** SLEPc Solution Results ***")
Print("******************************")
Print("Number of iterations of the method: %d" % its)
Print("Solution method                   : %s" % eps_type)
Print("Number of requested eigenvalues   : %d" % nev)
Print("Stopping condition                : tol=%.4g, maxit=%d" % (tol, maxit))
Print("Number of converged eigenpairs    : %d" % nconv)


if comm.rank == 0:

    Draw2(mesh, gfu, nconv, periodic_timer=True)
    _, polyData  = gfuActor2(mesh, gfu, nconv)

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


