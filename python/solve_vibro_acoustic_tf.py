from mpi4py import MPI
import sys
import ngsolve
from netgen.stl import *
from netgen.meshing import *


#from meshing.AcousticMesh import *
from pde.eigenfrequencies import *
from pde.boundarysourcesolver import *
from pde.preconditioning import fe_preconditioning
from pde.eigensystemsolver import solveEigensystem, solveEigenmodes

ngsolve.MPI_Init()

count = 30
#comm = MPI_Init()
comm = MPI.COMM_WORLD
if comm.rank == 0:
    path = sys.argv[2]
    mesh = ngsolve.Mesh(path)
    ngmesh = mesh.ngmesh
    ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
ngsmesh = ngsolve.Mesh(ngmesh)

print(path)

SetVisualization(clipping=True, clipnormal=tuple([0., 0., -1.]))

count = 10

ngsolve.SetNumThreads(8)
with ngsolve.TaskManager():
    eigenmodes, lams = solveEigenmodes(ngsmesh, steel, 1,"solid", 15, "slepc_gd")

print(lams)
Draw(eigenmodes)

air_fes = H1(ngsmesh, definedon="air", dirichlet=ngsmesh.Boundaries("solid|fixed"), order=2, complex=True)
n = specialcf.normal(3)
E = eigenmodes.MDComponent(6)
g = BoundaryFromVolumeCF(E)
roh = rohForAir(lams[6])
print("roh : ", roh)
gfu = solveAcousticBoundaryValue(air_fes, n, g, roh, ngsmesh.Boundaries("solid|fixed"))
print("done.")
Draw(gfu)

