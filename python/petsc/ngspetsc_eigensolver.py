import sys

from mpi4py import MPI
comm = MPI.COMM_WORLD

import ngsolve
import netgen
from ngsPETSc import EigenSolver

from Visualization.VtkNGSolve import gfuActor
from Visualization.qt.drawutils import Draw, Draw2

from elasticity.eigenfrequencies import build_elasticity_system_on_fes, steel


if comm.rank == 0:
    path      = sys.argv[1]
    mesh      = ngsolve.Mesh(path).ngmesh
    mesh.Distribute(comm)
else:
    mesh = netgen.meshing.Mesh.Receive(comm)
mesh = ngsolve.Mesh(mesh)

fes = ngsolve.VectorH1(mesh, order=2, complex=True)
a, b = build_elasticity_system_on_fes(steel, fes)
a.Assemble()
b.Assemble()

solverParameters={"eps_type":"krylovschur", "pc_type":"bddc", "st_type":"sinvert"}
solver = EigenSolver((a, b), fes, 20)
solver.solve()

for i in range(20):
    print(solver.eigenValue(i))
eigenMode, _ = solver.eigenFunction(7)

if comm.rank == 0:
    actor  = gfuActor(mesh, eigenMode, 1, 3)
    actor.select_function("eigenmode0")
    Draw(actor, periodic_timer=True)
