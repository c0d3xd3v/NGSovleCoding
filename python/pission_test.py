# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0
from mpi4py import MPI
from ngsolve import *
from netgen.geom2d import unit_square

import petsc4py.PETSc as psc
import ngsolve.ngs2petsc as n2p


# generate a triangular mesh of mesh-size 0.2
comm = MPI.COMM_WORLD
if comm.rank == 0:
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
    mesh.ngmesh.Distribute(comm)
else:
    mesh = Mesh(netgen.meshing.Mesh.Receive(comm))

# H1-conforming finite element space
fes = H1(mesh, order=3, dirichlet="left|bottom", complex=True)
# define trial- and test-functions
u, v = fes.TnT()

# the right hand side
f = LinearForm(fes)
f += 1 * v * dx

# the bilinear-form
a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx

a.Assemble()
f.Assemble()

psc_mat = n2p.CreatePETScMatrix(a.mat, fes.FreeDofs())
vecmap = n2p.VectorMapping(fes.ParallelDofs(), fes.FreeDofs())
psc_f, psc_u = psc_mat.createVecs()
vecmap.N2P(f.vec, psc_f)

pc = psc.PC()
pc.create()
pc.setType(pc.Type.GAMG)

ksp = psc.KSP().create()
ksp.setType(ksp.Type.GMRES)
ksp.setPC(pc)
def myKSPMonitor(arg0, arg1, arg2):
    print(arg1, arg2)
ksp.setMonitor(myKSPMonitor)

ksp.setOperators(psc_mat)
ksp.setTolerances(rtol=1e-6, atol=0, divtol=1e16, max_it=400)

gfu = GridFunction(fes)
vecmap.N2P(f.vec, psc_f)
ksp.solve(psc_f, psc_u)
vecmap.P2N(psc_u, gfu.vec)

Draw (gfu)
