from mpi4py import MPI
import ngsolve, sys, time
from ngsolve import *
import netgen.csg
import numpy as np

import petsc4py.PETSc as psc
import slepc4py.SLEPc as spc

from Visualization.qt.VtkFEMMeshWidget import showTest
from Visualization.VtkNGSolve import gfuActor
from Visualization.qt.drawutils import Draw, Draw2

from elasticity.eigenfrequencies import build_elasticity_system_on_fes, steel

comm = MPI.COMM_WORLD
rank = comm.rank
if comm.rank == 0:
    path      = sys.argv[1]
    mesh      = ngsolve.Mesh(path)
    ngmesh = mesh.ngmesh #netgen.csg.unit_cube.GenerateMesh(maxh=0.05)
    ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
mesh = Mesh(ngmesh)

fes = VectorH1(mesh, order=2, complex=True)
u,v = fes.TnT()
a, b = build_elasticity_system_on_fes(steel, fes)
a.Assemble()
b.Assemble()
################################################################################
def ngs2petscMat(fes, a):
    pardofs = fes.ParallelDofs()
    globnums, nglob = pardofs.EnumerateGlobally()
    locmat = a.mat.local_mat
    val,col,ind = locmat.CSR()
    ind = np.array(ind, dtype='int32')

    apsc_loc = psc.Mat().createAIJ(size=(locmat.height, locmat.width), csr=(ind,col,val), comm=MPI.COMM_SELF)

    IS = psc.IS().createBlock(bsize=1, indices=globnums, comm=comm)
    lgmap = psc.LGMap().create(bsize=1, indices=globnums, comm=comm)

    mat = psc.Mat().createPython(size=nglob, comm=comm)
    mat.setType(psc.Mat.Type.IS)
    mat.setLGMap(lgmap)
    mat.setISLocalMat(apsc_loc)
    mat.assemble()

    return mat, IS, lgmap
################################################################################

mat, IS, _ = ngs2petscMat(fes, a)
mat_b, IS, _ = ngs2petscMat(fes, b)

pc = psc.PC()
pc.create()
pc.setType(pc.Type.BDDC)

ksp = psc.KSP()
ksp.create()

ksp.setType(ksp.Type.CG)
ksp.setTolerances(rtol=1E-3)
ksp.setPC( pc )

F = spc.ST().create()
F.setType(F.Type.SINVERT)
F.setKSP( ksp )

start = time.time()
eps = spc.EPS()
eps.create()
eps.setST(F)
eps.setOperators(mat_b, mat)
eps.setType(spc.EPS.Type.KRYLOVSCHUR)
eps.setProblemType(spc.EPS.ProblemType.GHEP)
eps.setDimensions(10, spc.DECIDE)
eps.setFromOptions()
eps.solve()

end = time.time()
Print = psc.Sys.Print
its = eps.getIterationNumber()
eps_type = eps.getType()
nev, ncv, mpd = eps.getDimensions()
tol, maxit = eps.getTolerances()
nconv = eps.getConverged()

gfu = ngsolve.GridFunction(fes, multidim=nconv)
if nconv > 0:
    # Create the results vectors
    vr, wr = mat.getVecs()
    vi, wi = mat.getVecs()

    for i in range(nconv):
        real = mat.createVecLeft()
        imag = mat.createVecLeft()
        k = eps.getEigenpair(i, real, imag)
        v1loc = vr.getSubVector(IS)
        for j in range(len(gfu.vecs[i])):
            gfu.vecs[i].FV()[j] = real.getArray()[j]
        #print(len(gfu.vecs[i]), len(v1loc.getArray()))
        #print(k)

if comm.rank == 0:
    '''
    Print()
    Print("******************************")
    Print("*** SLEPc Solution Results ***")
    Print("******************************")
    Print()
    print("full mesh ", mesh.nv)
    Print("Number of iterations of the method: %d" % its)
    Print("Solution method: %s" % eps_type)
    Print("Number of requested eigenvalues: %d" % nev)
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    Print("Number of converged eigenpairs %d" % nconv)
    print(rank, mesh.nv, mesh.nface)
    print ("rank ", comm.rank,"global dofs =", fes.ndofglobal, ", local dofs =", fes.ndof, ", sum of local dofs =", comm.allreduce(fes.ndof))
    '''
    actor  = gfuActor(mesh, gfu, 10)
    actor.select_function("eigenmode6")
    Draw(actor, periodic_timer=True)
