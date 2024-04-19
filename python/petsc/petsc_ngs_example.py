import time
from mpi4py import MPI
#from ngsolve import *
import ngsolve
import numpy as np
from netgen.csg import unit_cube
from ngsolve.krylovspace import CGSolver

import petsc4py
import slepc4py

import netgen.csg as csg
from netgen.meshing import *
from netgen.libngpy._meshing import Point3d
from netgen.geom2d import unit_square
from netgen.occ import *
from netgen.geom2d import SplineGeometry

import petsc4py.PETSc as psc
import ngsolve.ngs2petsc as n2p

def NgsSolve(a, f, fes):
    gfu = GridFunction(fes)
    inv = CGSolver(a.mat, freedofs=fes.FreeDofs(), printing=False, maxiter=4000, tol=1e-16)
    gfu.vec.data = inv * f.vec
    return gfu

def NgsSolveDirect(a, f, fes):
    gfu = GridFunction(fes)
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec
    return gfu

def KrylovSolve(a, f, fes):
    psc_mat = n2p.CreatePETScMatrix(a.mat, fes.FreeDofs())
    vecmap = n2p.VectorMapping(fes.ParallelDofs(), fes.FreeDofs())
    psc_f, psc_u = psc_mat.createVecs()
    vecmap.N2P(f.vec, psc_f)

    ksp = psc.KSP(petsc_options={
    "ksp_type": "gmres",
    #lu
    #  mkl_pardiso
    #  mkl_cpardiso
    #  klu
    #  umfpack
    #  mumps
    #  superlu
    #cholesky
    #  mumps
    #  cholmod
    "pc_type": "gamg",
    "pc_factor_mat_solver_type": "mkl_pardiso"
    })
    ksp.create()
    ksp.setOperators(psc_mat)
    ksp.setTolerances(rtol=1e-10, atol=0, divtol=1e16, max_it=500)
    ksp.solve(psc_f, psc_u)

    gfu = ngsolve.GridFunction(fes)
    vecmap.P2N(psc_u, gfu.vec)
    return gfu

#if __name__ == "__main__":

petsc4py.init()
slepc4py.init()


comm = MPI.COMM_WORLD
if comm.rank == 0:
    '''
    air = Circle((0.5, 0.5), 0.6).Face()
    air.edges.name = 'outer'
    scatterer = MoveTo(0.7, 0.3).Rectangle(0.05, 0.4).Face()
    scatterer.edges.name = 'scat'
    geo = OCCGeometry(air - scatterer, dim=2)
    ngmesh = geo.GenerateMesh(maxh=0.005)
    mesh = Mesh(ngmesh)
    ngmesh = mesh.ngmesh.Distribute(comm)
    '''
    '''
    #ngmesh = unit_cube.GenerateMesh(maxh=0.05).Distribute(comm)
    brick = csg.OrthoBrick(Point3d(-0.5, -0.5, -0.01), Point3d(0.5, 0.5, 0.01))
    brick2 = csg.OrthoBrick(Point3d(0.25, -0.25, -0.005), Point3d(0.25, 0.25, 0.005))
    brick.bc('outer')
    brick2.bc('outer')
    bb = brick - brick2
    bb.bc('outer')
    geo = netgen.csg.CSGeometry()
    geo.Add(bb)
    ngmesh = geo.GenerateMesh(maxh=0.005).Distribute(comm)
    '''

    '''
    air = csg.Sphere((0., 0., 0.,), 1.0)
    air.bc('outer')
    brick = csg.OrthoBrick(Point3d(-0.25, -0.125, -0.01), Point3d(0.25, 0.125, 0.01))
    geo = netgen.csg.CSGeometry()
    bb = air - brick
    geo.Add(bb)
    ngmesh = geo.GenerateMesh(maxh=0.025).Distribute(comm)
    '''
    ngmesh = ngsolve.Mesh('../data/aoustic-ngs-example-3d.vol').ngmesh
    ngmesh.Distribute(comm)

else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)

mesh = ngsolve.Mesh(ngmesh)

ngsolve.Draw(mesh)


omega = 40.
pulse = ngsolve.CoefficientFunction(5e4*ngsolve.exp(-(40**2)*((ngsolve.x-0.0)*(ngsolve.x-0.0) + (ngsolve.y-0.0)*(ngsolve.y-0.0) + (ngsolve.z-0.15)*(ngsolve.z-0.15))))
ngsolve.Draw(pulse, mesh, "pulse")

#fes = H1(mesh, order=3, complex=True, dirichlet=".*")
fes = ngsolve.H1(mesh, order=1, complex=True)

u,v = fes.TnT()
#a = BilinearForm(grad(u)*grad(v)*dx+u*v*ds).Assemble()
#f = LinearForm(1*v*dx).Assemble()

a = ngsolve.BilinearForm(fes)
a += ngsolve.grad(u)*ngsolve.grad(v)*ngsolve.dx - omega**2*u*v*ngsolve.dx
a += -omega*1j*u*v * ngsolve.ds('outer')
a.Assemble()

f = ngsolve.LinearForm(pulse * v * ngsolve.dx).Assemble();

ngs_cg_start_t = time.time()
#gfu = NgsSolve(a, f, fes)
#gfu = NgsSolveDirect(a, f, fes)
ngs_cg_durration = time.time() - ngs_cg_start_t

petsc_cg_start_t = time.time()
gfu = KrylovSolve(a, f, fes)
petsc_cg_durration = time.time() - petsc_cg_start_t

ngsolve.Draw(gfu)

print("ngs_cg_durration   : ", ngs_cg_durration)
print("petsc_cg_durration : ", petsc_cg_durration)

