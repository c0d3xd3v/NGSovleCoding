from mpi4py import MPI
import ngsolve, sys, time
import netgen
import numpy as np
import vtk

import petsc4py.PETSc as psc
import slepc4py.SLEPc as spc
import ngsolve.ngs2petsc as n2p

from Visualization.VtkNGSolve import gfuActor, gfuActor2
from Visualization.qt.drawutils import Draw, Draw2

from elasticity.eigenfrequencies import build_elasticity_system_on_fes, steel

ksp_itrate=0

def SLEPcLOBPCG(a, b, fes):
    mat = n2p.CreatePETScMatrix(a.mat, fes.FreeDofs())
    mat_b = n2p.CreatePETScMatrix(b.mat, fes.FreeDofs())

    pc = psc.PC()
    pc.create()
    pc.setType(pc.Type.NONE)
    pc.setGAMGType(psc.PC.GAMGType.AGG)
    pc.setGAMGLevels(15)
    pc.setFactorSolverType(psc.Mat.SolverType.MKL_PARDISO)

    ksp = psc.KSP()
    ksp.create()
    ksp.setType(ksp.Type.GMRES)
    ksp.setGMRESRestart(20)
    ksp.setTolerances(rtol=1e-3, atol=0., max_it=100)
    ksp.setPC(pc)

    def myKSPMonitor(arg0, arg1, arg2):
        global ksp_itrate
        #print(arg1, arg2, arg3, arg4)
        ksp_itrate = arg1
        #if arg2 < 1e-2:
        print(arg1, arg2)

        #print(arg0.getConvergenceTest())
    ksp.setMonitor(myKSPMonitor)
    print(ksp.getConvergenceTest())

    st = spc.ST().create()
    st.setType(spc.ST.Type.PRECOND)
    st.setShift(40000.)
    st.setKSP(ksp)

    eps = spc.EPS()
    eps.create()

    #eps.setLOBPCGBlockSize(40)
    #eps.setLOBPCGLocking(True)
    eps.setType(spc.EPS.Type.LOBPCG)
    eps.setInterval(1000000., 1000000000.)
    eps.setTolerances(tol=1e-8, max_it=125)
    eps.setProblemType(spc.EPS.ProblemType.GHEP)
    eps.setDimensions(30, spc.DECIDE, spc.DECIDE)

    eps.setST(st)
    print("set operators ")
    eps.setOperators(mat, mat_b)

    def myMonitor(arg0, arg1, arg2, arg3, arg4):
        global ksp_itrate
        #print(arg1, arg2, arg3, arg4)
        print(ksp_itrate, arg1, arg2, arg3)
        print("------------------")
        #print(arg0.getConvergenceTest())
    eps.setMonitor(myMonitor)

    eps.setConvergenceTest(eps.Conv.NORM)
    eps.setWhichEigenpairs(spc.EPS.Which.SMALLEST_REAL)
    #eps.setTarget(100000)

    eps.solve()

    nconv = eps.getConverged()
    gfu = ngsolve.GridFunction(fes, multidim=nconv)
    vecMap = n2p.VectorMapping(fes.ParallelDofs(), fes.FreeDofs())
    for k in range(nconv):
        eignModePETScReal = mat.createVecLeft()
        eignModePETScImag = mat.createVecLeft()
        eps.getEigenvector(k, eignModePETScReal, eignModePETScImag)
        lam = eps.getEigenvalue(k)
        vecMap.P2N(eignModePETScReal, gfu.vecs[k])
        print(lam)

    return gfu, eps

def SLEPcGD(a, b, fes):
    max_it = 100
    req_ep = 14

    mat = n2p.CreatePETScMatrix(a.mat, fes.FreeDofs())
    mat_b = n2p.CreatePETScMatrix(b.mat, fes.FreeDofs())

    pc = psc.PC()
    pc.create()
    pc.setType(pc.Type.LU)
    pc.setFactorSolverType(psc.Mat.SolverType.MKL_PARDISO)

    ksp = psc.KSP()
    ksp.create()
    ksp.setType(ksp.Type.PREONLY)
    ksp.setPC(pc)
    def myKSPMonitor(arg0, arg1, arg2):
        #global ksp_itrate
        #print(arg1, arg2, arg3, arg4)
        #ksp_itrate = arg1
        #if arg2 < 1e-2:
        print(arg1, arg2)

        #print(arg0.getConvergenceTest())
    ksp.setMonitor(myKSPMonitor)

    st = spc.ST().create()
    st.setType(spc.ST.Type.PRECOND)
    st.setShift(40000.)
    st.setKSP(ksp)

    eps = spc.EPS()
    eps.create()

    eps.setType(spc.EPS.Type.GD)
    eps.setTolerances(tol=1e-10, max_it=max_it)
    eps.setProblemType(spc.EPS.ProblemType.GHEP)
    eps.setDimensions(req_ep, spc.DECIDE, spc.DECIDE)

    eps.setST(st)
    eps.setOperators(mat, mat_b)

    def myMonitor(arg0, arg1, arg2, arg3, arg4):
        #print(arg1, arg2, arg3, arg4)
        print("iteraions :", arg1, "/", max_it, " convered ep :",  arg2, "/", req_ep)
        #print(arg0.getConvergenceTest())
    eps.setMonitor(myMonitor)

    eps.setConvergenceTest(eps.Conv.NORM)
    eps.setWhichEigenpairs(spc.EPS.Which.SMALLEST_REAL)
    eps.setTarget(100000)

    eps.solve()

    nconv = eps.getConverged()
    gfu = ngsolve.GridFunction(fes, multidim=nconv)
    vecMap = n2p.VectorMapping(fes.ParallelDofs(), fes.FreeDofs())
    for k in range(nconv):
        eignModePETScReal = mat.createVecLeft()
        eignModePETScImag = mat.createVecLeft()
        eps.getEigenvector(k, eignModePETScReal, eignModePETScImag)
        lam = eps.getEigenvalue(k)
        vecMap.P2N(eignModePETScReal, gfu.vecs[k])
        print(lam)

    return gfu, eps

def SLEPcKrylovSchur(a, b, fes):
    mat = n2p.CreatePETScMatrix(a.mat, fes.FreeDofs())
    mat_b = n2p.CreatePETScMatrix(b.mat, fes.FreeDofs())

    st = spc.ST().create()
    st.setType(spc.ST.Type.SINVERT)
    st.setShift(4000.)

    eps = spc.EPS()
    eps.create()

    eps.setType(spc.EPS.Type.KRYLOVSCHUR)
    eps.setTolerances(tol=1e-4, max_it=4000)
    eps.setProblemType(spc.EPS.ProblemType.GHEP)
    eps.setDimensions(20, spc.DECIDE, spc.DECIDE)

    eps.setST(st)
    eps.setOperators(mat, mat_b)

    eps.solve()

    nconv = eps.getConverged()
    gfu = ngsolve.GridFunction(fes, multidim=nconv)
    vecMap = n2p.VectorMapping(fes.ParallelDofs(), fes.FreeDofs())
    for k in range(nconv):
        eignModePETScReal = mat.createVecLeft()
        eignModePETScImag = mat.createVecLeft()
        eps.getEigenvector(k, eignModePETScReal, eignModePETScImag)
        lam = eps.getEigenvalue(k)
        vecMap.P2N(eignModePETScReal, gfu.vecs[k])
        print(lam)

    return gfu, eps

comm = MPI.COMM_WORLD
if comm.rank == 0:
    path      = sys.argv[1]
    mesh      = ngsolve.Mesh(path)
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


