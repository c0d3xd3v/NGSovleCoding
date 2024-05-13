from mpi4py import MPI
import ngsolve

import petsc4py.PETSc as psc
import slepc4py.SLEPc as spc
import ngsolve.ngs2petsc as n2p

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

def SLEPcGD(a, b, fes, count = 14):
    max_it = 100
    req_ep = count

    a.Assemble()
    b.Assemble()

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
    lams = []
    for k in range(nconv):
        eignModePETScReal = mat.createVecLeft()
        eignModePETScImag = mat.createVecLeft()
        eps.getEigenvector(k, eignModePETScReal, eignModePETScImag)
        lams.append(eps.getEigenvalue(k))
        vecMap.P2N(eignModePETScReal, gfu.vecs[k])
        print(lams[k])

    return gfu, lams

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
