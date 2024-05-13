from ngsolve import *
import ngsolve

import petsc4py.PETSc as psc
import ngsolve.ngs2petsc as n2p


def KrylovSolve(a, f, fes):
    psc_mat = n2p.CreatePETScMatrix(a, fes.FreeDofs())
    vecmap = n2p.VectorMapping(fes.ParallelDofs(), fes.FreeDofs())
    psc_f, psc_u = psc_mat.createVecs()
    vecmap.N2P(f, psc_f)

    ksp = psc.KSP(petsc_options={
    "ksp_type": "preonly",
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
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mkl_pardiso"
    })
    ksp.create()
    def myKSPMonitor(arg0, arg1, arg2):
        #global ksp_itrate
        #print(arg1, arg2, arg3, arg4)
        #ksp_itrate = arg1
        #if arg2 < 1e-2:
        print(arg1, arg2)
    #print(arg0.getConvergenceTest())
    ksp.setMonitor(myKSPMonitor)

    ksp.setOperators(psc_mat)
    ksp.setTolerances(rtol=1e-10, atol=0, divtol=1e16, max_it=500)
    ksp.solve(psc_f, psc_u)

    gfu = ngsolve.GridFunction(fes)
    vecmap.P2N(psc_u, gfu.vec)
    return gfu