import ngsolve
from meshing.tools import *
from elasticity.eigenfrequencies import *


def solveHelmholtz(fes, gfu, roh=0.001): 
    u = fes.TrialFunction()
    v = fes.TestFunction()

    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    a += -roh**2 * u * v * dx
    a += -roh*1j * u * v * ds("outer")

    # bddc, h1amg, multigrid, local, direct
    precond = 'h1amg'
    pre = ngsolve.Preconditioner(a, precond)
    a.Assemble()

    lams = ngsolve.krylovspace.EigenValues_Preconditioner(mat=a.mat, pre=pre)
    print(lams)

    g = CoefficientFunction((0.))
    # the right hand side
    f = LinearForm(fes)
    f += g * v * dx
    f.Assemble()

    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec

    #gfu.vec.data += ngsolve.solvers.GMRes(a.mat, res, freedofs=fes.FreeDofs(), tol=1e-5, maxsteps=3000, restart=100, printrates=False)
    #solver = ngsolve.CGSolver(mat=a.mat, pre=pre, complex=True, conjugate=True, maxsteps=20000)
    #gfu.vec.data = solver * res
    #gfu.vec.data += ngsolve.solvers.CG(a.mat, res, freedofs=fes.FreeDofs(), tol=1e-9, maxsteps=300000, printrates=True)
    # inverse : sparsecholesky, mumps, umfpack
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="mumps") * res
    return gfu
