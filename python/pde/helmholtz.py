import ngsolve
from meshing.tools import *
from elasticity.eigenfrequencies import *


def solveHelmholtz(fes, gfu, roh=0.001): 
    u = fes.TrialFunction()
    v = fes.TestFunction()

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a += -roh**2 * u * v * dx
    a += -roh*1j * u * v * ds("outer")

    # bddc, h1amg, multigrid, local
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

    #gfu.vec.data += ngsolve.solvers.GMRes(a.mat, res, freedofs=fes.FreeDofs(), tol=1e-6, maxsteps=3000, restart=100, printrates=True)
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res
    return gfu
