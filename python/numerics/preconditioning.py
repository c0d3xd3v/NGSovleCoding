import ngsolve
from ngsolve.la import EigenValues_Preconditioner

def fe_preconditioning(fes, a, b, precond):
    # bddc, h1amg, multigrid, local, direct
    # precond = 'h1amg'
    pre = None
    if precond == 'identity':
        pre = ngsolve.IdentityMatrix(fes.ndof, complex=True)
    elif precond == 'gamg':
        pre = ngsolve.Preconditioner(a,"petsc", pctype="gamg", levels=10)
    else:
        pre = ngsolve.Preconditioner(a, precond)
        #pre = ngsolve.Preconditioner(a, type="direct", inverse = "masterinverse")

    a.Assemble()
    b.Assemble()

    #jac = a.mat.CreateBlockSmoother(solid_fes.CreateSmoothingBlocks())
    #preJpoint = a.mat.CreateSmoother(solid_fes.FreeDofs())
    lams = ngsolve.krylovspace.EigenValues_Preconditioner(mat=a.mat, pre=pre)
    #l0 = min(lams)
    #l1 = max(lams)
    kapa = -1.
    #if l0 != 0.:
    #    kapa = l1/l0

    return a, b, pre, kapa
