import ngsolve
from ngsolve.la import EigenValues_Preconditioner

def fe_preconditioning(a, b, precond):
    # bddc, h1amg, multigrid, local
    # precond = 'h1amg'
    pre = None
    if precond == 'identity':
        pre = ngsolve.IdentityMatrix(solid_fes.ndof, complex=True)
    else:
        pre = ngsolve.Preconditioner(a, precond)

    a.Assemble()
    b.Assemble()

    #jac = a.mat.CreateBlockSmoother(solid_fes.CreateSmoothingBlocks())
    #preJpoint = a.mat.CreateSmoother(solid_fes.FreeDofs())
    lams = ngsolve.krylovspace.EigenValues_Preconditioner(mat=a.mat, pre=pre)
    l0 = min(lams)
    l1 = max(lams)
    kapa = -1.
    if l0 != 0.:
        kapa = l1/l0

    return a, b, pre, kapa
