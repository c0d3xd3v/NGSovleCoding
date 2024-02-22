import ngsolve

def solveEigensystem(fes, a, b, count, solver_name, pre):
    gfu = ngsolve.GridFunction(fes, multidim=count)
    lams = []
    if solver_name == "arnoldi":
        lams = ngsolve.ArnoldiSolver(a.mat, b.mat, fes.FreeDofs(), list(gfu.vecs), 4000)
    elif solver_name == "pinvit":
        lams, evecs = ngsolve.solvers.PINVIT(a.mat, b.mat, pre, num=count, maxit=300, printrates=True, GramSchmidt=True)
        for i in range(len(evecs)):
            gfu.vecs[i].data = evecs[i]
    elif solver_name == "lobpcg":
        lams, ev = ngsolve.solvers.LOBPCG(a.mat, b.mat, pre, num=count, maxit=300, printrates=True)
        for i in range(count):
            gfu.vec.data[i] = ev[i][0]
    return gfu, lams
