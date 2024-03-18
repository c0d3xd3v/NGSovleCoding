import ngsolve

def solveEigensystem(fes, a, b, count, solver_name, pre):
    gfu = ngsolve.GridFunction(fes, multidim=count)
    lams = []
    if solver_name == "arnoldi":
        # inverse : 'sparsecholesky', 'pardiso', 'pardisospd', 'mumps', 'masterinverse', 'umfpack'
        lams = ngsolve.ArnoldiSolver(a.mat, b.mat, fes.FreeDofs(), list(gfu.vecs), 4000, inverse="mumps")
    elif solver_name == "pinvit":
        lams, evecs = ngsolve.solvers.PINVIT(a.mat, b.mat, pre, num=count, maxit=300, printrates=False, GramSchmidt=True)
        for i in range(len(evecs)):
            gfu.vecs[i].data = evecs[i]
    elif solver_name == "lobpcg":
        lams, ev = ngsolve.solvers.LOBPCG(a.mat, b.mat, pre, num=count, maxit=300, printrates=False)
        for i in range(count):
            gfu.vecs[i].data = ev[i]
    return gfu, lams
