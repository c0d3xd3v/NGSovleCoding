import ngsolve
#import ngsolve.solvers as solvers

def solveEigensystem(fes, a, b):
    gfu = ngsolve.GridFunction(fes, multidim=count)
    lams = ngsolve.ArnoldiSolver(a.mat, b.mat, fes.FreeDofs(), list(gfu.vecs), 4000)
    #lams, evecs = ngsolve.solvers.PINVIT(a.mat, b.mat, pre, num=count, maxit=int(count/2), printrates=True, GramSchmidt=True)
    #for i in range(len(evecs)):
    #    gfu.vecs[i].data = evecs[i]
    #e, ev = solvers.LOBPCG(a.mat, b.mat, pre, num=count, maxit=300, printrates=True)
    #for i in range(count):
    #    gfu.vec.data[i] = ev[i][0]
    return gfu
