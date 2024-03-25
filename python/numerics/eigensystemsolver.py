from math import sqrt

import ngsolve

from ngsolve.la import InnerProduct, MultiVector
from ngsolve import Projector, Norm, Matrix, Vector, IdentityMatrix

def LOBPCG(mata, matm, pre, num=1, maxit=20, initial=None, printrates=True):
    """Knyazev's cg-like extension of PINVIT"""
    import scipy.linalg

    r = mata.CreateRowVector()

    if initial:
        num=len(initial)
        uvecs = initial
    else:
        uvecs = MultiVector(r, num)

    vecs = MultiVector(r, 3*num)
    for v in vecs:
        r.SetRandom()
        v.data = pre * r

    if initial:
         vecs[0:num] = uvecs

    lams = Vector(num * [1])
    lams0 = Vector(num * [1])
    err0 = 0.
    errsum = 0.
    for i in range(maxit):
        uvecs.data = mata * vecs[0:num] - (matm * vecs[0:num]).Scale (lams)
        vecs[2*num:3*num] = pre * uvecs

        vecs.Orthogonalize(matm)

        asmall = InnerProduct (vecs, mata * vecs)
        msmall = InnerProduct (vecs, matm * vecs)

        ev,evec = scipy.linalg.eigh(a=asmall, b=msmall)
        lams = Vector(ev[0:num])
        err = Norm(lams0 - lams)

        checksum = 1
        if printrates and err0 != 0. and i > 1:
            es0 = errsum/(i-1)
            errsum += (err0- err)/err
            checksum = (es0 - errsum/i)
            print (i, "%2E" % abs(checksum), flush=True)

        uvecs[:] = vecs * Matrix(evec[:,0:num])
        vecs[num:2*num] = vecs[0:num]
        vecs[0:num] = uvecs

        if abs(checksum) < 1.e-7:
            break

        lams0 = Vector(lams)
        err0 = err

    return lams, uvecs


def solveEigensystem(fes, a, b, count, solver_name, pre):
    gfu = ngsolve.GridFunction(fes, multidim=count)
    lams = []
    if solver_name == "arnoldi":
        # inverse : 'sparsecholesky', 'pardiso', 'pardisospd', 'mumps', 'masterinverse', 'umfpack'
        lams = ngsolve.ArnoldiSolver(a.mat, b.mat, fes.FreeDofs(), list(gfu.vecs), 4000, inverse="pardiso")
    elif solver_name == "pinvit":
        lams, evecs = ngsolve.solvers.PINVIT(a.mat, b.mat, pre, num=count, maxit=300, printrates=False, GramSchmidt=True)
        for i in range(len(evecs)):
            gfu.vecs[i].data = evecs[i]
    elif solver_name == "lobpcg":
        lams, ev = LOBPCG(a.mat, b.mat, pre, num=count, maxit=3000, printrates=True)
        for i in range(count):
            gfu.vecs[i].data = ev[i]
    return gfu, lams

