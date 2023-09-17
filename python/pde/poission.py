import ngsolve
from meshing.tools import *
from elasticity.eigenfrequencies import *


def solvePoission(fes, gfu):
    # define trial- and test-functions
    u = fes.TrialFunction()
    v = fes.TestFunction()

    # the right hand side
    f = LinearForm(fes)
    f += 0. * v * dx
    # the bilinear-form 
    a = BilinearForm(fes, symmetric=True)
    a += grad(u)*grad(v)*dx

    a.Assemble()
    f.Assemble()

    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs()) * res

    return gfu
