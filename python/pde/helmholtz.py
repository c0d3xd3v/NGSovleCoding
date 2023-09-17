import ngsolve
from meshing.tools import *
from elasticity.eigenfrequencies import *


def solveHelmholtz(fes, gfu, roh=0.001): 
    u = fes.TrialFunction()
    v = fes.TestFunction()

    a = BilinearForm(fes)
    a += grad(u)*grad(v)*dx
    a += -roh**2*u*v*dx
    a += -roh*1j*u*v * ds("outer")
    a.Assemble()

    g = CoefficientFunction((0.))
     # the right hand side
    f = LinearForm(fes)
    f += g * v * dx

    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs()) * res

    #gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    
    return gfu