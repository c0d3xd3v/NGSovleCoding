import ngsolve
import sys
import os
from ngsolve import VectorH1


# A class container for (linear) material parameters
# ######################################################################################################################
class Material:
    E = 0
    nu = 0
    mu = 0
    lam = 0
    rho = 0
    type = "linear"

    def __init__(self, type="linear", E=0.0, nu=0.0, rho=0.0):
        self.E = E
        self.nu = nu
        self.type = type
        if abs(E) > 10e-15 and abs(nu) > 10e-15 and abs(nu + 1) > 10e-15 and abs(nu - 1 / 2) > 10e-15 and abs(rho) > 10e-15:
            self.mu = E / 2 / (1 + nu)
            self.lam = E * nu / ((1 + nu) * (1 - 2 * nu))
            self.rho = rho


# defining default materials:
# type=Material(type of model (e.g. "linear"), E (=young-module),nu (=poisson number), rho (=density))
aluminium = Material("linear", 69 * (10 ** 6), 0.33, 2700 * (10 ** (-9)))
steel = Material("linear", 220 * (10 ** 6), 0.28, 7700 * (10 ** (-9)))

# ######################################################################################################################
def build_elasticity_system_on_fes(_material, _fes):
    #_fes: VectorH1 = VectorH1(_mesh, order=2, dirichlet=_dirichlet, complex=True)
    _u, _v = _fes.TrialFunction(), _fes.TestFunction()
    #_a = ngsolve.BilinearForm(_fes, symmetric=True, eliminate_internal=True)
    _a = ngsolve.BilinearForm(_fes, symmetric=True)
    _a += ngsolve.SymbolicBFI(2 * _material.mu
                              * ngsolve.InnerProduct(1.0 / 2.0 * (ngsolve.grad(_u) + ngsolve.grad(_u).trans),
                                                     1.0 / 2.0 * (ngsolve.grad(_v) + ngsolve.grad(_v).trans))
                              + _material.lam * ngsolve.div(_u) * ngsolve.div(_v))
    #_b = ngsolve.BilinearForm(_fes, symmetric=True, eliminate_internal=True)
    _b = ngsolve.BilinearForm(_fes, symmetric=True)
    _b += ngsolve.SymbolicBFI(_material.rho * _u * _v)

    return _a, _b