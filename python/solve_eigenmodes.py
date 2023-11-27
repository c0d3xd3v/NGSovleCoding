import ngsolve

from VibroAcoustic import *
from elasticity.eigenfrequencies import *
from acoustic.boundarysourcesolver import *

path = "/home/kai/Development/github/NGSovleCoding/data/ridex-Body004.vol"

#SetVisualization(clipping=True, clipnormal=tuple([0., 0., -1.]))

mesh = ngsolve.Mesh(path)

solid_fes = VectorH1(mesh, order=2, complex=True)
eigenmodes, lams = solveElasticityEigenmodes(solid_fes, 10, (1+10000j), steel)
#print(len(lams))
Draw(eigenmodes)
