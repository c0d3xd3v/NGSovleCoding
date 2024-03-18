import ngsolve
from netgen.stl import *
from netgen.meshing import *


from meshing.AcousticMesh import *
from elasticity.eigenfrequencies import *
from acoustic.boundarysourcesolver import *
from numerics.preconditioning import fe_preconditioning
from numerics.eigensystemsolver import solveEigensystem

ngsolve.MPI_Init()

path = 'mesh-test.vol'

SetVisualization(clipping=True, clipnormal=tuple([0., 0., -1.]))

count = 10
ngsmesh = ngsolve.Mesh("mesh-test.vol")
solid_fes = VectorH1(ngsmesh, definedon="solid", order=2, complex=True)
a, b = build_elasticity_system_on_fes(steel, solid_fes)
a, b, pre, kapa = fe_preconditioning(solid_fes, a, b, 'h1amg')
eigenmodes, lams = solveEigensystem(solid_fes, a, b, count, "arnoldi", pre)

print(lams)
Draw(eigenmodes)

air_fes = H1(ngsmesh, definedon="air", dirichlet=ngsmesh.Boundaries("solid|fixed"), order=2, complex=True)
n = specialcf.normal(3)
E = eigenmodes.MDComponent(6)
g = BoundaryFromVolumeCF(E)
roh = rohForAir(lams[6])
print("roh : ", roh)
gfu = solveAcousticBoundaryValue(air_fes, n, g, roh, ngsmesh.Boundaries("solid|fixed"))

Draw(gfu)

