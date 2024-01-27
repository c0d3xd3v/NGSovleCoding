import ngsolve
from netgen.stl import *
from netgen.meshing import *


from VibroAcoustic import *
from elasticity.eigenfrequencies import *
from acoustic.boundarysourcesolver import *

#ngsolve.MPI_Init()

path = '../build-C++-Imported_Kit-Debug/test.vol'

#SetVisualization(clipping=True, clipnormal=tuple([0., 0., -1.]))

#solidMesh_ = ngsolve.Mesh(path)
#solidMesh = solidMesh_.ngmesh
#solidMesh.SetMaterial(1, "solid")

#mesh = generateVibroAcousticDomain_(solidMesh, maxh=20.0)
#mesh.Save("mesh-test.vol")

ngsmesh = ngsolve.Mesh("mesh-test.vol")
solid_fes = VectorH1(ngsmesh, definedon="solid", order=2, complex=True)
eigenmodes, lams = solveElasticityEigenmodes(solid_fes, 10, (0.+0.j), steel)
print(lams)
Draw(eigenmodes)

air_fes = H1(ngsmesh, definedon="air", dirichlet=ngsmesh.Boundaries("solid|fixed"), order=2, complex=True)
n = specialcf.normal(3)
E = eigenmodes.MDComponent(7)
g = BoundaryFromVolumeCF(E)
roh = rohForAir()
gfu = solveAcousticBoundaryValue(air_fes, n, g, roh, ngsmesh.Boundaries("solid|fixed"))

Draw(gfu)

# VTKOutput object
#vtk = VTKOutput(ma=solidMesh_,
#                coefs=[E.real, E.imag],
#                names = ["eigenmode0", "eigenmode1"],
#                filename="mode0",
#                subdivision=1,
#                legacy=False)
#vtk.Do()

# VTKOutput object
#vtk = VTKOutput(ma=ngmesh,
#                coefs=[gfu.real, gfu.imag],
#                names = ["pressurefield0Real", "pressurefield0Imag"],
#                filename="pressure0",
#                subdivision=1,
#                legacy=False)
#vtk.Do()
