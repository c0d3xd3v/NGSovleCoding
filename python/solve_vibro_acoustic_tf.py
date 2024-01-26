import ngsolve

from VibroAcoustic import * 
from elasticity.eigenfrequencies import *
from acoustic.boundarysourcesolver import *


#ngsolve.MPI_Init()

path = '../build-C++-Imported_Kit-Debug/test.vol' # "/home/kai/Development/github/NGSovleCoding/data/tuningfork.stl"

SetVisualization(clipping=True, clipnormal=tuple([0., 0., -1.]))

solidMesh_ = ngsolve.Mesh(path)
solidMesh = solidMesh_.ngmesh
solidMesh.SetMaterial(1, "solid")

#solid_fes = VectorH1(solidMesh_, definedon="solid", order=2, complex=True)
#eigenmodes, lams = solveElasticityEigenmodes(solid_fes, 10, (0.+0.j), steel)
#print(len(lams))
Draw(solidMesh_)

ngmesh = generateVibroAcousticDomain_(solidMesh, maxh=20.0)
Draw(ngmesh)

#air_fes = H1(ngmesh, definedon="air", dirichlet=ngmesh.Boundaries("solid|fixed"), order=2, complex=True)
#n = specialcf.normal(3)
#E = eigenmodes.MDComponent(7)
#g = BoundaryFromVolumeCF(E)
#roh = rohForAir()
#gfu = solveAcousticBoundaryValue(air_fes, n, g, roh, ngmesh.Boundaries("solid|fixed"))

#Draw(gfu)

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

