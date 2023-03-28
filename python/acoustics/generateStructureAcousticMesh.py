import ngsolve
from meshing.tools import *
from elasticity.eigenfrequencies import *

solidGeometry     = STLGeometry("../../test_data/tuningfork.stl")
solidMesh         = solidGeometry.GenerateMesh()
solidMeshCenter   = computeMeshCenter(solidMesh)
bndRadius         = meshBoundingRadius(solidMesh, solidMeshCenter)

ambientSphere     = CSGeometry()
ambientSphere.Add(Sphere(Pnt(solidMeshCenter), 2.5*bndRadius))
ambientSphereMesh = ambientSphere.GenerateMesh()

pmlSphere     = CSGeometry()
pmlSphere.Add(Sphere(Pnt(solidMeshCenter), 3.5*bndRadius))
pmlSphereMesh = pmlSphere.GenerateMesh()


# create an empty mesh
mesh = netgen.meshing.Mesh()

fd_solid   = mesh.Add(FaceDescriptor(bc=1,domin=3, domout=2,surfnr=1))
fd_outside = mesh.Add(FaceDescriptor(bc=2,domin=2, domout=1,surfnr=2))
fd_pml     = mesh.Add(FaceDescriptor(bc=3,domin=1, surfnr=3))

addDomain(mesh, solidMesh, fd_solid)
addDomain(mesh, ambientSphereMesh, fd_outside)
addDomain(mesh, pmlSphereMesh, fd_pml)

mesh.SetMaterial(3, "solid")
mesh.SetMaterial(2, "air")
mesh.SetMaterial(1, "pml")

mesh.GenerateVolumeMesh()
ngmesh = ngsolve.Mesh(mesh)

_fes: VectorH1 = VectorH1(ngmesh, definedon="solid", order=2, complex=True)
a, b = build_elasticity_system_on_fes(ngmesh, steel, _fes, 10)
u = ngsolve.GridFunction(_fes, multidim=20)
lams = ngsolve.ArnoldiSolver(a.mat, b.mat, _fes.FreeDofs(), list(u.vecs), 10)

fes2 = H1(ngmesh, definedon="air")
u2 = GridFunction(fes2, "u2")
u2.Set(0)

fes3 = H1(ngmesh, definedon="pml")
u3 = GridFunction(fes3, "u3")
u3.Set(1)


Draw(ngmesh)
Draw(u)
Draw(u2)
Draw(u3)