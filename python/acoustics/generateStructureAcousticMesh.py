import ngsolve
from meshing.tools import *


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

Draw(ngsolve.Mesh(mesh))
