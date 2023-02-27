import ngsolve
from meshing.tools import *


solidGeometry     = STLGeometry("../../test_data/tuningfork.stl")
solidMesh         = solidGeometry.GenerateMesh()
solidMeshCenter   = computeMeshCenter(solidMesh)
bndRadius         = meshBoundingRadius(solidMesh, solidMeshCenter)


ambientSphere     = CSGeometry()
ambientSphere.Add(Sphere(Pnt(solidMeshCenter), 2.5*bndRadius))
ambientSphereMesh = ambientSphere.GenerateMesh()


# create an empty mesh
mesh = netgen.meshing.Mesh()
mesh.SetMaterial(1, "solid")

fd_solid   = mesh.Add(FaceDescriptor(bc=1,domin=2,domout=1,surfnr=1))
fd_outside = mesh.Add(FaceDescriptor(bc=2,domin=1,surfnr=2))

addDomain(mesh, solidMesh, fd_solid)

Draw(ngsolve.Mesh(mesh))
