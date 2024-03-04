import sys
import ngsolve
from meshing.tools import *
from elasticity.eigenfrequencies import *


def generateVibroAcousticDomain(solidMesh, maxh):
    solidMeshCenter   = computeMeshCenter(solidMesh)
    bndRadius         = meshBoundingRadius(solidMesh, solidMeshCenter)

    ambientSphere     = CSGeometry()
    ambientSphere.Add(Sphere(Pnt(solidMeshCenter), 2.8*bndRadius))
    ambientSphereMesh = ambientSphere.GenerateMesh(maxh=maxh)

    # create an empty mesh
    mesh = netgen.meshing.Mesh()

    fd_fixed   = mesh.Add(FaceDescriptor(bc=2,domin=2, domout=1,surfnr=1))
    fd_solid   = mesh.Add(FaceDescriptor(bc=1,domin=2, domout=1,surfnr=1))
    fd_outside = mesh.Add(FaceDescriptor(bc=3,domin=1,surfnr=2))

    addDomainSave(mesh, solidMesh, fd_solid, fd_fixed)
    addDomain(mesh, ambientSphereMesh, fd_outside)

    mesh.SetMaterial(2, "solid")
    mesh.SetMaterial(1, "air")

    mesh.SetBCName(2, "outer")
    mesh.SetBCName(1, "fixed")
    mesh.SetBCName(0, "solid")

    mesh.GenerateVolumeMesh(maxh=maxh)
    ngmesh = ngsolve.Mesh(mesh)

    print(ngmesh.GetMaterials())
    print(ngmesh.GetBoundaries())
    print(mesh.GetNDomains())

    return mesh

if __name__ == "__main__":
    path = sys.argv[1]
    solidMesh = ngsolve.Mesh(path).ngmesh
    solidMesh.SetMaterial(1, "solid")
    mesh = generateVibroAcousticDomain(solidMesh, maxh=30.0)
    mesh.Save("mesh-test.vol")
