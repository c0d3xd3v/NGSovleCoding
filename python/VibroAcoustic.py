import ngsolve
from meshing.tools import *
from elasticity.eigenfrequencies import *

def generateVibroAcousticDomain(path):
    solidGeometry     = STLGeometry(path)
    solidMesh         = solidGeometry.GenerateMesh()
    solidMeshCenter   = computeMeshCenter(solidMesh)
    bndRadius         = meshBoundingRadius(solidMesh, solidMeshCenter)

    ambientSphere     = CSGeometry()
    ambientSphere.Add(Sphere(Pnt(solidMeshCenter), 2.5*bndRadius))
    ambientSphereMesh = ambientSphere.GenerateMesh()

    pmlSphere     = CSGeometry()
    pmlSphere.Add(Sphere(Pnt(solidMeshCenter), 4.5*bndRadius))
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

    mesh.SetBCName(2, "PmlOutInterface")
    mesh.SetBCName(1, "AirPmlInterface")
    mesh.SetBCName(0, "SolidAirInterface")

    #print(ngmesh.GetMaterials())
    #print(ngmesh.GetBoundaries())

    mesh.GenerateVolumeMesh(maxh=100.0)
    ngmesh = ngsolve.Mesh(mesh)

    #print(solidMesh.FaceDescriptors())
    #print(mesh.FaceDescriptors())
    return ngmesh

def addDomainSave(mesh, domainMesh, fd, fd_fixed):
    pmap = {}
    for e in domainMesh.Elements2D():
        for v in e.vertices:
            if v not in pmap:
                pmap[v] = mesh.Add(domainMesh[v])
    # copy surface elements from first mesh to new mesh
    # we have to map point-numbers:
    for e in domainMesh.Elements2D():
        #fd = domainMesh.FaceDescriptors()[e.index]
        if(e.index == 1):
            mesh.Add(Element2D(fd_fixed, [pmap[v] for v in e.vertices]))
        else:
            mesh.Add(Element2D(fd, [pmap[v] for v in e.vertices]))
    return mesh

def generateVibroAcousticDomain_(solidMesh, maxh):
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

    #print(ngmesh.GetMaterials())
    #print(ngmesh.GetBoundaries())
    #print(mesh.GetNDomains())
    #print(mesh)

    return ngmesh
