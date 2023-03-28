from ngsolve import *
from netgen.geom2d import SplineGeometry
from netgen.meshing import *
from elasticity.eigenfrequencies import *

def add2DMesh(input, container, fd_descriptor):
    pmap1 = {}
    for e in input.Elements2D():
        for v in e.vertices:
            if v not in pmap1:
                pmap1[v] = container.Add(input[v])
    # copy surface elements from first mesh to new mesh
    # we have to map point-numbers:
    for e in input.Elements2D():
        mesh.Add(Element2D(fd_descriptor, [pmap1[v] for v in e.vertices]))

def findSourceBoundaryElements(input):
    for el in input.Elements2D():
        print(el)

# Geometry
circlegeo = SplineGeometry()
circlegeo.AddCircle((0.0, 0.0), 0.8,  leftdomain=1)
circlegeo.AddRectangle((-0.05, -0.025), (0.025, 0.275), leftdomain=2, rightdomain=1)
circlegeo.SetDomainMaxH(2, 0.02) 

circlegeo.SetMaterial(2, "solid")
circlegeo.SetMaterial(1, "air")

rectgeo = SplineGeometry()
rectgeo.AddRectangle((-0.05, -0.025), (0.025, 0.275))


circlemesh = circlegeo.GenerateMesh()
#rectmesh = rectgeo.GenerateMesh(maxh=0.015)

# create an empty mesh
mesh = circlemesh #netgen.meshing.Mesh()

# a face-descriptor stores properties associated with a set of surface elements
# bc .. boundary condition marker,
# domin/domout .. domain-number in front/back of surface elements (0 = void),
# surfnr .. number of the surface described by the face-descriptor

#fd_inside = mesh.Add(FaceDescriptor(bc=1,  domin=1, domout=2, surfnr=1))
#fd_outside = mesh.Add(FaceDescriptor(bc=2, domin=2, domout=0, surfnr=2))

#add2DMesh(rectmesh, mesh, fd_inside)
#add2DMesh(circlemesh, mesh, fd_outside)

#mesh.SetMaterial(1, "solid")
#mesh.SetMaterial(2, "air")

ngmesh = ngsolve.Mesh(mesh)

fes1 = H1(ngmesh, definedon="solid")
u1 = GridFunction(fes1, "u1")
u1.Set((x*x+y*y))

fes2 = H1(ngmesh, definedon="air")
u2 = GridFunction(fes2, "u2")
u2.Set(0)


Draw(ngmesh)
Draw(u1)
Draw(u2)

