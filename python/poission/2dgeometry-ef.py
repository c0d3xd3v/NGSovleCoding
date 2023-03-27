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
circlegeo.AddCircle((0.0, 0.0), 0.8,  bc="outer")
circlegeo.AddRectangle((-0.05, -0.025), (0.025, 0.275), leftdomain=0, rightdomain=1, bc="source")

rectgeo = SplineGeometry()
rectgeo.AddRectangle((-0.05, -0.025), (0.025, 0.275))

circlemesh = circlegeo.GenerateMesh(maxh=0.075)
rectmesh = rectgeo.GenerateMesh(maxh=0.075)

# create an empty mesh
mesh = netgen.meshing.Mesh()

# a face-descriptor stores properties associated with a set of surface elements
# bc .. boundary condition marker,
# domin/domout .. domain-number in front/back of surface elements (0 = void),
# surfnr .. number of the surface described by the face-descriptor

fd_outside = mesh.Add(FaceDescriptor(bc=1,domin=1,surfnr=1))
fd_inside = mesh.Add(FaceDescriptor(bc=2,domin=2,domout=1,surfnr=2))

add2DMesh(circlemesh, mesh, fd_outside)
add2DMesh(rectmesh, mesh, fd_inside)

Draw(ngsolve.Mesh(mesh))

