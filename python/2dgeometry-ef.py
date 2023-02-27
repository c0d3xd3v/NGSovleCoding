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

mesh.Save("rectanglemesh.vol")

num = 20
shift = 10.0
dirichlet = [1]

mesh = ngsolve.Mesh(rectmesh)
a, b, fes = build_elasticity_system(mesh, steel, dirichlet, shift)
u_rect = ngsolve.GridFunction(fes, multidim=num)

lams = ngsolve.ArnoldiSolver(a.mat, b.mat, fes.FreeDofs(), u_rect.vecs, shift)

# ngsolve.Draw(u_rect)

# count = 0
#
# sourceVertexNumbers = []
# for el in circlemesh.Elements1D():
#     if circlemesh.GetBCName(el.edgenr - 1) == "source":
#         for v in el.vertices:
#             sourceVertexNumbers.append(v)
#
# sourceVertexNumbers = list(set(sourceVertexNumbers))
# print("sourceVertexNumbers : " + str(len(sourceVertexNumbers)))

# fp = u_rect.MDComponent(3)
#
# sum = 0
# for v in sourceVertexNumbers:
#     count = count + 1
#     mp = circlemesh[v]
#     fv = fp(ngsolve.Mesh(rectmesh)(mp[0], mp[1], mp[2]))
#     sum = sum + abs(sqrt(fv[0]*fv[0] + fv[1]*fv[1]))
#     print(sum)

# for el in circlemesh.Elements1D():
#      if circlemesh.GetBCName(el.edgenr - 1) == "source":
#          print("source vertices : ")
#          for v in el.vertices:
#              count = count + 1
#              mp = circlemesh[v]
#              fv = fp(ngsolve.Mesh(rectmesh)(mp[0], mp[1], mp[2]))
#              print(abs(fv[0]))
             # print(mp)

#print("count : " + str(count))

# for el in rectmesh.Elements1D():
#     for v in el.vertices:
#         mp = rectmesh[v]
#         print(mp)
#         print(fp(ngsolve.Mesh(rectmesh)(mp[0], mp[1], mp[2])))

cmesh = ngsolve.Mesh(circlemesh)
fes = VectorH1(cmesh, order=1, complex=True)
u, v = fes.TnT()

# Wavenumber
omega = 15
#
# Forms
a = BilinearForm(fes)
a += InnerProduct(grad(u), grad(v))*dx - omega**2*u*v*dx
a += -omega*1j*u*v*ds("outer")
a.Assemble()

F = u_rect.MDComponent(0)

f = LinearForm(fes)
f += -omega*v*F*dx
f.Assemble();

gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec

Draw(gfu)