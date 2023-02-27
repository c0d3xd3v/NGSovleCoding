from ngsolve import *
from netgen.geom2d import SplineGeometry
from netgen.meshing import *
from elasticity.eigenfrequencies import *

# Geometry
circlegeo = SplineGeometry()
circlegeo.AddCircle((0.0, 0.0), 0.8,  bc="outer")
circlegeo.AddRectangle((-0.05, -0.025), (0.025, 0.275), leftdomain=0, rightdomain=1, bc="source")

rectgeo = SplineGeometry()
rectgeo.AddRectangle((-0.05, -0.025), (0.025, 0.275), bc="source")

circlemesh = circlegeo.GenerateMesh(maxh=0.025)
rectmesh = rectgeo.GenerateMesh(maxh=0.025)




num = 20
shift = 10.0
dirichlet = [1]
mesh = ngsolve.Mesh(rectmesh)
a, b, fes = build_elasticity_system(mesh, steel, dirichlet, shift)
u_rect = ngsolve.GridFunction(fes, multidim=num)
lams = ngsolve.ArnoldiSolver(a.mat, b.mat, fes.FreeDofs(), u_rect.vecs, shift)




cmesh = ngsolve.Mesh(circlemesh)
fes = VectorH1(cmesh, order=1, complex=True)
u, v = fes.TnT()
# Wavenumber
omega = 20
# Forms
a = BilinearForm(fes)
a += InnerProduct(grad(u), grad(v))*dx - omega**2*u*v*dx
a += -omega*1j*u*v*ds("outer")
a.Assemble()

F = u_rect.MDComponent(0)

f = LinearForm(fes)
f += -omega*v*F*dx
f.Assemble()

gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec

Draw(gfu)

