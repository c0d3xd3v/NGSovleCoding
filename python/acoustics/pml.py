from ngsolve import *
from netgen.geom2d import SplineGeometry

# Geometry
geo = SplineGeometry()
geo.AddCircle((0.0, 0.0), 0.8,  bc="outer")
geo.AddRectangle((-0.05, -0.025), (0.025, 0.275),
                  leftdomain=0, rightdomain=1, bc="scat")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

fes = H1(mesh, order=5, complex=True)
u, v = fes.TnT()

# Wavenumber & source
omega = 15
pulse = 5e4*exp(-(40**2)*((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)))

# Forms
a = BilinearForm(fes)
a += grad(u)*grad(v)*dx - omega**2*u*v*dx
a += -omega*1j*u*v * ds("outer")
a.Assemble()

f = LinearForm(fes)
f += pulse * v * dx
f.Assemble();

count = 0
for v in f.vec:
    if v != 0:
        print(str(count) + " : " + str(v))
        count = count + 1

gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec
Draw(gfu)
