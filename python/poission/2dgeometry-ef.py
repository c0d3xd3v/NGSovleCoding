from ngsolve import *
from netgen.geom2d import SplineGeometry
from netgen.meshing import *
from elasticity.eigenfrequencies import *

# Geometry
circlegeo = SplineGeometry()
circlegeo.AddCircle((0.0, 0.0), 0.8,  leftdomain=1)
circlegeo.AddRectangle((-0.05, -0.025), (0.025, 0.275), leftdomain=2, rightdomain=1)
circlegeo.SetDomainMaxH(2, 0.02) 

circlegeo.SetMaterial(2, "solid")
circlegeo.SetMaterial(1, "air")

circlemesh = circlegeo.GenerateMesh()
ngmesh = ngsolve.Mesh(circlemesh)

fes1 = H1(ngmesh, definedon="solid")
u1 = GridFunction(fes1, "u1")
u1.Set((x*x+y*y))

fes2 = H1(ngmesh, definedon="air")
u2 = GridFunction(fes2, "u2")
u2.Set(0)


Draw(ngmesh)
Draw(u1)
Draw(u2)

