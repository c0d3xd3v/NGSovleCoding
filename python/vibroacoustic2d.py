from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry

from helmholtz import *

geo = SplineGeometry()
C = geo.AddCircle((0.0, 0.0), 2.8, 
                  leftdomain=1,
                  bc="outer")
side = 1.
geo.AddRectangle((-side*0.1, -side*0.0),
                 (side*0.1, side),
                 bcs=["fixed","solid","solid","solid"],
                 leftdomain=2, rightdomain=1)
geo.SetMaterial (1, "air")
geo.SetMaterial (2, "solid")
mesh = ngsolve.Mesh(geo.GenerateMesh(maxh=0.05))

solid_fes = VectorH1(mesh, definedon="solid", dirichlet="fixed", order=1, complex=True)
eigenmodes = solveElasticityEigenmodes(solid_fes)

fes = H1(mesh, definedon="air", dirichlet="fixed|solid", order=2, complex=True)
roh = 2.205
gfu = GridFunction(fes, name='gfu')
n = specialcf.normal(2)
E = eigenmodes.MDComponent(19)
g = roh*BoundaryFromVolumeCF(E)
gfu.Set(g*n, definedon=mesh.Boundaries("solid"))
gfu = solveHelmholtz(fes, gfu, roh)

print(mesh.GetMaterials())
print(mesh.GetBoundaries())

Draw(gfu)
Draw(eigenmodes)
