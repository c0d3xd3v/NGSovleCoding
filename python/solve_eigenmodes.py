import sys
import ngsolve

from numerics.preconditioning import fe_preconditioning
from numerics.eigensystemsolver import solveEigensystem
from elasticity.eigenfrequencies import build_elasticity_system_on_fes, steel


path = "/home/kai/Deckel.vol"
count = 30

mesh = ngsolve.Mesh(path)
solid_fes = ngsolve.VectorH1(mesh, order=2, complex=True)
a, b = build_elasticity_system_on_fes(steel, solid_fes, (0+0.0j))
a, b = fe_preconditioning(a, b, 'h1amg')
gfu = solveEigensystem(solid_fes, a, b)

ngsolve.Draw(gfu)

vertices = []
triangles = []
tetraedras = []
eigenmodes = []

for v in mesh.vertices:
    p = [v.point[0], v.point[1], v.point[2]]
    vertices.append(p)

for e in mesh.Elements(ngsolve.BND):
    v = e.vertices
    tri = [v[0].nr, v[1].nr, v[2].nr]
    triangles.append(tri)

for k in range(count):
    E = gfu.MDComponent(k)
    eigenmode = []
    for i, v in enumerate(vertices):
        x = mesh(v[0], v[1], v[2])
        e = E(x)
        eigenmode.append(e)
    eigenmodes.append(eigenmode)

print(len(vertices))
print(len(triangles))
print(len(eigenmodes))


eigenmodes = []
names = []
for i in range(count):
    E = gfu.MDComponent(i)
    eigenmodes.append(E.real)
    eigenmodes.append(E.imag)
    names.append("mode-" + str(lams[i]) + "-real")
    names.append("mode-" + str(lams[i]) + "-imag")

# VTKOutput object
vtk =  ngsolve.VTKOutput(ma=mesh,
                coefs=eigenmodes,
                names = names,
                filename="mode02",
                subdivision=0,
                legacy=False)
vtk.Do()
