import sys
import ngsolve

from numerics.preconditioning import fe_preconditioning
from numerics.eigensystemsolver import solveEigensystem
from elasticity.eigenfrequencies import build_elasticity_system_on_fes, steel

from Visualization.VtkNGSolve import gfuActor
from Visualization.qt.drawutils import Draw

path      = sys.argv[1]
mesh      = ngsolve.Mesh(path)
count     = 20

solid_fes       = ngsolve.VectorH1(mesh, order=2, complex=True)
a, b            = build_elasticity_system_on_fes(steel, solid_fes)
a, b, pre, kapa = fe_preconditioning(a, b, "h1amg")
gfu, lams       = solveEigensystem(solid_fes, a, b, count, "arnoldi", pre)

modeshapeActor  = gfuActor(mesh, gfu, count)
modeshapeActor.select_function("eigenmode12")
modeshapeActor.disableEdges()
Draw(modeshapeActor, periodic_timer=True)
