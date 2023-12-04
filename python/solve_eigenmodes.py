import sys
import ngsolve
from ngsolve import *
from VibroAcoustic import *
from elasticity.eigenfrequencies import *
from acoustic.boundarysourcesolver import *
import ngsolve.solvers as solvers
from ngsolve.la import EigenValues_Preconditioner


#path = sys.argv[1] # '../data/tuningfork_.vol' #sys.argv[1] #"/home/kai/Development/github/NGSovleCoding/data/ridex-Body004.vol"
path = '../build-C++-Imported_Kit-Debug/test.vol'
#SetVisualization(clipping=True, clipnormal=tuple([0., 0., -1.]))

mesh = ngsolve.Mesh(path)

solid_fes = VectorH1(mesh, order=2, complex=True)
a, b = build_elasticity_system_on_fes(steel, solid_fes, (0+0.0j))

# bddc, h1amg, multigrid, local
precond = 'h1amg'
pre = Preconditioner(a, precond)
#pre = IdentityMatrix(solid_fes.ndof, complex=True)

a.Assemble()
b.Assemble()

#jac = a.mat.CreateBlockSmoother(solid_fes.CreateSmoothingBlocks())
#preJpoint = a.mat.CreateSmoother(solid_fes.FreeDofs())
lams = ngsolve.krylovspace.EigenValues_Preconditioner(mat=a.mat, pre=pre)
#kapa = max(lams)/min(lams)
#print(kapa)

count = 15
gfu = GridFunction(solid_fes, multidim=count)

lams = ngsolve.ArnoldiSolver(a.mat, b.mat, solid_fes.FreeDofs(), list(gfu.vecs), 4000)
#lams, evecs = ngsolve.solvers.PINVIT(a.mat, b.mat, pre, num=count, maxit=int(count/2), printrates=True, GramSchmidt=True)
#for i in range(len(evecs)):
#    gfu.vecs[i].data = evecs[i]
#e, ev = solvers.LOBPCG(a.mat, b.mat, pre, num=count, maxit=300, printrates=True)
#for i in range(count):
#    gfu.vec.data[i] = ev[i][0]

#print(gfu.vecs.data[0])

E = gfu.MDComponent(10)
# VTKOutput object
vtk = VTKOutput(ma=mesh,
                coefs=[E.real, E.imag],
                names = ["eigenmode0", "eigenmode1"],
                filename="mode0",
                subdivision=1,
                legacy=False)
vtk.Do()

Draw(gfu)

ngmesh = mesh.ngmesh
print(lams)


