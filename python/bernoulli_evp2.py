from ngsolve import *
from ngsolve.meshes import Make1DMesh
from ngsolve.webgui import Draw

mesh = Make1DMesh(100)

fesw = H1(mesh, order=2, dirichlet="left|right")
fessigma = H1(mesh, order=2, dirichlet=" ")

w,v = fesw.TnT()
sigma,tau = fessigma.TnT()

a = BilinearForm(sigma*tau*dx).Assemble()
b = BilinearForm(grad(w)*grad(tau)*dx).Assemble()

#k = BilinearForm(grad(w)*grad(v)*dx).Assemble()

S = b.mat.T @ a.mat.Inverse() @ b.mat
print (S.GetOperatorInfo())

pre = b.mat.Inverse(fesw.FreeDofs())
pre = pre@pre   # precond for 4th order problem
evals, evecs = solvers.LOBPCG(S, b.mat, pre=pre, num=5, maxit=100)

Draw(evals)



