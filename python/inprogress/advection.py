from ngsolve import *
from netgen.geom2d import SplineGeometry

from netgen.occ import *
from pde.poission import solvePoission
'''
shape = Rectangle(2,2).Face().Move((-1,-1,0))
shape.edges.Min(X).name="left"
shape.edges.Max(X).name="right"
shape.edges.Min(Y).name="bottom"
shape.edges.Max(Y).name="top"
mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=0.05))
'''
periodic = SplineGeometry()
pnts = [ (-1,-1), (1,-1), (1,1), (-1,1) ]
pnums = [periodic.AppendPoint(*p) for p in pnts]

ltop = periodic.Append ( ["line", pnums[0], pnums[1]],bc="outer")
# This should be our master edge so we need to save its number.
lright = periodic.Append ( ["line", pnums[1], pnums[2]], bc="periodic")
periodic.Append ( ["line", pnums[2], pnums[3]], bc="outer")
# Minion boundaries must be defined in the same direction as master ones,
# this is why the the point numbers of this spline are defined in the reverse direction,
# leftdomain and rightdomain must therefore be switched as well!
# We use the master number as the copy argument to create a slave edge.
periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")
mesh = Mesh(periodic.GenerateMesh(maxh=0.05))

fes = H1(mesh, order=3)
u,v = fes.TnT()
time = 0.0
dt = 0.01

o0 = 0.025
o1 = 0.05
x0 = 0.5
y0 = 0.5
x1 = -0.5
y1 = -0.5
source0 = CF(exp(-0.5*((x/o0)**2 + (y/(5.*o0))**2)))
source1 = CF(exp(-0.5*(((x - x0)/o0)**2 + ((y - y0)/o0)**2)))
source2 = CF(exp(-0.5*(((x - x1)/o1)**2 + ((y - y1)/o1)**2)))

gfu = GridFunction(fes)
gfu.Set(source1 + source2)
scene = Draw(gfu,mesh,"u")

def TimeStepping(initial_cond = None, t0 = 0, tend = 1, nsamples = 30):
    if initial_cond:
        gfu.Set(initial_cond)
    cnt = 0; time = t0
    sample_int = int(floor(tend / dt / nsamples)+1)
    gfut = GridFunction(gfu.space,multidim=0)
    gfut.AddMultiDimComponent(gfu.vec)

    m = BilinearForm(fes, symmetric=False)
    m += u * v * dx
    m.Assemble()
    mstar = m.mat.CreateMatrix()

    fes_ = H1(mesh, order=3, dirichlet="outer")
    gfu_ = GridFunction(fes_)
    b = CoefficientFunction((2 * y * (1 - x * x), -2 * x * (1 - y * y)))
    while time < tend - 0.5 * dt:

        gfu_ = solvePoission(fes_, gfu_, gfu)
        a = BilinearForm(fes, symmetric=False)
        a += b * grad(u) * v * dx
        #a += 0.5*grad(u)*grad(v)*dx
        a.Assemble()
        mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
        # corresponds to M* = M + dt * A
        invmstar = mstar.Inverse(freedofs=fes.FreeDofs())

        res = - dt * a.mat * gfu.vec
        gfu.vec.data += invmstar * res

        print("\r",time,end="")
        if cnt % sample_int == 0:
            gfut.AddMultiDimComponent(gfu.vec)
        cnt += 1; time = cnt * dt
    return gfut

gfut = TimeStepping()
Draw(gfut, mesh, "gfut", interpolate_multidim=True, animate=True)
