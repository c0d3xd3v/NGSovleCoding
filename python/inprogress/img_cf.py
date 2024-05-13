from ngsolve import *
import ngsolve
from ngsolve import internal

from videoio.ngsvideoread import *


internal.visoptions.showsurfacesolution = True
internal.visoptions.vecfunction = True
internal.viewoptions.drawedges = False
internal.viewoptions.drawoutline = False

ngsvideo = NGSVideoRead()

dt = 0.0001
u,v = ngsvideo.ImageSpace.TnT()

a = BilinearForm(ngsvideo.ImageSpace, symmetric=False)
a += 0.5*grad(u)*grad(v)*dx
a.Assemble()

m = BilinearForm(ngsvideo.ImageSpace, symmetric=False)
m += u * v * dx
m.Assemble()
mstar = m.mat.CreateMatrix()

mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
# corresponds to M* = M + dt * A
invmstar = mstar.Inverse(freedofs=ngsvideo.ImageSpace.FreeDofs())

for k in range(ngsvideo.frame_count-1):
    for i in range(10):
        res = - dt * a.mat * ngsvideo.gfu.vecs[k]
        ngsvideo.gfu.vecs[k].data += invmstar * res

dgfu = GridFunction(ngsvideo.ImageSpace, multidim=ngsvideo.frame_count)
for k in range(1, ngsvideo.frame_count-1):
    dgfu.vecs[k-1].data = ngsvideo.gfu.vecs[k] - ngsvideo.gfu.vecs[k-1]

I = BilinearForm(ngsvideo.FlowSpace, symmetric=False)
U,V = ngsvideo.FlowSpace.TnT()
I += U*dgfu.MDComponent(2)*V*dx
I.Assemble()
iinv = I.mat.Inverse(freedofs=ngsvideo.FlowSpace.FreeDofs())
Igfu = GridFunction(ngsvideo.FlowSpace, multidim=ngsvideo.frame_count)
Igfu.vec.data = iinv*Igfu.vec

Draw(ngsvideo.gfu, ngsvideo.mesh, "f")
Draw(dgfu, ngsvideo.mesh, "df")
Draw(Igfu, ngsvideo.mesh, "Igfu")
'''
gradgfu = GridFunction(FlowSpace, multidim=0)
for k in range(frame_count-1):
    for i in range(10):
        res = - dt * a.mat * gfu.vecs[k]
        gfu.vecs[k].data += invmstar * res

    print(k)
    gfuk = GridFunction(ImageSpace)
    gfuk.vec.data = gfu.vecs[k]
    grad_ = GridFunction(FlowSpace)
    grad_.Set(-grad(gfuk))
    gradgfu.AddMultiDimComponent(grad_.vec)


dgfu = GridFunction(ImageSpace, multidim=frame_count)
dgradgfu = GridFunction(FlowSpace, multidim=frame_count)
for k in range(1, frame_count-1):
    dgfu.vecs[k-1].data = gfu.vecs[k] - gfu.vecs[k-1]
    dgradgfu.vecs[k-1].data = gradgfu.vecs[k] - gradgfu.vecs[k-1]
    dg_grad = GridFunction(FlowSpace)
    dg_grad.Set(dgfu.MDComponent(k-1)*(gradgfu.MDComponent(k)))
    dgradgfu.vecs[k - 1].data += dg_grad.vec

Draw(gradgfu, mesh, "gradgfu")
Draw(gfu, mesh, "gfu")
Draw(dgfu, mesh, "dgfu")
Draw(dgradgfu, mesh, "dgradgfu")
'''