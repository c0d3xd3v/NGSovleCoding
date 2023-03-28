from ngsolve import *
import netgen.geom2d as geom2d

geo = geom2d.SplineGeometry()
p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(0,0), (1,0), (1,1), (0,1)] ]
p5,p6 =  [ geo.AppendPoint(x,y) for x,y in [(2,0), (2,1)] ]
geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0)
geo.Append (["line", p2, p3], leftdomain=1, rightdomain=2)
geo.Append (["line", p3, p4], leftdomain=1, rightdomain=0)
geo.Append (["line", p4, p1], leftdomain=1, rightdomain=0)
geo.Append (["line", p2, p5], leftdomain=2, rightdomain=0)
geo.Append (["line", p5, p6], leftdomain=2, rightdomain=0)
geo.Append (["line", p6, p3], leftdomain=2, rightdomain=0)

mesh = geo.GenerateMesh(maxh=0.05)

Draw(ngsolve.Mesh(mesh))