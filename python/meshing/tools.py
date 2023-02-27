from ngsolve import *
from netgen.stl import *
from netgen.csg import *
from netgen.meshing import *


def computeMeshCenter(mesh):
    centerx = 0
    centery = 0
    centerz = 0
    count   = 0

    for e in mesh.Elements2D():
        for v in e.vertices:
            centerx = centerx + mesh[v][0]
            centery = centery + mesh[v][1]
            centerz = centerz + mesh[v][2]
            count = count + 1

    centerx = centerx / count
    centery = centery / count
    centerz = centerz / count

    return [centerx, centery, centerz]


def meshBoundingRadius(mesh, center):
    R = 0
    for e in mesh.Elements2D():
        for v in e.vertices:
                x = -center[0] + mesh[v][0]
                y = -center[1] + mesh[v][1]
                z = -center[2] + mesh[v][2]
                tmp = sqrt(x*x + y*y + z*z)
                if tmp > R:
                    R = tmp
    return R


def addDomain(mesh, domainMesh, fd):
    pmap = {}
    for e in domainMesh.Elements2D():
        for v in e.vertices:
            if v not in pmap:
                pmap[v] = mesh.Add(domainMesh[v])
    # copy surface elements from first mesh to new mesh
    # we have to map point-numbers:
    for e in domainMesh.Elements2D():
        mesh.Add(Element2D(fd, [pmap[v] for v in e.vertices]))
    return mesh
