import numpy as np
import ngsolve
from Visualization.vtkhelper import iglToVtkPolydata, addScalarCellData
from Visualization.VtkModeshapeActor import VtkModeShapeActor

def gfuActor(mesh, gfu, count):
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
            eigenmode.append(E.real(x))
        eigenmodes.append(eigenmode)

    polyData = iglToVtkPolydata(triangles, vertices)
    for mid in range(count):
        emin = np.min(eigenmodes[mid])
        emax = np.max(eigenmodes[mid])
        for k in range(len(eigenmodes[mid])):
            eigenmodes[mid][k] = eigenmodes[mid][k]/(emax - emin)
        name = "eigenmode" + str(mid)
        polyData = addScalarCellData(polyData, eigenmodes[mid], 3, name)

    modeshapeActor = VtkModeShapeActor()
    modeshapeActor.setDataset(polyData)
    return modeshapeActor
