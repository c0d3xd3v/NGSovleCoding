import time
import numpy as np
import ngsolve
from Visualization.vtkhelper import iglToVtkPolydata, addScalarCellData
from Visualization.VtkFEMActor import VtkFEMActor

def gfuActor(mesh, gfu, count, dim=3):
    vertices = []
    triangles = []
    tetraedras = []
    eigenmodes = []
    import netgen

    #print(np.array(mesh.ngmesh.Points()))
    sf_ = np.array(mesh.ngmesh.Elements2D())
    print(np.array(sf_[:]).shape)
    #print(np.array(mesh.ngmesh.Elements3D()))

    for v in mesh.vertices:
        p = [v.point[0], v.point[1], v.point[2]]
        vertices.append(p)

    for e in mesh.Elements(ngsolve.BND):
        v = e.vertices
        tri = [v[0].nr, v[1].nr, v[2].nr]
        triangles.append(tri)

    if dim == 1:
        E = gfu
        eigenmode = []
        for i, v in enumerate(vertices):
            x = mesh(v[0], v[1], v[2])
            eigenmode.append(E.real(x))
        eigenmodes.append(eigenmode)
    else:
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
        polyData = addScalarCellData(polyData, eigenmodes[mid], dim, name)

    modeshapeActor = VtkFEMActor()
    modeshapeActor.setDataset(polyData)
    return modeshapeActor, polyData

def gfuActor2(mesh, gfu, count, dim=3):
    eigenmodes = [0]*count
    import netgen

    time_start = time.time()
    vertices = [ [p[0], p[1], p[2]] for p in mesh.ngmesh.Points() ]
    print("point copy : ", time.time() - time_start)

    time_start = time.time()
    triangles2 = [(t[0][0:3] - 1).tolist() for t in np.array(mesh.ngmesh.Elements2D())]
    print("triangle copy : ", time.time() - time_start)

    time_start = time.time()
    polyData = iglToVtkPolydata(triangles2, vertices)
    print("polydata copy : ", time.time() - time_start)

    time_start = time.time()
    meshpoints = [mesh(v[0], v[1], v[2]) for v in vertices]
    print("meshpoints generated : ", time.time() - time_start)

    if dim == 1:
        E = gfu
        eigenmode = [E.real(x) for x in meshpoints]
        eigenmodes[0] = eigenmode
    else:
        for k in range(count):
            E = gfu.MDComponent(k)
            name = "eigenmode" + str(k)
            time_start = time.time()
            eigenmodes[k] = [ E.real(x) for x in meshpoints ]
            print(name + " extract : ", time.time() - time_start)
            polyData = addScalarCellData(polyData, eigenmodes[k], dim, name)

    modeshapeActor = VtkFEMActor()
    modeshapeActor.setDataset(polyData)
    return modeshapeActor, polyData
