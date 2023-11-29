import sys
import numpy as np
import ngsolve
import netgen.meshing as nm

import wildmeshing as wm

from meshing.tools import *

path = sys.argv[1] # "/home/kai/Development/github/NGSovleCoding/data/tuningfork.stl"

tetra = wm.Tetrahedralizer(
#    coarsen=False,
#    stop_quality=10,
#    max_its=20,
#    stage=2,
#    stop_p=-1,
#    epsilon=0.0001,
#    edge_length_r=1.0/20.0
    )

tetra.load_mesh(path)
tetra.set_log_level(16)
tetra.tetrahedralize()

tetmesh = tetra.get_tet_mesh()
surface = tetra.get_tracked_surfaces()

surf = tetra.get_surface_triangles()

points = tetmesh[0]
tets = tetmesh[1]
tris = surface[1][0]

mesh = nm.Mesh()
#mesh.SetMaterial(1, "default")
fds = mesh.Add(FaceDescriptor(bc=0, domin=0, surfnr=0))
#mesh = addSurfaceFromLists(fds, mesh, points, tris)
#mesh = addVolumeFromLists(0, mesh, points, tets)

pmap = {}
for i, p in enumerate(points):
    mp = MeshPoint(Point3d(p[0], p[1], p[2]))
    pmap[i] = mesh.Add(mp)

for i, tri in enumerate(tris):
    T = Element2D(fds, [pmap[v] for v in tri])
    mesh.Add(T)

points = tetmesh[0]
index = 0
for i, tet in enumerate(tets):
    vindices = [pmap[v] for v in tet]
    tmp = vindices[2]
    vindices[2] = vindices[3]
    vindices[3] = tmp
    T = Element3D(index, vindices)
    mesh.Add(T)

mesh.Update()

path, name = os.path.split(path)
name = os.path.basename(name).split('.')[0]
print(path + "/" + name + ".vol")
mesh.Save(path + "/" + name + "_.vol")
print(path + "/" + name + ".msh")
tetra.save(path + "/" + name + ".msh",
            floodfill=True,
            manifold_surface=True,
#            use_input_for_wn=True,
            correct_surface_orientation=True,
#            all_mesh=False
)



tris = tetra.get_surface_triangles()
