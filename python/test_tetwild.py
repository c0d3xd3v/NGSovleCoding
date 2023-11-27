import sys

import ngsolve
import netgen.meshing as nm

import wildmeshing as wm

from meshing.tools import *

path = sys.argv[1] # "/home/kai/Development/github/NGSovleCoding/data/tuningfork.stl"

tetra = wm.Tetrahedralizer(
#    coarsen=False,
    stop_quality=9,
#    max_its=20,
#    stage=2,
#    stop_p=-1,
    epsilon=0.00125,
    edge_length_r=1.0/25.0
    )

tetra.load_mesh(path)
tetra.set_log_level(16)
tetra.tetrahedralize()

tetmesh = tetra.get_tet_mesh()
surface = tetra.get_tracked_surfaces()
points = surface[0][0]
tris = surface[1][0]

mesh = nm.Mesh()
fds = mesh.Add(FaceDescriptor(bc=1, domin=0, surfnr=0))
mesh = addSurfaceFromLists(fds, mesh, points, tris)

points = tetmesh[0]
tets = tetmesh[1]
mesh = addVolumeFromLists(0, mesh, points, tets)

mesh.Update()

path, name = os.path.split(path)
name = os.path.basename(name).split('.')[0]
print(path + "/" + name + ".vol")
mesh.Save(path + "/" + name + ".vol")
print(path + "/" + name + ".msh")
tetra.save(path + "/" + name + ".msh",
#            floodfill=True,
#            manifold_surface=True,
#            use_input_for_wn=True,
#            correct_surface_orientation=True,
#            all_mesh=False
)
