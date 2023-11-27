import ngsolve
import netgen.meshing as nm
import wildmeshing as wm

from meshing.tools import *

path = "/home/kai/Development/github/NGSovleCoding/data/tw/tuningfork.msh"
tetra = wm.Tetrahedralizer()

tetra.load_mesh(path)

tetmesh = tetra.get_tet_mesh()
surface = tetra.get_tracked_surfaces()
mesh = nm.Mesh()

points = surface[0][0]
tris = surface[1][0]
fds = mesh.Add(FaceDescriptor(bc=0, domin=0, surfnr=0))
mesh = addSurfaceFromLists(fds, mesh, points, tris)

points = tetmesh[0]
tets = tetmesh[1]
mesh = addVolumeFromLists(0, mesh, points, tets)

mesh.Update()
mesh.Save("/home/kai/Development/github/NGSovleCoding/data/tw/tw_tuningfork.vol")
#tetra.save("/home/kai/Development/github/NGSovleCoding/data/tw/tw_tuningfork")
