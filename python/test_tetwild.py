import os
import sys
import igl
import numpy as np

import pytetwild
import netgen.meshing as nm

from meshing.tools import *
from Visualization.vtkhelper import iglToVtkPolydata


path = sys.argv[1]

file_path, file = os.path.split(path)[0], os.path.split(path)[1]
file_name, file_ext = os.path.splitext(file)[0], os.path.splitext(file)[1]

wrapper = pytetwild.FTetWildWrapper(stop_energy=10, ideal_edge_length_rel=0.005, eps_rel=0.0005)
sv, sf = igl.read_triangle_mesh(path)
wrapper.loadMeshGeometry(sv, sf)
wrapper.tetrahedralize()
tris, tets, nods = wrapper.getSurfaceIndices()

wrapper.save(file_path+"/" + file_name)
surface_mesh_path = file_path + "/" + file_name + "_surface.obj"
igl.write_triangle_mesh(surface_mesh_path, nods, tris)
vtkPolydata = iglToVtkPolydata(tris, nods)

mesh = nm.Mesh()

pmap = {}
for i, p in enumerate(nods):
    mp = MeshPoint(Point3d(p[0], p[1], p[2]))
    pmap[i] = mesh.Add(mp)

fds = mesh.Add(FaceDescriptor(bc=1, domin=0, surfnr=0))
mesh = addSurfaceFromLists(fds, mesh, pmap, tris)
mesh = addVolumeFromLists(0, mesh, pmap, tets)
mesh.Update()
mesh.Save(file_path+"/" + file_name + ".vol")
