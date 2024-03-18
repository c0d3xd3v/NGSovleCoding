import ipyparallel as ipp
from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry
import netgen

import time

from Visualization.qt.VtkFEMMeshWidget import showTest
from PySide6 import QtCore, QtWidgets
from PySide6.QtGui import QIcon, QAction
from PySide6.QtWidgets import QToolBar, QPushButton

def mpi_example():

    import netgen.gui
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    import numpy as np
    from netgen.csg import unit_cube
    from ngsolve.krylovspace import CGSolver
    #import petsc4py.PETSc as psc
    import ngsolve
    import sys

    from numerics.preconditioning import fe_preconditioning
    from numerics.eigensystemsolver import solveEigensystem
    from elasticity.eigenfrequencies import build_elasticity_system_on_fes, steel


    if comm.rank == 0:
        ngmesh = unit_cube.GenerateMesh(maxh=0.1)
        ngmesh.Distribute(comm)
    else:
        ngmesh = netgen.meshing.Mesh.Receive(comm)
    mesh = ngsolve.Mesh(ngmesh)

    count           = 20
    solid_fes       = ngsolve.VectorH1(mesh, order=2, complex=True)
    a, b            = build_elasticity_system_on_fes(steel, solid_fes)
    a, b, pre, kapa = fe_preconditioning(solid_fes, a, b, "local")
    gfu, lams       = solveEigensystem(solid_fes, a, b, count, "arnoldi", pre)

    return  mesh.nv #f"Hello World from rank {comm.Get_rank()}. total ranks={comm.Get_size()}"

def startMPI():
    # request an MPI cluster with 4 engines
    with ipp.Cluster(engines='mpi', n=5) as rc:
        # get a broadcast_view on the cluster which is best
        # suited for MPI style computation
        view = rc.broadcast_view()
        # run the mpi_example function on all engines in parallel
        start = time.time()
        r = view.apply_sync(mpi_example)
        print(time.time() - start)
        # Retrieve and print the result from the engines
        print(r)
    # at this point, the cluster processes have been shutdown

app = QtWidgets.QApplication(sys.argv)

button = QPushButton("start mpi")
button.clicked.connect(startMPI)
button.show()

app.exec()
