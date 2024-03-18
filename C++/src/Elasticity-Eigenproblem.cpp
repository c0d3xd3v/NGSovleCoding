#include <limits>

#include "pde/elasticitysystem.h"
#include "pde/eigensystemsolver.h"
#include "pde/solutionmapper.h"

#include <vtk/QQuickVTKRenderItem.h>

#include <ngstd.hpp>
#include <nginterface_v2.hpp>
#include <meshaccess.hpp>

#include <mpi.h>

int main(int argc, char** argv)
{
    std::cout << argv[1] << ", " << steel.lam << ", " << steel.mu << ", " << steel.rho << std::endl;

    MPI_Init(&argc, &argv);

    netgen::NgMPI_Comm comm(MPI_COMM_WORLD);

    netgen::Ngx_Mesh mesh(argv[1], comm);
    std::cout << mesh.GetCommunicator().Size() << std::endl;
    std::shared_ptr<MeshAccess> ma = make_shared<MeshAccess>(mesh.GetMesh());
    shared_ptr<netgen::Mesh> ngmesh = ma->GetNetgenMesh();

    netgen::Mesh *_m = ngmesh.get();
    //_m->Distribute();
    //_m->ReceiveParallelMesh();

/*
    Flags surf_flags;
    std::shared_ptr<HDivHighOrderSurfaceFESpace> surf_fes =
            std::make_shared<HDivHighOrderSurfaceFESpace>(ma, surf_flags);
*/

    Material &material = steel;
    const double precision = 1.0e-10;
    Array<double> dirbnd;
    ElasticityFESetup elasticityFes(ma, material, dirbnd);

    EigenSystemSolver ess(precision,
                          elasticityFes.getBfaMatrix(),
                          elasticityFes.getBfmMatrix(),
                          elasticityFes.getFes(), 4000);


    //Array<shared_ptr<BaseVector>> &evecs = ess.getEvecs();

    /*
    std::cout << std::endl << "evecs.Size : " << evecs.Size() << std::endl;
    for(int j = 0; j < evecs.Size(); j++)
    {
        std::complex<double> tmp = evecs[0]->FV<std::complex<double>>()[j];
        std::cout << tmp << std::endl;
    }
    */

    //SolutionMapper sm(ma);
    //sm.saveVTK(elasticityFes.getFes(), ess.getEvecs(), ess.getLams());
    //MPI_Finalize();
    return 0;
}
