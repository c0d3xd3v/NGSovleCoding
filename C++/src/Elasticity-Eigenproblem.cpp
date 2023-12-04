#include <limits>

#include "pde/elasticitysystem.h"
#include "pde/eigensystemsolver.h"
#include "pde/solutionmapper.h"

#include <vtk/QQuickVTKRenderItem.h>

int main(int argc, char** argv)
{
    Material &material = steel;
    const double precision = 1.0e-10;
    Array<double> dirbnd;
    std::shared_ptr<MeshAccess> ma = make_shared<MeshAccess>("test.vol");
/*
    Flags surf_flags;
    std::shared_ptr<HDivHighOrderSurfaceFESpace> surf_fes =
            std::make_shared<HDivHighOrderSurfaceFESpace>(ma, surf_flags);
*/
    ElasticityFESetup elasticityFes(ma, material, dirbnd);

    EigenSystemSolver ess(precision,
                          elasticityFes.getBfaMatrix(),
                          elasticityFes.getBfmMatrix(),
                          elasticityFes.getFes(), 4000);

    Array<shared_ptr<BaseVector>> &evecs = ess.getEvecs();

    /*
    std::cout << std::endl << "evecs.Size : " << evecs.Size() << std::endl;
    for(int j = 0; j < evecs.Size(); j++)
    {
        std::complex<double> tmp = evecs[0]->FV<std::complex<double>>()[j];
        std::cout << tmp << std::endl;
    }
    */

    SolutionMapper sm(ma);
    sm.saveVTK(elasticityFes.getFes(), ess.getEvecs(), ess.getLams());

    return 0;
}
