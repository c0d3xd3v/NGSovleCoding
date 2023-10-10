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
    std::shared_ptr<MeshAccess> ma = make_shared<MeshAccess>(
                "test.vol");
/*
    Flags surf_flags;
    std::shared_ptr<HDivHighOrderSurfaceFESpace> surf_fes =
            std::make_shared<HDivHighOrderSurfaceFESpace>(ma, surf_flags);
*/
    ElasticityFESetup elasticityFes(ma, material, dirbnd, 10.0);

    EigenSystemSolver ess(precision,
                          elasticityFes.getBfaMatrix(),
                          elasticityFes.getBfmMatrix(),
                          elasticityFes.getFes(),
                          elasticityFes.getShift());

    Array<shared_ptr<BaseVector>> &evecs = ess.getEvecs();
    for(int j = 0; j < evecs.Size(); j++)
    {
        float tmp = std::abs(evecs[0]->FV<std::complex<double> >()[j]);
        std::cout << tmp << std::endl;
    }
    /*
    Array<float> mean;
    shared_ptr<BaseVector> mean_vec = CreateBaseVector(evecs[0]->Size(), true, 1);
    mean.SetSize(evecs[0]->Size());

    for(int i = 0; i < mean.Size(); i++)
    {
        mean[i] = 0.0f;
        for(int j = 0; j < evecs.Size(); j++)
        {
            mean[i] += std::abs(evecs[j]->FV<std::complex<double>>()[i]);
        }
        mean[i] /= evecs.Size();
    }

    std::vector<float> var;
    for(int i = 0; i < mean.Size(); i++)
    {
        float s = 0.0f;
        for(int j = 0; j < evecs.Size(); j++)
        {
             float tmp = std::abs(evecs[j]->FV<std::complex<double> >()[i]) - mean[i];
             tmp *= tmp;
             s += tmp;
        }
        s /= evecs.Size();
        var.push_back(std::sqrt(s));
        mean_vec->FV<std::complex<double>>()[i] = std::sqrt(s);
        std::cout << std::sqrt(s) << std::endl;
    }

    Array<shared_ptr<BaseVector>> arr_mean;
    arr_mean.Append(mean_vec);

    SolutionMapper sm(ma);
    //sm.mapSolutionToMesh(elasticityFes.getFes(), arr_mean);
    sm.mapSolutionToMesh(elasticityFes.getFes(), ess.getEvecs());

    //mean_vec->Print(std::cout);
    std::cout << " - > " << mean_vec->Size() << ", " << evecs[0]->Size() << std::endl;
    */
    return 0;
}