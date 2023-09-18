#ifndef ELASTICITYSYSTEM_H
#define ELASTICITYSYSTEM_H

#undef NETGEN_PYTHON
#include <comp.hpp>
using namespace ngcomp;

#include "material.h"

class ElasticityFESetup
{
private:
    double shift;

    std::shared_ptr<VectorH1FESpace> fes;
    std::shared_ptr<T_BilinearFormSymmetric<Complex>> bfa;
    std::shared_ptr<T_BilinearFormSymmetric<Complex>> bfm;

public:
    ElasticityFESetup(std::shared_ptr<MeshAccess> ma,
                  Material &material,
                  Array<double> &dirbnd,
                  double shift = 10.0);
    double getShift();
    std::shared_ptr<VectorH1FESpace> &getFes();
    std::shared_ptr<BaseMatrix> getBfaMatrix();
    std::shared_ptr<BaseMatrix> getBfmMatrix();
};
#endif // ELASTICITYSYSTEM_H
