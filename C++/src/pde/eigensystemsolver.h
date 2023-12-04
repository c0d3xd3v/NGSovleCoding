#ifndef EIGENSYSTEMSOLVER_H
#define EIGENSYSTEMSOLVER_H

#undef NETGEN_PYTHON
#include <comp.hpp>
using namespace ngcomp;

class EigenSystemSolver
{
private:
    Array<Complex> lam;
    Array<shared_ptr<BaseVector>> evecs;

public:
    EigenSystemSolver(double precision,
                      std::shared_ptr<T_BilinearFormSymmetric<Complex> > a,
                      std::shared_ptr<T_BilinearFormSymmetric<Complex> > m,
                      std::shared_ptr<VectorH1FESpace> fes,
                      double shift = 4000);

    Array<Complex> &getLams();

    Array<shared_ptr<BaseVector>> &getEvecs();
};


#endif // EIGENSYSTEMSOLVER_H
