#include "eigensystemsolver.h"

EigenSystemSolver::EigenSystemSolver(
        double precision,
        std::shared_ptr<BaseMatrix> a,
        std::shared_ptr<BaseMatrix> m,
        std::shared_ptr<VectorH1FESpace> fes,
        double shift)
{
    EigenSystem es(*a);
    es.SetPrecision(precision);
    es.Calc();
    es.PrintEigenValues(std::cout);
    int num = 20;//es.NumEigenValues();
    std::cout << "NumEigenValues : " << num << std::endl;
    lam = Array<Complex>(num);
    evecs = Array<shared_ptr<BaseVector>>(num);
    Arnoldi<Complex> arnoldi(a, m, fes->GetFreeDofs());
    arnoldi.SetShift(shift);
    arnoldi.Calc(2*num+1, lam, num, evecs, nullptr);
}

Array<Complex> &EigenSystemSolver::getLams()
{
    return lam;
}

Array<shared_ptr<BaseVector>> &EigenSystemSolver::getEvecs()
{
    return evecs;
}
