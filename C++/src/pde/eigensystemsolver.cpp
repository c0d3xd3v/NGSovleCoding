#include "eigensystemsolver.h"

#include <comp.hpp>

EigenSystemSolver::EigenSystemSolver(
        double precision,
        std::shared_ptr<T_BilinearFormSymmetric<Complex>> a,
        std::shared_ptr<T_BilinearFormSymmetric<Complex>> m,
        std::shared_ptr<VectorH1FESpace> fes,
        double shift)
{
    LocalHeap lh(100000);

    Flags flags;
    //flags.SetFlag("type", "amgh1");

    //GetPreconditionerClasses().Print(std::cout);
    // multigrid, direct, local, bddc, bddcc, bddcrc, h1amg
    auto creator = GetPreconditionerClasses().GetPreconditioner("h1amg");
    shared_ptr<Preconditioner> pre = creator->creatorbf(a, flags, "h1amg");

    a->Assemble(lh);
    m->Assemble(lh);

    EigenSystem es(*(a->GetMatrixPtr()));
    es.SetPrecision(precision);
    es.SetPrecond(pre->GetMatrix());
    es.Calc();
    //es.PrintEigenValues(std::cout);

    int num = 15;
    lam = Array<Complex>(num);
    evecs = Array<shared_ptr<BaseVector>>(num);
    Arnoldi<Complex> arnoldi(a->GetMatrixPtr(), m->GetMatrixPtr(), fes->GetFreeDofs());
    arnoldi.SetShift(shift);
    arnoldi.Calc(3*num+1, lam, num, evecs);
    for(int i = 0; i < num; i++)
    {
        std::cout << lam[i] << std::endl;
    }
}

Array<Complex> &EigenSystemSolver::getLams()
{
    return lam;
}

Array<shared_ptr<BaseVector>> &EigenSystemSolver::getEvecs()
{
    return evecs;
}
