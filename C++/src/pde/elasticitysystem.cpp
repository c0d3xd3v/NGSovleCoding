#include "elasticitysystem.h"

ElasticityFESetup::ElasticityFESetup(
        std::shared_ptr<MeshAccess> ma,
        Material &material,
        Array<double> &dirbnd,
        double shift)
    : shift(shift)
{
    Flags flags_fes;
    flags_fes.SetFlag("order", 2);
    flags_fes.SetFlag("complex", true);
    //flags_fes.SetFlag("dirichlet", dirbnd);

    fes = make_shared<VectorH1FESpace>(ma, flags_fes);
    fes->Update();
    fes->FinalizeUpdate();

    const ProxyNode & u = fes->GetTrialFunction();
    const ProxyNode & v = fes->GetTestFunction();
    shared_ptr<CoefficientFunction> divu = u->Operator("div");
    shared_ptr<CoefficientFunction> divv = v->Operator("div");

    Flags flags_bfa;
    bfa = make_shared<T_BilinearFormSymmetric<Complex>>(fes, "a", fes->GetFlags());
    bfa->AddIntegrator(make_shared<SymbolicBilinearFormIntegrator>(
                           2.0 * material.mu *
                           InnerProduct(0.5 * (u->Deriv() + TransposeCF(u->Deriv())),
                                        0.5 * (v->Deriv() + TransposeCF(v->Deriv())))
                           + material.lam * InnerProduct(divu, divv)
                           + shift * material.rho * InnerProduct(u, v)
                           , VOL, VOL));

    bfm = make_shared<T_BilinearFormSymmetric<Complex>> (fes, "m", fes->GetFlags());
    bfm->AddIntegrator(make_shared<SymbolicBilinearFormIntegrator>(
                           material.rho * InnerProduct(u, v), VOL, VOL));

    LocalHeap lh(10000000);

    bfa->Assemble(lh);
    bfm->Assemble(lh);
}

double ElasticityFESetup::getShift()
{
    return shift;
}

std::shared_ptr<VectorH1FESpace> &ElasticityFESetup::getFes()
{
    return fes;
}

std::shared_ptr<BaseMatrix> ElasticityFESetup::getBfaMatrix()
{
    return bfa->GetMatrixPtr();
}

std::shared_ptr<BaseMatrix> ElasticityFESetup::getBfmMatrix()
{
    return bfm->GetMatrixPtr();
}
