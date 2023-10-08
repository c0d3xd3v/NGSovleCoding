#include "solutionmapper.h"

SolutionMapper::SolutionMapper(std::shared_ptr<MeshAccess> ma)
    : ma(ma)
{
    ElementRange bndelements = ma->Elements(BND);
    for(int i = 0; i < bndelements.Size(); i++)
    {
        Ngs_Element el = ma->GetElement(bndelements[i]);
        for(int vid : el.Vertices())
            surfaceTriangles.push_back(vid);
    }

    std::set<unsigned int> tmp = std::set<unsigned int>(surfaceTriangles.begin(), surfaceTriangles.end());
    surfaceVertices = std::vector<unsigned int>(tmp.begin(), tmp.end());
    fv = std::vector<std::complex<double>[3]>(surfaceVertices.size());
    points = std::vector<float>(surfaceVertices.size()*3);
    fvr = std::vector<float>(surfaceVertices.size()*3);
    f_mag = std::vector<float>(surfaceVertices.size());

    for(int c = 0; c < surfaceVertices.size(); c++)
    {
        Vec<3> mp = ma->GetPoint<3>(surfaceVertices[c]);
        points[c*3 + 0] = mp[0];
        points[c*3 + 1] = mp[1];
        points[c*3 + 2] = mp[2];
    }
}

void SolutionMapper::mapSolutionToMesh(std::shared_ptr<VectorH1FESpace> fes,
                                       Array<shared_ptr<BaseVector> > &evecs)
{
    int num = evecs.Size();
    Flags gfu_flags;
    gfu_flags.SetFlag("multidim", num);
    std::shared_ptr<GridFunction> gfu = CreateGridFunction(fes, "gfu", gfu_flags);
    gfu->Update();
    for(int i = 0; i < num; i++)
        gfu->GetVector(i).Set(1.0, (*evecs[i]));
    gfu->GetMeshAccess()->SelectMesh();

    LocalHeapMem<100000> _lh("viscf::GetValue");
    for(int c = 0; c < surfaceVertices.size(); c++)
    {
        HeapReset hr(_lh);
        std::shared_ptr<ngfem::DifferentialOperator> evaluator = fes->GetEvaluator(VOL);
        IntegrationPoint ip;
        int elnr = -1;

        Array<int> elems;
        ma->GetVertexElements(surfaceVertices[c], elems);
        elnr = ma->FindElementOfPoint(ma->GetPoint<3>(surfaceVertices[c]), ip, true);

        ElementId ei(VOL, elnr);
        const FiniteElement & fel = fes->GetFE(ei, _lh);

        Array<int> dnums(fel.GetNDof(), _lh);
        fes->GetDofNrs(ei, dnums);
        auto & trafo = fes->GetMeshAccess()->GetTrafo(ei, _lh);

        Vector<Complex> elvec(fel.GetNDof()*fes->GetDimension());
        Vector<Complex> values(evaluator->Dim());

        gfu->GetVector().Set(1.0, (*evecs[4]));
        gfu->GetElementVector(dnums, elvec);

        evaluator->Apply(fel, trafo(ip, _lh), elvec, values, _lh);

        //std::cout << surfaceVertices[c] << " : " << values << std::endl;
        std::complex<double> tmp = 0.0;
        for(int i = 0; i < evaluator->Dim(); i++)
            tmp += values[i]*values[i];
        f_mag[c] = tmp.real() + tmp.imag();
    }

    float f_mag_max = f_mag[0];
    float f_mag_min = f_mag[0];

    for(int i = 0; i < f_mag.size(); i++)
        if(f_mag[i] < f_mag_min) f_mag_min = f_mag[i];
    if(f_mag_min < 0.0)
        for(int i = 0; i < f_mag.size(); i++)
            f_mag[i] = (f_mag[i] - f_mag_min);

    for(int i = 0; i < f_mag.size(); i++)
        if(f_mag[i] > f_mag_max) f_mag_max = f_mag[i];

    for(int i = 0; i < f_mag.size(); i++)
        f_mag[i] = f_mag[i]/f_mag_max;

}

std::vector<unsigned int> &SolutionMapper::getTriangles()
{
    return surfaceTriangles;
}

std::vector<float> &SolutionMapper::getVertices()
{
    return points;
}

std::vector<float> &SolutionMapper::getF_mag()
{
    return f_mag;
}
