#ifndef SOLUTIONMAPPER_H
#define SOLUTIONMAPPER_H
#undef NETGEN_PYTHON
#include <comp.hpp>
using namespace ngcomp;

class SolutionMapper
{
private:
    std::shared_ptr<MeshAccess> ma;
    std::vector<unsigned int> surfaceVertices;
    std::vector<float> points;
    std::vector<unsigned int> surfaceTriangles;
    std::vector<std::complex<double>[3]> fv;
    std::vector<float> fvr;
    std::vector<float> f_mag;

public:
    SolutionMapper(std::shared_ptr<MeshAccess> ma);

    void mapSolutionToMesh(std::shared_ptr<VectorH1FESpace> fes, Array<shared_ptr<BaseVector>> &evecs);
    std::vector<unsigned int> &getTriangles();
    std::vector<float> &getVertices();
    std::vector<float> &getF_mag();
};
#endif // SOLUTIONMAPPER_H
