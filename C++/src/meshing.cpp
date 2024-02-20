#include <Meshing/tettools.h>
#include <igl/read_triangle_mesh.h>
#include <floattetwild/ftetwildwrapper.h>


int main(int argc, char** argv)
{
    std::cout << "EIGEN_WORLD_VERSION : " << EIGEN_WORLD_VERSION << std::endl;
    std::cout << "EIGEN_MAJOR_VERSION : " << EIGEN_MAJOR_VERSION << std::endl;
    std::cout << "EIGEN_WORLD_VERSION : " << EIGEN_MINOR_VERSION << std::endl;

    const std::string oldLocale = std::setlocale(LC_NUMERIC, nullptr);
    std::setlocale(LC_NUMERIC, "C");

    std::string path = "../../data/ridex-Body004.obj";

    Eigen::MatrixXf V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(path, V, F);

    FTetWildWrapper* ftetwildWrapper = new FTetWildWrapper(40., 0.0005);
    ftetwildWrapper->loadMeshGeometry(V, F);
    ftetwildWrapper->tetrahedralize();

    Eigen::MatrixXf mNodes;
    Eigen::MatrixXi mTris;
    Eigen::MatrixXi mTets;
    ftetwildWrapper->getSurfaceIndices(mTris, mTets, mNodes);
    netgen::Mesh *ngMesh = generateNGMesh(mNodes, mTris, mTets);
    ngMesh->Save("../../data/ridex-Body004.vol");

    std::setlocale(LC_NUMERIC, oldLocale.c_str());

    return 0;
}
