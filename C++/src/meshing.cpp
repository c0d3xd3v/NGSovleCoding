#include <igl/writeOBJ.h>

#include "Meshing/tettools.h"
#include "Ui/instantrendering.h"

int main (int argc, char ** argv)
{
    std::string path = argv[1]; //"/home/kai/Development/github/FloatTetwild/build/ridex-Body004.obj_.msh";
    std::cout << "path : " << path << std::endl;

    Eigen::MatrixXf nodes;
    Eigen::MatrixXi tris;
    Eigen::MatrixXi tets;
    loadMSH(path, nodes, tris, tets);
    netgen::Mesh *mesh = generateNGMesh(nodes, tris, tets);

    std::cout << "CheckVolumeMesh : " << mesh->CheckVolumeMesh() << std::endl;
    std::cout << "points : " << nodes.rows() << std::endl;
    std::cout << "tets   : " << tets.rows() << std::endl;
    std::cout << "tris   : " << tris.rows() << std::endl;

    igl::writeOBJ("test.obj", nodes, tris);
    mesh->Save("test.vol");

    renderTriangleMesh(nodes, tris);

    return 0;
}
