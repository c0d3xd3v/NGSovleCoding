#include <iostream>
#include <fstream>
#include <igl/readOBJ.h>
#include "Ui/instantrendering.h"

int main (int argc, char ** argv)
{
    Eigen::MatrixXd nodes;
    Eigen::MatrixXi tris;
    std::string path = "/home/kai/Development/github/NGSovleCoding/data/skull2.obj";

    std::filebuf fb;
    if(fb.open (path,std::ios::in))
    {
        std::istream is(&fb);
        igl::readOBJ(path, nodes, tris);
    }

    std::cout << "points : " << nodes.rows() << std::endl;
    std::cout << "tris   : " << tris.rows() << std::endl;

    renderTriangleMesh(nodes, tris);

    return 0;
}
