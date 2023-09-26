#ifndef TETTOOLS_H
#define TETTOOLS_H

#include <Eigen/Dense>
#include <netgen/nglib.h>
#include <netgen/meshing/meshing.hpp>

void loadMSH(std::string &path, Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets);
void getSurfaceTriangles(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tets, Eigen::MatrixXi &F2);
netgen::Mesh *generateNGMesh(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets);

#endif // TETTOOLS_H
