#ifndef TETTOOLS_H
#define TETTOOLS_H

#include <Eigen/Dense>
#include <nglib.h>
#include <meshing/meshing.hpp>

void loadMSH(std::string &path, Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets);
void getSurfaceTriangles(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tets, Eigen::MatrixXi &F2);
netgen::Mesh *generateNGMesh(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets);

#endif // TETTOOLS_H
