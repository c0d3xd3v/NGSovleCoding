#include "tettools.h"

#include <Meshing/MshLoader.h>
#include <floattetwild/MeshImprovement.h>

/*
    #include <tetrahedralize.hpp>
    double stop_quality = 10;
    int max_its = 100;
    int stage = 2;
    int stop_p = -1;
    double epsilon = 0.001;
    double edge_length_r = 0.025;
    bool skip_simplify = false;
    bool coarsen = false;
    wildmeshing_binding::Tetrahedralizer tetra(stop_quality, max_its, stage, stop_p,
                                               epsilon, edge_length_r,
                                               skip_simplify, coarsen);
*/

void getSurfaceTriangles(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tets, Eigen::MatrixXi &F2)
{
    floatTetWild::Mesh ftMesh;

    for(unsigned int i = 0; i < nodes.rows(); i++)
        ftMesh.tet_vertices.push_back(floatTetWild::MeshVertex(nodes.row(i)));

    for(unsigned int i = 0; i < tets.rows(); i++)
        ftMesh.tets.push_back(floatTetWild::MeshTet(tets.row(i)));

    ftMesh.is_input_all_inserted = true;

    Eigen::MatrixXd *V = Eigen::internal::aligned_new<Eigen::MatrixXd>(1);
    Eigen::MatrixXi *F = Eigen::internal::aligned_new<Eigen::MatrixXi>(1);
    floatTetWild::manifold_surface(ftMesh, *V, *F);

    Eigen::VectorXi tags(V->rows());
    tags.setZero(V->rows());
    for(int i = 0; i < V->rows(); i++)
    {
        Eigen::Vector3d vi = V->row(i);
        for(int j = 0; j < nodes.rows(); j++)
        {
            Eigen::Vector3d vj = nodes.row(j);
            if((vi - vj).norm() < 0.000001)
            {
                tags(i) = j;
                break;
            }
        }
    }

    F2.setZero(F->rows(), F->cols());
    for(int i = 0; i < F->rows(); i++)
        for(int j = 0; j < F->cols(); j++)
            F2(i,j) = tags((*F)(i,j));

    Eigen::internal::aligned_free(V);
    Eigen::internal::aligned_free(F);
}

netgen::Mesh *generateNGMesh(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets)
{
    netgen::Mesh *mesh = new netgen::Mesh();

    for(unsigned int i = 0; i < nodes.rows(); i++)
    {
        netgen::Point3d p(
            nodes.row(i)[0],
            nodes.row(i)[1],
            nodes.row(i)[2]);
        mesh->AddPoint(p);
    }

    netgen::FaceDescriptor fd(0, 0, 1, 0);
    int si = mesh->AddFaceDescriptor(fd);

    for(unsigned int i = 0; i < tris.rows(); i++)
    {
        netgen::Element2d el(
            tris.row(i)[0]+1,
            tris.row(i)[1]+1,
            tris.row(i)[2]+1);
        el.SetIndex(si);
        mesh->AddSurfaceElement(el);
    }

    for(unsigned int i = 0; i < tets.rows(); i++)
    {
        netgen::Element el(netgen::ELEMENT_TYPE::TET);
        el.PNum(1) = tets.row(i)[0]+1;
        el.PNum(3) = tets.row(i)[1]+1;
        el.PNum(2) = tets.row(i)[2]+1;
        el.PNum(4) = tets.row(i)[3]+1;
        el.SetIndex(0);
        mesh->AddVolumeElement(el);
    }

    return mesh;
}

void loadMSH(string &path, Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets)
{
    PyMesh::MshLoader mshLoader(path);
    nodes = Eigen::MatrixXd(mshLoader.get_nodes());
    tets = Eigen::MatrixXi(mshLoader.get_elements());

    nodes.resize(3, nodes.rows()/3);
    nodes.transposeInPlace();

    tets.resize(4, tets.rows()/4);
    tets.transposeInPlace();

    getSurfaceTriangles(nodes, tets,  tris);
}
