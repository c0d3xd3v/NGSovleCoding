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

void getSurfaceTriangles(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tets, Eigen::MatrixXi &F2)
{
    floatTetWild::Mesh ftMesh;

    for(unsigned int i = 0; i < nodes.rows(); i++)
        ftMesh.tet_vertices.push_back(floatTetWild::MeshVertex(nodes.row(i).cast<double>()));

    for(unsigned int i = 0; i < tets.rows(); i++)
        ftMesh.tets.push_back(floatTetWild::MeshTet(tets.row(i)));

    ftMesh.is_input_all_inserted = true;

    floatTetWild::get_boundary_surface_indices(ftMesh, F2);
}

netgen::Mesh *generateNGMesh(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets)
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

    netgen::FaceDescriptor fd(1, 1, 1, 1);
    fd.SetBCName("solid");
    fd.SetBCProperty(1);
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
    mesh->SetMaterial(1, "solid");
    return mesh;
}

void loadMSH(std::string  &path, Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris, Eigen::MatrixXi &tets)
{
    PyMesh::MshLoader mshLoader(path);
    nodes = Eigen::MatrixXf(mshLoader.get_nodes().cast<float>());
    tets = Eigen::MatrixXi(mshLoader.get_elements());

    nodes.resize(3, nodes.rows()/3);
    nodes.transposeInPlace();

    tets.resize(4, tets.rows()/4);
    tets.transposeInPlace();

    getSurfaceTriangles(nodes, tets,  tris);
}
