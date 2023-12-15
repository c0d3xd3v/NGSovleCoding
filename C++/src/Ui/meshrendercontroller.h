#ifndef MESHRENDERCONTROLLER_H
#define MESHRENDERCONTROLLER_H


#include <QObject>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkActor.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkTriangle.h>
#include <vtk/vtkPolyData.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkPolyDataMapper.h>
#include <vtk/vtkTriangleFilter.h>
#include <vtk/vtkPolyDataNormals.h>

#include "colormaphelper.h"
#include "Meshing/graph.h"
#include <floattetwild/ftetwildwrapper.h>

class MeshRenderController : public QObject, public GraphCompareCondition
{
    Q_OBJECT
private:
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    FTetWildWrapper* ftetwildWrapper;
    Eigen::MatrixXf nodes;
    Eigen::MatrixXi tris;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkCellArray> cells;
    vtkSmartPointer<vtkPolyData> mesh;
    vtkSmartPointer<vtkTriangleFilter> triangleFilter;
    vtkSmartPointer<vtkPolyDataNormals> polynormals;
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;

    //Eigen::SparseMatrix<int> adjacency;
    Graph *meshGraph;

    //std::list<vtkIdType> findCellNeighbours(vtkIdType id, Eigen::Vector3d &N, std::list<vtkIdType> &visited_non_neighbors);
    void initGrouppingCellData(vtkSmartPointer<vtkPolyData> polydata);
    void buildMeshGraph();
    bool compare(vtkIdType f1, vtkIdType f2);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MeshRenderController(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris);
    vtkActor *getActor();
    vtkSmartPointer<vtkPolyData> getPolydata();
    void selectCell(vtkIdType vtkId);
    void setColormap(ColorMap cm);
    void domeshing(double stop_energy, double rel_edge_length, double rel_eps);
    void renderSurfaceWithEdges();
    void renderSurface();
};

#endif // MESHRENDERCONTROLLER_H
