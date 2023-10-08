#ifndef MESHRENDERCONTROLLER_H
#define MESHRENDERCONTROLLER_H

#include <QObject>

#include <Eigen/Dense>
#include <Eigen/Sparse>

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

class MeshRenderController : public QObject, public GraphCompareCondition
{
    Q_OBJECT
private:
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> cells;
    vtkNew<vtkPolyData> mesh;
    vtkNew<vtkTriangleFilter> triangleFilter;
    vtkNew<vtkPolyDataNormals> polynormals;
    vtkNew<vtkPolyDataMapper> mapper;
    vtkNew<vtkActor> actor;

    //Eigen::SparseMatrix<int> adjacency;
    Graph *meshGraph;

    //std::list<vtkIdType> findCellNeighbours(vtkIdType id, Eigen::Vector3d &N, std::list<vtkIdType> &visited_non_neighbors);
    void initGrouppingCellData(vtkSmartPointer<vtkPolyData> polydata);
    void buildMeshGraph();
    bool compare(vtkIdType f1, vtkIdType f2);

public:
    MeshRenderController(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris);
    vtkActor *getActor();
    vtkSmartPointer<vtkPolyData> getPolydata();
    void selectCell(vtkIdType vtkId);
    void setColormap(ColorMap cm);
};

#endif // MESHRENDERCONTROLLER_H
