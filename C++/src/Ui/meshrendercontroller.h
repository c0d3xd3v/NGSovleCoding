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

class MeshRenderController : public QObject, public GraphCompareCondition
{
    Q_OBJECT
private:
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
    MeshRenderController(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris);
    vtkActor *getActor();
    vtkSmartPointer<vtkPolyData> getPolydata();
    void selectCell(vtkIdType vtkId);
    void setColormap(ColorMap cm);
};

#endif // MESHRENDERCONTROLLER_H
