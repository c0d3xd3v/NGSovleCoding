#ifndef RENDERCONTROLLER_H
#define RENDERCONTROLLER_H

#include <QObject>
#include <Eigen/Dense>

#include <vtk/vtkPoints.h>
#include <vtk/vtkTriangle.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkPolyData.h>
#include <vtk/vtkPolyDataMapper.h>
#include <vtk/vtkActor.h>
#include <vtk/vtkTriangleFilter.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkPolyDataNormals.h>

#include "colormaphelper.h"

class RenderController : QObject
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

    std::list<vtkIdType> findCellNeighbours(vtkIdType id, Eigen::Vector3d &N);

public:
    RenderController(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris);
    vtkActor *getActor();
    vtkSmartPointer<vtkPolyData> getPolydata();
    void selectCell(vtkIdType vtkId);
    void setColormap(ColorMap cm);
};

#endif // RENDERCONTROLLER_H
