#ifndef VTKMESHWRAPPER_H
#define VTKMESHWRAPPER_H

#include <Eigen/Dense>

#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataNormals.h>

#include "Ui/colormaphelper.h"
#include "Meshing/graph.h"

class VtkMeshWrapper : public GraphCompareCondition
{
private:
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkCellArray> cells;
    vtkSmartPointer<vtkPolyData> mesh;
    vtkSmartPointer<vtkTriangleFilter> triangleFilter;
    vtkSmartPointer<vtkPolyDataNormals> polynormals;
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;

    Graph *meshGraph;

    void initGrouppingCellData(vtkSmartPointer<vtkPolyData> polydata);
    void buildMeshGraph();
    bool compare(vtkIdType f1, vtkIdType f2);

public:
    VtkMeshWrapper(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris);
    vtkActor *getActor();
    vtkSmartPointer<vtkPolyData> getPolydata();
    void selectCell(vtkIdType vtkId);
    void setColormap(ColorMap cm);
    void renderSurfaceWithEdges();
    bool getRenderSurfaceWidthEdges();
    void renderSurface();
    bool getRenderSurface();
    void visibilityOff();
    void visibilityOn();
};

#endif // VTKMESHWRAPPER_H
