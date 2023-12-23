#ifndef MESHINGCONTROLLER_H
#define MESHINGCONTROLLER_H


#include <QObject>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <floattetwild/ftetwildwrapper.h>

#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkActor.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkTriangle.h>
#include <vtk/vtkPolyData.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkCellArray.h>

#include "colormaphelper.h"
#include "Meshing/graph.h"

#include "vtk/vtkmeshwrapper.h"

class MeshingController : public QObject
{
    Q_OBJECT
private:
    VtkMeshWrapper* currentMesh;
    VtkMeshWrapper* sourceMesh;
    VtkMeshWrapper* resultMesh;

    FTetWildWrapper* ftetwildWrapper;
    Eigen::MatrixXf nodes;
    Eigen::MatrixXi tris;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MeshingController(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris);
    void domeshing(double stop_energy, double rel_edge_length, double rel_eps);
    void selectCell(vtkIdType vtkId);
    void setColormap(ColorMap cm);
    vtkActor *getActor();
    vtkSmartPointer<vtkPolyData> getPolydata();
    void renderSurfaceWithEdges();
    void renderSurface();
    void setCurrentMeshToSource();
    void setCurrentMeshToResult();
    void hideCurrentMesh();
    void showCurrentMesh();
    bool getRenderSurfaceWidthEdges();
};

#endif // MESHINGCONTROLLER_H
