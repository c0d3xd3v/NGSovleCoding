#include <Eigen/Dense>
#include <igl/read_triangle_mesh.h>
#include <floattetwild/ftetwildwrapper.h>

#include "meshingcontroller.h"


MeshingController::MeshingController(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris) :
    nodes(nodes), tris(tris)
{
    sourceMesh = new VtkMeshWrapper(nodes, tris);
    resultMesh = nullptr;
    currentMesh = sourceMesh;
}

void MeshingController::domeshing(double stop_energy, double rel_edge_length, double rel_eps)
{
    FTetWildWrapper* ftetwildWrapper = new FTetWildWrapper(stop_energy, rel_edge_length, rel_eps);
    ftetwildWrapper->loadMeshGeometry(nodes, tris);
    ftetwildWrapper->tetrahedralize();
    //ftetwildWrapper->save();
    Eigen::MatrixXi tris_;
    Eigen::MatrixXi tets_;
    Eigen::MatrixXf nodes_;
    ftetwildWrapper->getSurfaceIndices(tris_, tets_, nodes_);
    std::cout << "tris : " << tris_.rows() << std::endl;
    std::cout << "tets : " << tets_.rows() << std::endl;
    std::cout << "nods : " << nodes_.rows() << std::endl;
    if(resultMesh != nullptr)
        delete resultMesh;
    resultMesh = new VtkMeshWrapper(nodes_, tris_);
    resultMesh->renderSurfaceWithEdges();
    sourceMesh->visibilityOff();
    delete ftetwildWrapper;
}

void MeshingController::selectCell(vtkIdType vtkId)
{
    currentMesh->selectCell(vtkId);
}

void MeshingController::setColormap(ColorMap cm)
{
    currentMesh->setColormap(cm);
}

vtkActor *MeshingController::getActor()
{
    if(currentMesh != nullptr)
        return currentMesh->getActor();
    else
        return nullptr;
}

vtkSmartPointer<vtkPolyData> MeshingController::getPolydata()
{
    return currentMesh->getPolydata();
}

void MeshingController::renderSurfaceWithEdges()
{
    if(currentMesh != nullptr)
        currentMesh->renderSurfaceWithEdges();
}

void MeshingController::renderSurface()
{
    if(currentMesh != nullptr)
        currentMesh->renderSurface();
}

void  MeshingController::setCurrentMeshToSource()
{
    currentMesh = sourceMesh;
}

void  MeshingController::setCurrentMeshToResult()
{
    currentMesh = resultMesh;
}

void MeshingController::hideCurrentMesh()
{
    if(currentMesh != nullptr)
        currentMesh->visibilityOff();
}

void MeshingController::showCurrentMesh()
{
    if(currentMesh != nullptr)
        currentMesh->visibilityOn();
}

bool MeshingController::getRenderSurfaceWidthEdges()
{
    if(currentMesh != nullptr)
        return currentMesh->getRenderSurfaceWidthEdges();
    else return false;
}
