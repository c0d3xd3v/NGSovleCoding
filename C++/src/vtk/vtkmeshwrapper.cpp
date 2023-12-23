#include <Eigen/Dense>

#include <vtk/vtkFloatArray.h>
#include <vtk/vtkOrientationMarkerWidget.h>
#include <vtk/vtkCellData.h>
#include <vtk/vtkLookupTable.h>
#include <vtk/vtkProperty.h>
#include <vtk/vtkCell.h>
#include <vtk/vtkNamedColors.h>


#include "vtkmeshwrapper.h"


VtkMeshWrapper::VtkMeshWrapper(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris)
{
    points = vtkSmartPointer<vtkPoints>::New();
    cells = vtkSmartPointer<vtkCellArray>::New();
    mesh = vtkSmartPointer<vtkPolyData>::New();
    triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    polynormals = vtkSmartPointer<vtkPolyDataNormals>::New();
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    actor = vtkSmartPointer<vtkActor>::New();

    for(unsigned int i = 0; i < nodes.rows(); i++)
    {
        points->InsertNextPoint(
            nodes.row(i)[0],
            nodes.row(i)[1],
            nodes.row(i)[2]
            );
    }

    for(unsigned int i = 0; i < tris.rows(); i++)
    {
        vtkIdType v1 =  tris.row(i)[0];
        vtkIdType v2 =  tris.row(i)[1];
        vtkIdType v3 =  tris.row(i)[2];

        vtkNew<vtkTriangle> triangle;
        triangle->GetPointIds()->SetId(0, v1);
        triangle->GetPointIds()->SetId(1, v2);
        triangle->GetPointIds()->SetId(2, v3);

        cells->InsertNextCell(triangle);
    }

    mesh->SetPoints(points);
    mesh->SetPolys(cells);

    triangleFilter->SetInputData(mesh);
    triangleFilter->Update();

    polynormals->SetInputData(triangleFilter->GetOutput());
    polynormals->ComputeCellNormalsOn();
    polynormals->ConsistencyOn();
    polynormals->Update();

    initGrouppingCellData(mesh);

    meshGraph = nullptr;
    buildMeshGraph();

    mapper->SetInputData(mesh);
    actor->SetMapper(mapper);
}

void VtkMeshWrapper::selectCell(vtkIdType vtkId)
{
    meshGraph->reset();
    std::vector<vtkIdType> gids = meshGraph->BFS(vtkId);

    vtkCellData *cellData = mesh->GetCellData();
    vtkDataArray *scalars = cellData->GetScalars();
    int  selectColor = (scalars->GetTuple1(vtkId) == 1000) ? 0 : 1000;

    for (auto it2 = gids.begin(); it2 != gids.end(); ++it2)
        scalars->SetTuple1(*it2, selectColor);

    mesh->Modified();
    mapper->UpdateDataObject();
}

void VtkMeshWrapper::initGrouppingCellData(vtkSmartPointer<vtkPolyData> polydata)
{
    vtkNew<vtkFloatArray> tgroupping;
    tgroupping->SetNumberOfComponents(1);
    tgroupping->SetName("groupping");

    int n = polydata->GetNumberOfCells();
    for(int i = 0; i < n; i++)
        tgroupping->InsertNextTuple1(0);

    polydata->GetCellData()->SetScalars(tgroupping);
}

void VtkMeshWrapper::buildMeshGraph()
{
    vtkSmartPointer<vtkPolyData> meshData = triangleFilter->GetOutput();
    vtkIdType nbrOfCells = meshData->GetNumberOfCells();

    vtkIdType n = mesh->GetNumberOfCells();
    if(meshGraph != nullptr) delete meshGraph;
    meshGraph = new Graph(n, this);
    meshGraph->reset();

    for(vtkIdType currCellId = 0; currCellId < nbrOfCells; currCellId++)
    {
        vtkSmartPointer<vtkIdList> currCellPointIds = vtkSmartPointer<vtkIdList>::New();
        meshData->GetCellPoints(currCellId, currCellPointIds);

        for(vtkIdType i = 0; i < currCellPointIds->GetNumberOfIds(); i++)
        {
            vtkIdType currCellPointId = currCellPointIds->GetId(i);
            vtkNew<vtkIdList> idList;
            idList->InsertNextId(currCellPointId);
            // Get the neighbors of the cell.
            vtkNew<vtkIdList> currPointCellNeighborIds;
            meshData->GetCellNeighbors(currCellId, idList, currPointCellNeighborIds);
            for (vtkIdType j = 0; j < currPointCellNeighborIds->GetNumberOfIds(); j++)
            {
                vtkIdType cellId = currPointCellNeighborIds->GetId(j);
                meshGraph->addEdge(currCellId, cellId);
            }
        }
    }
}

bool VtkMeshWrapper::compare(vtkIdType f1, vtkIdType f2)
{
    vtkDataArray* normalDataFloat = polynormals->GetOutput()->GetCellData()->GetNormals();
    Eigen::Vector3d n1(normalDataFloat->GetTuple3(f1));
    Eigen::Vector3d n2(normalDataFloat->GetTuple3(f2));
    double threshold = 0.99;
    return (n1.dot(n2) >= threshold);
}

void VtkMeshWrapper::setColormap(ColorMap cm)
{
    mapper->SetLookupTable(cm.ctf);
}

vtkActor *VtkMeshWrapper::getActor()
{
    return actor;
}

vtkSmartPointer<vtkPolyData> VtkMeshWrapper::getPolydata()
{
    return mesh;
}

void VtkMeshWrapper::renderSurfaceWithEdges()
{
    actor->GetProperty()->EdgeVisibilityOn();
}

bool VtkMeshWrapper::getRenderSurfaceWidthEdges()
{
    return actor->GetProperty()->GetEdgeVisibility();
}

void VtkMeshWrapper::renderSurface()
{
    actor->GetProperty()->EdgeVisibilityOff();
}

bool VtkMeshWrapper::getRenderSurface()
{
    return !actor->GetProperty()->GetEdgeVisibility();
}

void VtkMeshWrapper::visibilityOff()
{
    actor->VisibilityOff();
}

void VtkMeshWrapper::visibilityOn()
{
    actor->VisibilityOn();
}

