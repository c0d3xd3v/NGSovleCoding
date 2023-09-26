#include <stack>

#include <vtk/vtkFloatArray.h>
#include <vtk/vtkOrientationMarkerWidget.h>
#include <vtk/vtkCellData.h>
#include <vtk/vtkLookupTable.h>
#include <vtk/vtkProperty.h>
#include "rendercontroller.h"

#include <Eigen/Dense>

void initGrouppingCellData(vtkSmartPointer<vtkPolyData> polydata)
{
    vtkNew<vtkFloatArray> tgroupping;
    tgroupping->SetNumberOfComponents(1);
    tgroupping->SetName("groupping");

    int n = polydata->GetNumberOfCells();
    for(int i = 0; i < n; i++)
        tgroupping->InsertNextTuple1(0);

    polydata->GetCellData()->SetScalars(tgroupping);
}

RenderController::RenderController(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris)
{
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
        vtkNew<vtkTriangle> triangle;
        triangle->GetPointIds()->SetId(0, tris.row(i)[0]);
        triangle->GetPointIds()->SetId(1, tris.row(i)[1]);
        triangle->GetPointIds()->SetId(2, tris.row(i)[2]);
        cells->InsertNextCell(triangle);
    }

    mesh->SetPoints(points);
    mesh->SetPolys(cells);

    triangleFilter->SetInputData(mesh);
    triangleFilter->Update();

    polynormals->SetInputData(triangleFilter->GetOutput());
    polynormals->ComputeCellNormalsOn();
    polynormals->Update();

    initGrouppingCellData(mesh);

    mapper->SetInputData(mesh);
    actor->SetMapper(mapper);
    //actor->GetProperty()->EdgeVisibilityOn();
}

void RenderController::selectCell(vtkIdType vtkId)
{
    vtkCellData *cellData = mesh->GetCellData();    
    vtkDataArray* normalDataFloat = polynormals->GetOutput()->GetCellData()->GetNormals();

    std::list<vtkIdType> all_neighbors;
    std::stack<vtkIdType> not_validated;
    not_validated.push(vtkId);

    while(not_validated.size() > 0)
    {
        vtkIdType currId = not_validated.top();
        not_validated.pop();

        Eigen::Vector3d N(normalDataFloat->GetTuple3(currId));
        std::list<vtkIdType> neighbors = findCellNeighbours(currId, N);

        int count_added = 0;
        for (auto it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
        {
            if(std::find(
                    all_neighbors.begin(),
                    all_neighbors.end(),
                    *it1) == all_neighbors.end())
            {
                all_neighbors.push_back(*it1);
                not_validated.push(*it1);
                count_added++;
            }
        }
        all_neighbors.unique();

        std::cout << "stack size : " << not_validated.size() << std::endl;
        std::cout << "added : " << count_added << " : " << all_neighbors.size() << std::endl;
    }

    if(all_neighbors.size() == 0)
        all_neighbors.push_back(vtkId);

    vtkDataArray *scalars = cellData->GetScalars();
    int  selectColor = 1000;
    if(scalars->GetTuple1(vtkId) == selectColor)
        selectColor = 0;
    else
        selectColor = 1000;
    for (auto it2 = all_neighbors.begin(); it2 != all_neighbors.end(); ++it2)
        scalars->SetTuple1(*it2, selectColor);

    mesh->Modified();
    mapper->UpdateDataObject();
}

std::list<vtkIdType> RenderController::findCellNeighbours(vtkIdType id, Eigen::Vector3d &N)
{
    std::list<vtkIdType> neighbors;
    vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
    vtkDataArray* normalDataFloat = polynormals->GetOutput()->GetCellData()->GetNormals();
    triangleFilter->Update();
    triangleFilter->GetOutput()->GetCellPoints(id, cellPointIds);
    for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
    {
        vtkNew<vtkIdList> idList;
        idList->InsertNextId(cellPointIds->GetId(i));

        // Get the neighbors of the cell.
        vtkNew<vtkIdList> neighborCellIds;

        triangleFilter->GetOutput()->GetCellNeighbors(id, idList,
                                                      neighborCellIds);

        for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++)
        {
            vtkIdType _id = neighborCellIds->GetId(j);
            Eigen::Vector3d N2(normalDataFloat->GetTuple3(_id));
            double S = N.dot(N2);
            if(S >=0.995)
            {
                neighbors.push_back(neighborCellIds->GetId(j));
            }
            else std::cout << "discard ... " << std::endl;
        }
    }

    return neighbors;
}

void RenderController::setColormap(ColorMap cm)
{
    mapper->SetLookupTable(cm.ctf);
}

vtkActor *RenderController::getActor()
{
    return actor;
}

vtkSmartPointer<vtkPolyData> RenderController::getPolydata()
{
    return mesh;
}
