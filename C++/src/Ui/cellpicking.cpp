#include <QDebug>

#include <vtk/vtkCellData.h>
#include <vtk/vtkPolyDataMapper.h>
#include <vtk/vtkActor.h>

#include "cellpicking.h"

vtkStandardNewMacro(MouseInteractorStyle)

    MouseInteractorStyle::MouseInteractorStyle() : rc(nullptr)
{
    //selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    //selectedActor = vtkSmartPointer<vtkActor>::New();
}

void MouseInteractorStyle::OnLeftButtonDown()
{
    int* pos = this->GetInteractor()->GetEventPosition();

    vtkNew<vtkCellPicker> picker;
    picker->SetTolerance(0.0005);

    // Pick from this location.
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
    std::cout << "Cell id is: " << picker->GetCellId() << std::endl;
    std::stringstream ss;
    ss << picker->GetActor();

    qDebug() << ss.str().c_str();

    if (picker->GetCellId() != -1)
    {
        if(rc != nullptr)
            rc->selectCell(picker->GetCellId());
    }
    // Forward events.
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}
