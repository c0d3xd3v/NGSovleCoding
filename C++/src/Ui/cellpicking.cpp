#include "cellpicking.h"

vtkStandardNewMacro(MouseInteractorStyle)

MouseInteractorStyle::MouseInteractorStyle() :
    meshRenderController(nullptr)
{}

void MouseInteractorStyle::OnLeftButtonDown()
{
    int* pos = this->GetInteractor()->GetEventPosition();

    vtkNew<vtkCellPicker> picker;
    picker->SetTolerance(0.0005);
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    if (picker->GetCellId() != -1)
    {
        if(meshRenderController != nullptr)
        {
            typedef std::map<std::string, MeshRenderController*>::iterator
                Iterator;
            std::stringstream ss;
            ss << picker->GetActor();
            Iterator itr = meshRenderController->find(ss.str());
            if(itr != meshRenderController->end())
            {
                MeshRenderController* mrc = itr->second;
                mrc->selectCell(picker->GetCellId());
            }
        }
    }
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void MouseInteractorStyle::setMeshRenderController(std::map<std::string, MeshRenderController *> *newMeshRenderController)
{
    meshRenderController = newMeshRenderController;
}
