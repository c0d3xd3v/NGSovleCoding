#ifndef CELLPICKING_H
#define CELLPICKING_H

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <vtkCellPicker.h>
#include <vtkExtractSelection.h>

#include "meshingcontroller.h"

class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
private:
    std::map<std::string, MeshingController *> *meshRenderController;
public:
    static MouseInteractorStyle* New();
    MouseInteractorStyle();
    virtual void OnLeftButtonDown() override;
    void setMeshRenderController(std::map<std::string, MeshingController *> *newMeshRenderController);
};

#endif // CELLPICKING_H
