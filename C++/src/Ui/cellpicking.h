#ifndef CELLPICKING_H
#define CELLPICKING_H

#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkInteractorStyleTrackballCamera.h>

#include <vtk/vtkIdTypeArray.h>
#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkPolyData.h>

#include <vtk/vtkCellPicker.h>
#include <vtk/vtkExtractSelection.h>

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
