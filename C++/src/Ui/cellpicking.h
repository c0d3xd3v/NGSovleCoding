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

#include "meshrendercontroller.h"

class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static MouseInteractorStyle* New();

    MouseInteractorStyle();

    virtual void OnLeftButtonDown() override;

    vtkSmartPointer<vtkPolyData> Data;
    MeshRenderController *rc;
    //vtkSmartPointer<vtkDataSetMapper> selectedMapper;
    //vtkSmartPointer<vtkActor> selectedActor;
};

#endif // CELLPICKING_H
