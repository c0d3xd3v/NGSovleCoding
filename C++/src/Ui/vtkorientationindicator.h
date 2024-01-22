#ifndef VTKORIENTATIONINDICATOR_H
#define VTKORIENTATIONINDICATOR_H

#include <vtkAxesActor.h>
#include <vtkNamedColors.h>
#include <vtkPropAssembly.h>
#include <vtkAnnotatedCubeActor.h>

vtkSmartPointer<vtkAxesActor> MakeAxesActor(std::array<double, 3>& scale,
                                   std::array<std::string, 3> const& xyzLabels);

vtkSmartPointer<vtkAnnotatedCubeActor>
MakeAnnotatedCubeActor(std::array<std::string, 6> const& cubeLabels,
                       vtkNamedColors* colors);

vtkSmartPointer<vtkPropAssembly> MakeCubeActor(std::string const& labelSelector,
                                      vtkNamedColors* colors);


#endif // VTKORIENTATIONINDICATOR_H
