#ifndef VTKORIENTATIONINDICATOR_H
#define VTKORIENTATIONINDICATOR_H

#include <vtk/vtkAxesActor.h>
#include <vtk/vtkNamedColors.h>
#include <vtk/vtkPropAssembly.h>
#include <vtk/vtkAnnotatedCubeActor.h>

vtkNew<vtkAxesActor> MakeAxesActor(std::array<double, 3>& scale,
                                   std::array<std::string, 3> const& xyzLabels);

vtkNew<vtkAnnotatedCubeActor>
MakeAnnotatedCubeActor(std::array<std::string, 6> const& cubeLabels,
                       vtkNamedColors* colors);

vtkNew<vtkPropAssembly> MakeCubeActor(std::string const& labelSelector,
                                      vtkNamedColors* colors);


#endif // VTKORIENTATIONINDICATOR_H
