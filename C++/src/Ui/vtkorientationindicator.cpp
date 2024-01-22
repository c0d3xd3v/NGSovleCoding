#include "vtkorientationindicator.h"


#include <vtkProperty.h>
#include <vtkConeSource.h>
#include <vtkTextProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkCaptionActor2D.h>


vtkSmartPointer<vtkAxesActor> MakeAxesActor(std::array<double, 3>& scale/*,
                                   std::array<std::string, 3> const& xyzLabels*/
                                            )
{
    vtkNew<vtkAxesActor> axes;
    axes->SetScale(scale.data());
    axes->SetShaftTypeToCylinder();

    axes->SetXAxisLabelText("");
    axes->SetYAxisLabelText("");
    axes->SetZAxisLabelText("");

    axes->SetCylinderRadius(1.2 * axes->GetCylinderRadius());
    axes->SetConeRadius(1.5 * axes->GetConeRadius());
    axes->SetSphereRadius(1.5 * axes->GetSphereRadius());
    /*
    auto tprop = axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty();
    tprop->ItalicOn();
    tprop->ShadowOn();
    tprop->SetFontFamilyToTimes();
    // Use the same text properties on the other two axes.
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy(tprop);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy(tprop);
    */
    return axes;
}

vtkSmartPointer<vtkAnnotatedCubeActor>
MakeAnnotatedCubeActor(std::array<std::string, 6> const& cubeLabels,
                       vtkNamedColors* colors)
{
    // A cube with labeled faces.
    vtkNew<vtkAnnotatedCubeActor> cube;
    cube->SetXPlusFaceText(cubeLabels[0].c_str());
    cube->SetXMinusFaceText(cubeLabels[1].c_str());
    cube->SetYPlusFaceText(cubeLabels[2].c_str());
    cube->SetYMinusFaceText(cubeLabels[3].c_str());
    cube->SetZPlusFaceText(cubeLabels[4].c_str());
    cube->SetZMinusFaceText(cubeLabels[5].c_str());
    cube->SetFaceTextScale(0.5);
    cube->GetCubeProperty()->SetColor(colors->GetColor3d("Gainsboro").GetData());

    cube->GetTextEdgesProperty()->SetColor(
        colors->GetColor3d("LightSlateGray").GetData());

    // Change the vector text colors.
    cube->GetXPlusFaceProperty()->SetColor(
        colors->GetColor3d("Tomato").GetData());
    cube->GetXMinusFaceProperty()->SetColor(
        colors->GetColor3d("Tomato").GetData());
    cube->GetYPlusFaceProperty()->SetColor(
        colors->GetColor3d("DeepSkyBlue").GetData());
    cube->GetYMinusFaceProperty()->SetColor(
        colors->GetColor3d("DeepSkyBlue").GetData());
    cube->GetZPlusFaceProperty()->SetColor(
        colors->GetColor3d("SeaGreen").GetData());
    cube->GetZMinusFaceProperty()->SetColor(
        colors->GetColor3d("SeaGreen").GetData());
    return cube;
}

vtkSmartPointer<vtkPropAssembly> MakeCubeActor(std::string const& labelSelector,
                                      vtkNamedColors* colors)
{
    std::array<std::string, 3> xyzLabels;
    std::array<std::string, 6> cubeLabels;
    std::array<double, 3> scale;
    if (labelSelector == "sal")
    {
        // xyzLabels = std::array<std::string,3>{"S", "A", "L"};
        xyzLabels = std::array<std::string, 3>{"+X", "+Y", "+Z"};
        cubeLabels = std::array<std::string, 6>{"S", "I", "A", "P", "L", "R"};
        scale = std::array<double, 3>{1.5, 1.5, 1.5};
    }
    else if (labelSelector == "rsp")
    {
        // xyzLabels = std::array<std::string, 3>{"R", "S", "P"};
        xyzLabels = std::array<std::string, 3>{"+X", "+Y", "+Z"};
        cubeLabels = std::array<std::string, 6>{"R", "L", "S", "I", "P", "A"};
        scale = std::array<double, 3>{1.5, 1.5, 1.5};
    }
    else if (labelSelector == "lsa")
    {
        // xyzLabels = std::array<std::string, 3>{"L", "S", "A"};
        xyzLabels = std::array<std::string, 3>{"+X", "+Y", "+Z"};
        cubeLabels = std::array<std::string, 6>{"L", "R", "S", "I", "A", "P"};
        scale = std::array<double, 3>{1.5, 1.5, 1.5};
    }
    else
    {
        xyzLabels = std::array<std::string, 3>{"+X", "+Y", "+Z"};
        cubeLabels = std::array<std::string, 6>{"+X", "-X", "+Y", "-Y", "+Z", "-Z"};
        scale = std::array<double, 3>{1.5, 1.5, 1.5};
    }

    // We are combining a vtkAxesActor and a vtkAnnotatedCubeActor
    // into a vtkPropAssembly
    auto cube = MakeAnnotatedCubeActor(cubeLabels, colors);

    auto axes = MakeAxesActor(scale/*, xyzLabels*/);

    // Combine orientation markers into one with an assembly.
    vtkNew<vtkPropAssembly> assembly;
    assembly->AddPart(axes);
    assembly->AddPart(cube);
    return assembly;
}
