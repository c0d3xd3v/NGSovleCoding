project(ngsolve-apps)
cmake_minimum_required(VERSION 3.1)
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)
cmake_policy(SET CMP0072 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. 
)

find_package(Qt5 REQUIRED COMPONENTS Quick Core Gui Widgets Qml)
find_package(VTK REQUIRED COMPONENTS
    GUISupportQtQuick
    CommonColor
    CommonCore
    CommonMath
    CommonTransforms
    FiltersCore
    FiltersGeneral
    IOImage
    ImagingCore
    ImagingGeneral
    ImagingMorphological
    InteractionStyle
    InteractionWidgets
    RenderingAnnotation
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2
    cli11
    jsoncpp)

include_directories(src)

add_executable(Elasticity-Eigenproblem
    Elasticity-Eigenproblem.cpp
    material.cpp
    elasticitysystem.cpp
    eigensystemsolver.cpp
    solutionmapper.cpp
)

target_link_libraries(Elasticity-Eigenproblem
    ngcomp
    Qt5::Gui
    Qt5::Quick
    ${VTK_LIBRARIES})



add_executable(qvtktest src/main.cpp
    src/Ui/main.qml
    src/Ui/vtksceneeventprocessor.h
    src/Ui/vtksceneeventprocessor.cpp
    src/Ui/vtkorientationindicator.h
    src/Ui/vtkorientationindicator.cpp
    src/Ui/resources.qrc)

message("vtk libs : ${VTK_LIBRARIES}")

target_link_libraries(qvtktest
    ngcomp
    ${VTK_LIBRARIES}
    Qt5::Gui
    Qt5::Quick)





