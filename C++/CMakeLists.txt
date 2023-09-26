project(ngsolve-apps)
cmake_minimum_required(VERSION 3.1)
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)
cmake_policy(SET CMP0072 NEW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. 
)

find_package(Qt5 REQUIRED COMPONENTS Quick Gui)
find_package(VTK REQUIRED COMPONENTS
    GUISupportQtQuick
    CommonColor
    CommonCore
    CommonMath
    CommonTransforms
    FiltersCore
    FiltersGeneral
    FiltersExtraction
    InteractionStyle
    InteractionWidgets
    RenderingAnnotation
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2)
find_package(Eigen3)
find_package(libigl)

set(TETWILD_INCLUDE_DIR /home/kai/Development/github/fTetWild/build/include)
set(TETWILD_LIB_DIR /home/kai/Development/libs/fTetwild/lib/)

include_directories(src)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${VTK_INCLUDE_DIR})
include_directories(${TETWILD_INCLUDE_DIR})
include_directories(/home/kai/Development/libs/fTetwild/include/geogram1)
link_directories(${TETWILD_LIB_DIR})

message("vtk libs : ${VTK_LIBRARIES}")

add_executable(Elasticity-Eigenproblem
    Elasticity-Eigenproblem.cpp
    material.cpp
    elasticitysystem.cpp
    eigensystemsolver.cpp
    solutionmapper.cpp
    src/Meshing/IOUtils.cpp
    src/Meshing/IOUtils.h
    src/Meshing/MshLoader.cpp
    src/Meshing/MshLoader.h
    src/Meshing/MshSaver.cpp
    src/Meshing/MshSaver.h
)

target_link_libraries(Elasticity-Eigenproblem
    ngcomp
    Qt5::Gui
    Qt5::Quick
    ${VTK_LIBRARIES}
    Eigen3::Eigen)


add_executable(meshing
    meshing.cpp
    src/Ui/vtksceneeventprocessor.h
    src/Ui/vtksceneeventprocessor.cpp
    src/Ui/vtkorientationindicator.h
    src/Ui/vtkorientationindicator.cpp
    src/Ui/instantrendering.h
    src/Ui/instantrendering.cpp
    src/Ui/rendercontroller.h
    src/Ui/rendercontroller.cpp
    src/Ui/cellpicking.h
    src/Ui/cellpicking.cpp
    src/Ui/resources.qrc
    src/Meshing/IOUtils.cpp
    src/Meshing/IOUtils.h
    src/Meshing/MshLoader.cpp
    src/Meshing/MshLoader.h
    src/Meshing/MshSaver.cpp
    src/Meshing/MshSaver.h
    src/Meshing/tettools.h
    src/Meshing/tettools.cpp
    src/tinyxml2/tinyxml2.cpp
    src/tinyxml2/tinyxml2.h
    src/Ui/colormaphelper.h
    src/Ui/colormaphelper.cpp
)

target_link_libraries(meshing
    Qt5::Gui
    Qt5::Quick
    ${VTK_LIBRARIES}
    Eigen3::Eigen
    igl::igl_core
    ngcomp
    gmp
    FloatTetwild
    predicates
    fmt
    tbb_static
    geogram
    FastEnvelope
    IndirectPredicates
)

add_executable(qvtktest src/main.cpp
    src/Ui/main.qml
    src/Ui/vtksceneeventprocessor.h
    src/Ui/vtksceneeventprocessor.cpp
    src/Ui/vtkorientationindicator.h
    src/Ui/vtkorientationindicator.cpp
    src/Ui/instantrendering.h
    src/Ui/instantrendering.cpp
    src/Ui/rendercontroller.h
    src/Ui/rendercontroller.cpp
    src/tinyxml2/tinyxml2.cpp
    src/tinyxml2/tinyxml2.h
    src/Ui/colormaphelper.h
    src/Ui/colormaphelper.cpp
    src/Ui/cellpicking.h
    src/Ui/cellpicking.cpp
    src/Ui/resources.qrc)

target_link_libraries(qvtktest
    Qt5::Gui
    Qt5::Quick
    ${VTK_LIBRARIES}
    Eigen3::Eigen
    igl::igl_core)







