cmake_minimum_required(VERSION 3.1)

project(ngsolve-apps)

cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)
cmake_policy(SET CMP0072 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

if(WIN32)
    find_package(GMP REQUIRED)
else()
    set(GMP_LIBRARIES gmp)
endif()
find_package(NGSolve CONFIG REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Qt REQUIRED COMPONENTS Quick QuickControls2 Widgets Gui)
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
#include_directories(/home/kai/Development/github/FloatTetwild/3rdparty/libigl/external/eigen)
#find_package(Eigen3)
#find_package(libigl)


set(TETWILD_ROOT $ENV{TETWILD_DIR})
set(TETWILD_INCLUDE_DIR ${TETWILD_ROOT}/include/)
set(TETWILD_LIB_DIR ${TETWILD_ROOT}/lib/)

include_directories(${TETWILD_ROOT}/../3rdparty/include)
include_directories($ENV{IGL_INCLUDE_DIR})

include_directories(src)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${VTK_INCLUDE_DIR})
include_directories(${TETWILD_INCLUDE_DIR})
include_directories(${TETWILD_INCLUDE_DIR}/geogram1)
link_directories(${TETWILD_LIB_DIR})
include_directories(${TETWILD_ROOT}/../3rdparty/lib)
add_definitions(-DEIGEN_DONT_VECTORIZE)

message("vtk libs : ${VTK_LIBRARIES}")

set(SRC_PDE
   src/pde/material.cpp
   src/pde/elasticitysystem.cpp
   src/pde/eigensystemsolver.cpp
   src/pde/solutionmapper.cpp)

# add_executable(Elasticity-Eigenproblem
#    src/Elasticity-Eigenproblem.cpp
#    ${SRC_PDE}
#    src/Meshing/IOUtils.cpp
#    src/Meshing/IOUtils.h
#    src/Meshing/MshLoader.cpp
#    src/Meshing/MshLoader.h
#    src/Meshing/MshSaver.cpp
#    src/Meshing/MshSaver.h
#    src/pde/modemeanmaxvariation.h
#    src/pde/modemeanmaxvariation.cpp
# )

# target_link_libraries(Elasticity-Eigenproblem
#    ngcomp
#    Qt5::Gui
#    Qt5::Quick
#    Qt5::Widgets
#    ${VTK_LIBRARIES}
#    Eigen3::Eigen)


add_executable(meshing
   src/meshing.cpp
   src/Meshing/IOUtils.cpp
   src/Meshing/IOUtils.h
   src/Meshing/MshLoader.cpp
   src/Meshing/MshLoader.h
   src/Meshing/MshSaver.cpp
   src/Meshing/MshSaver.h
   src/Meshing/tettools.h
   src/Meshing/tettools.cpp
)

target_link_libraries(meshing
   ngcomp
   ${GMP_LIBRARIES}
   OpenMP::OpenMP_CXX
   FloatTetwild
   predicates
   fmt
   tbb_static
   geogram
   FastEnvelope
   IndirectPredicates
)

add_executable(qvtktest
    src/main.cpp
    src/Ui/vtksceneeventprocessor.h
    src/Ui/vtksceneeventprocessor.cpp
    src/Ui/vtkorientationindicator.h
    src/Ui/vtkorientationindicator.cpp
    src/Ui/meshingcontroller.h
    src/Ui/meshingcontroller.cpp
    src/Ui/colormaphelper.h
    src/Ui/colormaphelper.cpp
    src/Ui/cellpicking.h
    src/Ui/cellpicking.cpp
    src/Ui/resources.qrc
    src/Ui/vtkqtcontroller.h
    src/Ui/vtkqtcontroller.cpp
    src/Ui/qmlpanemeshinterface.h
    src/Ui/qmlpanemeshinterface.cpp
    src/tinyxml2/tinyxml2.cpp
    src/tinyxml2/tinyxml2.h
    src/Meshing/graph.h
    src/Meshing/graph.cpp
    src/Meshing/IOUtils.cpp
    src/Meshing/IOUtils.h
    src/Meshing/MshLoader.cpp
    src/Meshing/MshLoader.h
    src/Meshing/MshSaver.cpp
    src/Meshing/MshSaver.h
    src/Meshing/tettools.h
    src/Meshing/tettools.cpp
    src/Ui/meshingworker.h
    src/Ui/meshingworker.cpp
    src/vtk/vtkmeshwrapper.h
    src/vtk/vtkmeshwrapper.cpp
)

target_link_libraries(qvtktest
    Qt::Gui
    Qt::Quick
    Qt::Widgets
    ${VTK_LIBRARIES}
    ${GMP_LIBRARIES}
    OpenMP::OpenMP_CXX
    ngstd
#    wildmeshing.so
    FloatTetwild
    predicates
    fmt
    tbb_static
    geogram
    FastEnvelope
    IndirectPredicates
)
