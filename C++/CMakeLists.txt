project(ngsolve-apps)
cmake_minimum_required(VERSION 3.1)
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)
cmake_policy(SET CMP0072 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. 
)

add_executable(Elasticity-Eigenproblem
    Elasticity-Eigenproblem.cpp
    material.cpp
    elasticitysystem.cpp
    eigensystemsolver.cpp
    solutionmapper.cpp)

target_link_libraries(Elasticity-Eigenproblem ngcomp)





