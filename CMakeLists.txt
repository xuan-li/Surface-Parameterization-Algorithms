cmake_minimum_required(VERSION 3.0)

project(Parameterization)

set(CMAKE_CXX_STANDARD 14)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install CACHE PATH "cmake install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if(MSVC)
  set(CMAKE_DEBUG_POSTFIX "d")
else()
  set(CMAKE_DEBUG_POSTFIX "")
endif()

# if there are some customized FindXXX modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake" CACHE STRING "Modules for CMake" FORCE)

###################### user-selected option ####################
option(WITH_OPENMP "Enable OpenMP support?" ON)
if(WITH_OPENMP)
 find_package(OpenMP REQUIRED)
 if(OPENMP_FOUND)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 endif()
 add_definitions(-DWITH_OPENMP)
endif()

####################### header-only dependencies ###################

set(EIGEN3_INCLUDE_DIR "EIGEN3-NOT-FOUND" CACHE PATH "Eigen3 include dir")
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/external/LBFGS/include)

########### libigl ##########
set(LIBIGL_ROOT external/libigl)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)

find_package(LIBIGL REQUIRED)

####################### precompiled dependencies ####################

find_package(OpenMesh REQUIRED)
include_directories(${OpenMesh_INCLUDE_DIRS})
##################### import files ##################
add_definitions(-D_USE_MATH_DEFINES)
add_subdirectory(src)