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

####################### eigen ###################

include(FetchContent)
FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2
    URL_MD5 b9e98a200d2455f06db9c661c5610496
)
FetchContent_GetProperties(eigen)
FetchContent_Populate(eigen)
include_directories(${eigen_SOURCE_DIR})


####################### lbfgs ###################

FetchContent_Declare(
    lbfgs
    GIT_REPOSITORY https://github.com/yixuan/LBFGSpp.git
    GIT_TAG 7fea82aab31607fc51c2ffd751a21c031b9a1061
)
FetchContent_GetProperties(lbfgs)
FetchContent_Populate(lbfgs)
include_directories(${lbfgs_SOURCE_DIR}/include)

############### libigl #######################

option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_GLFW              "Use GLFW"           ON)
option(LIBIGL_WITH_IMGUI             "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG 237ffa20b1e2fd92e0be7862e9246c3774e4fc14
)
FetchContent_MakeAvailable(libigl)


############## OpenMesh #################

FetchContent_Declare(
    openmesh
    GIT_REPOSITORY https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
    GIT_TAG d50cad4640d6d3657c4d0188fbf27ff38e4bfdca
)
FetchContent_MakeAvailable(openmesh)
include_directories(${openmesh_SOURCE_DIR}/src)


# ####################### precompiled dependencies ####################

# find_package(OpenMesh REQUIRED)
# include_directories(${OpenMesh_INCLUDE_DIRS})
##################### import files ##################
add_definitions(-D_USE_MATH_DEFINES)
add_subdirectory(src)