#
# - Try to find OpenMesh
#
# Once done this will define:
#
#  OPENMESH_FOUND - system has OpenMesh
#  OpenMesh_INCLUDE_DIR - OpenMesh include directory
#  OpenMesh_LIBRARIES - Link these to use OpenMesh
#  OPENMESH_DEFINITIONS - Compiler switches required for using OpenMesh
#
find_path(OpenMesh_DIR 
    NAMES include/OpenMesh/Core/IO/Options.hh 
    PATHS ENV OpenMesh_DIR
) 

find_path(OpenMesh_INCLUDE_DIRS 
          NAMES OpenMesh/Core/IO/Options.hh
          PATHS /usr
                /usr/local
                ${OpenMesh_DIR}      
          PATH_SUFFIXES include
         )
find_library(OPENMESH_CORE_LIB NAME OpenMeshCore
    PATHS ${OpenMesh_DIR}/lib NO_DEFAULT_PATH)
find_library(OPENMESH_CORED_LIB NAME OpenMeshCored
    PATHS ${OpenMesh_DIR}/lib NO_DEFAULT_PATH)
find_library(OPENMESH_TOOLS_LIB NAME OpenMeshTools
    PATHS ${OpenMesh_DIR}/lib NO_DEFAULT_PATH)
find_library(OPENMESH_TOOLSD_LIB NAME OpenMeshToolsd
    PATHS ${OpenMesh_DIR}/lib NO_DEFAULT_PATH)

set(OpenMesh_LIBRARIES
  debug ${OPENMESH_CORED_LIB} ${OPENMESH_TOOLSD_LIB}
  optimized ${OPENMESH_CORE_LIB} ${OPENMESH_TOOLS_LIB}
)

message(${OpenMesh_LIBRARIES})

IF(OpenMesh_INCLUDE_DIRS AND OpenMesh_LIBRARIES)
   SET(OPENMESH_FOUND TRUE)
ENDIF(OpenMesh_INCLUDE_DIRS AND OpenMesh_LIBRARIES)

IF(OPENMESH_FOUND)
  IF(NOT OpenMesh_FIND_QUIETLY)
    MESSAGE(STATUS "Found OpenMesh: ${OpenMesh_INCLUDE_DIRS}")
  ENDIF(NOT OpenMesh_FIND_QUIETLY)
ELSE(OPENMESH_FOUND)
  IF(OpenMesh_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find OpenMesh")
  ENDIF(OpenMesh_FIND_REQUIRED)
ENDIF(OPENMESH_FOUND)

