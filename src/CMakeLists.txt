include_directories(Mesh Utilities BoundaryFirstFlattening OrbifoldEmbedding Viewer)

file(GLOB MESH_HEADERS Mesh/*.h)
file(GLOB MESH_SOURCES Mesh/*.cpp)
add_library(Mesh ${MESH_HEADERS} ${MESH_SOURCES})
target_link_libraries(Mesh ${OpenMesh_LIBRARIES})

file(GLOB UTILITIES_HEADERS Utilities/*.h)
file(GLOB UTILITIES_SOURCES Utilities/*.cpp)
add_library(Utilities ${UTILITIES_HEADERS} ${UTILITIES_SOURCES})
target_link_libraries(Utilities Mesh)

file(GLOB BFF_HEADERS BoundaryFirstFlattening/*.h)
file(GLOB BFF_SOURCES BoundaryFirstFlattening/*.cpp)
add_library(BoundaryFirstFlattening ${BFF_HEADERS} ${BFF_SOURCES})
target_link_libraries(BoundaryFirstFlattening Mesh igl::core)

file(GLOB OE_HEADERS OrbifoldEmbedding/*.h)
file(GLOB OE_SOURCES OrbifoldEmbedding/*.cpp)
add_library(OrbifoldEmbedding ${OE_HEADERS} ${OE_SOURCES})
target_link_libraries(OrbifoldEmbedding Mesh Utilities)


file(GLOB VIEWER_HEADERS Viewer/*.h)
file(GLOB VIEWER_SOURCES Viewer/*.cpp)
add_executable(Viewer ${VIEWER_HEADERS} ${VIEWER_SOURCES})
target_link_libraries(Viewer Mesh Utilities OrbifoldEmbedding BoundaryFirstFlattening igl::png igl::core igl::opengl igl::glfw igl::imgui OpenMeshCore OpenMeshTools)