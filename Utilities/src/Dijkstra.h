#ifndef DIJKSTRA
#define DIJKSTRA

#include <MeshDefinition.h>

void DijkstraShortestDist(SurfaceMesh *m_mesh, OpenMesh::VertexHandle src, OpenMesh::VPropHandleT<double> &dist, OpenMesh::VPropHandleT<OpenMesh::VertexHandle> &parent);


#endif