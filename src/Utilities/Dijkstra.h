#ifndef DIJKSTRA
#define DIJKSTRA

#include <MeshDefinition.h>
#include <list>

#ifndef INF
#define INF 0x3f3f3f3f
#endif

// This function is modified from my previous research codes.
void DijkstraShortestDist(
	SurfaceMesh &mesh, 
	OpenMesh::VertexHandle src, // root node
	OpenMesh::VPropHandleT<double> &dist, // property to store distance from the root
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> &parent //property to store parent
); 


#endif