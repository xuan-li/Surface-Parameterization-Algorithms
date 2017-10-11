#ifndef DIJKSTRA
#define DIJKSTRA

#include <MeshDefinition.h>

void DijkstraShortestDist(
	SurfaceMesh &mesh, 
	OpenMesh::VertexHandle src, // root node
	OpenMesh::VPropHandleT<double> &dist, // property to store distance from the root
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> &parent //property to store parent
); 


#endif