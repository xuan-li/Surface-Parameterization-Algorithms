#ifndef MESH_SLICER
#define MESH_SLICER

#include <MeshDefinition.h>
#include "Dijkstra.h"

class MeshSlicer {
public:
	MeshSlicer(
		std::function<void(OpenMesh::VertexHandle, OpenMesh::VertexHandle)> f1,
		std::function<void(OpenMesh::FaceHandle, OpenMesh::FaceHandle)> f2,
		std::function<void(OpenMesh::EdgeHandle, OpenMesh::EdgeHandle)> f3 = [](OpenMesh::EdgeHandle e1, OpenMesh::EdgeHandle e2) {},
		std::function<void(OpenMesh::HalfedgeHandle, OpenMesh::HalfedgeHandle)> f4 = [](OpenMesh::HalfedgeHandle h1, OpenMesh::HalfedgeHandle h2) {} );
	
	void SliceMeshToDisk(SurfaceMesh &src, SurfaceMesh &dst);

protected:
	std::function<void(OpenMesh::VertexHandle, OpenMesh::VertexHandle)> TransferVertexData;
	std::function<void(OpenMesh::FaceHandle, OpenMesh::FaceHandle)> TransferFaceData;
	std::function<void(OpenMesh::EdgeHandle, OpenMesh::EdgeHandle)> TransferEdgeData;
	std::function<void(OpenMesh::HalfedgeHandle, OpenMesh::HalfedgeHandle)> TransferHalfedgeData;

	OpenMesh::VPropHandleT<int> wedge;

	void FindCutGraphSphere();
	void FindCutGraphNonSphere();
	void ConstructWedge();
	
};


#endif // !MESH_SLICER
