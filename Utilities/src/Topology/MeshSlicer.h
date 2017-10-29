#ifndef MESH_SLICER
#define MESH_SLICER

#include "MeshDefinition.h"
#include "..\Dijkstra.h"

class MeshSlicer {
public:
	MeshSlicer(
		SurfaceMesh &mesh,
		std::function<void(OpenMesh::VertexHandle, OpenMesh::VertexHandle)> f1 = [](OpenMesh::VertexHandle v1, OpenMesh::VertexHandle v2) {},
		std::function<void(OpenMesh::FaceHandle, OpenMesh::FaceHandle)> f2 = [](OpenMesh::FaceHandle f1, OpenMesh::FaceHandle f2) {},
		std::function<void(OpenMesh::EdgeHandle, OpenMesh::EdgeHandle)> f3 = [](OpenMesh::EdgeHandle e1, OpenMesh::EdgeHandle e2) {},
		std::function<void(OpenMesh::HalfedgeHandle, OpenMesh::HalfedgeHandle)> f4 = [](OpenMesh::HalfedgeHandle h1, OpenMesh::HalfedgeHandle h2) {} );
	
	SurfaceMesh SliceMeshToDisk();
	std::vector<OpenMesh::VertexHandle> SplitTo(OpenMesh::VertexHandle v);
	std::vector<OpenMesh::VertexHandle> GetLongestPath();

protected:
	std::function<void(OpenMesh::VertexHandle, OpenMesh::VertexHandle)> TransferVertexData;
	std::function<void(OpenMesh::FaceHandle, OpenMesh::FaceHandle)> TransferFaceData;
	std::function<void(OpenMesh::EdgeHandle, OpenMesh::EdgeHandle)> TransferEdgeData;
	std::function<void(OpenMesh::HalfedgeHandle, OpenMesh::HalfedgeHandle)> TransferHalfedgeData;

	SurfaceMesh &mesh_;
	OpenMesh::VertexHandle base_point_;
	OpenMesh::HPropHandleT<int> wedge_;
	OpenMesh::EPropHandleT<bool> on_cut_;
	OpenMesh::VPropHandleT<std::vector<OpenMesh::VertexHandle>> split_to_;
	std::vector<OpenMesh::VertexHandle> longest_path_;

	void FindAndMarkCutGraphSphere();
	void FindAndMarkCutGraphNonSphere();
	void ConstructWedge();
	SurfaceMesh SliceAccordingToWedge();
};


#endif // !MESH_SLICER
