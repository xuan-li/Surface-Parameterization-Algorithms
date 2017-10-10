#include "MeshSlicer.h"



MeshSlicer::MeshSlicer(std::function<void(OpenMesh::VertexHandle, OpenMesh::VertexHandle)> f1, std::function<void(OpenMesh::FaceHandle, OpenMesh::FaceHandle)> f2, std::function<void(OpenMesh::EdgeHandle, OpenMesh::EdgeHandle)> f3, std::function<void(OpenMesh::HalfedgeHandle, OpenMesh::HalfedgeHandle)> f4)
{
	TransferVertexData = f1;
	TransferFaceData = f2;
	TransferEdgeData = f3;
	TransferHalfedgeData = f4;
}

void MeshSlicer::SliceMeshToDisk(SurfaceMesh & src, SurfaceMesh & dst)
{
}

void MeshSlicer::FindCutGraphSphere()
{
}

void MeshSlicer::FindCutGraphNonSphere()
{

}

void MeshSlicer::ConstructWedge()
{

}
