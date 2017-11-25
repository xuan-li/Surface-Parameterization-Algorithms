#ifndef BFF_INITIALIZER_H_
#define BFF_INITIALIZER_H_

#include <MeshDefinition.h>
#include <MeshDefinition.h>
#include <Topology/MeshSlicer.h>
#include <map>
#include <Eigen\Core>

class BFFInitializer {
public:
	BFFInitializer(SurfaceMesh &mesh);
	void Initiate(SurfaceMesh &sliced_mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag);
	std::vector<OpenMesh::VertexHandle> GetConeVertices() { return cone_vertices_; }
	OpenMesh::VPropHandleT<std::vector<OpenMesh::VertexHandle>> split_to() { return split_to_; }
protected:
	SurfaceMesh &mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vertices_;
	OpenMesh::VPropHandleT<bool> cone_flag_;
	OpenMesh::VPropHandleT<double> cone_angle_;
	OpenMesh::EPropHandleT<bool> slice_flag_;

	OpenMesh::VPropHandleT<std::vector<OpenMesh::VertexHandle>> split_to_;

protected:
	void CutMesh(SurfaceMesh &sliced_mesh);

};


#endif