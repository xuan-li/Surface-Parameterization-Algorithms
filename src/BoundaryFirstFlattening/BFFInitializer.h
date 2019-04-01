#ifndef BFF_INITIALIZER_H_
#define BFF_INITIALIZER_H_

#include <MeshDefinition.h>
#include <MeshDefinition.h>
#include <MeshSlicer.h>
#include <map>
#include <Eigen/Core>

#ifndef PI
#define PI 3.141592653
#endif

// This class is a warper of class MeshSlicer to serve for BFF algorithm.
// The first step is to cut the mesh to disk.
// And then it will set singularity flags and target curvatures according to user's input.
// The constructor receive three properties: 
//		cone_flag: whether a vertex is a cone point ,
//		cone_angle: the target sum angles around cones,
//		slice_flag: whether a edge is on slice.
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

	//this property stores the each vertex is splitted to what vertices.
	OpenMesh::VPropHandleT<std::vector<OpenMesh::VertexHandle>> split_to_;

protected:
	// Cut the mesh into a disk.
	void CutMesh(SurfaceMesh &sliced_mesh);

};


#endif