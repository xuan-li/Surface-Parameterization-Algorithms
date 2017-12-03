#ifndef MESH_MARKER_H_
#define MESH_MARKER_H_

#include <MeshDefinition.h>
#include "Dijkstra.h"
#include <Eigen\Core>
#include "StringParser.h"

#ifndef PI
#define PI 3.14159254
#endif // !PI

// This class is to store user-defined cones, cone angles and slices.
class MeshMarker {
public:
	void SetObject(SurfaceMesh &mesh) { p_mesh_ = &mesh; };
	void ResetMarker();
	void ComputeAndSetSlice(OpenMesh::VertexHandle v0, OpenMesh::VertexHandle v1);
	void SetSingularity(OpenMesh::VertexHandle v, double cone_angle);
	OpenMesh::VPropHandleT<bool> GetSingularityFlag() { return singularity_; }
	OpenMesh::EPropHandleT<bool> GetSliceFlag() { return slice_; }
	OpenMesh::VPropHandleT<double> GetConeAngleFlag() { return cone_angle_; }
	void GenerateMatrix(Eigen::MatrixXd &P, Eigen::MatrixXd &EP1, Eigen::MatrixXd &EP2);
	void LoadFromFile(std::string filename);
	void SaveToFile(std::string filename);
protected:
	SurfaceMesh *p_mesh_;
	OpenMesh::VPropHandleT<bool> singularity_;
	OpenMesh::VPropHandleT<double> cone_angle_;
	OpenMesh::EPropHandleT<bool> slice_;
	int n_edges_;
	int n_vertices_;
};

#endif // !MESH_MARKER_H_
