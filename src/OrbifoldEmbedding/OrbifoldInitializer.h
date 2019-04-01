#ifndef ORBIFOLD_INITIALIZER_H
#define ORBIFOLD_INITIALIZER_H

#include <MeshDefinition.h>
#include <MeshSlicer.h>
#include <map>
#include <Eigen/Core>
#include <EuclideanGeometry2D.h>
#include <HyperbolicGeometry.h>

#ifndef PI
#define PI 3.141592654
#endif // !PI


typedef std::complex<double> Complex;


// This class is a warper of MeshSlicer to serve for Orbifold algorithms.
// The first step is to cut the mesh to disk.
// And then it will set singularity flags and target curvatures according to user's input.
// The process receives three properties: 
//		cone_flag: whether a vertex is a cone point ,
//		cone_angle: the target sum angles around cones,
//		slice_flag: whether a edge is on slice.
class OrbifoldInitializer {

public:
	OrbifoldInitializer(SurfaceMesh &mesh);
	void Initiate(SurfaceMesh &sliced_mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag);
	
	// Compute Isometries needed for orbifold requirements.
	// Euclidean isometries are stored using homogenous matrx.
	// Hyperbolic isometries are stored using functions.
	void ComputeEuclideanTransformations(SurfaceMesh &sliced_mesh, OpenMesh::VPropHandleT<Eigen::Matrix3d> &vtx_transit);
	void ComputeHyperbolicTransformations(SurfaceMesh &sliced_mesh, OpenMesh::VPropHandleT<std::function<Complex(Complex const)>> &vtx_transit);

	std::vector<OpenMesh::VertexHandle> GetConeVertices() { return cone_vertices_; }

	std::vector<std::vector<OpenMesh::VertexHandle>> GetSegments() { return segments_vts_; }
	
protected:
	SurfaceMesh &mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vertices_;
	std::vector<std::vector<OpenMesh::VertexHandle>> segments_vts_;

	OpenMesh::VPropHandleT<bool> cone_flag_;
	OpenMesh::VPropHandleT<double> cone_angle_;
	OpenMesh::EPropHandleT<bool> slice_flag_;

protected:
	// Cut the mesh into a disk.
	void CutMesh(SurfaceMesh &sliced_mesh);

	// Cut the boundary into segments according to cones.
	void CutBoundaryToSegments(SurfaceMesh &sliced_mesh);

	// Initiate the coordinates of cones.
	void InitiateEConeCoords(SurfaceMesh &sliced_mesh);
	void InitiateEConeCoordsType1(SurfaceMesh &sliced_mesh);
	void InitiateEConeCoordsType2(SurfaceMesh &sliced_mesh);
	void InitiateEConeCoordsType3(SurfaceMesh &sliced_mesh, bool up );
	void InitiateHConeCoords(SurfaceMesh &sliced_mesh);

	double HyperbolicCosineLaw(double a, double b, double c);
	
};
#endif // !ORBIFOLD_MESH_SLICER_H

