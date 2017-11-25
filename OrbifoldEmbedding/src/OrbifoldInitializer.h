#ifndef ORBIFOLD_INITIALIZER_H
#define ORBIFOLD_INITIALIZER_H

#include <MeshDefinition.h>
#include <Topology/MeshSlicer.h>
#include <map>
#include <Eigen\Core>
#include <Geometry\EuclideanGeometry2D.h>
#include <Geometry\HyperbolicGeometry.h>

#ifndef PI
#define PI 3.141592654
#endif // !PI


typedef std::complex<double> Complex;

class OrbifoldInitializer {

public:
	OrbifoldInitializer(SurfaceMesh &mesh);
	void Initiate(SurfaceMesh &sliced_mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag);
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
	void CutMesh(SurfaceMesh &sliced_mesh);
	void CutBoundaryToSegments(SurfaceMesh &sliced_mesh);

	void InitiateEConeCoords(SurfaceMesh &sliced_mesh);
	void InitiateEConeCoordsType1(SurfaceMesh &sliced_mesh);
	void InitiateEConeCoordsType2(SurfaceMesh &sliced_mesh);
	void InitiateEConeCoordsType3(SurfaceMesh &sliced_mesh, bool up );

	void InitiateHConeCoords(SurfaceMesh &sliced_mesh);
	double HyperbolicCosineLaw(double a, double b, double c);
	
};
#endif // !ORBIFOLD_MESH_SLICER_H

