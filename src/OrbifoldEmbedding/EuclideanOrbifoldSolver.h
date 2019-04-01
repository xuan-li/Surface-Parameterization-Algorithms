#ifndef EUCLIDEAN_ORBIFOLD_H_
#define EUCLIDEAN_ORBIFOLD_H_

#include <MeshDefinition.h>
#include "OrbifoldInitializer.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#ifndef PI
#define PI 3.141592653
#endif


// This class is the implementation of paper: Orbifold Tutte Embeddings.
// The system is harmonic system plus rotation constraints.
// So the problem is linear.
class EuclideanOrbifoldSolver
{
public:
	EuclideanOrbifoldSolver(SurfaceMesh &mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag);
	SurfaceMesh Compute();
	std::vector<OpenMesh::VertexHandle> ConeVertices() { return cone_vts_; }
protected:
	SurfaceMesh &mesh_;
	SurfaceMesh sliced_mesh_;

	OpenMesh::VPropHandleT<bool> cone_flag_;
	OpenMesh::VPropHandleT<double> cone_angle_;
	OpenMesh::EPropHandleT<bool> slice_flag_;

	std::vector<OpenMesh::VertexHandle> cone_vts_;
	std::vector<std::vector<OpenMesh::VertexHandle>> segments_vts_;
	OpenMesh::VPropHandleT<Eigen::Matrix3d> vtx_transit_;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> vtx_rotation_center_;
	Eigen::SparseMatrix<double> A_;
	Eigen::VectorXd b_;
	Eigen::VectorXd X_;
	

protected:

	void InitOrbifold();


	double CosineLaw(double a, double b, double c);
	void ComputeCornerAngles();
	void ComputeHalfedgeWeights();

	void ConstructSparseSystem();
	void SolveLinearSystem();
	
	

};

#endif
