#ifndef HYPERBOLIC_ORBIFOLD_SOLVER_H_
#define HYPERBOLIC_ORBIFOLD_SOLVER_H_

#include <MeshDefinition.h>
#include "OrbifoldInitializer.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <Geometry/HyperbolicGeometry.h>

#include <LBFGS.h>

#ifndef PI
#define PI 3.141592653
#endif


// This is the implementation of paper Hyperbolic Orbifold Embeddings.
// The problem is non linear.
class HyperbolicOrbifoldSolver
{
public:
	HyperbolicOrbifoldSolver(SurfaceMesh &mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::EPropHandleT<bool> slice_flag);
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
	OpenMesh::VPropHandleT<std::function<Complex(Complex const)>> vtx_transit_;
	
	int n_cones_;

	double max_error = 1e-4;

protected:

	void InitOrbifold();
	
	double CosineLaw(double a, double b, double c);
	double AngleCosineLaw(double a, double b, double c);
	void ComputeCornerAngles();
	void ComputeHalfedgeWeights();
	
	void InitiateBoundaryData();
	void InitMap();

	void ComputeEdgeLength();
	double ComputeGradient();
	OpenMesh::Vec2d ComputeGradientOfDistance2(Complex src, Complex dst);
	OpenMesh::Vec2d ComputeGradient(OpenMesh::VertexHandle v);
	
	double ComputeEnergy();
	Eigen::VectorXd GetCoordsVector();
	Eigen::VectorXd GetGradientVector();
	void SetCoords(const Eigen::VectorXd &uv_vector);
	double OptimizationLoop(double step_length, double error);
	
	// Normalize boundary such that it satisfies orbifold requirement.
	void Normalize();
};

#endif // !HYPERBOLIC_ORBIFOLD_SOLVER_H_
