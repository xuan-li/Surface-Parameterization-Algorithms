#ifndef HYPERBOLIC_ORBIFOLD_SOLVER_H_
#define HYPERBOLIC_ORBIFOLD_SOLVER_H_

#include <MeshDefinition.h>
#include "OrbifoldMeshSlicer.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <Geometry/HyperbolicGeometry.h>

#include <LBFGS.h>

#ifndef PI
#define PI 3.141592653
#endif

class HyperbolicOrbifoldSolver
{
public:
	HyperbolicOrbifoldSolver(SurfaceMesh &mesh);
	SurfaceMesh Compute();

protected:
	SurfaceMesh &mesh_;
	SurfaceMesh sliced_mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vts_;
	std::vector<std::vector<OpenMesh::VertexHandle>> segments_vts_;
	OpenMesh::VPropHandleT<std::function<Complex(Complex const)>> vtx_transit_;
	
	int n_cones_;

protected:
	void CutToDist(int n_cones);

	void InitOrbifold();
	void InitType1();

	double CosineLaw(double a, double b, double c);
	double AngleCosineLaw(double a, double b, double c);
	void ComputeCornerAngles();
	void ComputeHalfedgeWeights();
	void InitMap();

	void ComputeEdgeLength();
	double ComputeGradient();
	OpenMesh::Vec2d ComputeGradientOfDistance2(Complex src, Complex dst);
	OpenMesh::Vec2d ComputeGradient(OpenMesh::VertexHandle v);
	OpenMesh::Vec2d ComputeIntrinsicGradient(OpenMesh::VertexHandle v);
	
	double ComputeEnergy();
	Eigen::VectorXd GetCoordsVector();
	Eigen::VectorXd GetGradientVector();
	void SetCoords(const Eigen::VectorXd &uv_vector);
	double OptimizationLoop(double step_length, double error);
	
	void Normalize();
};

#endif // !HYPERBOLIC_ORBIFOLD_SOLVER_H_
