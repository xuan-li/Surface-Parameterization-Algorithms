#ifndef HYPERBOLIC_ORBIFOLD_SOLVER_H_
#define HYPERBOLIC_ORBIFOLD_SOLVER_H_

#include <MeshDefinition.h>
#include "OrbifoldMeshSlicer.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <Geometry/HyperbolicGeometry.h>

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
	void ComputeGradient();
	
	
};

#endif // !HYPERBOLIC_ORBIFOLD_SOLVER_H_
