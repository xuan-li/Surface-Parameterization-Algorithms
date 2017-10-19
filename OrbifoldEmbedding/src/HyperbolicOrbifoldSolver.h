#ifndef HYPERBOLIC_ORBIFOLD_SOLVER_H_
#define HYPERBOLIC_ORBIFOLD_SOLVER_H_

#include <MeshDefinition.h>
#include "OrbifoldMeshSlicer.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

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
	OpenMesh::VPropHandleT<Eigen::Matrix2d> vtx_transit_;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> vtx_rotation_center_;

protected:
	void CutToDist();

	void InitOrbifold();
	void InitType1();

	double CosineLaw(double a, double b, double c);
	void ComputeCornerAngles();
	void ComputeHalfedgeWeights();
	
};

#endif // !HYPERBOLIC_ORBIFOLD_SOLVER_H_
