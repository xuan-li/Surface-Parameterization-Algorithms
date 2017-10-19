#ifndef EUCLIDEAN_ORBIFOLD_H_
#define EUCLIDEAN_ORBIFOLD_H_

#include <MeshDefinition.h>
#include "OrbifoldMeshSlicer.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

#ifndef PI
#define PI 3.141592653
#endif



class EuclideanOrbifoldSolver
{
public:
	EuclideanOrbifoldSolver(SurfaceMesh &mesh);
	SurfaceMesh Compute();
protected:
	SurfaceMesh &mesh_;
	SurfaceMesh sliced_mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vts_;
	std::vector<std::vector<OpenMesh::VertexHandle>> segments_vts_;
	OpenMesh::VPropHandleT<Eigen::Matrix2d> vtx_transit_;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> vtx_rotation_center_;
	Eigen::SparseMatrix<double> A_;
	Eigen::VectorXd b_;
	Eigen::VectorXd X_;
	

protected:
	void CutToDist();
	
	void InitOrbifold();
	void InitType1();

	double CosineLaw(double a, double b, double c);
	void ComputeCornerAngles();
	void ComputeHalfedgeWeights();

	void ConstructSparseSystem();
	void SolveLinearSystem();
	
	

};

#endif
