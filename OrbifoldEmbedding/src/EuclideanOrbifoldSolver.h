#ifndef EUCLIDEAN_ORBIFOLD_H_
#define EUCLIDEAN_ORBIFOLD_H_

#include <MeshDefinition.h>
#include "OrbifoldMeshSlicer.h"
#include "TransformationBuilder.h"
#include <Eigen/Core>


class EuclideanOrbifoldSolver
{
public:
	EuclideanOrbifoldSolver(SurfaceMesh &mesh);
	SurfaceMesh Compute();
protected:
	SurfaceMesh &mesh_;
	SurfaceMesh sliced_mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vts;
	std::vector<std::vector<OpenMesh::VertexHandle>> segments_vts_;
	std::map<int, Eigen::Matrix2d> transforms_;

protected:
	void CutToDist();
	
	void InitOrbifold();
	void InitType1();

	double CosineLaw(double a, double b, double c);
	void ComputeCornerAngles();
	void ComputeHalfedgeWeights();


};

#endif
