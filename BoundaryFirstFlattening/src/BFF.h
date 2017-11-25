#ifndef BOUNDARY_FIRST_FLATTENING
#define BOUNDARY_FIRST_FLATTENING

#include <MeshDefinition.h>
#include <Topology\MeshSlicer.h>
#include <Eigen\Sparse>
#include <Eigen\Dense>
#include <Eigen\IterativeLinearSolvers>
#include "BFFInitializer.h"

#ifndef PI
#define PI 3.141592653
#endif

class BFFSolver {
public:
	BFFSolver(SurfaceMesh &mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag);
	SurfaceMesh Compute(int mode = 0);
	std::vector<OpenMesh::VertexHandle>  ConeVertices() { return cone_vts_; }

protected:
	SurfaceMesh &mesh_;
	SurfaceMesh sliced_mesh_;

	std::vector<OpenMesh::VertexHandle> cone_vts_;

	OpenMesh::VPropHandleT<bool> cone_flag_;
	OpenMesh::VPropHandleT<double> cone_angle_;
	OpenMesh::EPropHandleT<bool> slice_flag_;

	OpenMesh::VPropHandleT<std::vector<OpenMesh::VertexHandle>> split_to_;

	Eigen::SparseMatrix<double> Delta_;

	int n_boundary_;
	int n_interior_;
	int n_cones_;

protected:
	void Init();

	double CosineLaw(double a, double b, double c);
	
	// Compute mesh data
	void ComputeCornerAngles(SurfaceMesh &mesh, Eigen::VectorXd &l = Eigen::VectorXd());
	void ComputeHalfedgeWeights(SurfaceMesh &mesh);
	void ComputeVertexCurvatures(SurfaceMesh &mesh, Eigen::VectorXd &l = Eigen::VectorXd());
	
	void ComputeConformalFactors();

	void ComputeLaplacian(SurfaceMesh &mesh, bool mode = false);
	
	void ReindexVertices(SurfaceMesh &mesh);

	void BoundaryUToTargetK(bool free_boundary = false);

	Eigen::VectorXd BoundaryUToTargetK(Eigen::VectorXd &u);

	Eigen::VectorXd BoundaryTargetKToU(Eigen::VectorXd &k);

	void IntegrateBoundaryCurve();
	
	void ExtendToInteriorHilbert();

	void ExtendToInteriorHarmonic();

	void ComputeHarmonicMatrix();

	void NormalizeUV();
};

#endif // !BOUNDARY_FIRST_FLATTENING
