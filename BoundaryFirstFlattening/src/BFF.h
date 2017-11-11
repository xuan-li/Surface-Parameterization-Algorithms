#ifndef BOUNDARY_FIRST_FLATTENING
#define BOUNDARY_FIRST_FLATTENING

#include <MeshDefinition.h>
#include <Topology\MeshSlicer.h>
#include <Eigen\Sparse>
#include <Eigen\Dense>
#include <Eigen\IterativeLinearSolvers>

#ifndef PI
#define PI 3.141592653
#endif

class BFFSolver {
public:
	BFFSolver(SurfaceMesh &mesh);
	SurfaceMesh Compute();
	std::vector<OpenMesh::VertexHandle>  ConeVertices() { return cone_vts_; }

protected:
	SurfaceMesh &mesh_;
	SurfaceMesh sliced_mesh_;
	std::vector<OpenMesh::VertexHandle> slice_path_;
	std::vector<OpenMesh::VertexHandle> cone_vts_;
	std::vector<OpenMesh::VertexHandle> sliced_cone_vts_;
	MeshSlicer slicer_;
	OpenMesh::VPropHandleT<bool> on_slice_;


	Eigen::SparseMatrix<double> Delta_;

	int n_boundary_;
	int n_interior_;
	int n_cones_;

	Eigen::VectorXd target_k_;
	
protected:
	void Init();
	void InitType1();
	void Slice();

	double CosineLaw(double a, double b, double c);
	void ComputeCornerAngles(SurfaceMesh &mesh, Eigen::VectorXd &l = Eigen::VectorXd());
	void ComputeVertexCurvatures(SurfaceMesh &mesh, Eigen::VectorXd &l = Eigen::VectorXd());
	
	void ComputeLaplacian(SurfaceMesh &mesh, bool mode = false);
	void ComputeConformalFactors();


	void ComputeHalfedgeWeights(SurfaceMesh &mesh);

	void ReindexVertices(SurfaceMesh &mesh);

	void BoundaryUToTargetK();

	Eigen::VectorXd BoundaryUToTargetK(Eigen::VectorXd &u);

	Eigen::VectorXd BoundaryTargetKToU(Eigen::VectorXd &k);

	void IntegrateBoundaryCurve();
	
	void ExtendToInteriorHilbert();

	void ExtendToInteriorHarmonic();

	void ComputeHarmonicMatrix();

	void NormalizeUV();

	void ComputeOrbifoldBoundaryData();
	
};

#endif // !BOUNDARY_FIRST_FLATTENING
