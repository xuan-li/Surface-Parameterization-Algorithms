#ifndef BOUNDARY_FIRST_FLATTENING
#define BOUNDARY_FIRST_FLATTENING

#include <MeshDefinition.h>
#include <Topology\MeshSlicer.h>
#include <Eigen\Sparse>
#include <Eigen\Dense>

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

	Eigen::VectorXd target_k_;
	
protected:
	void Init();
	void InitType1();
	void Slice();

	double CosineLaw(double a, double b, double c);
	void ComputeCornerAngles(SurfaceMesh &mesh);
	void ComputeVertexCurvatures(SurfaceMesh &mesh);
	
	void ComputeLaplacian(SurfaceMesh &mesh, bool mode = 0);
	void ComputeConformalFactors();


	void ComputeHalfedgeWeights(SurfaceMesh &mesh);

	void ReindexVertices(SurfaceMesh &mesh);

	void BoundaryUToTargetK();

	void IntegrateBoundaryCurve();
	
	void ExtendToInterior();

	void ComputeHarmonicMatrix();

	void NormalizeUV();
	
};

#endif // !BOUNDARY_FIRST_FLATTENING
