#ifndef BOUNDARY_FIRST_FLATTENING
#define BOUNDARY_FIRST_FLATTENING

#include <MeshDefinition.h>
#include <MeshSlicer.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "BFFInitializer.h"

#include <igl/active_set.h>

#ifndef PI
#define PI 3.141592653
#endif

// This class is the implementation of paper Boundary First Flattening.
// 
// BFF algorithm parameterize the surface according to boundary data.
// We need one of two kinds of boundary as input: conformal factors or target geodesic curvature.
// These two kinds of data can be converted to each other.
// Known boundary data is determined by user and processed in the initializer for BFF.
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

	// Cut the mesh into disk, and set all kinds of data and flags.
	void Init();

	double CosineLaw(double a, double b, double c);
	
	// Compute mesh data
	void ComputeCornerAngles(SurfaceMesh &mesh, Eigen::VectorXd l = Eigen::VectorXd());
	void ComputeHalfedgeWeights(SurfaceMesh &mesh);
	void ComputeVertexCurvatures(SurfaceMesh &mesh, Eigen::VectorXd l = Eigen::VectorXd());

	// Compute cotangent Laplacian operator.
	void ComputeLaplacian(SurfaceMesh &mesh, bool mode = false);
	
	// Seperate inner vertices and boundary vertices.
	void ReindexVertices(SurfaceMesh &mesh);

	// Convert boundary target geodesic curvature into conformal factors.
	void BoundaryTargetKKnown();

	// To minimize area distorsion, we preserve the length of boundary, i.e. u_B = 0.
	// And then we convert conformal factors u_B to target curvature.
	void FreeBoundary();

	// To compute a global parameterization with fixed cone angles.
	void GlobalParameterization();


	// The operator that convert boundary conformal factors to target curvatures.
	Eigen::VectorXd BoundaryUToTargetK(Eigen::VectorXd &u);

	// The operator that convert boundary target curvatures to conformal factors.
	Eigen::VectorXd BoundaryTargetKToU(Eigen::VectorXd &k);

	// Integrate boundary data into a closed loop.
	void IntegrateBoundaryCurve();
	
	// Given boundary's embedding, we use harmonic map to get one component.
	// And minimize conformal energy use hilbert transform over the other component.
	void ExtendToInteriorHilbert();

	// Use harmonic map on both components.
	void ExtendToInteriorHarmonic();

	void ComputeHarmonicMatrix();

	// Normalize uvs s.t. they all fall in unit circle.
	void NormalizeUV();
};

#endif // !BOUNDARY_FIRST_FLATTENING
