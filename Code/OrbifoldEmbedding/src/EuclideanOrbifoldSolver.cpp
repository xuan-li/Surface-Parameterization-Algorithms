#include "EuclideanOrbifoldSolver.h"
#include <list>

EuclideanOrbifoldSolver::EuclideanOrbifoldSolver(SurfaceMesh & mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag)
	:mesh_(mesh), cone_flag_(cone_flag), cone_angle_(cone_angle), slice_flag_(slice_flag)
{
	
}

SurfaceMesh EuclideanOrbifoldSolver::Compute()
{
	if (mesh_.n_vertices() > 10) {
		InitOrbifold();
		ComputeHalfedgeWeights();
		ConstructSparseSystem();
		SolveLinearSystem();
	}
	return sliced_mesh_;
}


void EuclideanOrbifoldSolver::InitOrbifold()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	OrbifoldInitializer initializer(mesh_);
	initializer.Initiate(mesh,cone_flag_, cone_angle_, slice_flag_);
	initializer.ComputeEuclideanTransformations(sliced_mesh_, vtx_transit_);
	cone_vts_ = initializer.GetConeVertices();
	segments_vts_ = initializer.GetSegments();

	std::cout << "Cone coordinates:\n";
	for (int i = 0; i < cone_vts_.size(); ++i) {
		Vec2d uv = mesh.texcoord2D(cone_vts_[i]);
		std::cout << uv[0] << "\t" << uv[1] << std::endl;
	}
	
}


double EuclideanOrbifoldSolver::CosineLaw(double a, double b, double c)
{
	double cs = (a * a + b * b - c * c) / (2 * a * b);
	assert(-1 <= cs && cs <= 1);
	return acos(cs);
}

void EuclideanOrbifoldSolver::ComputeCornerAngles()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); ++fiter) {
		FaceHandle f = *fiter;
		std::vector<HalfedgeHandle> he;

		for (auto fhiter = mesh.fh_iter(f); fhiter.is_valid(); ++fhiter) {
			he.push_back(*fhiter);
		}
		double l[3];
		l[0] = mesh.calc_edge_length(mesh.edge_handle(he[0]));
		l[1] = mesh.calc_edge_length(mesh.edge_handle(he[1]));
		l[2] = mesh.calc_edge_length(mesh.edge_handle(he[2]));
		for (int i = 0; i < 3; ++i) {
			double cs = CosineLaw(l[i], l[(i + 1) % 3], l[(i + 2) % 3]);
			mesh.data(he[i]).set_angle(cs);
		}
	}
}

void EuclideanOrbifoldSolver::ComputeHalfedgeWeights()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	ComputeCornerAngles();
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		EdgeHandle e = *eiter;
		HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
		HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);
		HalfedgeHandle h0_next;
		if(!mesh.is_boundary(h0))
			h0_next = mesh.next_halfedge_handle(h0);
		HalfedgeHandle h1_next;
		if(!mesh.is_boundary(h1))
			h1_next = mesh.next_halfedge_handle(h1);
		
		double weight = 0.;

		if (h0_next.is_valid())
			weight += 1. / tan(mesh.data(h0_next).angle());
		if (h1_next.is_valid())
			weight += 1. / tan(mesh.data(h1_next).angle());
		weight *= 0.5;
		if (weight < 0)
			weight = 0.01;
		mesh.data(h0).set_weight(weight);
		mesh.data(h1).set_weight(weight);
	}

}

void EuclideanOrbifoldSolver::ConstructSparseSystem()
{
	using namespace OpenMesh;
	using namespace Eigen;
	SurfaceMesh &mesh = sliced_mesh_;

	A_.resize(2 * mesh.n_vertices(), 2 * mesh.n_vertices());
	A_.setZero();
	b_.resize(2 * mesh.n_vertices());
	std::vector<Eigen::Triplet<double> > A_coefficients;

	// interior vertex satisfies normal harmonic condition
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh.data(v).is_singularity()) {
			auto uv = mesh.texcoord2D(v);
			b_(2 * v.idx()) = uv[0];
			b_(2 * v.idx() + 1) = uv[1];
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * v.idx(), 1.));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * v.idx() + 1, 1.));
		}
		else if (!mesh.is_boundary(v)) {
			b_(2 * v.idx()) = 0;
			b_(2 * v.idx() + 1) = 0;
			double s_w = 0;
			for (SurfaceMesh::VertexOHalfedgeIter vohiter = mesh.voh_iter(v); vohiter.is_valid(); ++vohiter) {
				HalfedgeHandle h = *vohiter;
				VertexHandle neighbor = mesh.to_vertex_handle(h);
				double n_w = mesh.data(h).weight();
				s_w += n_w;
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * neighbor.idx(), -n_w));
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * neighbor.idx() + 1, -n_w));
			}
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * v.idx(), s_w));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * v.idx() + 1, s_w));
		}
	}


	// handle boundary vts;
	for (int i = 0; i < segments_vts_.size() / 2; ++i) {
		for (auto it = segments_vts_[i].begin(); it != segments_vts_[i].end(); ++it) {
			VertexHandle v = *it;
			if (mesh.data(v).is_singularity()) continue;
			b_(2 * v.idx()) = 0;
			b_(2 * v.idx() + 1) = 0;
			double s_w = 0.;
			for (auto vviter = mesh.vv_iter(v); vviter.is_valid(); ++vviter) {
				VertexHandle neighbor = *vviter;
				HalfedgeHandle h = mesh.find_halfedge(v, neighbor);
				double n_w = mesh.data(h).weight();
				s_w += n_w;
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * neighbor.idx(), -n_w));
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * neighbor.idx() + 1, -n_w));
			}
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * v.idx(), s_w));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * v.idx() + 1, s_w));


			auto equiv = mesh.data(v).equivalent_vertex();
			Matrix2d coeff_equiv;
			coeff_equiv.setZero();
			Matrix3d T = mesh.property(vtx_transit_, equiv); // from equiv to v
			Matrix2d rotation_matrix = T.block(0,0,2,2);
			//std::cout << rotation_matrix << std::endl;
			for (auto vviter = mesh.vv_iter(equiv); vviter.is_valid(); ++vviter) {
				VertexHandle neighbor = *vviter;
				HalfedgeHandle h = mesh.find_halfedge(equiv, neighbor);
				double n_w = mesh.data(h).weight();
				auto coeff_equiv_neighbor = n_w * rotation_matrix;
				coeff_equiv += coeff_equiv_neighbor;
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * neighbor.idx(), -coeff_equiv_neighbor(0,0)));
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * neighbor.idx() + 1, -coeff_equiv_neighbor(0,1)));
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * neighbor.idx(), -coeff_equiv_neighbor(1,0)));
				A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * neighbor.idx() + 1, -coeff_equiv_neighbor(1,1)));
			}
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * equiv.idx(), coeff_equiv(0, 0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx(), 2 * equiv.idx() + 1, coeff_equiv(0, 1)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * equiv.idx(), coeff_equiv(1, 0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * v.idx() + 1, 2 * equiv.idx() + 1, coeff_equiv(1, 1)));

			
			b_.segment(2 * equiv.idx(), 2) = - T.block(0,2, 2, 1);
					
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx(), 2 * equiv.idx(), rotation_matrix(0,0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx(), 2 * equiv.idx() + 1, rotation_matrix(0, 1)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx() + 1, 2 * equiv.idx(), rotation_matrix(1, 0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx() + 1, 2 * equiv.idx()+1, rotation_matrix(1, 1)));

			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx(), 2 * v.idx(), -1.));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx() + 1, 2 * v.idx() + 1, - 1.));
			
		}
	}

	A_.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
}

void EuclideanOrbifoldSolver::SolveLinearSystem()
{
	using namespace OpenMesh;
	using namespace Eigen;
	SurfaceMesh &mesh = sliced_mesh_;

	//Eigen::SparseQR<Eigen::SparseMatrix<double>, COLAMDOrdering<int>> solver;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A_);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::VectorXd x = solver.solve(b_);
	std::cout << "Error:" << (A_ * x - b_).norm() << std::endl;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec2d uv(x(2 * v.idx()), x(2 * v.idx() + 1));
		mesh.set_texcoord2D(v, uv);
	}
	/*Eigen::MatrixXd B = A_.toDense();
	int n_boundary = mesh.GetBoundaries().front().size();
	Eigen::MatrixXd C(n_boundary * 2, mesh.n_vertices() * 2);
	int i = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh.is_boundary(v)) {
			C.block(2 * i, 0, 2, mesh.n_vertices() * 2) = B.block(2 * v.idx(), 0, 2, mesh.n_vertices() * 2);
			++i;
		}
	}
	std::cout << C.transpose() << std::endl;*/
}
