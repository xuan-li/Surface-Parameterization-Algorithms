#include "EuclideanOrbifoldSolver.h"
#include <list>

EuclideanOrbifoldSolver::EuclideanOrbifoldSolver(SurfaceMesh & mesh)
	:mesh_(mesh)
{

}

SurfaceMesh EuclideanOrbifoldSolver::Compute()
{
	InitOrbifold();
	ComputeHalfedgeWeights();
	ConstructSparseSystem();
	SolveLinearSystem();
	return sliced_mesh_;
}

void EuclideanOrbifoldSolver::CutToDist()
{
	OrbifoldMeshSlicer slicer(mesh_);
	sliced_mesh_ = slicer.CutAndSelectSingularities();
	cone_vts_.push_back(slicer.GetConeVertices().front());

}


void EuclideanOrbifoldSolver::InitOrbifold()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	cone_vts_.clear();
	CutToDist();
	cone_vts_.erase(++cone_vts_.begin(), cone_vts_.end());
	mesh.RequestBoundary();
	auto boundary = mesh.GetBoundaries().front();
	std::list<HalfedgeHandle> boundary_list(boundary.begin(), boundary.end());
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		HalfedgeHandle h = *it;
		if (mesh.from_vertex_handle(h) == cone_vts_[0]) {
			boundary_list.insert(boundary_list.end(), boundary_list.begin(), it);
			boundary_list.erase(boundary_list.begin(), it);
			break;
		}
		
	}

	segments_vts_.clear();
	std::vector<OpenMesh::VertexHandle> segment;
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		HalfedgeHandle h = *it;
		
		if (mesh.data(mesh.to_vertex_handle(h)).is_singularity()) {
			cone_vts_.push_back(mesh.to_vertex_handle(h));
			segments_vts_.push_back(segment);
			segment.clear();
			continue;
		}
		segment.push_back(mesh.to_vertex_handle(h));
	}
	cone_vts_.pop_back();

	InitType1();
}

void EuclideanOrbifoldSolver::InitType1()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	OpenMesh::Vec2d v0(0, 0);
	OpenMesh::Vec2d v1(0, 1);
	OpenMesh::Vec2d v2(1, 1);
	OpenMesh::Vec2d v3(1, 0);
	sliced_mesh_.set_texcoord2D(cone_vts_[0], v0);
	sliced_mesh_.set_texcoord2D(cone_vts_[1], v1);
	sliced_mesh_.set_texcoord2D(cone_vts_[2], v2);
	sliced_mesh_.set_texcoord2D(cone_vts_[3], v3);
	double thetas[4] = {-PI /2, PI/2 , -PI/2, PI/2};
	Eigen::Matrix2d transforms[4];
	VertexHandle rotate_center[4] = { cone_vts_[0], cone_vts_[2] , cone_vts_[2] , cone_vts_[0] };
	if (!vtx_transit_.is_valid()) {
		mesh.add_property(vtx_transit_);
	}

	if (!vtx_rotation_center_.is_valid()) {
		mesh.add_property(vtx_rotation_center_);
	}

	for (int i = 0; i < 4; ++i) {
		transforms[i] << cos(thetas[i]), -sin(thetas[i]),
						 sin(thetas[i]), cos(thetas[i]);
	}
	
	for (int i = 0; i < segments_vts_.size(); ++i) {
		for (auto it = segments_vts_[i].begin(); it != segments_vts_[i].end(); ++it) {
			VertexHandle v = *it;
			mesh.property(vtx_rotation_center_, v) = rotate_center[i];
			mesh.property(vtx_transit_, v) = transforms[i];
		}
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
	for (int i = 0; i < 2; ++i) {
		for (auto it = segments_vts_[i].begin(); it != segments_vts_[i].end(); ++it) {
			VertexHandle v = *it;
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
			Matrix2d rotation_matrix = mesh.property(vtx_transit_, equiv);
			std::cout << rotation_matrix << std::endl;
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

			
			b_(2 * equiv.idx()) = 0;
			b_(2 * equiv.idx() + 1) = 0;
			int rotate_center_idx = mesh.property(vtx_rotation_center_, equiv).idx();
			
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx(), 2 * rotate_center_idx, 1. - rotation_matrix(0, 0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx(), 2 * rotate_center_idx + 1, -rotation_matrix(0, 1)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx() + 1, 2 * rotate_center_idx, -rotation_matrix(1, 0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx() + 1, 2 * rotate_center_idx + 1, 1. - rotation_matrix(1, 1)));

			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx(), 2 * equiv.idx(),  rotation_matrix(0, 0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx(), 2 * equiv.idx() + 1, rotation_matrix(0, 1)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx() + 1, 2 * equiv.idx(), rotation_matrix(1, 0)));
			A_coefficients.push_back(Eigen::Triplet<double>(2 * equiv.idx() + 1, 2 * equiv.idx() + 1, rotation_matrix(1, 1)));

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
	//std::cout << A_.diagonal()<<std::endl << std::endl;
	//std::cout << b_ << std::endl << std::endl;
	//std::cout << x << std::endl << std::endl;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec2d uv(x(2 * v.idx()), x(2 * v.idx() + 1));
		mesh.set_texcoord2D(v, uv);
	}
}
