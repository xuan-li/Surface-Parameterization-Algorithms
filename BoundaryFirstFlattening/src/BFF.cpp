#include "BFF.h"

BFFSolver::BFFSolver(SurfaceMesh & mesh)
	:mesh_(mesh), slicer_(mesh)
{

}

SurfaceMesh BFFSolver::Compute()
{
	Init();
	
	ComputeConformalFactors();

	Slice();

	ComputeVertexCurvatures(sliced_mesh_);

	//ComputeOrbifoldBoundaryData();
	
	BoundaryUToTargetK();

	IntegrateBoundaryCurve();

	//ExtendToInteriorHarmonic();
	ExtendToInteriorHilbert();

	NormalizeUV();

	return sliced_mesh_;
}

void BFFSolver::Init()
{
	using namespace OpenMesh;
	slicer_.FindAndMarkCutGraphSphere();
	slice_path_ = slicer_.GetLongestPath();
	if (!on_slice_.is_valid()) {
		mesh_.add_property(on_slice_);
	}
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		mesh_.property(on_slice_, *viter) = false;
	}
	for (auto it = slice_path_.begin(); it != slice_path_.end(); ++it) {
		mesh_.property(on_slice_, *it) = true;
	}
	InitType1();
}

void BFFSolver::InitType1()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = mesh_;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.data(v).set_target_curvature(0);
		mesh.data(v).set_singularity(false);
	}
	mesh.data(slice_path_[0]).set_target_curvature(PI);
	mesh.data(slice_path_[0]).set_singularity(true);
	mesh.data(slice_path_[0]).set_order(4);
	/*mesh.data(slice_path_[slice_path_.size() / 2]).set_target_curvature(PI);
	mesh.data(slice_path_[slice_path_.size() / 2]).set_singularity(true);
	mesh.data(slice_path_[slice_path_.size() / 2]).set_order(2);*/
	for (auto vviter = mesh.vv_iter(slice_path_[slice_path_.size() / 3]); vviter.is_valid(); ++vviter) {
		VertexHandle neighbor = *vviter;
		int center_index = slice_path_.size() / 3;
		if (neighbor != slice_path_[center_index + 1] && neighbor != slice_path_[center_index - 1]) {
			HalfedgeHandle h = mesh.find_halfedge(slice_path_[center_index], neighbor);
			EdgeHandle e = mesh.edge_handle(h);
			slicer_.AddOnCutEdge(e);
			mesh.data(neighbor).set_singularity(true);
			mesh.data(neighbor).set_target_curvature(PI);
			mesh.data(neighbor).set_order(2);
			break;
		}
	}
	for (auto vviter = mesh.vv_iter(slice_path_[2 * slice_path_.size() / 3]); vviter.is_valid(); ++vviter) {
		VertexHandle neighbor = *vviter;
		int center_index = 2 * slice_path_.size() / 3;
		if (neighbor != slice_path_[center_index + 1] && neighbor != slice_path_[center_index - 1]) {
			HalfedgeHandle h = mesh.find_halfedge(slice_path_[center_index], neighbor);
			EdgeHandle e = mesh.edge_handle(h);
			slicer_.AddOnCutEdge(e);
			mesh.data(neighbor).set_singularity(true);
			mesh.data(neighbor).set_target_curvature(PI);
			mesh.data(neighbor).set_order(2);
			break;
		}
	}
	mesh.data(slice_path_.back()).set_target_curvature(PI);
	mesh.data(slice_path_.back()).set_singularity(true);
	mesh.data(slice_path_.back()).set_order(4);

	
	
}

void BFFSolver::Slice()
{
	using namespace OpenMesh;
	slicer_.ConstructWedge();
	slicer_.SliceAccordingToWedge(sliced_mesh_);
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		auto equiv_vts = slicer_.SplitTo(v);
		if (equiv_vts.size() == 2) {
			sliced_mesh_.data(equiv_vts[0]).set_equivalent_vertex(equiv_vts[1]);
			sliced_mesh_.data(equiv_vts[1]).set_equivalent_vertex(equiv_vts[0]);
		}
		for (auto it = equiv_vts.begin(); it != equiv_vts.end(); ++it) {
			VertexHandle nv = *it;
			double u = mesh_.data(v).u();
			sliced_mesh_.data(nv).set_u(u);
			sliced_mesh_.data(nv).set_singularity(mesh_.data(v).is_singularity());
		}
	}

	sliced_mesh_.RequestBoundary();
	auto boundary = sliced_mesh_.GetBoundaries().front();
	//std::list<HalfedgeHandle> boundary_list(boundary.begin(), boundary.end());
	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		VertexHandle v = sliced_mesh_.to_vertex_handle(*it);
		if (sliced_mesh_.data(v).is_singularity())
			cone_vts_.push_back(v);
	}

}

double BFFSolver::CosineLaw(double a, double b, double c)
{
	double cs = (a * a + b * b - c * c) / (2 * a * b);
	if (-1 > cs)
		return PI;
	else if (cs > 1)
		return 0;
	else
		return acos(cs);
}

void BFFSolver::ComputeCornerAngles(SurfaceMesh &mesh, Eigen::VectorXd &L)
{
	using namespace OpenMesh;
	bool with_length = false;
	if (L.size() == mesh.n_edges())
		with_length = true;
	for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); ++fiter) {
		FaceHandle f = *fiter;
		std::vector<HalfedgeHandle> he;

		for (auto fhiter = mesh.fh_iter(f); fhiter.is_valid(); ++fhiter) {
			he.push_back(*fhiter);
		}
		double l[3];
		if (with_length) {
			l[0] = L(mesh.edge_handle(he[0]).idx());
			l[1] = L(mesh.edge_handle(he[1]).idx());
			l[2] = L(mesh.edge_handle(he[2]).idx());
		}
		else {
			l[0] = mesh.calc_edge_length(mesh.edge_handle(he[0]));
			l[1] = mesh.calc_edge_length(mesh.edge_handle(he[1]));
			l[2] = mesh.calc_edge_length(mesh.edge_handle(he[2]));
		}
		for (int i = 0; i < 3; ++i) {
			double cs = CosineLaw(l[i], l[(i + 1) % 3], l[(i + 2) % 3]);
			mesh.data(he[i]).set_angle(cs);
		}
	}
}

void BFFSolver::ComputeVertexCurvatures(SurfaceMesh &mesh, Eigen::VectorXd &L)
{
	using namespace OpenMesh;
	ComputeCornerAngles(mesh, L);
	double sum = 0.0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		double angle_sum = 0.0;
		for (auto vihiter = mesh.vih_iter(v); vihiter.is_valid(); ++vihiter) {
			HalfedgeHandle h = *vihiter;
			if (mesh.is_boundary(h)) continue;
			angle_sum += mesh.data(h).angle();
		}
		if(!mesh.is_boundary(v))
			mesh.data(v).set_curvature(2 * PI - angle_sum);
		else
			mesh.data(v).set_curvature(PI - angle_sum);
		sum += mesh.data(v).curvature();
	}
	std::cout << "Total Curvature: " << sum/ PI << " pi."<<  std::endl;
}

void BFFSolver::ComputeLaplacian(SurfaceMesh & mesh, bool mode)
{
	using namespace OpenMesh;
	using namespace Eigen;
	
	ReindexVertices(mesh);
	ComputeHalfedgeWeights(mesh);

	Delta_.resize(mesh.n_vertices(), mesh.n_vertices());
	Delta_.setZero();
	std::vector<Eigen::Triplet<double>> A_coefficients;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mode && (viter == mesh.vertices_begin())){
			for (auto viter1 = mesh.vertices_begin(); viter1 != mesh.vertices_end(); ++viter1) {
				VertexHandle v1 = *viter1;
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v1).reindex(), 1.));
			}
			continue;
		}
		double s_w = 0;
		for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(v); vviter.is_valid(); ++vviter) {
			VertexHandle neighbor = *vviter;
			HalfedgeHandle h = mesh.find_halfedge(v, neighbor);
			double n_w = mesh.data(h).weight();
			
			s_w += n_w;
			A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(neighbor).reindex(), -n_w));
		}
		A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v).reindex(), s_w));
	}
	Delta_.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
}

void BFFSolver::ComputeConformalFactors()
{
	using namespace Eigen;
	SurfaceMesh &mesh = mesh_;
	using namespace OpenMesh;
	
	
	// Using Cherrier Formula
	
	// Construct Sparse system;
	ComputeLaplacian(mesh, true);
	ComputeVertexCurvatures(mesh_);
	Eigen::VectorXd b(mesh.n_vertices());
		
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		b(mesh.data(v).reindex()) = mesh.data(v).curvature() - mesh.data(v).target_curvature();
	}
	b(mesh.data(*mesh.vertices_begin()).reindex()) = 0;
	
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
	solver.compute(Delta_);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::VectorXd u = solver.solve(b);
	u.setZero();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.data(v).set_u(u(mesh.data(v).reindex()));
	}
	

}

void BFFSolver::ComputeHalfedgeWeights(SurfaceMesh &mesh)
{
	using namespace OpenMesh;
	ComputeCornerAngles(mesh);
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		EdgeHandle e = *eiter;
		HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
		HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);
		HalfedgeHandle h0_next;
		if (!mesh.is_boundary(h0))
			h0_next = mesh.next_halfedge_handle(h0);
		HalfedgeHandle h1_next;
		if (!mesh.is_boundary(h1))
			h1_next = mesh.next_halfedge_handle(h1);

		double weight = 0.;

		if (h0_next.is_valid())
			weight += 1. / tan(mesh.data(h0_next).angle());
		if (h1_next.is_valid())
			weight += 1. / tan(mesh.data(h1_next).angle());
		weight *= 0.5;
		//if (weight < 0)
			//weight = 0.01;
		mesh.data(h0).set_weight(weight);
		mesh.data(h1).set_weight(weight);
	}
}

void BFFSolver::ReindexVertices(SurfaceMesh & mesh)
{
	using namespace OpenMesh;
	int n_interior = 0;
	int n_boundary = 0;

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (!mesh.is_boundary(v)) {
			mesh.data(v).set_reindex(n_interior);
			++n_interior;
		}
	}

	
	if (mesh.GetBoundaries().size() == 0) return;
	auto boundary = mesh.GetBoundaries().front();
	std::list<HalfedgeHandle> boundary_list(boundary.begin(), boundary.end());
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		if (!mesh.data(mesh.from_vertex_handle(*it)).equivalent_vertex().is_valid()) {
			boundary_list.insert(boundary_list.end(), boundary_list.begin(), it);
			boundary_list.erase(boundary_list.begin(), it);
			break;
		}
	}
	
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		mesh.data(mesh.from_vertex_handle(*it)).set_reindex(n_interior+n_boundary);
		++n_boundary;
	}



	n_boundary_ = n_boundary;
	n_interior_ = n_interior;
}



void BFFSolver::BoundaryUToTargetK()
{
	using namespace Eigen;
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	
	VectorXd u(mesh.n_vertices());
	ReindexVertices(mesh);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		u(mesh.data(v).reindex()) = mesh.data(v).u();
	}
	
	Eigen::VectorXd u_B = u.segment(n_interior_, n_boundary_);

	Eigen::VectorXd target_k = BoundaryUToTargetK(u_B);

	auto boundary = sliced_mesh_.GetBoundaries().front();
	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		VertexHandle v = sliced_mesh_.to_vertex_handle(*it);
		//if (!mesh.data(v).is_singularity())
			mesh.data(v).set_target_curvature(target_k(mesh.data(v).reindex() - n_interior_));
		//else {
		//	mesh.data(v).set_target_curvature(PI / 2);
		//	target_k(mesh.data(v).reindex() - n_interior_) = PI / 2;
		//}
	}

	std::cout << "Singularities' curvature:" << std::endl;
	for (auto it = cone_vts_.begin(); it != cone_vts_.end(); ++it) {
		VertexHandle v = *it;
		std::cout << mesh.data(v).target_curvature() / PI << "pi" << std::endl;
	}

	Eigen::VectorXd b(mesh.n_vertices());
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		b(mesh.data(v).reindex()) = mesh.is_boundary(v) ? mesh.data(v).curvature() - mesh.data(v).target_curvature() : mesh.data(v).curvature();
	}


	/*Eigen::VectorXd new_u = BoundaryTargetKToU(target_k);
	std::cout << new_u << std::endl << std::endl << u_B << std::endl;

	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		VertexHandle v = sliced_mesh_.to_vertex_handle(*it);
			mesh.data(v).set_u(new_u(mesh.data(v).reindex() - n_interior_));
	}*/


}


Eigen::VectorXd BFFSolver::BoundaryUToTargetK(Eigen::VectorXd & u)
{
	using namespace Eigen;
	using namespace OpenMesh;
	
	SurfaceMesh &mesh = sliced_mesh_;

	ComputeLaplacian(mesh);
	VectorXd omega(mesh.n_vertices());
	VectorXd k(n_boundary_);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		omega(mesh.data(v).reindex()) = mesh.data(v).curvature();
	}

	k = omega.segment(n_interior_, n_boundary_);
	SparseLU<SparseMatrix<double>> solver;
	solver.compute(Delta_.block(0, 0, n_interior_, n_interior_));
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	Eigen::SparseMatrix<double> A_IB = Delta_.block(0, n_interior_, n_interior_, n_boundary_);
	Eigen::SparseMatrix<double> A_BB = Delta_.block(n_interior_, n_interior_, n_boundary_, n_boundary_);
	Eigen::SparseMatrix<double> A_BI = Delta_.block(n_interior_, 0, n_boundary_, n_interior_);
	VectorXd omega_I = omega.segment(0, n_interior_);
	Eigen::VectorXd to_inverse = omega_I - A_IB * u;
	Eigen::VectorXd inverse = solver.solve(to_inverse);
	assert((Delta_.block(0, 0, n_interior_, n_interior_) * inverse - to_inverse).norm() < 1e-7);
	Eigen::VectorXd h = A_BI * inverse + A_BB * u;

	Eigen::VectorXd target_k = k - h;

	std::cout << "Sum of target total curvature:" << target_k.sum() / PI << " pi" << std::endl;

	return target_k;

}

Eigen::VectorXd BFFSolver::BoundaryTargetKToU(Eigen::VectorXd & target_k)
{
	using namespace Eigen;
	using namespace OpenMesh;

	SurfaceMesh &mesh = sliced_mesh_;
	ComputeLaplacian(mesh, true);

	VectorXd omega(mesh.n_vertices());
	VectorXd h(mesh.n_vertices());

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (!mesh.is_boundary(v)) {
			omega(mesh.data(v).reindex()) = mesh.data(v).curvature();
			h(mesh.data(v).reindex()) = 0.;
		}
		else {
			omega(mesh.data(v).reindex()) = 0;
			h(mesh.data(v).reindex()) = mesh.data(v).curvature() - target_k(mesh.data(v).reindex() - n_interior_);
		}

	}

	omega(mesh.data(*mesh.vertices_begin()).reindex()) = 0;


	SparseLU<SparseMatrix<double>> solver;
	solver.compute(Delta_);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	} 

	Eigen::VectorXd u = solver.solve(omega  + h);

	return u.segment(n_interior_, n_boundary_);


}

void BFFSolver::IntegrateBoundaryCurve()
{
	using namespace OpenMesh;
	using namespace Eigen;
	SurfaceMesh &mesh = sliced_mesh_;

	HPropHandleT<int> reindex;
	mesh.add_property(reindex);

	VPropHandleT<double> cumulative_angle;
	VPropHandleT<Vec2d> tangent;
	mesh.add_property(cumulative_angle);
	mesh.add_property(tangent);

	auto boundary = mesh.GetBoundaries().front();
	std::list<HalfedgeHandle> boundary_list(boundary.begin(), boundary.end());
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		if (mesh.data(mesh.from_vertex_handle(*it)).is_singularity()) {
			boundary_list.insert(boundary_list.end(), boundary_list.begin(), it);
			boundary_list.erase(boundary_list.begin(), it);
			break;
		}
	}
	boundary = std::vector<HalfedgeHandle>(boundary_list.begin(), boundary_list.end());
	// Set edge length and Construct tanget matrix and normalization matrix;
	int index = 0;
	for (auto it = boundary.begin(); it != boundary.end(); ++it, ++index) {
		HalfedgeHandle h = *it;
		mesh.property(reindex, h) = index;
	}

	std::vector<int> oppo_relation(boundary.size());

	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		HalfedgeHandle h = *it;
		HalfedgeHandle oppo = mesh.opposite_halfedge_handle(mesh.data(mesh.opposite_halfedge_handle(h)).original_opposition());
		oppo_relation[mesh.property(reindex, h)] = mesh.property(reindex, oppo);
		oppo_relation[mesh.property(reindex, oppo)] = mesh.property(reindex, h);
	}

	Eigen::MatrixXd T(2, n_boundary_);
	Eigen::MatrixXd N_inverse(n_boundary_, n_boundary_);
	T.setZero();
	N_inverse.setZero();

	Eigen::VectorXd L_star(n_boundary_);


	mesh.property(cumulative_angle, mesh.from_vertex_handle(boundary.front())) = 0;
	int i = 0;
	for (auto it = boundary.begin(); it != boundary.end(); ++it, ++i) {
		VertexHandle v0 = mesh.from_vertex_handle(*it);
		VertexHandle v1 = mesh.to_vertex_handle(*it);
		EdgeHandle e = mesh.edge_handle(*it);
		double u0 = mesh.data(v0).u();
		double u1 = mesh.data(v1).u();
		double l = mesh.calc_edge_length(e);
		double l_star = exp(0.5 * (u0 + u1)) * l;
		double angle = mesh.property(cumulative_angle, v0);
		mesh.property(cumulative_angle, v1) = angle + mesh.data(v0).target_curvature();
		mesh.property(tangent, v0) = Vec2d(cos(angle), sin(angle));

		T(0,i) = cos(angle);
		T(1, i) = sin(angle);
		L_star(i) = l_star;

	}
	
	N_inverse(0, 0) = 0.5 * (L_star(L_star.size() - 1) + L_star[0]);
	for (int i = 1; i < boundary.size(); ++i) {
		N_inverse(i, i) = 0.5 *(L_star(i - 1) + L_star(i));
	}

	Eigen::MatrixXd newT(T.rows(), T.cols()/2);
	Eigen::VectorXd L_star_half(L_star.size() / 2);
	Eigen::MatrixXd N_inverse_half(N_inverse.rows() / 2, N_inverse.cols() / 2);
	
	index = 0;
	for (int i = 0; i < boundary.size(); ++i) {
		if (i > oppo_relation[i]) continue;
		newT.col(index) = T.col(i) + T.col(oppo_relation[i]);
		L_star_half(index) = L_star(i);
		N_inverse_half(index, index) = N_inverse(i, i);
		++index;
	}

	Eigen::VectorXd L_normalized_half = L_star_half - N_inverse_half * newT.transpose() * (newT * N_inverse_half* newT.transpose()).inverse() * newT * L_star_half;
	Eigen::VectorXd L_normalized(L_star.size());
	
	index = 0;
	for (int i = 0; i < boundary.size(); ++i) {
		if (i > oppo_relation[i]) continue;
		L_normalized(i) = L_normalized_half(index);
		L_normalized(oppo_relation[i]) = L_normalized_half(index);
		++index;
	}

	L_normalized = L_star - N_inverse * T.transpose() * (T * N_inverse* T.transpose()).inverse() * T * L_star;

	i = 0;
	mesh.set_texcoord2D(mesh.from_vertex_handle(boundary.front()), Vec2d(0, 0));
	for (auto it = boundary.begin(); it != boundary.end(); ++it, ++i) {
		VertexHandle v0 = mesh.from_vertex_handle(*it);
		VertexHandle v1 = mesh.to_vertex_handle(*it);
		Vec2d coord = mesh.texcoord2D(v0) + mesh.property(tangent, v0) * L_normalized(i);
		mesh.set_texcoord2D(v1, coord);
	}
	
}

void BFFSolver::ExtendToInteriorHilbert()
{
	using namespace Eigen;
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	ComputeHarmonicMatrix();
	
	Eigen::VectorXd a_boundary(mesh.n_vertices());
	a_boundary.setZero();

	auto boundary = mesh.GetBoundaries().front();
	for (int i = 0; i < boundary.size(); ++i) {
		VertexHandle v = mesh.to_vertex_handle(boundary[i]);
		Vec2d coord = mesh.texcoord2D(v);
		a_boundary(v.idx()) = coord[0];
	}

	SparseLU<SparseMatrix<double>> solver;
	solver.compute(Delta_);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::VectorXd a = solver.solve(a_boundary);

	ComputeLaplacian(mesh);
	Eigen::VectorXd h(mesh.n_vertices());
	h.setZero();
	for (int i = 0; i < boundary.size(); ++i) {
		VertexHandle v = mesh.to_vertex_handle(boundary[i]);
		VertexHandle v_prev = mesh.from_vertex_handle(boundary[i]);
		VertexHandle v_next = mesh.to_vertex_handle(boundary[(i + 1) % boundary.size()]);
		h(mesh.data(v).reindex()) = 0.5 * (a(v_prev.idx()) - a(v_next.idx()));
	}

	solver.compute(Delta_);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::VectorXd b = solver.solve(h);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle  v = *viter;
		mesh.set_texcoord2D(v, Vec2d(a(v.idx()), b(mesh.data(v).reindex())));
	}
}

void BFFSolver::ExtendToInteriorHarmonic()
{
	using namespace Eigen;
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	ComputeHarmonicMatrix();

	Eigen::MatrixXd uv_boundary(mesh.n_vertices(), 2);
	uv_boundary.setZero();

	auto boundary = mesh.GetBoundaries().front();
	for (int i = 0; i < boundary.size(); ++i) {
		VertexHandle v = mesh.to_vertex_handle(boundary[i]);
		Vec2d coord = mesh.texcoord2D(v);
		uv_boundary(v.idx(), 0) = coord[0];
		uv_boundary(v.idx(), 1) = coord[1];
	}

	SparseLU<SparseMatrix<double>> solver;
	solver.compute(Delta_);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::MatrixXd uv = solver.solve(uv_boundary);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle  v = *viter;
		mesh.set_texcoord2D(v, Vec2d(uv(v.idx(), 0), uv(v.idx(), 1)));
	}
}

void BFFSolver::ComputeHarmonicMatrix()
{
	using namespace OpenMesh;
	using namespace Eigen;
	SurfaceMesh &mesh = sliced_mesh_;

	Delta_.resize(mesh.n_vertices(), mesh.n_vertices());
	Delta_.setZero();
	std::vector<Eigen::Triplet<double>> A_coefficients;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh.is_boundary(v)) {
			A_coefficients.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.));
			continue;
		}
		double s_w = 0;
		for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(v); vviter.is_valid(); ++vviter) {
			VertexHandle neighbor = *vviter;
			HalfedgeHandle h = mesh.find_halfedge(v, neighbor);
			double n_w = mesh.data(h).weight();
			s_w += n_w;
			A_coefficients.push_back(Eigen::Triplet<double>(v.idx(), neighbor.idx(), -n_w));
		}
		A_coefficients.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), s_w));
	}
	Delta_.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
}

void BFFSolver::NormalizeUV()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	Vec2d s(0, 0);
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		s += mesh.texcoord2D(v);
	}
	s /= mesh.n_vertices();
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.set_texcoord2D(v, mesh.texcoord2D(v) -s);
	}

	double scale = 0;
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		scale = mesh.texcoord2D(v).norm() > scale ? mesh.texcoord2D(v).norm() : scale;
	}
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.set_texcoord2D(v, mesh.texcoord2D(v) / scale);
	}
}

void BFFSolver::ComputeOrbifoldBoundaryData()
{
	using namespace OpenMesh;
	using namespace Eigen;
	SurfaceMesh &mesh = sliced_mesh_;
	ReindexVertices(mesh);
	
	Eigen::SparseMatrix<double> A(mesh.n_vertices() + 2 + (n_boundary_-2)/2 + 2, mesh.n_vertices() + n_boundary_);
	A.setZero();
	std::vector<Eigen::Triplet<double>> A_coefficients;
	Eigen::VectorXd b(mesh.n_vertices() + 2 + (n_boundary_ - 2) / 2 + 2);
	b.setZero();

	ComputeHalfedgeWeights(mesh);
		

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;

		if (mesh.is_boundary(v)) {
			VertexHandle equiv = mesh.data(v).equivalent_vertex();
			if (equiv.is_valid() && mesh.data(v).reindex() > mesh.data(equiv).reindex()) continue;
			
			if (equiv.is_valid()) {
				double s_w = 0;
				for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(v); vviter.is_valid(); ++vviter) {
					VertexHandle neighbor = *vviter;
					HalfedgeHandle h = mesh.find_halfedge(v, neighbor);
					double n_w = mesh.data(h).weight();
					s_w += n_w;
					A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(neighbor).reindex(), -n_w));
				}
				for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(equiv); vviter.is_valid(); ++vviter) {
					VertexHandle neighbor = *vviter;
					HalfedgeHandle h = mesh.find_halfedge(equiv, neighbor);
					double n_w = mesh.data(h).weight();
					s_w += n_w;
					A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(neighbor).reindex(), -n_w));
				}

				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v).reindex(), s_w));
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v).reindex() + n_boundary_, 1));
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(equiv).reindex() + n_boundary_, 1));
				b(mesh.data(v).reindex()) = mesh.data(v).curvature() + mesh.data(equiv).curvature();

				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(equiv).reindex(), mesh.data(v).reindex(), 1));
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(equiv).reindex(), mesh.data(equiv).reindex(), -1));

				if (mesh.data(v).is_singularity()) {
					assert(mesh.data(v).reindex() + n_boundary_ < A.rows());
					//assert(mesh.data(equiv).reindex() + n_boundary_ < A.rows());
					A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex() + n_boundary_, mesh.data(v).reindex() + n_boundary_, 1));
					b(mesh.data(v).reindex() + n_boundary_) = PI / 2;
					A_coefficients.push_back(Eigen::Triplet<double>(A.rows() - 2, mesh.data(equiv).reindex() + n_boundary_, 1));
					b(A.rows() - 2) = PI / 2;
					
				}
				else {
					int index = mesh.data(v).reindex() + n_boundary_;
					int equiv_index = mesh.data(equiv).reindex();
					assert(mesh.data(v).reindex() + n_boundary_ < A.rows());
					A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex() + n_boundary_, mesh.data(v).reindex() + n_boundary_, 1));
					A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex() + n_boundary_, mesh.data(equiv).reindex() + n_boundary_, 1));
				}
			}
			else{
				double s_w = 0;
				for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(v); vviter.is_valid(); ++vviter) {
					VertexHandle neighbor = *vviter;
					HalfedgeHandle h = mesh.find_halfedge(v, neighbor);
					double n_w = mesh.data(h).weight();
					s_w += n_w;
					A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(neighbor).reindex(), -n_w));
				}
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v).reindex(), s_w));
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v).reindex() + n_boundary_, 1));
				b(mesh.data(v).reindex()) = mesh.data(v).curvature();
				assert(mesh.data(v).reindex() + n_boundary_ < A.rows());
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex() + n_boundary_, mesh.data(v).reindex() + n_boundary_, 1));
				b(mesh.data(v).reindex() + n_boundary_) = PI / 2;
			}
		}
		else // inner vertices
		{
			double s_w = 0;
			for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(v); vviter.is_valid(); ++vviter) {
				VertexHandle neighbor = *vviter;
				HalfedgeHandle h = mesh.find_halfedge(v, neighbor);
				double n_w = mesh.data(h).weight();
				s_w += n_w;
				A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(neighbor).reindex(), -n_w));
			}
			A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v).reindex(), s_w));
			A_coefficients.push_back(Eigen::Triplet<double>(mesh.data(v).reindex(), mesh.data(v).reindex() + n_boundary_, 1));
			b(mesh.data(v).reindex()) = mesh.data(v).curvature();
		}
	}

	for (int i = 0; i < mesh.n_vertices(); ++i) {
		A_coefficients.push_back(Triplet<double>(A.rows() - 1, i, 1));
	}

	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());

	SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
	solver.compute(A);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::VectorXd uk = solver.solve(b);
	VectorXd u_B = uk.segment(n_interior_, n_boundary_);
	//VectorXd target_k = BoundaryUToTargetK(u_B);
	VectorXd target_k = uk.segment(mesh.n_vertices(), n_boundary_);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh.is_boundary(v)) {
			mesh.data(v).set_u(u_B(mesh.data(v).reindex() - n_interior_));
			mesh.data(v).set_target_curvature(target_k(mesh.data(v).reindex() - n_interior_));
		}
	}

	

	std::cout << "Singularities' curvature:" << std::endl;
	for (auto it = cone_vts_.begin(); it != cone_vts_.end(); ++it) {
		VertexHandle v = *it;
		std::cout << mesh.data(v).target_curvature() / PI << "pi" << std::endl;
	}

	
	
}

