#include "BFF.h"

BFFSolver::BFFSolver(SurfaceMesh & mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag)
	:mesh_(mesh), cone_flag_(cone_flag), cone_angle_(cone_angle), slice_flag_(slice_flag)
{

}


SurfaceMesh BFFSolver::Compute(int mode)
{
	/*
	mode 0: BFF known k with hilbert extension
	mode 1: BFF known k with harmonic extension
	mode 2: BFF free boundary with hilbert extension
	mode 3: BFF free boundary with harmonic extension
	mode 4: BFF cone parameterization with hilbert extension
	mode 5: BFF cone parameterization with harmonic extension
	*/

	Init();

	ComputeVertexCurvatures(sliced_mesh_);
	
	if (mode == 2 || mode == 3)
		FreeBoundary();
	else if (mode == 0 || mode == 1)
		BoundaryTargetKKnown();
	else if (mode == 4 || mode == 5)
		GlobalParameterization();

	IntegrateBoundaryCurve();

	if(mode == 0 || mode == 2 || mode == 4)
		ExtendToInteriorHilbert();
	else
		ExtendToInteriorHarmonic();

	NormalizeUV();

	return sliced_mesh_;
}

void BFFSolver::Init()
{
	BFFInitializer initializer(mesh_);
	initializer.Initiate(sliced_mesh_, cone_flag_, cone_angle_, slice_flag_);
	cone_vts_ = initializer.GetConeVertices();
	split_to_ = initializer.split_to();
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

	mesh.RequestBoundary();
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

void BFFSolver::BoundaryTargetKKnown()
{
	using namespace Eigen;
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	VectorXd target_k(mesh.n_vertices());
	target_k.setZero();

	ReindexVertices(mesh);

	std::cout << "Singularities' curvature:" << std::endl;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh.data(v).is_singularity()) {
			target_k(mesh.data(v).reindex()) = mesh.data(v).target_curvature();
			std::cout << mesh.data(v).target_curvature() / PI << "pi" << std::endl;
		}
	}
	VectorXd k_B = target_k.segment(n_interior_, n_boundary_);
	VectorXd u = BoundaryTargetKToU(k_B);

	mesh.RequestBoundary();
	auto boundary = mesh.GetBoundaries().front();
	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		VertexHandle v = mesh.to_vertex_handle(*it);
		mesh.data(v).set_u(u(mesh.data(v).reindex() - n_interior_));
	}

}



void BFFSolver::FreeBoundary()
{
	using namespace Eigen;
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	
	VectorXd u(mesh.n_vertices());
	ReindexVertices(mesh);

	u.setZero();

	Eigen::VectorXd u_B = u.segment(n_interior_, n_boundary_);

	Eigen::VectorXd target_k = BoundaryUToTargetK(u_B);

	auto boundary = sliced_mesh_.GetBoundaries().front();
	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		VertexHandle v = sliced_mesh_.to_vertex_handle(*it);
		mesh.data(v).set_target_curvature(target_k(mesh.data(v).reindex() - n_interior_));
		mesh.data(v).set_u(0);
	}

	std::cout << "Singularities' curvature:" << std::endl;
	for (auto it = cone_vts_.begin(); it != cone_vts_.end(); ++it) {
		VertexHandle v = *it;
		std::cout << mesh.data(v).target_curvature() / PI << "pi" << std::endl;
	}

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.data(v).reindex();
	}


}

void BFFSolver::GlobalParameterization()
{
	using namespace Eigen;
	SurfaceMesh &mesh = mesh_;
	using namespace OpenMesh;
	// Using Cherrier Formula
	// Construct Sparse system;
	ComputeLaplacian(mesh, true);
	ComputeVertexCurvatures(mesh);
	Eigen::VectorXd b(mesh.n_vertices());
	
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (!mesh.property(cone_flag_, v)) {
			b(mesh.data(v).reindex()) = mesh.data(v).curvature();
		}
		else {
			if (mesh.is_boundary(v))
				b(mesh.data(v).reindex()) = mesh.data(v).curvature() - (PI - mesh.property(cone_angle_, v));
			else
				b(mesh.data(v).reindex()) = mesh.data(v).curvature() - (2 * PI - mesh.property(cone_angle_, v));
		}
	}
	
	b(mesh.data(*mesh.vertices_begin()).reindex()) = 0;
	
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
	solver.compute(Delta_);
	if (solver.info() != Eigen::Success){
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::VectorXd u = solver.solve(b);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.data(v).set_u(u(mesh.data(v).reindex()));
	}

	ReindexVertices(sliced_mesh_);
	u.resize(sliced_mesh_.n_vertices());
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		auto verts = mesh.property(split_to_, v);
		for (auto it = verts.begin(); it != verts.end(); ++it) {
			sliced_mesh_.data(*it).set_u(mesh.data(v).u());
			u(sliced_mesh_.data(*it).reindex()) = mesh.data(v).u();
		}
	}

	
	// Convert u to k
	VectorXd uB = u.segment(n_interior_, n_boundary_);
	Eigen::VectorXd target_k = BoundaryUToTargetK(uB);
	auto boundary = sliced_mesh_.GetBoundaries().front();
	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		VertexHandle v = sliced_mesh_.to_vertex_handle(*it);
		sliced_mesh_.data(v).set_target_curvature(target_k(sliced_mesh_.data(v).reindex() - n_interior_));
	}

	std::cout << "Singularities' curvature:" << std::endl;
	for (auto it = cone_vts_.begin(); it != cone_vts_.end(); ++it) {
		VertexHandle v = *it;
		std::cout << sliced_mesh_.data(v).target_curvature() / PI << "pi" << std::endl;
	}

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
			h(mesh.data(v).reindex()) = (mesh.data(v).curvature() - target_k(mesh.data(v).reindex() - n_interior_));
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

	VPropHandleT<double> cumulative_angle;
	VPropHandleT<Vec2d> tangent;
	mesh.add_property(cumulative_angle);
	mesh.add_property(tangent);

	Eigen::MatrixXd T(2, n_boundary_); // Tangent vector

	T.setZero();

	Eigen::VectorXd L_star(n_boundary_); // edge length after conformal map
	Eigen::VectorXd L(n_boundary_); // mesh euclidean length.
	Eigen::VectorXd L_dual(n_boundary_);
	auto boundary = mesh.GetBoundaries().front();
	mesh.property(cumulative_angle, mesh.from_vertex_handle(boundary.front())) = mesh.data(mesh.from_vertex_handle(boundary.front())).target_curvature();
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
		mesh.property(cumulative_angle, v1) = angle + mesh.data(v1).target_curvature();
		mesh.property(tangent, v0) = Vec2d(cos(angle), sin(angle));

		T(0, i) = cos(angle);
		T(1, i) = sin(angle);
		L_star(i) = l_star;
		L(i) = l;
	}

	for (int i = 0; i < n_boundary_; ++i) {
		L_dual((i + 1) % n_boundary_) = 1 / L_star(i);
	}


	// Reindex boundary halfedges
	HPropHandleT<int> reindex;
	mesh.add_property(reindex);

	int index = 0;
	for (auto it = boundary.begin(); it != boundary.end(); ++it, ++index) {
		HalfedgeHandle h = *it;
		mesh.property(reindex, h) = index;
	}

	// Equivalent edges produced by cut should be of the same length.
	std::vector<int> oppo_relation(boundary.size());
	
	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		HalfedgeHandle h = *it;
		HalfedgeHandle oppo_inner = mesh.data(mesh.opposite_halfedge_handle(h)).original_opposition();
		if (oppo_inner.is_valid()) {
			HalfedgeHandle oppo = mesh.opposite_halfedge_handle(oppo_inner);
			oppo_relation[mesh.property(reindex, h)] = mesh.property(reindex, oppo);
			oppo_relation[mesh.property(reindex, oppo)] = mesh.property(reindex, h);
		}
		else {
			oppo_relation[mesh.property(reindex, h)] = -1;
		}
	}

	std::vector<HalfedgeHandle> valid_halfedges;
	for (int i = 0; i < oppo_relation.size(); ++i) {
		if (i > oppo_relation[i] && oppo_relation[i] != -1) continue;
		valid_halfedges.push_back(boundary[i]);
	}


	int n_valid_h = valid_halfedges.size();
	
	Eigen::MatrixXd N_inverse(n_valid_h, n_valid_h);
	N_inverse.setZero();
	for (int i = 1; i < n_valid_h; ++i) {
		HalfedgeHandle h = valid_halfedges[i];
		int index = mesh.property(reindex, h);
		int oppo_index = oppo_relation[index];
		N_inverse(i, i) = 1. / L_dual(index);
		if(oppo_index >= 0)
			N_inverse(i, i) += 1. / L_dual(oppo_index);
	}

	Eigen::MatrixXd newT(T.rows(), n_valid_h);
	Eigen::VectorXd L_star_valid(n_valid_h);
	
	for (int i = 0; i < n_valid_h; ++i) {
		HalfedgeHandle h = valid_halfedges[i];
		int index = mesh.property(reindex, h);
		int oppo_index = oppo_relation[index];
		newT.col(i) = T.col(index);
		if (oppo_index >= 0)
			newT.col(i) += T.col(oppo_index);
		L_star_valid(i) = L_star(index);
	}

	// Use quadratic programming to get optimal solutions.
	Eigen::SparseMatrix<double> Q = (N_inverse * N_inverse).sparseView();
	Eigen::VectorXd B = - 2 * N_inverse.diagonal();
	
	Eigen::SparseMatrix<double> A_eq(2, Q.rows());
	A_eq = newT.sparseView();
	Eigen::VectorXd B_eq(2); B_eq.setZero();

	Eigen::SparseMatrix<double> A_ieq = (-Eigen::MatrixXd::Identity(Q.rows(), Q.cols())).sparseView();
	Eigen::VectorXd B_ieq = -Eigen::VectorXd::Constant(A_ieq.rows(), -1e-7);

	Eigen::VectorXd lx = Eigen::VectorXd::Constant(Q.cols(), -1000);
	Eigen::VectorXd ux = Eigen::VectorXd::Constant(Q.cols(), 1000);

	Eigen::VectorXd L_normalized_valid(L_star_valid.size());
	L_normalized_valid = L_star_valid;
	igl::active_set_params as;
	igl::active_set(Q, B, Eigen::VectorXi(), Eigen::VectorXd(), A_eq, B_eq, A_ieq, B_ieq, lx, ux, as, L_normalized_valid);
	Eigen::VectorXd L_normalized(L_star.size());
	

	for (int i = 0; i < n_valid_h; ++i) {
		HalfedgeHandle h = valid_halfedges[i];
		int index = mesh.property(reindex, h);
		int oppo_index = oppo_relation[index];
		L_normalized(index) = L_normalized_valid[i];
		if (oppo_index >= 0)
			L_normalized(oppo_index) = L_normalized_valid(i);
	}
	
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
		h(mesh.data(v).reindex()) = -0.5 * (a(v_next.idx()) - a(v_prev.idx()));
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


