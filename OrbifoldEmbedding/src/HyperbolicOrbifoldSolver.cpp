#include "HyperbolicOrbifoldSolver.h"



HyperbolicOrbifoldSolver::HyperbolicOrbifoldSolver(SurfaceMesh & mesh): mesh_(mesh)
{

}

SurfaceMesh HyperbolicOrbifoldSolver::Compute()
{
	using namespace Eigen;
	using namespace LBFGSpp;

	if (mesh_.n_vertices() > 10) {
		InitOrbifold();
		ComputeCornerAngles();
		ComputeHalfedgeWeights();
		InitMap();
		Normalize();
		this->ComputeGradient();
		std::function<double(const VectorXd&, VectorXd&)> fun=
			[this](const VectorXd& x, VectorXd& grad) ->double 
		{
			this->SetCoords(x);
			this->Normalize();
			this->ComputeGradient();
			grad = this->GetGradientVector();
			return this->ComputeEnergy();
		};

		LBFGSParam<double> param;
		param.epsilon = 1e-7;
		param.max_iterations = 2000;
		//param.min_step = 1e-5;
		LBFGSSolver<double> solver(param);
		VectorXd x = GetCoordsVector();
		
		double fx;
		int niter = solver.minimize(fun, x, fx);

		std::cout << niter << " iterations" << std::endl;
		std::cout << "f(x) = " << fx << std::endl;

		//double energy = OptimizationLoop(0.01, 1e-4);
		//std::cout << "Final Energy:" << energy << std::endl;
	}
	return sliced_mesh_;
}

void HyperbolicOrbifoldSolver::CutToDist(int n_cones)
{
	n_cones_ = n_cones;
	OrbifoldMeshSlicer slicer(mesh_);
	sliced_mesh_ = slicer.CutAndSelectSingularities(n_cones);
	cone_vts_ = slicer.GetConeVertices();
	segments_vts_ = slicer.GetSegments();
}

void HyperbolicOrbifoldSolver::InitOrbifold()
{
	InitType1();
}

void HyperbolicOrbifoldSolver::InitType1()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	CutToDist(7);

	double radius = AngleCosineLaw(2 * PI / n_cones_, PI / 4., PI / 4.);
	double ratio = (exp(radius) - 1) / (exp(radius) + 1);

	Vec2d p1(cos(PI / 2 + (1 + n_cones_ / 2) * 2 * PI / n_cones_)*ratio, sin(PI / 2 + (1 + n_cones_ / 2) * 2 * PI / n_cones_)*ratio);
	Vec2d pk(cos(PI / 2 + (n_cones_ / 2) * 2 * PI / n_cones_)*ratio, sin(PI / 2 + (n_cones_ / 2) * 2 * PI / n_cones_)*ratio);
	double edge_length = AngleCosineLaw(PI / 4, PI / 4, 2 * PI / n_cones_);
	double target_radis = (exp(edge_length/2) - 1) / (exp(edge_length/2) + 1);
	Complex p1_source(p1[0], p1[1]);
	Complex pk_source(pk[0], pk[1]);
	Complex p1_target(target_radis, 0);
	Complex pk_target(-target_radis, 0);

	if (!vtx_transit_.is_valid()) {
		mesh.add_property(vtx_transit_);
	}
	auto transformation = ComputeMobiusTransformation(p1_source, pk_source, p1_target, pk_target);

	for (int i = 0; i < n_cones_; ++i) {
		Vec2d uv = Vec2d(cos(PI / 2 + (i + 1 + n_cones_ / 2) * 2 * PI / n_cones_)*ratio, sin(PI / 2 + (i + 1 + n_cones_ / 2) * 2 * PI / n_cones_)*ratio);
		Complex uv_complex = transformation(Complex(uv[0], uv[1]));
		uv = Vec2d(uv_complex.real(), uv_complex.imag());
		mesh.set_texcoord2D(cone_vts_[i], uv);
		uv_complex = std::conj(uv_complex);
		uv = Vec2d(uv_complex.real(), uv_complex.imag());
		if(i != 0 && i != n_cones_ - 1)
			mesh.set_texcoord2D(cone_vts_[2 * n_cones_ - 2 - i], uv);
	}
	
	std::cout << "Cone coordinates:\n";
	for (int i = 0; i < cone_vts_.size(); ++i) {
		Vec2d uv = mesh.texcoord2D(cone_vts_[i]);
		std::cout << uv[0] << "\t" << uv[1] << std::endl;
	}

	for (int i = 0; i < segments_vts_.size(); ++i) {
		Vec2d s0_uv = mesh.texcoord2D(cone_vts_[i]);
		Vec2d s1_uv = mesh.texcoord2D(cone_vts_[(i + 1) % cone_vts_.size()]);
		Vec2d t0_uv = mesh.texcoord2D(cone_vts_[(2 * n_cones_ - 2 - i) % cone_vts_.size()]);
		Vec2d t1_uv = mesh.texcoord2D(cone_vts_[(2 * n_cones_ - 2 - i - 1) % cone_vts_.size()]);

		Complex s0(s0_uv[0], s0_uv[1]);
		Complex s1(s1_uv[0], s1_uv[1]);
		Complex t0(t0_uv[0], t0_uv[1]);
		Complex t1(t1_uv[0], t1_uv[1]);
		/*std::cout << s0 << "\t" << s1 << std::endl;
		std::cout << t0 << "\t" << t1 << std::endl << std::endl;*/

		auto transformation = ComputeMobiusTransformation(s0, s1, t0, t1);
		assert(abs(transformation(s1) - t1) < 1e-6);
		// initiate the boundary's coordinates along the way.
		Vec2d segment_start = mesh.texcoord2D(cone_vts_[i]);
		Vec2d segment_end = mesh.texcoord2D(cone_vts_[(i + 1) % cone_vts_.size()]);
		Vec2d segment_vector = segment_end - segment_start;
		double segment_length = (segment_vector).norm();
		double interval_length = segment_length / (segments_vts_[i].size() + 1);
		int j = 1;
		for (auto it = segments_vts_[i].begin(); it != segments_vts_[i].end(); ++it, ++j) {
			mesh.property(vtx_transit_, *it) = transformation;
			mesh.set_texcoord2D(*it, segment_start + segment_vector * j * interval_length / segment_length);
		}
	}

}

double HyperbolicOrbifoldSolver::CosineLaw(double a, double b, double c)
{
	//double cs = (cosh(a) * cosh(b) - cosh(c)) / (sinh(a) * sinh(b));
	double cs = (a * a + b * b - c * c) / (2 * a * b);
	assert(-1 <= cs && cs <= 1);
	return acos(cs);

}

double HyperbolicOrbifoldSolver::AngleCosineLaw(double a, double b, double c)
{
	double l = (cos(c) + cos(a)*cos(b)) / (sin(a) * sin(b));
	return acosh(l);
}

void HyperbolicOrbifoldSolver::ComputeCornerAngles()
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

void HyperbolicOrbifoldSolver::ComputeHalfedgeWeights()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	ComputeCornerAngles();
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
		if (weight < 0)
			weight = 0.01;
		mesh.data(h0).set_weight(weight);
		mesh.data(h1).set_weight(weight);
	}

}

void HyperbolicOrbifoldSolver::ComputeEdgeLength()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		EdgeHandle e = *eiter;
		HalfedgeHandle h = mesh.halfedge_handle(e, 0);
		VertexHandle v0 = mesh.from_vertex_handle(h);
		VertexHandle v1 = mesh.to_vertex_handle(h);
		Vec2d v0_uv = mesh.texcoord2D(v0);
		Vec2d v1_uv = mesh.texcoord2D(v1);
		Complex v0_complex(v0_uv[0], v0_uv[1]);
		Complex v1_complex(v1_uv[0], v1_uv[1]);
		mesh.data(e).set_length(HyperbolicDistance(v0_complex, v1_complex));
	}
}

double HyperbolicOrbifoldSolver::ComputeGradient()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	double max_gradient_norm = 0.0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec2d gradient = ComputeGradient(v);
		mesh.data(v).set_gradient(gradient);
		if (gradient.norm() > max_gradient_norm) {
			max_gradient_norm = gradient.norm();
		}
	}
	return max_gradient_norm;
}

OpenMesh::Vec2d HyperbolicOrbifoldSolver::ComputeGradientOfDistance2(Complex src_complex, Complex dst_complex)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	
	Vec2d src_uv(src_complex.real(), src_complex.imag());
	Vec2d dst_uv(dst_complex.real(), dst_complex.imag());
	Vec2d diff = src_uv - dst_uv;
	double distance = HyperbolicDistance(src_complex, dst_complex);

	double f = 1 + 2 * pow(diff.norm(), 2) / ((1 - pow(src_uv.norm(), 2))*(1 - pow(dst_uv.norm(), 2)));
	double darccosh = 1 / sqrt(f * f - 1);
	double constant_df = 4. / ((1 - pow(dst_uv.norm(), 2))*(1 - pow(src_uv.norm(), 2)));
	Vec2d vector_df = diff + pow(diff.norm(), 2) * src_uv / (1 - pow(src_uv.norm(), 2));
	return distance * darccosh * constant_df * vector_df;
}

OpenMesh::Vec2d HyperbolicOrbifoldSolver::ComputeGradient(OpenMesh::VertexHandle v)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	if (mesh.data(v).is_singularity()) {
		return Vec2d(0, 0);
	}
	Vec2d v_uv = mesh.texcoord2D(v);
	Complex v_complex(v_uv[0], v_uv[1]);

	Vec2d gradient(0, 0);

	for (SurfaceMesh::VertexOHalfedgeIter vohiter = mesh.voh_iter(v); vohiter.is_valid(); ++vohiter) {
		HalfedgeHandle h = *vohiter;
		VertexHandle neighbor = mesh.to_vertex_handle(h);
		auto neighbor_uv = mesh.texcoord2D(neighbor);
		Complex neighbor_complex(neighbor_uv[0], neighbor_uv[1]);
		double n_w = mesh.data(h).weight();
		assert(n_w > 0);
		gradient += n_w * ComputeGradientOfDistance2(v_complex, neighbor_complex);
	}

	if (!mesh.is_boundary(v)) return gradient;

	VertexHandle equiv = mesh.data(v).equivalent_vertex();

	for (SurfaceMesh::VertexOHalfedgeIter vohiter = mesh.voh_iter(equiv); vohiter.is_valid(); ++vohiter) {
		HalfedgeHandle h = *vohiter;
		VertexHandle neighbor = mesh.to_vertex_handle(h);
		auto neighbor_uv = mesh.texcoord2D(neighbor);
		Complex neighbor_complex(neighbor_uv[0], neighbor_uv[1]);
		std::function<Complex(Complex const)> transformation = mesh.property(vtx_transit_, equiv);
		neighbor_complex = transformation(neighbor_complex);
		double n_w = mesh.data(h).weight();
		gradient += n_w * ComputeGradientOfDistance2(v_complex, neighbor_complex);
	}
	double metric_factor = pow(1 - pow(v_uv.norm(), 2), 2) / 4.;
	return gradient * metric_factor;
}

OpenMesh::Vec2d HyperbolicOrbifoldSolver::ComputeIntrinsicGradient(OpenMesh::VertexHandle v)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	if (mesh.data(v).is_singularity()) {
		return Vec2d(0, 0);
	}
	Vec2d v_uv = mesh.texcoord2D(v);
	Complex v_complex(v_uv[0], v_uv[1]);

	Complex gradient(0, 0);
	
	for (SurfaceMesh::VertexOHalfedgeIter vohiter = mesh.voh_iter(v); vohiter.is_valid(); ++vohiter) {
		HalfedgeHandle h = *vohiter;
		VertexHandle neighbor = mesh.to_vertex_handle(h);
		auto neighbor_uv = mesh.texcoord2D(neighbor);
		Complex neighbor_complex(neighbor_uv[0], neighbor_uv[1]);
		double n_w = mesh.data(h).weight();
		assert(n_w > 0);
		gradient -= n_w * InverseExponentialMap(v_complex, neighbor_complex);
	}

	if (!mesh.is_boundary(v)) return Vec2d(gradient.real(), gradient.imag());

	VertexHandle equiv = mesh.data(v).equivalent_vertex();
	
	for (SurfaceMesh::VertexOHalfedgeIter vohiter = mesh.voh_iter(equiv); vohiter.is_valid(); ++vohiter) {
		HalfedgeHandle h = *vohiter;
		VertexHandle neighbor = mesh.to_vertex_handle(h);
		auto neighbor_uv = mesh.texcoord2D(neighbor);
		Complex neighbor_complex(neighbor_uv[0], neighbor_uv[1]);
		std::function<Complex(Complex const)> transformation = mesh.property(vtx_transit_, equiv);
		neighbor_complex = transformation(neighbor_complex);
		double n_w = mesh.data(h).weight();
		gradient -= n_w * InverseExponentialMap(v_complex, neighbor_complex);
	}

	return Vec2d(gradient.real(), gradient.imag());
}

double HyperbolicOrbifoldSolver::ComputeEnergy()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	double energy = 0;
	for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); ++hiter) {
		HalfedgeHandle h = *hiter;
		VertexHandle v = mesh.from_vertex_handle(h);
		VertexHandle tv = mesh.to_vertex_handle(h);
		Vec2d v_uv = mesh.texcoord2D(v);
		Vec2d tv_uv = mesh.texcoord2D(tv);
		Complex v_complex(v_uv[0], v_uv[1]);
		Complex tv_complex(tv_uv[0], tv_uv[1]);
		double weight = mesh.data(h).weight();
		energy += weight * pow(HyperbolicDistance(v_complex, tv_complex), 2);
	}
	return energy * 0.5;
}

Eigen::VectorXd HyperbolicOrbifoldSolver::GetCoordsVector()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	Eigen::VectorXd uv_vector(mesh.n_vertices() * 2);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec2d uv = mesh.texcoord2D(v);
		uv_vector(v.idx() * 2) = uv[0];
		uv_vector(v.idx() * 2 + 1) = uv[1];
	}
	return uv_vector;
}

Eigen::VectorXd HyperbolicOrbifoldSolver::GetGradientVector()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	Eigen::VectorXd gradient_vector(mesh.n_vertices() * 2);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec2d grad = mesh.data(v).gradient();
		gradient_vector(v.idx() * 2) = grad[0];
		gradient_vector(v.idx() * 2 + 1) = grad[1];
	}
	return gradient_vector;
}

void HyperbolicOrbifoldSolver::SetCoords(const Eigen::VectorXd & uv_vector)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec2d uv(uv_vector(v.idx() * 2), uv_vector(v.idx() * 2 + 1));
		mesh.set_texcoord2D(v, uv);
	}
}

double HyperbolicOrbifoldSolver::OptimizationLoop(double step_length, double error)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	//double energy_prev = 1e10;
	//double energy = ComputeEnergy();
	int epoch = 0;
	double gradient_length = 1e10;
	while (gradient_length > error) {
		gradient_length = ComputeGradient();
		if (epoch % 20) {
			std::cout << "Max Gradient Length:" << gradient_length << "\tEnergy:" << ComputeEnergy() << std::endl;
		}
		
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
			VertexHandle v = *viter;
			Vec2d uv = mesh.texcoord2D(v);
			Vec2d gradient = mesh.data(v).gradient();
			
			Complex uv_complex(uv[0], uv[1]);
			uv -= step_length * gradient;
			assert(uv.norm() < 1);
			mesh.set_texcoord2D(v, uv);
		}

		Normalize();
		//energy_prev = energy;
		//energy = ComputeEnergy();
		epoch++;
	}
	return gradient_length;
}

void HyperbolicOrbifoldSolver::Normalize()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	for (auto it = segments_vts_.begin(); it != segments_vts_.end(); ++it) {
		for (auto viter = (*it).begin(); viter != (*it).end(); ++viter) {
			VertexHandle v = *viter;
			VertexHandle equiv = mesh.data(v).equivalent_vertex();
			Vec2d equiv_uv = mesh.texcoord2D(equiv);
			Complex equiv_complex(equiv_uv[0], equiv_uv[1]);
			Complex v_complex = mesh.property(vtx_transit_, equiv)(equiv_complex);
			Vec2d v_uv(v_complex.real(), v_complex.imag());
			mesh.set_texcoord2D(v, v_uv);
		}
	}
}

void HyperbolicOrbifoldSolver::InitMap()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	using namespace Eigen;
	
	SparseMatrix<double> A(mesh.n_vertices() * 2, mesh.n_vertices() * 2);
	A.setZero();
	VectorXd b(mesh.n_vertices() * 2);

	std::vector<Triplet<double>> A_coefficients;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		int idx = v.idx();
		if (mesh.is_boundary(v)) {
			A_coefficients.push_back(Triplet<double>(2 * idx, 2 * idx, 1.));
			A_coefficients.push_back(Triplet<double>(2 * idx + 1, 2 * idx + 1, 1.));
			auto uv = mesh.texcoord2D(v);
			b(2 * idx) = uv[0];
			b(2 * idx + 1) = uv[1];
		}
		else {
			b(2 * idx) = 0;
			b(2 * idx + 1) = 0;
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
	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
	
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}
	Eigen::VectorXd x = solver.solve(b);
	//std::cout << "Error:" << (A * x - b).norm() << std::endl;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh.is_boundary(v)) continue;
		Vec2d uv(x(2 * v.idx()), x(2 * v.idx() + 1));
		mesh.set_texcoord2D(v, uv);
	}
}


