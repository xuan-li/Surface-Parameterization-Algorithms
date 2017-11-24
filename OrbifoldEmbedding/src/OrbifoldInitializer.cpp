#include "OrbifoldInitializer.h"
#include <list>

OrbifoldInitializer::OrbifoldInitializer(SurfaceMesh & mesh):mesh_(mesh)
{

}

void OrbifoldInitializer::Initiate(SurfaceMesh &sliced_mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = mesh_;
	cone_flag_ = cone_flag;
	cone_angle_ = cone_angle;
	slice_flag_ = slice_flag;
	CutMesh(sliced_mesh);
	CutBoundaryToSegments(sliced_mesh);
	
}

void OrbifoldInitializer::ComputeEuclideanTransformations(SurfaceMesh & sliced_mesh, OpenMesh::VPropHandleT<Eigen::Matrix3d>& vtx_transit)
{
	if (!vtx_transit.is_valid())
		sliced_mesh.add_property(vtx_transit);
	using namespace OpenMesh;
	InitiateEConeCoords(sliced_mesh);

	for (auto it = segments_vts_.begin(); it != segments_vts_.end(); ++it) {
		auto seg = *it;
		VertexHandle vs0 = seg.front();
		VertexHandle vs1 = seg.back();
		VertexHandle vt0 = sliced_mesh.data(vs0).equivalent_vertex();
		VertexHandle vt1 = sliced_mesh.data(vs1).equivalent_vertex();
		if (!vt0.is_valid())
			vt0 = vs0;
		if (!vt1.is_valid())
			vt1 = vs1;
		Complex s0(sliced_mesh.texcoord2D(vs0)[0], sliced_mesh.texcoord2D(vs0)[1]);
		Complex s1(sliced_mesh.texcoord2D(vs1)[0], sliced_mesh.texcoord2D(vs1)[1]);
		Complex t0(sliced_mesh.texcoord2D(vt0)[0], sliced_mesh.texcoord2D(vt0)[1]);
		Complex t1(sliced_mesh.texcoord2D(vt1)[0], sliced_mesh.texcoord2D(vt1)[1]);
		Eigen::MatrixXd T = ComputeHomogeousRigidTransformation(s0, s1, t0, t1);
		for (auto viter = seg.begin(); viter != seg.end(); ++viter) {
			VertexHandle v = *viter;
			if (sliced_mesh.data(v).is_singularity()) continue;
			sliced_mesh.property(vtx_transit, v) = T;
		}
	}
	
}

void OrbifoldInitializer::ComputeHyperbolicTransformations(SurfaceMesh & sliced_mesh, OpenMesh::VPropHandleT<std::function<Complex(Complex const)>>& vtx_transit)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh;
	if (!vtx_transit.is_valid())
		sliced_mesh.add_property(vtx_transit);

	InitiateHConeCoords(sliced_mesh);
	
	for (auto it = segments_vts_.begin(); it != segments_vts_.end(); ++it) {
		auto seg = *it;
		VertexHandle vs0 = seg.front();
		VertexHandle vs1 = seg.back();
		VertexHandle vt0 = sliced_mesh.data(vs0).equivalent_vertex();
		VertexHandle vt1 = sliced_mesh.data(vs1).equivalent_vertex();
		if (!vt0.is_valid())
			vt0 = vs0;
		if (!vt1.is_valid())
			vt1 = vs1;
		Complex s0(sliced_mesh.texcoord2D(vs0)[0], sliced_mesh.texcoord2D(vs0)[1]);
		Complex s1(sliced_mesh.texcoord2D(vs1)[0], sliced_mesh.texcoord2D(vs1)[1]);
		Complex t0(sliced_mesh.texcoord2D(vt0)[0], sliced_mesh.texcoord2D(vt0)[1]);
		Complex t1(sliced_mesh.texcoord2D(vt1)[0], sliced_mesh.texcoord2D(vt1)[1]);
	

		auto transformation = ComputeMobiusTransformation(s0, s1, t0, t1);
		assert(abs(transformation(s1) - t1) < 1e-6);

		for (auto vit = seg.begin(); vit != seg.end(); ++vit) {
			if(!mesh.data(*vit).is_singularity())
				mesh.property(vtx_transit, *vit) = transformation;
		}
	}
}

void OrbifoldInitializer::InitiateEConeCoords(SurfaceMesh & sliced_mesh)
{
	if (abs(sliced_mesh.data(cone_vertices_.front()).angle_sum() - 0.5 * PI) < 1e-5 &&
		cone_vertices_.size() == 4) {
		InitiateEConeCoordsType1(sliced_mesh);
	}
}

void OrbifoldInitializer::InitiateHConeCoords(SurfaceMesh & sliced_mesh)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh;
	int n_cones = (cone_vertices_.size()+2) / 2;

	double radius = HyperbolicCosineLaw(2 * PI / n_cones, PI / 4., PI / 4.);
	double ratio = (exp(radius) - 1) / (exp(radius) + 1);

	Vec2d p1(cos(PI / 2 + (n_cones / 2) * 2 * PI / n_cones)*ratio, sin(PI / 2 + (n_cones / 2) * 2 * PI / n_cones)*ratio);
	Vec2d pk(cos(PI / 2 + (1 + n_cones / 2) * 2 * PI / n_cones)*ratio, sin(PI / 2 + (1 + n_cones / 2) * 2 * PI / n_cones)*ratio);
	double edge_length = HyperbolicCosineLaw(PI / 4, PI / 4, 2 * PI / n_cones);
	double target_radis = (exp(edge_length / 2) - 1) / (exp(edge_length / 2) + 1);
	Complex p1_source(p1[0], p1[1]);
	Complex pk_source(pk[0], pk[1]);
	Complex p1_target(-target_radis, 0);
	Complex pk_target(target_radis, 0);

	auto transformation = ComputeMobiusTransformation(p1_source, pk_source, p1_target, pk_target);
	for (int i = 0; i < n_cones; ++i) {
		Vec2d uv = Vec2d(cos(PI / 2 + (-i + n_cones / 2) * 2 * PI / n_cones)*ratio, sin(PI / 2 + (-i + n_cones / 2) * 2 * PI / n_cones)*ratio);
		Complex uv_complex = transformation(Complex(uv[0], uv[1]));
		uv = Vec2d(uv_complex.real(), uv_complex.imag());
		mesh.set_texcoord2D(cone_vertices_[i], uv);
		uv_complex = std::conj(uv_complex);
		uv = Vec2d(uv_complex.real(), uv_complex.imag());
		if (i != 0 && i != n_cones - 1)
			mesh.set_texcoord2D(mesh.data(cone_vertices_[i]).equivalent_vertex(), uv);
	}
	std::cout << "Cone coordinates:\n";
	for (int i = 0; i < cone_vertices_.size(); ++i) {
		Vec2d uv = mesh.texcoord2D(cone_vertices_[i]);
		std::cout << uv[0] << "\t" << uv[1] << std::endl;
	}
}

double OrbifoldInitializer::HyperbolicCosineLaw(double a, double b, double c)
{
	double l = (cos(c) + cos(a)*cos(b)) / (sin(a) * sin(b));
	return acosh(l);
}


void OrbifoldInitializer::CutMesh(SurfaceMesh & sliced_mesh)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = mesh_;

	MeshSlicer slicer(mesh);
	slicer.ResetFlags();
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		EdgeHandle e = *eiter;
		if (mesh.property(slice_flag_, e)) {
			slicer.AddOnCutEdge(e);
		}
	}

	slicer.ConstructWedge();
	slicer.SliceAccordingToWedge(sliced_mesh);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		auto verts = slicer.SplitTo(v);
		if (verts.size() == 2) {
			sliced_mesh.data(verts[0]).set_equivalent_vertex(verts[1]);
			sliced_mesh.data(verts[1]).set_equivalent_vertex(verts[0]);
		}
		if (mesh.property(cone_flag_, v)) {
			for (auto it = verts.begin(); it != verts.end(); ++it) {
				sliced_mesh.data(*it).set_singularity(true);
				sliced_mesh.data(*it).set_angle_sum(mesh.property(cone_angle_, v) / verts.size());
			}
		}
	}
	
}

void OrbifoldInitializer::CutBoundaryToSegments(SurfaceMesh &sliced_mesh)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh;

	mesh.RequestBoundary();
	auto boundary = mesh.GetBoundaries().front();

	std::list<HalfedgeHandle> boundary_list(boundary.begin(), boundary.end());
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		HalfedgeHandle h = *it;
		if (mesh.data(mesh.from_vertex_handle(h)).is_singularity()  && !mesh.data(mesh.from_vertex_handle(h)).equivalent_vertex().is_valid()) {
			boundary_list.insert(boundary_list.end(), boundary_list.begin(), it);
			boundary_list.erase(boundary_list.begin(), it);
			break;
		}
	}

	segments_vts_.clear();
	std::vector<OpenMesh::VertexHandle> segment;
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		HalfedgeHandle h = *it;
		segment.push_back(mesh.from_vertex_handle(h));
		if (mesh.data(mesh.from_vertex_handle(h)).is_singularity()) {
			cone_vertices_.push_back(mesh.from_vertex_handle(h));
		}
		if (mesh.data(mesh.to_vertex_handle(h)).is_singularity()) {
			segment.push_back(mesh.to_vertex_handle(h));
			segments_vts_.push_back(segment);
			segment.clear();
		}
	}
}

void OrbifoldInitializer::InitiateEConeCoordsType1(SurfaceMesh & sliced_mesh)
{
	using namespace OpenMesh;
	OpenMesh::Vec2d v0(0, 0);
	OpenMesh::Vec2d v1(0, 1);
	OpenMesh::Vec2d v2(1, 1);
	OpenMesh::Vec2d v3(1, 0);
	sliced_mesh.set_texcoord2D(cone_vertices_[0], v0);
	sliced_mesh.set_texcoord2D(cone_vertices_[1], v1);
	sliced_mesh.set_texcoord2D(cone_vertices_[2], v2);
	sliced_mesh.set_texcoord2D(cone_vertices_[3], v3);
}

void OrbifoldInitializer::InitiateEConeCoordsType2(SurfaceMesh & sliced_mesh)
{
}

void OrbifoldInitializer::InitiateEConeCoordsType3(SurfaceMesh & sliced_mesh)
{
}

void OrbifoldInitializer::InitiateEConeCoordsType4(SurfaceMesh & sliced_mesh)
{
}
