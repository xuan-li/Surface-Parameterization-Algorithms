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
	return sliced_mesh_;
}

void EuclideanOrbifoldSolver::CutToDist()
{
	OrbifoldMeshSlicer slicer(mesh_);
	sliced_mesh_ = slicer.CutAndSelectSingularities();
	cone_vts = slicer.GetConeVertices();

}


void EuclideanOrbifoldSolver::InitOrbifold()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;
	CutToDist();

	mesh.RequestBoundary();
	auto boundary = mesh.GetBoundaries().front();
	std::list<HalfedgeHandle> boundary_list(boundary.begin(), boundary.end());
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		HalfedgeHandle h = *it;
		if (mesh.from_vertex_handle(h) == cone_vts[0]) {
			boundary_list.insert(boundary_list.end(), boundary_list.begin(), it);
			boundary_list.erase(boundary_list.begin(), it);
			break;
		}
		
	}
	std::vector<OpenMesh::VertexHandle> segment;
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++i) {
		HalfedgeHandle h = *it;
		if (mesh.data(mesh.to_vertex_handle(h)).is_singularity()) {
			segments_vts_.push_back(segment);
			segment.push_back(mesh.to_vertex_handle(h));
			segment.clear();
		}
		segment.push_back(mesh.from_vertex_handle(h));
	}

	InitType1();
}

void EuclideanOrbifoldSolver::InitType1()
{
	
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
		HalfedgeHandle h0_next = mesh.next_halfedge_handle(h0);
		HalfedgeHandle h1_next = mesh.next_halfedge_handle(h1);

		double weight = 0.5 * (1 / tan(mesh.data(h0_next).angle()) + 1 / tan(mesh.data(h1_next).angle()));
		mesh.data(h0).set_weight(weight);
		mesh.data(h1).set_weight(weight);

	}
}
