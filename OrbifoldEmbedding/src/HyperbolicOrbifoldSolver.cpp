#include "HyperbolicOrbifoldSolver.h"


HyperbolicOrbifoldSolver::HyperbolicOrbifoldSolver(SurfaceMesh & mesh): mesh_(mesh)
{

}

SurfaceMesh HyperbolicOrbifoldSolver::Compute()
{
	return SurfaceMesh();
}

void HyperbolicOrbifoldSolver::CutToDist()
{
	OrbifoldMeshSlicer slicer(mesh_);
	sliced_mesh_ = slicer.CutAndSelectSingularities();
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

	CutToDist();


}

double HyperbolicOrbifoldSolver::CosineLaw(double a, double b, double c)
{
	double cs = (cosh(a) * cosh(b) - cosh(c)) / (sinh(a) * sinh(b));
	return acos(cs);

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


