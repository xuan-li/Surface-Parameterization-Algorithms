#include "BFFInitializer.h"

BFFInitializer::BFFInitializer(SurfaceMesh & mesh)
: mesh_(mesh)
{

}

void BFFInitializer::Initiate(SurfaceMesh & sliced_mesh, OpenMesh::VPropHandleT<bool> cone_flag, OpenMesh::VPropHandleT<double> cone_angle, OpenMesh::EPropHandleT<bool> slice_flag)
{
	cone_flag_ = cone_flag;
	cone_angle_ = cone_angle;
	slice_flag_ = slice_flag;
	CutMesh(sliced_mesh);
}

void BFFInitializer::CutMesh(SurfaceMesh & sliced_mesh)
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

	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		EdgeHandle e = *eiter;
		if (mesh.is_boundary(e)) continue;
		HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
		HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);
		HalfedgeHandle h0_to = slicer.ConvertTo(h0);
		HalfedgeHandle h1_to = slicer.ConvertTo(h1);
		sliced_mesh.data(h0_to).set_original_opposition(h1_to);
		sliced_mesh.data(h1_to).set_original_opposition(h0_to);
	}

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle  v = *viter;
		if (mesh.property(cone_flag_, v)) {
			auto verts = slicer.SplitTo(v);
			for (auto it = verts.begin(); it != verts.end(); ++it) {
				sliced_mesh.data(*it).set_singularity(true);
			}
		}
	}


	split_to_ = slicer.split_to();

	sliced_mesh.RequestBoundary();
	auto boundary = sliced_mesh.GetBoundaries().front();
	for (auto it = boundary.begin(); it != boundary.end(); ++it) {
		HalfedgeHandle h = *it;
		if (sliced_mesh.data(sliced_mesh.from_vertex_handle(h)).is_singularity()) {
			cone_vertices_.push_back(sliced_mesh.from_vertex_handle(h));
		}
	}

}

