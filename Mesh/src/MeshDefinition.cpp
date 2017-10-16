#include "MeshDefinition.h"

SurfaceMesh::SurfaceMesh()
{

}

void SurfaceMesh::RequestBoundary()
{
	OpenMesh::HPropHandleT<bool> touched;
	this->add_property(touched);
	boundaries.clear();
	bool has_boundary = false;
	for (HalfedgeIter hiter = halfedges_begin(); hiter != halfedges_end(); ++hiter) {
		HalfedgeHandle h = *hiter;
		property(touched, h) = false;
		if (is_boundary(h)) {
			has_boundary = true;
		}
	}
	if (!has_boundary) return;
	bool stop;
	do {
		stop = true;
		HalfedgeHandle hs, hc;
		for (HalfedgeIter hiter = halfedges_begin(); hiter != halfedges_end(); ++hiter) {
			HalfedgeHandle h = *hiter;
			if (is_boundary(h) && !property(touched, h)) {
				stop = false;
				hs = h;
				hc = next_halfedge_handle(h);
				break;
			}
		}
		if (hs.is_valid()) {
			std::vector<HalfedgeHandle> boundary;
			boundary.push_back(hs);
			property(touched, hs) = true;
			while (hc != hs) {
				boundary.push_back(hc);
				property(touched, hc) = true;
				hc = next_halfedge_handle(hc);
			}
			boundaries.push_back(boundary);
		}
	} while (!stop);
}

std::vector<std::vector<OpenMesh::HalfedgeHandle>> SurfaceMesh::GetBoundaries()
{
	return boundaries;
}

void NormalizeMesh(SurfaceMesh & mesh)
{
	using namespace OpenMesh;
	Vec3d s(0, 0, 0);
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		s += mesh.point(v);
	}
	s /= mesh.n_vertices();
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.set_point(v, mesh.point(v) - s);
	}

	double scale = 0;
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		scale = mesh.point(v).norm() > scale ? mesh.point(v).norm() : scale;
	}
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		mesh.set_point(v, mesh.point(v) / scale);
	}
}
