#include "MeshDefinition.h"

SurfaceMesh::SurfaceMesh()
{

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
