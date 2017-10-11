#include "Dijkstra.h"
#include <list>

#ifndef INF
#define INF 0x3f3f3f3f
#endif

void DijkstraShortestDist(SurfaceMesh & mesh, OpenMesh::VertexHandle src, OpenMesh::VPropHandleT<double>& dist, OpenMesh::VPropHandleT<OpenMesh::VertexHandle>& parent)
{
	// check whether dist property is registered;
	if (!(dist.is_valid())) {
		mesh.add_property(dist);
	}

	if (!(parent.is_valid())) {
		mesh.add_property(parent);
	}


	OpenMesh::VPropHandleT<bool> touched;
	mesh.add_property(touched);


	printf("Calculate distence from vertex %d\r\n", src.idx());
	std::list<OpenMesh::VertexHandle> S;
	std::list<OpenMesh::VertexHandle> U;

	// initiate dist of all vertices as INF;
	for (SurfaceMesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		OpenMesh::VertexHandle v = *viter;
		mesh.property(dist, v) = INF;
		mesh.property(touched, v) = false;
	}

	// set dist of root as 0
	mesh.property(dist, src) = 0;
	S.push_back(src);
	mesh.property(touched, src) = true;

	// update the neighbors of the root
	for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(src); vviter.is_valid(); ++vviter) {
		OpenMesh::VertexHandle v = *vviter;
		mesh.property(dist, v) = mesh.calc_edge_length(mesh.edge_handle(mesh.find_halfedge(src, v)));
		mesh.property(parent, v) = src;
		U.push_back(v);
	}

	while (U.size() > 0)
	{
		OpenMesh::VertexHandle v_to_update;
		double shortest = 2 * INF;
		std::list<OpenMesh::VertexHandle>::iterator iv;
		for (std::list<OpenMesh::VertexHandle>::iterator it = U.begin(); it != U.end(); ++it)
		{
			OpenMesh::VertexHandle v = *it;
			if (mesh.property(dist, v) < shortest) {
				shortest = mesh.property(dist, v);
				v_to_update = v;
				iv = it;
			}
		}
		S.push_back(v_to_update);
		U.erase(iv);
		mesh.property(touched, v_to_update) = true;
		for (SurfaceMesh::VertexVertexIter vviter = mesh.vv_iter(v_to_update); vviter.is_valid(); ++vviter) {
			OpenMesh::VertexHandle v = *vviter;
			if (!mesh.property(touched, v)) {
				if (mesh.property(dist, v_to_update) + mesh.calc_edge_length(mesh.edge_handle(mesh.find_halfedge(v_to_update, v))) < mesh.property(dist, v)) {
					mesh.property(parent, v) = v_to_update;
					mesh.property(dist, v) = mesh.property(dist, v_to_update) + mesh.calc_edge_length(mesh.edge_handle(mesh.find_halfedge(v_to_update, v)));
					U.push_back(v);
				}
			}
		}
	}
}
