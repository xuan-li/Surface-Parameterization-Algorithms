#include "MeshSlicer.h"
#include <list>


MeshSlicer::MeshSlicer(
	SurfaceMesh &mesh, 
	std::function<void(OpenMesh::VertexHandle, OpenMesh::VertexHandle)> f1,
	std::function<void(OpenMesh::FaceHandle, OpenMesh::FaceHandle)> f2,
	std::function<void(OpenMesh::EdgeHandle, OpenMesh::EdgeHandle)> f3,
	std::function<void(OpenMesh::HalfedgeHandle, OpenMesh::HalfedgeHandle)> f4
): mesh_(mesh)
{
	TransferVertexData = f1;
	TransferFaceData = f2;
	TransferEdgeData = f3;
	TransferHalfedgeData = f4;

	mesh_.add_property(wedge_);
	mesh_.add_property(on_cut_);

	for (auto eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		OpenMesh::EdgeHandle e = *eiter;
		mesh_.property(on_cut_, e);
	}
	for (auto hiter = mesh_.halfedges_begin(); hiter != mesh_.halfedges_end(); ++hiter) {
		OpenMesh::HalfedgeHandle h = *hiter;
		mesh_.property(wedge_, h);
	}
}

SurfaceMesh MeshSlicer::SliceMeshToDisk()
{
	int euler = mesh_.n_vertices() - mesh_.n_edges() + mesh_.n_faces();
	if (euler == 2)
		FindAndMarkCutGraphSphere();
	else
		FindAndMarkCutGraphNonSphere();
	ConstructWedge();
	return SliceAccordingToWedge();
}

void MeshSlicer::FindAndMarkCutGraphSphere()
{
	using namespace OpenMesh;
	base_point_ = *(mesh_.vertices_begin());
	OpenMesh::VPropHandleT<double> dist;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> parent;
	DijkstraShortestDist(mesh_, base_point_, dist, parent);
	
	// Find the biggest dist
	double max_dist = 0;
	VertexHandle end;
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh_.property(dist, v) > max_dist) {
			max_dist = mesh_.property(dist, v);
			end = v;
		}
	}

	std::vector<VertexHandle> longest_path_vertices;

	VertexHandle cptr = end;
	longest_path_vertices.push_back(end);
	while (mesh_.property(parent, cptr).is_valid()) {
		cptr = mesh_.property(parent, cptr);
		longest_path_vertices.push_back(cptr);
	}

	std::reverse(longest_path_vertices.begin(), longest_path_vertices.end());
	
	for (int i = 0; i < longest_path_vertices.size() - 1; ++i) {
		EdgeHandle e = mesh_.edge_handle(mesh_.find_halfedge(longest_path_vertices[i], longest_path_vertices[i+1]));
		mesh_.property(on_cut_, e) = true;
	}
}

void MeshSlicer::FindAndMarkCutGraphNonSphere()
{

}

void MeshSlicer::ConstructWedge()
{
	using namespace OpenMesh;

	for (SurfaceMesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		std::list<HalfedgeHandle> halfedges_around;
		for (SurfaceMesh::VertexIHalfedgeCWIter vihiter = mesh_.vih_cwiter(v); vihiter.is_valid(); ++vihiter) {
			HalfedgeHandle h = *vihiter;
			if (mesh_.is_boundary(h)) continue;
			halfedges_around.push_back(h);
		}
		for (std::list<HalfedgeHandle>::iterator it = halfedges_around.begin(); it != halfedges_around.end(); ++it) {
			HalfedgeHandle h = *it;
			if (mesh_.is_boundary(v)) {
				if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h))) {
					halfedges_around.insert(halfedges_around.end(), halfedges_around.begin(), it);
					halfedges_around.erase(halfedges_around.begin(), it);
					break;
				}
			}
			else {
				if (mesh_.property(on_cut_,mesh_.edge_handle(h))) {
					halfedges_around.insert(halfedges_around.end(), halfedges_around.begin(), it);
					halfedges_around.erase(halfedges_around.begin(), it);
					break;
				}
			}
		}

		/*construct wedge around vertex v*/
		int w = -1;
		for (std::list<HalfedgeHandle>::iterator it = halfedges_around.begin(); it != halfedges_around.end(); ++it) {
			HalfedgeHandle h = *it;
			if (mesh_.property(on_cut_, mesh_.edge_handle(h))
				|| mesh_.is_boundary(mesh_.edge_handle(h))) {
				w++;
			}
			mesh_.property(wedge_, h) = w;
		}
	}
}

SurfaceMesh MeshSlicer::SliceAccordingToWedge()
{
	using namespace OpenMesh;

	HPropHandleT<VertexHandle> new_end; // store the new end of an halfedge after being cutted.
	mesh_.add_property(new_end);

	SurfaceMesh new_mesh;

	int max_vid = mesh_.n_vertices() - 1;

	for (SurfaceMesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		// triverse around in-halfedges, add new vertices if needed and bind halfedges to their new ends.
		for (int i = 0; i <= mesh_.property(max_wedge_, v); ++i) {
			VertexHandle nv = new_mesh.add_vertex(mesh_.point(v));
			new_mesh->Idb(nv) = mesh_.Id(v);
			TransferVertexProperties(v, nv);
		}
	}

	/*construct face*/
	for (SurfaceMesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		FaceHandle f = *fiter;
		std::vector<VertexHandle> verts;
		std::vector<VertexHandle> old_verts;
		for (SurfaceMesh::FaceHalfedgeIter fhiter = mesh_.fh_iter(f); fhiter.is_valid(); ++fhiter) {
			HalfedgeHandle h = *fhiter;
			VertexHandle v = mesh_.to_vertex_handle(h);
			old_verts.push_back(v);
			int new_vid = mesh_.Id(v) + max_vid * mesh_.Wedge(h);
			verts.push_back(new_mesh->IdVertex(new_vid));
		}
		FaceHandle new_f = new_mesh->CreateFace(verts, mesh_.Id(f));
		TransferFaceProperties(f, new_f);
		for (int i = 0; i < 3; ++i) {
			HalfedgeHandle he = new_mesh->find_halfedge(verts[i], verts[(i + 1) % 3]);
			HalfedgeHandle ohe = mesh_.find_halfedge(old_verts[i], old_verts[(i + 1) % 3]);
			TransferHalfedgeProperties(ohe, he);
			EdgeHandle e = new_mesh->edge_handle(he);
			EdgeHandle oe = mesh_.edge_handle(ohe);
			TransferEdgeProperties(oe, e);
		}
	}
	TransferMeshProperties(m_mesh, new_mesh);
	new_mesh->RequestBoundary();
	return new_mesh;
}
