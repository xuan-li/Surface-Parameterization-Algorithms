#include "MeshSlicer.h"
#include <list>
#include <map>

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
	mesh_.add_property(split_to_);




}

void MeshSlicer::SliceMeshToDisk(SurfaceMesh &slice_mesh)
{
	int euler = mesh_.n_vertices() - mesh_.n_edges() + mesh_.n_faces();
	if (euler == 2)
		FindAndMarkCutGraphSphere();
	else
		FindAndMarkCutGraphNonSphere();
	ConstructWedge();
	SliceAccordingToWedge(slice_mesh);
}

std::vector<OpenMesh::VertexHandle> MeshSlicer::SplitTo(OpenMesh::VertexHandle v)
{
	return mesh_.property(split_to_,v);
}

std::vector<OpenMesh::VertexHandle> MeshSlicer::GetLongestPath()
{
	return longest_path_;
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
	
	longest_path_ = longest_path_vertices;
	for (auto eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		OpenMesh::EdgeHandle e = *eiter;
		mesh_.property(on_cut_, e) = false;
	}
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
	for (auto hiter = mesh_.halfedges_begin(); hiter != mesh_.halfedges_end(); ++hiter) {
		OpenMesh::HalfedgeHandle h = *hiter;
		mesh_.property(wedge_, h) = 0;
	}
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
		int w = 0;
		std::list<HalfedgeHandle>::iterator it = halfedges_around.begin();
		mesh_.property(wedge_, *it) = 0;
		++it;
		for (; it != halfedges_around.end(); ++it) {
			HalfedgeHandle h = *it;
			if (mesh_.property(on_cut_, mesh_.edge_handle(h))
				|| mesh_.is_boundary(mesh_.edge_handle(h))) {
				w++;
			}
			mesh_.property(wedge_, h) = w;
		}
	}
}

void MeshSlicer::SliceAccordingToWedge(SurfaceMesh &new_mesh)
{
	using namespace OpenMesh;
	HPropHandleT<VertexHandle> new_end; // store the new end of an halfedge after being cutted.
	mesh_.add_property(new_end);

	HPropHandleT<HalfedgeHandle> halfedge_split_to;
	mesh_.add_property(halfedge_split_to);

	int max_vid = mesh_.n_vertices() - 1;

	for (SurfaceMesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		// triverse around in-halfedges, add new vertices if needed and bind halfedges to their new ends.
		std::map<int, VertexHandle> whether_new_vert_exists;
		for (SurfaceMesh::VertexIHalfedgeIter vihiter = mesh_.vih_iter(v); vihiter.is_valid(); ++vihiter) {
			HalfedgeHandle h = *vihiter;
			int wedge = mesh_.property(wedge_, h);
			if (whether_new_vert_exists[wedge].is_valid())
				mesh_.property(new_end, h) = whether_new_vert_exists[wedge];
			else {
				VertexHandle new_vertex = new_mesh.add_vertex(mesh_.point(v));
				mesh_.property(split_to_, v).push_back(new_vertex);
				mesh_.property(new_end, h) = new_vertex;
				whether_new_vert_exists[wedge] = new_vertex;
				TransferVertexData(v, new_vertex);
			}
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
			VertexHandle nv = mesh_.property(new_end, h);
			verts.push_back(nv);
		}
		FaceHandle new_f = new_mesh.add_face(verts);
		TransferFaceData(f, new_f);
		for (int i = 0; i < 3; ++i) {
			HalfedgeHandle he = new_mesh.find_halfedge(verts[i], verts[(i + 1) % 3]);
			HalfedgeHandle ohe = mesh_.find_halfedge(old_verts[i], old_verts[(i + 1) % 3]);
			mesh_.property(halfedge_split_to, ohe) = he;
			TransferHalfedgeData(ohe, he);
			EdgeHandle e = new_mesh.edge_handle(he);
			EdgeHandle oe = mesh_.edge_handle(ohe);
			TransferEdgeData(oe, e);
		}
	}
	new_mesh.RequestBoundary();
	for (auto eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		EdgeHandle e = *eiter;
		HalfedgeHandle h0 = mesh_.halfedge_handle(e, 0);
		HalfedgeHandle h1 = mesh_.halfedge_handle(e, 1);
		HalfedgeHandle new_h0 = mesh_.property(halfedge_split_to, h0);
		HalfedgeHandle new_h1 = mesh_.property(halfedge_split_to, h1);
		if (mesh_.is_boundary(e)) {
			assert(new_h0.is_valid());
			assert(new_h1.is_valid());
		}
		new_mesh.data(new_h0).set_original_opposition(new_h1);
		new_mesh.data(new_h1).set_original_opposition(new_h0);
	}
}