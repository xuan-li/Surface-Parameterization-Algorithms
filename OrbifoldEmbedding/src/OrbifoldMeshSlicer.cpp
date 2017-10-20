#include "OrbifoldMeshSlicer.h"
#include <list>

OrbifoldMeshSlicer::OrbifoldMeshSlicer(SurfaceMesh & mesh):mesh_(mesh)
{

}

SurfaceMesh OrbifoldMeshSlicer::CutAndSelectSingularities(int n_cones)
{
	MeshSlicer slicer(mesh_);
	SurfaceMesh sliced_mesh = slicer.SliceMeshToDisk();
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		OpenMesh::VertexHandle v = *viter;
		sliced_mesh.data(v).set_singularity(false);
		auto equivlent_splition = slicer.SplitTo(v);
		//assert(equivlent_splition.size() <= 2);
		if (equivlent_splition.size() == 2) {
			auto v1 = equivlent_splition[0];
			auto v2 = equivlent_splition[1];
			sliced_mesh.data(v1).set_equivalent_vertex(v2);
			sliced_mesh.data(v2).set_equivalent_vertex(v1);
		}
	}
	auto longest_path = slicer.GetLongestPath();


	std::vector<OpenMesh::VertexHandle> cone_vtx_original;
	cone_vtx_original.push_back(longest_path.front());
	for (int i = 1; i < n_cones - 1; ++i) {
		cone_vtx_original.push_back(longest_path[i*longest_path.size() / n_cones]);
	}
	cone_vtx_original.push_back(longest_path.back());

	for (auto it = cone_vtx_original.begin(); it != cone_vtx_original.end(); ++it) {
		auto vts = slicer.SplitTo(*it);
		for (auto it2 = vts.begin(); it2 != vts.end(); ++it2)
			cone_vertices_.push_back(*it2);
	}

	for (auto it = cone_vertices_.begin(); it != cone_vertices_.end(); ++it) {
		OpenMesh::VertexHandle v = *it;
		sliced_mesh.data(v).set_singularity(true);
	}

	sliced_mesh_ = sliced_mesh;

	CutBoundaryToSegments();

	return sliced_mesh;
}

void OrbifoldMeshSlicer::CutBoundaryToSegments()
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = sliced_mesh_;

	slice_root_vtx_ = cone_vertices_.front();
	cone_vertices_.erase(++cone_vertices_.begin(), cone_vertices_.end());
	assert(cone_vertices_.size() == 1);
	
	mesh.RequestBoundary();
	auto boundary = mesh.GetBoundaries().front();
	std::list<HalfedgeHandle> boundary_list(boundary.begin(), boundary.end());
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		HalfedgeHandle h = *it;
		if (mesh.from_vertex_handle(h) == slice_root_vtx_) {
			boundary_list.insert(boundary_list.end(), boundary_list.begin(), it);
			boundary_list.erase(boundary_list.begin(), it);
			break;
		}
	}

	segments_vts_.clear();
	std::vector<OpenMesh::VertexHandle> segment;
	for (auto it = boundary_list.begin(); it != boundary_list.end(); ++it) {
		HalfedgeHandle h = *it;
		if (mesh.data(mesh.to_vertex_handle(h)).is_singularity()) {
			cone_vertices_.push_back(mesh.to_vertex_handle(h));
			segments_vts_.push_back(segment);
			segment.clear();
			continue;
		}
		segment.push_back(mesh.to_vertex_handle(h));
	}
	cone_vertices_.pop_back();

}
