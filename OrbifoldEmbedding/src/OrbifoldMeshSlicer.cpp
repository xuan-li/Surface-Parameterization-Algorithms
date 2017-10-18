#include "OrbifoldMeshSlicer.h"

OrbifoldMeshSlicer::OrbifoldMeshSlicer(SurfaceMesh & mesh):mesh_(mesh)
{

}

SurfaceMesh OrbifoldMeshSlicer::CutAndSelectSingularities()
{
	MeshSlicer slicer(mesh_);
	SurfaceMesh sliced_mesh = slicer.SliceMeshToDisk();
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		OpenMesh::VertexHandle v = *viter;
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
	OpenMesh::VertexHandle origin = longest_path.front();
	OpenMesh::VertexHandle middle = longest_path[longest_path.size() / 2];
	OpenMesh::VertexHandle end = longest_path.back();

	cone_vertices_.push_back(slicer.SplitTo(origin).front());
	cone_vertices_.push_back(slicer.SplitTo(middle).front());
	cone_vertices_.push_back(slicer.SplitTo(middle).back());
	cone_vertices_.push_back(slicer.SplitTo(end).front());

	for (auto it = cone_vertices_.begin(); it != cone_vertices_.end(); ++it) {
		OpenMesh::VertexHandle v = *it;
		sliced_mesh.data(v).set_singularity(true);
	}

	return sliced_mesh;
}
