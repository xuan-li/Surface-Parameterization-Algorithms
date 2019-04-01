#include "MeshMarker.h"

void MeshMarker::ResetMarker()
{
	SurfaceMesh &mesh = *p_mesh_;
	using namespace OpenMesh;
	if (!singularity_.is_valid())
		mesh.add_property(singularity_);
	if (!slice_.is_valid())
		mesh.add_property(slice_);
	if (!cone_angle_.is_valid())
		mesh.add_property(cone_angle_);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		OpenMesh::VertexHandle v = *viter;
		mesh.property(singularity_, v) = false;
		mesh.property(cone_angle_, v) = 2 *PI;
	}
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		OpenMesh::EdgeHandle e = *eiter;
		mesh.property(slice_, e) = false;
	}
	n_edges_ = 0;
	n_vertices_ = 0;
}

// Compute the shortest path between two vertex and set slice flag.
void MeshMarker::ComputeAndSetSlice(OpenMesh::VertexHandle v0, OpenMesh::VertexHandle v1)
{
	SurfaceMesh &mesh = *p_mesh_;
	using namespace OpenMesh;
	VPropHandleT<double> dist;
	VPropHandleT<VertexHandle> parent;
	
	DijkstraShortestDist(mesh, v0, dist, parent);
	std::vector<VertexHandle> slice_vertices;

	VertexHandle cptr = v1;
	slice_vertices.push_back(v1);
	while (mesh.property(parent, cptr).is_valid()) {
		cptr = mesh.property(parent, cptr);
		slice_vertices.push_back(cptr);
	}

	std::reverse(slice_vertices.begin(), slice_vertices.end());

	for (int i = 0; i < slice_vertices.size() - 1; ++i) {
		EdgeHandle e = mesh.edge_handle(mesh.find_halfedge(slice_vertices[i], slice_vertices[i + 1]));
		if (!mesh.property(slice_, e)) {
			mesh.property(slice_, e) = true;
			++n_edges_;
		}
	}

}

// The cone angle is  cone_angle * PI;
void MeshMarker::SetSingularity(OpenMesh::VertexHandle v, double cone_angle)
{
	SurfaceMesh &mesh = *p_mesh_;
	using namespace OpenMesh;
	if (!mesh.property(singularity_, v)) {
		mesh.property(singularity_, v) = true;
		mesh.property(cone_angle_, v) = cone_angle * PI;
		++n_vertices_;
	}
	else {
		mesh.property(cone_angle_, v) = cone_angle * PI;
	}
}

void MeshMarker::GenerateMatrix(Eigen::MatrixXd & P, Eigen::MatrixXd & EP1, Eigen::MatrixXd & EP2)
{
	SurfaceMesh &mesh = *p_mesh_;
	using namespace OpenMesh;

	P.resize(n_vertices_, 3);
	EP1.resize(n_edges_, 3);
	EP2.resize(n_edges_, 3);

	int index = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		OpenMesh::VertexHandle v = *viter;
		if (mesh.property(singularity_, v)) {
			auto p = mesh.point(v);
			P.row(index) = Eigen::Vector3d(p[0], p[1], p[2]);
			++index;
		}
	}
	index = 0;
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		OpenMesh::EdgeHandle e = *eiter;
		if (mesh.property(slice_, e)) {
			HalfedgeHandle h = mesh.halfedge_handle(e, 0);
			VertexHandle v1 = mesh.from_vertex_handle(h);
			VertexHandle v2 = mesh.to_vertex_handle(h);
			auto p1 = mesh.point(v1);
			auto p2 = mesh.point(v2);
			EP1.row(index) = Eigen::Vector3d(p1[0], p1[1], p1[2]);
			EP2.row(index) = Eigen::Vector3d(p2[0], p2[1], p2[2]);
			++index;
		}
	}
}

void MeshMarker::LoadFromFile(std::string filename)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = *p_mesh_;
	ResetMarker();
	std::ifstream f(filename);
	std::string line;
	if (f.is_open())
	{
		while (!f.eof())
		{
			std::getline(f, line);
			std::vector<std::string> parser;
			Split(" ", line, parser);
			if (parser.size() == 0) continue;
			if (parser[0] == "v") {
				VertexHandle v = mesh.vertex_handle(stoi(parser[1]));
				mesh.property(singularity_, v) = true;
				mesh.property(cone_angle_, v) = atof(parser[2].c_str());
				++n_vertices_;
			}
			else if (parser[0] == "e") {
				EdgeHandle e = mesh.edge_handle(stoi(parser[1]));
				mesh.property(slice_, e) = true;
				++n_edges_;
			}
		}
		f.close();
	}
}

void MeshMarker::SaveToFile(std::string filename)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = *p_mesh_;

	std::ofstream f(filename);

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		if (mesh.property(singularity_, v)) {
			f << "v " << v.idx() << " " << std::setprecision(9) << mesh.property(cone_angle_, v) << "\n";
		}
	}
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); ++eiter) {
		EdgeHandle e = *eiter;
		if (mesh.property(slice_, e)) {
			f << "e " << e.idx() << "\n";
		}
	}

	f.close();
	
}
