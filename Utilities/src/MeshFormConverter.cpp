#include "MeshFormConverter.h"

typedef SurfaceMesh Mesh;

void OpenMeshToMatrix(SurfaceMesh & mesh, Eigen::MatrixXd & V, Eigen::MatrixXi & F)
{
	int nv = mesh.n_vertices();
	int nf = mesh.n_faces();

	V.setZero();
	V.resize(nv, 3);

	/*load vertex data*/
	for (Mesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		Mesh::VertexHandle v = *viter;

		int idx = v.idx();
		Mesh::Point point = mesh.point(v);

		for (int i = 0; i < 3; i++) {
			V(idx, i) = point[i];
		}
	}

	F.resize(nf, 3);
	F.setConstant(-1);
	for (Mesh::FaceIter fiter = mesh.faces_begin(); fiter != mesh.faces_end(); ++fiter) {
		Mesh::FaceHandle f = *fiter;
		int i = 0;
		for (Mesh::FaceVertexIter fviter = mesh.fv_iter(f); fviter.is_valid(); ++fviter) {
			Mesh::VertexHandle v = *fviter;
			F(f.idx(), i) = v.idx();
			i++;
		}
	}
}



void OpenMeshCoordToMatrix(SurfaceMesh & mesh, Eigen::MatrixXd &UV)
{
	int nv = mesh.n_vertices();
	int nf = mesh.n_faces();

	UV.setZero();
	UV.resize(nv, 2);
	/*load vertex data*/
	for (Mesh::VertexIter viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		Mesh::VertexHandle v = *viter;
		int idx = v.idx();
		Mesh::TexCoord2D point = mesh.texcoord2D(v);

		for (int i = 0; i < 2; i++) {
			UV(idx, i) = point[i];
		}
	}
}

void MatrixToOpenMesh(Eigen::MatrixXd & V, Eigen::MatrixXi & F, SurfaceMesh & mesh)
{

}

void HalfedgesToMatrix(SurfaceMesh & mesh, std::vector<OpenMesh::HalfedgeHandle> halfedges, Eigen::MatrixXd & P1, Eigen::MatrixXd & P2)
{
	int n_halfedges = halfedges.size();
	P1.resize(n_halfedges, 3);
	P2.resize(n_halfedges, 3);
	for (int i = 0; i < n_halfedges; i++) {
		OpenMesh::HalfedgeHandle h = halfedges[i];
		OpenMesh::VertexHandle v1 = mesh.from_vertex_handle(h);
		OpenMesh::VertexHandle v2 = mesh.to_vertex_handle(h);
		OpenMesh::Vec3d p1 = mesh.point(v1);
		OpenMesh::Vec3d p2 = mesh.point(v2);
		P1.row(i) = Eigen::RowVector3d(p1[0], p1[1], p1[2]);
		P2.row(i) = Eigen::RowVector3d(p2[0], p2[1], p2[2]);
	}
}
