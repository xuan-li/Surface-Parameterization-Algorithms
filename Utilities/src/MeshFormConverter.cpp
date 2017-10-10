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
	UV.resize(nv, 3);
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
