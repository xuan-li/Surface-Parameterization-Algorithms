#ifndef MESH_FORM_CONVERTER
#define MESH_FORM_CONVERTER

#include <Eigen\Dense>
#include <MeshDefinition.h>

void OpenMeshToMatrix(SurfaceMesh &mesh, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
void OpenMeshToMatrix(SurfaceMesh &mesh, Eigen::MatrixXd &V, Eigen::MatrixXd &V_normal, Eigen::MatrixXi &F, Eigen::MatrixXd &F_normal);

void OpenMeshCoordToMatrix(SurfaceMesh &mesh, Eigen::MatrixXd &UV);

void MatrixToOpenMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, SurfaceMesh &mesh);

void HalfedgesToMatrix(SurfaceMesh &mesh, std::vector<OpenMesh::HalfedgeHandle> halfedges, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2);

#endif // !MESH_FORM_CONVERTER
