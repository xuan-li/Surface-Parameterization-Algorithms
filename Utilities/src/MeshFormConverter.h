#ifndef MESH_FORM_CONVERTER
#define MESH_FORM_CONVERTER

#include <Eigen\Dense>
#include <MeshDefinition.h>

void OpenMeshToMatrix(SurfaceMesh &mesh, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

void MatrixToOpenMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, SurfaceMesh &mesh);


#endif // !MESH_FORM_CONVERTER
