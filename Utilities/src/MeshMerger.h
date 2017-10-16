#ifndef MESH_MERGER_H
#define MESH_MERGER_H

#include <Eigen\Core>
void MergeMeshMatrix(Eigen::MatrixXd & V1, Eigen::MatrixXi & T1, Eigen::MatrixXd &TC1,
	Eigen::MatrixXd & V2, Eigen::MatrixXi & T2, Eigen::MatrixXd &TC2,
	Eigen::MatrixXd & OV, Eigen::MatrixXi & OT, Eigen::MatrixXd &OTC);


#endif // !MESH_MATRIX_MERGER_H
