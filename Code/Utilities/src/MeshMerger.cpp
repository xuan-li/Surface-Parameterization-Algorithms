#include "MeshMerger.h"

void MergeMeshMatrix(Eigen::MatrixXd  V1, Eigen::MatrixXi T1, Eigen::MatrixXd TC1,
	Eigen::MatrixXd  V2, Eigen::MatrixXi  T2, Eigen::MatrixXd TC2,
	Eigen::MatrixXd & OV, Eigen::MatrixXi & OT, Eigen::MatrixXd &OTC)
{
	using namespace Eigen;
	if (V1.rows() == 0) {
		OV = V2;
		OT = T2;
		OTC = TC2;
	}
	else if (V2.rows() == 0) {
		OV = V1;
		OT = T1;
		OTC = TC1;
	}

	MatrixXd bigV(V1.rows() + V2.rows(), 3);
	MatrixXi bigT(T1.rows() + T2.rows(), 3);
	MatrixXd bigTC(TC1.rows() + TC2.rows(), 3);
	bigV.block(0, 0, V1.rows(), 3) = V1;
	bigT.block(0, 0, T1.rows(), 3) = T1;
	bigTC.block(0, 0, TC1.rows(), 3) = TC1;
	bigV.block(V1.rows(), 0, V2.rows(), 3) = V2;
	bigT.block(T1.rows(), 0, T2.rows(), 3) = T2 + Eigen::MatrixXi::Constant(T2.rows(), T2.cols(), V1.rows());
	bigTC.block(TC1.rows(), 0, TC2.rows(), 3) = TC2;
	OV = bigV;
	OT = bigT;
	OTC = bigTC;

}
