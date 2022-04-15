#include "LineCylinder.h"
#include <Eigen/Geometry>

bool LineCylinders(const Eigen::MatrixXd & P1, const Eigen::MatrixXd & P2, const double & radius, const Eigen::MatrixXd & C, const int res, const bool colorPerVertex, Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & TC)
{
	using namespace Eigen;
	V.resize(2 * res*P1.rows(), 3);
	T.resize(2 * res*P1.rows(), 3);
	RowVector3d ZAxis; ZAxis << 0.0, 0.0, 1.0;
	RowVector3d YAxis; YAxis << 0.0, 1.0, 0.0;

	int NewColorSize = (colorPerVertex ? V.rows() : T.rows());
	TC.resize(NewColorSize, 3);

	MatrixXd PlanePattern(res, 2);
	for (int i = 0; i < res; i++) {
		std::complex<double> CurrRoot = exp(2 * M_PI*std::complex<double>(0, 1)*(double)i / (double)res);
		PlanePattern.row(i) << CurrRoot.real(), CurrRoot.imag();
	}


	for (int i = 0; i < P1.rows(); i++) {
		RowVector3d NormAxis = (P2.row(i) - P1.row(i)).normalized();
		RowVector3d PlaneAxis1 = NormAxis.cross(ZAxis);
		if (PlaneAxis1.norm() < 10e-2)
			PlaneAxis1 = NormAxis.cross(YAxis).normalized();
		else
			PlaneAxis1 = PlaneAxis1.normalized();
		RowVector3d PlaneAxis2 = NormAxis.cross(PlaneAxis1).normalized();
		for (int j = 0; j < res; j++) {
			int v1 = 2 * res*i + 2 * j;
			int v2 = 2 * res*i + 2 * j + 1;
			int v3 = 2 * res*i + 2 * ((j + 1) % res);
			int v4 = 2 * res*i + 2 * ((j + 1) % res) + 1;
			V.row(v1) << P1.row(i) + (PlaneAxis1*PlanePattern(j, 0) + PlaneAxis2*PlanePattern(j, 1))*radius;
			V.row(v2) << P2.row(i) + (PlaneAxis1*PlanePattern(j, 0) + PlaneAxis2*PlanePattern(j, 1))*radius;

			if (colorPerVertex) {
				TC.row(v1) << C.row(i);
				TC.row(v2) << C.row(i);
			}


			T.row(2 * res*i + 2 * j) << v3, v2, v1;
			T.row(2 * res*i + 2 * j + 1) << v4, v2, v3;

			if (!colorPerVertex) {
				TC.row(2 * res*i + 2 * j) << C.row(i);
				TC.row(2 * res*i + 2 * j + 1) << C.row(i);
			}
		}
	}
	return true;
}
