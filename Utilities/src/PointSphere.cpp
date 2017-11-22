#include "PointSphere.h"

bool PointSpheres(const Eigen::MatrixXd & points, const double & radius, const Eigen::MatrixXd & colors, const int res, const bool colorPerVertex, Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & TC)
{
	using namespace Eigen;
	V.resize(res*res*points.rows(), 3);
	T.resize(2 * (res - 1)*res*points.rows(), 3);
	TC.resize((colorPerVertex ? V.rows() : T.rows()), 3);

	for (int i = 0; i < points.rows(); i++) {
		RowVector3d center = points.row(i);

		//creating vertices
		for (int j = 0; j < res; j++) {
			double z = center(2) + radius*cos(M_PI*(double)j / (double(res - 1)));
			for (int k = 0; k < res; k++) {
				double x = center(0) + radius*sin(M_PI*(double)j / (double(res - 1)))*cos(2 * M_PI*(double)k / (double(res - 1)));
				double y = center(1) + radius*sin(M_PI*(double)j / (double(res - 1)))*sin(2 * M_PI*(double)k / (double(res - 1)));
				V.row((res*res)*i + j*res + k) << x, y, z;
				if (colorPerVertex)
					TC.row((res*res)*i + j*res + k) << colors.row(i);
			}
		}


		//creating faces
		for (int j = 0; j < res - 1; j++) {
			for (int k = 0; k < res; k++) {
				int v1 = (res*res)*i + j*res + k;
				int v2 = (res*res)*i + (j + 1)*res + k;
				int v3 = (res*res)*i + (j + 1)*res + (k + 1) % res;
				int v4 = (res*res)*i + j*res + (k + 1) % res;
				T.row(2 * (((res - 1)*res)*i + res*j + k)) << v1, v2, v3;
				T.row(2 * (((res - 1)*res)*i + res*j + k) + 1) << v4, v1, v3;
				if (!colorPerVertex) {
					TC.row(2 * (((res - 1)*res)*i + res*j + k)) << colors.row(i);
					TC.row(2 * (((res - 1)*res)*i + res*j + k) + 1) << colors.row(i);
				}
			}
		}
	}

	return true;
}

