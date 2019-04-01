#ifndef LINE_CYLINDER
#define LINE_CYLINDER

#include <Eigen/Core>

// This function is copied from libhedra library:
//		https://github.com/avaxman/libhedra/blob/master/include/hedra/line_cylinders.h	
bool LineCylinders(const Eigen::MatrixXd& P1,
	const Eigen::MatrixXd& P2,
	const double& radius,
	const Eigen::MatrixXd& C,
	const int res,
	const bool colorPerVertex,
	Eigen::MatrixXd& V,
	Eigen::MatrixXi& T,
	Eigen::MatrixXd& TC);

#endif // !LINE_CYLINDER
