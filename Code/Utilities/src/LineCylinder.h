#ifndef LINE_CYLINDER
#define LINE_CYLINDER

#include <Eigen/Core>

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
