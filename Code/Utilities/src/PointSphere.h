#ifndef POINTSPHERE_H
#define POINTSPHERE_H

#include <Eigen/Core>

bool PointSpheres(const Eigen::MatrixXd& points,
	const double& radius,
	const Eigen::MatrixXd& colors,
	const int res,
	const bool colorPerVertex,
	Eigen::MatrixXd& V,
	Eigen::MatrixXi& T,
	Eigen::MatrixXd& TC);


#endif // !POINTSPHERE
