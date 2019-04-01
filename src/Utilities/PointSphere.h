#ifndef POINTSPHERE_H
#define POINTSPHERE_H

#include <Eigen/Core>

// This function is copied from libhedra library:
//		https://github.com/avaxman/libhedra/blob/master/include/hedra/point_spheres.h
bool PointSpheres(const Eigen::MatrixXd& points,
	const double& radius,
	const Eigen::MatrixXd& colors,
	const int res,
	const bool colorPerVertex,
	Eigen::MatrixXd& V,
	Eigen::MatrixXi& T,
	Eigen::MatrixXd& TC);


#endif // !POINTSPHERE
