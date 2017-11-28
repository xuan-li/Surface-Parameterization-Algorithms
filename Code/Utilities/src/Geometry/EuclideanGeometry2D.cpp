#include "EuclideanGeometry2D.h"

std::function<Complex(Complex)> ComputeRigidTransformation(Complex const s0, Complex const s1, Complex const t0, Complex const t1)
{
	Complex source_direction = s1 - s0;
	Complex target_direction = t1 - t0;
	Complex rotate_factor = target_direction / source_direction;
	auto result = [=](Complex p0)->Complex {return (p0 - s0)*rotate_factor + t0; };
	assert(std::abs(result(s1) - t1) < 1e-7);
	return [=](Complex p0)->Complex {return (p0 - s0)*rotate_factor + t0; };
}

Eigen::Matrix3d ComputeHomogeousRigidTransformation(Complex const s0, Complex const s1, Complex const t0, Complex const t1)
{
	Complex rotation = (t1 - t0) / (s1 - s0);
	Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d Ts = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d Tt = Eigen::Matrix3d::Identity();
	R(0, 0) = rotation.real();
	R(0, 1) = -rotation.imag();
	R(1, 0) = rotation.imag();
	R(1, 1) = rotation.real();
	Ts(0, 2) = -s0.real();
	Ts(1, 2) = -s0.imag();
	Ts(2, 2) = 1;
	Tt(0, 2) = t0.real();
	Tt(1, 2) = t0.imag();
	Tt(2, 2) = 1;
	Eigen::Matrix3d T = Tt * R * Ts;

	// test
	Eigen::Vector3d s(s1.real(), s1.imag(), 1);
	Eigen::Vector3d t(t1.real(), t1.imag(), 1);
	assert((T*s - t).norm() < 1e-5);
	return T;
}
