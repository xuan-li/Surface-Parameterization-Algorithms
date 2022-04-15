#ifndef EULIDEAN_GEOMETRY_2D_H_
#define EULIDEAN_GEOMETRY_2D_H_


#include <complex>
#include <iostream>
#include <Eigen/Core>
#include <functional>

typedef std::complex<double> Complex;

std::function<Complex(Complex)> ComputeRigidTransformation(Complex const s0, Complex const s1, Complex const t0, Complex const t1);

Eigen::Matrix3d ComputeHomogeousRigidTransformation(Complex const s0, Complex const s1, Complex const t0, Complex const t1);


#endif