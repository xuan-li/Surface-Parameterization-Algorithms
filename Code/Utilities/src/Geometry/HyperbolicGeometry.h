#ifndef HYPERBOLIC_GEOMETRY_H_
#define HYPERBOLIC_GEOMETRY_H_
#include <functional>
#include <iostream>
#include <complex>
#include "Circle.h"

typedef std::complex<double> Complex;

/*
These codes are used to handle geometry on Poincare disk
*/



std::function<Complex(Complex const)> ComputeMobiusTransformation(Complex const s0, Complex const s1, Complex const t0, Complex const t1);

double HyperbolicDistance(Complex p0, Complex p1);

Complex InverseExponentialMap(Complex p0, Complex p1);

#endif // !HYPERBOLIC_GEOEMETRY_H_
