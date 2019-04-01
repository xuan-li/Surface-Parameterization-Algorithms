#ifndef CIRCLE_2D_H_
#define CIRCLE_2D_H_
#endif // !CIRCLE_2D_H_

#include <complex>
#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

class Circle {
public:
	Circle(std::complex<double> center, double radius);
	std::complex<double> center() { return center_; }
	double radius() { return radius_; }
	void set_center(std::complex<double> center) { center_ = center; }
	void set_radius(double radius) { radius_ = radius; }
	std::complex<double> Tangent(std::complex<double> point);
	double PolarAngle(std::complex<double> point);
protected:
	std::complex<double> center_;
	double radius_;
};

Circle CircleFromPoints(std::vector<std::complex<double>> points);

std::complex<double> ReflectByCircle(std::complex<double> point, Circle circle);