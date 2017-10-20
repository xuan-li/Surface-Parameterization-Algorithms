#include "Circle.h"

Circle::Circle(std::complex<double> center, double radius):center_(center), radius_(radius)
{

}

std::complex<double> Circle::Tangent(std::complex<double> point)
{
	std::complex<double> diff = point - center_;
	diff /= abs(diff);
	return std::complex<double>(-diff.imag(), diff.real());
}

double Circle::PolarAngle(std::complex<double> point)
{
	std::complex<double> diff = point - center_;
	diff /= abs(diff);
	double t = diff.imag() / diff.real();
	double angle = atan(t);
	if (diff.imag() < 0) {
		angle += 3.141592653;
	}
	return angle;
}

Circle CircleFromPoints(std::vector<std::complex<double>> points)
{
	assert(points.size() >= 3);
	std::complex<double> p1 = points[0];
	std::complex<double> p2 = points[1];
	std::complex<double> p3 = points[2];

	double offset = p2.real() * p2.real() + p2.imag() * p2.imag();
	double bc = (p1.real() * p1.real() + p1.imag() * p1.imag() - offset) / 2.0;
	double cd = (offset - p3.real() * p3.real() - p3.imag() * p3.imag()) / 2.0;
	double det = (p1.real() - p2.real()) * (p2.imag() - p3.imag()) - (p2.real() - p3.real())* (p1.imag() - p2.imag());

	assert(abs(det) < 1e-7);

	double idet = 1 / det;

	double centerx = (bc * (p2.imag() - p3.imag()) - cd * (p1.imag() - p2.imag())) * idet;
	double centery = (cd * (p1.real() - p2.real()) - bc * (p2.real() - p3.real())) * idet;
	std::complex<double> center(centerx, centery);
	double radius = abs(p2 - center);

	return Circle(center, radius);
}

std::complex<double> ReflectByCircle(std::complex<double> point, Circle circle)
{
	double dist_original = abs(point - circle.center());
	double dist_target = circle.radius() * circle.radius() / dist_original;
	return circle.center() + (point - circle.center()) * dist_target / dist_original;
}
