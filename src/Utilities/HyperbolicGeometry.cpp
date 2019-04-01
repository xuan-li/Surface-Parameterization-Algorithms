#include "HyperbolicGeometry.h"
#include <assert.h>

std::function<Complex(Complex const)> ComputeMobiusTransformation(Complex const s0, Complex const s1, Complex const t0, Complex const t1)
{
	auto s0_to_zero = [=](Complex const c)->Complex {
		return (c - s0) / (Complex(1,0) - std::conj(s0)*c);
	};

	auto t0_to_zero = [=](Complex const c)->Complex {
		return (c - t0) / (Complex(1, 0) - std::conj(t0)*c);
	};

	auto zero_to_t0 = [=](Complex const c)->Complex {
		return (c + t0) / (Complex(1, 0) + std::conj(t0)*c);
	};

	Complex new_s1 = s0_to_zero(s1);
	Complex new_t1 = t0_to_zero(t1);
	
	Complex rotation_coeff = new_t1 / new_s1;

	auto result = [=](Complex const c)->Complex {
		Complex resualt = (c - s0) / (Complex(1, 0) - std::conj(s0)*c);
		resualt *= rotation_coeff;
		resualt = (resualt + t0) / (Complex(1, 0) + std::conj(t0)*resualt);
		return resualt;
	};
	auto test = result(s1);
	assert(std::abs(result(s1) - t1) < 1e-6);
	return result;
}

double HyperbolicDistance(Complex p0, Complex p1)
{
	double dist  = acosh(
		1 + 2 * abs(p0 - p1) * abs(p0 - p1) /
		((1. - abs(p0) * abs(p0)) * (1 - abs(p1) * abs(p1)))
	);
	return dist;
}

Complex InverseExponentialMap(Complex p0, Complex p1)
{
	Circle unit(Complex(0, 0), 1.0);
	Complex reflection = ReflectByCircle(p0, unit);
	Circle geodesic_circle = CircleFromPoints(std::vector<Complex>({ p0, p1, reflection }));
	double tangent_length = pow(1 - abs(p0) * abs(p0),2) / 4;
	double angle0 = geodesic_circle.PolarAngle(p0);
	double angle1 = geodesic_circle.PolarAngle(p1);
	Complex tangent = geodesic_circle.Tangent(p0);
	if (angle1 < angle0)
		tangent = -tangent;
	double dist = HyperbolicDistance(p0, p1);
	double ratio = pow(1 - pow(abs(p0), 2), 2) / 4.0;
	return tangent * dist * ratio;
}
