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
