#ifndef SLERP_HPP
#define SLERP_HPP
#include <Eigen/Dense>
#include <vector>
#include "exp_and_log_map.hpp"

// spherical linear interpolation (SLERP)
// R1 and R2 are 3x3 rotation matrices,
// where the function computes the geodesic line
// from R1 to R2 where lamda : [0, 1]
// R(0) = R1 and R(1) = R2
// the function returns the interpolated rotation matrices
// from lamda = 0 to lamda = 1
// based on the number of steps, default is 25
// the resulting motion is in constant angular velocity
// along the surface of the hypersphere
void slerp_R(const Mat3f &__R1, const Mat3f &__R2, int __n_steps, std::vector<Mat3f> &Rs_interp__) {
	if (Rs_interp__.size())
		Rs_interp__.clear();
	Vec3f k;
	log_map_from_SO3_to_so3(__R1.inverse()*__R2, k);
	for (int i=1; i<__n_steps; i++) {	
		Mat3f R;
		float lambda = 1.0f * static_cast<float>(i) / static_cast<float>(__n_steps-1);
		exp_map_from_so3_to_SO3(lambda*k, R);
		Rs_interp__.push_back(__R1*R);
	}
}

#endif
