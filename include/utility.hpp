#ifndef UTILITY_HPP
#define UTILITY_HPP
#include <Eigen/Dense>

typedef Eigen::Vector3f Vec3f;
typedef Eigen::Matrix3f Mat3f;

const float _PI = 3.14159265359f;
const float _PI_OVER_180 = 0.01745329251f;
const float _180_OVER_PI = 57.2957795131f;

inline Vec3f deg2rad(const Vec3f &__deg) {
	return {__deg[0]*_PI_OVER_180, __deg[1]*_PI_OVER_180, __deg[2]*_PI_OVER_180};
}

inline Vec3f rad2deg(const Vec3f &__rad) {
	return {__rad[0]*_180_OVER_PI, __rad[1]*_180_OVER_PI, __rad[2]*_180_OVER_PI};
}


inline void skew_symmetric(const Vec3f &__k, Mat3f &K__) {
	K__ << 0, -__k(2), __k(1),
				 __k(2), 0, -__k(0),
				-__k(1), __k(0), 0;
}

void vee_operator(const Mat3f &__K, Vec3f &k__) {
	if (-__K != __K.transpose())
		throw std::invalid_argument("NotSymmetricMatrix");
	k__ << -__K(1,2), __K(0,2), -__K(0,1);
}

void R_from_euler_angles_rad(const Vec3f &__euler_angles_rad, Mat3f &R__) {
	float roll = __euler_angles_rad(0);
	float pitch = __euler_angles_rad(1);
	float yaw = __euler_angles_rad(2);
	float c_x = cos(roll);
	float c_y = cos(pitch);
	float c_z = cos(yaw);
	float s_x = sin(roll);
	float s_y = sin(pitch);
	float s_z = sin(yaw);
	Mat3f R_x {{1.f, 0.f, 0.f},
						 {0.f, c_x, -s_x},
						 {0.f, s_x, c_x}};
	Mat3f R_y {{c_y, 0, s_y},
						 {0.f, 1.f, 0.f},
						 {-s_y, 0.f, c_y}};
	Mat3f R_z {{c_z, -s_z, 0.f},
						 {s_z, c_z, 0.f},
						 {0.f, 0.f , 1.f}};
	R__ = R_z * R_y * R_x;
}


#endif
