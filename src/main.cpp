//#include "exp_and_log_map.hpp"
#include "cpp_utility.hpp"
#include "slerp.hpp"

int main(int argc, char *argv[]) {
	Vec3f k1(1,2,3), k2, k3;
	Mat3f R;
	k1 = k1.normalized();
	exp_map_from_so3_to_SO3(k1, R);
	log_map_from_SO3_to_so3(R, k2);
	std::cout<<"k1 :\n"<<k1<<"\n";
	std::cout<<"R :\n"<<R<<"\n";
	std::cout<<"k2 :\n"<<k2<<"\n";
	Quat q;
	exp_map_from_so3_to_SU2(k1, q);
	log_map_from_SU2_to_so3(q, k3);
	std::cout<<"k1 :\n"<<k1<<"\n";
	std::cout<<"q :\n"<<q<<"\n";
	std::cout<<"k3 :\n"<<k3<<"\n";
	std::cout<<"k1-k3 :\n"<<k1-k3<<"\n";

	Vec3f euler_angles_deg1 {125.f, 100.f, 40.f};
	Vec3f euler_angles_deg2 {75.f, 150.f, 250.f};
	Vec3f euler_angles_rad1 = deg2rad(euler_angles_deg1);
	Vec3f euler_angles_rad2 = deg2rad(euler_angles_deg2);
	std::cout<<"euler_angles_deg1 : "<<euler_angles_deg1<<"\n";
	std::cout<<"euler_angles_deg2 : "<<euler_angles_deg2<<"\n";
	std::cout<<"euler_angles_rad1 : "<<euler_angles_rad1<<"\n";
	std::cout<<"euler_angles_rad2 : "<<euler_angles_rad2<<"\n";
	Mat3f R1, R2;
	R_from_euler_angles_rad(euler_angles_rad1, R1);
	R_from_euler_angles_rad(euler_angles_rad2, R2);
	int n_steps = 30;
	std::vector<Mat3f> Rs_interp;
	slerp_R(R1, R2, n_steps, Rs_interp);
	std::cout<<"R1 :\n"<<R1<<"\n";
	std::cout<<"R2 :\n"<<R2<<"\n";
	print_vector(Rs_interp);
	return 0;
}


