from exp_and_log_map import *
from metrics import *

# spherical linear interpolation (SLERP)
# R1 and R2 are 3x3 rotation matrices,
# where the function computes the geodesic line
# from R1 to R2 where lamda : [0, 1]
# R(0) = R1 and R(1) = R2
# the function returns the interpolated rotation matrices
# from lamda = 0 to lamda = 1
# based on the number of steps, default is 25
# the resulting motion is in constant angular velocity
# along the surface of the hypersphere
def slerp(R1: np.ndarray, R2: np.ndarray, n_steps: int = 25):
    if not isinstance(R1, np.ndarray) and not isinstance(R2, np.ndarray):
        raise NotNumpyArray
    if R1.shape[0] != R1.shape[1] and R2.shape[0] != R2.shape[1]:
        raise NotSquareMatrix
    if R1.shape[0] != 3 and R2.shape[0] != 3:
        raise Not3by3Matrix
    if not isinstance(n_steps, int):
        raise NotInt
    if n_steps < 0:
        raise NotPositive
    Rs_interp = []
    for i, lamda in enumerate(np.linspace(0.0, 1.0, n_steps)):
        if i == 0:
            continue
        R1_inv = np.linalg.inv(R1)
        k_log, theta_log, R_log = log_map_from_SO3_to_so3(R1_inv @ R2)
        k, theta, R = exp_map_from_so3_to_SO3(lamda*k_log)
        Rs_interp.append(R1 @ R)
    return Rs_interp

# euler_angles in rad : [roll, pitch, yaw]
def R_from_euler_angles_rad(euler_angles_rad: np.ndarray):
    if not isinstance(euler_angles_rad, np.ndarray):
        raise NotNumpyArray
    if euler_angles_rad.ndim != 1:
        raise NotOneDimArray
    if euler_angles_rad.shape[0] != 3:
        raise Not3Vector
    roll = euler_angles_rad[0]
    pitch = euler_angles_rad[1]
    yaw = euler_angles_rad[2]
    R_x = np.array([[1,0,0],[0,np.cos(roll),-np.sin(roll)],[0,np.sin(roll),np.cos(roll)]])
    R_y = np.array([[np.cos(pitch),0,np.sin(pitch)],[0,1,0],[-np.sin(pitch),0,np.cos(pitch)]])
    R_z = np.array([[np.cos(yaw), -np.sin(yaw), 0],[np.sin(yaw),np.cos(yaw),0],[0,0,1]])
    return R_z @ R_y @ R_x

def main():
    euler_angles_deg1 = np.array([0,0,45])
    euler_angles_deg2 = np.array([45,45,45])
    R1 = R_from_euler_angles_rad(deg2rad(euler_angles_deg1))
    R2 = R_from_euler_angles_rad(deg2rad(euler_angles_deg2))
    k_inv1, theta_inv1, R_log1 = log_map_from_SO3_to_so3(R1)
    k_inv2, theta_inv2, R_log2 = log_map_from_SO3_to_so3(R2)
    print(f'k1 :\n{k_inv1}\nk2 :\n{k_inv2}')
    print(f'theta1 : {rad2deg(theta_inv1)}\t theta2 : {rad2deg(theta_inv2)}')
    print(f'R_log1 :\n{R_log1}\nR_log2 :\n{R_log2}')
    slerp(R1, R2, 3)

if __name__=='__main__':
    main()


