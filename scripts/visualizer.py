from slerp import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def draw_sphere(n_pts: int = 50):
    xs, ys, zs = [], [], []
    for i in range(n_pts):
        for j in range(n_pts):
            xs.append(np.sin(np.pi*i/n_pts)*np.cos(2*np.pi*j/n_pts))
            ys.append(np.sin(np.pi*i/n_pts)*np.sin(2*np.pi*j/n_pts))
            zs.append(np.cos(np.pi*i/n_pts))
    return xs, ys, zs

def visualize_qs(qs_interp: list, draw_3d: bool):
    f, ax = plt.subplots(1,4, figsize=(20,6))
    f.suptitle(f'Spherical Linear Interpolation from q1 to q2 into {len(qs_interp)} points in Axis-Angle Representation')
    idxs = [*range(len(qs_interp))]
    #ks_norm1, ks_norm2, ks_x1, ks_y1, ks_z1, ks_x2,ks_y2, ks_z2 = [], [], [], [], [], [], [], [] 
    thetas, ks_x, ks_y, ks_z = [], [], [], []
    ks_norm, ks_norm_x, ks_norm_y, ks_norm_z = [], [], [], []
    #thetas_1, thetas_2 = [], []
    #qs_r1, qs_x1, qs_y1, qs_z1 = [], [], [], []
    #qs_r2, qs_x2, qs_y2, qs_z2 = [], [], [], []
    #for q1, q2 in zip(qs_interp, qs_interp_alter):
    for q_interp in qs_interp:
        k = log_map_from_SU2_to_so3(q_interp)
        k_norm = np.linalg.norm(k)
        thetas.append(rad2deg(k_norm))
        ks_x.append(k[0])
        ks_y.append(k[1])
        ks_z.append(k[2])
        ks_norm.append(k_norm)
        ks_norm_x.append(k[0]/k_norm)
        ks_norm_y.append(k[1]/k_norm)
        ks_norm_z.append(k[2]/k_norm)
        '''
        k1 = log_map_from_SU2_to_so3(q1)
        k2 = log_map_from_SU2_to_so3(q2)
        k1_norm = np.linalg.norm(k1)
        k2_norm = np.linalg.norm(k2)
        ks_norm1.append(k1_norm)
        ks_norm2.append(k2_norm)
        ks_x1.append(k1[0])
        ks_y1.append(k1[1])
        ks_z1.append(k1[2])
        ks_x2.append(k2[0])
        ks_y2.append(k2[1])
        ks_z2.append(k2[2])
        thetas_1.append(rad2deg(np.linalg.norm(k1)))
        thetas_2.append(rad2deg(np.linalg.norm(k2)))
        qs_r1.append(q1[0])
        qs_x1.append(q1[1])
        qs_y1.append(q1[2])
        qs_z1.append(q1[3])
        qs_r2.append(q2[0])
        qs_x2.append(q2[1])
        qs_y2.append(q2[2])
        qs_z2.append(q2[3])
        '''
    ax[0].scatter(idxs, thetas, alpha=0.5, s=8, marker='o', c='r')
    #ax[0].scatter(idxs, thetas_2, alpha=0.5, s=8, marker='+', c='r')
    #ax[0].legend(('q1', 'q0'))
    ax[0].set_title('axis angle in degrees of quaternion')
   
    ax[1].scatter(idxs, ks_norm, alpha=0.5, s=8, marker='o', c='orange')
    ax[1].scatter(idxs, ks_x, alpha=0.5, s=8, marker='o', c='r')
    ax[1].scatter(idxs, ks_y, alpha=0.5, s=8, marker='o', c='g')
    ax[1].scatter(idxs, ks_z, alpha=0.5, s=8, marker='o', c='b')
    '''
    ax[1].scatter(idxs, ks_norm2, alpha=0.5, s=8, marker='+', c='orange')
    ax[1].scatter(idxs, ks_x2, alpha=0.5, s=8, marker='+', c='r')
    ax[1].scatter(idxs, ks_y2, alpha=0.5, s=8, marker='+', c='g')
    ax[1].scatter(idxs, ks_z2, alpha=0.5, s=8, marker='+', c='b')
    '''
    ax[1].legend(('Frobenius norm of axis vector of quaternion\n equals to angle in radians', 'component in x axis', 'component in y axis', 'component in z axis'), fontsize=10)
    ax[1].set_title('axis vector of quaternion') 
    '''
    ax[2].scatter(idxs, qs_r1, alpha=0.5, s=8, marker='o', c='r')
    ax[2].scatter(idxs, qs_r2, alpha=0.5, s=8, marker='+', c='b')
    ax[2].scatter(idxs, qs_x1, alpha=0.5, s=8, marker='o', c='r')
    ax[2].scatter(idxs, qs_y1, alpha=0.5, s=8, marker='o', c='g')
    ax[2].scatter(idxs, qs_z1, alpha=0.5, s=8, marker='o', c='b')
    ax[2].scatter(idxs, qs_x2, alpha=0.5, s=8, marker='+', c='r')
    ax[2].scatter(idxs, qs_y2, alpha=0.5, s=8, marker='+', c='g')
    ax[2].scatter(idxs, qs_z2, alpha=0.5, s=8, marker='+', c='b')
    ax[2].legend(('q1_r','q2_r','q1_x','q1_y','q1_z','q2_x','q2_y','q2_z'))
    ax[2].set_title('quaternion components')
    '''
def visualize_Rs(Rs_interp: list, draw_3d: bool):
    f, ax = plt.subplots(1,4,figsize=(20,6))
    plt.subplots_adjust(left=0.03,bottom=0.05,right=0.97,top=0.88,wspace=0.5,hspace=0.4)
    f.suptitle(f'Spherical Linear Interpolation from R1 to R2 into {len(Rs_interp)} points in Axis-Angle Representation')

    idxs, thetas, ks_x, ks_y, ks_z = [*range(len(Rs_interp))], [], [], [], []
    ks_norm, ks_norm_x, ks_norm_y, ks_norm_z = [], [], [], []
    idxs_metrics = [*range(len(Rs_interp)-1)]
    geodesic_metrics, hyperbolic_metrics, frobenius_metrics = compute_metrics(Rs_interp[0], Rs_interp)
    for R_interp in Rs_interp:
        k_inv, theta_inv, R_log = log_map_from_SO3_to_so3(R_interp)
        thetas.append(rad2deg(theta_inv))
        ks_x.append(k_inv[0])
        ks_y.append(k_inv[1])
        ks_z.append(k_inv[2])
        k_norm = np.linalg.norm(k_inv)
        ks_norm.append(k_norm)
        ks_norm_x.append(k_inv[0]/k_norm)
        ks_norm_y.append(k_inv[1]/k_norm)
        ks_norm_z.append(k_inv[2]/k_norm)
    ax[0].scatter(idxs, thetas)
    ax[0].set_title('axis angle in degrees')
    ax[1].plot(idxs, ks_norm, c='orange', alpha=0.5)
    ax[1].plot(idxs, ks_x, c='r', alpha=0.5)
    ax[1].plot(idxs, ks_y, c='g', alpha=0.5)
    ax[1].plot(idxs, ks_z, c='b', alpha=0.5)
    ax[1].scatter(idxs, ks_norm, c='orange', alpha=0.5)
    ax[1].scatter(idxs, ks_x, c='r', s=8, alpha=0.5)
    ax[1].scatter(idxs, ks_y, c='g', s=8, alpha=0.5)
    ax[1].scatter(idxs, ks_z, c='b', s=8, alpha=0.5)
    ax[1].legend(('Frobenius norm of axis vector,\nequals to angle in radians','component in x axis','component in y axis','component in z axis'), fontsize=10)
    ax[1].set_title('axis vector')
    ax[2].plot(idxs_metrics, geodesic_metrics, c='r')
    ax[2].plot(idxs_metrics, hyperbolic_metrics, c='g')
    ax[2].plot(idxs_metrics, frobenius_metrics, c='b')
    ax[2].scatter(idxs_metrics, geodesic_metrics, s=8,  c='r')
    ax[2].scatter(idxs_metrics, hyperbolic_metrics, s=8,  c='g')
    ax[2].scatter(idxs_metrics, frobenius_metrics, s=8,  c='b')
    ax[2].legend(('geodesic metric', 'hyperbolic metric', 'frobenius metric'))
    ax[2].set_title('various metrics used\nbetween R1 and Ri')
    ax[3].plot(idxs_metrics, geodesic_metrics, c='r')
    ax[3].scatter(idxs_metrics, geodesic_metrics, s=8,  c='r')
    ax[3].set_yscale('log')
    ax[3].set_title('geodesic metric plotted \nin log scale for accuracy')
    plt.savefig('slerp_0.png')

    if draw_3d:
        sphere_x, sphere_y, sphere_z = draw_sphere(100)
        f3d = plt.figure(figsize=(10,10))
        f3d.suptitle('Axis unit vectors of SLERP from R1 to R2 shown alongside the unit sphere')
        ax3d = plt.axes(projection='3d')
        plt.subplots_adjust()
        ax3d.scatter3D(ks_x, ks_y, ks_z, marker='^', s=10, c='r', alpha=1.0)
        ax3d.scatter3D(ks_norm_x,ks_norm_y,ks_norm_z, marker='+', s=10, c='g', alpha=1.0)
        ax3d.scatter3D(sphere_x, sphere_y, sphere_z, marker='.', s=3, c='c', alpha=0.3)
        ax3d.legend(('axis-angle vector', 'normalized axis-angle vector', 'unit sphere'))
        for (k_x, k_y, k_z, k_norm) in zip(ks_x, ks_y, ks_z, ks_norm):
            ax3d.plot3D([0,k_x],[0,k_y],[0,k_z], c='r', linewidth=1, alpha=0.5)
            ax3d.plot3D([0,k_x/k_norm],[0,k_y/k_norm],[0,k_z/k_norm], c='g', linewidth=1, alpha=0.5)
        plt.savefig('slerp_1.png')

def main():
    #euler_angles_deg1 = np.array([125,100,40])
    #euler_angles_deg2 = np.array([75,150,150])

    euler_angles_deg1 = np.array([125,100,40])
    euler_angles_deg2 = np.array([75,150,250])
    euler_angles_rad1 = deg2rad(euler_angles_deg1)
    euler_angles_rad2 = deg2rad(euler_angles_deg2)
    print(f'euler_angles_rad1 :\n{euler_angles_rad1}')
    print(f'euler_angles_rad2 :\n{euler_angles_rad2}')

    R1 = R_from_euler_angles_rad(euler_angles_rad1)
    R2 = R_from_euler_angles_rad(euler_angles_rad2)
    print(f'R1 :\n{R1}')
    print(f'R2 :\n{R2}')
    Rs_interp = slerp_R(R1, R2, 30)
    q1 = q_from_euler_angles_rad(euler_angles_rad1)
    q2 = q_from_euler_angles_rad(euler_angles_rad2)
    qs_interp_direct = slerp_q_direct(q1, q2, 30)
    qs_interp_exp_and_log = slerp_q_exp_and_log(q1, q2, 30) 

    draw_3d = True
    visualize_Rs(Rs_interp, draw_3d)
    visualize_qs(qs_interp_direct, draw_3d)
    plt.show()

if __name__=='__main__':
    main()

