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

def visualize(Rs_interp: list, draw_3d: bool):
    f, ax = plt.subplots(1,4,figsize=(20,6))
    plt.subplots_adjust(left=0.03,bottom=0.05,right=0.97,top=0.88,wspace=0.5,hspace=0.4)
    f.suptitle(f'Spherical Linear Interpolation from R1 to R2 into {len(Rs_interp)} points in Axis-Angle Representation')

    idxs, thetas, ks_x, ks_y, ks_z = [*range(len(Rs_interp))], [], [], [], []
    ks_norm, ks_norm_x, ks_norm_y, ks_norm_z = [], [], [], []
    idxs_metrics = [*range(len(Rs_interp)-1)]
    print(f'Rs_interp[0] :\n{Rs_interp[0]}')
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
    plt.show()

def main():
    euler_angles_deg1 = np.array([125,100,40])
    euler_angles_deg2 = np.array([75,150,150])
    R1 = R_from_euler_angles_rad(deg2rad(euler_angles_deg1))
    R2 = R_from_euler_angles_rad(deg2rad(euler_angles_deg2))
    k_inv1, theta_inv1, R_log1 = log_map_from_SO3_to_so3(R1)
    k_inv2, theta_inv2, R_log2 = log_map_from_SO3_to_so3(R2)
    Rs_interp = slerp(R1, R2, 30)
    Rs_interp[0] = R1
    #print(f'Rs_interp :\n{Rs_interp}')
    #print(f'k1 :\n{k_inv1}\nk2 :\n{k_inv2}')
    #print(f'theta1 : {rad2deg(theta_inv1)}\t theta2 : {rad2deg(theta_inv2)}')
    #print(f'R_log1 :\n{R_log1}\nR_log2 :\n{R_log2}')

    draw_3d = True
    visualize(Rs_interp, draw_3d)

if __name__=='__main__':
    main()

