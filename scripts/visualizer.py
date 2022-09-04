from slerp import *
from conversion import *
import matplotlib as mpl
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

def visualize(init_conditions: tuple, Rs_interp: list, qs_interp: list, draw_3d: bool, idx: int):
    ### MPL Formatting
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'
    mpl.rcParams['legend.fontsize'] = 12
    #mpl.rcParams['lines.linewidth'] = 1
    #mpl.rcParams['lines.linestyle'] = 'dashdot'
    
    euler_angles_deg1 = init_conditions[0]
    euler_angles_deg2 = init_conditions[1]
    R1 = init_conditions[2]
    R2 = init_conditions[3]
    q1 = init_conditions[4]
    q2 = init_conditions[5]
    f, ax = plt.subplots(1,4,figsize=(20,6))
    plt.subplots_adjust(left=0.03,bottom=0.15,right=0.97,top=0.7,wspace=0.5,hspace=0.4)
    suptitle = f'Comparison of Spherical Linear Interpolation from A to B between using rotation matrices and quaternions to {len(Rs_interp)} points in euler angles and axis-angle representation, '+'where\n' 
    suptitle += r'$A : euler\ angles_{deg} = [%.1f, %.1f, %.1f]^T$'%(euler_angles_deg1[0], euler_angles_deg1[1], euler_angles_deg1[2])+'\n'
    suptitle += r'$B : euler\ angles_{deg} = [%.1f, %.1f, %.1f]^T$'%(euler_angles_deg2[0], euler_angles_deg2[1], euler_angles_deg2[2])+'\n'
    suptitle += r'$R_A = \begin{bmatrix} %.3f & %.3f & %.3f \\ %.3f & %.3f & %.3f \\ %.3f & %.3f & %.3f \end{bmatrix}$'%(R1[0,0],R1[0,1],R1[0,2],R1[1,0],R1[1,1],R1[1,2],R1[2,0],R1[2,1],R1[2,2])
    suptitle += r'$R_B = \begin{bmatrix} %.3f & %.3f & %.3f \\ %.3f & %.3f & %.3f \\ %.3f & %.3f & %.3f \end{bmatrix}$'%(R1[0,0],R1[0,1],R1[0,2],R1[1,0],R1[1,1],R1[1,2],R1[2,0],R1[2,1],R1[2,2])
    suptitle += r'$q_A = [%.3f, %.3f, %.3f, %.3f]^T$'%(q1[0], q1[1], q1[2], q1[3])
    suptitle += r'$q_B = [%.3f, %.3f, %.3f, %.3f]^T$'%(q2[0], q2[1], q2[2], q2[3])+'\n'
    f.suptitle(suptitle)

    ### Data Extraction
    idxs, idxs_metrics = [*range(len(Rs_interp))], [*range(len(Rs_interp)-1)]
    R_rolls, R_pitchs, R_yaws = [], [], []
    q_rolls, q_pitchs, q_yaws = [], [], []
    R_thetas, R_ks_x, R_ks_y, R_ks_z = [], [], [], []
    q_thetas, q_ks_x, q_ks_y, q_ks_z = [], [], [], []
    R_ks_norm, R_ks_norm_x, R_ks_norm_y, R_ks_norm_z = [], [], [], []
    q_ks_norm, q_ks_norm_x, q_ks_norm_y, q_ks_norm_z = [], [], [], []

    #geodesic_metrics, hyperbolic_metrics, frobenius_metrics = compute_metrics(Rs_interp[0], Rs_interp)
    for R_interp in Rs_interp:
        euler_angles_rad = euler_angles_rad_from_R(R_interp)
        R_rolls.append(rad2deg(euler_angles_rad[0]))
        R_pitchs.append(rad2deg(euler_angles_rad[1]))
        R_yaws.append(rad2deg(euler_angles_rad[2]))
        
        k = log_map_from_SO3_to_so3(R_interp)
        theta = np.linalg.norm(k)
        R_thetas.append(rad2deg(theta))
        R_ks_x.append(k[0])
        R_ks_y.append(k[1])
        R_ks_z.append(k[2])
        k_norm = np.linalg.norm(k)
        R_ks_norm.append(k_norm)
        R_ks_norm_x.append(k[0]/k_norm)
        R_ks_norm_y.append(k[1]/k_norm)
        R_ks_norm_z.append(k[2]/k_norm)
    for q_interp in qs_interp:
        euler_angles_rad = euler_angles_rad_from_q(q_interp)
        q_rolls.append(rad2deg(euler_angles_rad[0]))
        q_pitchs.append(rad2deg(euler_angles_rad[1]))
        q_yaws.append(rad2deg(euler_angles_rad[2]))

        k = log_map_from_SU2_to_so3(q_interp)
        theta = np.linalg.norm(k)
        q_thetas.append(rad2deg(theta))
        q_ks_x.append(k[0])
        q_ks_y.append(k[1])
        q_ks_z.append(k[2])
        q_ks_norm.append(k_norm)
        q_ks_norm_x.append(k[0]/k_norm)
        q_ks_norm_y.append(k[1]/k_norm)
        q_ks_norm_z.append(k[2]/k_norm)

    ax[0].scatter(idxs, R_rolls, c='r', marker='o', s=20, alpha=0.5)
    ax[0].scatter(idxs, R_pitchs, c='g', marker='o', s=20, alpha=0.5)
    ax[0].scatter(idxs, R_yaws, c='b', marker='o', s=20, alpha=0.5)
    ax[0].scatter(idxs, q_rolls, c='m', marker='^', s=20, alpha=0.5)
    ax[0].scatter(idxs, q_pitchs, c='lime', marker='^', s=20, alpha=0.5)
    ax[0].scatter(idxs, q_yaws, c='c', marker='^', s=20, alpha=0.5)
    ax[0].legend(('R_roll','R_pitch','R_yaw','q_roll','q_pitch','q_yaw'))
    ax[0].set_title('euler angles in degrees')

    ax[1].scatter(idxs, R_thetas, c='r', marker='o', s=20, alpha=0.5)
    ax[1].scatter(idxs, q_thetas, c='m', marker='^', s=20, alpha=0.5)
    ax[1].legend(('R_theta', 'q_theta'))
    ax[1].set_title('axis angle in degrees')
    
    ax[2].scatter(idxs, R_ks_norm, c='orange', marker='o', s=20, alpha=0.5)
    ax[2].scatter(idxs, R_ks_x, c='r', marker='o', s=20, alpha=0.5)
    ax[2].scatter(idxs, R_ks_y, c='g', marker='o', s=20, alpha=0.5)
    ax[2].scatter(idxs, R_ks_z, c='b', marker='o', s=20, alpha=0.5)
    ax[2].legend(('norm of axis vector','R_k_x','R_k_y','R_k_z'))
    ax[2].set_title('axis vector')
    plt.savefig(f'R_q_slerp_comparison_{idx}.png')

    '''
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
    '''
    if draw_3d:
        sphere_x, sphere_y, sphere_z = draw_sphere(100)
        f3d = plt.figure(figsize=(12,12))
        f3d.suptitle(suptitle)
        #f3d.suptitle('Axis unit vectors of SLERP from R1 to R2 shown alongside the unit sphere')
        ax3d = plt.axes(projection='3d')
        plt.subplots_adjust()
        ax3d.scatter3D(R_ks_x, R_ks_y, R_ks_z, marker='^', s=20, c='r', alpha=1.0)
        ax3d.scatter3D(R_ks_norm_x,R_ks_norm_y, R_ks_norm_z, marker='+', s=20, c='g', alpha=1.0)
        ax3d.scatter3D(q_ks_x, q_ks_y, q_ks_z, marker='^', s=20, c='m', alpha=1.0)
        ax3d.scatter3D(q_ks_norm_x,q_ks_norm_y, q_ks_norm_z, marker='+', s=20, c='lime', alpha=1.0)
        ax3d.scatter3D(sphere_x, sphere_y, sphere_z, marker='.', s=3, c='c', alpha=0.3)
        ax3d.legend(('R_ks', 'normalized R_ks', 'q_ks', 'normalized q_ks', 'unit sphere'))
        for (k_x, k_y, k_z, k_norm) in zip(R_ks_x, R_ks_y, R_ks_z, R_ks_norm):
            ax3d.plot3D([0,k_x],[0,k_y],[0,k_z], c='r', linewidth=1, alpha=0.5)
            ax3d.plot3D([0,k_x/k_norm],[0,k_y/k_norm],[0,k_z/k_norm], c='g', linewidth=1, alpha=0.5)
        for (k_x, k_y, k_z, k_norm) in zip(q_ks_x, q_ks_y, q_ks_z, q_ks_norm):
            ax3d.plot3D([0,k_x],[0,k_y],[0,k_z], c='m', linewidth=1, alpha=0.5)
            ax3d.plot3D([0,k_x/k_norm],[0,k_y/k_norm],[0,k_z/k_norm], c='lime', linewidth=1, alpha=0.5)
            
        plt.savefig(f'R_q_slerp_comparison_3d_{idx}.png')

def compute_and_visualize(euler_angles_deg1: np.ndarray, euler_angles_deg2: np.ndarray, draw_3d: bool, idx: int):
    if not isinstance(euler_angles_deg1, np.ndarray) or not isinstance(euler_angles_deg2, np.ndarray):
        raise NotNumpyArray
    if euler_angles_deg1.ndim != 1 or euler_angles_deg2.ndim != 1:
        raise NotOneDimArray
    if euler_angles_deg1.shape[0] != 3 or euler_angles_deg2.shape[0] != 3:
        raise Not3Vector
    euler_angles_rad1 = deg2rad(euler_angles_deg1)
    euler_angles_rad2 = deg2rad(euler_angles_deg2)

    R1 = R_from_euler_angles_rad(euler_angles_rad1)
    R2 = R_from_euler_angles_rad(euler_angles_rad2)
    Rs_interp = slerp_R(R1, R2, 30)
    
    q1 = q_from_euler_angles_rad(euler_angles_rad1)
    q2 = q_from_euler_angles_rad(euler_angles_rad2)
    qs_interp_direct = slerp_q_direct(q1, q2, 30)
    qs_interp_exp_and_log = slerp_q_exp_and_log(q1, q2, 30) 

    init_conditions = (euler_angles_deg1, euler_angles_deg2, R1, R2, q1, q2)
    visualize(init_conditions, Rs_interp, qs_interp_direct, draw_3d, idx)

def main():
    draw_3d = True
    euler_angles_deg1_s = [np.array([125,100,40]), np.array([125,100,40])]
    euler_angles_deg2_s = [np.array([75,150,150]), np.array([75,150,250])]
    for (i, (euler_angles_deg1, euler_angles_deg2)) in enumerate(zip(euler_angles_deg1_s, euler_angles_deg2_s)):
        compute_and_visualize(euler_angles_deg1, euler_angles_deg2, draw_3d, i)

    '''
    euler_angles_rad1 = deg2rad(euler_angles_deg1)
    euler_angles_rad2 = deg2rad(euler_angles_deg2)

    R1 = R_from_euler_angles_rad(euler_angles_rad1)
    R2 = R_from_euler_angles_rad(euler_angles_rad2)
    Rs_interp = slerp_R(R1, R2, 30)
    
    q1 = q_from_euler_angles_rad(euler_angles_rad1)
    q2 = q_from_euler_angles_rad(euler_angles_rad2)
    qs_interp_direct = slerp_q_direct(q1, q2, 30)
    qs_interp_exp_and_log = slerp_q_exp_and_log(q1, q2, 30) 

    init_conditions = (euler_angles_deg1, euler_angles_deg2, R1, R2, q1, q2)
    visualize(init_conditions, Rs_interp, qs_interp_direct, draw_3d)
    ''' 
    plt.show()

if __name__=='__main__':
    main()

