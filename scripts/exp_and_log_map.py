import numpy as np
import matplotlib.pyplot as plt

#import ctypes
#def ref_count(variable):
#    return ctypes.c_long.from_address(id(variable)).value

class Error(Exception):
    pass
class NonSquareMatrix(Error):
    pass
class NonSymmetricMatrix(Error):
    pass

class NotNumpyArray(Error):
    pass
class NotOneDimArray(Error):
    pass
class NotTwoDimArray(Error):
    pass

class Not3Vector(Error):
    pass
class Not3by3Matrix(Error):
    pass

class NotPositive(Error):
    pass
class NotInt(Error):
    pass
class NotFloat(Error):
    pass


def cross_product(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if not isinstance(a, np.ndarray) or not isinstance(b, np.ndarray):
        raise NotNumpyArray
    if a.ndim != 1 or b.ndim != 1:
        raise NotOneDimArray
    if a.shape[0] != 3 or b.shape[0] != 3:
        raise Not3Vector
    return np.array([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])

def trace(A: np.ndarray) -> float:
    if not isinstance(A, np.ndarray):
        raise NotNumpyArray
    if A.shape[0] != A.shape[1]:
        raise NonSquareMatrix
    trace = 0
    for i in range(A.shape[0]):
        trace += A[i][i]
    return trace

def vee_operator(A: np.ndarray) -> np.ndarray:
    if not isinstance(A, np.ndarray):
        raise NotNumpyArray
    if A.shape[0] != A.shape[1]:
        raise NonSquareMatrix
    if (-A != A.T).all():
        raise NonSymmetricMatrix
    return np.array([-A[1,2], A[0,2], -A[0,1]])


def skew_symmetric(v: np.ndarray) -> np.ndarray:
    if not isinstance(v, np.ndarray):
        raise NotNumpyArray
    if v.ndim != 1:
        raise NotOneDimArray
    return np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])

def skew_symmetric_inplace(v: np.ndarray, out: np.ndarray):
    if not isinstance(v, np.ndarray):
        raise NotNumpyArray
    if v.ndim != 1:
        raise NotOneDimArray
    if not isinstance(out, np.ndarray):
        raise NotNumpyArray
    if out.shape != (3,3):
        raise Not3by3Matrix
    out[0,1] = -v[2]
    out[0,2] = v[1]
    out[1,0] = v[2]
    out[1,2] = -v[0]
    out[2,0] = -v[1]
    out[2,1] = v[0]

# there are four vectors present : 
# w : the vector to rotate by angle theta about the axis of rotation k
# w_des : the desired final vector after the rotation
# k : axis of rotation
# theta : angle of rotation
def R_from_so3_to_SO3(w: np.ndarray, w_des: np.ndarray) -> np.ndarray:
    if not isinstance(w, np.ndarray) and not isinstance(w_des, np.ndarray):
        raise NotNumpyArray
    if w.ndim != 1 and not w_des.ndim != 1:
        raise NotOneDimArray
    cp = cross_product(w, w_des)
    # axis of rotation vector k
    k = cp / np.linalg.norm(cp)
    # angle of rotation theta (rad)
    theta = np.arccos(w.dot(w_des) / (np.linalg.norm(w) * np.linalg.norm(w_des)))
    K = skew_symmetric(k) 
    # R = I_3 + sin(theta)*K + (1-cos(theta))*K^2
    # where the skew_symmetric and symmetric decomposition of R is :
    # R = (R - R.T)/2 + (R + R.T)/2
    # any square matrix can be written uniquely as a skew-symmetric matrix
    # sin(theta)*K  : (R-R.T)/2 -> skew_symmetric
    # I_3 + (1-cos(theta))*K^2 : (R+R.T)/2 -> symmetric
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * K@K
    return k, theta, R

def exp_map_from_so3_to_SO3(k: np.ndarray, theta: float):
    if not isinstance(k, np.ndarray):
        raise NotNumpyArray
    if not isinstance(theta, float):
        raise NotFloat
    if k.ndim != 1:
        raise NotOneDimArray
    if k.shape[0] != 3:
        raise Not3Vector
    if theta < 0 or theta > 2*np.pi:
        theta = theta % (2*np.pi)
    K = skew_symmetric(k)
    # R = I_3 + sin(theta)*K + (1-cos(theta))*K^2
    # where the skew_symmetric and symmetric decomposition of R is :
    # R = (R - R.T)/2 + (R + R.T)/2
    # any square matrix can be written uniquely as a skew-symmetric matrix
    # sin(theta)*K  : (R-R.T)/2 -> skew_symmetric
    # I_3 + (1-cos(theta))*K^2 : (R+R.T)/2 -> symmetric
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * K@K
    return k, theta, R

def log_map_from_SO3_to_so3(R: np.ndarray) -> np.ndarray:
    if not isinstance(R, np.ndarray):
        raise NotNumpyArray
    if R.ndim != 2:
        raise NotTwoDimArray
    if R.shape[0] != R.shape[1]:
        raise NotSquareMatrix
    if R.shape[0] != 3:
        raise Not3by3Matrix
    tr = trace(R)
    # angle of rotation theta (rad)
    theta = np.arccos((trace(R) - 1) / 2)
    #print('theta :\n',theta)
    R_log = (1 / (2*np.sin(theta))) * (R - R.T)
    # axis of rotation vector k
    k = vee_operator(R_log)
    return k, theta, R_log

def main():
    w_s = [np.ones(3)*np.sqrt(1/3), np.array([np.sqrt(1/2), 0, np.sqrt(1/2)]), np.array([0, np.sqrt(1/2), np.sqrt(1/2)]), np.array([0,0,1])]
    w_des_s = [np.array([1,0,0]), np.array([1,0,0]), np.array([1,0,0]), np.array([1,0,0])]

    all_test_pass = True
    print_res = False

    for i, (w, w_des) in enumerate(zip(w_s, w_des_s)):
        test_pass = True
        print(f'\nTest {i} :\n')
        print('original vector w :\n',w)
        print('desired vector w_des :\n',w_des)
        
        k, theta, R = R_from_so3_to_SO3(w, w_des)
        if not np.allclose(R@w, w_des):
            test_pass = False
            all_test_pass = False
        if print_res:
            print('rotation matrix R obtained from Rodrigues formula,\n or exponential map of the original vector w to SO3 is :\n', R)
            if test_pass:
                print('---PASS---')
                print('w_new = R @ w == w_des is within tolerance of rtol=1e-05, atol=1e-08,\nwhere\n')
            else:
                print('---WRONG---')
                print('w_new = R @ w == w_des is NOT within tolerance of rtol=1e-05, atol=1e-08,\nwhere\n')
            print('w_new = R @ w is : ',R@w)
            print('w_des is : ',w_des)

        k_log, theta_log, R_log = log_map_from_SO3_to_so3(R)
        test_pass = True
        if not np.allclose(k, k_log):
            test_pass = False
            all_test_pass = False
        if print_res:
            print('logarithm map of the rotation matrix R to so3 :\n',R_log)
            if test_pass:
                print('---PASS---')
                print('the axis of rotation from the exponential map and logarithm map are within tolerance of rtol=1e-05, atol=1e-08,\nwhere\n')
            else:
                print('---WRONG---')
                print('the axis of rotation from the exponential map and logarithm map are NOT within tolerance of rtol=1e-05, atol=1e-08,\nwhere\n')
            print('axis of rotation from exponential map is :\n',k)
            print('axis of rotation from logarithm map is :\n',k_log)

        test_pass = True
        if not np.allclose(theta, theta_log):
            test_pass = False
            all_test_pass = False
        if print_res:
            if test_pass:
                print('---PASS---')
                print('the angle of rotation from the exponential map and logarithm map are within tolerance of rol=1e-05, atol=1e-08,\nwhere\n')
            else:
                print('---WRONG---')
                print('the angle of rotation from the exponential map and logarithm map are NOT within tolerance of rol=1e-05, atol=1e-08,\nwhere\n')
            print('angle of rotation is :\n',theta)
            print('angle of rotation is :\n',theta_log)


    if all_test_pass:
        print('---ALL_PASS---\n')
    '''  
    q = q_from_so3_to_SU2(w, w_des)
    print('quaternion q obtained from Rodrigues formula,\n or exponential map of the original vector w to SU2 is :\n',q)
    '''

    '''
    try:
        skew_symmetric(v, v_skew)
        print(id(v_skew))
    except NotNumpyArray:
        print('use numpy array')
    except NotOneDimArray:
        print('use one dim vector')
    print(v_skew)
    '''

if __name__=='__main__':
    main()