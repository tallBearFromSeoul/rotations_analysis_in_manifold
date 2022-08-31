import numpy as np
import matplotlib.pyplot as plt
from utility import *

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
class DivisionByZero(Error):
    pass

# R = I_3 + sin(theta)*K + (1-cos(theta))*K^2
    # where the skew_symmetric and symmetric decomposition of R is :
    # R = (R - R.T)/2 + (R + R.T)/2
    # any square matrix can be written uniquely as a skew-symmetric matrix
    # sin(theta)*K  : (R-R.T)/2 -> skew_symmetric
    # I_3 + (1-cos(theta))*K^2 : (R+R.T)/2 -> symmetric
def exp_map_from_so3_to_SO3(k: np.ndarray):
    if not isinstance(k, np.ndarray):
        raise NotNumpyArray
    if k.ndim != 1:
        raise NotOneDimArray
    if k.shape[0] != 3:
        raise Not3Vector
    K = skew_symmetric(k)
    theta = np.linalg.norm(k)
    if theta == 0:
        raise DivisionByZero
    R = np.eye(3) + np.sin(theta)/theta * K + (1 - np.cos(theta))/(theta**2) * K@K
    return k, theta, R

def log_map_from_SO3_to_so3(R: np.ndarray):
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
    R_log = (theta / (2*np.sin(theta))) * (R - R.T)
    # axis of rotation vector k
    k = vee_operator(R_log)
    return k, theta, R_log

