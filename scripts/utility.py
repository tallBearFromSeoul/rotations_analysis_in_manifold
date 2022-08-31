import numpy as np

def deg2rad(deg: float):
    return deg*np.pi/180.0
def rad2deg(rad: float):
    return rad*180.0/np.pi

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


