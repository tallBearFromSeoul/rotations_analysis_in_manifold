
def rodrigues_k_vec(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if not isinstance(a, np.ndarray) or not isinstance(b, np.ndarray):
        raise NotNumpyArray
    if a.ndim != 1 or b.ndim != 1:
        raise NotOneDimArray
    if a.shape[0] != 3 or b.shape[0] != 3:
        raise Not3Vector
    cross = cross_product(a, b)
    l1_norm = cross.sum()
    return cross / l1_norm

def rodrigues_rotation_mat(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    k_vec = rodrigues_k_vec(a, b)
    K = skew_symmetric(k_vec)
    return np.eye(3) + np.sin(theta)*K + (1-cos(theta))*K@K

def q_from_so3_to_SU2(w: np.ndarray, w_des: np.ndarray) -> np.ndarray:
    if not isinstance(w, np.ndarray) and not isinstance(w_des, np.ndarray):
        raise NotNumpyArray
    if w.ndim != 1 and not w_des.ndim != 1:
        raise NotOneDimArray
    cp = cross_product(w, w_des)
    theta = np.arccos(w.dot(w_des) / np.linalg.norm(w) * np.linalg.norm(w_des))
    if (np.allclose(w.sum(),0)):
        return np.array([1,0,0,0])
    else:
        return np.array([np.cos(theta/2), np.sin(theta/2)*w[0]/theta, np.sin(theta/2)*w[1]/theta, np.sin(theta/2)*w[2]/theta])



