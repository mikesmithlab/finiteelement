import numpy as np

def get_reduced_force_vector(F, disp_ids):
    F_red = np.copy(F)
    F_red = np.delete(F_red, disp_ids, 0)
    return F_red