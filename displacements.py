import numpy as np

def restrained_DoF(displacements):
    """
    Find locations of node ids associated with restrained
    DoF. N.B returned vector has x and y ids merged so 
    node id 2 in y direction has value of 4.
    """
    dispx_ids = displacements[:,:2]
    dispx_ids = 2*dispx_ids[~np.isnan(dispx_ids).any(axis=1)][:,0]-2
    dispy_ids = np.delete(displacements, 1, axis=1)
    dispy_ids = 2*dispy_ids[~np.isnan(dispy_ids).any(axis=1)][:,0]-1
    disp_ids =  np.int_(np.concatenate((dispx_ids,dispy_ids))) 
    return disp_ids

def calc_global_disp(U, disp_ids, disp):
    """
    Construct global displacement vector
    """ 
    num_nodes = np.shape(disp)[0]
    UG=np.zeros((num_nodes*2,1))
    #Add specified displacements 
    mask = np.zeros((num_nodes*2), dtype=bool)
    mask[disp_ids] = True
    UG[mask] = disp[:,1:].flatten()[disp_ids].reshape((np.size(disp_ids),1))
    UG[~mask] = U
    return UG

def calc_new_node_coords(node_coords, UG, scale=1):
    node_coords[:,:2]=node_coords[:,:2]+scale*UG.reshape((-1,2))
    return node_coords