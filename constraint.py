import sys
import matplotlib.pyplot as plt
import numpy as np
from visualise import show_msh
from DataAnalysis.select_data import get_pts
from gmsh import read_gmsh

    
def set_constraint(node_ids, node_coords, element_nodes, DoF=None, value=None):
    """
    DoF - 'x', 'y'
    value - a numerical value of displacement or force.
    """
    if DoF == 'x':
        DoF = 0
    elif DoF == 'y':
        DoF = 1
    else:
        raise Exception

    
    fig, ax = plt.subplots()
    ax, node_pts = show_msh(ax, node_coords, element_nodes)

    constrained_node_ids, _ = get_pts(ax, node_pts, type='Rectangle',enter_closes=True)

    displacement_constraints = np.empty((msh_params['node']['total_num'],3))*np.nan   
    displacement_constraints[:,0] = node_ids
    if constrained_node_ids is not None:
        for node_id in constrained_node_ids:
            displacement_constraints[node_id][DoF+1]=value
    return displacement_constraints

def combine_constraints(constraint, total=None, add=False):
    """
    Add constraints to combined vector. 
    If position is already occupied the value can be added if add=True
    or ignored if add=False
    """
    if total is None:
        total = np.empty(np.shape(constraint))*np.nan
        total[:,0] = constraint[:,0]
    for index, val in enumerate(constraint):
        #x constraints
        if np.isnan(total[index, 1]):
            total[index, 1] = val[1]
        elif add:
            #site already has non nan value. add to existing number if add=True
            total[index, 1] = total[index, 1] + val[1]
        
        #y constrainsts
        if np.isnan(total[index, 2]):
            total[index, 2] = val[2]
        elif add:
            #site already has non nan value. add to existing number if add=True
            total[index, 2] = total[index, 2] + val[2]

    return total



if __name__ == '__main__':
    #mesh_dir = "C:\\Users\\ppzmis\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_dir = "C:\\Users\\mikei\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_filename = 'rect.msh'
    constraint_type = 'Displacement'
    #constraint_type = 'Force'

    node_ids, node_coords, element_ids, element_nodes, msh_params = read_gmsh(mesh_dir + mesh_filename)

    if constraint_type == 'Displacement':
        values = set_constraint(node_ids, node_coords, element_nodes, DoF='x',value=0)
        values2 = set_constraint(node_ids, node_coords, element_nodes, DoF='y',value=0)
        constraints = combine_constraints(values)
        constraints = combine_constraints(values2, total=constraints)
    else:
        values = set_constraint(node_ids, node_coords, element_nodes, DoF='y',value=-10000)
        #values2 = set_constraint(node_ids, node_coords, element_nodes, DoF='y',value=-10000)
        constraints = combine_constraints(values)
        #constraints = combine_constraints(values2, total=constraints)
      

    np.savetxt(mesh_dir + mesh_filename[:-4] + '_' + constraint_type + '.constraint', constraints)

     

    
    
    


    #Support X


    #Support Y

   
   
