
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from DataAnalysis.select_data import get_pts
from gmsh import read_gmsh









        
    
def show_msh(ax, node_coords, element_nodes):
    #Plot node pts
    node_pts = ax.scatter(node_coords[:,0],node_coords[:,1],s=5,c='b')

    #Plot gray patch for mesh
    for n, element in enumerate(element_nodes):
        x = np.array([node_coords[element[0]-1][0],node_coords[element[1]-1][0],node_coords[element[2]-1][0],node_coords[element[3]-1][0]])
        y = np.array([node_coords[element[0]-1][1],node_coords[element[1]-1][1],node_coords[element[2]-1][1],node_coords[element[3]-1][1]])
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),facecolor='gray',alpha=0.1))
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),fill=None,edgecolor='gray',alpha=0.8))
    
   
    
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    return ax, node_pts


    



    
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

    subplot_kw = dict(autoscale_on=True)
    fig, ax = plt.subplots(subplot_kw=subplot_kw)
    ax, node_pts = show_msh(ax, node_coords, element_nodes)

    constrained_node_ids, _ = get_pts(ax, node_pts, type='Rectangle',enter_closes=True)

    displacement_constraints = np.empty((msh_params['node']['total_num'],3))*np.nan   
    displacement_constraints[:,0] = node_ids
    if constrained_node_ids is not None:
        for node_id in constrained_node_ids:
            displacement_constraints[node_id-1][DoF+1]=value
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
    mesh_dir = "C:\\Users\\ppzmis\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    #mesh_dir = "C:\\Users\\mikei\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_filename = 'rect.msh'
    #constraint_type = 'Displacement'
    constraint_type = 'Force'


    node_ids, node_coords, element_ids, element_nodes, msh_params = read_gmsh(mesh_dir + mesh_filename)
    

    if constraint_type == 'Displacement':
        values = set_constraint(node_ids, node_coords, element_nodes, DoF='x',value=0)
        values2 = set_constraint(node_ids, node_coords, element_nodes, DoF='y',value=0)
        constraints = combine_constraints(values)
        constraints = combine_constraints(values2, total=constraints)
    else:
        values = set_constraint(node_ids, node_coords, element_nodes, DoF='x',value=0)
        values2 = set_constraint(node_ids, node_coords, element_nodes, DoF='y',value=-20)
        constraints = combine_constraints(values)
        constraints = combine_constraints(values2, total=constraints)
      

    np.savetxt(mesh_dir + mesh_filename[:-4] + '_' + constraint_type + '.constraint', constraints)

     

    
    
    


    #Support X


    #Support Y

   
   
