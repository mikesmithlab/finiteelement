
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from DataAnalysis.select_data import get_pts
from gmsh import read_gmsh









        
    
def show_msh(ax, node_coords, element_nodes):
    x = node_coords[:,0]
    y = node_coords[:,1]
    
    #self.ax.set_aspect('equal',adjustable='box')

    node_pts = ax.scatter(x,y,s=5,c='r')
    for n, element in enumerate(element_nodes):
        x = np.array([node_coords[element[0]-1][0],node_coords[element[1]-1][0],node_coords[element[2]-1][0],node_coords[element[3]-1][0]])
        y = np.array([node_coords[element[0]-1][1],node_coords[element[1]-1][1],node_coords[element[2]-1][1],node_coords[element[3]-1][1]])
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),facecolor='gray',alpha=0.2))
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),fill=None,edgecolor='gray',alpha=1))
    
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    return ax, node_pts


    



    
def set_constraint(node_coords, element_nodes, DoF=None, displacement=None):
    """
    DoF - 0 = x, 1=y
    """
    subplot_kw = dict(autoscale_on=True)
    fig, ax = plt.subplots(subplot_kw=subplot_kw)
    ax, node_pts = show_msh(ax, node_coords, element_nodes)

    constrained_node_ids, _ = get_pts(ax, node_pts, type='Rectangle',enter_closes=True)

    displacement_constraints = np.empty((msh_params['node']['total_num'],3))*np.nan   
    if constrained_node_ids is not None:
        for node_id in constrained_node_ids:
            displacement_constraints[node_id-1][0]=node_id
        for node_id in constrained_node_ids:
            displacement_constraints[node_id-1][DoF+1]=displacement
    return displacement_constraints

def combine_constraints(constraint1, constraint2):
    pass

if __name__ == '__main__':
    mesh_dir = "C:\\Users\\ppzmis\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    #mesh_dir = "C:\\Users\\mikei\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_filename = 'pdms_stamp.msh'
    node_ids, node_coords, element_ids, element_nodes, msh_params = read_gmsh(mesh_dir + mesh_filename)
    
    values = set_constraint(node_coords, element_nodes, DoF=0,displacement=0)
    print(values[0:300,:])
    values2 = set_constraint(node_coords, element_nodes, DoF=1,displacement=0)
    print(values2)
    new_values = combine_constraints(values, values2)
    print(new_values)

    #Restraints are indicated by a mask. 
    # pin, mask = [1, 1]
    # roller free in x [0, 1]
    # roller free in y [1, 0]

    

    
    
    


    #Support X


    #Support Y

   
   
