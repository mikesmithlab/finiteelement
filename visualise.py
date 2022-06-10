import matplotlib.pyplot as plt
import matplotlib.patches as patch
from gmsh import read_gmsh
import numpy as np
from displacements import calc_new_node_coords

def show_msh(ax, node_coords, element_nodes, show_nodes=True, title=''):
    #Plot node pts
    if show_nodes:
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

def plot_deformation(scale=2):