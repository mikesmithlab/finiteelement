import matplotlib.pyplot as plt
import matplotlib.patches as patch
from gmsh import read_gmsh
import numpy as np
from displacements import calc_new_node_coords

def show_msh(ax, node_coords, element_nodes, show_nodes=True, color='gray',alpha=0.2,wire_alpha=0.8, title='',margin=0.05):
    #Plot node pts
    if show_nodes:
        node_pts = ax.scatter(node_coords[:,0],node_coords[:,1],s=5,c='b')
    else:
        node_pts = None

    #Plot gray patch for mesh
    for n, element in enumerate(element_nodes):
        x = np.array([node_coords[element[0]-1][0],node_coords[element[1]-1][0],node_coords[element[2]-1][0],node_coords[element[3]-1][0]])
        y = np.array([node_coords[element[0]-1][1],node_coords[element[1]-1][1],node_coords[element[2]-1][1],node_coords[element[3]-1][1]])
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),facecolor=color,alpha=alpha))
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),fill=None,edgecolor='gray',alpha=wire_alpha))
    
    minX = node_coords[:,0].min()
    maxX = node_coords[:,0].max()
    minY = node_coords[:,1].min()
    maxY = node_coords[:,1].max()

    ax.set_title(title)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')

    dx = margin*(maxX - minX)
    dy = margin*(maxY - minY)

    ax.set_xlim([minX - dx, dx + maxX])
    ax.set_ylim([minY - dy, dy + maxY])

    return ax, node_pts

def plot_deformation(node_coords, UG, element_nodes, scale=2):
    subplot_kw = dict(autoscale_on=True)
    fig, ax = plt.subplots(subplot_kw=subplot_kw)

    #Plot undeformed structure
    ax, node_pts = show_msh(ax, node_coords, element_nodes, alpha=0.05,wire_alpha=0.3,show_nodes=False)
    
    #Calc new node coords
    new_node_coords = calc_new_node_coords(node_coords,UG,scale=scale)

    #Plot deformed structure
    ax, node_pts = show_msh(ax, new_node_coords, element_nodes, alpha=0.6,wire_alpha=0.9, color='red',show_nodes=False)
    ax.axis('equal')
    plt.show()