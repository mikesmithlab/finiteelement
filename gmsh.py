import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from DataAnalysis.select_data import get_pts

def read_gmsh(mesh_file,nodes_in_element=4):
    """
    read_gmsh takes a .msh file created using GMSH as an input which specifies
    a starting mesh. 

    Inputs: 
    mesh_file   :   full path to .msh file
    nodes_in_element    : kwarg defining how many nodes make up an element
    
    It returns:
    node_ids : A 1D numpy array of nodal ids [0,1,2,3...]
    node_coords : A 2D numpy array of nodal coords  [[x,y,z],[x2,y2,z2]...]
    element_ids : A 1D numpy array of element ids  [ 0,1,2,3...]
    element_nodes   : A 2D numpy array of element node ids  [[1,4,5,7], [3,4,6,7]...]
    """

    with open(mesh_dir + mesh_filename, 'r') as f:
        msh = f.read()

    node_data = msh.split('$Nodes')[1].split('$EndNodes')[0].split('\n')[1:-1]
    element_data = msh.split('$Elements')[1].split('$EndElements')[0].split('\n')[1:-1]
    msh_params = {}
    msh_params['node'] = {}
    msh_params['element'] = {}
    
    #Extract node data
    msh_params['node']['total_num'] = int(node_data[0].split(' ')[1])
    msh_params['node']['id_min'] = int(node_data[0].split(' ')[2])
    msh_params['node']['id_max'] = int(node_data[0].split(' ')[3])

    node_ids = []
    node_coords = []
    for i in range(1,len(node_data)):
        num_bits = len(node_data[i].rstrip(' ').split(' '))

        if num_bits == 1:
            node_ids.append(int(node_data[i]))
        elif num_bits == 3:
            node_coords.append([float(coord) for coord in node_data[i].rstrip(' ').split(' ')])

    node_ids = np.array(node_ids)
    node_coords = np.array(node_coords)

    #Extract element data
    msh_params['element']['total_num'] = int(element_data[0].split(' ')[1])
    msh_params['element']['id_min'] = int(element_data[0].split(' ')[2])
    msh_params['element']['id_max'] = int(element_data[0].split(' ')[3])

    element_ids = []
    element_nodes = []
    for i in range(1,len(element_data)):
        num_bits = len(element_data[i].rstrip(' ').split(' '))

        if num_bits == (nodes_in_element + 1):
            single_element = element_data[i].rstrip(' ').split(' ')
            element_ids.append(int(single_element[0]))
            element_nodes.append([int(node_id) for node_id in single_element[1:]])

    element_ids = np.array(element_ids)
    element_nodes = np.array(element_nodes)

    return node_ids, node_coords, element_ids, element_nodes, msh_params

def show_msh(node_ids, node_coords, element_ids, element_nodes):
    x = node_coords[:,0]
    y = node_coords[:,1]

    fig, ax = plt.subplots()
    fig.gca().set_aspect('equal',adjustable='box')

    node_pts = ax.scatter(x,y,s=5,c='r')
    for n, element in enumerate(element_nodes):
        x = np.array([node_coords[element[0]-1][0],node_coords[element[1]-1][0],node_coords[element[2]-1][0],node_coords[element[3]-1][0]])
        y = np.array([node_coords[element[0]-1][1],node_coords[element[1]-1][1],node_coords[element[2]-1][1],node_coords[element[3]-1][1]])
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),facecolor='gray',alpha=0.2))
        ax.add_patch(patch.Polygon(xy=list(zip(x,y)),fill=None,edgecolor='gray',alpha=1))
    
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    return ax, node_pts


def set_constraint(ax, node_pts,type='Lasso'):
    get_pts(ax, node_pts,type=type,enter_closes=True)

   


if __name__ == '__main__':
    mesh_dir = "C:\\Users\\mikei\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_filename = 'rect.msh'

    node_ids, node_coords, element_ids, element_nodes, msh_params = read_gmsh(mesh_dir + mesh_filename)

    ax, node_pts = show_msh(node_ids, node_coords, element_ids, element_nodes)
    
    constrained_node_id, value = set_constraint(ax, node_pts)

    





