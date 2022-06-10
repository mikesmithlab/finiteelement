from gmsh import read_gmsh
import numpy as np
from stiffness_matrix import build_KP, get_stiffness_matrix
from displacements import restrained_DoF, calc_global_disp
from forces import get_reduced_force_vector
from visualise import plot_deformation

def text_output(disp_ids, FG):
    print(disp_ids)
    print('Reactions')
    print('----------------------------------')
    for id in disp_ids:
        if id%2!=0:
            Ry = FG[id]
            print('Ry at node {id} = {Ry} kN'.format(id=id,Ry=np.round(Ry/1000,2)))
        else:
            Rx = FG[id]
            print('Rx at node {id} = {Rx} kN'.format(id=id,Rx=np.round(Rx/1000,2)))
        

    

if __name__ == '__main__':
    #mesh_dir = "C:\\Users\\ppzmis\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_dir = "C:\\Users\\mikei\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_filename = 'rect.msh'
    
    """
    Fundamental constants / parameter setting
    """
    E = 200E9 # Young's Modulus  N/m^2
    nu = 0.3 # Poisson's ratio

    #Plane Stress Material Matrix
    C = (E/(1-nu**2))*np.array([[1, nu, 0],[nu, 1, 0],[0, 0, (1-nu)/2]])
    t= 0.1 # Thickness beam must be small

    #Plane Strain Material Matrix
    #C = ((E*(1 - nu))) / ((1 + nu) * (1 - 2*nu)) * np.array([[1, nu / (1 - nu), 0],[nu / (1 - nu), 1, 0],[0, 0, (1 - 2*nu)/(2*(1-nu))]])
    #t = 1

    #Gauss Parameters
    alpha = [1, 1]
    sp = [-0.5773502692, 0.5773502692] # Sampling points

    """End parameter setting"""

    node_ids, node_coords, element_ids, element_nodes, msh_params = read_gmsh(mesh_dir + mesh_filename)

    disp = np.loadtxt(mesh_dir + mesh_filename[:-4] + '_Displacement.constraint')
    force = np.nan_to_num(np.loadtxt(mesh_dir + mesh_filename[:-4] + '_Force.constraint'), copy=False, nan=0.0)
    print(force)
    #Build the force vector
    F = np.reshape(force[:,1:],(2*np.shape(force)[0],1))
    
    #Global Stiffness Matrix
    KP = build_KP(node_coords, element_nodes, alpha=alpha, sp=sp, t=t, C=C)

    #Find restrained DoF
    disp_ids = restrained_DoF(disp)
   

    #Extract Structure Stiffness Matrix - Delete rows and columns
    KS = get_stiffness_matrix(KP, disp_ids)

    #Construct reduced force vector
    F_red = get_reduced_force_vector(F, disp_ids)

    #Solve for structure
    U = KS.I * F_red

    #Get global displacement vector
    UG=calc_global_disp(U, disp_ids, disp)

    #Global force vector
    FG=np.matmul(KP,UG)

    #Display deflection
    plot_deformation(node_coords, UG, element_nodes,scale=1)

    #Text output
    text_output(disp_ids, FG)

    

