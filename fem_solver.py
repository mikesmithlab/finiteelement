from gmsh import read_gmsh
import numpy as np

def calculate_ke(alpha, sp, x, y):
        KE = np.zeros([8,8])

        #Sum across all 'r' sample points
        for i, r in enumerate(sp):
            for j, s in enumerate(sp):
                #partial derivatives of shape functions
                dh1dr = 0.25*(1+s)
                dh2dr = -0.25*(1-s)
                dh3dr = -0.25*(1-s)
                dh4dr = 0.25*(1+s)

                dh1ds = 0.25*(1+r)
                dh2ds = -0.25*(1-r)
                dh3ds = -0.25*(1-r)
                dh4ds = 0.25*(1+r)

                #Components of Jacobian Matrix
                dxdr = x[0]*dh1dr + x[1]*dh2dr + x[2]*dh3dr + x[3]*dh4dr
                dxds = x[0]*dh1ds + x[1]*dh2ds + x[2]*dh3ds + x[3]*dh4ds
                dydr = y[0]*dh1dr + y[1]*dh2dr + y[2]*dh3dr + y[3]*dh4dr
                dyds = y[0]*dh1ds + y[1]*dh2ds + y[2]*dh3ds + y[3]*dh4ds

                #Build Jacobian and inverse
                J = np.matrix([[dxdr, dydr],[dxds, dyds]])
                invJ = J.I
                detJ = np.linalg.det(J)

                #Compile matrices of partial derivs of H1 and H2
                dH1 = np.matrix([[dh1dr, 0,dh2dr, 0,dh3dr, 0,dh4dr, 0],[dh1ds, 0,dh2ds, 0,dh3ds, 0,dh4ds, 0]])
                dH2 = np.matrix([[0, dh1dr, 0,dh2dr, 0,dh3dr, 0,dh4dr],[0,dh1ds, 0,dh2ds, 0,dh3ds, 0,dh4ds]])

                #Calculate the Strain Displacement Matrix
                B = np.matrix([[1,0],[0,0],[0,1]])*invJ*dH1 + np.matrix([[0,0],[0,1],[1,0]])*invJ*dH2

                KE = KE + alpha[i]*alpha[j]*t*B.T*C*B*detJ
        return KE


def add_localKE(KE, ke):
    pass


def build_KE(elements):
    num_elements = np.shape(elements)[0]
    KE = np.zeros([num_elements, num_elements])

    for n, element in enumerate(elements):

        node1 = element[0]
        node2 = element[1]
        node3 = element[2]
        node4 = element[3]

        




if __name__ == '__main__':
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

    







    mesh_dir = "C:\\Users\\ppzmis\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    #mesh_dir = "C:\\Users\\mikei\\OneDrive - The University of Nottingham\\Documents\\FEM\\MeshFiles\\"
    mesh_filename = 'rect.msh'

    node_ids, node_coords, element_ids, element_nodes, msh_params = read_gmsh(mesh_dir + mesh_filename)

    disp = np.loadtxt(mesh_dir + mesh_filename[:-4] + '_Displacement.constraint')
    force = np.nan_to_num(np.loadtxt(mesh_dir + mesh_filename[:-4] + '_Force.constraint'), copy=False, nan=0.0)



    #Assign Point Loads
    # He uses P for magnitude and pointLoadAxis to indicate direction
    
    #Build the global force vector
    forceVector = np.reshape(force[:,1:],(2*np.shape(force)[0],1))
    
    #Element Stiffness Matrix
    