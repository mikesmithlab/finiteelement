import numpy as np

def calculate_ke(alpha, sp, x, y, t=None, C=None):
    """
    calculate_ke calculates the local stiffness matrix for a single element
    alpha
    sp
    x, y are nodal coordinates
    """
    ke = np.zeros([8,8])

    #Sum across all 'r' sample points
    for i, r in enumerate(sp):
        for j, s in enumerate(sp):
            #partial derivatives of shape functions
            dh1dr = 0.25*(1+s)
            dh2dr = -0.25*(1+s)
            dh3dr = -0.25*(1-s)
            dh4dr = 0.25*(1-s)

            dh1ds = 0.25*(1+r)
            dh2ds = 0.25*(1-r)
            dh3ds = -0.25*(1-r)
            dh4ds = -0.25*(1+r)

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

            ke = ke + alpha[i]*alpha[j]*t*B.T*C*B*detJ
    return ke


def add_ke2KP(KP, element, nodes, alpha=None, sp=None, t=None, C=None):
    """
    add a single local stiffness matrix to the global stiffness matrix
    """
    node1 = element[0]
    node2 = element[1]
    node3 = element[2]
    node4 = element[3]

    #x coords
    x = np.array([
                nodes[node1 - 1,0],
                nodes[node2 - 1,0],
                nodes[node3 - 1,0],
                nodes[node4 - 1,0]
                ])

    #y coords
    y = np.array([
                nodes[node1 - 1,1],
                nodes[node2 - 1,1],
                nodes[node3 - 1,1],
                nodes[node4 - 1,1]
                ])

    ke = calculate_ke(alpha, sp, x, y, t=t, C=C)
    
    #Indices for primary stiffness matrix KE
    i1x = 2*node1 - 2
    i1y = 2*node1 - 1
    i2x = 2*node2 - 2
    i2y = 2*node2 - 1
    i3x = 2*node3 - 2
    i3y = 2*node3 - 1
    i4x = 2*node4 - 2
    i4y = 2*node4 - 1

    indices = [i1x,i1y,i2x,i2y,i3x,i3y,i4x,i4y]
    indexArray = np.ix_(indices, indices)
    KP[indexArray]=KP[indexArray]+ke

    return KP


def build_KP(nodes, elements, alpha=None, sp=None, t=None, C=None):
    """
    build_KE calculates the global stiffness matrix of the entire structure
    nodes
    elements
    """
    num_nodes = np.shape(nodes)[0]
    KP = np.zeros([2*num_nodes, 2*num_nodes])

    for n, element in enumerate(elements):
        KP = add_ke2KP(KP, element, nodes, alpha=alpha, sp=sp, t=t, C=C)

    return KP
    
def get_stiffness_matrix(KP, disp_ids):
    """
    Deletes these rows and columns from global stiffness matrix KP
    """
    KS = np.delete(KP, disp_ids, axis=0)
    KS = np.delete(KS, disp_ids, axis=1)
    KS = np.matrix(KS)
    return KS