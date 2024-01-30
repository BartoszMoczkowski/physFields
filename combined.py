import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as pth

elem_type = "triangular"    #choose between triangular and rectangular elements

object_width = 20       #input width of the object
object_height = 10      #input height of the object

k1 = 237     #thermal conductivity of the left half of the object
k2 = 150     #thermal conductivity of the right half of the object

#no. of elements in a row/column + 1 (i.e if you want four elements in a row set s_x to 5)
s_x = 16    #has to be even for a convenient split in half
s_y = 16

elem_width = object_width / s_x
elem_height = object_height / s_y

#dictionary to store boundary conditions: left field - boundary type, right field - value of temperature or heat flux
boundary_conditions = {"upper" : ["Neumann", 0], "lower" : ["Neumann", 0], 
                       "left" : ["Dirichlet", 400], "right" : ["Dirichlet", 300]}

#creates a list of all the nodes with their x, y coords and saves them to the specified matrix
def generate_nodes(mat):
    x = np.linspace(0, object_width, s_x + 1)
    y = np.linspace(0, object_height, s_y + 1)

    for xi in x:
        for yi in y:
            mat.append((xi, yi))

#makes a list of elements used in the method, in the form of a list of node indexes
def generate_elems(mat, type, node_mat):
    if type == "rectangular":
        for i in range(0, len(node_mat)-s_x-1):
            check = i % (s_x + 1)   
            if check != s_x:
                mat.append([i, i+1, i+s_x+1, i+s_x+2])
    else:
        for i in range(0, len(node_mat) - 2 - s_x, 2):
            check = i % (2 * (s_x + 1))
            if check <= s_x:
                if check != 0:
                    mat.append([i + s_x, i + s_x + 1, i])
                    mat.append([i + s_x, i, i - 1])
                if check != s_x:
                    mat.append([i + s_x + 1, i + s_x + 2, i])
                    mat.append([i + s_x + 2, i + 1, i])
        
            else:
                mat.append([i + s_x, i + s_x + 1, i])
                mat.append([i + s_x, i, i - 1])
                mat.append([i + s_x + 1, i + s_x + 2, i])
                mat.append([i + s_x + 2, i + 1, i])

#makes an iteration of the K matrix for the specific rectangular element
def rect_mat_maker(node1, node2, node3, node4, K):
    loc_mat = [[1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2),1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2),-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2),-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2)],
               [1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2),1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2),-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2),-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2)],
               [-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2),-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2),1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2),1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2)],
               [-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2),-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2),1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2),1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2)]]
    
    n_mat = [node1, node2, node3, node4]
    
    for i in range(4):
        for j in range(4):
            K[n_mat[i]][n_mat[j]] += loc_mat[i][j]
    return K

#makes an iteration of the K matrix for a specific triangular element
def tri_mat_maker(node1, node2, node3, node_mat, K):
    elem_area = np.sqrt(elem_width ** 2 + elem_height ** 2)
    A = 1 / 4 / elem_area
    b1 = node_mat[node2][1] - node_mat[node3][1]
    b2 = node_mat[node3][1] - node_mat[node1][1]
    b3 = node_mat[node1][1] - node_mat[node2][1]
    c1 = node_mat[node3][0] - node_mat[node2][0]
    c2 = node_mat[node1][0] - node_mat[node3][0]
    c3 = node_mat[node2][0] - node_mat[node1][0]

    loc_mat = [[A*(b1**2+c1**2), A*(b1*b2+c1*c2), A*(b1*b3+c1*c3)], 
               [A*(b1*b2+c1*c2), A*(b2**2+c2**2), A*(b2*b3+c2*c3)],
               [A*(b1*b3+c1*c3), A*(b2*b3+c2*c3), A*(b3**2+c3**2)]]
    
    n_mat = [node1, node2, node3]
    
    for i in range(3):
        for j in range(3):
            K[n_mat[i]][n_mat[j]] += loc_mat[i][j]
    return K

#adjusts the F and K matrix to account for all the boundary conditions accordingly
def boundary_check(bnd_str, F, K, n_mat_length):
    match bnd_str:
        case "upper":
            ran = range(0, s_x + 1)
        case "lower":
            ran = range(-s_x - 1, -0)
        case "left":
            ran = [i for i in range(n_mat_length) if i % (s_x + 1) == 0]
        case "right":
            ran = [i for i in range(n_mat_length) if i % (s_x + 1) == s_x]
        case _:
            print("Wrong boundary condition")
    
    if boundary_conditions[bnd_str][0] == "Dirichlet":
        for i in ran:
            F[i] = boundary_conditions[bnd_str][1]
            K[i] = np.zeros(K[i].shape)
            K[i][i] = 1
    elif boundary_conditions[bnd_str][0] == "Neumann":
        for i in ran:
            F[i] = F[i] + boundary_conditions[bnd_str][1]
    return F, K

def corner_tweak(F):
    if boundary_conditions["upper"][1] == "Dirichlet" and boundary_conditions["left"] == "Dirichlet":
        F[0] = (F[1] + F[s_x + 1]) / 2
    if boundary_conditions["upper"][1] == "Dirichlet" and boundary_conditions["right"] == "Dirichlet":
        F[s_x] = (F[2 * s_x + 1] + F[s_x - 1]) / 2
    if boundary_conditions["lower"][1] == "Dirichlet" and boundary_conditions["left"] == "Dirichlet":
        F[-s_x - 1] = (F[-s_x] + F[-2 * s_x - 2]) / 2
    if boundary_conditions["lower"][1] == "Dirichlet" and boundary_conditions["right"] == "Dirichlet":
        F[-1] = (F[-1] + F[-s_x + 1]) / 2
    return F

#executes all the functions in the required order and returns the temperature vector
def main():
    nodes = []
    elems = []

    generate_nodes(nodes)
    generate_elems(elems, elem_type, nodes)

    K = np.zeros((len(nodes), len(nodes)))

    if elem_type == "rectangular":
        for e in elems:
            K = rect_mat_maker(e[0], e[1], e[2], e[3], K)
    else:
        for e in elems:
            K = tri_mat_maker(e[0], e[1], e[2], nodes, K)

    F = np.zeros(len(nodes))

    for i, row in enumerate(K):
        for j in range(len(row)):
            j_mod = i % (s_y + 1)
            if j_mod < (s_y / 2):
                K[i][j] = K[i][j] * k1
            elif j_mod > (s_y / 2):
                K[i][j] = K[i][j] * k2
            else:
                K[i][j] = K[i][j] * (k1 + k2) / 2
    
    n_mat_length = len(nodes)
    F, K = boundary_check("upper", F, K, n_mat_length)
    F, K = boundary_check("lower", F, K, n_mat_length)
    F, K = boundary_check("left", F, K, n_mat_length)
    F, K = boundary_check("right", F, K, n_mat_length)

    F = corner_tweak(F)

    T = np.linalg.solve(K, F)
    return T

T = main()
plt.matshow(T.reshape(17, 17))
plt.show()
