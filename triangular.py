import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as pth
import pandas as pd

object_width = 20       #input width of the object
object_height = 10      #input height of the object

k1 = 150     #thermal conductivity of the left half of the object
k2 = 237     #thermal conductivity of the right half of the object

mesh_density = 4    #has to be even for a convenient split in half

elem_width = object_width / mesh_density
elem_height = object_height / mesh_density
elem_area = np.sqrt(elem_width ** 2 + elem_height **2)
A = 1 / 4 / elem_area

#dictionary to store boundary conditions: left field - boundary type, right field - value of temperature or heat flux
boundary_conditions = {"upper" : ["Dirichlet", 370], "lower" : ["Dirichlet", 320], 
                       "left" : ["Dirichlet", 400], "right" : ["Dirichlet", 300]}

nodes = []

x = np.linspace(0, object_width, mesh_density + 1)
y = np.linspace(0, object_height, mesh_density + 1)

for xi in x:
    for yi in y:
        nodes.append((xi, yi))

elements = []       #format: [element_type, node1, node2, node3]

for i in range(0, len(nodes) - 2 - mesh_density, 2):
    check = i % (2 * (mesh_density + 1))
    if check <= mesh_density:
        if check != 0:
            elements.append(["Type 3", i + mesh_density, i + mesh_density + 1, i])
            elements.append(["Type 4", i + mesh_density, i, i - 1])
        if check != mesh_density:
            elements.append(["Type 1", i + mesh_density + 1, i + mesh_density + 2, i])
            elements.append(["Type 2", i + mesh_density + 2, i + 1, i])
        
    else:
        elements.append(["Type 3", i + mesh_density, i + mesh_density + 1, i])
        elements.append(["Type 4", i + mesh_density, i, i - 1])
        elements.append(["Type 1", i + mesh_density + 1, i + mesh_density + 2, i])
        elements.append(["Type 2", i + mesh_density + 2, i + 1, i])

#good up to here 

K = np.zeros((len(nodes), len(nodes)))

def loc_mat_maker(node1, node2, node3):
    b1 = nodes[node2][1] - nodes[node3][1]
    b2 = nodes[node3][1] - nodes[node1][1]
    b3 = nodes[node1][1] - nodes[node2][1]
    c1 = nodes[node3][0] - nodes[node2][0]
    c2 = nodes[node1][0] - nodes[node3][0]
    c3 = nodes[node2][0] - nodes[node1][0]

    loc_mat = [[A*(b1**2+c1**2), A*(b1*b2+c1*c2), A*(b1*b3+c1*c3)], 
               [A*(b1*b2+c1*c2), A*(b2**2+c2**2), A*(b2*b3+c2*c3)],
               [A*(b1*b3+c1*c3), A*(b2*b3+c2*c3), A*(b3**2+c3**2)]]
    
    n_mat = [node1, node2, node3]
    
    for i in range(3):
        for j in range(3):
            K[n_mat[i]][n_mat[j]] += loc_mat[i][j]

for elem in elements:
    loc_mat_maker(elem[1], elem[2], elem[3])

F = np.zeros(len(nodes))

def boundary_check(bnd_str):
    match bnd_str:
        case "upper":
            ran = range(0, mesh_density + 1)
        case "lower":
            ran = range(-mesh_density - 1, -0)
        case "left":
            ran = [i for i in range(len(nodes)) if i % (mesh_density + 1) == 0]
        case "right":
            ran = [i for i in range(len(nodes)) if i % (mesh_density + 1) == mesh_density]
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
    F[0] = (F[1] + F[mesh_density + 1]) / 2
    F[mesh_density] = (F[2 * mesh_density + 1] + F[mesh_density - 1]) / 2
    F[-mesh_density - 1] = (F[-mesh_density] + F[-2 * mesh_density - 2]) / 2
    F[-1] = (F[-1] + F[-mesh_density + 1]) / 2

for i, row in enumerate(K):
    for j in range(len(row)):
        j_mod = j % (mesh_density + 1)
        if j_mod < (mesh_density / 2):
            K[i][j] = K[i][j] * k1
        elif j_mod > (mesh_density / 2):
            K[i][j] = K[i][j] * k2
        else:
            K[i][j] = K[i][j] * (k1 + k2) / 2

#for the purposees of viewing the K matrix in a more readabl format
"""
df = pd.DataFrame (K)

filepath = 'K_mat.xlsx'

df.to_excel(filepath, index=False)
"""

boundary_check("upper")
boundary_check("lower")
boundary_check("left")
boundary_check("right")

T = np.linalg.solve(K, F)
print(T)
