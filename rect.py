import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as pth

object_width = 20       #input width of the object
object_height = 10      #input height of the object

k1 = 237     #thermal conductivity of the left half of the object
k2 = 150     #thermal conductivity of the right half of the object

s_x = 4    #has to be even for a convenient split in half
s_y = 4

elem_width = object_width / s_x
elem_height = object_height / s_y
elem_area = elem_width*elem_height
A = 1 / 4 / elem_area

#dictionary to store boundary conditions: left field - boundary type, right field - value of temperature or heat flux
boundary_conditions = {"upper" : ["Dirichlet", 370], "lower" : ["Dirichlet", 320], 
                       "left" : ["Dirichlet", 400], "right" : ["Dirichlet", 300]}


nodes = []

x = np.linspace(0, object_width, s_x + 1)
y = np.linspace(0, object_height, s_y + 1)

for xi in x:
    for yi in y:
        nodes.append((xi, yi))

elements = []       #format: [node1, node2, node3, node4]

for i in range(0, len(nodes)-s_x-1):
    check = i % (s_x+1)   
    if check != s_x:
        elements.append([i, i+1, i+s_x+1, i+s_x+2])

print(elements)

#good up to here 

K = np.zeros((len(nodes), len(nodes)))

def loc_mat_maker(node1, node2, node3, node4):

    loc_mat = [[1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2),1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2),-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2),-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2)],
               [1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2),1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2),-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2),-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2)],
               [-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2),-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2),1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2),1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2)],
               [-1/(6*elem_width*elem_height)*(2*elem_width**2-elem_height**2),-1/(6*elem_width*elem_height)*(elem_width**2+elem_height**2),1/(6*elem_width*elem_height)*(elem_width**2-2*elem_height**2),1/(3*elem_width*elem_height)*(elem_width**2+elem_height**2)]]
    
    n_mat = [node1, node2, node3, node4]
    
    for i in range(4):
        for j in range(4):
            K[n_mat[i]][n_mat[j]] += loc_mat[i][j]

for elem in elements:
    loc_mat_maker(elem[0], elem[1], elem[2], elem[3])

F = np.zeros(len(nodes))

def boundary_check(bnd_str):
    match bnd_str:
        case "upper":
            ran = range(0, s_x + 1)
        case "lower":
            ran = range(-s_x - 1, -0)
        case "left":
            ran = [i for i in range(len(nodes)) if i % (s_y + 1) == 0]
        case "right":
            ran = [i for i in range(len(nodes)) if i % (s_y + 1) == s_y]
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
    F[0] = (F[1] + F[s_y + 1]) / 2
    F[s_y] = (F[2 * s_y + 1] + F[s_y - 1]) / 2
    F[-s_y - 1] = (F[-s_y] + F[-2 * s_y - 2]) / 2
    F[-1] = (F[-1] + F[-s_y + 1]) / 2

for i, row in enumerate(K):
    for j in range(len(row)):
        j_mod = j % (s_x + 1)
        if j_mod < (s_x / 2):
            K[i][j] = K[i][j] * k1
        elif j_mod > (s_x / 2):
            K[i][j] = K[i][j] * k2
        else:
            K[i][j] = K[i][j] * (k1 + k2) / 2

import pandas as pd

## convert your array into a dataframe
df = pd.DataFrame (K)

## save to xlsx file

filepath = 'K_mat.xlsx'

df.to_excel(filepath, index=False)

boundary_check("upper")
boundary_check("lower")
boundary_check("left")
boundary_check("right")

print(F)


T = np.linalg.solve(K, F)
print(T)
