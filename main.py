import numpy as np
import scipy.linalg as lin
import networkx as nx
from io import StringIO

alpha = 0
beta = -1

def lin_polyene_mat(n_orbital):
    mat = np.zeros([n_orbital,n_orbital])
    for i in range(n_orbital-1):
        mat[i,i] = alpha
        mat[i,i+1] = beta
    mat += mat.T
    return mat

def cyc_polyene_mat(n_orbital):
    mat = np.zeros([n_orbital,n_orbital])
    for i in range(n_orbital-1):
        mat[i,i] = alpha
        mat[i,i+1] = beta
    mat[0,n_orbital-1] = beta
    mat += mat.T
    return mat

def degen_check(e_values):
    output = ""
    degen_dict = {}
    for i in range(len(e_values)): # Using a sliding window to check the sorted eigenvalues for repeats
        e_value1 = round(e_values[i],5) # first value of sliding window

        if i == len(e_values)-1: # if reach end of e_values list, make a dictionary item if there isnt any key 
            degen_dict.setdefault(e_value1,1)
            break 

        e_value2 = round(e_values[i+1],5) # second value of sliding window

        if e_value1 == e_value2:
            degen_dict[e_value1] = degen_dict.setdefault(e_value1,1) + 1 # if values match, add one to key. if key does not exist, add the key and set the default value to 1 and add 1
        else:
            degen_dict.setdefault(e_value1,1) # if values does not match do nothing. If key does not exist, set the default value as 1.
    
    output = "Energy      | Degeneracies\n" \
             "--------------------------\n"
    for key, value in reversed(degen_dict.items()): # generate a string to display the energies and degeneracies
        if key*-1 >= 0:
            output += f'α+{str(abs(-1*key))+'β':<9} | {value}\n'
        else:
            output += f'α-{str(abs(-1*key))+'β':<9} | {value}\n'
    print(output)
    return output

def save_output(text): # save output to text
    file_name = "output.txt"
    f = open(file_name, 'w', encoding="utf-8")
    f.write(text)
    f.close()


def tetrahedron_mat():
    mat = np.ones([4,4])
    mat = mat*beta
    for i in range(4):
        mat[i,i] = alpha
    return mat

def cube_mat():
    mat = np.zeros([8,8])
    for i in range(8):
        mat[i,i] = alpha
    for i in range(7):
        mat[i,i+1] = beta
    mat[0,7] = mat[0,5] = mat[1,4] = mat[2,7] = mat[3,6] = beta
    mat += mat.T
    return mat

def dodecahedron_mat():
    dodeca_graph = nx.dodecahedral_graph()
    mat = nx.adjacency_matrix(dodeca_graph)
    mat = mat.toarray()
    mat = np.array(mat)*beta
    return mat

def octahedron_mat():
    octa_graph = nx.octahedral_graph()
    mat = nx.adjacency_matrix(octa_graph)
    mat = mat.toarray()
    mat = np.array(mat)*beta
    return mat

def icosahedron_mat():
    icosa_graph = nx.icosahedral_graph()
    mat = nx.adjacency_matrix(icosa_graph)
    mat = mat.toarray()
    mat = np.array(mat)*beta
    return mat

def buckyball_mat():
    with open("buckyball.txt", 'r') as f:
        buckyball_text = f.read()
        mat = np.genfromtxt(StringIO(buckyball_text), delimiter=" ")
        mat = mat*beta
        return mat
    
def solve_mat(mat):
    evals, evects = lin.eigh(mat)
    output = degen_check(evals)
    save_output(output) 
    # print(evects)



print("Select a structure to determine the molecular orbital energies:")
print('''1. linear polyene
2. cyclic polyene
3. tetrahedron
4. cube
5. dodecahedron
6. octahedron
7. icosahedron
8. buckministerfullerene

Type the selection number of the structure:''')

selection = ""
orbital_num_selection = ""
input_success = False
function_list = [lin_polyene_mat, cyc_polyene_mat, tetrahedron_mat, cube_mat, 
                 dodecahedron_mat, octahedron_mat, icosahedron_mat, buckyball_mat]

while not input_success:
    selection = input()
    if not selection.isdigit():
        print("Please enter a valid number!")
    elif 0 < int(selection) <= 8:
        selection = int(selection)
        input_success = True
    else:
        print("Please enter a valid number!")

if selection == 1 or selection == 2: # check if linear or cyclic polyenes is selected and then ask for the number of orbitals in the system
    input_success = False
    print("Please enter the number of atoms/orbitals in the polyene")
    while not orbital_num_selection.isdigit():
        orbital_num_selection = input()
        if not orbital_num_selection.isdigit():
            print("Please enter a valid integer!")
    solve_mat(function_list[selection-1](int(orbital_num_selection)))
else:
    solve_mat(function_list[selection-1]())



