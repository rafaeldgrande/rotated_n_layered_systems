
import numpy as np
import sys

a = float(sys.argv[1])  #2.46  # Parametro de rede
N = int(sys.argv[2])  # index n
M = int(sys.argv[3])  # index m
arq_1st_layer = sys.argv[4]  # File for first layer
arq_2nd_layer = sys.argv[5]  # File for second layer

reinforce_distance_layers = True  # If True, adds a small offset to the second layer to ensure separation
d_interlayer = 6.0 # Distance between layers in angstroms (distance between highest values of z of the first layer and lowest values of z of the second layer)

'''Usage:
python generate_twisted_2L_WSe2.py a N M first_layer_file second_layer_file

File for each layer should be like this:
W x y z
Se x y z
Se x y z

'''

def inside_supercell(x, y, L_supercell):  # Verifica se o vetor posicao ta dentro da supercelula
    # (x ,y) = alfa*L1 + beta*L2, onde L1 = L(1, 0) e L2 = L(0.5, sqrt(3)/2)
    alfa, beta = (x-y/np.sqrt(3))/L_supercell, 2*y/(L_supercell*np.sqrt(3))
    if 0 <= alfa < 1 - delta and 0 <= beta < 1 - delta:
        return True
    else:
        return False
    
def get_prim_cell_coordinates(file_coord):
    arq = open(file_coord, 'r')

    coordinates = []
    elements = []
    for line in arq:
        line_split = line.split()
        if len(line_split) == 4:
            element = line_split[0]
            x, y, z = map(float, line_split[1:])
            coordinates.append(np.array([x, y, z]))
            elements.append(element)
    arq.close()

    return np.array(coordinates), elements

def get_element_indices(elements_list):
    unique_elements = []
    for elem in elements_list:
        if elem not in unique_elements:
            unique_elements.append(elem)
    element_to_index = {elem: idx+1 for idx, elem in enumerate(unique_elements)}
    indexed_list = [element_to_index[elem] for elem in elements_list]
    return element_to_index, indexed_list
    

theta = np.arccos((N**2 + 4*N*M + M**2)/(2*(N**2 + N*M + M**2)))
L_supercell = a*np.sqrt(N**2+N*M+M**2)
L = L_supercell/abs(M-N)  # Periodicidade do padrao de Moire

Expected_number_prim_cells = (N**2 + N*M + M**2)  # Numero esperado de celulas primarias na supercelula
Expected_number_atoms_per_layer = Expected_number_prim_cells * 3  # Numero esperado de atomos por camada (3 atomos por celula primaria)
print(f"Predicted number of atoms per layer 3 * (N**2 + N*M + M**2): {Expected_number_atoms_per_layer}")

delta = 1e-3*a/L_supercell

# Cell definition in LAMMPS
xhi = L_supercell
yhi = L_supercell*np.sin(np.pi/3)
xy = L_supercell*np.cos(np.pi/3)


# Hexagonal lattice vectors. Used in quantum espresso
a1 = a*np.array([np.sqrt(3)/2, 0.5, 0])
a2 = a*np.array([np.sqrt(3)/2, -0.5, 0])

# Base - supercelula
A1 = L_supercell*np.array([1, 0, 0])
A2 = L_supercell*np.array([0.5, np.sqrt(3)/2.0, 0])

print(f"Rotation angle (theta): {np.degrees(theta):.3f} degrees")
print(f"Supercell size: {L_supercell} Angstroms")
print(f"Supercel lattice vectors A1: {A1}, A2: {A2}")


# Reading first and second layer coordinates
print(f"Reading coordinates from files {arq_1st_layer} and {arq_2nd_layer}")
coords_1st_layer, elements_1st_layer = get_prim_cell_coordinates(arq_1st_layer)
coords_2nd_layer, elements_2nd_layer = get_prim_cell_coordinates(arq_2nd_layer)

# print zmin and zmax for each layer
zmin_1st = np.min(coords_1st_layer[:, 2])
zmax_1st = np.max(coords_1st_layer[:, 2])
zmin_2nd = np.min(coords_2nd_layer[:, 2])
zmax_2nd = np.max(coords_2nd_layer[:, 2])
print(f"Layer 1: zmin = {zmin_1st}, zmax = {zmax_1st}")
print(f"Layer 2: zmin = {zmin_2nd}, zmax = {zmax_2nd}")
print(f"Original interlayer distance: {zmin_2nd - zmax_1st} Angstroms")

print(f"Total atoms in first layer: {len(coords_1st_layer)}")
print(f"Total atoms in second layer: {len(coords_2nd_layer)}")

if reinforce_distance_layers:
    print(f"Reinforcing interlayer distance by {d_interlayer} Angstroms")
    coords_2nd_layer[:, 2] += (zmax_1st - zmin_2nd + d_interlayer) 


rmax = L_supercell * np.sqrt(3)  # generous overestimate of radius
max_ind = int(np.ceil(rmax / np.linalg.norm(a1)))

def rotate_coordinates(coords, rot_angle, elements):
    moire_coords = []
    moire_elements = []
    for i in range(-max_ind, max_ind):
        for j in range(-max_ind, max_ind):
            for idx_atom, atomic_coord in enumerate(coords):  # loop over coords:
                x, y, z = atomic_coord + i*a1 + j*a2
                
                # Rotating coordinates
                x_rot = x*np.cos(rot_angle/2.0) + y*np.sin(rot_angle/2.0)
                y_rot = -x*np.sin(rot_angle/2.0) + y*np.cos(rot_angle/2.0)

                if inside_supercell(x_rot, y_rot, L_supercell):
                    moire_coords.append(np.array([x_rot, y_rot, z]))
                    moire_elements.append(elements[idx_atom])

    return np.array(moire_coords), moire_elements


# rotating the coordinates of the first and second layers
coords_1st_twisted_layer, elements_1st_twisted_layer = rotate_coordinates(coords_1st_layer, theta, elements_1st_layer)
coords_2nd_twisted_layer, elements_2nd_twisted_layer = rotate_coordinates(coords_2nd_layer, -theta, elements_2nd_layer)

print(f"Total atoms in moire first layer: {len(coords_1st_twisted_layer)}")
print(f"Total atoms in moire second layer: {len(coords_2nd_twisted_layer)}")

if len(coords_1st_twisted_layer) != Expected_number_atoms_per_layer:
    print("WARNING: Number of atoms in first twisted layer does not match expected number.")
if len(coords_2nd_twisted_layer) != Expected_number_atoms_per_layer:
    print("WARNING: Number of atoms in second twisted layer does not match expected number.")

Total_atoms = len(coords_1st_twisted_layer) + len(coords_2nd_twisted_layer)

# mapping elements to indices
element_to_index_1st, indexed_elements_1st = get_element_indices(elements_1st_twisted_layer)
element_to_index_2nd, indexed_elements_2nd = get_element_indices(elements_2nd_twisted_layer)

indexed_elements_2nd = list(np.array(indexed_elements_2nd) + max(indexed_elements_1st))

Nelements = max(indexed_elements_2nd)

print("Writing coordinates to lammps files...")

topo_lammps = open(str(N)+'_'+str(M)+'_twisted.top', 'w')
xyz_file = open(str(N)+'_'+str(M)+'_twisted.xyz', 'w')
xyz_file.write(f"{Total_atoms}\n\n")

topo_lammps.write('\n\n# indices '+str(N)+' '+str(M)+'\n')
topo_lammps.write('# '+str(L_supercell)+' angstrons - supercell size\n')
topo_lammps.write('# '+str((180.0/np.pi)*theta)+' degrees - rotation angle\n')
topo_lammps.write(f'\n\n{Total_atoms} atoms\n')
topo_lammps.write(f'{Nelements}  atom types \n\n')

topo_lammps.write('0.0 '+str(xhi)+' xlo xhi \n')
topo_lammps.write('0.0 '+str(yhi)+' ylo yhi \n')
topo_lammps.write('-50.0 50.0 zlo zhi \n')
topo_lammps.write(str(xy)+' 0.0 0.0 xy xz yz \n\n')
topo_lammps.write('Atoms \n\n')

# 1st layer atoms
atom_index = 0
mol_index = 1  
for i, (coord, elem_index) in enumerate(zip(coords_1st_twisted_layer, indexed_elements_1st)):
    x, y, z = coord
    atom_index += 1
    # atom_index molecule_index atom_type x y z
    # topo_lammps.write(f"{i+1} 1 {elem_index} {x} {y} {z}\n")
    topo_lammps.write(f"{atom_index:6d} {mol_index:4d} {elem_index:4d} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    
    xyz_file.write(f"{elements_1st_twisted_layer[i]} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    
# 2nd layer atoms
mol_index = 2
for i, (coord, elem_index) in enumerate(zip(coords_2nd_twisted_layer, indexed_elements_2nd)):
    x, y, z = coord
    atom_index += 1
    # atom_index molecule_index atom_type x y z
    # topo_lammps.write(f"{i+1 + len(coords_1st_twisted_layer)} 2 {elem_index} {x} {y} {z}\n")
    topo_lammps.write(f"{atom_index:6d} {mol_index:4d} {elem_index:4d} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    xyz_file.write(f"{elements_2nd_twisted_layer[i]} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    
xyz_file.close()
topo_lammps.close()

print(f"Coordinates written to {str(N)}_{str(M)}_twisted.top")