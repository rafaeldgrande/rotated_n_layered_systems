
import numpy as np
import sys

# Tolerance for boundary checking - slightly larger to handle floating point errors
delta = 1e-2

a = float(sys.argv[1])  #2.46  # Parametro de rede
N = int(sys.argv[2])  # index n
M = int(sys.argv[3])  # index m
# arq_1st_layer = sys.argv[4]  # File for first layer
# arq_2nd_layer = sys.argv[5]  # File for second layer
d_Se_W = 1.6
d_W_W = 6.5  # Distance between W atoms in different layers (angstroms)

reinforce_distance_layers = True  # If True, adds a small offset to the second layer to ensure separation
d_interlayer = 6.0 # Distance between layers in angstroms (distance between highest values of z of the first layer and lowest values of z of the second layer)


def rotate_coordinates(coords, rot_angle, elements):
    moire_coords = []
    moire_elements = []
    non_moire_coords = []
    non_moire_elements = []
    for i in range(-max_ind, max_ind):
        for j in range(-max_ind, max_ind):
            origin_cell = i*a1 + j*a2
            x_origin_rot = origin_cell[0]*np.cos(rot_angle/2.0) + origin_cell[1]*np.sin(rot_angle/2.0)
            y_origin_rot = -origin_cell[0]*np.sin(rot_angle/2.0) + origin_cell[1]*np.cos(rot_angle/2.0)
            cell_inside_supercell = inside_supercell(x_origin_rot, y_origin_rot, L_supercell)
            for idx_atom, atomic_coord in enumerate(coords):  # loop over coords:
                x, y, z = atomic_coord + origin_cell
                
                # Rotating coordinates
                x_rot = x*np.cos(rot_angle/2.0) + y*np.sin(rot_angle/2.0)
                y_rot = -x*np.sin(rot_angle/2.0) + y*np.cos(rot_angle/2.0)
                
                if cell_inside_supercell:
                    moire_coords.append(np.array([x_rot, y_rot, z]))
                    moire_elements.append(elements[idx_atom])
                else:
                    non_moire_coords.append(np.array([x_rot, y_rot, z]))
                    non_moire_elements.append(elements[idx_atom])

    return np.array(moire_coords), moire_elements, np.array(non_moire_coords), non_moire_elements

def rotate_coordinates_integer_moire(
    coords,
    elements,
    origins,
    rot_angle
):
    """
    Integer-based construction of a twisted layer using precomputed
    primitive-cell origins inside the moiré supercell.

    Parameters
    ----------
    coords : ndarray (Nat, 3)
        Primitive-cell atomic coordinates
    elements : list
        Atomic species for primitive cell
    origins : ndarray (Ncell, 3)
        Primitive-cell origins inside moiré supercell
    rot_angle : float
        Rotation angle (full angle; rotation applied as ±angle/2 externally)

    Returns
    -------
    moire_coords : ndarray (Nat * Ncell, 3)
        Twisted-layer atomic coordinates
    moire_elements : list
        Atomic species list (same ordering)
    """
    cos_t = np.cos(rot_angle / 2.0)
    sin_t = np.sin(rot_angle / 2.0)

    moire_coords = []
    moire_elements = []

    for R in origins:
        for coord, elem in zip(coords, elements):
            x, y, z = coord + R

            # Rotation
            x_rot =  cos_t * x + sin_t * y
            y_rot = -sin_t * x + cos_t * y

            moire_coords.append([x_rot, y_rot, z])
            moire_elements.append(elem)

    return np.array(moire_coords), moire_elements

def cartesian_to_lammps(xc, yc, A_LAMMPS, B_LAMMPS):
    lx = A_LAMMPS[0]
    bx = B_LAMMPS[0]
    by = B_LAMMPS[1]

    v = yc / by
    u = (xc - v * bx) / lx

    # u %= 1.0
    # v %= 1.0

    x_lmp = u * lx + v * bx
    y_lmp = v * by

    return x_lmp, y_lmp

# def inside_supercell(x, y, L_supercell):  # Verifica se o vetor posicao ta dentro da supercelula
#     # (x ,y) = alfa*L1 + beta*L2, onde L1 = L(1, 0) e L2 = L(0.5, sqrt(3)/2)
#     x_red, y_red = x / L_supercell, y / L_supercell
#     alfa, beta = x_red-y_red/np.sqrt(3), 2*y_red/np.sqrt(3)
#     if 0 < alfa <= 1 and 0 < beta <= 1:
#         return True
#     else:
#         return False

def moire_primitive_origins(N, M, a1, a2):
    """
    Fully integer-based enumeration of primitive-cell origins
    inside a commensurate moiré supercell.

    Returns exactly N^2 + N*M + M^2 cells.
    """
    # Supercell vectors in integer (a1, a2) basis
    A1 = np.array([N, M])
    A2 = np.array([-M, N + M])

    det = A1[0]*A2[1] - A1[1]*A2[0]
    assert det == N**2 + N*M + M**2

    origins = []
    ij_indices = []

    # Bounding box in integer space
    imin = min(0, A1[0], A2[0], A1[0]+A2[0])
    imax = max(0, A1[0], A2[0], A1[0]+A2[0])
    jmin = min(0, A1[1], A2[1], A1[1]+A2[1])
    jmax = max(0, A1[1], A2[1], A1[1]+A2[1])

    for i in range(imin, imax):
        for j in range(jmin, jmax):
            # Solve for (alpha, beta) using integer arithmetic
            # [i] = alpha*A1 + beta*A2
            # [j]
            alpha_num =  i*A2[1] - j*A2[0]
            beta_num  = -i*A1[1] + j*A1[0]

            # Check inclusion: 0 <= alpha < 1 and 0 <= beta < 1
            if (0 <= alpha_num < det) and (0 <= beta_num < det):
                R = i*a1 + j*a2
                origins.append(R)
                ij_indices.append((i, j))

    return np.array(origins), ij_indices



    
def inside_supercell(x, y, L_supercell, tol=1e-8):
    # Reduced coordinates in the supercell basis
    x_red = x / L_supercell
    y_red = y / L_supercell

    alfa = x_red - y_red / np.sqrt(3)
    beta = 2.0 * y_red / np.sqrt(3)
    
    if abs(alfa - 1.0) < 1e-6 or abs(beta - 1.0) < 1e-6:
        print("Boundary cell:", alfa, beta)

    return (
        (alfa >= -tol) and (alfa < 1.0 - tol) and
        (beta >= -tol) and (beta < 1.0 - tol)
    )
    
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
L_supercell = np.sqrt(N**2+N*M+M**2)
L = L_supercell/abs(M-N)  # Periodicidade do padrao de Moire

Expected_number_prim_cells = (N**2 + N*M + M**2)  # Numero esperado de celulas primarias na supercelula
Expected_number_atoms_per_layer = Expected_number_prim_cells * 3  # Numero esperado de atomos por camada (3 atomos por celula primaria)
print(f"Expected number of atoms per layer 3 * (N**2 + N*M + M**2): {Expected_number_atoms_per_layer}")
print(f"Expected number of primitive cells in supercell (N**2 + N*M + M**2): {Expected_number_prim_cells}")
print(f"Supercell lattice parameter in lattice units {L_supercell:.4f} ")
print(f"Supercell lattice parameter in angstroms {L_supercell*a:.2f} ")



# # Cell definition in LAMMPS
# xhi = L_supercell
# yhi = L_supercell*np.sin(np.pi/3)
# xy = L_supercell*np.cos(np.pi/3)


# Hexagonal lattice vectors. Used in quantum espresso
# a1 = a*np.array([np.sqrt(3)/2, 0.5, 0])
# a2 = a*np.array([np.sqrt(3)/2, -0.5, 0])

# a1 = a*np.array([np.sqrt(3)/2,  1/2, 0])
# a2 = a*np.array([np.sqrt(3)/2, -1/2, 0])

a1 = np.array([1, 0, 0])
a2 = np.array([1/2, np.sqrt(3)/2, 0])

tau_W = np.array([0, 0, 0])
tau_Se1 = (a1 + a2) / 3.0 + np.array([0, 0, d_Se_W/a])
tau_Se2 = (a1 + a2) / 3.0 + np.array([0, 0, -d_Se_W/a])

coords_1st_layer = np.array([
    tau_W,
    tau_Se1,
    tau_Se2
])

coords_2nd_layer = np.array([
    tau_W + np.array([0, 0, d_W_W/a]),
    tau_Se1 + np.array([0, 0, d_W_W/a]),
    tau_Se2 + np.array([0, 0, d_W_W/a])
])

elements_1st_layer = ['W', 'Se', 'Se']
elements_2nd_layer = ['W', 'Se', 'Se']

# # Base - supercelula
# A1 = L_supercell*np.array([1, 0, 0])
# A2 = L_supercell*np.array([0.5, np.sqrt(3)/2.0, 0])


# True moiré supercell vectors
A1 = N*a1 + M*a2
A2 = -M*a1 + (N+M)*a2

A1_cart = a * A1
A2_cart = a * A2

factor = 1.0

xlo, xhi = 0.0, A1_cart[0] * factor  # Slightly larger to avoid boundary issues
ylo, yhi = 0.0, A2_cart[1] * factor
xy = A2_cart[0] * factor

lx = xhi - xlo
ly = yhi - ylo
b = np.sqrt(ly**2 + xy**2)
bx = xy  # = b * cos(gamma)
by = np.sqrt(b**2 - bx**2)  # = b * sin(gamma)
# A_LAMMPS = np.array([lx, 0, 0])
# B_LAMMPS = np.array([bx, by, 0])

A_LAMMPS = A1_cart
B_LAMMPS = A2_cart

xlo, xhi = 0.0, A_LAMMPS[0]
ylo, yhi = 0.0, B_LAMMPS[1]
xy = B_LAMMPS[0]


print(f"Rotation angle (theta): {np.degrees(theta):.3f} degrees")
print(f"Supercel lattice vectors (in lattice units) A1: {A1}, A2: {A2}")

print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
origins, ij_list = moire_primitive_origins(N, M, a1, a2)
print("Number of primitive cells:", len(origins))
print(f"Expected number of primitive cells in supercell (N**2 + N*M + M**2): {Expected_number_prim_cells}")
assert len(origins) == N**2 + N*M + M**2
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")


# Reading first and second layer coordinates
# print(f"Reading coordinates from files {arq_1st_layer} and {arq_2nd_layer}")
# coords_1st_layer, elements_1st_layer = get_prim_cell_coordinates(arq_1st_layer)
# coords_2nd_layer, elements_2nd_layer = get_prim_cell_coordinates(arq_2nd_layer)

# print zmin and zmax for each layer
# zmin_1st = np.min(coords_1st_layer[:, 2])
# zmax_1st = np.max(coords_1st_layer[:, 2])
# zmin_2nd = np.min(coords_2nd_layer[:, 2])
# zmax_2nd = np.max(coords_2nd_layer[:, 2])
# print(f"Layer 1: zmin = {zmin_1st}, zmax = {zmax_1st}")
# print(f"Layer 2: zmin = {zmin_2nd}, zmax = {zmax_2nd}")
# print(f"Original interlayer distance: {zmin_2nd - zmax_1st} Angstroms")

print(f"Total atoms in first layer: {len(coords_1st_layer)}")
print(f"Total atoms in second layer: {len(coords_2nd_layer)}")

# if reinforce_distance_layers:
#     print(f"Reinforcing interlayer distance by {d_interlayer} Angstroms")
#     coords_2nd_layer[:, 2] += (zmax_1st - zmin_2nd + d_interlayer) 


rmax = L_supercell * np.sqrt(3)  # overestimate of radius
max_ind = int(np.ceil(rmax / np.linalg.norm(a1)))

print("max_ind:", max_ind)

# rotating the coordinates of the first and second layers
# coords_1st_twisted_layer, elements_1st_twisted_layer, non_moire_coords_1st, non_moire_elements_1st = rotate_coordinates(coords_1st_layer, theta, elements_1st_layer)
# coords_2nd_twisted_layer, elements_2nd_twisted_layer, non_moire_coords_2nd, non_moire_elements_2nd = rotate_coordinates(coords_2nd_layer, -theta, elements_2nd_layer)

coords_1st_twisted_layer, elements_1st_twisted_layer = \
    rotate_coordinates_integer_moire(
        coords_1st_layer,
        elements_1st_layer,
        origins,
        +theta
    )

coords_2nd_twisted_layer, elements_2nd_twisted_layer = \
    rotate_coordinates_integer_moire(
        coords_2nd_layer,
        elements_2nd_layer,
        origins,
        -theta
    )
    
coords_1st_layer_cart = []
coords_2nd_layer_cart = []

for coord in coords_1st_twisted_layer:
    coords_1st_layer_cart.append(coord * a)
    
for coord in coords_2nd_twisted_layer:
    coords_2nd_layer_cart.append(coord * a)

print(f"Total atoms in moire first layer (before dedup): {len(coords_1st_twisted_layer)}")
print(f"Total atoms in moire second layer (before dedup): {len(coords_2nd_twisted_layer)}")

print(f"Total atoms in moire first layer: {len(coords_1st_twisted_layer)}")
print(f"Total atoms in moire second layer: {len(coords_2nd_twisted_layer)}")

if len(coords_1st_twisted_layer) != Expected_number_atoms_per_layer:
    print("WARNING: Number of atoms in first twisted layer does not match expected number.")
    print(f"Expected: {Expected_number_atoms_per_layer}, Found: {len(coords_1st_twisted_layer)}")
if len(coords_2nd_twisted_layer) != Expected_number_atoms_per_layer:
    print("WARNING: Number of atoms in second twisted layer does not match expected number.")
    print(f"Expected: {Expected_number_atoms_per_layer}, Found: {len(coords_2nd_twisted_layer)}")

Total_atoms = len(coords_1st_twisted_layer) + len(coords_2nd_twisted_layer)

# mapping elements to indices
element_to_index_1st, indexed_elements_1st = get_element_indices(elements_1st_twisted_layer)
element_to_index_2nd, indexed_elements_2nd = get_element_indices(elements_2nd_twisted_layer)

indexed_elements_2nd = list(np.array(indexed_elements_2nd) + max(indexed_elements_1st))

Nelements = max(indexed_elements_2nd)

# rotate everything to make A1 to be alligned along x
angle_to_rotate = -np.arctan2(A1_cart[1], A1_cart[0])

def rotate_z(coord, angle):
    x, y, z = coord
    x_new = x * np.cos(angle) - y * np.sin(angle)
    y_new = x * np.sin(angle) + y * np.cos(angle)
    return np.array([x_new, y_new, z])

coords_1st_layer_cart = [rotate_z(coord, angle_to_rotate) for coord in coords_1st_layer_cart]
coords_2nd_layer_cart = [rotate_z(coord, angle_to_rotate) for coord in coords_2nd_layer_cart]
A1_cart = rotate_z(A1_cart, angle_to_rotate)
A2_cart = rotate_z(A2_cart, angle_to_rotate)

# Check coordinates before wrapping
coords_1st_array = np.array(coords_1st_layer_cart)
coords_2nd_array = np.array(coords_2nd_layer_cart)
print(f"\nBefore wrapping - Layer 1:")
print(f"  x range: [{coords_1st_array[:,0].min():.4f}, {coords_1st_array[:,0].max():.4f}]")
print(f"  y range: [{coords_1st_array[:,1].min():.4f}, {coords_1st_array[:,1].max():.4f}]")
print(f"\nBefore wrapping - Layer 2:")
print(f"  x range: [{coords_2nd_array[:,0].min():.4f}, {coords_2nd_array[:,0].max():.4f}]")
print(f"  y range: [{coords_2nd_array[:,1].min():.4f}, {coords_2nd_array[:,1].max():.4f}]")

# Update LAMMPS vectors after rotation
A_LAMMPS = A1_cart
B_LAMMPS = A2_cart

xlo, ylo = 0.0, 0.0
xhi = A1_cart[0]
xy = A2_cart[0]
yhi = np.sqrt(np.linalg.norm(A2_cart[:2])**2 - xy**2)

# Filter atoms to keep only those inside the parallelogram defined by A1_cart and A2_cart
def is_inside_parallelogram(coord, A1, A2, tol=1e-6):
    """Check if a point is inside the parallelogram defined by A1 and A2"""
    x, y, z = coord
    # Solve for fractional coordinates: (x,y) = s*A1 + t*A2
    # This is a 2x2 linear system
    det = A1[0]*A2[1] - A1[1]*A2[0]
    if abs(det) < 1e-10:
        return False
    
    s = (x*A2[1] - y*A2[0]) / det
    t = (A1[0]*y - A1[1]*x) / det
    
    # Inside if 0 <= s < 1 and 0 <= t < 1
    return (-tol <= s < 1.0 - tol) and (-tol <= t < 1.0 - tol)

print(f"\nFiltering atoms inside parallelogram:")
print(f"  A1_cart: {A1_cart[:2]}")
print(f"  A2_cart: {A2_cart[:2]}")

coords_1st_filtered = []
elements_1st_filtered = []
for coord, elem in zip(coords_1st_layer_cart, elements_1st_twisted_layer):
    if is_inside_parallelogram(coord, A1_cart, A2_cart):
        coords_1st_filtered.append(coord)
        elements_1st_filtered.append(elem)

coords_2nd_filtered = []
elements_2nd_filtered = []
for coord, elem in zip(coords_2nd_layer_cart, elements_2nd_twisted_layer):
    if is_inside_parallelogram(coord, A1_cart, A2_cart):
        coords_2nd_filtered.append(coord)
        elements_2nd_filtered.append(elem)

coords_1st_layer_cart = coords_1st_filtered
coords_2nd_layer_cart = coords_2nd_filtered
elements_1st_twisted_layer = elements_1st_filtered
elements_2nd_twisted_layer = elements_2nd_filtered

print(f"After filtering:")
print(f"  Layer 1: {len(coords_1st_layer_cart)} atoms (expected: {Expected_number_atoms_per_layer})")
print(f"  Layer 2: {len(coords_2nd_layer_cart)} atoms (expected: {Expected_number_atoms_per_layer})")

# Update indexed elements after filtering
element_to_index_1st, indexed_elements_1st = get_element_indices(elements_1st_twisted_layer)
element_to_index_2nd, indexed_elements_2nd = get_element_indices(elements_2nd_twisted_layer)
indexed_elements_2nd = list(np.array(indexed_elements_2nd) + max(indexed_elements_1st))

# Diagnostic: Check for duplicate or very close atoms
def check_min_distances(coords, layer_name, threshold=0.1):
    """Check minimum distance between atoms"""
    coords_array = np.array(coords)
    min_dist = float('inf')
    closest_pair = None
    
    for i in range(len(coords_array)):
        for j in range(i+1, len(coords_array)):
            dist = np.linalg.norm(coords_array[i][:2] - coords_array[j][:2])
            if dist < min_dist:
                min_dist = dist
                closest_pair = (i, j)
    
    print(f"\n{layer_name} - Minimum distance: {min_dist:.4f} Å")
    if min_dist < threshold:
        i, j = closest_pair
        print(f"  WARNING: Atoms {i} and {j} are very close!")
        print(f"  Atom {i}: {coords_array[i][:2]}")
        print(f"  Atom {j}: {coords_array[j][:2]}")
        
        # Check if they're near boundaries
        if coords_array[i][0] < 1.0 or coords_array[i][0] > xhi - 1.0:
            print(f"  Atom {i} is near x boundary")
        if coords_array[i][1] < 1.0 or coords_array[i][1] > yhi - 1.0:
            print(f"  Atom {i} is near y boundary")
        if coords_array[j][0] < 1.0 or coords_array[j][0] > xhi - 1.0:
            print(f"  Atom {j} is near x boundary")
        if coords_array[j][1] < 1.0 or coords_array[j][1] > yhi - 1.0:
            print(f"  Atom {j} is near y boundary")
    
    return min_dist

min_dist_1 = check_min_distances(coords_1st_layer_cart, "Layer 1", threshold=1.5)
min_dist_2 = check_min_distances(coords_2nd_layer_cart, "Layer 2", threshold=1.5)

# Expected nearest neighbor distance in WSe2
expected_nn = a  # This is the lattice parameter
print(f"\nExpected nearest neighbor distance: {expected_nn:.4f} Å")
print(f"Box dimensions: xhi={xhi:.4f}, yhi={yhi:.4f}, xy={xy:.4f}")

import matplotlib.pyplot as plt
plt.figure()

for i, (coord, elem_index) in enumerate(zip(coords_1st_layer_cart, indexed_elements_1st)):
    # x, y, z = coord * a
    xc, yc, zc = coord 
    plt.plot(xc, yc, 'ro', markersize=2, alpha=1)
    
# for i, (coord, elem_index) in enumerate(zip(coords_2nd_layer_cart, indexed_elements_2nd)):
#     # x, y, z = coord * a
#     xc, yc, zc = coord
#     plt.plot(xc, yc, 'bo', markersize=2, alpha=1)

# Draw lattice vectors A1 and A2
plt.arrow(0, 0, A1_cart[0], A1_cart[1], head_width=0.5, head_length=0.5, fc='green', ec='green', linewidth=2, label='A1')
plt.arrow(0, 0, A2_cart[0], A2_cart[1], head_width=0.5, head_length=0.5, fc='orange', ec='orange', linewidth=2, label='A2')
plt.plot([A1_cart[0], A1_cart[0]+A2_cart[0]], [A1_cart[1], A1_cart[1]+A2_cart[1]], 'k--', linewidth=2)
plt.plot([A2_cart[0], A1_cart[0]+A2_cart[0]], [A2_cart[1], A1_cart[1]+A2_cart[1]], 'k--', linewidth=2)

plt.arrow(0, 0, A_LAMMPS[0], A_LAMMPS[1], head_width=0.5, head_length=0.5, fc='red', ec='red', linewidth=2, label='A1 LAMMPS')
plt.arrow(0, 0, B_LAMMPS[0], B_LAMMPS[1], head_width=0.5, head_length=0.5, fc='blue', ec='blue', linewidth=2, label='A2 LAMMPS')
plt.plot([A_LAMMPS[0], A_LAMMPS[0]+B_LAMMPS[0]], [A_LAMMPS[1], A_LAMMPS[1]+B_LAMMPS[1]], 'k--', linewidth=2)
plt.plot([B_LAMMPS[0], A_LAMMPS[0]+B_LAMMPS[0]], [B_LAMMPS[1], A_LAMMPS[1]+B_LAMMPS[1]], 'k--', linewidth=2)



for IC in [-1, 0, 1]:
    for JC in [-1, 0, 1]:
        if IC == 0 and JC == 0:
            continue
        print(f"Translating by IC={IC}, JC={JC}")
        for i, (coord, elem_index) in enumerate(zip(coords_1st_layer_cart, indexed_elements_1st)):
            # x, y, z = coord * a
            xc, yc, zc = coord 
            xc += A_LAMMPS[0] * IC + B_LAMMPS[0] * JC
            yc += A_LAMMPS[1] * IC + B_LAMMPS[1] * JC
            plt.plot(xc, yc, 'k*', markersize=2, alpha=1)
            
        # for i, (coord, elem_index) in enumerate(zip(coords_2nd_layer_cart, indexed_elements_2nd)):
        #     # x, y, z = coord * a
        #     xc, yc, zc = coord
        #     xc += A_LAMMPS[0] * IC + B_LAMMPS[0] * JC
        #     yc += A_LAMMPS[1] * IC + B_LAMMPS[1] * JC
        #     plt.plot(xc, yc, 'b*', markersize=2, alpha=0.2)


plt.axis('equal')
plt.legend()
plt.savefig('1st_layer_twisted.png', dpi=300)
# plt.show()




print("Writing coordinates to lammps files...")

topo_lammps = open(str(N)+'_'+str(M)+'_twisted.top', 'w')
xyz_file = open(str(N)+'_'+str(M)+'_twisted.xyz', 'w')
xyz_file.write(f"{Total_atoms}\n\n")

topo_lammps.write('\n\n# indices '+str(N)+' '+str(M)+'\n')
topo_lammps.write('# '+str(L_supercell * a)+' angstrons - supercell size\n')
topo_lammps.write('# '+str((180.0/np.pi)*theta)+' degrees - rotation angle\n')
topo_lammps.write(f'\n\n{Total_atoms} atoms\n')
topo_lammps.write(f'{Nelements}  atom types \n\n')

topo_lammps.write(f"{xlo:.8f} {xhi:.8f} xlo xhi\n")
topo_lammps.write(f"{ylo:.8f} {yhi:.8f} ylo yhi\n")
topo_lammps.write("-50.0 50.0 zlo zhi\n")
topo_lammps.write(f"{xy:.8f} 0.0 0.0 xy xz yz\n\n")
topo_lammps.write('Atoms \n\n')

# 1st layer atoms
atom_index = 0
mol_index = 1  
for i, (coord, elem_index) in enumerate(zip(coords_1st_layer_cart, indexed_elements_1st)):
    # x, y, z = coord * a
    x, y, z = coord
    atom_index += 1
    # atom_index molecule_index atom_type x y z
    # topo_lammps.write(f"{i+1} 1 {elem_index} {x} {y} {z}\n")
    topo_lammps.write(f"{atom_index:6d} {mol_index:4d} {elem_index:4d} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    xyz_file.write(f"{elements_1st_twisted_layer[i]} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    
# 2nd layer atoms
mol_index = 2
for i, (coord, elem_index) in enumerate(zip(coords_2nd_layer_cart, indexed_elements_2nd)):
    # x, y, z = coord * a
    x, y, z = coord
    atom_index += 1
    # atom_index molecule_index atom_type x y z
    # topo_lammps.write(f"{i+1 + len(coords_1st_twisted_layer)} 2 {elem_index} {x} {y} {z}\n")
    topo_lammps.write(f"{atom_index:6d} {mol_index:4d} {elem_index:4d} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    xyz_file.write(f"{elements_2nd_twisted_layer[i]} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    
    
xyz_file.close()
topo_lammps.close()

print(f"Coordinates written to {str(N)}_{str(M)}_twisted.top")