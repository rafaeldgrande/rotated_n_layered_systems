

from __future__ import division
import numpy as np
import argparse
import sys

def rotate2d(v, theta_rad):
    """Rotate a 2D vector or (N,2) array by theta_rad around the origin."""
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    R = np.array([[c, -s],
                  [s,  c]])
    return R @ v
def expected_atom_counts(sc, natoms_per_cell_L1, natoms_per_cell_L2):
    """
    Compute expected number of primitive cells and atoms in the supercell.
    Uses the integer index determinant: N = |det([[n1,p1],[n2,p2]])|
    """
    N_L1 = abs(sc['n1']*sc['p2'] - sc['n2']*sc['p1'])
    N_L2 = abs(sc['m1']*sc['q2'] - sc['m2']*sc['q1'])

    print(f"Layer 1: {N_L1} primitive cells × {natoms_per_cell_L1} atoms = {N_L1*natoms_per_cell_L1} atoms")
    print(f"Layer 2: {N_L2} primitive cells × {natoms_per_cell_L2} atoms = {N_L2*natoms_per_cell_L2} atoms")
    print(f"Total:   {N_L1*natoms_per_cell_L1 + N_L2*natoms_per_cell_L2} atoms")

    return N_L1 * natoms_per_cell_L1, N_L2 * natoms_per_cell_L2

def dentro_supercelula(r_pos, V1_lammps, V2_lammps, delta=1e-6):
    """Check if 2D position r_pos is inside the supercell
    spanned by V1_lammps and V2_lammps."""
    M = np.array([[V1_lammps[0], V2_lammps[0]],
                  [V1_lammps[1], V2_lammps[1]]])
    alpha, beta = np.linalg.solve(M, r_pos)
    return (-delta <= alpha < 1 - delta) and (-delta <= beta < 1 - delta)

def hex_vecs(a, theta_rad=0.0):
    """2D primitive vectors for a hexagonal lattice rotated by theta_rad.

    Convention:  e1 = a*(cos t, sin t)
                 e2 = a*(cos(t+60), sin(t+60))
    giving       e1 = a*(1, 0)  and  e2 = a*(1/2, sqrt(3)/2)  at t=0.
    """
    c,   s   = np.cos(theta_rad),             np.sin(theta_rad)
    c60, s60 = np.cos(theta_rad + np.pi/3),   np.sin(theta_rad + np.pi/3)
    return a * np.array([c, s]), a * np.array([c60, s60])

parser = argparse.ArgumentParser(
    description='Build a twisted TMD heterobilayer supercell for LAMMPS.')
parser.add_argument('--a1',    type=float, default=3.31,
                    help='Lattice constant layer 1 (A) [WSe2 default]')
parser.add_argument('--a2',    type=float, default=3.11,
                    help='Lattice constant layer 2 (A) [MoS2 default]')
parser.add_argument('--theta', type=float, default=1.0,
                    help='Twist angle (degrees)')
parser.add_argument('--h1',    type=float, default=1.6,
                    help='Chalcogen-metal height layer 1 (A)')
parser.add_argument('--h2',    type=float, default=1.5,
                    help='Chalcogen-metal height layer 2 (A)')
parser.add_argument('--d_int', type=float, default=7.5,
                    help='Metal-metal interlayer distance (A)')
parser.add_argument('--tol',   type=float, default=0.05,
                    help='Commensurability tolerance (A)')
parser.add_argument('--N_max', type=int,   default=100,
                    help='Max integer search range for supercell indices')
parser.add_argument('--output', type=str,  default=None,
                    help='Output .dat filename (default: auto-generated)')
parser.add_argument('--m1_metal',  type=float, default=0.0,
                    help='Atomic mass of L1 metal (0 = omit Masses section)')
parser.add_argument('--m1_chalc',  type=float, default=0.0,
                    help='Atomic mass of L1 chalcogen')
parser.add_argument('--m2_metal',  type=float, default=0.0,
                    help='Atomic mass of L2 metal')
parser.add_argument('--m2_chalc',  type=float, default=0.0,
                    help='Atomic mass of L2 chalcogen')
args = parser.parse_args()

a1     = args.a1
a2     = args.a2
theta  = args.theta
d_MCalc1 = args.h1
d_MCalc2 = args.h2
d_MM   = args.d_int

elements_L1 = ['W', 'Se', 'Se']
elements_L1_indexes = [1, 2, 3]
elements_L2 = ['Mo', 'S', 'S']
elements_L2_indexes = [4, 5, 6]

# in plane lattice vectors L1
a1L1 = a1*np.array([np.sqrt(3)/2, 0.5, 0])
a2L1 = a1*np.array([np.sqrt(3)/2, -0.5, 0])

# in plane lattice vectors L2
a1L2 = a2*np.array([np.sqrt(3)/2, 0.5, 0])
a2L2 = a2*np.array([np.sqrt(3)/2, -0.5, 0])

# # motifs layer 1
# tau1 = a1*np.array([0, 0, 0])
# tau2 = a1*np.array([1/np.sqrt(3), 0, d_MCalc1/a1])
# tau3 = a1*np.array([1/np.sqrt(3), 0, -d_MCalc1/a1])
# TAUS_L1 = [tau1, tau2, tau3]

# # motifs layer 2
# tau4 = a2*np.array([0, 0, d_MCalc2/a2])
# tau5 = a2*np.array([1/np.sqrt(3), 0, (d_MM + d_MCalc2)/a2])
# tau6 = a2*np.array([1/np.sqrt(3), 0, (d_MM - d_MCalc2)/a2])
# TAUS_L2 = [tau4, tau5, tau6]



# finding supercell lattice vectors
def moire_period(a1, a2, theta_rad):
    denom = np.sqrt(a1**2 + a2**2 - 2.0*a1*a2*np.cos(theta_rad))
    return a1 * a2 / denom if denom > 1e-12 else np.inf

def find_supercell(a1, a2, theta_deg, tol=0.05, N_max=60):
    """
    Find two shortest linearly independent commensurate supercell vectors.

    For each candidate layer-1 vector L = n1*e1 + n2*e2, find the nearest
    layer-2 integer combination m1*f1+m2*f2; accept if residual < tol.

    Returns dict with keys:
        L1, L2        – Cartesian 2D supercell vectors (layer-1 unrotated frame)
        n1,n2,m1,m2   – layer-1 / layer-2 integer indices for L1
        p1,p2,q1,q2   – layer-1 / layer-2 integer indices for L2
                        (already corrected for any orientation flip)
        res1, res2    – commensurability residuals (Angstrom)
    """
    theta_rad = np.radians(theta_deg)
    e1, e2 = hex_vecs(a1, 0.0)
    f1, f2 = hex_vecs(a2, theta_rad)

    F    = np.array([[f1[0], f2[0]], [f1[1], f2[1]]])
    Finv = np.linalg.inv(F)

    accepted = []
    a_min = 0.3 * min(a1, a2)

    for n1 in range(-N_max, N_max + 1):
        for n2 in range(-N_max, N_max + 1):
            if n1 == 0 and n2 == 0:
                continue
            L   = n1*e1 + n2*e2
            mag = np.linalg.norm(L)
            if mag < a_min:
                continue
            mf  = Finv @ L
            m1r = int(np.round(mf[0]))
            m2r = int(np.round(mf[1]))
            res = np.linalg.norm(L - (m1r*f1 + m2r*f2))
            if res < tol:
                accepted.append({'L': L.copy(), 'n1': n1, 'n2': n2,
                                 'm1': m1r, 'm2': m2r, 'res': res, 'mag': mag})

    if not accepted:
        print(f'ERROR: no commensurate vector found (tol={tol} A, N_max={N_max}).')
        print('Try increasing --tol or --N_max.')
        sys.exit(1)

    accepted.sort(key=lambda d: d['mag'])
    c1 = accepted[0]
    L1 = c1['L']

    c2 = None
    for c in accepted[1:]:
        cross = abs(L1[0]*c['L'][1] - L1[1]*c['L'][0])
        if cross > 1.0:
            if c2 is None or c['mag'] < c2['mag']:
                c2 = c

    if c2 is None:
        print('ERROR: could not find two linearly independent supercell vectors.')
        sys.exit(1)

    L2 = c2['L']
    # Ensure positive orientation (L1 x L2 > 0) so yhi > 0 for LAMMPS
    flipped = (L1[0]*L2[1] - L1[1]*L2[0]) < 0
    if flipped:
        L2 = -L2
    sign = -1 if flipped else 1   # propagate sign to integer indices

    return {
        'L1': L1, 'L2': L2,
        'n1': c1['n1'],       'n2': c1['n2'],
        'm1': c1['m1'],       'm2': c1['m2'],       'res1': c1['res'],
        'p1': sign*c2['n1'],  'p2': sign*c2['n2'],
        'q1': sign*c2['m1'],  'q2': sign*c2['m2'],  'res2': c2['res'],
    }
     
    
sc = find_supercell(a1, a2, theta, tol=args.tol, N_max=args.N_max)
theta_rot = -np.arctan2(sc['L1'][1], sc['L1'][0])

V1_lammps = rotate2d(sc['L1'], theta_rot)
V2_lammps = rotate2d(sc['L2'], theta_rot)

# Lattice vectors: hex_vecs frame → LAMMPS frame
_e1_L1, _e2_L1 = hex_vecs(a1, 0.0)
_e1_L2, _e2_L2 = hex_vecs(a2, np.radians(theta))

a1L1 = np.append(rotate2d(_e1_L1, theta_rot), 0.0)
a2L1 = np.append(rotate2d(_e2_L1, theta_rot), 0.0)
a1L2 = np.append(rotate2d(_e1_L2, theta_rot), 0.0)
a2L2 = np.append(rotate2d(_e2_L2, theta_rot), 0.0)

# Chalcogen 2D position: fractional (1/3, 1/3) in hex_vecs = hollow site of M triangle
_tau_X_L1 = rotate2d((1/3)*_e1_L1 + (1/3)*_e2_L1, theta_rot)
_tau_X_L2 = rotate2d((1/3)*_e1_L2 + (1/3)*_e2_L2, theta_rot)

tau1 = np.array([0.0,           0.0,           0.0      ])
tau2 = np.array([_tau_X_L1[0], _tau_X_L1[1],  d_MCalc1 ])
tau3 = np.array([_tau_X_L1[0], _tau_X_L1[1], -d_MCalc1 ])
TAUS_L1 = [tau1, tau2, tau3]

tau4 = np.array([0.0,           0.0,           d_MM            ])  # was d_MCalc2 ← bug
tau5 = np.array([_tau_X_L2[0], _tau_X_L2[1],  d_MM + d_MCalc2 ])
tau6 = np.array([_tau_X_L2[0], _tau_X_L2[1],  d_MM - d_MCalc2 ])
TAUS_L2 = [tau4, tau5, tau6]

delta = 1e-3*(a1+a2)/(2*np.linalg.norm(V1_lammps))

N1 = round(np.linalg.norm(V1_lammps)/np.linalg.norm(a1))
N2 = round(np.linalg.norm(V1_lammps)/np.linalg.norm(a2))

Nrep = 3 * max(N1, N2)

expected_L1, expected_L2 = expected_atom_counts(sc, len(TAUS_L1), len(TAUS_L2))

# gen_pos: everything is already in LAMMPS frame, no extra rotation needed
def gen_pos(Nrep, TAUS_LIST, V1_lammps, V2_lammps, elements_indexes, a1L, a2L):
    atom_positions = []
    for i in range(-Nrep, Nrep + 1):
        if i % 10 == 0:
            print(f"Processing: {i+Nrep}/{2*Nrep+1} ({100*(i+Nrep)/(2*Nrep+1):.2f} %) ")
        for j in range(-Nrep, Nrep + 1):
            r_origin = i*a1L + j*a2L
            if dentro_supercelula(r_origin[:2], V1_lammps, V2_lammps):
                for itau, tau in enumerate(TAUS_LIST):
                    r_pos = r_origin + tau
                    element_index = elements_indexes[itau]
                    linha = f'{element_index}   {r_pos[0]:.9f}   {r_pos[1]:.9f}   {r_pos[2]:.9f}\n'
                    atom_positions.append(linha)
    return atom_positions

def remove_boundary_duplicates(atom_list, V1, V2, n_basis=3,
                                boundary_tol=0.01, pos_tol=0.1):
    """
    Remove unit cells that were double-counted at the supercell boundary.

    When the commensurability residual is nonzero, some origins that should
    sit exactly at fractional coord 1 land just below it and pass the
    containment test twice.  For each origin near a boundary, this function
    applies +-V1/V2 translations and checks whether the shifted position
    coincides with another origin already in the list.  When a duplicate
    pair is found, the one whose fractional coordinate is closest to 1
    (the spurious extra copy) is removed.
    """
    V1 = np.array(V1[:2])
    V2 = np.array(V2[:2])
    M    = np.array([[V1[0], V2[0]], [V1[1], V2[1]]])
    Minv = np.linalg.inv(M)

    n_cells = len(atom_list) // n_basis
    origins = np.array([
        [float(atom_list[k * n_basis].split()[1]),
         float(atom_list[k * n_basis].split()[2])]
        for k in range(n_cells)
    ])
    fracs = (Minv @ origins.T).T          # shape (n_cells, 2)

    translations = [n1 * V1 + n2 * V2
                    for n1 in (-1, 0, 1)
                    for n2 in (-1, 0, 1)
                    if not (n1 == 0 and n2 == 0)]

    remove = set()
    for i in range(n_cells):
        if i in remove:
            continue
        fi = fracs[i]
        if not any(abs(fi[d]) < boundary_tol or abs(fi[d] - 1.0) < boundary_tol
                   for d in range(2)):
            continue
        for T in translations:
            dists = np.linalg.norm(origins - (origins[i] + T), axis=1)
            for j in np.where(dists < pos_tol)[0]:
                j = int(j)
                if j == i or j in remove:
                    continue
                fj = fracs[j]
                # Remove whichever has a fractional coord closest to 1
                score_i = max(fi[0] % 1.0, fi[1] % 1.0)
                score_j = max(fj[0] % 1.0, fj[1] % 1.0)
                remove.add(i if score_i > score_j else j)
                break

    result = []
    for k in range(n_cells):
        if k not in remove:
            result.extend(atom_list[k * n_basis:(k + 1) * n_basis])
    print(f"  remove_boundary_duplicates: removed {len(remove)} duplicate cell(s)")
    return result


atoms_L1 = gen_pos(Nrep, TAUS_L1, V1_lammps, V2_lammps, elements_L1_indexes, a1L1, a2L1)
atoms_L2 = gen_pos(Nrep, TAUS_L2, V1_lammps, V2_lammps, elements_L2_indexes, a1L2, a2L2)
atoms_L2 = remove_boundary_duplicates(atoms_L2, V1_lammps, V2_lammps, n_basis=len(TAUS_L2))

assert len(atoms_L1) == expected_L1, f"L1: expected {expected_L1}, got {len(atoms_L1)}"
assert len(atoms_L2) == expected_L2, f"L2: expected {expected_L2}, got {len(atoms_L2)}"

print(f'Got {len(atoms_L1)} atoms in L1 and {len(atoms_L2)} atoms in L2')


# LAMMPS data file
# V1_lammps is aligned with x after theta_rot, so V1[1]~0; V2 gives xy tilt and yhi
Total_atomos = len(atoms_L1) + len(atoms_L2)
xhi = V1_lammps[0]
yhi = V2_lammps[1]
xy  = V2_lammps[0]

out_file = args.output if args.output else f"WSe2_MoS2_theta_{theta:.2f}_twisted.dat"
with open(out_file, 'w') as topo_lammps:
    topo_lammps.write(f'# WSe2/MoS2 heterobilayer, theta={theta:.4f} deg\n')
    topo_lammps.write(f'# supercell size: {np.linalg.norm(V1_lammps):.6f} angstroms\n')
    topo_lammps.write(f'\n  {Total_atomos} atoms\n')
    topo_lammps.write('6  atom types \n\n')

    topo_lammps.write(f'0.0 {xhi:.9f} xlo xhi \n')
    topo_lammps.write(f'0.0 {yhi:.9f} ylo yhi \n')
    topo_lammps.write('-50.0 50.0 zlo zhi \n')
    topo_lammps.write(f'{xy:.9f} 0.0 0.0 xy xz yz \n\n')
    if args.m1_metal > 0:
        topo_lammps.write('Masses\n\n')
        topo_lammps.write(f'1 {args.m1_metal}\n')
        topo_lammps.write(f'2 {args.m1_chalc}\n')
        topo_lammps.write(f'3 {args.m1_chalc}\n')
        topo_lammps.write(f'4 {args.m2_metal}\n')
        topo_lammps.write(f'5 {args.m2_chalc}\n')
        topo_lammps.write(f'6 {args.m2_chalc}\n')
        topo_lammps.write('\n')
    topo_lammps.write('Atoms \n\n')

    ind_molecula = 1
    for i, linha_str in enumerate(atoms_L1):
        fields = linha_str.split()
        atom_type, x, y, z = fields[0], fields[1], fields[2], fields[3]
        topo_lammps.write(f"{i+1} {ind_molecula} {atom_type} {x} {y} {z}\n")

    ind_molecula = 2
    offset = len(atoms_L1)
    for i, linha_str in enumerate(atoms_L2):
        fields = linha_str.split()
        atom_type, x, y, z = fields[0], fields[1], fields[2], fields[3]
        topo_lammps.write(f"{i+1+offset} {ind_molecula} {atom_type} {x} {y} {z}\n")