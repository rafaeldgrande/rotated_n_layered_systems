
# Generates twisted bilayer MoSe2-WSe2 (AA registry at theta=0) LAMMPS data files.
# AA is pinned at the rotation pivot (origin): Mo1 of MoSe2 (L1) above W2 of WSe2 (L2).
# L1 rotated by +theta/2, L2 by -theta/2.
# Adapted from generate_twisted_2L_MoS2.py.

from __future__ import division
import numpy as np
import sys
import os

def dentro_supercelula(x, y, L_supercell):
    alfa, beta = (x - y/np.sqrt(3))/L_supercell, 2*y/(L_supercell*np.sqrt(3))
    if 0 <= alfa < 1 - delta and 0 <= beta < 1 - delta:
        return True
    return False

N, M = int(sys.argv[1]), int(sys.argv[2])

# Equilibrium parameters from AB_3R bilayer relaxation (sw/mod + KC)
a      = 3.2977    # in-plane lattice constant (Å)
d_MoW  = 6.4598   # Mo-W interlayer distance (Å)
d_Mo_Se = 1.6558  # Mo-Se vertical offset in MoSe2 L1 (Å)
d_W_Se  = 1.6443  # W-Se vertical offset in WSe2 L2 (Å)

# Atom types
# L1 (MoSe2): 1=Mo1, 2=Se2(top), 3=Se1(bot)
# L2 (WSe2):  4=W2,  5=Se4(top), 6=Se3(bot)

outdir = f'{N}_{M}'
os.makedirs(outdir, exist_ok=True)

arq_out = open(f'{outdir}/t2L_MoSe2WSe2_{N}_{M}.xyz', 'w')
topo_lammps = open(f'{outdir}/twisted_MoSe2WSe2_initial.top', 'w')

ATOMS_L1, ATOMS_L2 = [], []
ATOMS_L1_TYPES, ATOMS_L2_TYPES = [], []

theta = np.arccos((N**2 + 4*N*M + M**2) / (2*(N**2 + N*M + M**2)))
L_supercell = a * np.sqrt(N**2 + N*M + M**2)

delta = 1e-3 * a / L_supercell

xhi = L_supercell
yhi = L_supercell * np.sin(np.pi/3)
xy  = L_supercell * np.cos(np.pi/3)

teste_ind = 3 * max(N, M)

# Primitive lattice vectors (same orientation as MoS2 script)
a1 = a * np.array([np.sqrt(3)/2,  0.5, 0])
a2 = a * np.array([np.sqrt(3)/2, -0.5, 0])

# Supercell lattice vectors
A1 = L_supercell * np.array([1, 0, 0])
A2 = L_supercell * np.array([0.5, np.sqrt(3)/2, 0])

# Basis for L1 (MoSe2): AA pinned at origin → Mo at (0,0,0)
tau1 = a * np.array([0,             0, 0])
tau2 = a * np.array([1/np.sqrt(3),  0,  d_Mo_Se/a])
tau3 = a * np.array([1/np.sqrt(3),  0, -d_Mo_Se/a])

TAUS_L1 = [tau1, tau2, tau3]

# Basis for L2 (WSe2): AA registry at origin → W at (0,0,d_MoW)
tau4 = a * np.array([0,             0,  d_MoW/a])
tau5 = a * np.array([1/np.sqrt(3),  0, (d_MoW + d_W_Se)/a])
tau6 = a * np.array([1/np.sqrt(3),  0, (d_MoW - d_W_Se)/a])

TAUS_L2 = [tau4, tau5, tau6]

# Layer 1: rotate by +theta/2
for i in range(-teste_ind, teste_ind):
    for j in range(-teste_ind, teste_ind):
        for itau, tau in enumerate(TAUS_L1):
            x, y, z = i*a1 + j*a2 + tau
            x_rot =  x*np.cos(theta/2) + y*np.sin(theta/2)
            y_rot = -x*np.sin(theta/2) + y*np.cos(theta/2)
            if dentro_supercelula(x_rot, y_rot, L_supercell):
                element = str(itau + 1)   # 1=Mo1, 2=Se2, 3=Se1
                linha = f'{element}   {x_rot}   {y_rot}   {z}\n'
                ATOMS_L1.append(linha)
                ATOMS_L1_TYPES.append(element)

# Layer 2: rotate by -theta/2
for i in range(-teste_ind, teste_ind):
    for j in range(-teste_ind, teste_ind):
        for itau, tau in enumerate(TAUS_L2):
            x, y, z = i*a1 + j*a2 + tau
            x_rot =  x*np.cos(-theta/2) + y*np.sin(-theta/2)
            y_rot = -x*np.sin(-theta/2) + y*np.cos(-theta/2)
            if dentro_supercelula(x_rot, y_rot, L_supercell):
                element = str(itau + 4)   # 4=W2, 5=Se4, 6=Se3
                linha = f'{element}   {x_rot}   {y_rot}   {z}\n'
                ATOMS_L2.append(linha)
                ATOMS_L2_TYPES.append(element)

Total_atomos = len(ATOMS_L1) + len(ATOMS_L2)

# XYZ output
arq_out.write(f'{Total_atomos}\n\n')
for l in ATOMS_L1:
    arq_out.write(l)
for l in ATOMS_L2:
    arq_out.write(l)
arq_out.close()

# LAMMPS topology file
topo_lammps.write(f'\n\n# indices {N} {M}\n')
topo_lammps.write(f'# {L_supercell} angstrons - supercell size\n')
topo_lammps.write(f'# {(180.0/np.pi)*theta} degrees - rotation angle\n')
topo_lammps.write(f'\n\n  {Total_atomos} atoms\n')
topo_lammps.write('6  atom types \n\n')

topo_lammps.write(f'0.0 {xhi} xlo xhi \n')
topo_lammps.write(f'0.0 {yhi} ylo yhi \n')
topo_lammps.write('-50.0 50.0 zlo zhi \n')
topo_lammps.write(f'{xy} 0.0 0.0 xy xz yz \n\n')
topo_lammps.write('Atoms \n\n')

ind_molecula = 1
for i, linha in enumerate(ATOMS_L1):
    parts = linha.split()
    atype, x, y, z = ATOMS_L1_TYPES[i], parts[1], parts[2], parts[3]
    topo_lammps.write(f'{i+1} {ind_molecula} {atype} {x} {y} {z}\n')

ind_molecula = 2
for i, linha in enumerate(ATOMS_L2):
    parts = linha.split()
    atype, x, y, z = ATOMS_L2_TYPES[i], parts[1], parts[2], parts[3]
    topo_lammps.write(f'{i+1+len(ATOMS_L1)} {ind_molecula} {atype} {x} {y} {z}\n')

topo_lammps.close()

print(f'N={N} M={M}  theta={np.degrees(theta):.4f} deg  L={L_supercell:.4f} A  atoms={Total_atomos}')
