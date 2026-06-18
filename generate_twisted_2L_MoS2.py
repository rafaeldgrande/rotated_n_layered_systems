
# Generates twisted bilayer MoS2 (AA registry at theta=0) LAMMPS data files.
# AA is pinned at the rotation pivot (origin), so the AA high-symmetry point
# sits at (0,0) for every twist angle.
# Adapted from generate_twisted_2L_WSe2_new.py.

from __future__ import division
import numpy as np
import sys

def dentro_supercelula(x, y, L_supercell):  # Verifica se o vetor posicao ta dentro da supercelula
    # (x ,y) = alfa*L1 + beta*L2, onde L1 = L(1, 0) e L2 = L(0.5, sqrt(3)/2)
    alfa, beta = (x-y/np.sqrt(3))/L_supercell, 2*y/(L_supercell*np.sqrt(3))
    if 0 <= alfa < 1 - delta and 0 <= beta < 1 - delta:
        return True
    else:
        return False

N, M = int(sys.argv[1]), int(sys.argv[2])  # Indices supercelula from command line

arq_out = open('t2LMoS2_'+str(N)+'_'+str(M)+'.xyz', 'w')

topo_lammps = open(f"{N}_{M}_twisted.top", 'w')

Arq_out = []
ATOMS_L1, ATOMS_L2 = [], []
ATOMS_L1_TYPES = []
ATOMS_L2_TYPES = []

# Equilibrium values from SW+KC vc-relax of untwisted bilayer MoS2, AB_3R stacking
# (see /Users/rdelgrande/work/TEMP_DATA/moire_classical_FF/MoS2-MoS2/AB_3R/bilayer_MoS2_AB_3R_relaxed.dat)
a = 3.11702       # Parametro de rede (relaxed, AB_3R)
d_MoMo = 6.08719  # Distancia interplanar Mo-Mo (relaxed, AB_3R)
d_Mo_S = 1.55883  # Mo-S out-of-plane offset (relaxed, AB_3R)

elements = ['Mo', 'S', 'S']

teste_ind = 3*max(N, M)

theta = np.arccos((N**2 + 4*N*M + M**2)/(2*(N**2 + N*M + M**2)))
L_supercell = a*np.sqrt(N**2 + N*M + M**2)
L = L_supercell/abs(M-N)  # Periodicidade do padrao de Moire

delta = 1e-3*a/L_supercell

# Celula lammps
xhi = L_supercell
yhi = L_supercell*np.sin(np.pi/3)
xy = L_supercell*np.cos(np.pi/3)


# Base - celula MoS2
a1 = a*np.array([np.sqrt(3)/2, 0.5, 0])
a2 = a*np.array([np.sqrt(3)/2, -0.5, 0])

# Base - supercelula
A1 = L_supercell*np.array([1, 0, 0])
A2 = L_supercell*np.array([0.5, np.sqrt(3)/2.0, 0])

# Primeira camada (registro de referencia)
tau1 = a*np.array([0, 0, 0])
tau2 = a*np.array([1/np.sqrt(3), 0, d_Mo_S/a])
tau3 = a*np.array([1/np.sqrt(3), 0, -d_Mo_S/a])

TAUS_L1 = [tau1, tau2, tau3]

# Segunda camada - registro AA: mesma base em-plano da camada 1 (Mo2 sobre Mo1),
# deslocada apenas em z. O ponto de rede i=j=0 fica exatamente no pivo de
# rotacao (origem), entao o registro AA permanece fixo na origem para
# qualquer angulo de torcao theta.
tau4 = a*np.array([0, 0, d_MoMo/a])
tau5 = a*np.array([1/np.sqrt(3), 0, (d_MoMo + d_Mo_S)/a])
tau6 = a*np.array([1/np.sqrt(3), 0, (d_MoMo - d_Mo_S)/a])

TAUS_L2 = [tau4, tau5, tau6]

for i in range(-teste_ind, teste_ind):
    for j in range(-teste_ind, teste_ind):
        if i % 10 == 0:  # Print progress every 10 iterations
            print(f"Processing i: {i}/{teste_ind}")

        for itau, tau in enumerate(TAUS_L1):

            x , y, z = i*a1 + j*a2 + tau

            # Rotacionando
            x_rot = x*np.cos(theta/2.0) + y*np.sin(theta/2.0)
            y_rot = -x*np.sin(theta/2.0) + y*np.cos(theta/2.0)

            if dentro_supercelula(x_rot, y_rot, L_supercell) == True:

                element = str(itau + 1)  # 1: Mo, 2: S, 3: S

                linha = f'{element}   {x_rot}   {y_rot}   {z}\n'
                Arq_out.append(linha)
                ATOMS_L1.append(linha)
                ATOMS_L1_TYPES.append(element)

# Segunda camada

for i in range(-teste_ind, teste_ind):
    for j in range(-teste_ind, teste_ind):

        for itau, tau in enumerate(TAUS_L2):

            x , y, z = i*a1 + j*a2 + tau

            # Rotacionando
            x_rot = x*np.cos(-theta/2.0) + y*np.sin(-theta/2.0)
            y_rot = -x*np.sin(-theta/2.0) + y*np.cos(-theta/2.0)

            if dentro_supercelula(x_rot, y_rot, L_supercell) == True:

                element = str(itau + 4)  # 4: Mo, 5: S, 6: S

                linha = f'{element}   {x_rot}   {y_rot}   {z}\n'
                Arq_out.append(linha)
                ATOMS_L2.append(linha)
                ATOMS_L2_TYPES.append(element)

# Arquivo de saida

Total_atomos = len(ATOMS_L1) + len(ATOMS_L2)

arq_out.write(str(Total_atomos)+'\n\n')

for i in range(len(ATOMS_L1)):
    arq_out.write(ATOMS_L1[i])
for i in range(len(ATOMS_L2)):
    arq_out.write(ATOMS_L2[i])

arq_out.close()

# Celula unitaria

arq = open('cell_basis', 'w')

arq.write(str(L_supercell)+' 0.0   0.0 \n')
arq.write(str(L_supercell/2.0) + '       '+ str(L_supercell*np.sqrt(3)/2) +'0.0 \n')
arq.write('0 0 20.0')

arq.close()

# Arquivo de topologia lammps

topo_lammps.write('\n\n# indices '+str(N)+' '+str(M)+'\n')
topo_lammps.write('# '+str(L_supercell)+' angstrons - supercell size\n')
topo_lammps.write('# '+str((180.0/np.pi)*theta)+' degrees - rotation angle\n')
topo_lammps.write(f'\n\n  {Total_atomos} atoms\n')
topo_lammps.write('6  atom types \n\n')

topo_lammps.write('0.0 '+str(xhi)+' xlo xhi \n')
topo_lammps.write('0.0 '+str(yhi)+' ylo yhi \n')
topo_lammps.write('-50.0 50.0 zlo zhi \n')
topo_lammps.write(str(xy)+' 0.0 0.0 xy xz yz \n\n')
topo_lammps.write('Atoms \n\n')

ind_molecula = 1
for i in range(len(ATOMS_L1)):
    linha = ATOMS_L1[i].split()
    atom_type = ATOMS_L1_TYPES[i]
    x, y, z = linha[1], linha[2], linha[3]
    # atom_index mol_index atom_type x y z
    topo_lammps.write(f"{i+1} {ind_molecula} {atom_type} {x} {y} {z}\n")

ind_molecula = 2
for i in range(len(ATOMS_L2)):
    linha = ATOMS_L2[i].split()
    atom_type = ATOMS_L2_TYPES[i]
    x, y, z = linha[1], linha[2], linha[3]
    # atom_index mol_index atom_type x y z
    topo_lammps.write(f"{i+1+len(ATOMS_L1_TYPES)} {ind_molecula} {atom_type} {x} {y} {z}\n")

topo_lammps.close()

print('done')
