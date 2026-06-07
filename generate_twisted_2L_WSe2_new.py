
# Fazendo uma gambiarra para fazer 4LG

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

# print(sys.argv[1])
# print(type(sys.argv[1]))
N, M = int(sys.argv[1]), int(sys.argv[2])  # Indices supercelula from command line
# dir N_M
# N, M = int(sys.argv[1].split('_')[0]), int(sys.argv[1].split('_')[1])

arq_visualizacao = open('visual_t6LG_'+str(N)+'_'+str(M)+'.xyz', 'w')
arq_out = open('t6LG_'+str(N)+'_'+str(M)+'.xyz', 'w')

# topo_lammps = open('t6LG_'+str(N)+'_'+str(M)+'.top', 'w')
topo_lammps = open(f"{N}_{M}_twisted.top", 'w')

Arq_out = []
ATOMS_L1, ATOMS_L2 = [], []
ATOMS_L1_TYPES = []
ATOMS_L2_TYPES = []

a = 3.29  # Parametro de rede
d_WW = 7.5 #3.5  # Distancia interplanar
d_W_Se = 1.66

elements = ['W', 'Se', 'Se']


teste_ind = 3*max(N, M)

theta = np.arccos((N**2 + 4*N*M + M**2)/(2*(N**2 + N*M + M**2)))
L_supercell = a*np.sqrt(N**2 + N*M + M**2)
L = L_supercell/abs(M-N)  # Periodicidade do padrao de Moire

delta = 1e-3*a/L_supercell

# Celula lammps
xhi = L_supercell
yhi = L_supercell*np.sin(np.pi/3)
xy = L_supercell*np.cos(np.pi/3)


# Base - celula grafeno
a1 = a*np.array([np.sqrt(3)/2, 0.5, 0])
a2 = a*np.array([np.sqrt(3)/2, -0.5, 0])

# Base - supercelula

A1 = L_supercell*np.array([1, 0, 0])
A2 = L_supercell*np.array([0.5, np.sqrt(3)/2.0, 0])

# Primeira camada

# Motivos - stacking ABA
tau1 = a*np.array([0, 0, 0])
tau2 = a*np.array([1/np.sqrt(3), 0, d_W_Se/a])
tau3 = a*np.array([1/np.sqrt(3), 0, -d_W_Se/a])

TAUS_L1 = [tau1, tau2, tau3]

# Motivos - stacking ABA
tau4 = a*np.array([0, 0, d_WW/a])
tau5 = a*np.array([1/np.sqrt(3), 0, (d_WW + d_W_Se)/a])
tau6 = a*np.array([1/np.sqrt(3), 0, (d_WW - d_W_Se)/a])

TAUS_L2 = [tau4, tau5, tau6]
Camada1 = []
Camada2 = []

#for i in range(-M, N):
#    for j in range(M, N+M):

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

                element = str(itau + 1)  # 1: W, 2: Se, 3: Se

                linha = f'{element}   {x_rot}   {y_rot}   {z}\n'
                Camada1.append(linha)
                Arq_out.append(linha)
                ATOMS_L1.append(linha)
                ATOMS_L1_TYPES.append(element)

            else:

                linha = f'7   {x_rot}   {y_rot}   {z}\n'
                Camada1.append(linha)

# Segunda camada

for i in range(-teste_ind, teste_ind):
    for j in range(-teste_ind, teste_ind):

        for itau, tau in enumerate(TAUS_L2):

            x , y, z = i*a1 + j*a2 + tau

            # Rotacionando
            x_rot = x*np.cos(-theta/2.0) + y*np.sin(-theta/2.0)
            y_rot = -x*np.sin(-theta/2.0) + y*np.cos(-theta/2.0)

            if dentro_supercelula(x_rot, y_rot, L_supercell) == True:
                
                element = str(itau + 4)  # 4: W, 5: Se, 6: Se
                
                linha = f'{element}   {x_rot}   {y_rot}   {z}\n'
                Camada2.append(linha)
                Arq_out.append(linha)
                ATOMS_L2.append(linha)
                ATOMS_L2_TYPES.append(element)
            else:

                linha = f'7   {x_rot}   {y_rot}   {z}\n'
                Camada2.append(linha)


# Arquivo para visualizacao

Total_atomos_vis = len(Camada1) + len(Camada2)

arq_visualizacao.write(str(Total_atomos_vis)+'\n\n')

for i in range(len(Camada1)):
    arq_visualizacao.write(Camada1[i])

for i in range(len(Camada2)):
    arq_visualizacao.write(Camada2[i])

arq_visualizacao.close()

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
