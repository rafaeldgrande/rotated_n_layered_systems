
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
plt.rc('lines', linewidth=1.0)
plt.rc('text', usetex=True)

diretorios = ['11_10', '27_26', '26_25', '20_19', '16_15', '28_27',
              '21_20', '3_2', '4_3', '5_4', '6_5', '7_6', '8_7', '9_8', '14_13', '15_14',
              '17_16', '18_17', '19_18', '22_21', '23_22', '24_23', '25_24', '2_1']

E_AB = -7.4379866430036925

#E_per_atom = [E_AB]
#Cos_angle = [1.0]
#Angle = [0.0]

E_per_atom = []
Cos_angle = []
Angle = []

for dir in diretorios:

    command = os.popen("grep 'Energy initial' "+dir+"/log -A 1 | awk {'print $3'} | tail -1")
    E = float(command.read().split('\n')[0])

    command = os.popen("grep atoms "+dir+"/*top | awk {'print $1'} | tail -1")
    Tot_atoms = float(command.read().split(':')[1])

    E_per_atom.append(E/Tot_atoms)

    command = os.popen("grep degrees "+dir+"/*top | awk {'print $2'}")
    theta = float(command.read().split('\n')[0])*np.pi/180.0
    Cos_angle.append(np.cos(theta))
    Angle.append(theta*180/np.pi)

    print(dir, E/Tot_atoms, theta)


plt.figure(1)
theta_c = 1.1
#plt.plot([np.cos(theta_c*np.pi/180), np.cos(theta_c*np.pi/180)], [min(E_per_atom), max(E_per_atom)], 'k--')

plt.plot(Cos_angle, E_per_atom, 'ro')
plt.xlabel(r'$cos(\theta)$')
plt.ylabel(r'$E/atom (eV)$')

plt.figure(2)
plt.plot(Angle, E_per_atom, 'ro')
plt.xlabel(r'$\theta$', fontsize=20)
plt.ylabel(r'$E/atom (eV)$', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

#plt.xticks([0.0, 1.0, 1.3, 2.0, 3.0])
plt.grid(axis='x')

#for theta_c in [1.1, 1.2, 1.3]:
#    plt.plot([theta_c, theta_c], [min(E_per_atom), max(E_per_atom)], 'k--')

plt.show()
