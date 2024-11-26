import numpy as np
import matplotlib.pyplot as plt

import datetime
current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

#fixed well size and seperation of 1
#vary the number of wells and record the energy for each eigenstate

num_well = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
E1 = np.array([27700.24, 27727.97, 27737.22, 27741.85, 27744.63, 27746.48, 27747.80, 27748.89, 27749.56, 27750.18])
E2 = np.array([27708.83, 27727.97, 27737.22, 27741.85, 27744.63, 27746.48, E1[6], E1[7], E1[8], E1[9]]) #can copy E1 for any past 6 
E3 = np.array([27721.34, 27736.45, 27737.25, 27741.88, 27744.63, 27746.48, E1[6], E1[7], E1[8], E1[9]])
E4 = np.array([27733.13, 27736.49, 27745.73, 27741.88, 27744.63, 27746.48, E1[6], E1[7], E1[8], E1[9]])
E5 = np.array([27741.88, 27748.26, 27745.73, 27750.32, 27744.66, 27746.51, E1[6], E1[7], E1[8], E1[9]])

Es = [E1, E2, E3, E4, E5]

colors = plt.cm.cool(np.linspace(0, 1, len(Es))) #'winter', 'cool', and 'brg' are good

for i, E in enumerate(Es):
    plt.scatter(num_well, E, color=colors[i], label=f'Eigenstate {i+1}')

plt.xlabel('Number of Wells')
plt.ylabel('Energy')
plt.legend(loc='best')
plt.savefig(f'Outputs/Energy/{current_time}_energy_num_analysis.png', dpi=800)
#plt.show()
plt.close()

#now fix number of wells at 2, size at 1
#vary the seperation to investigate the effect of tunneling on energy

""" doing it manually takes forever and doesn't give a great picture, removed functionality for more continuous calculation
seperations = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 3.5, 4, 4.5, 5])

E1 = np.array([27700.16, 24875.48, 24875.49, 27727.97, 27727.97, 25971.74, 25971.71, 27737.25, 27737.22, 26463.63, 27741.85, 26742.92, 27744.63, 26922.92, 27746.51])
E3 = np.array([27708.72, 24883.98, 24884.02, 27736.45, 27736.45, 25980.27, 25980.17, 27745.81, 27745.71, 26472.13, 27750.34, 26751.41, 27753.12, 26931.38, 27755.07])
E5 = np.array([27722.58, 24896.33, 24896.00, 27748.26, 27748.22, 25992.05, 25991.94, 27757.59, 27757.46, 26483.90, 27762.09, 26763.15, 27764.87, 26943.11, 27766.85])

Es = [E1, E3, E5]

colors = plt.cm.cool(np.linspace(0, 1, len(Es))) #'winter', 'cool', and 'brg' are good

for i, E in enumerate(Es):
    plt.scatter(seperations, E, color=colors[i], label=f'Eigenstate {2*(i+1)-1}')

plt.xlabel('Well seperation')
plt.ylabel('Energy')
plt.legend(loc='best')
plt.show()
"""

from main import solve_well
seperations = np.linspace(0, 5.1, 50)

Es = []

for i, x in enumerate(seperations):
    value, vector = solve_well(1, x, 2, -25)
    Es.append(value[0]) #ground state energy

plt.scatter(seperations, Es, c='black', s=10)

plt.xlabel('Well seperation')
plt.ylabel('Energy')
plt.savefig(f'Outputs/Energy/{current_time}_energy_separation_analysis.png', dpi=800)
#plt.show()
plt.close()