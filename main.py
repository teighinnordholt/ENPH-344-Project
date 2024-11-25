import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

def V(x, V0, a, L, n_well):

    #returning the potential for a point x
    #we construct the well to be symmetric so we only consider positive x

    #getting getting total length
    well_length = n_well * a + (n_well-1) * L

    #returning "infinity" beyond ends
    if np.abs(x) >= well_length/2 + L/2:
        return 1e10

    #case for even number of wells
    if n_well % 2 == 0:

        #list of integers to index wells
        ns = np.arange(1, n_well//2 + 1)

        #negative edges defined as where the well falls going in the positive direction, positive edges where it rises
        neg_edges = L/2 + (ns-1) * L + (ns-1)*a
        pos_edges = L/2 + a + (ns-1) * L + (ns-1)*a

    #case for odd number of wells
    else:

        #list of integers to index wells
        ns = np.arange(1, n_well//2 + 2)

        #negative edges defined as where the well falls going in the positive direction, positive edges where it rises
        neg_edges = a/2 + (ns-1) * a + ns * L; neg_edges = np.append([0], neg_edges[:-1]) #to account for first one being half the size, make the list for the rest the append them to the first 
        pos_edges = a/2 + (ns-1) * a + (ns-1) * L

    #checking if our point is between a negative and positive edge for each well in the positive quadrant
    for i, _ in enumerate(pos_edges):
        if np.abs(x) > neg_edges[i] and np.abs(x) < pos_edges[i]:
            return V0

    #potential is zero if not in a well or outside the region
    return 0

def numeric_solver(x, V):

    N = len(x)
    del_x = (x[-1] - x[0]) / (N - 1)

    #selecting nice units
    h_bar = 1
    m = 1

    H = np.zeros((N, N))

    #constructing the hamiltonian matrix
    for i in range(1, N - 1):
        H[i, i] = (h_bar**2 / (m * del_x**2)) * 2 + V[i]  #diagonals
        H[i, i - 1] = H[i, i + 1] = -h_bar**2 / (2 * m * del_x**2)  #off-diagonals

    #boundary conditions
    H[0, 0] = H[N - 1, N - 1] = (h_bar**2 / (m * del_x**2)) * 2 + V[0]

    #solving eigenvalues, vectors
    eigenvalues, eigenvectors = eigh(H)

    #sort by energy
    sorted_indices = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]

    #normalize wavefunctions
    for i in range(len(eigenvectors)):
        eigenvectors[:, i] /= np.linalg.norm(eigenvectors[:, i]) * np.sqrt(del_x)

    return eigenvalues, eigenvectors

#list of parameter
width = 1 #well size
dist = 3 #well seperation
num = 6 #number of wells
V0 = -25 #potentials
length = int(num * width + (num) * dist)

#number of discrete points, dynamic to length based on factor found to give 'sharp' wells
N = 500*length

#create x array and compute potentials
xs = np.linspace((-3*length/2)//1, (3*length/2), N)
Vs = np.array([V(x, V0, width, dist, num) for x in xs])

#solve potential
solution = numeric_solver(xs, Vs)
eigenvalues = solution[0]
eigenvectors = solution[1]

#number of eigenstates to plot
num_eigenstate = 6
colors = plt.cm.cool(np.linspace(0, 1, num_eigenstate)) #'winter', 'cool', and 'brg' are good

plt.plot(xs, Vs, label='Potential V(x)', c='black')

offset = 0

for i in range(num_eigenstate):

    if i > 0:
        offset += 1.25*(max(eigenvectors[:, i-1]) - min(eigenvectors[:,i]))
        plt.axhline(y=offset, color=colors[i], linestyle='--', linewidth=0.5)

    plt.plot(xs, eigenvectors[:, i]+offset, c=colors[i],
                label=f'Eigenstate {i+1}, E={eigenvalues[i]:.2f}')

plt.xlim((-1.25*length/2, 1.25*length/2))
plt.ylim((1.1*V0, -0.3*V0+offset))

plt.title(f'Well size = {width:.2f}, Well Seperation = {dist:.2f}', loc='left', x=-0)

plt.xlabel('Position')
plt.ylabel('Energy')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(f'Outputs/num{num}_size{width}_sep{dist}_v{abs(V0)}.png', dpi=800)
#plt.show()