import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

def V(x, V0, a, L, n_well):

    #returning the potential for a point x
    #we construct the well to be symmetric so we only consider positive x
    
    #getting getting total length
    well_length = n_well * a + (n_well-1) * L

    #returning "infinity" beyond ends
    if np.abs(x) >= well_length/2:
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

#list of parameters to test
widths = np.array([1,2]) #well widths
dists = np.array([1,1]) #distances between wells
nums = np.array([3,2]) #number of wells
V0s = np.array([-25,-25]) #potentials
lengths = np.array(nums * widths + (nums-1) * dists)
num_params = len(nums)

#number of discrete points
N = 1000

#create x array and compute potential for each set of parameters
Xs = np.array([np.linspace(-3*lengths[i]/2, 3*lengths[i]/2, N) for i in range(num_params)])
Vs = np.array([[V(x, V0s[i], widths[i], dists[i], nums[i]) for x in Xs[i]] for i in range(num_params)])

#solve each set of parameters
solutions = [numeric_solver(xs, vs) for xs, vs in zip(Xs, Vs)]
eigenvalues = [solution[0] for solution in solutions]
eigenvectors = [solution[1] for solution in solutions]

#plotting all potentials together
fig, axs = plt.subplots(num_params, 1, sharex=True)

for i in range(num_params):

    axs[i].plot(Xs[i], Vs[i], label="Potential V(x)", c='black')
    
    for j in range(nums[i]):
        axs[i].plot(Xs[i], eigenvectors[i][:, j], 
                    label=f"Eigenstate {j+1}, E={eigenvalues[i][j]:.2f}")

    axs[i].set_xlim((-1.1*max(lengths)/2, 1.1*max(lengths)/2))
    axs[i].set_ylim((1.1*V0s[i], -0.3*V0s[i]))
    axs[i].set_title(f"Potential Well Configuration {i + 1}", loc='left', x=-0)
    axs[i].set_ylabel("Energy")
    axs[i].legend(loc='lower right')

axs[-1].set_xlabel("Position")
plt.tight_layout()
plt.show()

#plotting all eigenstates of one potential seperately
for i in range(num_params):
    fig, axs = plt.subplots(nums[i], 1, sharex=True)
    
    for j in range(nums[i]):

        axs[j].plot(Xs[i], Vs[i], c='black')
        axs[j].plot(Xs[i], eigenvectors[i][:, j])
        
        axs[j].set_xlim((-1.1 * lengths[i] / 2, 1.1 * lengths[i] / 2))
        axs[j].set_ylim((1.1 * V0s[i], -0.3 * V0s[i]))
        
        axs[j].set_title(f"Eigenstate Energy {eigenvalues[i][j]}", loc='left')
        axs[j].set_ylabel("Energy")
    
    axs[-1].set_xlabel("Position")

    plt.tight_layout()
    plt.show()