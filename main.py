import numpy as np
import matplotlib.pyplot as plt

def V(x, V0, a, L, n_well):

    #returning the potential for a point x
    #we construct the well to be symmetric so we only consider positive x
    
    #getting getting total length
    well_length = n_well * a + (n_well-1) * L

    #returning "infinity" beyond ends
    if np.abs(x) >= well_length/2:
        return 10

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
        
well_width = 5
well_dist = 2
well_num = 6
well_length = well_num * well_width + (well_num-1) * well_dist

V0 = -10

xs = np.linspace(-3*well_length/2, 3*well_length/2, 1000)
Vs = [V(x, V0, well_width, well_dist, well_num) for x in xs]

plt.plot(xs, Vs)
plt.show()