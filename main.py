import numpy as np
import matplotlib.pyplot as plt

def V(x, V0, a, L, n_well):

    #returning the potential for a point x
    #we construct the well to be symmetric so we only consider positive x

    #case for even number of wells
    if n_well % 2 == 0:

        #getting getting total length
        well_length = n_well * a + (n_well-1) * L

        #returning infinity beyond ends
        if np.abs(x) >= well_length/2:
            return np.inf

        #list of integers to index wells
        ns = np.arange(1, n_well//2+1)
        
        #negative edges defined as where the well falls going in the positive direction, positive edges where it rises
        neg_edges = L/2 + (ns-1) * L + (ns-1)*a
        pos_edges = L/2 + a + (ns-1) * L + (ns-1)*a

        #checking if our point is between a negative and positive edge for each well in the positive quadrant
        for i in range(n_well//2):
            if np.abs(x) > neg_edges[i] and np.abs(x) < pos_edges[i]:
                return V0

        #potential is zero if not in a well or outside the region
        return 0

    #case for odd number of wells
    else:
        print('odd')

well_width = 2
well_dist = 1
well_num = 4
well_length = well_num * well_width + (well_num-1) * well_dist

V0 = -10

xs = np.linspace(-3*well_length/2, 3*well_length/2, 1000)
Vs = [V(x, V0, well_width, well_dist, well_num) for x in xs]

plt.plot(xs, Vs)
plt.show()