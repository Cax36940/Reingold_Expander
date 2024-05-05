import USTCON as st
from utils_graph import *
import numpy as np

def test_exp(l = 3, N = 20):

    G = random_connected_graph(N)
    rotGreg = st.rotRegularize(G, N)

    D = 67*67 # Can't really be changed because of H
    D16 = D**16
    n = 8

    # Test our high level algorithm
    Gexp_simple = st.rotGexp_simple(rotGreg, st.rotH_complete, l, D16, D, n)

    # Test low-level Reingold algorithm
    Gexp_reingold = st.rotGexp_Reingold(rotGreg, st.rotH_complete, D, l)


    vertex = np.random.randint(0,N*N)
    for i in range(l):
        vertex *= D16
        vertex += np.random.randint(0, 10000) # function doesn't suport int above int64 limit
    
    edge = np.random.randint(0, 10000)

    new_vertex1, new_edge1 = Gexp_simple(vertex, edge)
    new_vertex2, new_edge2 = Gexp_reingold(vertex, edge)

    # print(vertex, edge)
    print(new_vertex1, new_edge1)
    print(new_vertex2, new_edge2)

    is_equal = (new_vertex1 == new_vertex2) and (new_edge1 == new_edge2)

    return is_equal

def test_exp_n_times(nb_iter):

    N = 8
    D = 67*67 # Can't really be changed because of H
    D16 = D**16
    n = 8
    l = 3    
    
    average = 0
    for i in range(nb_iter):
        G = random_connected_graph(N)
        rotGreg = st.rotRegularize(G, N)

        # Test our high level algorithm
        Gexp_simple = st.rotGexp_simple(rotGreg, st.rotH_complete, l, D16, D, n)

        # Test low-level Reingold algorithm
        Gexp_reingold = st.rotGexp_Reingold(rotGreg, st.rotH_complete, D, l)


        vertex = np.random.randint(0,N*N)
        for i in range(l):
            vertex *= D16
            vertex += np.random.randint(0, 10000)

        edge = np.random.randint(0, 10000)

        new_vertex1, new_edge1 = Gexp_simple(vertex, edge)
        new_vertex2, new_edge2 = Gexp_reingold(vertex, edge)

        if (new_vertex1 == new_vertex2) and (new_edge1 == new_edge2) :
            average += 1
    print(average/nb_iter)

if __name__ == "__main__":
    print(test_exp())
    # test_exp_n_times(100)

