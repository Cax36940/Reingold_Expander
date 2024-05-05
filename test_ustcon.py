import USTCON as st
from utils_graph import *
import numpy as np

def test_exp():
    
    N = 8
    G = random_connected_graph(N)
    rotGreg = st.rotRegularize(G, N)

    D = 67*67 # Can't really be changed because of H
    D16 = D**16
    n = 8
    
    l = 1
    Gexp_test = st.rotMainTransform(rotGreg, st.rotH_complete, D16, D, n)
    Gexp_simple = st.rotGexp_simple(rotGreg, st.rotH_complete, l, D16, D, n)
    Gexp_reingold = st.rotGexp_Reingold(rotGreg, st.rotH_complete, D, l)


    vertex = np.random.randint(0,N*N)
    for i in range(l):
        vertex *= D16
        vertex += np.random.randint(0, 10000) # function doesn't suport int above int64 limit
    
    edge = np.random.randint(0, 10000)

    new_vertex0, new_edge0 = Gexp_test(vertex, edge)
    new_vertex1, new_edge1 = Gexp_simple(vertex, edge)
    new_vertex2, new_edge2 = Gexp_reingold(vertex, edge)

    print(vertex, edge)
    print(new_vertex0, new_edge0)
    print(new_vertex1, new_edge1)
    print(new_vertex2, new_edge2)

    is_equal = (new_vertex1 == new_vertex2) and (new_edge1 == new_edge2)

    return is_equal



if __name__ == "__main__":
    print(test_exp())

