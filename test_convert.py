import USTCON as st
from utils_graph import *


##### Test Convertion functions #####

def test_matrix_to_rot():
    """
        Test adjacencyMatrix_to_rotMap function on cycle graph
    """
    N = 8

    G = cycle_graph(N)

    rotG = adjacencyMatrix_to_rotMap(G)

    test = True

    test = test and (rotG(0, 0) == (1, 0))
    test = test and (rotG(0, 1) == (N-1, 0))
    test = test and (rotG(1, 0) == (0, 0))
    test = test and (rotG(1, 1) == (2, 0))
    test = test and (rotG(N-2, 0) == (N-3, 1))
    test = test and (rotG(N-2, 1) == (N-1, 1))
    test = test and (rotG(N-1, 0) == (0, 1))
    test = test and (rotG(N-1, 1) == (N-2, 1))
    
    for v in range(2, N-2):
        for i in [0, 1]:
            if i == 0 :
                test = test and (rotG(v, i) == ((v-1)%N, 1))
            else:
                test = test and (rotG(v, i) == ((v+1)%N, 0))

    return test


def test_rot_to_matrix():
    """
        Test rotMap_to_adjacencyMatrix function on cycle graph
    """
    N = 8

    rotG = st.rotH_simple(N)

    G = rotMap_to_adjacencyMatrix(rotG, N, 2)
    trueG = cycle_graph(N)

    test = True
    for i in range(N):
        for j in range(N):
            test = test and (G[i][j] == trueG[i][j])

    print(G)
    print(trueG)
    
    return test


def test_convert():
    """
    Convert a graph to its rotation map and back to matrix
    """
    N = 8
    D = 3

    G = random_regular_graph(N, D)

    # Convert to rotation map
    rotG = st.adjacencyMatrix_to_rotMap(G)

    # Convert back to ajacency matrix
    testG = st.rotMap_to_adjacencyMatrix(rotG, N, D)

    # Compare the initial graph with the resulting graph
    test = True
    for i in range(N):
        for j in range(N):
            test = test and (G[i][j] == testG[i][j])

    print(G)
    print(testG)
    
    return test


if __name__ == "__main__":
    print(test_convert())
    print(test_matrix_to_rot())
    print(test_rot_to_matrix())

