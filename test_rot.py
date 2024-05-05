from matrix_op import *
import USTCON as st
from utils_graph import *

def GraphEquality(Ga, Gb, N):
    for i in range(N):
        for j in range(N):
            if Ga[i][j] != Gb[i][j]:
                print(i,j)
                return False
    return True

def test_regularize_rot():
    """
    Compare the regularisation obtained using matrix or rotation map
    """
    N = 8
    G = random_connected_graph(N)

    trueGreg = regularize(G)
    
    rotGreg = st.rotRegularize(G, N)
    Greg = rotMap_to_adjacencyMatrix(rotGreg, N*N, 3)

    return GraphEquality(Greg, trueGreg, N*N)

def test_zig_zag_rot():
    """
    Compare the zigzag product obtained using matrix or rotation map
    """
    N = 8
    D = 3
    d = 2 # We are using cycle graph here
    G = random_regular_graph(N, D)
    H = cycle_graph(D)

    trueGoH = zigzagProductMatrix(G, H)

    rotG = adjacencyMatrix_to_rotMap(G)
    rotH = adjacencyMatrix_to_rotMap(H)

    rotGoH = st.rotZigzagProduct(rotG, rotH, D, d)
    GoH = rotMap_to_adjacencyMatrix(rotGoH, N*D, d*d)

    return GraphEquality(GoH, trueGoH, N*D)

def test_power_rot():
    """
    Compare the power obtained using matrix or rotation map
    """
    N = 8
    D = 3
    n = 5
    
    G = random_regular_graph(N, D)

    truepowG = graphPowerMatrix(G, n)

    rotG = adjacencyMatrix_to_rotMap(G)

    rotpowG = st.rotGraphPower(rotG, D, n)
    powG = rotMap_to_adjacencyMatrix(rotpowG, N, D**n)

    return GraphEquality(powG, truepowG, N)


def test_main_transform_rot():
    """
    Compare the main transform using matrix or rotation map
    """
    N = 8
    D = 3
    n = 5
    d = 2

    G = random_regular_graph(N, D)
    H = cycle_graph(3)

    trueGt = mainTransformMatrix(G, H, n)

    rotG = adjacencyMatrix_to_rotMap(G)
    rotH = adjacencyMatrix_to_rotMap(H)

    rotGt = st.rotMainTransform(rotG, rotH, D, d, n)
    Gt = rotMap_to_adjacencyMatrix(rotGt, N*D, (d*d)**n)

    return GraphEquality(Gt, trueGt, N*D)


if __name__ == "__main__":
    print(test_regularize_rot())
    print(test_power_rot())
    print(test_zig_zag_rot())
    print(test_main_transform_rot())