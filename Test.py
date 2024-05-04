
from reingold import *
import USTCON as st

def test_main_transform_step(N, D) :

    # H is the complete graph on D vertices
    H = complete_graph(D)
    print(H)
    # print(H)

    # G is a random 3-regular graph
    G = random_graph(N, D)
    print(G)

    lamb = secondEV(G)

    # computes the zigzag product between G and H
    GoH = zigzagProductMatrix(G, H)

    # 8th power of the resulting graph
    GoH8 = graphPowerMatrix(GoH, 8)

    lamb2 = secondEV(GoH8)

    print("λ(G)² : " + str(lamb**2))
    print("λ(GoH⁸) : " + str(lamb2))

# test_main_transform_step(8, 3)

def test_zig_zag_product():
    G = random_graph(6,4)
    #nx.draw(nx.from_numpy_array(G))
    print(G)
    H = np.array([[0.,0.,1.,0.],[0.,0.,0.,1.],[1.,0.,0.,0.],[0.,1.,0.,0.]])
    #nx.draw(nx.from_numpy_array(H))
    GoH = zigzagProductMatrix(G,H)
    print(GoH)

    nx.draw(nx.from_numpy_array(GoH))
    plt.show()

def zigzag_eigen(lamb, alph):
    a = (1 - alph**2)
    b = sqrt(a*a * lamb*lamb + 4*alph**2)
    return 1/2*(a * lamb + b)

def test_eigen_value_zigzag_product():
    N = 8
    D = 3
    # H is the complete graph on D vertices
    H = complete_graph(D)
    print(H)
    # print(H)

    # G is a random 3-regular graph
    G = random_graph(N, D)
    print(G)

    lamb = secondEV(G)
    alpha = secondEV(H)
    # computes the zigzag product between G and H
    GoH = zigzagProductMatrix(G, H)

    lamb = secondEV(GoH)

    lambth = zigzag_eigen(lamb, alpha)

    print("eigen value with zigzag product :\t",lamb)
    print("eigen value with theoritical formula :\t",lambth)

def test_eigen_value_power():
    N = 8
    D = 3

    # G is a random 3-regular graph
    G = random_graph(N, D)
    print(G)

    lamb = secondEV(G)
    # computes the zigzag product between G and H
    G8 = graphPowerMatrix(G, 8)

    lamb8 = secondEV(G8)

    lambth = lamb ** 8

    print("eigen value with 8th power : \t\t", lamb8)
    print("eigen value with theoritical formula : \t", lambth)

def test_rotMap():
    G = np.array([[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]]) #cycle de 4 sommets

    f = st.rot(G)

    for i in range(4):
        for j in range(2):
            print((i,j + 1),"->",f(i,j + 1))

def test_rotMap_to_Amatrix():
    G = np.array([[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]]) #cycle de 4 sommets

    f = st.rot(G)

    H = st.RotMap_to_adjacenceMatrix(f,2,len(G))

    print(G)
    print(H)

test_rotMap_to_Amatrix()
#test_eigen_value_zigzag_product()
# test_eigen_value_power()
# test_zig_zag_product()