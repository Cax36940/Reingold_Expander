from matrix_op import *
from utils import *

def test_eigen_value_power():
    N = 8
    D = 3

    # G is a random 3-regular graph
    G = random_regular_graph(N, D)
    # print(G)

    lamb = secondEV(G)
    # computes the zigzag product between G and H
    G8 = graphPowerMatrix(G, 8)

    lamb8 = secondEV(G8)

    lambth = lamb ** 8

    print("eigen value with 8th power : \t\t", lamb8)
    print("eigen value with theoritical formula : \t", lambth)

def zigzag_eigen(lamb, alph):
    a = (1 - alph**2)
    b = sqrt(a**2 * lamb**2 + 4*alph**2)
    return 1/2*(a * lamb + b)

def test_eigen_value_zigzag_product():
    N = 8
    D = 3
    # H is the complete graph on D vertices
    H = complete_graph(D)
    # print(H)

    # G is a random 3-regular graph
    G = random_regular_graph(N, D)
    # print(G)

    lamb = secondEV(G)
    alpha = secondEV(H)
    # computes the zigzag product between G and H
    GoH = zigzagProductMatrix(G, H)

    lamb = secondEV(GoH)

    lambth = zigzag_eigen(lamb, alpha)

    print("eigen value with zigzag product :\t",lamb)
    print("eigen value with theoritical formula :\t",lambth)

    return lambth >= lamb 

def test_eigen_value_main_transform() :
    N = 8
    D = 3

    # H is the complete graph on D vertices
    H = complete_graph(D)
    # print(H)

    # G is a random 3-regular graph
    G = random_regular_graph(N, D)
    # print(G)

    lamb = secondEV(G)

    # computes the zigzag product between G and H
    GoH = zigzagProductMatrix(G, H)

    # 8th power of the resulting graph
    GoH8 = graphPowerMatrix(GoH, 8)

    lamb2 = secondEV(GoH8)

    print("λ(G)² : " + str(lamb**2))
    print("λ(GoH⁸) : " + str(lamb2))
    return lamb**2 >= lamb2

def test_n_times_powering(n):
    N = 8
    D = 3
    H = complete_graph(D)

    accurate = 0
    for i in range(n):
        G = random_regular_graph(N,D)
        lamb = secondEV(G)
        G8 = graphPowerMatrix(G, 8)
        lamb8 = secondEV(G8)
        lambth = lamb ** 8
        if abs(lamb8 - lambth) < 10E-10 :# We try to avoid floating point error
            accurate +=1 
    print(accurate/n)


def test_n_times_zigzag_product(n):
    accurate = 0

    N = 8
    D = 3
    
    H = complete_graph(D)

    for i in range(n):
        G = random_regular_graph(N, D)

        lamb = secondEV(G)
        alpha = secondEV(H)

        GoH = zigzagProductMatrix(G, H)

        lamb = secondEV(GoH)

        lambth = zigzag_eigen(lamb, alpha)

        if lambth >= lamb :
            accurate += 1
    print(accurate/n)

def test_n_times_main_transformation(n):
    accurate = 0
    N = 8
    D = 3

    for i in range(n):
        H = complete_graph(D)

        G = random_regular_graph(N, D)

        lamb = secondEV(G)

        GoH = zigzagProductMatrix(G, H)

        GoH8 = graphPowerMatrix(GoH, 8)

        lamb2 = secondEV(GoH8)
        if (lamb**2 >= lamb2):
            accurate += 1

    print(accurate/n)


if __name__ == "__main__":
    test_eigen_value_zigzag_product()
    test_eigen_value_power()
    test_eigen_value_main_transform()
    test_n_times_powering(1000)
    test_n_times_zigzag_product(1000)
    test_n_times_main_transformation(1000)