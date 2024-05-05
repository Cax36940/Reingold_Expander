from matrix_op import *
from utils import *
import matplotlib.pyplot as plt


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


def test_eigen_main_transform_graph():
    D = 3
    H = complete_graph(D)
    nb_iter = 100



    X = [x for x in range(4, 51, 2)]
    lambG_list = []
    lamb2G_list = []
    lambGoH8_list = []
    for N in X:
        sum_lambG = 0
        sum_lamb2G = 0
        sum_lambGoH8 = 0
        for i in range(nb_iter):
            print(N, i)
            G = random_regular_graph(N, D)
            GoH = zigzagProductMatrix(G, H)
            GoH8 = graphPowerMatrix(GoH, 8)

            lambG = secondEV(G)
            lamb2G = lambG * lambG
            lambGoH8 = secondEV(GoH8)

            sum_lambG += lambG
            sum_lamb2G += lamb2G
            sum_lambGoH8 += lambGoH8

        lambG_list.append(sum_lambG/nb_iter)
        lamb2G_list.append(sum_lamb2G/nb_iter)
        lambGoH8_list.append(sum_lambGoH8/nb_iter)
    
    plt.plot(X, lambG_list, label='λ(G)')
    plt.plot(X, lamb2G_list, label='λ(G)²')
    plt.plot(X, lambGoH8_list, label='λ(GoH⁸)')

    plt.xlabel('Number of Vertices', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=20)
    plt.show()





if __name__ == "__main__":
    test_eigen_value_zigzag_product()
    test_eigen_value_power()
    test_eigen_value_main_transform()
    # test_n_times_powering(1000)
    # test_n_times_zigzag_product(1000)
    # test_n_times_main_transformation(1000)
    test_eigen_main_transform_graph()