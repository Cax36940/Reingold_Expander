import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt ##### Useful Functions #####

def random_graph(N, D):
    """
    Returns a random D-regular connected graph on N vertices

    N : int
    D : int

    return G : graph
    """
    # Generate a random D-regular graph
    G = nx.random_regular_graph(D, N)
    while not nx.is_connected(G): # Verify if it is connected
        G = nx.random_regular_graph(D, N)
    return nx.to_numpy_array(G)


def graphPowerRotationMap(G, n):
    """
        Computes G^n

        G : graph
        n : int

        return : graph G'
    """

def graphPowerMatrix(G, n):
    """
        Computes G^n

        G : graph (int[][])
        n : int

        return : graph G'
    """
    return np.linalg.matrix_power(G, n)

def zigzagProduct(G, H):
    """
        Computes G Ⓩ H, the zigzag product between G and H

        G : graph
        H : graph

        return : G' : graph
    """


def neighbours(G, v) :
    """
    Return the neighbours of vertex v in the graph G

    G : graph (int[][])
    v : int

    return : neigh : int[]

    """
    neigh = []
    for i in range(len(G)) :
        for j in range(int(G[v][i])):
            neigh.append(i)
    
    return neigh

def zigzagProductMatrix(G, H):
    """
        Computes G Ⓩ H, the zigzag product between G and H

        G : graph (int[][])
        H : graph (int[][])

        return : G' : graph (int[][])
    """

    degG = int(sum(G[0]))

    if degG != len(H):
        print("Graph H don't have the right size, must be equal to deg of G")
        return [[]]

    newN = int(len(G) * degG)

    GoH = np.zeros((newN, newN))

    for i in range(newN):
        iG = int(i // degG)
        iH = int(i % degG)
        neigh_iG = neighbours(G, iG) # Neighbours of i in G
        neigh_iH = neighbours(H, iH) # Neighbours of i in H (short edge)

        jG_list = [neigh_iG[k] for k in neigh_iH] # (long edge)

        for jG in jG_list:
            for jH in neighbours(H, neighbours(G, jG).index(iG)): # (short edge)
                GoH[i][jG * degG + jH] += 1
                GoH[jG * degG + jH][i] += 1
    
    return GoH / 2 # Every edge is counted 2 times here


def secondEV(G):
    """
        Computes the second largest eigenvalue of G

        G : graph (int[][])

        return : λ : float
    """
    eigen = np.linalg.eigvals(G) # eigenvalues of G
    eigen = map(abs, eigen) # take absolute value
    eigen = sorted(eigen) # sort
    return eigen[-2]/sum(G[0]) # return second greatest normalized value

def pathLength(N):
    """
        Computes the max length of a log path in a (N, D, 1/2) expander

        N : int

        return : Delta : int
    """
    return 4 * int(np.ceil(N))

def maxPower(N, D):
    """
        Computes the max number of iteration for the Reingold Algorithm

        N : int
        D : int

        return : l : int
    """ 
    return int(np.ceil(-np.log2(-np.log2(1-1/(D*N*N)))))


def rotationMapConstructor(G):
    """
        Returns the rotationMap function for the graph G

        G : graph (int[][])
    
        return : (int, int) -> ()
    """

##### USTCON #####

# Notation N or D instead of [|1, N|] or [|1, D|]

def rotH_simple(vertex, edge, nb_vertices = 3) :
    """
    Rotation map for cycle on nb_vertices vertices
    Uses pi-permutation (0, 1) -> (1, 0)

    vertex : int    in [0, nb_vertices]
    edge :   int    in [0, 1]
    nb_vertices : int (default value = 3)
    
    return (vertex', edge') : (int, int)
    
    """
    if edge == 0 :
        return (vertex-1)%nb_vertices, 1
    return (vertex-1)%nb_vertices, 0

def rotH_complete(vertex, edge):
    """
        Rotation map of (D^16, D, 1/2) graph H

        vertex : int[16]    in D^16
        edge : int       in D

        return : (vertex', edge') : (int[16], int[2])
    """

    q = 67
    d = 31
    x = edge // q
    y = edge % q

    # convert vertex form int[16] in D^16 to int[32] in q^32
    a = [0 for i in range(d + 1)]

    for i in range(32):
        if i % 2 :
            a[i] = vertex[i//2] % q
        else :
            a[i] = vertex[i//2] // q
        
    b = [(y * (x**i))%q for i in range(d + 1)]

    c = [(a[i] + b[i])%q for i in range(d + 1)]

    new_vertex = [0 for i in range(16)]

    for i in range(16):
        new_vertex[i] = c[2*i] * q + c[2*i+1]

    return [new_vertex, x*q + ((-y)%q)]


def rotGreg(vertex, edge, G):
    """
        Rotation map of N² vertex regular graph based on G

        G :      graph (int[][])
        vertex : (int, int) in N x N
        edge :   int
        
        return : (vertex', edge') : ((int, int), int)
    """
    N = len(G)
    if edge == 1:
        return (vertex[0],(vertex[1]+1)%N), 2
    if edge == 2:
        return (vertex[0],(vertex[1]-1)%N), 1
    if G[vertex[0]][vertex[1]] :
        return (vertex[1],vertex[0]), edge

    return vertex, edge

#TODO Test une itération de Reingold avec une H = 3-cycle
#       on pourrait faire le test avec plusieurs graphes de degré 3
def rotGexp_complete(vertex, edge, G, rotH, D):
    """
        Rotation map of (N² x (D^16)^L, D^16, 1/2) expander graph

        global const D
        global const N

        dependancies :
            rotH
            rotGreg

        vertex : [(int, int), int[L][16]]   in N² x (D^16)^L
        edge : int[16]                      in D^16



        return : (vertex', edge') : ([(int, int), int[L][16]], int[16])
    """

    N = len(G)
    # l = maxPower(N, D)
    l = 1
    # Allocate variables
    a = [[0 for i in range(16)] for k in range(l + 1)]
    I = 0
    j = [0 for i in range(l + 1)]

    # Initialise variables
    for i in range(l):
        for k in range(16):
            a[i][k] = vertex[i][k]
    for k in range(16):  
        a[l][k] = edge[k]

    I = l

    for i in range(1, l + 1):
        j[i] = 1
    

    while 1 :

        # Step 1
        # print(a)
        new_a = rotH(a[I-1],a[I][j[I]-1])
        a[I-1] = new_a[0]
        a[I][j[I]-1] = new_a[1]

        # Step 2
        if j[I]%2:
            if I == 1 :
                if a[0][:-2] == [1 for i in range(15)]:
                    if a[0][15] <= 3 :
                        print((a[0][15], piPerm(a[0][15])))
                        a[0] = piPerm(a[0])
                j[I] += 1

        # Step 3
            elif I > 1 :
                j[I-1] = 1
                I -= 1
         
        # Step 4
        elif j[I] == 16 :
            #reverse a[I]
            for i in range(8):
                tmp = a[I][i]
                a[I][i] = a[I][15 - i]
                a[I][15 - i] = tmp
        
        # Step 5
            if I == l :
                break

        # Step 6
            elif I < l :
                j[I + 1] += 1
                I += 1
        else:
            j[I] += 1

    return path, vertex


def Aexp(s, t):
    """
        Return true if s and t are connected in Gexp

        global const D
        global const N

        computes Delta

        dependancies :
            rotGexp

        s : [(int, int), int[L][16]]   in N² x (D^16)^L
        t : [(int, int), int[L][16]]   in N² x (D^16)^L

        return : con : boolean
    """


def Acon(s, t):
    """
        Return true if s and t are connected in Gexp, implying connexity in G

        global const D
        global const N

        dependancies :
            Aexp

        s : (int, int)  in N x N
        t : (int, int)  in N x N

        return : con : boolean
    """




##### UTS / UXS #####

def piPerm(x):
    """
        Compute the pi permutation (2,1,...) on x

        x : int

        return : y : int
    """
    if x == 1 :
        return 2
    if x == 2 :
        return 1
    return x


def Ae2p(vertex, edge, N):
    """
        Compute the path, in a graph with pi-permutation, taken while walking on an edge in its expander transformation

        global const D

        dependancies :
            maxPower
            rotH
            piPerm

        vertex : int[L][16]
        edge : int[16]                      in D^16
        N : int

        (computations are similar to the ones of rotGexp)

        return : path : int[][16]     in (D^16)*
                 next_vertex : int[L][16]

    """
    path = []
    global D
    l = maxPower(N, D)
    l = 6
    # Allocate variables
    a = [[0 for i in range(16)] for k in range(l + 1)]
    I = 0
    j = [0 for i in range(l + 1)]

    # Initialise variables
    for i in range(l):
        for k in range(16):
            a[i][k] = vertex[i][k]
    for k in range(16):  
        a[l][k] = edge[k]

    I = l

    for i in range(1, l + 1):
        j[i] = 1
    

    while 1 :

        # Step 1
        # print(a)
        new_a = rotH(a[I-1],a[I][j[I]-1])
        a[I-1] = new_a[0]
        a[I][j[I]-1] = new_a[1]

        # Step 2
        if j[I]%2:
            if I == 1 :
                if a[0][:-2] == [1 for i in range(15)]:
                    if a[0][15] <= 3 :
                        print((a[0][15], piPerm(a[0][15])))
                        a[0] = piPerm(a[0])
                j[I] += 1

        # Step 3
            elif I > 1 :
                j[I-1] = 1
                I -= 1
         
        # Step 4
        elif j[I] == 16 :
            #reverse a[I]
            for i in range(8):
                tmp = a[I][i]
                a[I][i] = a[I][15 - i]
                a[I][15 - i] = tmp
        
        # Step 5
            if I == l :
                break

        # Step 6
            elif I < l :
                j[I + 1] += 1
                I += 1
        else:
            j[I] += 1

    return path, vertex





def Asrch(path, N, L): #not really what was defined in the paper
    """
        Compute path in a graph with pi-permutation, taken while walking on Delta an edges in its expander transformation

        dependancies :
            piPerm
            Ae2p

        path : path : int[L][16]     in (D^16)^Delta
        N : int

        calls Ae2p for each edge of the path

        return : path : int[][16]     in (D^16)*
    """
    entire_path = []
    vertex = [0 for i in range(L)]

    for edge in path :
        tmp_path, vertex = Ae2p(vertex, edge, N)
        entire_path = entire_path + list(tmp_path)
    

    entire_path = entire_path + [piPerm(x) for x in entire_path].reverse()
    return entire_path


def UTS(N, Dmax):
    """
        Compute UTS on N*Dmax vertex D^16 regular graph with pi permutation labelling

        N : int
        Dmax : int

        dependancies :
            Asrch

        calls Asrch on every path of length Delta of edges in D^16

        return : UTS : int[][16] in D^16
    """

    Delta = pathLength(N, Dmax)
    L = maxPower(N, Dmax)
    uts = []
    for i in range((D**16)**Delta):
        path = [(i//((D**16)**k))%D**16 for k in range(L)]
        uts = uts + Asrch(path, N, L)

    return uts


def toXS(TS):
    """
        Convert TS to XS considering (2, 1, ...) pi permutation

        global const D

        UTS : int[][16]         in (D^16)*

        return : XS : int[]
    """




    return 0


def UXS(N, Dmax):
    """
        Compute UXS on N vertex graph with max degre Dmax

        global const D

        N : int
        Dmax : int

        dependancies :
            UTS
            
            toXS

        calls UTS
        calls toXS

        return : UXS : int[] in Dmax*
    """
    return toXS(UTS(N, Dmax))
#TODO faire un graph, sa transformation en TS en XS
# print(maxPower(8, 67*67))
# Ae2p([[1 for i in range(16)] for k in range(6)], [1 for i in range(16)], 8)

def complete_graph(n):
    H = np.zeros((n,n))
    for i in range(n):
        H[i][(i+1)%n] = 1
        H[(i+1)%n][i] = 1
    
    return H

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


test_eigen_value_zigzag_product()
# test_eigen_value_power()
# test_zig_zag_product()