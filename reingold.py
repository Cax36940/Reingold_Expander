import networkx as nx
import numpy as np

##### Useful Functions #####

def random_graph(D, N):
    """
    Returns a random D-regular connected graph on N vertices

    D : int
    N : int

    return G : graph
    """
    # Generate a random D-regular graph
    G = nx.random_regular_graph(D, N)
    while not nx.is_connected(G): # Verify if it is connected
        G = nx.random_regular_graph(D, N)
    return nx.to_numpy_array(G)


def graphPower(G, n):
    """
        Computes G^n

        G : graph
        n : int

        return : graph G'
    """


def zigzagProduct(G, H):
    """
        Computes G Ⓩ H, the zigzag product between G and H

        G : graph
        H : graph

        return : G' : graph
    """


def secondEV(G):
    """
        Computes the second largest eigenvalue of G

        G : graph

        return : λ : float
    """

def pathLength(N, D):
    """
        Computes the max length of a log path in a (N, D, 1/2) expander

        N : int
        D : int

        return : Delta : int
    """

def maxPower(N, D):
    """
        Computes the max number of iteration for the Reingold Algorithm

        N : int
        D : int

        return : l : int
    """ 
    return int(np.ceil(-np.log2(-np.log2(1-1/(D*N*N)))))

##### USTCON #####

# Notation N or D instead of [|1, N|] or [|1, D|]

D = 67 * 67 # q = 67

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
        return [(vertex-1)%nb_vertices, 1]
    return [(vertex-1)%nb_vertices, 0]

def rotH(vertex, edge):
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


    invy = -y % q

    new_vertex = [0 for i in range(16)]

    for i in range(16):
        new_vertex[i] = c[2*i] * 67 + c[2*i+1]


    return [new_vertex, x*67 + ((-y)%67)]


def rotGreg(vertex, edge):
    """
        Rotation map of N² vertex D^16 regular graph based on G

        global graph G
        global const D
        global const N

        vertex : (int, int) in N x N
        edge : int[16]      in D^16

        return : (vertex', edge') : ((int, int), int[16])
    """
    return 0

#TODO Test une itération de Reingold avec une H = 3-cycle
#       on pourrait faire le test avec plusieurs graphes de degré 3
def rotGexp(vertex, edge):
    """
        Rotation map of (N² x (D^16)^L, D^16, 1/2) expander graph

        global const D
        global const N

        computes L as in 3.1

        dependancies :
            rotH
            rotGreg

        vertex : [(int, int), int[L][16]]   in N² x (D^16)^L
        edge : int[16]                      in D^16

        return : (vertex', edge') : ([(int, int), int[L][16]], int[16])
    """
    return 0


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
        Return true if s and t are connected in Greg, implying connexity in G

        global const D
        global const N

        computes L

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

def test_main_transformation_step(D, N) :
    g = random_graph(D, N)
    # print(g)
    eigen = np.linalg.eigvals(g)
    eigen = map(abs, eigen)
    eigen = sorted(eigen)
    # print(eigen)
    lamb = eigen[-2]/D
    print("Lambda value for a random (" + str(N) + ", " + str(D) + " lambda) graph : " + str(lamb))


test_main_transformation_step(3, 8)