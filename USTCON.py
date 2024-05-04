import numpy as np
import networkx as nx

def rot(G): #G : matrice d'adjacence, on suppose qu'il n'existe pas d'arêtes parallèles et que G est régulier
    def fct(v,i):
        sum = 0
        N = len(G)
        k = 0 
        while (k < N):
            sum += G[v][k] 
            if (sum == i):
                 w = k
                 break
            k += 1
        k = 0 
        sum = 0
        while (k < N):
            sum += G[k][w]
            if (v == k):
                j = sum
                break
            k += 1
        return (w,j)
    return fct
#à tester


def RotMap_to_adjacenceMatrix(f,D,N):
    AMatrix = np.array([[0]*N]*N)
    for u in range(N):  
        for i in range(D):
            (v,j) = f(u,i)
            AMatrix[u][v] = 1
    return AMatrix

'''
def getSecondMax(tab):
    max1, max2 = tab[0], tab[0]
    for i in range(len(tab)):
        if tab[i] >= max1:
            max1 = tab[i]
        elif (tab[i] >= max2):
            max2 = tab[i]
    return max2

def getSecondEigenValue(M,D):
    newM = normalize(M,D)
    eigenvalues, eigenvectors = np.linalg.eig(newM)
    return getSecondMax(eigenvalues)


def normalize(M,D):
    n,m = M.shape()
    for i in range(n):
        for j in range(m):
            M[i][j] = (1/D) * M[i][j]
    return M
'''


##### Useful Functions #####

def graphPower(G, n):
    """
        Computes G^n

        G : graph
        n : int

        return : graph G'
    """
    def NRot(node, path): #node in [0,N] and path = (x for x in [D])
        Nnode = node
        Npath = []
        for i in range(n):
            (Nnode, idx) = G(Nnode,path[i])
            Npath.append(idx)
        return (Nnode,Npath)
    return NRot


def zigzagProduct(G, H):
    """
        Computes G Ⓩ H, the zigzag product between G and H

        G : Rotation map of the graph
        H : Rotation map of the graph

        return : G' : graph
    """
    def Rot(v,a,i,j):
        aa, ii = H(a,i)
        w,bb = G(v,aa)
        b,jj = H(bb,j)
        return ((w,b),jj,ii)
    return Rot


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

