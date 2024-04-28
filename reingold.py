import numpy as np;

def rot(G): #G : matrice d'adjacence, on suppose qu'il n'existe pas d'arêtes parallèles et que G est régulier
    def fct(v,i):
        sum = 0
        N = len(G)
        k = 0 
        while (k < N):
            sum += G[v][k] 
            if (i == sum):
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




##### Useful Functions #####

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


def secondEV(G): #done alrdy
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

##### USTCON #####

# Notation N or D instead of [|1, N|] or [|1, D|]

def rotH(vertex, edge):
    """
        Rotation map of (D^16, D^2, 1/2) graph H

        global const D

        vertex : int[16]    in D^16
        edge : int[2]       in D^2

        return : (vertex', edge') : (int[16], int[2])
    """
    return 0


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



def Ae2p(vertex, edge, N):
    """
        Compute the path, in a graph with pi-permutation, taken while walking on an edge in its expander transformation

        global const D

        dependancies :
            rotH
            piPerm

        vertex : [(int, int), int[L][16]]   in N² x (D^16)^L        #(int, int) is not necessary
        edge : int[16]                      in D^16
        N : int

        (computations are similar to the ones of rotGexp)

        return : path : int[][16]     in (D^16)*

    """
    return 0


def Asrch(path, N): #not really what was defined in the paper
    """
        Compute path in a graph with pi-permutation, taken while walking on L an edges in its expander transformation


        dependancies :
            Ae2p


        path : path : int[L][16]     in (D^16)^L
        N : int

        calls Ae2p for each edge of the path


        return : path : int[][16]     in (D^16)*
    """




def UTS(N, Dmax):
    """
        Compute UTS on N*Dmax vertex D^16 regular graph with pi permutation labelling

        N : int
        Dmax : int

        dependancies :
            Asrch

        calls Asrch on every path of length L of edges in D^16

        return : UTS : int[][16] in D^16
    """



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
