from utils_graph import *
import math 
from utils import *
""" Notation N or D instead of [|1, N|] or [|1, D|] """

##### Constant Size Graph #####

def rotH_simple(nb_vertices = 3) :
    """
    Rotation map for cycle on nb_vertices vertices
    Uses pi-permutation (0, 1) -> (1, 0)

    This is a smaller graph for testing purpose

    vertex : int    in [0, nb_vertices]
    edge :   int    in [0, 1]
    nb_vertices : int (default value = 3)
    
    return (vertex', edge') : (int, int)
    
    """
    def rotH(vertex, edge):
        if edge == 0 :
            return ((vertex-1)%nb_vertices, 1)
        return ((vertex+1)%nb_vertices, 0)
    return rotH

def rotH_complete(vertex, edge):
    """
        Rotation map of (D^16, D, 1/2) graph H

        This is the real constant size graph that should be used

        vertex : int[16]    in D^16
        edge : int       in D

        return : (vertex', edge') : (int[16], int[2])
    """

    print(edge)
    q = 67
    d = 31
    x = edge // q
    y = edge % q

    # convert vertex from int in D^16 to int[32] in q^32
    a = int_to_list(vertex, 32, q)
    
    b = [(y * (x**i))%q for i in range(d + 1)]

    c = [(a[i] + b[i])%q for i in range(d + 1)]

    new_vertex = list_to_int(c, 32, q)
    new_edge = x*q + ((-y)%q)
    print(new_edge)
    return (new_vertex, new_edge)


##### Transformation Functions #####

def rotRegularize(G, N):
    """
        Rotation map of N² vertex regular graph based on G

        G :      graph (int[][])
        vertex : int in N x N
        edge :   int
        
        return : (vertex', edge') -> (int, int)
    """
    def rotGreg(vertex, edge):
        if edge == 0:
            return ((vertex // N) * N + (vertex + 1) % N, 1)
        if edge == 1:
            return ((vertex // N) * N + (vertex - 1) % N, 0)
        if G[vertex // N][vertex % N] :
            return ((vertex % N) * N + (vertex // N), edge)

        return vertex, edge
    return rotGreg

def rotGraphPower(G, D, n):
    """
        Computes the rotation map of G^n with D the degree of G

        rotG : graph
        D : int
        n : int

        return : Rotation map of the graph G^n
    """
    def NRot(v, a): #node in [0,N] and path = (x for x in [D])
        w = v
        b = 0
        for i in range(n):
            b *= D
            (w, edgeG) = G(w,(a//D**(n-i-1))%D)
            b += edgeG
        return (w,b)
    return NRot


def rotZigzagProduct(G, H, D, d):
    """
        Computes G Ⓩ H, the zigzag product between G and H (graph on D vertex of degree d)

        rotG : Rotation map of the graph G
        rotH : Rotation map of the graph H
        D : int
        d : int

        return : Rotation map of the graph GoH
    """
    def Rot(vertex, edge):
        v, a = vertex // D, vertex % D
        i, j = edge // d, edge % d
        aa, ii = H(a,i)
        w, bb = G(v,aa)
        b, jj = H(bb,j)
        return (w * D + b, jj * d + ii)
    return Rot


def rotMainTransform(rotG, rotH, D, d, n):
    """
        Computes the rotation map of the n-th power of the zigzag product of G and H

        rotG : Rotation map of the graph G
        rotH : Rotation map of the graph H
        D : int
        d : int
        n : int

        return : Rotation map
    """
    return rotGraphPower(rotZigzagProduct(rotG, rotH, D, d), d*d, n) 


def rotGexp_simple(rotG, rotH, l, D, d, n):
    """
        Computes the rotation map for l steps of Reingold main transformation

        rotG : Rotation map of the graph G
        rotH : Rotation map of the graph H
        l : int
        n : int

        return : Rotation map
    """
    
    rotGexp = rotMainTransform(rotG, rotH, D, d, n)
    
    for i in range(l-1):
        rotGexp = rotMainTransform(rotGexp, rotH, D, d, n)

    return rotGexp


def rotGexp_Reingold(rotG, rotH, D, l):
    """
        Rotation map of (N² x (D^16)^L, D^16, 1/2) expander graph

        This is the algorithm found at the end of Reingold research paper

        dependancies :
            rotH
            rotGreg

        vertex : int  in N² x (D^16)^L
        edge :   int  in D^16
        rotG : (int, int) -> (int, int)
        rotH : (int, int) -> (int, int)
        N : int
        D : int
        l : int

        return : rotGexp : Rotation map
    """

    D16 = D**16

    def rotGexp(vertex, edge):

        # Allocate variables
        v = 0
        a = [[0 for i in range(16)] for k in range(l + 1)]
        I = 0
        j = [0 for i in range(l + 1)]

        # Initialise variables
        for i in range(l-1, -1, -1):
            atmp = vertex % D16
            for k in range(15, -1, -1):
                a[i][k] = atmp % D
                atmp = atmp // D
            vertex = vertex // D16
        v = vertex
        
        for k in range(15, -1, -1):  
            a[l][k] = edge % D
            edge = edge // D
        print(a)
        I = l

        for i in range(1, l + 1):
            j[i] = 1
    
        while True :

            # Step 1
            aI_1 = list_to_int(a[I-1], 16, D)
            aI_1, a[I][j[I]-1] = rotH(aI_1,a[I][j[I]-1])
            a[I-1] = int_to_list(aI_1, 16, D)

            # Step 2
            if j[I]%2 == 1:
                if I == 1 :
                    # Convert a[0] from list to int
                    a0 = list_to_int(a[0], 16, D)
                    
                    v, a0 = rotG(v, a0)

                    a[0] = int_to_list(a0, 16, D)

                    j[I] += 1

                # Step 3
                elif I > 1 :
                    j[I-1] = 1
                    I -= 1
        
            # Step 4
            elif j[I] == 16 :
                a[I].reverse()
        
            # Step 5
                if I == l :
                    break
        
            # Step 6
                elif I < l :
                    j[I + 1] += 1
                    I += 1
            else:
                j[I] += 1
        

        new_vertex = v
        for i in range(l):
            new_vertex *= D16
            new_vertex += list_to_int(a[i], 16, D)
        print(a)
        new_edge = list_to_int(a[l], 16, D)

        return new_vertex, new_edge
    return rotGexp

##### USTCON #####

def Aexp(G, s, t):
    """
        Return true if s and t are connected in Gexp

        global const D
        global const N

        computes Delta

        dependancies :
            rotGexp

        s : [(int, int), int[L][16]]   in N² x (D^16)^L
        t : [(int, int), int[L][16]]   in N² x (D^16)^L
        G : a D-regular graph with N vertices

        return : con : boolean
    """
    rotG = adjacencyMatrix_to_rotMap(G)
    rotH = rotH_complete

    lmax = math.log(2,N) 
    
    def Aexprec(vertex,edge,l):
        if l >= lmax:
            return edge 
        else:
            (Nvertex,Nedge) = rotGexp_Reingold(vertex,edge,rotG,rotH,N,D)
            Aexprec(Nvertex, Nedge, l+1)
    l = 0
    for i in range(D):
        if (Aexprec(s,i,l) == t):
            return True
        l = 0
    return False
        



def Acon(G, s, t):
    """
        Return true if s and t are connected in Greg, implying connexity in G

        global const D
        global const N

        dependancies :
            Aexp

        s : (int, int)  in N x N
        t : (int, int)  in N x N

        return : con : boolean
    """
    lmax = math.log(2,N)
    pS = (s, [1 for i in range(lmax)])
    pT = (t, [1 for i in range(lmax)])
    return Aexp(G,pS,pT)
