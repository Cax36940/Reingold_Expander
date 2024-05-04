from utils_graph import *

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
            return (vertex-1)%nb_vertices, 1
        return (vertex+1)%nb_vertices, 0
    return rotH

def rotH_complete(vertex, edge):
    """
        Rotation map of (D^16, D, 1/2) graph H

        This is the real constant size graph that should be used

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


##### Transformation Functions #####

def rotRegularize(G):
    """
        Rotation map of N² vertex regular graph based on G

        G :      graph (int[][])
        vertex : int in N x N
        edge :   int
        
        return : (vertex', edge') -> (int, int)
    """
    def rotGreg(vertex, edge):
        N = len(G)
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


def rotGexp_simple(rotG, rotH, l):
    """
        Computes the rotation map for l steps of Reingold main transformation

        rotG : Rotation map of the graph G
        rotH : Rotation map of the graph H
        l : int

        return : Rotation map
    """
    # This algorithm can't be executed with the true value of l
    
    rotGexp = rotMainTransform(rotG, rotH, 8)
    
    for i in range(l-1):
        rotGexp = rotMainTransform(rotG, rotH)

    return rotGexp


def rotGexp_Reingold(vertex, edge, rotG, rotH, N, D):
    """
        Rotation map of (N² x (D^16)^L, D^16, 1/2) expander graph

        This is the algorithm found at the end of Reingold research paper

        global const D
        global const N

        dependancies :
            rotH
            rotGreg

        vertex : [(int, int), int[L][16]]   in N² x (D^16)^L
        edge : int[16]                      in D^16


        return : (vertex', edge') : ([(int, int), int[L][16]], int[16])
    """

    def rotGexp(vertex, edge):
        # This algorithm can't be executed with the true value of l
        # l = maxPower(N, D)

        l = 5

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
        
        return vertex, edge
    return rotGexp

##### USTCON #####

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

        dependancies :
            Aexp

        s : (int, int)  in N x N
        t : (int, int)  in N x N

        return : con : boolean
    """
