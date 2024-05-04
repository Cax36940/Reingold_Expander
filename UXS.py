"""
This file is not finished
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