import numpy as np
import matplotlib.pyplot as plt
from math import sqrt ##### Useful Functions #####
from utils_graph import *

def regularize(G):
    N = len(G[0])
    Greg = np.zeros((N*N, N*N))
    cycle = cycle_graph(N)

    # Replace every vertex by a cycle of length N
    for v in range(N):
        for i in range(N):
            for j in range(N):
                Greg[v*N + i][v*N + j] = cycle[i][j]
    
    # There is an edge between (v,w) and (w,v) if there is one in G
    # Else self loop
    for i in range(N):
        for j in range(N):
            if G[i][j] :
                Greg[i*N + j][j*N + i] = 1
            else :
                Greg[i*N + j][i*N + j] = 1
    return Greg

def graphPowerMatrix(G, n):
    """
        Computes G^n

        G : graph (int[][])
        n : int

        return : graph G'
    """
    return np.linalg.matrix_power(G, n)


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


def mainTransformMatrix(G, H, n):
    return graphPowerMatrix(zigzagProductMatrix(G, H), n)







#TODO faire un graph, sa transformation en TS en XS
# print(maxPower(8, 67*67))
# Ae2p([[1 for i in range(16)] for k in range(6)], [1 for i in range(16)], 8)

