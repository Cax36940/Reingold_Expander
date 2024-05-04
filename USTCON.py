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
            (v,j) = f(u,i+1)
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


