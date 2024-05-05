import networkx as nx
import numpy as np

##### Rotation Map and Adjacency Matrix #####

def adjacencyMatrix_to_rotMap(G):
    """
        Returns the rotation map function for regular graph G given its adjacency matrix
        Does not works with parallel edges
        Fix a labelling that depends on the adjacency matrix


        G : graph (int[][])

        return rotG : (int, int) -> (int, int)
    
    """
    def fct(v,i):
        sum = 0
        N = len(G)
        k = 0 
        w = 0
        while (k < N):
            sum += G[v][k] 
            if (sum > i):
                w = k
                break
            k += 1
        k = 0 
        sum = 0
        j = 0
        while (k < N):
            sum += G[w][k]
            if (v == k):
                j = sum - 1
                break
            k += 1
        return (w,int(j))
    return fct



def rotMap_to_adjacencyMatrix(f, N, D):
    """
        Returns the adjacency matrix of a graph, given its rotation map

        f : (int, int) -> (int, int)
        N : int
        D : int

        return : AMatrix : int[][]
    """
    AMatrix = np.array([[0]*N]*N)
    for u in range(N):
        for i in range(D):
            (v,j) = f(u,i)
            AMatrix[u][v] += 1
    return AMatrix


##### Graph Generation #####

def random_connected_graph(N):
    """
    Returns a random connected graph on N vertices

    N : int

    return G : graph (int[][])
    """
    # Generate a random D-regular graph
    G = nx.gnp_random_graph(N, 0.5)
    while not nx.is_connected(G): # Verify if it is connected
        G = nx.gnp_random_graph(N, 0.5)
    return nx.to_numpy_array(G)

def random_non_connected_graph(N):
    """
    Returns a random non-connected graph on N vertices

    N : int

    return G : graph (int[][])
    """
    # Generate a random D-regular graph
    G = nx.gnp_random_graph(N, 0.1)
    while nx.is_connected(G): # Verify if it is connected
        G = nx.gnp_random_graph(N, 0.1)
    return nx.to_numpy_array(G)


def random_regular_graph(N, D):
    """
    Returns a random D-regular connected non-bipartite graph on N vertices

    N : int
    D : int

    return G : graph (int[][])
    """
    # Generate a random D-regular graph
    G = nx.random_regular_graph(D, N)
    while not nx.is_connected(G) or nx.is_bipartite(G): # Verify if it is connected
        G = nx.random_regular_graph(D, N)
    return nx.to_numpy_array(G)

def cycle_graph(N):
    """
    Returns cycle graph on N vertices

    N : int

    return G : graph (int[][])
    """
    G = np.zeros((N,N))
    for i in range(N):
        G[i][(i+1)%N] = 1
        G[(i+1)%N][i] = 1
    return G

def complete_graph(N):
    """
    Returns complete graph on N vertices

    N : int

    return G : graph (int[][])
    """
    G = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j :
                G[i][j] = 0
            else :
                G[i][j] = 1
    return G


##### Other #####

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