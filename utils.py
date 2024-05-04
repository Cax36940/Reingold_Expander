import numpy as np

def deg(G):
    """
    Returns the degree of G

    G : graph (int[][])
    
    return : D : int
    """
    return int(sum(G[0]))


def pathLength(N):
    """
        Computes the max length of a log path in a (N, D, 1/2) expander

        N : int

        return : Delta : int
    """
    return 4 * int(np.ceil(N))

def maxPower(N, D):
    """
        Computes the max number of iteration for the Reingold Algorithm main transformation

        N : int
        D : int

        return : l : int
    """ 
    return int(np.ceil(-np.log2(-np.log2(1-1/(D*N*N)))))

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
