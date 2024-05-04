from matrix_op import *
import USTCON as st
from utils import *
from utils_graph import *



def test_zig_zag_product():
    """
        This example come the explanation of the zigzag product on Wikipedia
        Wikipedia : https://fr.wikipedia.org/wiki/Produit_zig-zag_de_graphes
    """
    G = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            if i != j and abs(i-j) != 3:
                G[i][j] = 1
    
    H = np.array([[0.,1.,0.,0.],[1.,0.,0.,0.],[0.,0.,0.,1.],[0.,0.,1.,0.]])

    GoH = zigzagProductMatrix(G,H)

    # Draw the graphs
    plt.figure(figsize=(15, 5))

    plt.subplot(131)
    G = nx.from_numpy_array(G)
    nx.draw(G, nx.circular_layout(G))
    plt.title('Graph G')

    plt.axvline(x=1.2, color='gray', linestyle='--')

    plt.subplot(132)
    H = nx.from_numpy_array(H)
    nx.draw(H, nx.circular_layout(H))
    plt.title('Graph H')

    plt.axvline(x=1.2, color='gray', linestyle='--')

    plt.subplot(133)
    GoH = nx.from_numpy_array(GoH)
    nx.draw(GoH, nx.circular_layout(GoH))
    plt.title('Zigzag Product Graph (G ◦ H)')

    plt.tight_layout()
    plt.show()



# def test_regularize():

if __name__ == "__main__":
    test_zig_zag_product()

