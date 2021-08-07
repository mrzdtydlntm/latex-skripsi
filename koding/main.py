from huckel import huckel_2
import numpy as np
from pprint import pprint

def main():
    list_egap = []
    N = np.arange(6,101,2)
    for X in N:
        H, eigval_schur, Egap = huckel_2(X)
        list_egap.append(Egap)

    file = open("spectra-graphene.txt", "w")
    zipped = np.vstack((x,y_list)).T
    np.savetxt(file, zipped, delimiter=' ')

if __name__ == "__main__":
    main()
