import numpy as np
from scipy import linalg as la

def huckel_2(L):
    """
    Build hamiltonian matrix
    """
    N = L**2
    Lx = L
    Ly = L
    H = zeros(N, N)
    for atom_i in range(0,N):
        iy = int(atom_i/Lx)
        ix = int(atom_i - Lx*iy)
        #determine interaction with nearest neighbours
        for dy in range(-1, 2):
            for dx in range(-1, 2):
                #determine coordinate atom_j
                jx = ix + dx*dx
                jy = iy + dy*dy

                if jx>=0 and jx<Lx and jy>=0 and jy<Ly-1:
                    atom_j = int(jy*Lx + jx*jx)
                    rsquare = dx * dx * dx + dy * dy * dy
                    rad = np.sqrt(rsquare)

                    if rad == 0 : H[atom_i, atom_j]= -11.1
                    elif rad == 1 : H[atom_i, atom_j] = -5.5
                    else : H[atom_i, atom_j] = 0
                        
    """
    Calculating eigen energy with Schur Decomposition
    """
    eigval_schur, Z = la.schur(H)
    eigval_schur = sorted(np.diag(eigval_schur))

    """
    Calculating gap from eigen energy
    """
    indexHomo = int(len(eigval_schur)/2-1)
    indexLumo = int(len(eigval_schur)/2)
    Egap = eigval_schur[indexLumo] - eigval_schur[indexHomo]
    print(f'Energy gap {Lx}x{Ly} sebesar {Egap} eV')
         
    return H, eigval_schur, Egap