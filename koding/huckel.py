import numpy as np
from scipy import linalg as la

def huckel(Lx, Ly):
    """
    Build hamiltonian matrix
    """
    N = Lx*Ly
    H = zeros(N, N)
    for atom_i in range(0,N):
        iy = int(atom_i/Lx)
        ix = int(atom_i - Lx*iy)
        #determine interaction with nearest neighbours
        for dy in range(-1, 2):
            for dx in range(-1, 2):
                #determine coordinate atom_j
                jx = ix + dx
                jy = iy + dy

                if jx>=0 and jx<Lx and jy>=0 and jy<Ly:
                    atom_j = int(jy*Lx + jx)
                    rsquare = dx * dx + dy * dy
                    rad = np.sqrt(rsquare)

                    if rad == 0 : H[atom_i, atom_j]= -11.1
                    elif rad == 1 : H[atom_i, atom_j] = -5.5
                    else : H[atom_i, atom_j] = 0

    """
    Calculating eigen energy with Schur Decomposition
    """
    # eigval_schur, Z = la.schur(H)
    # eigval_schur = sorted(np.diag(eigval_schur))

    """
    Calculating eigen energy and eigen number
    """
    eigval, eignum = la.eig(H)
    eigval = sorted(eigval) #sorting secara ascending

    """
    Calculating gap from eigen energy
    """
    indexHomo = int(len(eigval)/2-1) #nilai terkecil
    indexLumo = int(len(eigval)/2) #nilai terbesar
    Egap = eigval[indexLumo] - eigval[indexHomo] #dalam satuan elektronvolt

    """
    Calculating deltaE
    """
    dE = []
    for indexEnergyLumo in range(0, indexLumo):
        for indexEnergyHomo in range (0,indexHomo+1):
            dE.append(eigval[indexEnergyLumo + indexLumo] - eigval[indexHomo - indexEnergyHomo])

    """
    Calculating dipole moment
    """
    energyNumber = int(len(eignum))
    M = []
    for indexEnergyLumo in range(0, indexLumo):
        for indexEnergyHomo in range (0, indexHomo+1):
            mtot = 0
            for indexC in range(0, energyNumber):
                mtot += eignum[indexC][indexEnergyLumo + indexLumo] * eignum[indexC][indexHomo - indexEnergyHomo]
            M.append(mtot)
    Ms = sorted(M, reverse=True)

    """
    Calculating Absorbance Spectra
    """
    def dipole_moment(Ndiff, Gamma, E):
        result = 0.0
        for atom_i in range (0, Ndiff):
            disc = E - dE[atom_i]
            disc /= Gamma # gamma = width nya
            result+=(Ms[atom_i]/(1+disc**2))
        return result
    Npoints = 3000
    x = np.linspace (0, 4, Npoints)
    y_list = []
    for i in range(0,Npoints):
        y = dipole_moment(N, 0.01, x[i])
        y_list.append(y.real)

    return H, eigval, eignum, Egap, dE, Ms, x, y_list

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
                jx = ix + dx
                jy = iy + dy

                if jx>=0 and jx<Lx and jy>=0 and jy<Ly:
                    atom_j = int(jy*Lx + jx)
                    rsquare = dx * dx + dy * dy
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