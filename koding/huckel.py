import numpy as np
from scipy import linalg as la

def huckel(Lx, Ly):
    """
    Build hamiltonian matrix
    """
    N = Lx*Ly
    H = np.zeros(shape = (N, N))
    for atom_i in range(0,N):
        alpha = -10.7
        beta = -1.58
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

                    if rad == 0: 
                        H[atom_i][atom_j]= alpha
                    elif rad == 1:
                        if Lx%2 == 0 and abs(atom_i-atom_j)==Lx and atom_i%2 == 1: #jika panjang x genap, dan atomi-atomj = panjang x dan atom i ganjil
                            H[atom_i][atom_j]=0
                        elif Lx%2 == 1 and iy%2 == 0 and abs(atom_i-atom_j)==Lx and atom_i%2==1: #jika panjang x ganjil, iy genap, absolut dari atomi-atomj = panjang x dan atom i ganjil
                            H[atom_i][atom_j]=0
                        elif Lx%2 == 1 and iy%2 == 1 and abs(atom_i-atom_j)==Lx and atom_i%2==0: #jika panjang x ganjil, iy ganjil, atomi-atomj = panjang x dan atom i genap
                            H[atom_i][atom_j]=0    
                        else : 
                            H[atom_i][atom_j] = beta
                        # if abs(atom_i-atom_j)==Lx and atom_i%2==1: 
                            # H[atom_i][atom_j]=0
                        # else: 
                            # H[atom_i][atom_j] = -1.58
                    else: 
                        H[atom_i][atom_j] = 0

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
        Npoints = 1000
        x = np.linspace (0, 1, Npoints)
        y_list = []
        for i in range(0,Npoints):
            y = dipole_moment(N, 0.01, x[i])
            y_list.append(y.real)

    return H, eigval, eignum, Egap, dE, Ms, x, y_list

def huckel_2(L):
    """
    Build hamiltonian matrix
    """
    N = L*2
    H = np.zeros(shape = (N, N))
    for atom_i in range(0,N):
        alpha = -10.7
        beta = -1.58
        iy = int(atom_i/L)
        ix = int(atom_i - L*iy)
        #determine interaction with nearest neighbours
        for dy in range(-1, 2):
            for dx in range(-1, 2):
                #determine coordinate atom_j
                jx = ix + dx
                jy = iy + dy

                if jx>=0 and jx<L and jy>=0 and jy<2:
                    atom_j = int(jy*L + jx)
                    rsquare = dx * dx + dy * dy
                    rad = np.sqrt(rsquare)

                    if rad == 0: 
                        H[atom_i][atom_j]= alpha
                    elif rad == 1:
                        if L%2 == 0 and abs(atom_i-atom_j)==L and atom_i%2 == 1: #jika panjang x genap, dan atomi-atomj = panjang x dan atom i ganjil
                            H[atom_i][atom_j]=0
                        elif L%2 == 1 and iy%2 == 0 and abs(atom_i-atom_j)==L and atom_i%2==1: #jika panjang x ganjil, iy genap, absolut dari atomi-atomj = panjang x dan atom i ganjil
                            H[atom_i][atom_j]=0
                        elif L%2 == 1 and iy%2 == 1 and abs(atom_i-atom_j)==L and atom_i%2==0: #jika panjang x ganjil, iy ganjil, atomi-atomj = panjang x dan atom i genap
                            H[atom_i][atom_j]=0    
                        else : 
                            H[atom_i][atom_j] = beta
                        # if abs(atom_i-atom_j)==L and atom_i%2==1: 
                            # H[atom_i][atom_j]=0
                        # else: 
                            # H[atom_i][atom_j] = -1.58
                    else: 
                        H[atom_i][atom_j] = 0
                        
        """
        Calculating eigen energy and eigen number
        """
        eigval, eignum = la.eig(H)
        eigval = sorted(eigval)

        """
        Calculating gap from eigen energy
        """
        indexHomo = int(len(eigval)/2-1)
        indexLumo = int(len(eigval)/2)
        Egap = eigval[indexLumo] - eigval[indexHomo]

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
        l = []
        A = []
        for i in range(len(Ms)):
            for m in range(1,1000):
                el = 0.1 + m*0.0012
                l.append(el)
                a = 0
                a+=abs(Ms[i]**2/(dE[i]-(1.2/el-0.1*i)))
                A.append(a)
        
    return H, eigval, eignum, Egap, dE, Ms, l, A