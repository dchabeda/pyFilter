from fd_Hmat import *

def Hmatreal(psitot, potl, ksqr, ist, par):
    '''Calculate and diagonalize the Hamiltonian matrix'''
    
    # Calculate <psi_j|H|psi_i>
    H = np.zeros((ist['mstot'], ist['mstot']))
    hmat_file = open("hmat.dat", "w")
    #hmat_file.write("Calculating the H matrix:\n")
    for ims in range(ist['mstot']):
        # Calculate phi_i = |H|psi_i>
        psi = psitot[ims].copy()
        phi = psi.copy()
        phi, psi = hamiltonian(phi, psi, potl, ksqr)
        for jms in range(ims, ist['mstot']):
            # Calculate <psi_j|H|psi_i>
            H[jms,ims] = H[ims,jms] = np.sum(np.multiply(psitot[jms].conj(), phi)).real*par['dv']
    
    for i in range(ist['mstot']):
        for j in range(ist['mstot']):
            hmat_file.write(f"{H[i,j]: .6f} ")
        hmat_file.write("\n")
        
    evals, H = eigh(H) # Diagonalize the matrix
    
    tpsi = np.zeros(ist['mstot'])
    for jx in range(ist['nx']):
        for jy in range(ist['ny']):
            for jz in range(ist['nz']):
                for jms in range(ist['mstot']):
                    tpsi[jms] = psitot[jms][jx,jy,jz]
                    psitot[jms][jx,jy,jz] = 0
                for jms in range(ist['mstot']):
                    sum_tot = 0.0
                    for ims in range(ist['mstot']):
                        sum_tot += H[ims,jms]*tpsi[ims]
                    psitot[jms][jx,jy,jz] = sum_tot

    return psitot, evals