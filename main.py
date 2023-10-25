from fd_main import *

def main():
    '''Main driver for the Filter Diagonalization package'''
    
    print("Running FILTER DIAGONALIZATION (Python version!)\n")
    print(datetime.now(), '\n')
    print('*****************************************************************\n')
    Ti = time.time()
    print("1. INITIALIZING \n")
    ti = time.time()
    potl, ksqr, el, par, ist, atm_ids, coords, atm = init()
    tf = time.time()
    print(f"Done initializing job, wall run time {tf - ti:.4f} (sec)\n")
    
    
    print('\n*****************************************************************\n')
    print("2. CALCULATING ENERGY RANGE OF HAMILTONIAN\n")
    ti = time.time()
    get_energy_range(potl,ksqr,par,ist)
    tf = time.time()
    print(f"Done calculating energy range, wall run time {tf - ti:.4f} (sec)\n")
    
    
    print('\n*****************************************************************\n')
    print("3. FILTERING STATES\n")
    ti = time.time()
    
    print("Setting parameters for the Newton interpolation...")
    par['dt'] = (ist['nCheby'] / (2.5 * par['dE']))**2
    print(f"nCheby = {ist['nCheby']} dE = {par['dE']:.5f} dt = {par['dt']:.5f}")
    
    # Generate the Chebyshev expansion coefficients
    an, zn = coefficient(par, ist) 
    
    psitot = []
    for jns in range(ist['ns']):
        seed = int(time.time())+jns # random number generator seed based on current system time
        #print(f"seed = {seed}\n")
        psi = init_psi(par, ist, seed)
        for jms in range(ist['ms']):
            psitot.append(psi)
    
    psitot = np.array(psitot, dtype=np.float64)
    for jns in range(ist['ns']):
        psitot[jns*ist['ms']:jns*ist['ms']+ist['ms']] = filtering(psitot[jns*ist['ms']:jns*ist['ms']+ist['ms']], potl, ksqr, el, an, zn, ist, par, jns)
    
    # Write cube files of the filtered states
    for i in range(ist['ns']):
        for j in range(ist['ms']):
            write_cube(psitot[i*ist['ms']+j], f"psi-filt-{i}-{j}", atm_ids , coords, atm, par, ist)
    
    tf = time.time()
    print(f"Done filtering states, wall run time {tf - ti:.4f} (sec)\n")
    
    
    print('\n*****************************************************************\n')
    print("3. ORTHOGONALIZING STATES\n")
    ti = time.time()
    psitot = psitot.reshape((ist['ns']*ist['ms'], ist['ngrid'])).T # Cast the grid-based states into a matrix
    psitot = orth(psitot).T # orthogonalize the states
    psitot = psitot.reshape((psitot.shape[0], ist['nx'], ist['ny'], ist['nz'])) # Recast back onto the grid
    ist['mstot'] = psitot.shape[0]
    
    print(f"Number of orthogonal states, mstot = {ist['mstot']}\n")
    
    psitot = normalize_all(psitot, par['dv'], ist['mstot'])
    
    # Write cube files of the orthogonalized states
    for i in range(ist['mstot']):
            write_cube(psitot[i], f"psi-ortho-{i}", atm_ids , coords, atm, par, ist)
    
    tf = time.time()
    print(f"Done orthogonalizing states, wall run time {tf - ti:.4f} (sec)\n")
    
    
    print('\n*****************************************************************\n')
    print("4. DIAGONALIZING STATES\n")
    ti = time.time()
    
    psitot, evals = Hmatreal(psitot, potl, ksqr, ist, par)
    psitot = normalize_all(psitot, par['dv'], ist['mstot'])
    
    # Write cube files of the orthogonalized states
    for i in range(ist['mstot']):
        write_cube(psitot[i], f"psi-diag-{i}", atm_ids , coords, atm, par, ist)
    
    tf = time.time()
    print(f"Done diagonalizing states, wall run time {tf - ti:.4f} (sec)\n")
    
    
    print('\n*****************************************************************\n')
    print("5. CALCULATING EIGENSTATE VARIANCE\n")
    ti = time.time()
    sige = calc_sigma_E(psitot, potl, ksqr, ist, par);
    np.savetxt("eval.dat", np.array([evals, sige]).T , fmt="% .8f % .4e", delimiter="\n") # save the energies
    tf = time.time()
    print(f"Done calc sigma E, wall run time {tf - ti:.4f} (sec)\n")
    Tf = time.time()
    print('\n*****************************************************************')
    print("Done with program FILTER DIAGONALIZATION (Python version!)\n")
    print(f"Wall run time {timedelta(seconds=(round(Tf - Ti)))}") 
    print(datetime.now(), '\n')
    print('*****************************************************************')

main()