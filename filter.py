from fd_filter import *

# Filter functions

def filtering(psitot, potl, ksqr, el, an, zn, ist, par, jns):
    psi = psitot[0].copy() # Get the first
    print("filtering\t%.3e"%(np.min(psi)))
    phi = np.zeros((ist['nx'], ist['ny'], ist['nz']))
    psitot = filter_states(psi, phi, psitot, potl, ksqr, an, zn, ist, par)
    psitot = normalize_all(psitot, par['dv'], ist['ms'])
    
    evals = energy_all(psi, phi, psitot, potl, ksqr, ist, par)
    
    eval_file = open(f"eval-filt-{jns}.dat", "w")
    for jms in range(ist['ms']):
        eval_file.write(f"{jms} {evals[jms]} {el[jms]}\n")
    eval_file.close()

    return psitot
    
    
def filter_states(psin, psim1, psi0, potl, ksqr, an, zn, ist, par):
    
    for ie in range(ist['ms']):
        ncie = ist['nCheby']*ie
        psi0[ie] = an[ncie].real * psin
    #print("psi0\t%.3e"%(np.min(psi0)))
    for j in range(ist['nCheby']):
        psim1 = psin.copy()
        
        psim1, psin = hnorm(psim1, psin, potl, ksqr, zn[j-1], ist, par)
        for ie in range(ist['ms']):
            ncie = ist['nCheby']*ie
            psi0[ie] += an[ncie+j].real * psin.real
            #if ie == 0:
                #print("%d\t%.3e"%(j, np.min(psi0[ie])))
    print("min\t%.3e"%(np.min(psi0[0])))
    return psi0

def hnorm(psim1, psin, potl, ksqr, zm1, ist, par):
    psim1, psin = hamiltonian_norm(psim1, psin, potl, ksqr)
    #print("hnorm\t%.3e"%(np.min(psin)))
    psin = par['dE_1'] * psin.real - (2.0 + zm1 + par['Vmin'] * par['dE_1']) * psim1.real\
        + (par['dE_1'] * psin.imag - (2.0 + zm1 + par['Vmin'] * par['dE_1']) * psim1.imag)*1j

    return psim1, psin

def hamiltonian_norm(psi, phi, potl, ksqr):
    
    phi = kinetic(phi, ksqr)
    phi += np.multiply(potl, psi)
    
    return psi, phi
    
