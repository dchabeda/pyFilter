from fd_energy import *

# Energy functions
def energy(psi,phi,potl,ksqr,dv):
    phi = psi.copy()
    phi, psi = hamiltonian(phi, psi, potl, ksqr)
    ene = np.sum(np.multiply(psi.conj(), phi))*dv # expectation value of energy
    return ene.real
    

def energy_all(psi, phi, psi0, potl, ksqr, ist, par):
    evals = []
    
    for ie in range(ist['ms']):
        psi = psi0[ie]
        phi = psi.copy()
        phi, psi = hamiltonian(phi, psi, potl, ksqr)
        ene = np.sum(np.multiply(psi.conj(), phi))*par['dv']
        evals.append(ene.real)
        
    return evals

    
def get_energy_range(potl,ksqr,par,ist):
    #seed = int(time.time())
    seed = 874917403
    
    psi = init_psi(par, ist, seed)
    phi = psi.copy # Save a copy of the original wavefunction for expectation value
    
    ene_old = 0; ene = ene_old + 0.1
    ene_file = open("ene-init.dat", "w")
    for i in range(500):
        if np.abs((ene - ene_old)/ene) > 1.0e-4:
            psi, phi = hamiltonian(psi, phi, potl,ksqr)
            psi *= (norma := norm(psi, par['dv']))
            ene_old = ene
            ene = energy(psi, phi, potl, ksqr, par['dv'])
            ene_file.write(f"{i:2d} {ene_old:10.8f} {ene:10.8f} {norma:10.8f}\n")
            
    
    ene_file.close()
    par['dE'] = 1.1*(ene - par['Vmin'])
    par['dE_1'] = 4.0 / par['dE']
    
    print(f"Energy range calculated, dE = {par['dE']:.4f} Hartree {par['dE']*27.211:.5f} eV\n")


def calc_sigma_E(psitot, potl, ksqr, ist, par):
    sige = np.zeros(ist['mstot'], dtype=np.float64)
    
    for ims in range(ist['mstot']):
        psi = psitot[ims].copy()
        phi = psi.copy()
        phi, psi = hamiltonian(phi, psi, potl, ksqr)
        
        ene = np.sum(np.multiply(psitot[ims], phi).real)*par['dv'] # Compute <E>_i = <psi_i|H|psi_i>
        
        psi = phi.copy()
        phi, psi = hamiltonian(phi, psi, potl, ksqr)
        
        ene2 = np.sum(np.multiply(psitot[ims], phi).real)*par['dv']
        sige[ims] = np.sqrt(np.abs(ene2 - ene**2))
    
    return sige
    
def hamiltonian(phi, psi, potl, ksqr):
    psi = phi.copy()
    phi = kinetic(phi, ksqr) # Compute kinetic energy
    phi += np.multiply(potl, psi) # K.E. + potential is the Hamiltonian
    return phi, psi

    
def kinetic(psi, ksqr):
    psifft = fftn(psi, s=[psi.shape[0], psi.shape[1], psi.shape[2]]) # Forward FT the wavefunction
    np.multiply(ksqr, psifft, out=psifft) # Multiply by k^2
    psi = ifftn(psifft,s=[psifft.shape[0], psifft.shape[1], psifft.shape[2]]) # Inverse FT to get K.E.|Psi>
    
    return psi
    
