import numpy as np

def norm(psi, dv):
    integral = np.sum(np.abs(psi)**2) * dv
    return 1/np.sqrt(integral)

def normalize_all(psitot, dv, ms):
    for ie in range(ms):
        norma = np.sqrt((np.sum(np.abs(psitot[ie])**2) * dv))
        if norma != 0:
            #print("norma: ", norma)
            psitot[ie] /= norma
        else:
            print(f"norm is 0! values of psi ~ {psitot[ie][2,2,2]}")
    return psitot

