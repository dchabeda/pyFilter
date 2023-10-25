import numpy as np
import matplotlib.pyplot as plt

def filt_func(x, dt):
    # Define the filter function
    return np.sqrt(dt / np.pi) * np.exp(- dt * x**2)


def coefficient(par, ist):
    '''Compute the coefficients of the Newton interpolation for a Chebyshev expansion of the filter function (a Gaussian)'''
    print("\nDetermining Chebyshev coefficients...")
    
    # Get the energy grid
    El = np.loadtxt("Egrid.dat")
    #print('The energy grid: ', El)

    # Define parameters for readability
    nc = ist['nCheby']
    ms = ist['ms']
    dt = par['dt']
    res = 1.0; Smin = -2.0; Smax = 2.0
    scale = (Smax - Smin) / par['dE']
    
    # Allocate arrays for the coefficients
    samp = np.zeros(nc, dtype=np.float64)
    an = np.zeros(nc*ms, dtype=complex)
    
    samploc, rho = samp_points_ashkenazy(Smin, Smax, nc)
    
    samp = samploc.real
    
    for ie in range(ms):
        x = (samp[0]+2.0)/scale + par['Vmin']
        an[nc*ie+0] += filt_func(x - El[ie], dt)
        an[nc*ie+0] += 0.0j
    
        x = (samp[1]+2.0)/scale + par['Vmin']
        an[nc*ie+1] += (filt_func(x - El[ie], dt) - an[nc*ie].real)/(samp[1] - samp[0])
        an[nc*ie+1] += (-an[nc*ie+0].imag)/(samp[1] - samp[0])*1j
        
        for j in range(2,nc):
            res = 1.0; sumre = sumim = 0.0
            for i in range(1,j):
                res *= (samp[j] - samp[i-1])
                sumre += an[nc*ie+i].real * res
                sumim += an[nc*ie+i].imag * res
            
            res *= (samp[j] - samp[j-1])
            x = (samp[j]+2.0)/scale + par['Vmin']
            an[nc*ie+j] += (filt_func(x-El[ie], dt) - an[nc*ie+0].real - sumre) / res
            an[nc*ie+j] += ((-an[nc*ie+0].imag - sumim) / res)*1j

    coeff_file = open("coeff.dat", "w")
    for i in range(nc*ms):
        coeff_file.write(f"{i} {an[i].real:8.6g} {an[i].imag:8.6g}\n")
    
    check_function(an,samploc,ist,par,El);
    print("\tdone.\n")
    return an, samploc

def samp_points_ashkenazy(Smin, Smax, nc):
    print("\tDetermining sample points...")
    
    # Setting initialization parameters
    Srange = 0.0
    nc3 = nc*32
    minim = Smin * Srange
    maxim = Smax * Srange
    imfrac = int(nc3 / 8)
    dsre = (Smax - Smin)/(nc3 - 1 - imfrac)
    dsim = (maxim - minim)/(imfrac - 1)
    
    # Allocating arrays for sample points
    point = np.zeros(nc, dtype=complex)
    veca = np.zeros(nc3,dtype=np.float64)
    samp = np.zeros(nc3, dtype=complex)
    
    for j in range(nc3 - imfrac):
        samp[j] += Smin + j*dsre
    for j in range(nc3 - imfrac, nc3):
        samp[j] += (minim + (j - nc3 + imfrac) * dsim)*1j
    
    jrnd = 0
    point[0] = samp[jrnd]
    
    for k in range(nc3):
        Sdel = (samp[k].real - point[0].real)**2 + (samp[k].imag - point[0].imag)**2
        if Sdel < 1.0e-10: veca[k] = -1.0e30
        else: veca[k] = np.log(Sdel)
        
    for j in range(1,nc):
        kmax = 0; fkmax = -2.0e30
        for k in range(nc3):
            if veca[k] > fkmax:
                fkmax = veca[k]
                kmax = k
        point[j] = samp[kmax]
        for k in range(nc3):
            Sdel = (samp[k].real - point[j].real)**2 + (samp[k].imag - point[j].imag)**2
            if Sdel < 1.0e-10: veca[k] = -1.0e30
            else: veca[k] += np.log(Sdel)
    
    fk = 1.0
    for j in range(nc):
        fk *= point[j].real**2 + point[j].imag**2
    
    fk = np.sqrt(fk)
    fk = np.power(fk, 1.0/nc)
    
    zn_file = open("zn.dat", "w")
    for i in range(nc):
        zn_file.write(f"{i} {point[i].real:.8g} {point[i].imag:.8g}\n")
    zn_file.close()
    
    return point, fk
    

    
def check_function(an, samp, ist, par, El):
    dx = 0.01
    x_vals = []; y_vals = []
    func_file = open("func.dat", "w")
    x = -2.0 + 0.0*1j
    fig, ax = plt.subplots(figsize=(9,6))
    for el in El:
        while x.real < 2.0:
            xn = 1.0 + 0.0*1j
            f = an[0]
            for j in range(1,ist['nCheby']):
                xm1 = (x.real - samp[j-1].real) + (x.imag - samp[j-1].imag)*1j
                xn = (xm1.real*xn.real - xm1.imag*xn.imag) + (xm1.real*xn.imag + xm1.imag*xn.real)*1j
                ctmp = (an[j].real*xn.real - an[j].imag*xn.imag) + (an[j].real*xn.imag + an[j].imag*xn.real)*1j

                f += ctmp
            #print(ctmp)
            xunsc = ((x.real + 2.0) * par['dE'] / 4.0 + par['Vmin'])
            func_file.write(f"{xunsc:.6g} {f.real:.6g} {filt_func(xunsc - El[0], par['dt'])}\n")
            x_vals.append(xunsc); y_vals.append(f.real)
            x += dx
        func_file.close()
        x_vals = np.array(x_vals); y_vals = np.array(y_vals)
        ax.plot(x_vals+(el-El[0]), y_vals)
    ax.set_xlim(2*par['VBmin']-par['CBmax'],2*par['CBmax']-par['VBmin'])
    ax.set_xlabel("Energy (a.u.)", fontsize=14)
    ax.set_ylabel("Filter Function", fontsize=14)
    plt.savefig("func.png", dpi=120)
    


