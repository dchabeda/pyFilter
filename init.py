from fd_init import *

#Init Functions
def init():
    '''Initialize  job specs, configuration, grid, and potentials'''
    par = {}; ist = {}; atm = {}
    print("Reading input and initializing job specifications...\n")

    with open("input.par", "r") as f:
        init_file = f.readlines()
    
    atm['id'] = {'Cd':48, 'Se':34, 'Cs':55, 'Pb':82, 'I':53, 'Br':35, 'P1':1, 'P2':2, 'P3':3, 'P4':4, 'PC5':5, 'PC6':6}
    atm['name'] = {48:'Cd', 34:'Se', 55:'Cs', 82:'Pb', 53:'I', 35:'Br', 1:'P1', 2:'P2', 3:'P3', 4:'P4', 5:'PC5', 6:'PC6'}
    
    # Grab inputs from the init_file
    ist['nx'] = int(init_file[0].split('\n')[0])
    ist['ny'] = int(init_file[1].split('\n')[0])
    ist['nz'] = int(init_file[2].split('\n')[0])
    ist['ms'] = int(init_file[3].split('\n')[0])
    ist['ns'] = int(init_file[4].split('\n')[0])
    ist['nCheby'] = int(init_file[5].split('\n')[0])
    par["VBmin"] = float(init_file[6].split()[0])
    par["VBmax"] = float(init_file[6].split()[1])
    par["CBmin"] = float(init_file[7].split()[0])
    par["CBmax"] = float(init_file[7].split()[1])
    par['setdGridFlag'] = int(init_file[9].split()[0])
    if par['setdGridFlag']:
        par['dx'] = float(init_file[10].split()[0])
        par['dy'] = float(init_file[10].split()[0])
        par['dz'] = float(init_file[10].split()[0])
    # Additional inputs can be easily introduced (or ignored) here
    
    # Hardcoded parameters
    par['fermiEnergy'] = -0.180
    par['Ekinmax'] = 10.0
    ist['nx_1'] = 1/ist['nx']
    ist['ny_1'] = 1/ist['ny']
    ist['nz_1'] = 1/ist['nz']
    ist['ngrid'] = ist['nx']*ist['ny']*ist['nz']
    ist['natomtype'] = 20
    ist['mstot'] = ist['ns']*ist['ms'] # Total number of energies computed (number of filters * energy windows/filter)
    par['Vmin'] = 1.0e10;
    par['Vmax'] = -1.0e10;
    
    # Read the atom configuration
    natom, atm_ids, coords = read_conf(atm)
    atm_names = [atm['name'][i] for i in atm_ids]
    
    ist['natom'] = natom
    print(f"Number of atoms in the system, natom = {ist['natom']}")
    print(f"Number of grid points used: nx = {ist['nx']} ny = {ist['nx']} nz = {ist['nx']}")
    print(f"Number of filter cycles, ns = {ist['ns']}\nNumber of states per filter, ms = {ist['ms']}")
    print(f"Length of Newton interpolation used, ncheby = {ist['nCheby']}")
    print(f"VBmin = {par['VBmin']}\t, VBmax = {par['VBmax']}")
    print(f"CBmin = {par['CBmin']}\t, CBmax = {par['CBmax']}")
    print(f"Kinetic energy maximum = {par['Ekinmax']}")
    if par['setdGridFlag']:
          print(f"Grid spacing: dx = {par['dx']}\tdy = {par['dy']}\tdz = {par['dz']}")
    
  
    #print('The atomic configuration: ', coords)
    # Plot the atomic configuration
    
    print('\nVisualization of atomic configuration: \n')
    '''fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(projection='3d')
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(projection='3d')
    style = {48:['tab:blue',5], 34:['darkorange',8], 1:['darkorange',3], 2:['tab:blue',3]}
    
    for i in range(natom):
        ax.scatter(coords[i][0], coords[i][1], coords[i][2], color=style[atm_ids[i]][0], linewidth=style[atm_ids[i]][1])
    
    ax.set_xlabel('X (Bohr)')
    ax.set_ylabel('Y (Bohr)')
    ax.set_zlabel('Z (Bohr)')
    ax.set_zticks([])'''
    #ax.view_init(elev=3, azim=3, vertical_axis='y')
    
    
    # Get dot dimensions
    xd, yd, zd = get_dot_size(coords)
    print(f'The dot size in Bohr, xd = {xd:.2f} yd = {yd:.2f} zd = {zd:.2f}')
    xd = round(0.5 * xd + 5.0)
    yd = round(0.5 * yd + 5.0)
    zd = round(0.5 * zd + 5.0)
    print(f"Box (quadrant) dimensions: xd = {xd:.2f} yd = {yd:.2f} zd = {zd:.2f}")
    
    # Initial parameters for the potential and grid
    # x-direction
    par['xmin'] = -xd; par['xmax'] = xd
    mx = 1.0
    if par['setdGridFlag'] == 0:
        par['dx'] = (par['xmax'] - par['xmin'])/ ist['nx']
    par['dkx'] = (2*np.pi)/(ist['nx']*par['dx'])
    
    # y-direction
    par['ymin'] = -yd; par['ymax'] = yd
    my = 1.0
    if par['setdGridFlag'] == 0:
        par['dy'] = (par['ymax'] - par['ymin'])/ ist['ny']
    par['dky'] = (2*np.pi)/(ist['ny']*par['dy'])
    
    # z-direction
    par['zmin'] = -zd; par['zmax'] = zd
    mz = 1.0
    if par['setdGridFlag'] == 0:
        par['dz'] = (par['zmax'] - par['zmin'])/ ist['nz']
    par['dkz'] = (2*np.pi)/(ist['nz']*par['dz'])
    
    par['dv'] = par['dx']*par['dy']*par['dz']
    par['dr'] = np.sqrt(par['dx']**2 + par['dy']**2 + par['dz']**2)
    print(f"Grid point spacing: dx = {par['dx']:.4f} dy = {par['dy']:.4f} dz = {par['dz']:.4f} dv = {par['dv']:.4f} dr = {par['dr']:.4f}")

    #Initializing k^2
    ksqrx = np.zeros(ist['nx'])
    for j in range(1, int(ist['nx']/2) + 1):
        #print('ind1: ', j); print('ind2: ', ist['nx'] - j)
        ksqrx[j] = ksqrx[ist['nx']-j] = (0.5 * (j*par['dkx'])**2)#*(ist['nx_1']*ist['ny_1']*ist['nz_1'])/mx
        
    ksqry = np.zeros(ist['ny'])
    for j in range(1, int(ist['ny']/2) + 1):
        ksqry[j] = ksqry[ist['ny']-j] = (0.5 * (j*par['dky'])**2)#*(ist['nx_1']*ist['ny_1']*ist['nz_1'])/my

    ksqrz = np.zeros(ist['nz'])
    for j in range(1, int(ist['nz']/2) + 1):
        ksqrz[j] = ksqrz[ist['nz']-j] = (0.5 * (j*par['dkz'])**2)#*(ist['nx_1']*ist['ny_1']*ist['nz_1'])/mz
  
    par['Ekinmax'] #*= (ist['nx_1']*ist['ny_1']*ist['nz_1'])
    
    
    ksqr = np.zeros((ist['nx'],ist['ny'],ist['nz']))            
    for jx in range(ist['nx']):
        for jy in range(ist['ny']):
            for jz in range(ist['nz']):
                ksqr[jx,jy,jz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz]
                if ksqr[jx,jy,jz] > par['Ekinmax']: ksqr[jx,jy,jz] = par['Ekinmax']
    

    # Read the psuedopotentials
    atm_types = list(set(atm_ids))
    atm_name_set = [atm['name'][i] for i in atm_types]
    pots = read_pot(atm_name_set)
    
    # Initialize the grid
    vx = []; vy = []; vz = []
    if par['setdGridFlag'] == 0:
        for jx in range(ist['nx']): vx.append(par['xmin'] + jx*par['dx'])
        for jy in range(ist['ny']): vy.append(par['ymin'] + jy*par['dy'])
        for jz in range(ist['nz']): vz.append(par['zmin'] + jz*par['dz'])
    if par['setdGridFlag'] == 1:
        par['xmin'] = - ist['nx']*par['dx'] / 2
        par['ymin'] = - ist['ny']*par['dy'] / 2
        par['zmin'] = - ist['nz']*par['dz'] / 2
        for jx in range(ist['nx']): vx.append(par['xmin'] + jx*par['dx'])
        for jy in range(ist['ny']): vy.append(par['ymin'] + jy*par['dy'])
        for jz in range(ist['nz']): vz.append(par['zmin'] + jz*par['dz'])
    
        
    # Initialize the big potential at all grid points
    print("\nCalculating the local potential...")
    potl = np.zeros((ist['nx'], ist['ny'], ist['nz']))
    for jx in range(ist['nx']):
        for jy in range(ist['ny']):
            for jz in range(ist['nz']):
                sum_pot = 0.0
                for iatom in range(natom):
                    atm_name = atm['name'][atm_ids[iatom]]
                    dx = vx[jx] - coords[iatom,0]
                    dy = vy[jy] - coords[iatom,1]
                    dz = vz[jz] - coords[iatom,2]
                    dr_grid = np.sqrt(dx**2 + dy**2 + dz**2)
                    
                    sum_pot += interpolate(dr_grid, pots[atm_name])
                potl[jx,jy,jz] = sum_pot
                if par['Vmax'] < potl[jx,jy,jz]: par['Vmax'] = potl[jx,jy,jz]
                if par['Vmin'] > potl[jx,jy,jz]: par['Vmin'] = potl[jx,jy,jz]
    print("\tdone.")
    write_cube(potl, 'localPot', atm_ids, coords, atm, par, ist)
    
    par['dE'] = 0.5 * np.pi**2 / (mx*par['dx']**2) + 0.5 * np.pi**2 / (my*par['dy']**2) + 0.5 * np.pi**2 / (mz*par['dz']**2)
    
    print(f"Potential energy range: Vmin = {par['Vmin']:.5f} Vmax = {par['Vmax']:.5f} dV = {par['Vmax']-par['Vmin']:.5f}\n")
    print(f"Energy separation for kinetic energy, dT = {par['dE']:.5f}\n")
    
    
    #Setting the energy grid!
    egrid = []
    #e_range = (par['CBmax'] - par['CBmin']) + (par['VBmax'] - par['VBmin'])
    #msCB = int(ist['ms']*(par['CBmax'] - par['CBmin'])/e_range)
    msCB = int(ist['ms']/2)
    msVB = int(ist['ms'] - msCB)
    print(f"msVB = {msVB} msCB = {msCB}")
    if msCB < 1 or msVB < 1: print('Error initializing energy grid. msCB or msVB < 1!\n'); exit()
    delVB = (par['VBmax'] - par['VBmin'])/msVB
    for jx in range(msVB):
        egrid.append(round(par['VBmin'] + jx*delVB, 6))
    delCB = (par['CBmax'] - par['CBmin'])/msCB
    for jx in range(msCB):
        egrid.append(round(par['CBmin'] + jx*delCB,6))
    print(f"Spacing between states, VB = {delVB:.5f} CB = {delCB:.5f}")
    
    np.savetxt("Egrid.dat",egrid,fmt='%.6f', delimiter='\n') # save the energy grid file for later use
    
    return potl, ksqr, egrid, par, ist, atm_ids, coords, atm

def init_psi(par, ist, seed):
    '''Creates a random wavefunction on the grid. Each point is a uniform random value between -1 and 1'''
    np.random.seed(seed)
    
    psi = 1.0 - 2.0 * np.random.rand(ist['nx'],ist['ny'],ist['nz']) # Create a random wavefunc on [-1,1]
        
    norma = norm(psi, par['dv']) # Get the normalization constant for the wavefunction
    if norma == 0:
        print('Norm of 0! Exiting\n'); exit()
    
    psi *= norma
    
    return psi
  