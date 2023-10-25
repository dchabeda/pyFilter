import numpy as np

# Read functions
def read_conf(atm):
    '''Read the atom configuration and return the centered position of all atoms'''
    
    conf_file = open("conf.par", "r")
    coords = []
    atm_ids = []
    for i,line in enumerate(conf_file):
        if i == 0:
            natom = int(line.split()[0])
        else:
            if line.split() != []:
                atm_id = atm['id'][line.split()[0]] # get atom ids from the atm dict
                atm_ids.append(atm_id) 
                x = float(line.split()[1]); y = float(line.split()[2]); z = float(line.split()[3])
                coords.append([x, y, z])
    
    if int(len(coords)) != int(natom):
        print('length of conf does not match natom!')
        print('Exiting...\n'); exit()
    
    coords = np.array(coords, dtype=np.float64)
    
    com = np.mean(coords, axis=0)
    
    coords[:,0] -= com[0]
    coords[:,1] -= com[1]
    coords[:,2] -= com[2]
    
    return natom, atm_ids, coords
    
def read_pot(atm_names):
    pots = {}
    for atm_name in atm_names:
        pot = np.loadtxt(f"pot{atm_name}.par", dtype=np.float64) # load in the potential
        pot[:,1] -= pot[-1,1] # shift the potential so it decays to 0
        x_diff = pot[1,0] - pot[0,0]
        pots[atm_name] = (x_diff, pot)
    
    return pots