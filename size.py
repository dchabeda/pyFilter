import numpy as np

# Misc functions
def get_dot_size(coords):
    disx = 0.0; disy = 0.0; disz = 0.0
    for i in range(len(coords)-1):
        for j in range(i+1,len(coords)):
            dx = np.sqrt((coords[i,0] - coords[j,0])**2)
            dy = np.sqrt((coords[i,1] - coords[j,1])**2)
            dz = np.sqrt((coords[i,2] - coords[j,2])**2)
            if dx > disx:
                disx = dx
            if dy > disy:
                disy = dy
            if dz > disz:
                disz = dz
                
    return disx, disy, disz
 