
def interpolate(dr_grid, potatom):
    n = len(potatom[1])
    dr_pot = potatom[0]
    
    i = int(dr_grid/dr_pot)
    if (i > (n - 2)): return 0.0 #print('hit!');
    
    a = (potatom[1][i+1,1] - potatom[1][i,1])/ (potatom[1][i+1,0] - potatom[1][i,0])
    b = potatom[1][i,1] - potatom[1][i,0] * a
    
    return a * dr_grid + b
