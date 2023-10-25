
def write_cube(rho, filename, atm_ids, coords, atm, par, ist):
    
    natoms = len(atm_ids)
    
    cube_file = open(f'{filename}.cube', "w")
    cube_file.write("CUBE FILE\n")
    cube_file.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
    
    cube_file.write(f"{natoms} {par['xmin']:12.6f}{par['ymin']:12.6f}{par['zmin']:12.6f}\n")
    cube_file.write(f"{ist['nx']} {par['dx']} {0.0:12.6f} {0.0:12.6f}\n");
    cube_file.write(f"{ist['ny']} {0.0:12.6f} {par['dy']} {0.0:12.6f}\n");
    cube_file.write(f"{ist['nz']} {0.0:12.6f} {0.0:12.6f} {par['dz']}\n");
    
    for i in range(natoms):
        cube_file.write(f"{atm_ids[i]} {0.0:12.6f} {coords[i,0]:12.6f} {coords[i,1]:12.6f} {coords[i,2]:12.6f}\n")
    
    for jx in range(ist['nx']):
        for jy in range(ist['ny']):
            for jz in range(ist['nz']):
                cube_file.write(f"{rho[jx,jy,jz]:12.5f} ")
                if jz % 6 == 5:
                    cube_file.write('\n')
            cube_file.write('\n')
    cube_file.close()
    
