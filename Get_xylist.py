import sys
import numpy as np
def main(get_x=False, get_y=False,real_fibpos=False,  cell_size=1., xymin= -5., xymax=5.,z=0.04, dither_set=-1 ):
    get_x=bool(get_x)
    get_y=bool(get_y)
    real_fibpos=bool(real_fibpos)
    cell_size = float(cell_size)
    xymin = float(xymin)
    xymax = float(xymax)
    z = float(z)
    dither_set = int(dither_set)
    # print cell_size
    #  print xymin, xymax, cell_size, Nbin
    if real_fibpos: # use the real position of the fibers
        filename = 'some_dependencies/fiber_location_kpc_at_z%0.3f.txt' %z
        all_pos = np.loadtxt(filename)
        if get_x :
            col_idx = 0
        if get_y :
            col_idx = 1
        col_idx = col_idx + 2 * dither_set
        xylist= all_pos[:, col_idx]

    else:
        Nbin = int((xymax - xymin)/cell_size)+1
        xylist= np.linspace(xymin, xymax, Nbin)

    # print the results
    for i in xylist:
        print i,

if __name__ == "__main__":
    main(**dict(arg.split('=') for arg in sys.argv[1:])) # kwargs
