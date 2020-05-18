
# coding: utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM



def angdiam_fib(z, angsize=1.0, Om0=0.3, H0=70):
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    lumdist = cosmo.luminosity_distance(z).value # in Mpc
    diam = angsize/(60.*60.*360/(2*np.pi)) #radians
    return diam*lumdist*1000 /(1+z)**2 # in kpc




def main(z):
    z=float(z)
    arcsec_to_kpc = angdiam_fib(z,angsize=1.0)
# D. Law+2015 figure 6 
# https://iopscience.iop.org/article/10.1088/0004-6256/150/1/19 
    mm_to_arcsec = 16.6   # 1 mm equals to 16.6 arcsec
    # MaNGA 
    raw_file= 'ma134-56995-1.par'
    

    try:
        data=np.genfromtxt(raw_file, dtype=None, skip_header=69)
    except:
        print "Please download a MaNGA harness metrology file and put it under the same folder"
        print "The one we used is ma134-56995-1.par, you could find it here:"
        print "https://svn.sdss.org/public/repo/manga/mangacore/tags/v1_6_2/metrology/ma134/ma134-56995-1.par"


    n_fib = len(data)
    fib_x = np.zeros(n_fib)
    fib_y = np.zeros(n_fib)
    fib_type = np.empty(n_fib, dtype='S5')


    for i in range(n_fib):
        fib_x[i] = data[i][9]
        fib_y[i] = data[i][10]
        fib_type[i] = data[i][5]




# add dithering
    fib_x1 = fib_x - 0.075 # dithering 1
    fib_y1 = fib_y + 0.075/np.sqrt(3) # dithering 1

    fib_x2 = fib_x - 0.075 # dithering 1
    fib_y2 = fib_y - 0.075/np.sqrt(3) # dithering 1

# select the IFU fibers and abandon sky fibers
    IFU_fib_idx = np.bool_(fib_type!='SKY')
# original position
    IFU_fib_x = fib_x[IFU_fib_idx]
    IFU_fib_y = fib_y[IFU_fib_idx]
# dithering 1
    IFU_fib_x1 = fib_x1[IFU_fib_idx]
    IFU_fib_y1 = fib_y1[IFU_fib_idx]
# dithering 2
    IFU_fib_x2 = fib_x2[IFU_fib_idx]
    IFU_fib_y2 = fib_y2[IFU_fib_idx]


# In[21]:


    IFU_fib_x_kpc  = IFU_fib_x  * mm_to_arcsec * arcsec_to_kpc
    IFU_fib_y_kpc  = IFU_fib_y  * mm_to_arcsec * arcsec_to_kpc
    IFU_fib_x1_kpc = IFU_fib_x1 * mm_to_arcsec * arcsec_to_kpc
    IFU_fib_y1_kpc = IFU_fib_y1 * mm_to_arcsec * arcsec_to_kpc
    IFU_fib_x2_kpc = IFU_fib_x2 * mm_to_arcsec * arcsec_to_kpc
    IFU_fib_y2_kpc = IFU_fib_y2 * mm_to_arcsec * arcsec_to_kpc


# In[24]:


    header="""
    To know where to place the circular apertures on the simulated galaxy,
    we take a representative RSS file for the largest MaNGA IFU bundle of 127 fibres.
    All in unit of kpc.
    Taking 1 mm equals to 16.6 arcsec and 1 arcsec = 0.7913 kpc at redshift z = 0.04. 
    The first and second columns are the original x and y position of the fibers(0).
    The third and forth columns are the dithered x and y position of the fibers(dithered 1).
    The fifth and sixth columns are the dithered x and y position of the fibers(dithered 2).
    Check Law+2015 figure 6 & 9.
    1
    | \\
    |   0
    | /
    2
    """
    save_data = np.c_[IFU_fib_x_kpc, IFU_fib_y_kpc, IFU_fib_x1_kpc, IFU_fib_y1_kpc, IFU_fib_x2_kpc, IFU_fib_y2_kpc] 
    filename = 'fiber_location_kpc_at_z%0.3f.txt' %z
    np.savetxt(filename, save_data, fmt='%9.2f', header=header)

if __name__ == "__main__":
    main(sys.argv[1])
