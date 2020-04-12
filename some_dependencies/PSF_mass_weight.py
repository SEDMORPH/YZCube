
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from SDSS_PSF import SDSS_PSF_kernel
from astropy.convolution import convolve
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

def angdiam_fib(z, angsize=1.0, Om0=0.3, H0=70):
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    lumdist = cosmo.luminosity_distance(z).value # in Mpc
    diam = angsize/(60.*60.*360/(2*np.pi)) #radians
    return diam*lumdist*1000 /(1+z)**2 # in kpc


# output setting
x_min = -5.0# kpc, the minimum of the x-axis value of the output
x_max =  5.0
y_min = -5.0
y_max =  5.0





def main(z):
    z=float(z)
# paramters
    kpc_in_arcsec = angdiam_fib(z,angsize=1.0)# 1 arcsec is corresponding to 0.8 kpc at z=0.04
    fib_size = 1.0 * kpc_in_arcsec # as its radius is 1 arcsec
    img_res = 0.05 # kpc, image resolution



    for band_name in "ugriz":
# #### Caculate the PSF kernel
        comb_ker = SDSS_PSF_kernel(bandname=band_name, kpc_in_arcsec=kpc_in_arcsec, cell_size=img_res) 



# create the real image
        real_img_nx = int(np.ceil( (x_max - x_min)/ img_res ))+1
        real_img_ny = int(np.ceil( (y_max - y_min)/ img_res ))+1
        real_img = np.zeros([real_img_ny, real_img_nx])
        real_img[(real_img_ny-1)/2, (real_img_nx-1)/2] = 1.0
# convolve with PSF
        real_img_PSF = convolve(real_img, comb_ker)


# put the fiber on the positive part of the x-axis
        fib_x = np.arange(real_img_nx)* img_res
        fib_y = np.zeros_like(fib_x)
        rimg_x = np.repeat(np.arange(real_img_nx), real_img_ny).reshape([real_img_nx, real_img_ny]).transpose()*img_res+x_min
        rimg_y = np.repeat(np.arange(real_img_ny), real_img_nx).reshape([real_img_ny, real_img_nx])*img_res+y_min


# #### Now, do the observation. We need convert to all the stuff to the unit of output pixel size.

# In[13]:


        mass_weight = np.zeros_like(fib_x)
        for i in range(len(fib_x)):
            rlist = np.sqrt( (rimg_x - fib_x[i])**2 + (rimg_y - fib_y[i])**2   )
            temp_idx = np.bool_(rlist <= fib_size)
            mass_weight[i] = np.sum(real_img_PSF[temp_idx])
            if mass_weight[i] < 1e-6:
                break




        weight_col = fits.Column(name='PSF_weight', format='E', array=mass_weight )
        hdu = fits.BinTableHDU.from_columns(fits.ColDefs([weight_col]))
        hdu.name='weight'
        if not os.path.isdir("./PSF_mass_weight/"):
            os.mkdir("./PSF_mass_weight/")
        fits_name = "./PSF_mass_weight/"+band_name+"_band_PSF_mass_weight_res_"+str(img_res)+'kpc_at_z%.3f.fits' %z
        try:
            hdu.writeto(fits_name)
        except:
            os.remove(fits_name)
            hdu.writeto(fits_name)

if __name__ == "__main__":
    main(sys.argv[1])
