import numpy as np
from astropy.convolution import convolve, Gaussian2DKernel
def SDSS_PSF_kernel(bandID=None, bandname =None, cell_size=0.5, kpc_in_arcsec=None):
    """
    Calculate the PSF convolution kernel using SDSS data, 
    check the coments given by Melina in sedm2_sdssimage.pro
    
    Please give either a bandID or a bandname. 
    kpc_in_arcsec, how many kpc is 1 arcsec corresponding to 
    at the required redshift. If it is not given, assume z=0.04
    """ 
    
    #-- CCD pixel size in arcsec, http://classic.sdss.org/dr3/instruments/imager/index.html
    pixel_size = 0.396
    band_list = np.array(['u','g','r','i','z'])
    if bandID is not None:
        print "Check the band you are using: %s" %band_list[bandID]
    elif bandname is not None:
        if bandname in band_list:
            print "Using %s band" %bandname
            bandID = np.where(band_list == bandname)[0][0]
#             print bandID
        else:
            print "The specified band is not availiable, please use one of the following bands"
            print band_list
            sys.exit()
    else:
        print "please give either a bandID or a bandname"
        sys.exit()
    
    if kpc_in_arcsec is None:
        kpc_in_arcsec = 0.79134391 # 1 arcsec is corresponding to 0.8 kpc at z=0.04

    #-- PSF parameters (Milena)
    #-- Values from SDSS to construct a 2-gaussian psf (median of values for 500000 SDSS fields)
    psf_s12g = np.array([1.44527,1.36313,1.21196,1.114143,1.20288]) # in pixels
    psf_s22g = np.array([3.1865,3.06772,2.82102,2.75094,3.09981]) # in pixels
    psf_ratio =np.array([0.081989,0.081099,0.06811,0.059653,0.054539])
    psf_width = np.array([1.54254,1.44980,1.31878,1.25789,1.29391])/pixel_size # originally in arcsec

    psf_s12g_cell = psf_s12g * pixel_size * kpc_in_arcsec / cell_size
    # the size of kernel, in number of pixel in side length, cover 10 kpc
    ker_size = int(np.ceil(10.0/cell_size)) 
    if (ker_size%2) !=1:
        ker_size =  ker_size+1
    kernel1 = Gaussian2DKernel(psf_s12g_cell[bandID], x_size=ker_size, y_size=ker_size)

    psf_s22g_cell = psf_s22g * pixel_size * kpc_in_arcsec / cell_size
    kernel2 = Gaussian2DKernel(psf_s22g_cell[bandID], x_size=ker_size, y_size=ker_size)

    k1_max = kernel1.array.max()
    k2_max = kernel2.array.max()

#     print k1_max, k2_max

    k2_max_s = k1_max * psf_ratio[bandID] # what the k2_max should be
#     print k2_max_s

    kernel2_new = kernel2 * (1./ k2_max *k2_max_s)
#     print kernel2_new.array.max()/kernel1.array.max() , psf_ratio[bandID] 

    comb_ker = (kernel1+kernel2_new) * (1./ (kernel1.array.sum() + kernel2_new.array.sum()) )
#     print comb_ker.array.sum()[0]
    
    return comb_ker

#  def PSF_convovle(datacube, )
