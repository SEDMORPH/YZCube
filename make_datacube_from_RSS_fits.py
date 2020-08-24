
# coding: utf-8

# In[1]:


import os
import sys
import time
import numpy as np
from astropy.io import fits
from some_dependencies.read_cell_spec import read_cell_spec


# In[2]:


def ml_gridspec(flux, fib_x, fib_y, rlim, sigma, dim_out=[100,100],
                scale=None, ivar=None, maskvec=None,
               xomin=-15, xomax=15, yomin=-15, plot_final_image=True):
    """
    An upgrade version of ml_griddata that allows us to grid RSS input to a datacube if
    no 
    check ml_griddata for more information
    Oringal code by David Law:
    https://trac.sdss.org/browser/repo/manga/mangadrp/trunk/pro/math/ml_griddata.pro
    """


    # Dimensions
    # check whether flux is a set of a spectra
    nx_out = dim_out[1]
    ny_out = dim_out[0]
    f_shape = flux.shape
    if len(f_shape) == 1: # the input is a slice of RSS
        n_fib = len(flux)
        cube  = np.zeros([ny_out, nx_out]) 
    elif len(f_shape) == 2: # the input is the RSS
        n_fib = f_shape[0]
        n_wave = f_shape[1]
        cube = np.zeros([ny_out, nx_out, n_wave])
        
    all_dim = [ny_out, nx_out, n_fib]# Output X size, Y size, total samples
    

    # Default scale multiplier is 1.0
    if scale is None:
        scale = 1.0

    # Default inverse variance vector has everything
    # set to a constant of 1
    ivar = None
    if ivar is None:
        ivar = np.ones(n_fib)
#     print scale

    # Default mask vector has everything equal to 0
    maskvec = None
    if maskvec is None:
        maskvec = np.zeros(n_fib)


    # omit the mask part
    # maskimg=lonarr(dim_out)
    # if (keyword_set(dispvec)) then dispimg=fltarr(dim_out)
    # if (keyword_set(predispvec)) then predispimg=fltarr(dim_out)

    # Safety check that x,y,f have the same length
    nx=len(fib_x)
    ny=len(fib_y)
    nf=len(flux)
    ni=len(ivar)
    nm=len(maskvec)
    if ((nx != ny) or (nx != nf) or (nx != ni) or (nx != nm)) :
        print 'WARNING!  x,y,f,ivar,maskvec do not have the same length.'


    # X and Y output pixel coordinate arrays
    arr_xcoord = np.arange(nx_out) 
    arr_ycoord = np.arange(ny_out)
    

    # Calculate a 3d array of weights for all locations in the image
    # for all input elements.  Weight array combines a term describing
    # distance from the input element and a term describing the inverse
    # variance of the input element.
    arr_weights = np.zeros(all_dim)


    for i in range(n_fib):
        # Array defining radii in output grid away from the fiber location
        # Initialize to rlim+1 (i.e. something greater than rlim
        # radius criterion within which we actually care about the radii)
        arr_radius = np.zeros(dim_out)+ (rlim+1.0)
        # Figure out the region of influence in which we actually need to calculate radii
        xmin = int(np.max( [np.floor(fib_x[i] - rlim), 0]))
        xmax = int(np.min( [np.ceil( fib_x[i] + rlim), (nx_out-1) ]))
        xmaxp1= xmax+1
        ymin = int(np.max( [np.floor(fib_y[i] - rlim), 0]))
        ymax = int(np.min( [np.ceil( fib_y[i] + rlim), (ny_out-1) ]))
        ymaxp1= ymax+1

        # Calculate actual radii in this region of influence
        temp_nx = len(arr_xcoord[xmin:xmaxp1])
        temp_ny = len(arr_ycoord[ymin:ymaxp1])

        # python is row-first
        temp_x_2d = np.repeat(arr_xcoord[xmin:xmaxp1], temp_ny).reshape([temp_nx, temp_ny]).transpose()
        temp_y_2d = np.repeat(arr_ycoord[ymin:ymaxp1], temp_nx).reshape([temp_ny, temp_nx])
        # python is row-first
        arr_radius[ymin:ymaxp1,xmin:xmaxp1] = np.sqrt(                                (temp_x_2d-fib_x[i])**2 + (temp_y_2d-fib_y[i])**2 )
        tocalc = np.bool_(arr_radius <= rlim)
        # Weights are the exponential falloff of influence with
        # increasing distance from the input fiber location.  Things with ivar=0
        # will be given zero weight later on, keep them non-zero here
        # so that we can just track which input elements affect which
        # output pixels
        # arr_weights[:,:,i][tocalc] = np.exp(-0.5/sigma**2*arr_radius[tocalc]**2)
        arr_weights[:,:,i][tocalc] = np.exp(-0.5/sigma**2*arr_radius[tocalc]**2)



    # Figure out the normalization matrix- sum of arr_weights
    # Safety case for where there is only 1 exposure, so no 3rd dimension
    if n_fib == 1:
        matr_norm=arr_weights 
    # Sum over the 3rd dimension of arr_weights.  First make sure that
    # any input element that has zero inverse variance contributions nothing
    # to the normalization sum.  Do this by taking a logical AND between ivar
    # and 1 (which will give 0 where ivar=0, and 1 elsewhere) and
    # recasting this as a 3d array of the correct dimensions so that it
    # can simply be multiplied by the arr_weights.
    else:
        matr_norm = np.sum(arr_weights, axis=2)

    # Flag where the normalization matrix is zero; there is no good data here
    nodata=np.bool_(matr_norm == 0)
    # We don't want to divide by zero where there is no data; set the normalization
    # matrix to 1 in these cases
    matr_norm[nodata]=1.

        
    if len(f_shape) == 1: # the input is a just a slice of the RSS
        for i in range(n_fib):
            alpha= arr_weights[:,:,i] * (ivar[i] and 1) / matr_norm
            cube = cube+flux[i]*alpha
            
    elif len(f_shape) == 2: # the input is the RSS
        for i in range(n_fib):
            print "Processing fiber: %d/381...\r" %(i+1),
            alpha= arr_weights[:,:,i] * (ivar[i] and 1) / matr_norm
            # apply alpha to all wavelength. Need to make alpha to 3d at first
            d3_alpha = np.repeat(alpha[:, :, np.newaxis], n_wave, axis=2)
            cube = cube+flux[i]*d3_alpha

    cube = cube*scale
    print "Finished!"
    
    return cube


# In[3]:


def datacube_dir(fileseq, snap, cell_size=1.,
        faceon=False, style='', tauv=1.0, mu=0.3, cir_fib=False,
        with_metal=False, one_comp_dust=False, with_PSF=False):
    """
    find the dir of the datacube. The dir is where RSS are stored
    cir_fib: use the circular fiber
    """
    #where we store the output of the SEDmorph code
    str_snap = str(snap).zfill(3)
    dir_out = '/share/data/yz69/SEDMORPH/SimnsGadget3_output/'+fileseq+'/'

    fits_outstr = "_tauv%0.1f_mu%0.1f" %(tauv, mu)
    if faceon:
        fits_outstr=fits_outstr+'_fo'     
    if one_comp_dust:
        fits_outstr=fits_outstr+'_one_comp_dust'     
    if style.lower() == 'sedmorph':
        style=''
    elif style=="":
        pass
    else:
        style='_'+style
    if cir_fib:
        data_cube_dir = "DataCube"+fits_outstr+'_'+str_snap+style+'_cir_radius_%0.2f' %(cell_size)
    else:
        data_cube_dir = "DataCube"+fits_outstr+'_'+str_snap+style+'_size_%0.2f' %(cell_size)
    if with_metal:
        data_cube_dir=data_cube_dir+'_with_metal'
    if with_PSF:
        data_cube_dir=data_cube_dir+'_with_PSF'
    data_cube_dir = dir_out+data_cube_dir + '/'
    
    return data_cube_dir



def oldfits_removing_confirmation(fitsname, wait_time=5):
    if os.path.isfile(fitsname): 
        print "Fits File: %r already exits!!!" % fitsname
        print "Do you really want to remove this old one?"

        for i in range(wait_time, 0, -1):
            warning_str = "Press Ctrl+C in %ds to cancel overwriting and stop the script." %i
            sys.stdout.write(warning_str+"\r")
            time.sleep(1)
            sys.stdout.flush()
        sys.stdout.write("\r\n")
        sys.stdout.flush()
        print "Confimred, try to remove old fits file..."
        try:
            os.remove(fitsname)
            print "Old fits file removed, please write the new one."
        except OSError:
            print "Unable to remove fits file."


# In[6]:


def main(fileseq, snap, z =0.04, style='star_age',tauv=1.0, mu = 0.3,
        pixel_size_arcsec= 0.5,rlim_arcsec= 1.6, sigma_arcsec= 0.7, #check D. Law+16 
        faceon = False, one_comp_dust = False, with_metal=True, with_PSF=True,):

    snap=int(snap)
    z = float(z)
    tauv=float(tauv)
    mu = float(mu)
    faceon = bool(int(faceon))
    one_comp_dust = bool(int(one_comp_dust))
    with_metal = bool(int(with_metal))
    with_PSF = bool(int(with_PSF))



    #  pixel_size_arcsec= 0.5,rlim_arcsec= 1.6,sigma_arcsec= 0.7

    dir_in = '/share/data/common/SEDM2/major/'+fileseq+'/'

# In[7]:


# paramters
    kpc_in_arcsec = 0.79134391 # 1 arcsec is corresponding to 0.8 kpc at z=0.04
    #  arcsec_in_kpc = 1.0/kpc_in_arcsec
    fib_size = 0.79134391 # as its radius is 1 arcsec
    pixel_size = pixel_size_arcsec * kpc_in_arcsec # in side length in kpc
    rlim = rlim_arcsec * kpc_in_arcsec 
    sigma = sigma_arcsec * kpc_in_arcsec


# In[8]:


# output setting
    xomin = -5.0# kpc, the minimum of the x-axis value of the output
    xomax =  5.0
    yomin = -5.0
    yomax =  5.0


# In[9]:


    hdr = fits.Header()
    hdr['fileseq'] = (fileseq, "breif code for the simulation")
    hdr['snapshot'] = snap
    hdr['spectype'] = ( style, "type of the spectra")
    hdr['faceon'] = (faceon, "Whether rotate the galaxy to be faceon")
    hdr['Author'] ='Yirui Zheng (yz69@st-andrews.ac.uk)'
# hdr['code_ver'] = ("master-5a2dc4f", "The version of the SEDmorph code, branch-commit" ) 
    hdr['redshift'] = z
    hdr['pixel_size_arcsec'] = ( pixel_size_arcsec, "in arcsec")
    hdr['rlim_arcsec'] = (rlim_arcsec, "in arcsec")
    hdr['sigma_arcsec'] = (sigma_arcsec, "in arcsec")
    hdr['pixel_size'] = (pixel_size, "in kpc")
    hdr['rlim'] = (rlim, "in kpc")
    hdr['sigma'] = (sigma, "in kpc")
    hdr['xomin'] = (xomin, "in kpc")
    hdr['xomax'] = (xomax, "in kpc")
    hdr['yomin'] = (yomin, "in kpc")
    hdr['yomax'] = (yomax, "in kpc")
# hdr['gal_cen'] = (str(center), "xyz1,xyz2")
    hdr['spsmodel'] = ("BC03_xmiles", "Stellar Population Synthesis model")
    hdr['with_metal'] = (with_metal, "metallicity model")
    hdr['note'] = "the cell is marked by the offset on x and y aixes to the BH of the first galaxy"
    primary_hdu = fits.PrimaryHDU(header=hdr)


# In[10]:
# In[11]:


#write cell position info into a table
# out put dimensions
    nx_out = int(np.ceil( (xomax - xomin)/ pixel_size ))+1
    ny_out = int(np.ceil( (yomax - yomin)/ pixel_size ))+1
    dim_out = [ny_out, nx_out]

    to_x = np.repeat( np.arange(nx_out), ny_out).reshape([nx_out, ny_out]).transpose()
    to_y = np.repeat( np.arange(ny_out), nx_out).reshape([ny_out, nx_out])
    pixel_pos_x = to_x*pixel_size+xomin
    pixel_pos_y = to_y*pixel_size+yomin
# pos_hdu.name='cell_pos'
    x_col = fits.Column(name='x', format=str(nx_out)+'E', unit='kpc', array=pixel_pos_x)
    y_col = fits.Column(name='y', format=str(ny_out)+'E', unit='kpc', array=pixel_pos_y)
    pos_hdu = fits.BinTableHDU.from_columns(fits.ColDefs([x_col, y_col]))
    pos_hdu.name='pixel_pos'


# In[12]:


# read the position of all fibers (3 dithering sets)
    fib_pos_file = 'fiber_location_kpc_at_z%0.3f.txt' %z
    fib_pos_file = 'some_dependencies/fiber_location_kpc_at_z%0.3f.txt' %z
    #  fib_pos_file = 'fiber_location_kpc.txt'
    fib_pos_all = np.loadtxt(fib_pos_file)
    fib_x = np.concatenate((fib_pos_all[:,0], fib_pos_all[:,2], fib_pos_all[:,4]))
    fib_y = np.concatenate((fib_pos_all[:,1], fib_pos_all[:,3], fib_pos_all[:,5]))
    n_fib = len(fib_x)


# In[13]:


#write wavelength info into a table
    spec=read_cell_spec(fileseq, snap, cell_size=fib_size, 
                                cell_x_offset=fib_x[0], cell_y_offset=fib_y[0],
                                style=style,  with_metal=with_metal, cir_fib=True,with_PSF=with_PSF )
    wave = spec['WAVE'][0]
    wave_col = fits.Column(name='wave', format='E', unit='Angstorm', array=wave)
    wave_hdu = fits.BinTableHDU.from_columns(fits.ColDefs([wave_col]))
    wave_hdu.name='wavelength'


# In[14]:


# Scaling factor to be applied to the output.  Required
# as input flux units tend to be per fiber area, while
# output flux units should be per spaxel.
    scale = pixel_size**2 / (np.pi*fib_size**2)

    if (len(fib_x) != n_fib) or (len(fib_y) != n_fib):
        print "length of the fiber postion array does no match n_fib. Exit..."
        sys.exit()

# read the first fiber and initialise RSS 2-d array
    temp = read_cell_spec(fileseq, snap, cell_size=fib_size, 
                            cell_x_offset=fib_x[0], cell_y_offset=fib_y[0],
                            style=style,  with_metal=with_metal, cir_fib=True,with_PSF=with_PSF )
    wave = temp['WAVE'][0]

    RSS = np.zeros([n_fib, len(wave)])
    RSS_notau = np.zeros([n_fib, len(wave)])
    RSS[0]=temp['SPEC_TAU'][0]
    RSS_notau[0]=temp['SPEC_NOTAU'][0]
# read the other fibers
    for i in range(1, n_fib):
        temp = read_cell_spec(fileseq, snap, cell_size=fib_size, 
                                cell_x_offset=fib_x[i], cell_y_offset=fib_y[i],
                                style=style,  with_metal=with_metal, cir_fib=True,with_PSF=with_PSF )
        RSS[i]=temp['SPEC_TAU'][0]
        RSS_notau[i]=temp['SPEC_NOTAU'][0]

# shift the data to make all fiber position in the unit of pixels
    fib_x_pix = (fib_x - xomin)/ pixel_size
    fib_y_pix = (fib_y - yomin)/ pixel_size
    rlim_pix = rlim/pixel_size # convert to pixel unit
    sigma_pix = sigma/pixel_size # convert to pixel unit

# out put dimensions
    nx_out = int(np.ceil( (xomax - xomin)/ pixel_size ))+1
    ny_out = int(np.ceil( (yomax - yomin)/ pixel_size ))+1
    dim_out = [ny_out, nx_out]

    to_x = np.repeat( np.arange(nx_out), ny_out).reshape([nx_out, ny_out]).transpose()
    to_y = np.repeat( np.arange(ny_out), nx_out).reshape([ny_out, nx_out])
    pixel_pos_x = to_x*pixel_size+xomin
    pixel_pos_y = to_y*pixel_size+yomin
#     temp_idx = np.bool_(pixel_pos_x>=-5)*np.bool_(pixel_pos_x<=5)*np.bool_(pixel_pos_y>=-5)*np.bool_(pixel_pos_y<=5)

    cube = ml_gridspec(flux= RSS, fib_x=fib_x_pix, fib_y=fib_y_pix, rlim=rlim_pix,
                       scale=scale, sigma=sigma_pix, dim_out=dim_out, plot_final_image =False )
    cube = np.transpose(cube, (2,0,1))

    cube_notau = ml_gridspec(flux= RSS_notau, fib_x=fib_x_pix, fib_y=fib_y_pix, rlim=rlim_pix,
                       scale=scale, sigma=sigma_pix, dim_out=dim_out, plot_final_image =False )
    cube_notau = np.transpose(cube_notau, (2,0,1))


# In[15]:


# write to hdu
    hdu_list = [primary_hdu, wave_hdu, pos_hdu]
    cube_hdu = fits.ImageHDU( cube  ) 
    cube_hdu.name = 'Datacube_dustted'
    cube_hdu.header['unit'] = "Lsun/Angstrom"
    hdu_list.append(cube_hdu)
    cube_notau_hdu = fits.ImageHDU( cube_notau  ) 
    cube_notau_hdu.name = 'Datacube_nodust'
    cube_notau_hdu.header['unit'] = "Lsun/Angstrom"
    hdu_list.append(cube_notau_hdu)


# In[16]:


#Write all tables created above into same fits file
    datacube_hdul = fits.HDUList(hdu_list)
    datacube_hdul.info()


# In[17]:


    cube_dir=datacube_dir(fileseq, snap, cell_size=fib_size,
                         style=style,  with_metal=with_metal,
                         cir_fib=True,with_PSF=with_PSF )
    print cube_dir


# In[21]:


    fits_outstr = "_tauv%0.1f_mu%0.1f" %(tauv, mu)
    if faceon:
        fits_outstr=fits_outstr+'_fo'
    if one_comp_dust:
        fits_outstr=fits_outstr+'_one_comp_dust'
    datacube_fits = cube_dir+"DataCube_"+fileseq+fits_outstr+'_'+str(snap).zfill(3)+style+'.fits'
    print "save the datacube to file: %s" %datacube_fits
# datacube_fits = 'test.fits'
    oldfits_removing_confirmation(datacube_fits, wait_time=5)
    datacube_hdul.writeto(datacube_fits)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
            **dict(arg.split('=') for arg in sys.argv[3:])) # kwargs
