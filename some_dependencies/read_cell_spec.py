from astropy.io import fits
from decimal import Decimal, ROUND_HALF_UP

def read_cell_spec(fileseq, snap, cell_x_offset=0, cell_y_offset=0, cell_size=1.,
        faceon=False, style='', storage_style=None, tauv=1.0, mu=0.3, cir_fib=False,
        with_metal=False, one_comp_dust=False, with_PSF=False):
    """
    Read the cell spectra created by the SEDmorph code sedm2_cell_spec.pro

    storage_style:the storage arrangement of the cell_spec has been changed in Oct
                  -->Aug, the older one used in August
                  -->if none of the above is specificied , choose the up_to_date method
    cir_fib: use the circular fiber
    """
    #where we store the output of the SEDmorph code
    str_snap = str(snap).zfill(3)
    dir_out = '/share/data/yz69/SEDMORPH/SimnsGadget3_output/'+fileseq+'/'
    
    if storage_style == "Aug":
        dir_spec = dir_out
        cell_str = 'cell_'+'{0:+d}'.format(int(cell_x_offset)).zfill(3)+ \
                    '{0:+d}'.format(int(cell_y_offset)).zfill(3)
        cell_str = cell_str+'_size_'+'{0:0.1f}'.format(cell_size)
        if faceon:
            fits_file = dir_spec+cell_str+'_spec_tauv1.0_mu0.3_fo_'+str_snap+style+'.fits'
        else:
            fits_file = dir_spec+cell_str+'_spec_tauv1.0_mu0.3_'+str_snap+style+'.fits'
            
    else: # use the newest one, which should be more reasonable.
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



        # round up to deal with strange numbers like -0.985
        cell_size= float(Decimal(str(cell_size)).quantize(Decimal('1e-2'),  ROUND_HALF_UP))
        cell_x_offset= float(Decimal(str(cell_x_offset)).quantize(Decimal('1e-2'),  ROUND_HALF_UP))
        cell_y_offset= float(Decimal(str(cell_y_offset)).quantize(Decimal('1e-2'),  ROUND_HALF_UP))
        if cir_fib:
            data_cube_dir = "DataCube"+fits_outstr+'_'+str_snap+style+'_cir_radius_%0.2f' %(cell_size)
        else:
            data_cube_dir = "DataCube"+fits_outstr+'_'+str_snap+style+'_size_%0.2f' %(cell_size)
        if with_metal:
            data_cube_dir=data_cube_dir+'_with_metal'
        if with_PSF:
            data_cube_dir=data_cube_dir+'_with_PSF'
        data_cube_dir = data_cube_dir + '/'

        if cir_fib:
            cell_str = "cell_%+0.2f%+0.2f_cir_radius_%0.2f" %(cell_x_offset, cell_y_offset, cell_size)
        else:
            cell_str = "cell_%+0.2f%+0.2f_size_%0.2f" %(cell_x_offset, cell_y_offset, cell_size)
        fits_file = dir_out+ data_cube_dir+'spec_'+cell_str+'.fits'

    return fits.getdata(fits_file)
