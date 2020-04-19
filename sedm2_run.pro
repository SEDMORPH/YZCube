;+
; NAME:
;	SEDM2_RUN
;
; PURPOSE:
;	Builds images and spectra, and compile output files to run
;	through hyperion.
;
;
; CALLING SEQUENCE:
;
;	SEDM_RUN, [fileseq], snapID=snapID
;
;
; OPTIONAL INPUT:
;
;       FILESEQ: default '2xSab_00_vw'
;
; KEYWORD PARAMETERS:
;
;       SNAPID: default all
;	      REDSHIFT: default 0.04
;       TAUV: default 1.0
;       MU_D: default 0.3
;
;       SDSSIMAGE: =1 to create a mock SDSS image. Note that z!=0.
;       IMAGESIZE: For SDSSIMAGE: default 1; in arcmin
;       cell_spectra: create spectra for the cell descibed below:
;                    cell_x_offset: the offset on x-axis from the center of the first galaxy
;                    cell_y_offset: the offset on y-axis from the center of the first galaxy
;                    cell_size: the size of the cell, we now using square cells
;                    cir_fib: use a circular fiber when creating cell_spectra
;                    fib_radius   : the radius of the circular fiber, if unset, use (cell_size/2.0)
;                    arcsec: if this keyword is set, the parameters above are in arcsec.
;                            Otherwise they are in kpc. Use kpc by default
;                    with_PSF: if set, include the PSF effect when creating cell spectra
;       one_comp_dust: use 1-component dust instead of the default 2-component mag_AB_dust
;                      adopt tau_old for all stars, i.e. tau_young = tau_old
;
; EXAMPLE:
;       SEDM2_RUN, snapID=75,/skipgas
;
;
; DEPENDENCIES:
;        hyp_readsnap
;        sedm_readised.pro
;        sedm_buildsed.pro
;        sedm_mocksdssimage.pro
;        sedm_mocksdssimage_movie.pro
;        hyp_gassfh.pro
;	 sedm_spec.pro
;
; NOTES:
;       sedm2_codeunits.inc  contains code units, cosmology and conversions
;       sedm2_directories.inc contains directories
;
; MODIFICATION HISTORY:
; 	Written by:	Vivienne Wild July 2015, based on an older
; 	pre-hyperion version.
;
;-


PRO SEDM2_RUN, fileseq, snapID = snapID, redshift=redshift, tauv=tauv, mu_d = mu_d, $
    sdssimage=sdssimage,imagesize=imagesize, faceon=faceon, sfrmovie=sfrmovie, sdssmovie=sdssmovie, $
    hyperion=hyperion, gassfh = gassfh,  movieorientation=movieorientation,spectra=spectra,pca=pca, $
    centerslist=centerslist, spec_star_age=spec_star_age, with_metal=with_metal, $
    cell_spectra=cell_spectra, cell_x_offset=cell_x_offset, cell_y_offset=cell_y_offset, $
    cell_size=cell_size, cir_fib=cir_fib, fib_radius=fib_radius, arcsec=arcsec, $
    spec_style=spec_style, rtfaceon=rtfaceon, one_comp_dust=one_comp_dust, with_PSF=with_PSF, $
    plot_cell_spec=plot_cell_spec, sedm_mocksdssimage_movie=sedm_mocksdssimage_movie

  if n_elements(fileseq) eq 0 then begin
     print, 'please provide input file sequence'
     return
  endif

;fileseq = '2xSab_00_vw'

;;------------------------------------------------------------------
;;-- Fixed parameters and variables
;;------------------------------------------------------------------


  rotate_seed = 130950          ;this ensures we can reproduce the randomly selected orientation

  model_str = 'bc03'

  if n_elements(redshift) eq 0 then redshift = 0.04               ;redshift to put galaxies at
  if n_elements(tauv) eq 0 then tauv = 1.0                    ;effective optical depth in the V band
  if n_elements(mu_d) eq 0 then mu_d = 0.3                    ;fraction of optical depth in the ISM

  if n_elements(cell_x_offset) eq 0 then cell_x_offset = 0.0  ; if not specified, use the center of the first galaxy
  if n_elements(cell_y_offset) eq 0 then cell_y_offset = 0.0  ; if not specified, use the center of the first galaxy
  ; if n_elements(cell_size) eq 0 then cell_size = 1.0          ; default cell size, 1 kpc * 1kpc
  if n_elements(spec_style) eq 0 then spec_style = ''         ; SEDMoprh style by default


;;------------------------------------------------------------------
;; Directories and basic numbers
;;------------------------------------------------------------------
  @sedm2_directories.inc

  if file_test(dir_OUT) eq 0 then spawn, 'mkdir '+ dir_OUT

  nsnap_total = n_elements(file_search(dir_in+'snap_*.hdf5'))


;;------------------------------------------------------------------
;; Calculate centers of the two halo for all snapshots
;;------------------------------------------------------------------
  if KEYWORD_SET(centerslist) then SEDM2_centerslist, fileseq, indir=dir_in, outdir=dir_in

;;------------------------------------------------------------------
;; Create total optical spectra
;;------------------------------------------------------------------

  if keyword_Set(spec_star_age) then SEDM2_spec_star_age, dir_in, dir_out,tauv,mu_d, $
                                         snap = snapID, model_str=model_str, models_dir=dir_models, with_metal=with_metal

;;------------------------------------------------------------------
;; Create  optical spectra for a cell
;;------------------------------------------------------------------

 if keyword_Set(cell_spectra) then SEDM2_CELL_SPEC_STAR_AGE, dir_in, dir_out,tauv,mu_d,redshift,  cell_x_offset, cell_y_offset,$
    					                            cell_size=cell_size, cir_fib=cir_fib, fib_radius=fib_radius, arcsec=arcsec, $
                                          snap = snapID, style = spec_style, rtfaceon=rtfaceon, model_str=model_str, models_dir=dir_models,$
                                          one_comp_dust=one_comp_dust, with_metal=with_metal, with_PSF=with_PSF, dir_PSF_weight=dir_PSF_weight,$
                                          plot_cell_spec=plot_cell_spec


;;------------------------------------------------------------------
;; Create spectral indices
;;------------------------------------------------------------------

  if keyword_Set(pca) then SEDM2_pca, dir_in, dir_out, tauv, mu_d, dir_pca_data, one_comp_dust=one_comp_dust, style=spec_style, with_metal=with_metal;# snap = snapID





END
