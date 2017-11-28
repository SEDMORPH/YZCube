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
;	REDSHIFT: default 0.04
;       TAUV: default 1.0
;       MU_D: default 0.3
;
;       SDSSIMAGE: =1 to create a mock SDSS image. Note that z!=0. 
;       IMAGESIZE: For SDSSIMAGE: default 1; in arcmin
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
    centerslist=centerslist

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

;;------------------------------------------------------------------
;; Directories and basic numbers
;;------------------------------------------------------------------
  @sedm2_directories.inc

  if file_test(dir_OUT) eq 0 then spawn, 'mkdir '+ dir_OUT

  nsnap_total = n_elements(file_search(dir_in+'snap_*.hdf5'))

;;------------------------------------------------------------------
;; Create SFHs for all gas and star particles (stored, only run once)
;;------------------------------------------------------------------

;if n_elements(file_search(dir_in+'*_gassfh.sav')) lt nsnap_total

  if keyword_set(gassfh) then $ ;only compute if they are not all there
     SEDM2_GASSFH, fileseq, dir_IN, models_dir=dir_models, model_str = model_str

;;------------------------------------------------------------------
;; Create output files for Hyperion
;;------------------------------------------------------------------
  if keyword_set(hyperion) then SEDM2_HYP_OUTPUT, fileseq, dir_models, indir=dir_in, outdir=dir_out

;;------------------------------------------------------------------
;; Calculate centers of the two halo for all snapshots
;;------------------------------------------------------------------
  if KEYWORD_SET(centerslist) then SEDM2_centerslist, fileseq, indir=dir_in, outdir=dir_out


;;------------------------------------------------------------------
;; Create mock SDSS images
;;------------------------------------------------------------------
  
  if keyword_Set(sdssimage) then begin
     if redshift eq 0 then message, 'Redshift must be >0 to make an image'

     filterlist = 'sdss.lis'
     if n_elements(imagesize) eq 0 then imagesize=1. 

     SEDM2_sdssimage, dir_in, dir_out, dir_filters,dir_code,filterlist, redshift, tauv,mu_d, imagesize,$
                      snap = snapID, model_str=model_str, $
                      models_dir=dir_models,rotate_seed=rotate_seed, faceon=faceon

  endif
 

;;------------------------------------------------------------------
;; Create total optical spectra
;;------------------------------------------------------------------
  
  if keyword_Set(spectra) then SEDM2_spec, dir_in, dir_out,tauv,mu_d, $
                                           snap = snapID, model_str=model_str, models_dir=dir_models


;;------------------------------------------------------------------
;; Create spectral indices
;;------------------------------------------------------------------
  
  if keyword_Set(pca) then SEDM2_pca, dir_in, dir_out, tauv, mu_d;, snap = snapID

  
  
;;------------------------------------------------------------------
;; Create movies: this relies on image and SED files already being made
;;------------------------------------------------------------------

  ;;-- SFR movie
  if keyword_set(sfrmovie) then SEDM2_sfrmovie, dir_in, dir_out

  ;;-- SDSS image movie
  if keyword_set(sdssmovie) then begin
     if n_elements(imagesize) eq 0 then imagesize=1.
     if n_elements(movieorientation) eq 0 then movieorientation=0 ;face-on
     
     SEDM2_sdssmovie, dir_in, dir_out, redshift, tauv, mu_d,imagesize, movieorientation
  endif

  ;;-- optical spectra movies


  
END
