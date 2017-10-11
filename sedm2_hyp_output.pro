;+
; NAME:
;	SEDM2_HYP_OUTPUT
;
; PURPOSE:
;	This procedure reads in a Gadget file and outputs a fits file
;	containing the gas densities of the particles and their
;	positions. This can then be read into Hyperion for binning
;	onto a regular grid and further radiative transfer calculations 
;
;
; CALLING SEQUENCE:
;
;	SEDM2_HYP_OUPUT, Filename
;
;
; INPUTS:
;	Filename: name of input Gadget file. Output filename will be
;	          built from input filename.
;
;	
; KEYWORD PARAMETERS:
;	INDIR:	Directory where file is stored. Default is current
;	        directory. 
;
;	OUTDIR:	Directory where output file should be placed. Default
;	        is input directory. 
;
;
; EXAMPLE:
;       SEDM2_HYP_OUTPUT, 'Galaxy_model_Sc',                                             $
;                    INDIR='/Volumes/data/gadget_simulations/Galaxy_model_Sc/',     $
;                    OUTDIR='/Volumes/data/gadget_simulations/Galaxy_model_Sc/'
;
;
; DEPENDENCIES: 
;        sedm2_readsnap 
;        sedm2_readised.pro
;        sedm2_buildsed.pro
;		     
; 
; NOTES: 
;  The output structure contained in the fits file contains the
;  following parameters and units:
; 
;  x, y, z, position kpc
;  vx, vy, vz, velocity km/s
;  mass, mass  M_solar
;  sph_smooth, SPH smoothing length kpc
;  sfr, star formation rate  Msol/yr
;  id, original ID in Gadget file 
;  rho, density of gas M_sun /kpc^3 
;  nh, neutral gas fraction
;  ind_vel, placeholder for full code
;  
; MODIFICATION HISTORY:
; 	Written by:	Joe Llama April, 2014
;       Adapted by:     Vivienne Wild July 2016 from hyp_output.pro	
;-

PRO outputsources, dir, filename, source, name
fcode = '(f10.4)'
openw, lun, dir + filename + '_' + name + '.txt', /get_lun
for i = 0, n_elements(source.mass)-1 do $
   printf, lun, $
           string(source[i].x,format=fcode), ', ',  $
           string(source[i].y,format=fcode), ', ',  $
           string(source[i].z,format=fcode), ', ',  $
           string(source[i].vx,format=fcode), ', ', $
           string(source[i].vy,format=fcode), ', ', $
           string(source[i].vz,format=fcode), ', ', $
           string(source[i].mass,format='(I6.6)')
free_lun, lun
END
 

PRO SEDM2_HYP_OUTPUT, fileseq, models_dir,indir=indir, outdir=outdir

  if not (keyword_set(indir))  then cd, current=indir
  if not (keyword_set(outdir)) then outdir = indir

;make sure that indir has the '/' on the directory structure
  if (strmid(indir, 0, 1, /reverse_offset) ne '/') then $
     indir = indir+'/' 

  filename = file_search(indir+'*_???.hdf5',count=Nsnap) 

  for i=0,Nsnap-1 do $
     filename[i] = (strsplit(filename[i],'/',/extract,count=n))[n-1]
  print, 'SEDM2_HYP_OUTPUT:', string(Nsnap, format='(I03)') +' files found'
  
;;-- read in the SSPs: 1Msolar formed
  Z = '62'                      ;solar
  stem = models_dir + 'bc2003_hr_m'
  file_ssps = stem + Z + '_chab_ssp.ised'
  seds = SEDM_READISED(file_ssps, t, lambda)   
  age_ssp = t/1e9               ; in Gyr to match simns
  nlambda = n_elements(lambda)
  nage = n_elements(age_ssp)


  for i=0,Nsnap-1 do begin
     fname = (strsplit(filename[i],'.',/extract))[0]
     file_outgas   =  outdir + fname + '_gas.fits'
     file_outstars =  outdir + fname + '_sfh.fits'
     
     print, file_outgas, file_outstars

;;-- read Gadget file
     SEDM2_READSNAP, indir+filename[i], stars=stars, gas=gas, halo=halo,sfr=sfr, snap_time=snap_time,/getstars,/getgas,/gethalo

     if i eq 0 then oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value. 
     ind_oldstars = where(stars.id ge oldstars_minID,compl=ind_newstars, noldstars) ;ID the old stars
     if ind_newstars[0] ne -1 then nnewstars = n_elements(ind_newstars) else nnewstars = 0
     nstars = n_elements(stars)

;;!!! HERE!!! don't have potential, so can't get centers....

;;-- Run hyp_find_centers_merger to get the centers of the two galaxies;
;;-- NOTE: We have already run hyp_readsnap so call with halo specified 
;     SEDM2_FIND_CENTERS_MERGER, halo=halo, separation=separation

;;-- output the gas as a fits file 
     mwrfits, gas, file_outgas, /create

;;-- assign SSPs to old stars
     SEDM2_BUILDSED, age_ssp, stars, oldstars_minID, sfr, snap_time

;;-- Restore gassfh and starsfh
     restore, indir + fname + '_gassfh.sav'  

;;-- Create new array for all stars and add newstars SFH
     starsfh = fltarr(nstars,nage)
     if nnewstars ne 0 then starsfh[ind_newstars,*] = newstarsfh

;;-- Add old stars SFH (mass particles at single SSP ages)
     for j=0, nage-1 do begin
        ind = where(stars[ind_oldstars].ind_ssp eq j, count)
        if count eq 0  then continue
        starsfh[ind_oldstars[ind],j] = stars[ind_oldstars[ind]].mass
     endfor


;;-- output the SFHs for star particles
;;-- a table of nparticles x nage (221)
 
     mkhdr, hdr, starsfh, /extend
     ;; sxaddpar, hdr, 'X1_center',separation[0]
     ;; sxaddpar, hdr, 'Y1_center',separation[1]
     ;; sxaddpar, hdr, 'Z1_center',separation[2]
     ;; sxaddpar, hdr, 'X2_center',separation[3]
     ;; sxaddpar, hdr, 'Y2_center',separation[4]
     ;; sxaddpar, hdr, 'Z2_center',separation[5]
     ;; sxaddpar, hdr, 'Radial Distance',separation[6]


     mwrfits, starsfh,  file_outstars, hdr, /create
     mwrfits, gassfh,   file_outstars, /silent

     outputsources, outdir, fname, stars, 'stars'

  endfor


END


