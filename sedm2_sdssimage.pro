;+
; NAME:
;	SEDM2_SDSSIMAGE
;
; PURPOSE:
;	Builds sdss images including realistic sky, noise and psf
;
;
; CALLING SEQUENCE:
;
;
;
; OPTIONAL INPUT:
;
;
; KEYWORD PARAMETERS:
;
;
; EXAMPLE:
;
;
;
; DEPENDENCIES:
;               sedm2_getssps -> sedm2_readised
;               vw_func_rotate
;               sedm2_conv_filters
;               sedm2_sdssgrid
;               sedm2_readsnap
;
; NOTES:
;
;
; MODIFICATION HISTORY:
; 	Written by:	Vivienne Wild July 2016, based on SEDM_MOCKSDSSIMAGE
;
;-





;;------------------------------------------------------------------

function magtocounts, mag_AB, p
common image, b,aa,kk,air, gain, exptime

;; f/f0 = counts/exptime * 10^0.4*(aa + kk * airmass)
f_over_f0 = (2*b[p])*sinh(-mag_AB*(alog(10)/2.5)-alog(b[p]))

counts = exptime*f_over_f0*10d^(-0.4*(aa+kk*air)[p])

return,counts

end

;;------------------------------------------------------------------

function countstomag,counts, p
common image, b,aa,kk,air, gain, exptime

f_over_f0 = counts*10d^(0.4*(aa+kk*air)[p])/exptime
mag_AB = -(asinh(f_over_f0 / (2*b[p]))+alog(b[p]))/(alog(10)/2.5)

return, mag_AB

end

;;------------------------------------------------------------------

function countstojansky, counts,p
common image, b,aa,kk,air, gain, exptime

f_over_f0 = counts*10d^(0.4*(aa+kk*air)[p])/exptime
flux = 3631*f_over_f0

return, flux

end

;;------------------------------------------------------------------

FUNCTION JANSKYTOMAG, f_nu

;; f_nu_rest in Jansky
;; mag is AB
;; 1Jy = 10-23 erg s-1 Hz-1 cm-2
;f_nu_ergs = f_nu/1e23
;mag_AB = -alog10(f_nu_ergs)*2.5 - 48.6

mag_AB = -2.5*alog10(f_nu/3631.)

return, mag_AB


END



;;******************************************************************
;;******************************************************************


PRO SEDM2_SDSSIMAGE, dir_in, dir_out, dir_filters, dir_code, filterlist, redshift, tauv, mu_d,imagesize,faceon=faceon,$
                snap=snap_in,model_str=model_str,models_dir=models_dir, rotate_seed=rotate_seed



  @sedm2_codeunits.inc

;;------------------------------------------------------------------
;; Image parameters - stored in common block
;;------------------------------------------------------------------

  common image, b,aa,kk,air,gain, exptime

;;-- CCD pixel size in arcsec
  pixel_size = 0.396

;;-- PSF parameters (Milena)
;;-- Values from SDSS to construct a 2-gaussian psf (median of values for 500000 SDSS fields)
  psf_s12g = [1.44527,1.36313,1.21196,1.114143,1.20288] ;; in pixels
  psf_s22g = [3.1865,3.06772,2.82102,2.75094,3.09981] ;; in pixels
  psf_ratio =[0.081989,0.081099,0.06811,0.059653,0.054539]
  psf_width = [1.54254,1.44980,1.31878,1.25789,1.29391]/pixel_size ;; originally in arcsec

  b = [1.4,0.9,1.2,1.8,7.4]*10d^(-10)

;;-- median SDSS values from Jairo
  sky_level_pixel = janskytomag([ 1.20542,1.67145,  3.99916, 7.07624,  20.7627]*(3.631*10d^(-6)))+2.5*alog10(1/pixel_size^2)    ;nanomaggy -> jansky
  aa = [ -24.1270,  -24.5128,  -24.1073,  -23.6936,   -21.9465]
  kk = [ 0.506005, 0.181950,  0.103586,  0.0628950,  0.0557740]
  air = [ 1.19020,1.19027, 1.19027, 1.19022,1.19022]
  gain = [1.6,3.995,4.725,4.885,4.775]
  darkvar = [9.30250,1.44000,1.00000, 5.76000, 1.00000 ]

;;-- Values in magnitudes, converted from SDSS skyerr (Milena):
  skysig = [0.0047951106042681858, 0.0015403498079674194, 0.00098253076542588599, 0.00088190450872838123, 0.0010950415323111404]

;;-- Exposure time (seconds??)
  exptime = 53.907456

;;------------------------------------------------------------------
;; Input and output files
;;------------------------------------------------------------------

;;-- set up plotting file
  outstr = '_z'+string(redshift,form='(F0.3)')
  outstr = outstr+'_tauv'+string(tauv,form='(F0.1)')
  outstr = outstr+'_mu'+string(mu_d,form='(F0.1)')
  outstr = outstr+'_'+string(imagesize,form='(I0)')+'arcmin'

  if n_elements(snap_in) gt 0 then psfile = dir_out+'sdssimage'+outstr+'_'+string(snap_in,form='(I3.3)')+'.ps' $
  else psfile = dir_out+'sdssimage'+outstr+'.ps'

;;-- find snapshot filenames from gadget simulation files
  filename = file_search(dir_in+'/*.hdf5',count=nsnap)
  filename0 = filename[0]

  if n_elements(snap_in) gt 0 then begin

     SEDM2_READSNAP, filename[0], stars=stars, /getstars ;need to get minID of old stars from first file
     oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.

     tmp = where(strmatch(filename, '*'+string(snap_in,form='(I3.3)')+'.hdf5') eq 1)
     if tmp[0] eq -1 then message, 'Something wrong with identifying snapshot file: ', snap_in
     filename = filename[tmp]

  endif


  nsnap = n_elements(filename)

  print, 'SEDM2_SDSSIMAGE: number of snapshots to process',nsnap

;;------------------------------------------------------------------
;; Pre-computations
;;------------------------------------------------------------------

;;-- read SSPs !!!! NEED TO ADD DIFFERENT METALLICITIES
  ssps = SEDM2_GETSSPS(models_dir, model_str,'62')
  age_ssp = ssps.age
  nssps = n_elements(age_ssp)

;;-- convolve SSPs with filter functions
  readcol, dir_filters+filterlist, filters,ll_eff,form='(A,F)' ;list of filter files
  nband = n_elements(filters)
  ssps_lum = fltarr(nssps,nband)
  band_name = ['u','g','r','i','z']

  for i=0,nband-1 do begin
     readcol, dir_filters+filters[i],lam_fil,fil,/silent
     for j=0,nssps-1 do ssps_lum[j,i] = SEDM2_CONV_FILTER(ssps.lambda,ssps.seds[*,j],lam_fil,fil,redshift)
  endfor

;;-- set up 6 orientations + face-on
; leave for now - just do face-on
  if keyword_set(faceon) then begin
     rotate = VW_FUNC_ROTATE([0,0,0])
     norien = 1
  endif else rotate = SEDM2_ROTATE(fileseq,norien,seed=rotate_seed) ;norien is output

;;-- decide the nx and ny note the grid!!! The grids need to generate for every orientation seperately
  SEDM2_SDSSGRID, xx,yy,nx,ny,redshift,imagesize=imagesize
  

;;-- distance to the galaxy
  dist = lumdist(redshift,H0=hubparam*100,omega_m=omega_m,lambda=omega_l,/silent)*Mpc_in_cm

;;------------------------------------------------------------------
;;-- dust attenuation in each of the filters
;;------------------------------------------------------------------
;; Note the approximation here. We should attenuation the SSPs, then
;; apply the filter functions. But instead we attenuate the
;; photometry.
  tau_young = mu_d*tauv*( (5500./ll_eff)^0.7) + (1-mu_d)*tauv*( (5500./ll_eff)^1.3)
  tau_old = mu_d*tauv*( (5500./ll_eff)^0.7)

;;------------------------------------------------------------------
;;-- Loop over all snapshots
;;------------------------------------------------------------------

  ps1c, psfile
  time = systime(1)

  mass_young = dblarr(nsnap,norien)
  
  ;;------------yrzheng, for a moving box
  ;; read in the centers list for all snapshots
  temp =  file_search(dir_in+'*_???.hdf5',count=tot_nsnap)
  centerlist = fltarr(6, tot_nsnap)

  openr, centertxt, dir_in+'centers.txt', /get_lun
  readf, centertxt, centerlist
  free_lun, centertxt


  for i=0, nsnap-1 do begin ;Nsnap-1 do begin
     tmp = (strsplit(filename[i],'/',/extract,count=n))[n-1]
     filename_short = (strsplit(tmp,'.',/extract))[0]
     str_snap = (strsplit(filename_short,'_',/extract))[1] ;don't use i as could be only doing a single snapshot

     print, 'SEDM2_SDSSIMAGE building image for snapshot:'+str_snap

;;-- outfile (single snapshot, all components, all orientations)
;     outfile = dir_out+'sdssimage'+outstr+'_'+str_snap+'.sav'
     outfile_fits = dir_out+'sdssimage'+outstr+'_'+str_snap+'.fits'

;;-- output arrays for each component of the final image
     image_jansky = (image_conv = (image_noise = (image_noise_counts = fltarr(nx,ny,norien,nband))))
     image_jansky_dust = (image_conv_dust = (image_noise_dust = (image_noise_counts_dust = fltarr(nx,ny,norien,nband))))

;;-- read simulation files  - 9 secs
     SEDM2_READSNAP, filename[i], stars=stars, gas=gas, sfr=sfr, snap_time=snap_time,/getstars,/getgas

;     ntmp = n_elements(gas)
;     gas = gas[round(randomu(234985,ntmp/10.)*(ntmp-1))]
;     stars = stars[round(randomu(234985,ntmp/10.)*(ntmp-1))]

     if size(stars,/type) eq 8 then nstars = n_elements(stars) else nstars=0
     if size(gas,/type)   eq 8 then ngas   = n_elements(gas)   else ngas=0

     if n_elements(oldstars_minid) eq 0 and i eq 0 then oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.

     ind_oldstars = where(stars.id ge oldstars_minID,compl=ind_newstars, noldstars) ;ID the old stars
     if ind_newstars[0] ne -1 then nnewstars = n_elements(ind_newstars) else nnewstars = 0
     if nnewstars+noldstars ne nstars then stop

;;-- fill up ind_ssp star structures
     if i eq 0 then plot=1 else plot=0
     if noldstars gt 0 then SEDM2_BUILDSED, age_ssp, stars, oldstars_minid, sfr, snap_time,plot=plot

;;-- read gas particle and new star particle SFHs for this snapshot - 0.5 secs
     restore,  dir_in+filename_short+'_gassfh.sav'

;     gassfh = gassfh[round(randomu(234985,ntmp/10.)*(ntmp-1)),*]

;;------------------------------------------------------------------
;;-- loop over orientations
;;------------------------------------------------------------------
     for j=0,norien-1 do begin

;;------------------------------------------------------------------
;;-- loop through each pixel and sum luminosities in SSPs to create
;;   images in jansky. ~10 secs
;;------------------------------------------------------------------
        time = systime(1)
        
        ;;---use the center of the first halo as the center of the grid
        snapnum = round(float(str_snap))
        two_center = centerlist[*, snapnum]
        center = two_center[0:2]
        SEDM2_SDSSGRID, xx,yy,nx,ny,redshift, imagesize=imagesize, center =  center, rotate = rotate[*,*,j]
        print, 'grid details:', 'z=',redshift, 'starting point=',min(xx),'nx=',nx

        SEDM2_ZSUM, nx, ny, xx, yy, gas, stars, ind_newstars, ind_oldstars, gassfh, newstarsfh, rotate[*,*,j], dist, $ ;input
                    ssps_lum, age_ssp, ll_eff, tau_young, tau_old, $                             ;SSPs
                    flux_all_jansky, flux_all_jansky_dust, tmp                                        ;output

        mass_young[i,j] = tmp
        image_jansky[*,*,j,*] = flux_all_jansky
        image_jansky_dust[*,*,j,*] = flux_all_jansky_dust

        print, systime(1)-time
;;------------------------------------------------------------------
;;-- convolve with PSF 0.02sec
;;------------------------------------------------------------------
        for p=0,nband-1 do begin

           ;;Milena:
           ;; Accounting for the different amplitudes of the gaussian components:
           kernel1 = gaussian_function(replicate(psf_s12g[p],2),width=psf_width[p],maximum=1.0)
           kernel1 = kernel1/total(kernel1)
           kernel2 = gaussian_function(replicate(psf_s22g[p],2),width=psf_width[p],maximum=psf_ratio[p])
           kernel2 = kernel2/total(kernel2)
           kernel = (kernel1+kernel2)/total(kernel1+kernel2)

           image_conv[*,*,j,p] = convol(reform(image_jansky[*,*,j,p]),kernel)
           image_conv_dust[*,*,j,p] = convol(reform(image_jansky_dust[*,*,j,p]),kernel)

        endfor

;;------------------------------------------------------------------
;;-- add noise 0.02sec
;;------------------------------------------------------------------
        for p=0,nband-1 do begin

           mag_AB = JANSKYTOMAG(image_conv[*,*,j,p])                                  ;convert to AB magnitudes
           mag_AB_dust = JANSKYTOMAG(image_conv_dust[*,*,j,p])                                  ;convert to AB magnitudes

          if band_name[p] eq 'u' then begin
              mag_AB = mag_AB+0.04 ;add sdss zeropoint offsets
              mag_AB_dust = mag_AB_dust+0.04
           endif
           if band_name[p] eq 'z' then begin
              mag_AB = mag_AB-0.02
              mag_AB_dust = mag_AB_dust-0.02
           endif
           image_counts = magtocounts(mag_AB,p) > 0 ;counts/pixel: force to be >= zero
           image_counts_dust = magtocounts(mag_AB_dust,p) > 0 ;counts/pixel: force to be >= zero

           sky_counts =  magtocounts(sky_level_pixel[p],p)

            ;; error(counts) = sqrt([counts+sky]/gain + Npix*(dark_variance+skyErr))
           error_sigma = sqrt((image_counts+sky_counts)/gain[p]+1*(darkvar[p]+skySig[p]))
           error_sigma_dust = sqrt((image_counts_dust+sky_counts)/gain[p]+1*(darkvar[p]+skySig[p]))

           dim1 = (size(image_counts,/dim))[0]
           dim2 = (size(image_counts,/dim))[0]
           error_in_counts = randomn(seed,dim1,dim2) * error_sigma ;Gaussian random noise with sigma = error in counts
           error_in_counts_dust = randomn(seed,dim1,dim2) * error_sigma_dust ;Gaussian random noise with sigma = error in counts

           error_in_ff0 = (error_in_counts/exptime) * 10d^(0.4*(aa+kk*air)[p])
           error_in_jansky = 3631*error_in_ff0 ; error in flux : S = 3631*f/f0

           error_in_ff0_dust = (error_in_counts_dust/exptime) * 10d^(0.4*(aa+kk*air)[p])
           error_in_jansky_dust = 3631*error_in_ff0_dust ; error in flux : S = 3631*f/f0

           image_noise[*,*,j,p] = image_conv[*,*,j,p] + error_in_jansky
           image_noise_counts[*,*,j,p] = image_counts + error_in_counts
           image_noise_dust[*,*,j,p] = image_conv_dust[*,*,j,p] + error_in_jansky_dust
           image_noise_counts_dust[*,*,j,p] = image_counts_dust + error_in_counts_dust


        endfor


;;------------------------------------------------------------------
;;-- plots
;;------------------------------------------------------------------
        if j eq 0 then begin
           !p.multi=[0,2,3]


           ;; output from simulation, in jansky, r-band
           CGimage, bytscl(alog10(image_jansky[*,*,j,2]),/nan), /nointerp,multimargin=1,/keep_aspect_ratio

           ;; convolved with PSF, in jansky, r-band
           CGimage, bytscl(alog10(image_conv[*,*,j,2]),/nan), /nointerp,multimargin=1,/keep_aspect_ratio

           ;; add noise, in jansky, r-band
           CGimage, bytscl(alog10(image_noise[*,*,j,2]),/nan), /nointerp,multimargin=1,/keep_aspect_ratio


          if i eq 0 then begin

;;--  compare to real galaxy
;; can use this: img = mrdfits(imgname,/fscale)
;; image 1 has z=0.013

             ;;origin band setting function
             ;;gband = float((vw_readfits(dir_code+'sdssimages/stamp_1g.fits',hdr,/silent)-1000))            ;; counts = value - softbias ; this includes the sky background
             ;;rband = float((vw_readfits(dir_code+'sdssimages/stamp_1r.fits',hdr,/silent)-1000))
             ;;iband = float((vw_readfits(dir_code+'sdssimages/stamp_1i.fits',hdr,/silent)-1000))
             ;;-------------------

             ;;my bands ------------
             gband = float((mrdfits(dir_code+'sdssimages/stamp_1g.fits',/fscale)-1000))            ;; counts = value - softbias ; this includes the sky background
             rband = float((mrdfits(dir_code+'sdssimages/stamp_1r.fits',/fscale)-1000))
             iband = float((mrdfits(dir_code+'sdssimages/stamp_1i.fits',/fscale)-1000))
             ;;----------------------

              ;; it's crucial to get the sky level correct for the colour images
                                ;sky_r = 160 & sky_g = 81 & sky_i = 195  ;; image 6
              sky_r = 131 & sky_g = 72 & sky_i = 218  ;; image 1

              plot_histogram,image_noise_counts[*,*,j,2],nbins=100,min = min(rband),max=max(rband),xtitle='counts',title='dotted=real, full=mock'
              plot_histogram,rband,nbins=100,/over,linestyle=1



              ;;set minmax to be the same as the real image for better comparison
              CGimage, bytscl(alog10(image_noise_counts[*,*,j,2]),/nan,max = max(alog10(rband),/nan),min=min(alog10(rband),/nan)), /nointerp,multimargin=1,/keep_aspect_ratio
              CGimage, bytscl(alog10(rband),/nan,max = max(alog10(rband),/nan),min=min(alog10(rband),/nan)), /nointerp,multimargin=1,/keep_aspect_ratio

;;-- 3 colour comparison
              dim1 = (size(rband,/dim))[0]
              dim2 = (size(rband,/dim))[1]
              rgbimage = fltarr(dim1,dim2,3)
              rgbimage[*,*,2] = gband-sky_g
              rgbimage[*,*,1] = rband-sky_r
              rgbimage[*,*,0] = iband-sky_i

              tmp = [3,2,1]
;              rgbmock = reform(image_noise_counts[*,*,j,tmp]) ;this does not have sky added (just noise from sky)

              rgbimage_jans = rgbimage
              rgbimage_jans[*,*,0] = countstojansky(rgbimage[*,*,0],3)
              rgbimage_jans[*,*,1] = countstojansky(rgbimage[*,*,1],2)
              rgbimage_jans[*,*,2] = countstojansky(rgbimage[*,*,2],1)


              ;SEDM_RGBIM,  rgbmock_jans,[0,1,2]; these look very similar, but the counts version looses the very faint features at the edge, where we've set counts=1 above
              SEDM2_RGBIM,  reform(image_noise[*,*,j,*]),[3,2,1]
              SEDM2_RGBIM,  reform(image_noise_dust[*,*,j,*]),[3,2,1]
              SEDM2_RGBIM, rgbimage_jans, [0,1,2]

           endif else begin
              ;; colour image
              CGimage, bytscl(alog10(image_noise_counts[*,*,j,2]),/nan,max = max(alog10(rband),/nan),min=min(alog10(rband),/nan)), /nointerp,multimargin=1,/keep_aspect_ratio ;still scaled to rband image
              SEDM2_RGBIM,  reform(image_noise[*,*,j,*]),[3,2,1]
              SEDM2_RGBIM,  reform(image_noise_dust[*,*,j,*]),[3,2,1]
           endelse

        endif

     endfor                     ;orientation


;;------------------------------------------------------------------
;;-- output
;;------------------------------------------------------------------
     mwrfits, image_noise, outfile_fits,/create
     mwrfits, image_noise_dust, outfile_fits,/silent ;goes in next HDU


  endfor                        ;snapshot

  plot, mass_young[*,0],ytitle='SFR'

  ps2




END
