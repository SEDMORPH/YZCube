PRO SEDM2_SPEC_STAR_AGE, dir_in, dir_out, tauv,mu_d, $
                snap = snap_in, model_str=model_str, $
                models_dir=dir_models, with_metal=with_metal

;;------------------------------------------------------------------
;; Parameters
;;------------------------------------------------------------------
;  ll_min = 3700.                ;min wavelength
;  ll_max = 9100.                ;max wavelength
;  with_metal   if being set, use the meatllicity given in sedm2_codeunits.inc
;               otherwise, use the a unifiy Metallicity (m62, i.e. Z=0.02 )for all stars,

;;------------------------------------------------------------------
;; Input and output files
;;------------------------------------------------------------------

  @sedm2_codeunits.inc

;;-- set up plotting file
  outstr = '_tauv'+string(tauv,form='(F0.1)')
  outstr = outstr+'_mu'+string(mu_d,form='(F0.1)')

  if KEYWORD_SET(with_metal) then begin
    if n_elements(snap_in) gt 0 then psfile = dir_out+'spectra'+outstr+'_'+string(snap_in,form='(I3.3)')+'_star_age_with_metal.ps' $
    else psfile = dir_out+'spectra'+outstr+'_star_age_with_metal.ps'
  endif else begin
    if n_elements(snap_in) gt 0 then psfile = dir_out+'spectra'+outstr+'_'+string(snap_in,form='(I3.3)')+'_star_age.ps' $
    else psfile = dir_out+'spectra'+outstr+'_star_age.ps'
  endelse


;;-- find snapshot filenames from gadget simulation files
  filename = file_search(dir_in+'/*.hdf5',count=nsnap)
  filename0 = filename[0]

  if n_elements(snap_in) gt 0 then begin

     SEDM2_READSNAP, filename[0], stars=stars, /getstars ;need to get minID of old stars from first file
     ; oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.
     oldstars_minid = -1 ;Treat all stars like old stars
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
if KEYWORD_SET(with_metal) then begin
  Z_keys = Z_models.keys
  ; read in the SSPs for the first metallicity value, record all the results
  ssps = SEDM2_GETSSPS(dir_models, model_str,Z_keys[0])
  age_ssp = ssps.age
  nssps = n_elements(age_ssp)
  lambda = ssps.lambda;[ind]
  nlambda = n_elements(lambda)
  ssps_lum_metal = fltarr(nlambda, nssps, N_z)
  ssps_lum_metal[*,*,0] = ssps.seds;[ind,*]
  ssps = 0

  ; read in the SSPs for the resest metallicity value, record the Luminosity only
  for mm=1, N_z-1 do begin
    ssps = SEDM2_GETSSPS(dir_models, model_str,Z_keys[mm])
    ssps_lum_metal[*,*,mm] = ssps.seds;[ind,*]
    ssps = 0
  endfor

endif else begin
  ;use the a unifiy Metallicity (m62, i.e. Z=0.02 )for all stars
  ssps = SEDM2_GETSSPS(dir_models, model_str,'62')
  age_ssp = ssps.age
  nssps = n_elements(age_ssp)
;  ind = where(ssps.lambda ge ll_min and ssps.lambda le ll_max,nlambda)
  lambda = ssps.lambda;[ind]
  nlambda = n_elements(lambda)
  ssps_lum = ssps.seds;[ind,*]
  ssps = 0
endelse

;;------------------------------------------------------------------
;;-- dust attenuation
;;------------------------------------------------------------------
  tau_young = mu_d*tauv*( (5500./lambda)^0.7) + (1-mu_d)*tauv*( (5500./lambda)^1.3)
  tau_old = mu_d*tauv*( (5500./lambda)^0.7)

;;------------------------------------------------------------------
;;-- Loop over all snapshots
;;------------------------------------------------------------------

  ps1c, psfile
  time = systime(1)

  mass_young = dblarr(nsnap)

  for i=0, nsnap-1 do begin ;Nsnap-1 do begin

     tmp = (strsplit(filename[i],'/',/extract,count=n))[n-1]
     filename_short = (strsplit(tmp,'.',/extract))[0]
     str_snap = (strsplit(filename_short,'_',/extract))[1] ;don't use i as could be only doing a single snapshot

     print, 'SEDM2_SPEC_STAR_AGE building spectrum for snapshot:'+str_snap

;;-- outfile (single snapshot, all components, all orientations)
     if KEYWORD_SET(with_metal) then begin
       outfile_fits = dir_out+'spec'+outstr+'_'+str_snap+'_star_age_with_metal.fits'
     endif else begin
       outfile_fits = dir_out+'spec'+outstr+'_'+str_snap+'_star_age.fits'
     endelse

;;-- output arrays for each component of the final image
     ;; keep the same structure as the SEDmorph method output,
     ;; but put all flux into spec_os_*, the other arrays (spec_g_* and spec_ns_*) will be empty
     ;; may consider to have a new structure in the future
     spec_g_old = (spec_g_young = (spec_ns_old = (spec_ns_young = (spec_os_old = (spec_os_young = fltarr(nlambda))))))

;;-- read simulation files  - 9 secs
     ; SEDM2_READSNAP, filename[i], stars=stars, gas=gas, sfr=sfr, snap_time=snap_time,/getstars,/getgas
     SEDM2_READSNAP, filename[i], stars=stars, sfr=sfr, snap_time=snap_time,/getstars



     if size(stars,/type) eq 8 then nstars = n_elements(stars) else nstars=0
     ; if size(gas,/type)   eq 8 then ngas   = n_elements(gas)   else ngas=0

     ; if n_elements(oldstars_minid) eq 0 and i eq 0 then oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.
     if n_elements(oldstars_minid) eq 0 and i eq 0 then oldstars_minid = -1 ;Treat all stars as old stars

     ind_oldstars = where(stars.id ge oldstars_minID,compl=ind_newstars, noldstars) ;ID the old stars
     if ind_newstars[0] ne -1 then nnewstars = n_elements(ind_newstars) else nnewstars = 0
     if nnewstars+noldstars ne nstars then stop
     print, i, nnewstars, noldstars, nstars

;;-- fill up ind_ssp star structures
     if i eq 0 then plot=1 else plot=0
     if noldstars gt 0 then SEDM2_BUILDSED, age_ssp, stars, oldstars_minid, snap_time,sfr=sfr, plot=plot

;;-- read gas particle and new star particle SFHs for this snapshot - 0.5 secs
     ; restore,  dir_in+'gassfh_uniform_Z_5_March/'+filename_short+'_gassfh.sav'


;;-- loop over SSPs to build integrated spectra
    if KEYWORD_SET(with_metal) then begin

      print, "Using metal"
      star_metal_bin = uintarr(noldstars)
      sedm2_z_ind, stars.metal, Z_models.values, Z_ind
      star_metal_bin = Z_ind

      accum_nstar_metal = 0 ;accumulative count of stars
      print, "Z_keys   number of the stars in this bins"
      for mm=0, N_z-1 do begin
        ind_metal = where(star_metal_bin eq mm, mcount)
        if mcount gt 0 then begin
            print, Z_keys[mm], mcount
            accum_nstar_metal += mcount

            mm_stars = stars[ind_metal]
            for j=0,nssps-1 do begin
              ind = where(mm_stars.ind_ssp eq j,nn2)
              ; print, nn2
              if ind[0] ne -1 then begin
                 if age_ssp[j] le 0.01 then spec_os_young = spec_os_young+total(mm_stars[ind].mass)*ssps_lum_metal[*,j,mm] $
                 else spec_os_old = spec_os_old+total(mm_stars[ind].mass)*ssps_lum_metal[*,j,mm]
              endif ; have particle in the ssp
            endfor

        endif ; have particles in the meatllicity bin

      endfor ; for loop of metallicity
      ;check that all stars have been counted
      if accum_nstar_metal ne nstars then  begin
        print, "Not all stars are counted for the spectra, stop here"
        stop
      endif else begin
        print, "nstars          total star number in all metallicity bins"
        print, nstars, accum_nstar_metal
        print, "All stars have been counted."
      endelse

    endif else begin ; loop over meatllicity bins

      print, "Uni_metal"
      for j=0,nssps-1 do begin

        ; if noldstars gt 0 then begin
        ind = where(stars.ind_ssp eq j,nn2)
        if ind[0] ne -1 then begin
          if age_ssp[j] le 0.01 then spec_os_young = spec_os_young+total(stars[ind].mass)*ssps_lum[*,j] $
          else spec_os_old = spec_os_old+total(stars[ind].mass)*ssps_lum[*,j]
        endif
        ; endif
      endfor
    endelse ; use meatllicity


;;-- sum components and add dust
     spec_ns = spec_ns_young+spec_ns_old
     spec_ns_dust = spec_ns_young*exp(-tau_young)+spec_ns_old*exp(-tau_old)
     spec_g = spec_g_young+spec_g_old
     spec_g_dust = spec_g_young*exp(-tau_young)+spec_g_old*exp(-tau_old)
     spec_os = spec_os_young+spec_os_old
     spec_os_dust = spec_os_young*exp(-tau_young)+spec_os_old*exp(-tau_old)

     spec_notau = spec_ns+spec_g+spec_os
     spec_tau = spec_ns_dust+spec_g_dust+spec_os_dust

;;-- PCA components


;;-- save output file
     struc = {wave:lambda, spec_tau:spec_tau,spec_notau:spec_notau,spec_stars:spec_ns, spec_bulge:spec_os, spec_gas:spec_g}
     mwrfits, struc, outfile_fits,/create


;;------------------------------------------------------------------
;;-- some figures
;;------------------------------------------------------------------

     plot, lambda,spec_notau[*,0],/xlog,/ylog,/xs,xtitle='Wavelength [A]',ytitle=textoidl('Luminosity [L_\odot/A]'),title='Orien 0: Snapshot '+str_snap;+' N*='+string(nstars,form='(I0)')
     oplot, lambda,spec_os[*,0],color=cgcolor('red')
     oplot, lambda,spec_ns[*,0],color=cgcolor('cyan')
     oplot, lambda,spec_g[*,0],color=cgcolor('purple')

     xyouts, 0.8,0.9,'Total',/normal
     xyouts, 0.8,0.85,'Old stars',/normal,color=cgcolor('blue')
     xyouts, 0.8,0.75,'New Stars',/normal,color=cgcolor('cyan')
     xyouts, 0.8,0.7,'Gas',/normal,color=cgcolor('purple')

     plot, lambda,spec_tau[*,0],/xlog,/ylog,/xs,xtitle='Wavelength [A]',ytitle=textoidl('Luminosity [L_\odot/A]'),title='Orien 0: Snapshot '+str_snap+' with dust';+' N*='+string(nstars,form='(I0)')
     oplot, lambda,spec_os_dust[*,0],color=cgcolor('red')
     oplot, lambda,spec_ns_dust[*,0],color=cgcolor('cyan')
     oplot, lambda,spec_g_dust[*,0],color=cgcolor('purple')

     plot, lambda, -alog(spec_tau[*,0]/spec_notau[*,0]),/xlog,/xs,ytitle=textoidl('\tau_\lambda'),xtitle='Wavelength [A]'
     oplot, lambda, -alog(spec_g_dust[*,0]/spec_g[*,0]),color=cgcolor('purple')
     oplot, lambda, -alog(spec_ns_dust[*,0]/spec_ns[*,0]),color=cgcolor('cyan')
     oplot, lambda, -alog(spec_os_dust[*,0]/spec_os[*,0]),color=cgcolor('red')


     ;; shorteer wavelength range
     plot, lambda,spec_notau[*,0],/xs,xr=[3000,9000],xtitle='Wavelength [A]',ytitle=textoidl('Luminosity [L_\odot/A]'),title='Orien 0: Snapshot '+str_snap+' N*='+string(nstars,form='(I0)')
     oplot, lambda,spec_os[*,0],color=cgcolor('red')
     oplot, lambda,spec_ns[*,0],color=cgcolor('cyan')
     oplot, lambda,spec_g[*,0],color=cgcolor('purple')
     xyouts, 0.8,0.4,'Total',/normal
     xyouts, 0.8,0.35,'Old stars',/normal,color=cgcolor('red')
     xyouts, 0.8,0.25,'New Stars',/normal,color=cgcolor('cyan')
     xyouts, 0.8,0.2,'Gas',/normal,color=cgcolor('purple')


     plot, lambda,spec_tau[*,0],/xs,xr=[3000,9000],title='Orien 0: Snapshot '+str_snap+' with dust'+' N*='+string(nstars,form='(I0)')
     oplot, lambda,spec_os_dust[*,0],color=cgcolor('red')
     oplot, lambda,spec_ns_dust[*,0],color=cgcolor('cyan')
     oplot, lambda,spec_g_dust[*,0],color=cgcolor('purple')
     xyouts, 0.8,0.9,'Total',/normal
     xyouts, 0.8,0.85,'Old stars',/normal,color=cgcolor('red')
     xyouts, 0.8,0.75,'New Stars',/normal,color=cgcolor('cyan')
     xyouts, 0.8,0.7,'Gas',/normal,color=cgcolor('purple')

     plot, lambda, -alog(spec_tau[*,0]/spec_notau[*,0]),/xs,xr=[3000,9000],ytitle=textoidl('\tau_\lambda')
     oplot, [5500,5500],[0,5],linestyle=2
     oplot, lambda, -alog(spec_g_dust[*,0]/spec_g[*,0]),color=cgcolor('purple')
     oplot, lambda, -alog(spec_ns_dust[*,0]/spec_ns[*,0]),color=cgcolor('cyan')
     oplot, lambda, -alog(spec_os_dust[*,0]/spec_os[*,0]),color=cgcolor('red')


  endfor


END
