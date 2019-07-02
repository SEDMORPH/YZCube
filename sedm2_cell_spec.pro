
PRO SEDM2_CELL_SPEC, dir_in, dir_out, tauv,mu_d,redshift, cell_x_offset, cell_y_offset,cell_size, arcsec=arcsec,$
                snap = snap_in, style=style, model_str=model_str,$
                models_dir=dir_models, rtfaceon=rtfaceon, one_comp_dust=one_comp_dust, $
                with_metal=with_metal, cir_fib=cir_fib
;+
; create spectra for the cell descibed below:
;
; cell_x_offset: the offset on x-axis from the center of the first galaxy
; cell_y_offset: the offset on y-axis from the center of the first galaxy
; cell_size    : the size of the cell, side length of a square fiber or diameter of a circular fiber
; arcsec       : if this keyword is set, the parameters above are in arcsec.
;                Otherwise they are in kpc. Use kpc by default
;                We do not suggest use arcsec!!! As we want keep the spectra in rest-frame
; rtfaceon     : rotate to faceon, if switched on, rotate the disk to be a face on one.
; style        : the spectra style, choose from the following choices
;                  "" (empty string) --> SEDmorph method
;                  "star_age" --> star_age method
;                   also support _eagle and _eagle_minus, but these are not well tested yet.(20-Aug-2018)
; one_comp_dust: use tau_old for all stars, i.e. tau_young = tau_old
; cir_fib      : circular fiber
;_


;;------------------------------------------------------------------
;; Parameters
;;------------------------------------------------------------------
;  ll_min = 3700.                ;min wavelength
;  ll_max = 9100.                ;max wavelength

;;------------------------------------------------------------------
;; Input and output files
;;------------------------------------------------------------------
  @sedm2_codeunits.inc

  if NOT KEYWORD_SET(style) then style='' ;;SEDMoprh style by default
  ;; use add an underscore for non-SEDmorph style, so that we use file_style for the input and output files
  if (strlowcase(style) eq "sedmorph") || (style eq '') then begin
    style=''
    file_style=''
  endif else file_style='_'+style
;;-- set the output file name
  outstr = '_tauv'+string(tauv,form='(F0.1)')
  outstr = outstr+'_mu'+string(mu_d,form='(F0.1)')
  if KEYWORD_SET(rtfaceon) then outstr=outstr+'_fo'
  if KEYWORD_SET(one_comp_dust) then outstr=outstr+'_one_comp_dust'
;;-- set up plotting file
  cell_str = 'cell_'+string(cell_x_offset, form='(F+0.2)')+string(cell_y_offset, form='(F+0.2)')
  if KEYWORD_SET(arcsec) then begin
      ;; convert the parameters value from arcsec to kpc
      print, "Check the cell parameters in aresec"
      print, "x_offset | y_offset | cell_size"
      print, cell_x_offset, cell_y_offset,cell_size
      ; converting
      kpc_in_arcsec = angdiam_fib(redshift,angsize=1.0)
      cell_x_offset = cell_x_offset * kpc_in_arcsec
      cell_y_offset = cell_y_offset * kpc_in_arcsec
      cell_size = cell_size * kpc_in_arcsec
      print, "NOTE! The spectra is still in rest-frame."
      print, "NOTE! The spectra is still in rest-frame."
      print, "NOTE! The spectra is still in rest-frame."
      ; cell_str = cell_str+'_in_arcsec_z'+string(redshift,form='(F0.3)')
  endif

  if KEYWORD_SET(cir_fib) then begin
    cell_str = cell_str+'_cir_size_'+string(cell_size, form='(F0.2)' )
  endif else begin
    cell_str = cell_str+'_size_'+string(cell_size, form='(F0.2)' )
  endelse

  print, "Check the cell parameters in kpc"
  print, "x_offset | y_offset | cell_size"
  print, cell_x_offset, cell_y_offset,cell_size



  if n_elements(snap_in) le 0 then psfile = dir_out+cell_str+'_spectra'+outstr+string(file_style)+'.ps'

;;-- find snapshot filenames from gadget simulation files
  filename = file_search(dir_in+'/*.hdf5',count=nsnap)
  filename0 = filename[0]

  if n_elements(snap_in) gt 0 then begin

     SEDM2_READSNAP, filename[0], stars=stars, /getstars ;need to get minID of old stars from first file
     oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.
     ; for star_age method, treat all stars as old stars, so oldstars_minid should be set to -1
     if style eq 'star_age' then oldstars_minid=-1

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
    if style ne 'star_age' then begin
      print, "Meatllicity is only work for star_age method currently! --- 12-Mar-19"
      stop
    endif
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
  if KEYWORD_SET(one_comp_dust) then tau_young = tau_old

;;------------------------------------------------------------------
;;-- Loop over all snapshots
;;------------------------------------------------------------------

  if n_elements(snap_in) le 0 then ps1c, psfile
  time = systime(1)

  mass_young = dblarr(nsnap)

  for i=0, nsnap-1 do begin ;Nsnap-1 do begin

     tmp = (strsplit(filename[i],'/',/extract,count=n))[n-1]
     filename_short = (strsplit(tmp,'.',/extract))[0]
     str_snap = (strsplit(filename_short,'_',/extract))[1] ;don't use i as could be only doing a single snapshot

     print, 'SEDM2_SPEC building spectrum for snapshot:'+str_snap
 ;;-- make a directory for each data cube.
     if KEYWORD_SET(cir_fib) then begin
       data_cube_dir = "DataCube"+outstr+'_'+str_snap+string(file_style)+'_cir_size_'+string(cell_size, form='(F0.2)' )
     endif else begin
       data_cube_dir = "DataCube"+outstr+'_'+str_snap+string(file_style)+'_size_'+string(cell_size, form='(F0.2)' )
     endelse

     if KEYWORD_SET(with_metal) then data_cube_dir=data_cube_dir+'_with_metal'
     data_cube_dir = data_cube_dir + '/'
     if file_test(dir_out+data_cube_dir) eq 0 then spawn, 'mkdir '+ dir_out+data_cube_dir
     if n_elements(snap_in) gt 0 then psfile = dir_out+data_cube_dir+cell_str+'.ps' & ps1c, psfile
;;-- outfile (single snapshot, all components, all orientations)
     outfile_fits = dir_out+ data_cube_dir+'spec_'+cell_str+'.fits'
     print, outfile_fits
     print, psfile
     ; wait, 10
;;-- output arrays for each component of the final image
     spec_g_old = (spec_g_young = (spec_ns_old = (spec_ns_young = (spec_os_old = (spec_os_young = fltarr(nlambda))))))

;;-- read simulation files  - 9 secs
     SEDM2_READSNAP, filename[i], stars=stars, gas=gas, sfr=sfr, snap_time=snap_time,/getstars,/getgas


     if n_elements(oldstars_minid) eq 0 and i eq 0 then oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.
     if style eq 'star_age' then oldstars_minid=-1
 ;;-- select particles in the center part
      temp =  file_search(dir_in+'*_???.hdf5',count=tot_nsnap)
      centerlist = fltarr(6, tot_nsnap)
      ;;read the centers
      openr, centertxt, dir_in+'centers.txt', /get_lun
      readf, centertxt, centerlist
      free_lun, centertxt
      two_center = centerlist[*,uint(str_snap)]
      center = two_center[0:2]

      ;;rotate to faceon
      if KEYWORD_SET(rtfaceon) then begin
        ; outstr=outstr+'_fo' ;;fo for Face on
        theta1=0.
        phi1 = 0.
        ;; Calculate rotation matrix
        fileseq_ex = (strsplit(dir_in,'/',/extract))[-1] ;;ex for EXtracted
        ; print, fileseq_ex
        orbit = strmid(fileseq_ex, strlen(fileseq_ex)-2, 2 )
        if orbit eq '07' then theta1 = -109 & phi1 = -60
        pi = !CONST.pi
        theta1 *=pi/180
        phi1 *= pi/180
        rotate = vw_func_rotate([theta1,0,0]) ## vw_func_rotate([0,0, phi1])

        ;;do the rotation
        ; print, n_elements(gas.x)
        xyz_gas = transpose(rotate ## [[gas.x-center[0]],[gas.y-center[1]],[gas.z-center[2]]])

        gas.x = reform(xyz_gas[0,*]) + center[0]
        gas.y = reform(xyz_gas[1,*]) + center[1]
        gas.z = reform(xyz_gas[2,*]) + center[2]
        ; print, n_elements(gas.x)

        xyz_stars = transpose(rotate ## [[stars.x-center[0]],[stars.y-center[1]],[stars.z-center[2]]])
        stars.x = reform(xyz_stars[0,*]) + center[0]
        stars.y = reform(xyz_stars[1,*]) + center[1]
        stars.z = reform(xyz_stars[0,*]) + center[2]
        ; print, "check gas.x[0], after rotation", gas.x[0][0]

      endif



      center[0] += cell_x_offset
      center[1] += cell_y_offset
      print, 'oldstars_minid', oldstars_minid

      print,"cell_size", cell_size
      ;print, stars[0]
      ind_newstars = where(stars.id lt oldstars_minid, nnewstars,compl=ind_oldstars)
      newstars = stars[ind_newstars]
      print, "Check newstars:"
      get_center_ind, newstars, cen_part_ind, center=center, box_size=cell_size, cir_fib=cir_fib
      newstar_cen_ind = cen_part_ind
      print, "Check stars:"
      get_center_ind, stars, cen_part_ind, center=center, box_size=cell_size, cir_fib=cir_fib
      ;print, cen_part_ind
      stars = stars[cen_part_ind]
      print, "Check gas:"
      get_center_ind, gas, cen_part_ind, center=center, box_size=cell_size, cir_fib=cir_fib
      gas_cen_ind = cen_part_ind
      gas = gas[gas_cen_ind]
      ; print, n_elements(gas)



      if size(stars,/type) eq 8 then nstars = n_elements(stars) else nstars=0
      if size(gas,/type)   eq 8 then ngas   = n_elements(gas)   else ngas=0

     ind_newstars = where(stars.id lt oldstars_minid, nnewstars)
     ind_oldstars = where(stars.id ge oldstars_minID,compl=ind_newstars, noldstars) ;ID the old stars
     if ind_newstars[0] ne -1 then nnewstars = n_elements(ind_newstars) else nnewstars = 0
     if nnewstars+noldstars ne nstars then stop
     print, 'starnum', nnewstars, noldstars, nstars
     print, 'gasnum', ngas

;;-- fill up ind_ssp star structures
     if i eq 0 then plot=1 else plot=0
     if noldstars gt 0 then SEDM2_BUILDSED, age_ssp, stars, oldstars_minid, sfr, snap_time,plot=plot

 ;;-- read gas particle and new star particle SFHs for this snapshot - 0.5 secs
     if style ne "star_age" then begin
    	 savefile = dir_in+filename_short+'_gassfh'+file_style+'.sav'
	     restore, savefile
    	 gassfh=gassfh[gas_cen_ind,*]
    	 newstarsfh=newstarsfh[newstar_cen_ind, *]
     endif


;;-- loop over SSPs to build integrated spectra
     for j=0,nssps-1 do begin

       ; for star_age method, treat all stars as old stars, so oldstars_minid was set to be -1
       ;therefore, we need only calculate "old stars" for star_age method
       if style ne "star_age" then begin
          ;;-- gas
          if ngas gt 0 then begin
             if age_ssp[j] le 0.01 then spec_g_young = spec_g_young+total(gassfh[*,j])*ssps_lum[*,j] $
             else spec_g_old = spec_g_old+total(gassfh[*,j])*ssps_lum[*,j]
          endif

          ;;-- new stars
          if nnewstars gt 0 then begin
             if age_ssp[j] le 0.01 then spec_ns_young = spec_ns_young+total(newstarsfh[*,j])*ssps_lum[*,j] $
             else spec_ns_old = spec_ns_old+total(newstarsfh[*,j])*ssps_lum[*,j]
          endif
        endif ;;style

        ;;-- old stars
        if noldstars gt 0 then begin
          if style ne 'star_age' then begin
            print, "Meatllicity is only work for star_age method currently! --- 12-Mar-19"
            stop
          endif

          oldstars = stars[ind_oldstars]
          if KEYWORD_SET(with_metal) then begin

            print, "Using metal"
            ostars_metal_bin = uintarr(noldstars)
            sedm2_z_ind, oldstars.metal, Z_models.values, Z_ind
            ostars_metal_bin = Z_ind

            accum_nostars_metal = 0 ;accumulative count of oldstars
            print, "Z_keys   number of the oldstars in this bins"
            for mm=0, N_z-1 do begin
              ind_metal = where(ostars_metal_bin eq mm, mcount)
              if mcount gt 0 then begin
                  print, Z_keys[mm], mcount
                  accum_nostars_metal += mcount

                  mm_oldstars = oldstars[ind_metal]
                  for j=0,nssps-1 do begin
                    ind = where(mm_oldstars.ind_ssp eq j,nn2)
                    ; print, nn2
                    if ind[0] ne -1 then begin
                       if age_ssp[j] le 0.01 then spec_os_young = spec_os_young+total(mm_oldstars[ind].mass)*ssps_lum_metal[*,j,mm] $
                       else spec_os_old = spec_os_old+total(mm_oldstars[ind].mass)*ssps_lum_metal[*,j,mm]
                    endif ; have particle in the ssp
                  endfor

              endif ; have particles in the meatllicity bin

            endfor ; for loop of metallicity
            ;check that all oldstars have been counted
            if accum_nostars_metal ne noldstars then  begin
              print, "Not all oldstars are counted for the spectra, stop here"
              stop
            endif else begin
              print, "noldstars          total star number in all metallicity bins"
              print, noldstars, accum_nostars_metal
              print, "All oldstars have been counted."
            endelse

          endif else begin ; loop over meatllicity bins

            print, "Uni_metal"
            for j=0,nssps-1 do begin

              ; if noldstars gt 0 then begin
              ind = where(oldstars.ind_ssp eq j,nn2)
              if ind[0] ne -1 then begin
                if age_ssp[j] le 0.01 then spec_os_young = spec_os_young+total(oldstars[ind].mass)*ssps_lum[*,j] $
                else spec_os_old = spec_os_old+total(oldstars[ind].mass)*ssps_lum[*,j]
              endif
              ; endif
            endfor
          endelse ; use meatllicity



        ;    ind = where(stars[ind_oldstars].ind_ssp eq j,nn2)
        ;    if ind[0] ne -1 then begin
        ;       if age_ssp[j] le 0.01 then spec_os_young = spec_os_young+total(stars[ind_oldstars[ind]].mass)*ssps_lum[*,j] $
        ;       else spec_os_old = spec_os_old+total(stars[ind_oldstars[ind]].mass)*ssps_lum[*,j]
        ;    endif
        endif


     endfor


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
     struc = {wave:lambda, spec_tau:spec_tau,spec_notau:spec_notau,spec_stars:spec_ns, spec_ns_dust:spec_ns_dust, spec_bulge:spec_os, spec_os_dust:spec_os_dust, spec_gas:spec_g, spec_g_dust:spec_g_dust}
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

     plot, lambda,spec_tau[*,0],/xlog,/ylog,/xs,xtitle='Wavelength [A]',ytitle=textoidl('Luminosity [L_\odot/A]'),title='Orien 0: Snapshot '+str_snap+' with dust'+' (cell region)';+' N*='+string(nstars,form='(I0)')
     oplot, lambda,spec_os_dust[*,0],color=cgcolor('red')
     oplot, lambda,spec_ns_dust[*,0],color=cgcolor('cyan')
     oplot, lambda,spec_g_dust[*,0],color=cgcolor('purple')

     plot, lambda, -alog(spec_tau[*,0]/spec_notau[*,0]),/xlog,/xs,ytitle=textoidl('\tau_\lambda'),xtitle='Wavelength [A]'
     oplot, lambda, -alog(spec_g_dust[*,0]/spec_g[*,0]),color=cgcolor('purple')
     oplot, lambda, -alog(spec_ns_dust[*,0]/spec_ns[*,0]),color=cgcolor('cyan')
     oplot, lambda, -alog(spec_os_dust[*,0]/spec_os[*,0]),color=cgcolor('red')

     ;; shorteer wavelength range
     plot, lambda,spec_notau[*,0],/xs,xr=[3000,9000],title='Orien 0: Snapshot '+str_snap+' N*='+string(nstars,form='(I0)')+' (cell region)'
     oplot, lambda,spec_os[*,0],color=cgcolor('red')
     oplot, lambda,spec_ns[*,0],color=cgcolor('cyan')
     oplot, lambda,spec_g[*,0],color=cgcolor('purple')

     plot, lambda,spec_tau[*,0],/xs,xr=[3000,9000],title='Orien 0: Snapshot '+str_snap+' with dust'+' N*='+string(nstars,form='(I0)')+' (cell region)'
     oplot, lambda,spec_os_dust[*,0],color=cgcolor('red')
     oplot, lambda,spec_ns_dust[*,0],color=cgcolor('cyan')
     oplot, lambda,spec_g_dust[*,0],color=cgcolor('purple')
     print, 'gas_max',max(spec_g_dust[*,0])
     print, 'newstar_max',max(spec_ns_dust[*,0])
     ;print, spec_ns_dust[*,0]
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
