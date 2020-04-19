FUNCTION check_size,cell_size=cell_size, fib_radius=fib_radius, cir_fib=cir_fib
;+
; cell_size    : the size of the cell, side length of a square fiber or diameter of a circular fiber
; cir_fib      : use circular fiber
; fib_radius   : the radius of the circular fiber, if unset, use (cell_size/2.0)
;                if both fib_radius and cell_size is set, ignore cell_size when using cir_fib
;                ignore fib_radius when using sqaure fiber
;-

  if n_elements(cell_size) eq 0 then cs_set = 0 else cs_set = 1
  if n_elements(fib_radius) eq 0 then fr_set = 0 else fr_set = 1
  ; print, cell_size

  if (cs_set eq 1) && (fr_set eq 1) then begin
    print, "Both cell_size and fib_radius are set."
    if KEYWORD_SET(cir_fib) then begin
      print, "Please specify fib_radius only as we are now using a circular fiber."
    endif else begin
      print, "Please specify cell_size only as we are now using a sqaure fiber."
    endelse
    stop
  endif ; both params are set


  if (cs_set eq 1) AND (fr_set eq 0) then begin
    if KEYWORD_SET(cir_fib) then begin
      print, "We are now using a sqaure fiber, please specify the side length of the fiber by cell_size!"
      stop
    endif else begin
      print, "cell_size:", cell_size
    endelse
  endif

  if (cs_set eq 0) AND (fr_set eq 1) then begin
    if KEYWORD_SET(cir_fib) then begin
      print, "fib_radius:", fib_radius
    endif else begin
      print, "We are now using a sqaure fiber, please specify the side length of the fiber by cell_size!"
      stop
    endelse
  endif ; both params are set

  if (cs_set eq 0) AND (fr_set eq 0) then begin
    if KEYWORD_SET(cir_fib) then begin
      print,"Please specify the radius of the circular fiber by fib_radius."
    endif else begin
      print,"Please specify the side length of the square fiber by cell_size."
    endelse
    stop
  endif

END


FUNCTION params_in_kpc, cell_x_offset, cell_y_offset, cell_size=cell_size, fib_radius=fib_radius,$
  arcsec=arcsec, cir_fib=cir_fib, redshift=redshift
;+
; check whether the parameters are in kpc unit. If not, convert to kpc
;-
  if KEYWORD_SET(arcsec) then begin
      ;; convert the parameters value from arcsec to kpc
      print, "Check the cell parameters in aresec"
      if KEYWORD_SET(cir_fib) then begin
        print, "x_offset | y_offset | fiber radius"
        print, cell_x_offset, cell_y_offset, fib_radius
      endif else begin
        print, "x_offset | y_offset | cell_size"
        print, cell_x_offset, cell_y_offset,cell_size
      endelse
      ; converting
      kpc_in_arcsec = ANGDIAM_FIB( redshift, angsize=1.0)
      cell_x_offset = cell_x_offset * kpc_in_arcsec
      cell_y_offset = cell_y_offset * kpc_in_arcsec
      if KEYWORD_SET(cir_fib) then begin
        fib_radius = fib_radius * kpc_in_arcsec
      endif else begin
        cell_size = cell_size * kpc_in_arcsec
      endelse

      print, "NOTE! The spectra is still in rest-frame."
      print, "NOTE! The spectra is still in rest-frame."
      print, "NOTE! The spectra is still in rest-frame."
      ; cell_str = cell_str+'_in_arcsec_z'+string(redshift,form='(F0.3)')
  endif

  print, "Check the cell parameters in kpc"
  if KEYWORD_SET(cir_fib) then begin
    print, "x_offset | y_offset | fiber radius"
    print, cell_x_offset, cell_y_offset, fib_radius
    return, [cell_x_offset, cell_y_offset, fib_radius]
  endif else begin
    print, "x_offset | y_offset | cell_size"
    print, cell_x_offset, cell_y_offset,cell_size
    return, [cell_x_offset, cell_y_offset, cell_size]
  endelse

END



PRO SEDM2_CELL_SPEC_STAR_AGE, dir_in, dir_out, tauv,mu_d,redshift, cell_x_offset, cell_y_offset,$
    		cell_size=cell_size, cir_fib=cir_fib, fib_radius=fib_radius, arcsec=arcsec, $
                snap=snap_in, style=style, model_str=model_str,$
                models_dir=dir_models, rtfaceon=rtfaceon, one_comp_dust=one_comp_dust, $
                with_metal=with_metal, with_PSF=with_PSF, dir_PSF_weight=dir_PSF_weight, $
                plot_cell_spec=plot_cell_spec
;+
; create spectra for the cell descibed below:
;
; cell_x_offset: the offset on x-axis from the center of the first galaxy
; cell_y_offset: the offset on y-axis from the center of the first galaxy
; cell_size    : the size of the cell, side length of a square fiber or diameter of a circular fiber
; cir_fib      : use circular fiber
; fib_radius   : the radius of the circular fiber, if unset, use (cell_size/2.0)
;                if both fib_radius and cell_size is set, ignore cell_size when using cir_fib
;                ignore fib_radius when using sqaure fiber
; arcsec       : if this keyword is set, the parameters above are in arcsec.
;                Otherwise they are in kpc. Use kpc by default
;                We do not suggest use arcsec!!! As we want keep the spectra in rest-frame
; rtfaceon     : rotate to faceon, if switched on, rotate the disk to be a face on one.
; style        : the spectra style, choose from the following choices
;                  "" (empty string) --> SEDmorph method--> turn out to be useless, removed 
;                  "star_age" --> star_age method, --> the methods descriped in Y.Zheng+2020
;                   also support _eagle and _eagle_minus, but these are not well tested yet.(20-Aug-2018)
;                    _eagle and _eagle_minus --> turn out to be useless, removed
; one_comp_dust: use tau_old for all stars, i.e. tau_young = tau_old
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


  ;check the parameters, stop the program if the right size parameter is not set
  temp=check_size(cell_size=cell_size, fib_radius=fib_radius, cir_fib=cir_fib)

  if (strlowcase(style) ne "star_age")  then begin
      print, "During the devolpment of the code, we also include spectrum styles like" 
      print, "sedmorph, eagle, eagle_minus. They turn out to be unnecessary so removed." 
      print, "This public version only include the star_age method." 
      print, "The star_age method is the one described in Y.Zheng+2020." 
      stop
  endif else file_style='_'+style
;;-- set the output file name
  outstr = '_tauv'+string(tauv,form='(F0.1)')
  outstr = outstr+'_mu'+string(mu_d,form='(F0.1)')
  if KEYWORD_SET(rtfaceon) then outstr=outstr+'_fo'
  if KEYWORD_SET(one_comp_dust) then outstr=outstr+'_one_comp_dust'
;;-- set up plotting file
  cell_str = 'cell_'+string(cell_x_offset, form='(F+0.2)')+string(cell_y_offset, form='(F+0.2)')

  xys = params_in_kpc(cell_x_offset, cell_y_offset, cell_size=cell_size, fib_radius=fib_radius, arcsec=arcsec, cir_fib=cir_fib, redshift=redshift)
  cell_x_offset = xys[0]
  cell_y_offset = xys[1]
  if KEYWORD_SET(cir_fib) then begin
    fib_radius=xys[2]
    cell_str = cell_str+'_cir_radius_'+string(fib_radius, form='(F0.2)' )
  endif else begin
    cell_size=xys[2]
    cell_str = cell_str+'_size_'+string(cell_size, form='(F0.2)' )
  endelse


  if KEYWORD_SET(with_metal) then begin
    if style ne 'star_age' then begin
      print, "Meatllicity is only work for star_age method currently! --- 12-Mar-19"
      stop
    endif
  endif



  if n_elements(snap_in) le 0 then psfile = dir_out+cell_str+'_spectra'+outstr+string(file_style)+'.ps'

;;-- find snapshot filenames from gadget simulation files
  filename = file_search(dir_in+'/*.hdf5',count=nsnap)
  filename0 = filename[0]

  if n_elements(snap_in) gt 0 then begin

     SEDM2_READSNAP, filename[0], stars=stars, /getstars, /no_sfr_log ;need to get minID of old stars from first file
     oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.
     ; for star_age method, treat all stars as old stars, so oldstars_minid should be set to -1
     if style eq 'star_age' then oldstars_minid=-1

     tmp = where(strmatch(filename, '*'+string(snap_in,form='(I3.3)')+'.hdf5') eq 1)
     if tmp[0] eq -1 then message, 'Something wrong with identifying snapshot file: ', snap_in
     filename = filename[tmp]

  endif


  nsnap = n_elements(filename)

  print, 'SEDM2_CELL_SPEC_STAR_AGE: number of snapshots to process',nsnap

;;------------------------------------------------------------------
;; Pre-computations
;;------------------------------------------------------------------
;;-- read SSPs 
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
    ; 0 has been read above
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

  if KEYWORD_SET(plot_cell_spec) then begin
    if n_elements(snap_in) le 0 then ps1c, psfile
  endif
  time = systime(1)

  mass_young = dblarr(nsnap)

  for i=0, nsnap-1 do begin ;Nsnap-1 do begin

     tmp = (strsplit(filename[i],'/',/extract,count=n))[n-1]
     filename_short = (strsplit(tmp,'.',/extract))[0]
     str_snap = (strsplit(filename_short,'_',/extract))[1] ;don't use i as could be only doing a single snapshot

     print, 'SEDM2_CELL_SPEC_STAR_AGE building spectrum for snapshot:'+str_snap
 ;;-- make a directory for each data cube.
     if KEYWORD_SET(cir_fib) then begin
       data_cube_dir = "DataCube"+outstr+'_'+str_snap+string(file_style)+'_cir_radius_'+string(fib_radius, form='(F0.2)' )
     endif else begin
       data_cube_dir = "DataCube"+outstr+'_'+str_snap+string(file_style)+'_size_'+string(cell_size, form='(F0.2)' )
     endelse

     if KEYWORD_SET(with_metal) then data_cube_dir=data_cube_dir+'_with_metal'
     if KEYWORD_SET(with_PSF) then data_cube_dir=data_cube_dir+'_with_PSF'
     data_cube_dir = data_cube_dir + '/'
     if file_test(dir_out+data_cube_dir) eq 0 then spawn, 'mkdir '+ dir_out+data_cube_dir
     if KEYWORD_SET(plot_cell_spec) then begin
       if n_elements(snap_in) gt 0 then psfile = dir_out+data_cube_dir+cell_str+'.ps' & ps1c, psfile
     endif
;;-- outfile (single snapshot, all components, all orientations)
     outfile_fits = dir_out+ data_cube_dir+'spec_'+cell_str+'.fits'
     print, outfile_fits
     if KEYWORD_SET(plot_cell_spec) then print, psfile
     ; wait, 10
;;-- output arrays for each component of the final image
;;   NOTE: for star_age method, all stars are treated as old stars
;;   In other spectrum styles, there are 3 components, gas, new stars formed during the simulation
;;   and old stars that are already in place prior to start of simulation
     spec_g_old = (spec_g_young = (spec_ns_old = (spec_ns_young = (spec_os_old = (spec_os_young = fltarr(nlambda))))))

;;-- read simulation files  - 9 secs
     ; use no_sfr_log, unless you want to plot the SFR log file
     SEDM2_READSNAP, filename[i], stars=stars, snap_time=snap_time,/getstars,/no_sfr_log


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

      print, "Check stars:"
      get_center_ind, stars, cen_part_ind, mass_weight, center=center, cell_size=cell_size,fib_radius=fib_radius, cir_fib=cir_fib, with_PSF=with_PSF, redshift=redshift, dir_PSF_weight=dir_PSF_weight
      ;print, cen_part_ind
      stars.mass = stars.mass * mass_weight
      stars = stars[cen_part_ind]



     ind_oldstars = where(stars.id ge oldstars_minID,compl=ind_newstars, noldstars) ;ID the old stars
     if ind_newstars[0] ne -1 then nnewstars = n_elements(ind_newstars) else nnewstars = 0

;;-- fill up ind_ssp star structures
     if i eq 0 then plot=1 else plot=0
     if not KEYWORD_SET(plot_cell_spec) then plot=0
     if noldstars gt 0 then SEDM2_BUILDSED, age_ssp, stars, oldstars_minid, snap_time,plot=plot


    ;;-- old stars
    if noldstars gt 0 then begin

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
              for j=0,nssps-1 do begin ; loop over ssp
                ind = where(mm_oldstars.ind_ssp eq j,nn2)
                ; print, nn2
                if ind[0] ne -1 then begin
                   if age_ssp[j] le 0.01 then spec_os_young = spec_os_young+total(mm_oldstars[ind].mass)*ssps_lum_metal[*,j,mm] $
                   else spec_os_old = spec_os_old+total(mm_oldstars[ind].mass)*ssps_lum_metal[*,j,mm]
                endif ; have particle in the ssp
              endfor ; loop over ssp
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

    endif ; if there are old stars


;;-- sum components and add dust
;;-- here we treat all stars as old stars
     ; spec_ns = spec_ns_young+spec_ns_old
     ; spec_ns_dust = spec_ns_young*exp(-tau_young)+spec_ns_old*exp(-tau_old)
     ; spec_g = spec_g_young+spec_g_old
     ; spec_g_dust = spec_g_young*exp(-tau_young)+spec_g_old*exp(-tau_old)
     spec_os = spec_os_young+spec_os_old
     spec_os_dust = spec_os_young*exp(-tau_young)+spec_os_old*exp(-tau_old)

     spec_notau = spec_os
     spec_tau = spec_os_dust

;;-- PCA components


;;-- save output file
     ; struc = {wave:lambda, spec_tau:spec_tau,spec_notau:spec_notau,spec_stars:spec_ns, spec_ns_dust:spec_ns_dust, spec_bulge:spec_os, spec_os_dust:spec_os_dust, spec_gas:spec_g, spec_g_dust:spec_g_dust}
     struc = {wave:lambda, spec_tau:spec_tau,spec_notau:spec_notau}
     mwrfits, struc, outfile_fits,/create



     if KEYWORD_SET(plot_cell_spec) then begin
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
     endif


  endfor


END
