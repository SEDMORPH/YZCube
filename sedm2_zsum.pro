;; From: Flux Units and NICMOS
;C. J. Skinner
;February 23, 1996

;; lambda in microns, f_lam in erg/s/cm^2/AA, f_nu in janksy

FUNCTION fluxtojansky,f_lam,lambda

  ind = where(lambda gt 1000)
  if ind[0] ne -1 then message, 'FLUXTOJANKSY: lambda must be in microns!'

  beta = 3e-13                  ; F_lambda in erg/s/cm^2/AA
  ll = double(lambda)           ;just in case

  f_nu = f_lam*ll^2/beta

return, f_nu

END


PRO SEDM2_ZSUM, nx, ny, xx, yy, gas, stars, ind_newstars, ind_oldstars, gassfh, newstarsfh, rotate, dist, $
                ssps_lum, age_ssp, ll_eff, tau_young, tau_old, $
                flux_all_jansky, flux_all_jansky_dust, mass_young, vel=vel

  @sedm2_codeunits.inc
  age_young = 0.01              ;1e7 years - IDs young stars

  if size(gas,/type) eq 8 then ngas = n_elements(gas) else ngas=0
  if ind_newstars[0] ne -1 then nnewstars = n_elements(ind_newstars) else nnewstars = 0
  if ind_oldstars[0] ne -1 then noldstars = n_elements(ind_oldstars) else noldstars = 0

  ;;------------------------------------------------------------------
  ;;-- locate positions of all particles in the grid

  if ngas gt 0 then xyz_gas = transpose(rotate ## [[gas.x],[gas.y],[gas.z]]) ; transpose because hist_nd works on Ndim xPdata arrays
  if nnewstars gt 0 then xyz_newstars = transpose(rotate ## [[stars[ind_newstars].x],[stars[ind_newstars].y],[stars[ind_newstars].z]])
  if noldstars gt 0 then xyz_oldstars = transpose(rotate ## [[stars[ind_oldstars].x],[stars[ind_oldstars].y],[stars[ind_oldstars].z]])


  if Ngas gt 1 then $
     hist_g = hist_nd(xyz_gas[0:1,*],  nbins=[nx,ny], min=[min(xx),min(yy)], max = [max(xx), max(yy)], reverse_indices=RI_g)

  if Nnewstars gt 1 then $
     hist_ns = hist_nd(xyz_newstars[0:1,*],  nbins=[nx,ny], min=[min(xx),min(yy)], max = [max(xx), max(yy)], reverse_indices=RI_ns)

  if Noldstars gt 1 then $
     hist_os = hist_nd(xyz_oldstars[0:1,*],  nbins=[nx,ny], min=[min(xx),min(yy)], max = [max(xx), max(yy)], reverse_indices=RI_os)


  ;;------------------------------------------------------------------
  ;;-- useful stuff
  nssps = (size(ssps_lum,/dim))[0]
  nlambda = n_elements(ll_Eff)

  if (size(ssps_lum))[0] ge 2 then if (size(ssps_lum,/dim))[1] ne nlambda then stop

  ;; SSPs affected by BCs
  ind_young = where(age_ssp lt age_young,nn_young,compl=ind_old)
  nn_old = n_elements(ind_old)

  ;;------------------------------------------------------------------
  ;; loop over all pixels
  flux_all_jansky = (flux_all_jansky_dust = dblarr(nx,ny,nlambda))
  mass_young = 0d               ;mass of stars behind BCs (i.e. SFR in last 1e7 yrs)

  lum_g_out = (lum_ns_out = (lum_os_out = dblarr(nx,ny,nlambda))) ;for debugging
  sfr_g = dblarr(nx,ny)

  for kbin=0,nx-1 do begin
     for lbin=0,ny-1 do begin
        ind = kbin+long(nx)*lbin ;1D index for this bin
        lum_os = (lum_g = (lum_ns = (lum_young = (lum_old = dblarr(nlambda)))))

        ;;------------------------------------------------------------------
        ;; first check whether there are any particles in this x,y bin

        nexist_ns = (nexist_g = (nexist_os =  0))
        if Nnewstars gt 1 then if ri_ns[ind] ne ri_ns[ind+1] then nexist_ns=1
        if Ngas gt 1 then if ri_g[ind] ne ri_g[ind+1] then nexist_g=1
        if Noldstars gt 1 then if ri_os[ind] ne ri_os[ind+1] then nexist_os=1

        if nexist_ns+nexist_g+nexist_os eq 0 then continue



        ;;------------------------------------------------------------------
        ;;-- no velocity
        if not(keyword_Set(vel)) then begin

           ;; gas

           if nexist_g eq 1 then begin
              ind2 = ri_g[ri_g[ind]:ri_g[ind+1]-1] ;reverse indices relating to this x,y bin

              ;; verifying reverse indices - this is slow!
              ;indtmp = where(gas.x ge float(xx[kbin]) and gas.x lt float(xx[kbin+1]) and gas.y ge float(yy[lbin]) and gas.y lt float(yy[lbin+1]))
              ;if abs(n_elements(ind2) - n_elements(indtmp)) gt 1 then stop

              mass = total(gassfh[ind2,*],1) ;total mass in each ssp, [nssp]

              lum_g = lum_g + total(rebin(mass,nssps,nlambda)*ssps_lum,1) ;total(mass in ssp * lum of ssp)
              sfr_g[kbin,lbin] = total(gas[ind2].sfr)
              lum_young = lum_young + total(rebin(mass[ind_young],nn_young,nlambda)*ssps_lum[ind_young,*],1)
              lum_old = lum_old + total(rebin(mass[ind_old],nn_old,nlambda)*ssps_lum[ind_old,*],1)
              mass_young = mass_young+total(mass[ind_young]) ;for SFR

           endif

           ;; newstars

           if nexist_ns eq 1 then begin
              ind2 = ri_ns[ri_ns[ind]:ri_ns[ind+1]-1] ;reverse indices relating to this x,y bin

              mass = total(newstarsfh[ind2,*],1)      ;total mass in each ssp, [nssp]

              lum_ns = lum_ns + total(rebin(mass,nssps,nlambda)*ssps_lum,1)
              lum_young = lum_young + total(rebin(mass[ind_young],nn_young,nlambda)*ssps_lum[ind_young,*],1)
              lum_old = lum_old + total(rebin(mass[ind_old],nn_old,nlambda)*ssps_lum[ind_old,*],1)
              mass_young = mass_young+total(mass[ind_young]) ;for SFR
           endif

           ;; oldstars

           if nexist_os eq 1 then begin
              ind2 = ri_os[ri_os[ind]:ri_os[ind+1]-1] ;reverse indices relating to this x,y bin

              if n_elements(ind2) gt 1 then begin
                 ind_uniq = uniq(stars[ind_oldstars[ind2]].ind_ssp,sort(stars[ind_oldstars[ind2]].ind_ssp))
                 inds = stars[ind_oldstars[ind2[ind_uniq]]].ind_ssp ;unique ssp-ids of particles in this x,y bin
              endif else inds = stars[ind_oldstars[ind2]].ind_ssp

              for s=0,n_elements(inds)-1 do begin ;for each of these unique ssp-ids, add SSP to lum
                 ind3 = where(stars[ind_oldstars[ind2]].ind_ssp eq inds[s])
                 if age_ssp[inds[s]] le age_young then lum_young = lum_young + total(stars[ind_oldstars[ind2[ind3]]].mass)*ssps_lum[inds[s],*]
                 if age_ssp[inds[s]] gt age_young then lum_old = lum_old + total(stars[ind_oldstars[ind2[ind3]]].mass)*ssps_lum[inds[s],*]

                 lum_os = lum_os+total(stars[ind_oldstars[ind2[ind3]]].mass)*ssps_lum[inds[s],*]
              endfor
           endif

           ;;------------------------------------------------------------------
           ;;-- with velocity information
        endif else begin

           stop                 ;doesn't work below here


        endelse

        lum_all = lum_g+lum_ns+lum_os ;Lsol/AA/pixel
        lum_all_dust = lum_young*exp(-tau_young)+lum_old*exp(-tau_old)

        flux_all = lum_all*Lsol_in_erg/(4*!pi*dist^2) ; models are in Lsol/AA, convert to erg/s/cm^2/AA using redshift
        flux_all_dust = lum_all_dust*Lsol_in_erg/(4*!pi*dist^2)

        flux_all_jansky[kbin,lbin,*] = fluxtojansky(flux_all,ll_eff/1d4) ; convert to jansky/pixel
        flux_all_jansky_dust[kbin,lbin,*] = fluxtojansky(flux_all_dust,ll_eff/1d4)

        ;; for testing
        lum_g_out[kbin,lbin,*]  = lum_g
        lum_ns_out[kbin,lbin,*]  = lum_ns
        lum_os_out[kbin,lbin,*]  = lum_os

     endfor
  endfor

  ;;-- debugging
  ;; this shows that although there is gas throughout the disk, the
  ;; luminosity is totally driven by the SFR. The outer gas is not
  ;; forming stars, so the luminosity and gas particle distributions
  ;; are very different.

  ;; !p.multi=[0,1,2]
  ;; if total(hist_g) gt 0 then begin
  ;;    plot, total(hist_g,2)/total(hist_g)
  ;;    oplot, total(lum_g_out[*,*,2],2)/total(lum_g_out[*,*,2]),color=cgcolor('red')
  ;;    oplot, total(sfr_g,2)/total(sfr_g),color=cgcolor('blue')
  ;; endif
  ;; if total(hist_os) gt 0 then begin
  ;;    plot, total(hist_os,2)/total(hist_os)
  ;;    oplot, total(lum_os_out[*,*,2],2)/total(lum_os_out[*,*,2]),color=cgcolor('red')
  ;; endif


  return


END
