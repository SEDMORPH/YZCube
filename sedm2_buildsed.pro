;;** for the old stars that are already in place prior to start of simulation,
;;assign SSPs based on their formation time as set up by Natalia. This
;;involves filling the ind_ssp array in the structure "stars"

PRO SEDM2_BUILDSED, age_ssp, stars, oldstars_minID, sfr, snap_time, plot=plot

  Nssp = n_elements(age_ssp)

  ;;-- find the old stars in the array
  ind_oldstars = where(stars.id ge oldstars_minID, noldstars)

  if noldstars ne 0 then oldstarsfh = fltarr(noldstars,Nssp) else return

  ;;-- assign closest SSP to them. This might cause discreteness
  ;;   problems early on in simulation where the SSPs are much more
  ;;   finely spaced than the age of the stars. But it probably
  ;;   isn't worth doing top hats, as we'd have to then include them
  ;;   in the gassfh calculations.
  diff = abs(rebin(stars[ind_oldstars].age, Noldstars, Nssp,/sample) - rebin(transpose(age_ssp),Noldstars,Nssp,/sample) )
  minarr = min(diff,dimension=2,ind_ssp_1d)
  ind_ssp_2d = array_indices(size(diff,/dim),ind_ssp_1d,/dimensions) ;convert back to 2D indices
  stars[ind_oldstars].ind_ssp = reform(ind_ssp_2d[1,*])

  ;;-- plot the SFH
  if keyword_set(plot) then begin
     !p.multi=[0,1,2]
     hist = histogram(snap_time - stars[ind_oldstars].age,binsize=0.1,loc=x) ;formation time
     plot, sfr[0,*],sfr[1,*],xr=[0,13.7]
     oplot, x, hist*stars[0].mass/0.1d9,color=cgcolor('blue')
     plot, sfr[0,*],sfr[1,*],xr=[0,13.7],/ylog,yr=[0.01,max(sfr[1,*])]
     oplot, x, hist*stars[0].mass/0.1d9,color=cgcolor('blue')
  endif

END
