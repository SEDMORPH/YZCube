FUNCTION inside_the_box, xlist, ylist, N_part, xc, yc, box_size=box_size

  xl = where(xlist le (xc + box_size/2.0))
  xl_ind = intarr(N_part)
  xl_ind[xl] =1
  ;comment out print, "have xl?",max(xl_ind)

  xg = where(xlist ge (xc - box_size/2.0))
  xg_ind = intarr(N_part)
  xg_ind[xg]=1

  yl = where(ylist le (yc + box_size/2.0))
  yl_ind = intarr(N_part)
  yl_ind[yl]=1

  yg = where(ylist ge (yc - box_size/2.0))
  yg_ind = intarr(N_part)
  yg_ind[yg]=1

  inside_the_box =  xl_ind * xg_ind * yl_ind * yg_ind

  return, inside_the_box
END



PRO get_center_ind, particle, cen_part_ind,mass_weight, $
  center=center, cell_size=cell_size, fib_radius=fib_radius, $
  cir_fib=cir_fib, with_PSF=with_PSF, redshift=redshift, dir_PSF_weight=dir_PSF_weight
;+
; To get the indices that tells us wether the particle should be considered.
; If the PSF is not required, then simply select the particles inside the cell or the circle;
; set mass_weight to 1. as it is unnecessary.
; If the PSF effect is required, select the particles inside a radius of 4 kpc
; and calculate the mass_weight as well. The further the particle is from the fiber,
; the less weight the particle should have. If the particle is further than 4 kpc,
; its weight will be smaller then 1e-5, set the weight to 0 and,
; therefore, ignore the contribution from this particle.
;
; cell_size    : in kpc, the size of the cell, side length of a square fiber
; cir_fib      : use circular fiber
; fib_radius   : in kpc,  the radius of the circular fiber
;-

    xc=center[0]
    yc=center[1]
    N_part = n_elements(particle)
    xlist = particle.x
    ylist = particle.y

; square fiber
  if NOT KEYWORD_SET(cir_fib) then begin
    temp = inside_the_box(xlist, ylist, N_part, xc, yc, box_size=cell_size)
    print,"particle number in side the square fiber:" ,total(temp)
    cen_part_ind = where(temp eq 1)

  endif



  if KEYWORD_SET(cir_fib) then begin
    ; print, "Searching particles inside the circle"
    ; print, fib_radius
    if NOT KEYWORD_SET(with_PSF) then begin
      temp = inside_the_box(xlist, ylist, N_part, xc, yc, box_size=(fib_radius*2))
      rlist = fltarr(N_part)
      rlist += fib_radius + 1.0 ;Initialize to  cell_size+1, (i.e. something greater than cell_size)
      in_box_idx = where(temp eq 1)
      rlist[in_box_idx]=sqrt( (xlist[in_box_idx] - xc)^2 + (ylist[in_box_idx] -yc)^2  )
      ; print, rlist[in_box_idx]

      in_circle = where( rlist le fib_radius )
      in_circle_idx = intarr(N_part)
      in_circle_idx[in_circle] = 1
      temp = temp * in_circle_idx
      cen_part_ind = where(temp eq 1)
      print,"particle number in side the circular fiber(NO PSF effect):" ,total(temp)
    endif else begin
      temp = inside_the_box(xlist, ylist, N_part, xc, yc, box_size=(4.0*2))
      rlist = fltarr(N_part)
      rlist += 4.0 + 1.0 ;Initialize to  cell_size+1, (i.e. something greater than cell_size)
      in_box_idx = where(temp eq 1)
      rlist[in_box_idx]=sqrt( (xlist[in_box_idx] - xc)^2 + (ylist[in_box_idx] -yc)^2  )
      ; print, rlist[in_box_idx]

      in_circle = where( rlist le 4.0 )
      in_circle_idx = intarr(N_part)
      in_circle_idx[in_circle] = 1
      temp = temp * in_circle_idx
      print,"particle number in side the circular fiber(with PSF effect):" ,total(temp)
      cen_part_ind = where(temp eq 1)

      ; hyp_gassfh
    endelse

  endif


  if KEYWORD_SET(with_PSF) then begin
    if NOT KEYWORD_SET(cir_fib) then begin
      print, "PSF function only works for circular fiber currently! --- 11-July-19"
      stop
    endif
    if  n_elements(redshift) eq 0 then begin
      print, "Need to specify the redshift to choose the PSF."
      print, "As we are interested in data around 4000A and our PSF is not wavelength dependent, the PSF around the band that the 4000A break falls into is used at all other wavelength. For example, if the galaxy is at z=0.2, then the g-band PSF is used at all wavelength."
      stop
    endif
    print, "Calculate the mass weight with the PSF effect..."
    mass_weight = fltarr(N_part)
    ; then do the calculation
    bin_size = 0.05 ;kpc
    ; effective wavelength of ugriz, http://skyserver.sdss.org/dr1/en/proj/advanced/color/sdssfilters.asp
    wave_eff = [3543.0, 4770.0, 6231.0, 7625.0, 9134.0] 
    band_list = ['u', 'g', 'r', 'i', 'z']
    temp = min( abs(wave_eff - 4000.0*(1+redshift)), band_id )
    band_name = band_list[band_id] 
    print, "redshift: ", redshift, "      band_name: ", band_name
    filename = dir_PSF_weight+band_name+"_band_PSF_mass_weight_res_"+string(bin_size, form='(F0.2)')+'kpc.fits'
    data = mrdfits(filename, 1, hdr)
    arr_weight_bin = data.PSF_weight;read in the weight at different radius
    bin_idx = round(rlist/bin_size)
    mass_weight = arr_weight_bin[bin_idx]
  endif else begin
    mass_weight = 1.0
  endelse

END
