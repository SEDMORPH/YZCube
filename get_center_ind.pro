PRO get_center_ind, particle, cen_part_ind,mass_weight, $
  center=center, cell_size=cell_size, fib_radius=fib_radius, $
  cir_fib=cir_fib, with_PSF=with_PSF
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


  ;check the parameters, stop the program if the right size parameter is not set
  temp=check_size(cell_size=cell_size, fib_radius=fib_radius, cir_fib=cir_fib)


  ;;comment out print, center
  xc=center[0]
  yc=center[1]
  N_part = n_elements(particle)
  xlist = particle.x
  ylist = particle.y

  print, "cell_size:  ", cell_size

  xl = where(xlist le (xc + cell_size/2.0))
  ;;comment out print, xl
  xl_ind = intarr(N_part)
  ;xl_ind = 0
  xl_ind[xl] =1
  ;comment out print, "have xl?",max(xl_ind)

  xg = where(xlist ge (xc - cell_size/2.0))
  xg_ind = intarr(N_part)
  xg_ind[xg]=1
  ;comment out print, "have xg?",max(xg_ind)

  yl = where(ylist le (yc + cell_size/2.0))
  yl_ind = intarr(N_part)
  yl_ind[yl]=1
  ;comment out print, "have yl?",max(yl_ind)

  yg = where(ylist ge (yc - cell_size/2.0))
  yg_ind = intarr(N_part)
  yg_ind[yg]=1
  ;comment out print, "have yg?",max(yg_ind)


  temp = intarr(N_part)
  ; print, n_elements(temp)
  ; print, total(temp)
  temp = temp+1
  print,"particle number at the beginning        " ,total(temp)
  temp = temp * xl_ind
  print,"particle number after requiring x < xmax" ,total(temp)
  temp = temp * xg_ind
  print,"particle number after requiring x > xmin" ,total(temp)
  temp = temp * yl_ind
  print,"particle number after requiring y < ymax" ,total(temp)
  temp = temp * yg_ind
  ;temp = xl_ind * xg_ind * yl_ind * yg_ind
  print,"particle number after requiring y > ymin" ,total(temp)

  if KEYWORD_SET(cir_fib) then begin
    ; print, "Searching particles inside the circle"
    rlist = fltarr(N_part)
    rlist += cell_size + 1.0 ;Initialize to  cell_size+1, (i.e. something greater than cell_size)
    in_box_idx = where(temp eq 1)
    rlist[in_box_idx]=sqrt( (xlist[in_box_idx] - xc)^2 + (ylist[in_box_idx] -yc)^2  )
    ; print, rlist[in_box_idx]

    in_circle = where( rlist le (cell_size/2.0) )
    in_circle_idx = intarr(N_part)
    in_circle_idx[in_circle] = 1
    temp = temp * in_circle_idx
    print, "particle number after requiring distance < r" ,total(temp)

  endif

  cen_part_ind = where(temp eq 1)

  if KEYWORD_SET(with_PSF) then begin
    if NOT KEYWORD_SET(cir_fib) then begin
      print, "PSF function only works for circular fiber currently! --- 11-July-19"
      stop
    endif
    mass_weight = fltarr(N_part)
    ; then do the calculation
  endif else begin
    mass_weight = 1.0
  endelse

END
