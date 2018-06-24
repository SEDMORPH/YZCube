pro get_center_ind, particle, cen_part_ind, center=center, box_size=box_size
  ;;comment out print, center
  xc=center[0]
  yc=center[1]
  N_part = n_elements(particle)
  xlist = particle.x
  ylist = particle.y

  print, box_size

  xl = where(xlist le (xc + box_size/2.0))
  ;;comment out print, xl
  xl_ind = intarr(N_part)
  ;xl_ind = 0
  xl_ind[xl] =1
  ;comment out print, "have xl?",max(xl_ind)

  xg = where(xlist ge (xc - box_size/2.0))
  xg_ind = intarr(N_part)
  xg_ind[xg]=1
  ;comment out print, "have xg?",max(xg_ind)

  yl = where(ylist le (yc + box_size/2.0))
  yl_ind = intarr(N_part)
  yl_ind[yl]=1
  ;comment out print, "have yl?",max(yl_ind)

  yg = where(ylist ge (yc - box_size/2.0))
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

  cen_part_ind = where(temp eq 1)

END
