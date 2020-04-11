; +
;
; This procedure read in the 12 values Metallicity given by sedm2_readsnap,
; then compute Z for particles.
; We will find out the closest values in the Z_values for each particles
; and finally return the indices in Z_values list.
;
;
; MODIFICATION HISTORY:
; 	Written by:	Yirui Zheng, March 2018
; -


PRO sedm2_z_ind, Metallicity, Z_values, Z_ind
  
  Z = (total(Metallicity, 1) - Metallicity[6, *] - Metallicity[0, *] ) / total(Metallicity, 1)
  ; in some cases, the metallicity is extremely low,
  ; resulting Z values calculated by the method above to be a negative value
  ; in these cases, assgin the particle to the lowest Z bin. 
  nega_ind = where( Z le 0 )
  Z[nega_ind] = min(Z_values)

  n_Z = n_elements(Z_values)
  Npart = n_elements(Z)
  Z_ind = intarr(Npart)
  ;;-- find closest Z values to particles
  diff = abs( rebin(alog10(Z), Npart, n_Z, /sample) - rebin(transpose(alog10(Z_values)), Npart, n_Z, /sample) )
  ;--------- get the min of each row of the diff       ----------
  ;----------may cause Segmentation fault (core dumped)----------
  ; minarr = min(diff,dimension=2,ind_ssp_1d)
  ; ind_ssp_2d = array_indices(size(diff,/dim),ind_ssp_1d,/dimensions) ;convert back to 2D indices
  ; Z_ind = reform(ind_ssp_2d[1,*])
  ;--------- get the min of each row of the diff       ----------
  ;----------may cause Segmentation fault (core dumped)----------

  for i=0, Npart-1 do begin
      temp_arr = reform(diff[i,*])
      Z_ind[i] = where(temp_arr EQ min(temp_arr))
  endfor


END
