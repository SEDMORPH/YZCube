; +
;
; This procedure read in the 12 values Metallicity given by sedm2sedm2_readsnap,
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

  n_Z = n_elements(Z_values)
  Npart = n_elements(Z)
  ;;-- find closest Z values to particles
  diff = abs( rebin(alog10(Z), Npart, n_Z) - rebin(transpose(alog10(Z_values)), Npart, n_Z) )
  minarr = min(diff,dimension=2,ind_ssp_1d)
  ind_ssp_2d = array_indices(size(diff,/dim),ind_ssp_1d,/dimensions) ;convert back to 2D indices
  Z_ind = reform(ind_ssp_2d[1,*])

  ; help, z_ind

END
