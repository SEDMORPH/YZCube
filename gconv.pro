;;*** PLEASE NOTE THAT SPECTRA MUST BE LOGARITHMICALLY BINNED TO BE
;;    CONVOLVED TO A CONSTANT VELOCITY DISPERSION!!!!

;; to use e.g.  flux = gconv(flux, sigma_pix)
;; where     
;;     cspeed = 2.99792e5
;;     loglamtov = cspeed * alog(10)
;;     dw = alog10(wave_extend[1]) - alog10(wave_extend[0])

;;     vdisp_add = sqrt(vdisp^2 - data_disp^2) ; Deconvolve template resolution
;;     sigma_pix = vdisp_add / loglamtov / dw

;; and
;;    data_disp is velocity dispersion of model spectra. Use
;     75. for Stelib/BC03; 26. for IndoUS/BC03; 58. for MILES/BC03

function gconv, x, sigma, edge_wrap=edge_wrap, fwhm = fwhm, $
                edge_truncate = edge_truncate
   ;; Convolve the x-vector with a Gaussian kernel - the kernel size
   ;; is set to 4 times the sigma of the Gaussian.


   if (n_elements(fwhm) gt 0) then $
    sigma = fwhm/2.3548

   ;; special case for no smoothing
   if sigma eq 0 then return, x

   binfactor=1
   ksize=round(4.0*sigma[0]+1.0)*2
   xx = findgen(ksize)-ksize/2

   kernel=exp(-xx^2/(2*sigma[0]^2))
   kernel=kernel/total(kernel)

   sm = convol(x, kernel, edge_wrap=edge_wrap, edge_truncate=edge_truncate)

   return, sm
end 

