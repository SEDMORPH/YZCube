;;-- convolve SED with a filter function
;;
;; Using equation given here http://www.gama-survey.org/wiki/extra:filtercurves
;;
;; INPUT lambda (xx) and f_lambda (flux) in rest-frame : this code redshifts both
;; OUTPUT : f_lambda_obs 

;; exact copy of vwsc_conv_filter August 2015

FUNCTION sedm2_conv_filter,xx,flux,lam_fil,fil,z
  
  
  if (size(flux))[0] eq 1 then ngal = 1 else ngal = (size(flux,/dim))[1]

  if n_elements(z) eq 1 then begin
     x = xx*(1+z)                                           ;redshfit the SSP
     
     if (max(x) lt max(lam_fil)) or (min(x) gt min(lam_fil)) then begin
        print, 'VWSC_CONV_FILTER: warning - filter extends beyond wavelength array'
        out = 0.0 
     endif else begin
        ind = where(x gt min(lam_fil) and x lt max(lam_fil),npix) ;cut down to range in filter
        x = x[ind]
        ff = flux[ind,*]/(1+z)  ;(*wave_rest/wave_obs) : rest -> observed frame /AA
        linterp, lam_fil,fil,x,filr_interp
        
        dlam = x[1:npix-1]-x[0:npix-2]   
        dlam = [dlam,dlam[npix-2]]
        
        if ngal eq 1 then out = total(x*ff*dlam*filr_interp)/total(x*dlam*filr_interp)
        ;; SLOW: total(ff*rebin(dlam*filr_interp,npix,ngal),1)/total(dlam*filr_interp)
        if ngal gt 1 then out =  reform((ff ## transpose(x*dlam*filr_interp)) /total(x*dlam*filr_interp) )
     endelse

  endif else begin              ;multiple redshifts

     out = fltarr(n_elements(z),ngal) 

     for i=0,n_elements(z)-1 do begin
        x = xx*(1+z[i])                                           ;redshfit the SSP

        if (max(x) lt max(lam_fil)) or (min(x) gt min(lam_fil)) then out[i] = 0.0 else begin
           
           ind = where(x gt min(lam_fil) and x lt max(lam_fil),npix) ;cut down to range in filter
           x = x[ind]
           ff = flux[ind,*]/(1+z[i]) ;(*wave_rest/wave_obs) :  rest -> observed frame /AA

           linterp, lam_fil,fil,x,filr_interp
           
           dlam = x[1:npix-1]-x[0:npix-2]   
           dlam = [dlam,dlam[npix-2]]

           if ngal eq 1 then out[i] = total(x*ff*dlam*filr_interp)/total(x*dlam*filr_interp)
           ;; SLOW: transpose(total(ff*rebin(dlam*filr_interp,npix,ngal),1)/total(dlam*filr_interp))
           if ngal gt 1 then out[i,*] = transpose( reform((ff ## transpose(x*dlam*filr_interp)) /total(x*dlam*filr_interp) ))
        endelse

     endfor
  endelse

  
  return,  out
   
END



