;;** AIM: return either ages, or full ssp arrays, for a given model
;;   and metallicity
;;
;; INPUT: model = 'bc03' or 'cb08' at the moment
;;        Z = metallicity as a string 12-82 for cb08, 22-72 for bc03 (no interpolation)
;;            or as a value (interpolate - not yet implemented)
;;
;; OUTPUT: nlambda = nb wavelength bins
;;         nssps = nb ssps
;;         age_ssp = age array
;;
;; KEYWORDS: /novel = do not include third dimension to seds array
;;                    containing different velocities
;;     
;; RETURN: default is {seds[nlambda,nssps,{nvel}], age, lambda} 
;;         age is in Gyr
;; 
;; TO DO: interpolation over metallicity not yet implemented
;; 
;;******************************************************************
;;******************************************************************

FUNCTION SEDM2_GETSSPS, models_dir, model, Z, nlambda, nssps, age_ssp,  vel=vel


; if ((model eq 'CB08') or (model eq 'cb08')) then $
;       stem = !dataDIR+'bc03/models/CB08/cb2008_hr_miless_n' $
; else if ((model eq 'BC03') or (model eq 'bc03')) then $
        stem = models_dir+'bc2003_hr_m' ;$ 
; else message, 'model not recognised'

;;-- check string i.e. single metallicity
  if size(Z,/type) eq 7 then filename = stem+Z+'_chab_ssp.ised' $
  else message, 'interpolation between metallicities not yet implemented'

;;-- read in SSPs
  seds_zerovel = SEDM2_READISED(filename,t,lambda)   
  
  t = t/1e9                     ; in Gyr to match simns
  nlambda = n_elements(lambda)
  nssps = n_elements(t)
  age_ssp = t

;;-- if not needing velocity information
  if NOT(keyword_set(vel)) then ssp = {seds:seds_zerovel, age:t, lambda:lambda} $
  
  else begin
     print, 'I am wrong'
     stop

;;-- include a range of velocities for the particles
;;-- Information from velinfo.inc
nvel = 17   ; always odd number so zero is included 
velrange = 400  ; +/- this value is used
vel = findgen(nvel)        ;0 -> 6.0
vel = vel / (0.5*(nvel-1)) ;0 -> 2.0
vel = (vel -1.0)*velrange    ;-400 -> 400
dvel = vel[1]-vel[0]
      
;;-- currently only works where delta_lambda = 1AA, so remove the rest
;;   of array to ensure its not used
;;   note extra inelegant cut at 1000AA is because of 2 pesky bins at ~131AA
      ;; dl = lambda[1:nlambda-1] - lambda[0:nlambda-2]
      ;; ind = (where(dl gt 0.99 and dl lt 1.01 and lambda gt 1000,nn_tmp))[1:nn_tmp-2]
      ;; lambda = lambda[ind]
      ;; seds_zerovel = seds_zerovel[ind,*]
      ;; nlambda = n_elements(ind)
      ;; splog, 'cut down wavelength range because of velocities:', minmax(lambda)
      
      seds = fltarr(nlambda,nssps,nvel)
      for i=0,nvel-1 do begin
         ;; lambda in AA ; c in m/s ; vel in km/s -> dl in AA ;
         ;; 1AA/pixel -> dl in pixels as required for this interpolation method
         ;; delta_lambda = lambda * vel[i] *1e3 / !speedoflight 

         ;; ;; set sign here so that a blueshift for -tve velocity
         ;; arr_interp = findgen(nlambda)-delta_lambda                                           ;position of interpolates
         ;; seds_new = transpose(interpolate(transpose(seds_zerovel), arr_interp,missing=0.0))   ;; checked
         
         wave_new = lambda*(1+vel[i] *1e3 / !speedoflight) ;+tve vel is redshifted
         wave_shift = lambda-wave_new                      ;redshift = wave_shift is negative
         
;!!! THIS NEEDS MORE THINKING ABOUT
         
         seds[*,*,i] = seds_new ;now put shifted array into observed array
         
      endfor
      
      ssp = {seds:seds, age:t, lambda:lambda, vel:vel}
   endelse
  
  ;; check the ISED files, sometimes there are problems
  if size(t,/dim) ne (size(seds_zerovel,/dim))[1] then stop
  if size(lambda,/dim) ne (size(seds_zerovel,/dim))[0] then stop
  
  return, ssp

END
