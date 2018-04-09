;
;+
; NAME:
;       vw_integral.pro
;
; PURPOSE:
;       Routine to perform trapezoidal integration in X,Y between limits
;       xmin to xmax.
;
; CALLING SEQUENCE:
;	result = integral(x,y,xmin,xmax)
;
; INPUTS:
;	x,y - vectors to be integrated
;	xmin,xmax - vectors with lower and upper integral limits
; OUTPUTS:
;	the integrations between xmin and xmax are returned
;	as the function result
; RESTRICTIONS:
;	The values in XMIN must be greater than or equal to the minimum
;	of X. The values in XMAX must be less than or equal to the 
;	maximum of X. X must be in ascending order.
;
; HISTORY:
;	Version 1,  18/07/06
;-
;------------------------------------------------------------------

FUNCTION VW_INTEGRAL, X, F, wl_1, wl_2

;; Nb of spectra
if (SIZE(F))[0] eq 1 then ngal = 1 $
else ngal = (SIZE(F,/dim))[1]

npix = n_elements(X)

;; make sure continuum windows are in wave range
if min(X) gt wl_1 or max(X) lt wl_2 then return, lonarr(ngal)-1L

;; ensure that X  runs from low to high
X2 = double(X[SORT(X)])
F2 = double(F[SORT(X),*])

i_int = WHERE(X2 ge wl_1 and X2 le wl_2,nr)

;; whole bins
f_int = total(rebin(X2[i_int[1:nr-1]]-X2[i_int[0:nr-2]],nr-1,ngal)*(F2[i_int[1:nr-1],*]+F2[i_int[0:nr-2],*])/2.,1)

;; low wavelength part bins
if i_int[0] ne 0 then f_int = f_int + ((X2[i_int[0]]-wl_1)/(X2[i_int[0]]-X2[i_int[0]-1]))*(F2[i_int[0],*]+F2[i_int[0]-1,*])/2.

;; high wavelength part bins
if i_int[nr-1] ne npix -1 then f_int = f_int + ((wl_2-X2[i_int[nr-1]])/(X2[i_int[nr-1]+1]-X2[i_int[nr-1]]))*(F2[i_int[nr-1]+1,*]+F2[i_int[nr-1],*])/2.

return, f_int

END
