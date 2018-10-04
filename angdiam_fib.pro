;; ANGDIAM_FIB
;;
;; AIM: determine the proper physical size (not comoving) covered by a 3" fibre as a function
;; of redshift.
;;
;; INPUT: z = redshift, flt or fltarr
;; OUTPUT: proper size in kpc
;; OPTIONAL INPUT: angsize, default 3" SDSS fiber
;;                 Omega_m, Lambda0, H0 (default 0.3,0.7,70)
;; For inverse (kpc -> arcsec) see angdiam_kpc.pro
;;******************************************************************

FUNCTION ANGDIAM_FIB, z, angsize=angsize,Omega_m=Omega_m, Lambda0=Lambda0, H0=H0,silent=silent

if N_ELEMENTS(z) eq 0 then begin
    print,'USAGE: size(proper kpc)=ANGDIAM_FIB( z, angsize=angsize,Omega_m=Omega_m, Lambda0=Lambda0, H0=H0)'
    return,0
endif

;arcsec - default size of SDSS fiber
if NOT(KEYWORD_SET(angsize)) then angsize = 3

diam = angsize/(60.*60.*360/(2*!PI)) ;radians

if NOT(KEYWORD_SET(Omega_m)) then Omega_m=0.3
if NOT(KEYWORD_SET(Lambda0)) then Lambda0=0.7
if NOT(KEYWORD_SET(H0)) then H0=70

lumdist = lumdist(z,Omega_m=Omega_m, Lambda0=Lambda0, H0=H0,silent=silent)

;; angsize = (360/2pi *60*60) * dl * (1+z)^2/(1000*lumdist)
;; lumdist is in Mpc
;; dl is proper size in kpc
;; angsize is in arcsec

RETURN, diam*lumdist*1000 /(1+z)^2




END
