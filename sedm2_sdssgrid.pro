;; ANGDIAM
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

FUNCTION ANGDIAM, z, angsize,Omega_m=Omega_m, Lambda0=Lambda0, H0=H0,silent=silent

if N_ELEMENTS(z) eq 0 then begin
    print,'USAGE: size(proper kpc)=ANGDIAM( z, angsize,Omega_m=Omega_m, Lambda0=Lambda0, H0=H0)'
    return,0
endif

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
;;------------------------------------------------------------------
;; ANGDIAM_KPC
;;
;; AIM: determine the angular size of an obsject of proper physical
;; size (kpc) at a given redshift
;;
;; INPUT: z = redshift, flt or fltarr
;; OUTPUT: angular size in arcsec
;;
;; For inverse see angdiam_fib
;;******************************************************************

FUNCTION ANGDIAM_kpc, z, size,Omega_m=Omega_m, Lambda0=Lambda0, H0=H0,silent=silent

if N_ELEMENTS(z) eq 0 or n_elements(size) eq 0 then begin
    print,'USAGE: size(arcsec)=ANGDIAM_kpc( z, size[kpc],Omega_m=Omega_m, Lambda0=Lambda0, H0=H0)'
    return,0
endif

if NOT(KEYWORD_SET(Omega_m)) then Omega_m=0.3
if NOT(KEYWORD_SET(Lambda0)) then Lambda0=0.7
if NOT(KEYWORD_SET(H0)) then H0=70

lumdist = lumdist(z,Omega_m=Omega_m, Lambda0=Lambda0, H0=H0,silent=silent)

;; angsize = (360/2pi *60*60) * dl * (1+z)^2/(1000*lumdist)
;; lumdist is in Mpc
;; dl is proper size in kpc
;; angsize is in arcsec

diam = size*(1+z)^2 / lumdist / 1000d ;radians

RETURN, diam*(60d*60*360/(2*!PI))

END

;;******************************************************************
;;******************************************************************
;; output: xx,yy,nx,ny
;; input: redshift

PRO SEDM2_SDSSGRID,xx,yy,nx,ny,redshift, imagesize=imagesize, center=center, rotate=rotate
  
  @sedm2_codeunits.inc          ;H, Om, Ol

  if redshift eq 0 then begin
     print, 'sedm2_sdsssgrid does not work with z=0!'
     stop
  endif

  ;; pixel size
  size_pix_arcsec = 0.396             ;arcsec
  size_pix_kpc = angdiam(redshift, size_pix_arcsec,Omega_m=Omega_m, Lambda0=omega_l, H0=hubparam*100,/silent)
  
  ;; total box size required to be output. Needs to be ~500 for merger
  ;; simulations, but for now use much smaller size
;  size_box_kpc = 150.               ;kpc
;  size_box_arcsec = angdiam_kpc(redshift,size_box_kpc,Omega_m=Omega_m, Lambda0=omega_l, H0=hubparam)
 
  if not(keyword_set(imagesize)) then size_box_arcsec = 60   $ ;1 arcmin postage stamp
  else size_box_arcsec = 60*imagesize                          

  size_box_kpc = angdiam(redshift, size_box_arcsec,Omega_m=Omega_m, Lambda0=omega_l, H0=hubparam*100,/silent)

  print, 'SEDM2_SDSSGRID: size of box in arcsec and kpc: ',size_box_arcsec, size_box_kpc

  ;; set up square image grid
  nx = round(size_box_arcsec/size_pix_arcsec)
  xx = findgen(nx+1)*size_pix_kpc - size_box_kpc/2. ;gives cell edges, including final edge
  
  ny = nx
  yy = xx
  
  ;;----------set a moving box
  if not(KEYWORD_SET(center)) then center = [0, 0, 0]

  if (KEYWORD_SET(rotate)) then $
    center = transpose(rotate ## center)
    print, "rotated_center", center

  xx = xx + center[0]
  yy = yy + center[1]


END
