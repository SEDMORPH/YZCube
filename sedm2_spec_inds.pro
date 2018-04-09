;
;+
; NAME:
;       SPEC_INDS.PRO
;
; PURPOSE:
;       Calculate a variety of spectral indices from an input (flux) spectrum
;       and continuum for models (no pixel errors).
;
; CALLING SEQUENCE:
;
;	VAL = SPEC_INDS(J,X,F)
;
; INPUTS:
;       J = index indentification (or array)
;       X = wavelength scale
;       F = flux (or array of fluxes)
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       SED  = compute directly from SED (default)
;       LICK = compute on Lick system
;       FFN  = compute using fitting function
;
; OUTPUTS:
;       VAL = index value (or array of values)
;
; OPTIONAL OUTPUTS:
;
; INDEX IDs: 101 = B4000VN (Balogh et al.)
;            102 = D4000 (Bruzual et al.)
;
; REQUIREMENTS: spec_inds.inc
;-
;------------------------------------------------------------------

FUNCTION D4000, X, F,ngal,narrow=narrow,sdss=sdss,gc=gc,vw=vw

if keyword_set(gc) then begin
    print, 'Gorgas et al. D4000 not yet implimented'
    return, -1
endif


nbin = N_ELEMENTS(X)

if keyword_set(narrow) then begin
    ;; Balogh et al., 1999, ApJ, 527, 54
    wl_b1 = 3850.0  &  wl_b2 = 3950.0
    wl_r1 = 4000.0  &  wl_r2 = 4100.0
endif else if keyword_set(sdss) then begin
    ;;Stoughton et al.
    wl_b1 = 3751.0  &  wl_b2 = 3951.0
    wl_r1 = 4051.0  &  wl_r2 = 4251.0
endif else if keyword_set(vw) then begin
;    wl_b1 = 3850.0  &  wl_b2 = 3920.0
;    wl_r1 = 4000.0  &  wl_r2 = 4100.0
    wl_b1 = 3820.0  &  wl_b2 = 3920.0
    wl_r1 = 4050.0  &  wl_r2 = 4250.0
endif else begin
    ;; Bruzual, 1983, ApJ, 273, 105
    wl_b1 = 3750.0  &  wl_b2 = 3950.0
    wl_r1 = 4050.0  &  wl_r2 = 4250.0
endelse

if min(x) gt wl_b1 then begin
    print, 'spectrum not blue enough',min(x)
    return,-1
endif

if max(x) lt wl_r2 then begin
    print, 'spectrum not red enough'
    return,-1
endif

;; integrate in F-nu units!
if NOT(KEYWORD_SET(SDSS)) then F2 = F*rebin(X^2,nbin,ngal) else F2=F

;; integrate
f_r = vw_integral(X, F2, wl_r1, wl_r2)
f_b = vw_integral(X, F2, wl_b1, wl_b2)

;; This gives exactly the same results, but doesn't do arrays
;f_r = integral(X, F2, wl_r1,wl_r2)
;f_b = integral(X, F2, wl_b1,wl_b2)

;; normalise
f_r = f_r/(wl_r2-wl_r1)
f_b = f_b/(wl_b2-wl_b1)

return, f_r/f_b

END

;;------------------------------------------------------------------

;; Taken from Jarle's code lineindex.pro
FUNCTION LINEINDEX, X, F, info, ngal

;; make sure continuum windows are in wave range
if min(X) gt info.wl_b1 or max(X) lt info.wl_r2 then return, lonarr(ngal)-1L

npix = n_elements(X)  ;;pixels number on the CCD of the spectrometer, i.e. the number of wavelength

;; Find average red & blue continuum fluxes
f_blue_cont = vw_integral(X, F, info.wl_b1, info.wl_b2) / (info.wl_b2 - info.wl_b1)
f_red_cont  = vw_integral(X, F, info.wl_r1, info.wl_r2) / (info.wl_r2 - info.wl_r1)

f_blue_cont = TRANSPOSE(f_blue_cont)
f_red_cont = TRANSPOSE(f_red_cont)

wl_blue = (info.wl_b1 + info.wl_b2) / 2.0
wl_red  = (info.wl_r1 + info.wl_r2) / 2.0

;; Compute continuum for the index (different for BH than Lick & DTT)
lick = strpos(info.name, 'Lick_') & lick = lick[0] ne -1
dtt  = strpos(info.name, 'DTT_') & dtt = dtt[0] ne -1
bh   = strpos(info.name, 'BH_') & bh = bh[0] ne -1

if bh then f_cont = (f_red_cont + f_blue_cont) / 2.0 $
else begin
    m      = (f_red_cont - f_blue_cont) / (wl_red - wl_blue)
    b      = f_blue_cont
    f_cont = rebin(m,npix,ngal) * rebin((X - wl_blue),npix,ngal) + rebin(b,npix,ngal)
endelse


;; Compute equivalent width

wli = [info.wl_1, info.wl_2]
if info.unit eq 'A' then val = vw_integral(X, 1 - F/rebin(f_cont,npix,ngal), wli[0], wli[1])


;; Compute magnitudes

if info.unit eq 'mag' then begin

  nflux = vw_integral(X, F/rebin(f_cont,npix,ngal), wli[0], wli[1])
  val = -2.5 * alog10(nflux / (wli[1] - wli[0]))

endif

return, val

END


;;******************************************************************
;; MAIN PROGRAM
;;******************************************************************

FUNCTION SPEC_INDS, J, X, F, SED=SED, LICK=LICK, FFN=FFN

;;contains structure "index"
@spec_inds.inc

;; Number of spectra/galaxies
;;F[*,i] is the sed of i th galaxy
if (SIZE(F))[0] eq 1 then ngal = 1 $
else ngal = (SIZE(F,/dim))[1]

if N_ELEMENTS(F)/ngal ne N_ELEMENTS(X) then stop

;; set up output array
;; (in future this can be a structure with values, errors, IDs and
;; names)
nindices = n_elements(J)
val      = fltarr(nindices,ngal)
n        = 0L


;;------------------------------------------------------------------

for i = 0L, nindices -1 do begin

    if J[i] eq 100 then val[n,*] = d4000(X,F,ngal,/vw)

    ;; B4000VN (narrow D4000 Balogh et al. 1999)
    ;; checked this matches BC03 stoc burst values (dust atten)
    if J[i] eq 101 then val[n,*] = d4000(X,F,ngal,/narrow)

    ;; D4000 (Bruzual, 1983, ApJ, 273, 105)
    if J[i] eq 102 then val[n,*] = d4000(X,F,ngal)

    ;; D4000 (Stoughton et al.)
    if J[i] eq 103 then val[n,*] = d4000(X,F,ngal,/sdss)

    ;; D4000 (Gorgos et al. NOT IMPLEMENTED)
    if J[i] eq 104 then val[n,*] = d4000(X,F,ngal,/gc)

    ;; extended Lick indices
    if J[i] lt 38 then val[n,*] = lineindex(X,F,index[J[i]],ngal)

    n=n+1L

endfor

RETURN, val

END
