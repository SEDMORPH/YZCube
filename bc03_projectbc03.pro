;;*** PROJECTBC03
;
; NOTE: - probably better to use ../VO/projectvo_bc03.pro if storing
;         files as it keeps track of files
;
; AIM: - to calculate pc amplitudes for model spectra
;
; INPUT: specstr = structure containing:
;           specarr - array of fluxes
;           wave    - wave array VACUUM
;           data_disp - best guess dispersion of spectra (see below)
;        nrecon = no. of eigenspectra to project
;        runno  = eigenspectra runno to use
;
; OUTPUT: returns PC Amplitudes
;
; OPTIONAL OUTPUT: norm - the normalisation used for each spectrum
;                  flxarr - rebinned, convolved flux (before
;                           normalisation) for wave_extend wavelength range
;
; OPTIONAL INPUT: plotspec - binary array of whether to plot or not
;                      (will always plot 20 whether you like it
;                      or not!)
;                 vdisp    - use a different velocity dispersion to
;                            the eigenspectra
;
;
;
;
; NOTES:
; *** input data_disp is velocity dispersion of model spectra. Use
;     75. for Stelib/BC03; 26. for IndoUS/BC03; 58. for MILES/BC03
;
;*****************************************************************

FUNCTION bc03_projectbc03, specstr, nrecon, runno,$ ;inputs
                           plotspec=plotspec,silent=silent,maskem=maskem,$ ;keyword switches
                           usegappy=usegappy,usenormgappy=usenormgappy, $ ;keyword switches
                           norm=norm,vdisp=vdisp,flxarr=flxarr,$ ; optional outputs
                           chi2=chi2,wave_espec=wave_espec, wave_extend = wave_extend ;optional outputs

; @bc03dir.inc
; @voinfo.inc

if NOT(TAG_EXIST(specstr,'data_disp')) then begin
    print, 'Need to know data dispersion if input own array'
    print, 'SPECSTR = {specarr, wave, data_disp, outfile}'
    return,-1
endif else data_disp = specstr.data_disp

if (size(specstr.specarr))[0] eq 2 then ngal = (size(specstr.specarr,/dim))[1] else ngal = 1

if N_ELEMENTS(plotspec) ge 1 then begin
    if N_ELEMENTS(plotspec) eq 1  then begin
        if plotspec eq 1 then plotspec = lonarr(ngal)+1 else if plotspec eq 0 then plotspec = lonarr(ngal)
    endif
endif else begin
    plotspec = lonarr(ngal)    ;plot first 20 at least
    if ngal gt 20 then plotspec[0:19]=1. else plotspec=plotspec+1
    if ngal eq 1 then plotspec=0 ;if doing them 1 by 1 turn it off
endelse

if keyword_set(maskem) and not(keyword_set(usegappy)) and not(keyword_set(usenormgappy)) then begin
    if not(keyword_set(silent)) then print,'PROJECTBC03: for /maskem I am setting /normgappy!'
    usenormgappy=1
endif

;; if runno has been input as a string convert it
if size(runno,/type) eq 7 then begin
    if not(keyword_set(silent)) then splog,'converting runno to interger'
    runno=uint(runno)
endif

strrunno = string(runno,form='(I0)')

if not(keyword_set(silent)) then splog,'Projecting runno',runno
print,"dir_vopca: ",dir_vopca
print, "strrunno: ", strrunno
print, "csv: ", dir_vopca+'pcavo_info.csv'
;; get eigenspectra
infile = DIR_VOpca+'ESPEC/pcavo_espec_'+strrunno+'.sav'
restore, infile

if not(keyword_set(silent)) then splog, 'Number of eigenspectra being used:',nrecon

espec = espec[*,0:nrecon-1]
wave_espec = wave
minwave = min(wave_espec)
maxwave = max(wave_espec)
npix_espec = (size(espec))[1]

;; extend espectra wave array for convolution of incoming data
minwave_log = round(alog10(minwave-100)*10000.)/10000.
maxwave_log = round(alog10(maxwave+100)*10000.)/10000.
npix = round((maxwave_log-minwave_log)/0.0001 +1)
wave_extend = 10D^(round((0.0001D*dindgen(npix) + minwave_log)*10000)/10000.)

;; get extra info for this PCA run
;; collect useful information on eigenspectra
nlines = file_lines(dir_vopca+'pcavo_info.csv')
openr,1,dir_vopca+'pcavo_info.csv'
tags = ''
readf,1,tags
data = strarr(nlines-1)
readf,1,data
close,1
especinfo = csv_str(tags,data,/silent)
ind = where(especinfo.runno eq runno)
if ind[0] eq -1 then begin
    print,'PROJECTBC03:no info for this runno, returning'
    return,0
endif
especinfo = especinfo[ind]

if strmatch(especinfo.norm,'*mean flux*') then begin
    norm_w0 = minwave
    norm_w1 = maxwave
endif else begin
    norm_w0 = uint(strmid(especinfo.norm,10,4,/reverse_offset))
    norm_w1 = uint(strmid(especinfo.norm,5,4,/reverse_offset))
    if norm_w0 eq 0 or norm_w1 eq 0 then begin
        print,'PROJECTBC03: normalisation info wrong, returning'
        return,0
    endif
endelse
if strmatch(especinfo.extra,'*CaH masked*') then begin
    print,'Ability to mask CaH no longer available'
    return,-1
endif

if (min(wave_espec) ne norm_w0 or max(wave_espec) ne norm_w1) and keyword_set(usenormgappy) then begin
    print,'I think Normgappy only works for eigenspectra normalised over whole wave range'
    return,-1
endif

;; sort out input arrays
wave = specstr.wave
specarr = specstr.specarr
if n_elements(wave) eq 0 or $
  n_elements(specarr) eq 0 then begin
    splog,'SPECSTR not complete'
    return,-1
endif

;; Velocity dispersion at which the eigenspectra were made
;; if vdisp=0.0, no convolution of input spectra to PCA
if n_elements(vdisp) eq 0. then begin
    vdisp = especinfo.vdisp
    if not(keyword_set(silent)) then splog, 'vdisp of eigenspectra',vdisp
endif else if not(keyword_set(silent)) then splog, 'vdisp requested',vdisp
if not(keyword_set(silent)) then splog, 'vdisp of input spectra',data_disp


if vdisp ne 0.0 and vdisp gt data_disp then begin
    if not(keyword_set(silent)) then splog,'convolving spectra to match eigenspectra or requested vdisp'
    cspeed = 2.99792e5
    loglamtov = cspeed * alog(10)
    dw = alog10(wave_extend[1]) - alog10(wave_extend[0])

    vdisp_add = sqrt(vdisp^2 - data_disp^2) ; Deconvolve template resolution
    sigma_pix = vdisp_add / loglamtov / dw
    NOCONV=0
endif else NOCONV=1


if not(keyword_set(silent)) then splog, 'nb of galaxies to project:',ngal

;;------------------------------------------------------------------
;; PROJECT SPECTRA

if (where(plotspec eq 1))[0] ne -1 and ngal gt 1 then !p.multi=[0,1,5]

pcs = fltarr(nrecon,ngal)
norm = fltarr(ngal)

if ARG_PRESENT(flxarr) then flxarr = fltarr(n_elements(wave_extend),ngal)
if arg_present(chi2) then chi2 = fltarr(ngal)

for i=0L,ngal-1 do begin

    if not(keyword_set(silent)) then if i/5000. eq i/5000 then print,'PROJECTVO_BC03: on spectrum',i

    ;; rebin spectrum onto same as espec
    linterp, wave,specarr[*,i],wave_extend,flux

    ;; convolve with Gaussian
    if not(NOCONV) then flux = gconv(flux, sigma_pix)


    ;; mask "emission" lines - not sure where this should go
    if keyword_set(maskem) then begin
        masksize = 5.0
        error = fltarr(npix_espec)+1.

        ind = where(wave_espec gt 4102.9-masksize and wave_espec lt 4102.9+masksize)
        if ind[0] ne -1 then error[ind]=0.

        ind = where(wave_espec gt 3971.195-masksize and wave_espec lt 3971.195+masksize)
        if ind[0] ne -1 then error[ind]=0.

        ind = where(wave_espec gt 3890.151-masksize and wave_espec lt 3890.151+masksize)
        if ind[0] ne -1 then error[ind]=0.

        masksize = 2.5

        ind = where(wave_espec gt 3836.472-masksize and wave_espec lt 3836.472+masksize)
        if ind[0] ne -1 then error[ind]=0.

        ind = where(wave_espec gt 3798.976-masksize and wave_espec lt 3798.976+masksize)
        if ind[0] ne -1 then error[ind]=0.

;; in the real data we only mask this if an AGN. Actually it is strong
;; in starforming galaxies too - here mask it anyway
        ind = where(wave_espec gt 3869-5. and wave_espec lt 3869+5.)
        if ind[0] ne -1 then error[ind]=0.

    endif else error = fltarr(npix_espec)+1.



    if ARG_PRESENT(flxarr) then flxarr[*,i] = flux

    ind = where(wave_extend ge min(wave_espec) and wave_extend le max(wave_espec))
    flux = flux[ind]
    if n_elements(ind) ne n_elements(wave_espec) then stop



    ;; calculate normalisation
    ind = where(wave_espec ge norm_w0 and wave_espec le norm_w1 and error eq 1.)
    if ind[0] eq -1 then begin
        print, 'PROJECTBC03: something wrong with normalisation, returning'
        return,-1
    endif

    norm[i] = mean(flux[ind])
    flux2 = flux/norm[i]

    ;; calculate pcs
    if keyword_set(usegappy) then pcs[*,i] = vwpca_gappy(flux2,error,espec[*,0:nrecon-1],meanarr,/silent) $
    else if keyword_set(usenormgappy) then begin
        pcs[*,i] = vwpca_normgappy(flux,error,espec[*,0:nrecon-1],meanarr,norm=tmp)
        norm[i] = tmp
    endif else $
      pcs[*,i] = (flux2-meanarr) ## TRANSPOSE(espec[*,0:nrecon-1])


    if arg_present(chi2) then chi2[i] = total( (flux - (vwpca_reconstruct(pcs[*,i],espec[*,0:nrecon-1])+meanarr)*norm[i])^2)/float(npix_espec)

;plot to check
    if plotspec[i] eq 1 then begin
        if nrecon ge 3 then n=2 else n=nrecon-1
        plot, wave_espec,flux,/xstyle,xrange=minmax(wave_espec),title=string(pcs[0:n,i],form='('+string(n+1,form='(I0)')+'F7.2)'),/ynozero,psym=10
        oplot, wave_espec,(vwpca_reconstruct(pcs[*,i],espec[*,0:nrecon-1])+meanarr)*norm[i],color=cgcolor('red'),psym=10
    endif

endfor

return,pcs

END
