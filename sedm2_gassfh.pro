;+
; NAME:
;	SEDM2_GASSFH
;
; PURPOSE:
;       This procedure calculates the star formation history of gas
;       and new star particles, based on their SFR during the course of the simulation, and saves them
;       to an output file
;
; CALLING SEQUENCE:
;
;	SEDM2_GASSFH, Fileseq, Indir
;
; INPUTS:
;	Fileseq: Core name of Gadget files to be read in
;
; OPTIONAL INPUT:
;       Indir: Directory of Gadget files Default is current directory.
;
; KEYWORD PARAMETERS:
;	OUTDIR:	Directory where file is stored. Default it INDIR
;
;
;
; EXAMPLE:
;       SEDM2_GASSFH, Fileseq, Indir, outdir=outdir
;       sedm2_gassfh, 'Galaxy_model_Sab', 'Galaxy_model_Sab/'
;
; NOTES:

; 1) Normalisation of SFR. Each individual star particle does not have the total mass of a star
; particle when we integrate over the SFH, because the time at which
; it turns into a star particle is randomised, to some extent. The spectra
; are based on the SFH, so the total mass that we extract from the
; spectra will be equivalent to SUM(SFH_gas+SFH_stars). We see from
; the final ratio calculated here that we will underestimate the total
; stellar mass of new stars in the galaxy, compared to SUM(stars.mass)
; by about 10% in the first ~30 snapshots. Because the gas particles
; are also forming stars, we do not renormalise the star particles to
; have the correct stellar mass, as this would overestimate the number
; of very young stars in the system. This may cause problems for accurate mass comparisons in the future...
;
; 2) Matching SSPs to particles. The snapshots are spaced by,
; typically, 2e7 yrs. We do not have finer resolution for the SFR(t)
; of the gas particles. Rather than assign SSP=0 to all t=0 particles,
; have assumed a top hat SFR(t) and allowed all SSPs with
; t-delta_t/2<t<t+delta_t/2 to contribute, weighted by the time
; interval between SSPs.
;
; Note that structures return as -1 if no particles are present
;
; ~9 mins total for 2xSab_13 (low res)
; 20.5MB each snapshot
;
; MODIFICATION HISTORY:
; 	Written by:	Vivienne Wild July 2016 based on hyp_gassfh
;-


PRO SEDM2_GASSFH, fileseq, indir, outdir=outdir, quiet=quiet,  $
                model_str=model_str, models_dir=models_dir

  @sedm2_codeunits.inc          ; import code units from include file, needed for SimnUnitTime and hubparam

  if not keyword_set(outdir) then outdir = indir
  if not keyword_set(model_str) then model_str = 'bc03'
  if n_elements(indir) eq 0 then indir='./'
  if not keyword_set(outdir) then outdir = indir

  filename_long = file_search(indir+'*_???.hdf5',count=Nsnap)

  print, 'SEDM2_GASSFH found: '+string(Nsnap, format='(I03)')+ ' files'
  if Nsnap eq 0 then message, 'No files found: '+indir+'*_???.hdf5'

  filename = strarr(nsnap)
  for i=0,Nsnap-1 do $
     filename[i] = (strsplit(filename_long[i],'/',/extract,count=n))[n-1]

  psfile = outdir+'sedm_gassfh.ps' ;output sfh plots

;;-- get time between snapshots. Silly way to do it, but I
;;   don't see another safe way.
  file_id = H5F_OPEN(filename_long[0])
  AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
  Snap_Time0 = H5A_READ(AttrID)*SimnUnitTime/hubparam

  ;; comment out the orignal block and use a sligtly unsafer way to test gassfh for a subset >>>>>>
  ;; do remember to change it back
  ; file_id = H5F_OPEN(filename_long[150])
  ; AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
  ; Snap_Time150 = H5A_READ(AttrID)*SimnUnitTime/hubparam
  ; delta_age = (Snap_Time150 - Snap_Time0)/150*1d9 ;in yr

  file_id = H5F_OPEN(filename_long[1])
  AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
  Snap_Time1 = H5A_READ(AttrID)*SimnUnitTime/hubparam
  delta_age = (Snap_Time1 - Snap_Time0)/1.0*1d9 ;in yr
  ;;<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< test subset


;;-- read ssps from SSP files to get time steps
;; TO DO: added to evernote
;; BC_str,'62' need to be passed in from calling routine so we can
;; change the SSPs. Ideally need to be a single input too.
;;===========================% Keyword MODELS_DIR not allowed in call to: SEDM2_GETSSPS
  ;tmp = SEDM2_GETSSPS(model_str,'62', models_dir=models_dir) ;origin
  tmp= SEDM2_GETSSPS(models_dir, model_str,'62')
  ;repalce the origin line with that find in sdssimage
;;===========================yrzheng
  age_ssp = tmp.age & tmp=0     ;in Gyr
  nssps = n_elements(age_ssp)
  ;delta_age_ssp = age_ssp[findgen(nssps-1)+1]-age_ssp[findgen(nssps-1)] ; in Gyr
;;===============yrzheng, Jan, 2018
;;====Change the division of the delta_age_ssp
;;====Let it match the "closest bin" algorithm
  delta_age_ssp = fltarr( n_elements(age_ssp) )
  delta_age_ssp[0] = (age_ssp[1] - age_ssp[0])/2.0
  delta_age_ssp[1:-2] = (age_ssp[2:-1] - age_ssp[0:-3])/2.0
  delta_age_ssp[-1] = (age_ssp[-1] - age_ssp[-2])/2.0
  ;print, n_elements(delta_age_ssp)


  age_snap = dblarr(Nsnap)

  time = systime(1)
  masscum1 = 0d
  sfr_track = fltarr(nsnap)

;;------------------------------------------------------------------
;;-- For every snapshot calculate and store SFH of gas and star
;;   particles in form of SSP indices
  for i=0, Nsnap-1 do begin ;; !!! must be run from 0 as sfr[] is built in loop!!!

;;-- build output filename
     outfile =  outdir+(strsplit(filename[i],'.',/extract))[0] + '_gassfh.sav'

;;-- read simulation files
     SEDM2_READSNAP, filename_long[i], gas=gas, stars=stars, /getstars, /getgas, $
                     snap_time=snap_time, sfr=sfr_global

     age_snap[i] = snap_time    ; in Gyr, starts at ~10.7Gyr

     if i eq 0 then oldstars_minid = min(stars.id) ;minimum ID number of old stars. New stars will have ID < this value.

;;-- initialise output arrays
     if size(gas,/type) eq 8  then begin
        ngas = n_elements(gas)
        ; gassfh = fltarr(ngas,nssps)
        gassfh = fltarr(ngas,nssps, n_Z) ;;funny face (n_Z). Okay, n_Z is the how many metallicity models we have now
     endif else ngas=0

     if size(stars,/type) eq 8 then begin
        ind_newstars = where(stars.id lt oldstars_minid, nnewstars,compl=ind_oldstars)
        if ind_oldstars[0] ne -1 then noldstars = n_elements(ind_oldstars) else noldstars = 0

        if nnewstars ne 0 then newstarsfh = fltarr(nnewstars,nssps, n_Z)
;        if noldstars ne 0 then oldstarsfh = fltarr(noldstars,nssps)
     endif else begin
        nnewstars=0
        noldstars=0
     endelse
     print, 'SEDM2_GASSFH: snapshot # old and new stars', i, noldstars, nnewstars


     ;Below we require that all snapshots are equally spaced in time. Use a tolerance of 10^-3 of delta_age.
;; Joe - possible error in last snapshot file
;;     if i gt 0 and abs((age_snap[i]-age_snap[i-1])*1e9 - delta_age) gt delta_age/100. then message, 'Snapshots must be equally spaced in time'

;;-- initialise array to store sfr and metallicity for all snapshots
     if i eq 0 then begin
       sfr = fltarr(ngas+nnewstars, Nsnap)
       ind_metal_history = uintarr(ngas, Nsnap) ;;record the ind_metal_history
       ind_metal_history = ind_metal_history + 999 ;;use 999 as symbol for unassigned data
     endif


;;-- if generations = 2 in Gadget-2 allvars.h then one gas particle
;;   makes 2 gas particles. See email from Peter 2/5/2014. VW checked and compared to mass
;;   of gas particles and this is actually the correct way to recover the correct gas particle ID:
;;   stars_id[ind_double] = stars[ind_double].id + long64(2)^(long64(31))
;;   currently not using this, as its very complicated to keep track of non-unique IDs!
     if nnewstars gt 0 then begin
        ind_double = where(stars[ind_newstars].id le 0)
        if ind_double[0] ne -1 then message, 'This code only works for 1:1 conversion of gas:star particles, set generations=1 in allvar.h'
     endif


;;-- age associated with this, and all preceding, snapshots.
     gas_age = snap_time - age_snap[0:i] ; in Gyr, for Snapshot 0 -> current

;;-- find closest SSP to the ages of all preceeding snapshots
     diff = abs(rebin(gas_age, i+1, Nssps,/sample) - rebin(transpose(age_ssp),i+1,Nssps) ) ; all in Gyr
     minarr = min(diff,dimension=2,ind_ssp_1d)
     ind_ssp_2d = array_indices(size(diff,/dim),ind_ssp_1d,/dimensions) ;convert back to 2D indices
     ind_ssp = reform(ind_ssp_2d[1,*])


;;------------------------------------------------------------------
;;-- Assign SFH of gas particles
;;------------------------------------------------------------------

;; note that sfr of snap0 = 0.0! not sure what to do

     if ngas gt 0 then begin ; some gas particles

; this should never happen!
;        junk = where(gas.id le 0,count) ;check indices are sensible
;        if count gt 0 then message, 'Gas indices screwed up'

;;-- Store the instantaneous SFR of each gas particle at this snapshot, also metallicity
        ;; this array is then used in subsequent snapshots to build up a SFR(t) over the duration of the simulations
        sfr[gas.id-1,i] = gas.sfr ; particles are in arbitrary places in gas array, here we order from 0 to Ngas+Nstars for easy access
        sedm2_z_ind, gas.metal, Z_models.values, Z_ind
        ind_metal_history[gas.id-1, i] = Z_ind
        if nnewstars gt 0 then begin
          newstars=stars[ind_newstars]
          sedm2_z_ind, newstars.metal, Z_models.values, Z_ind
          ind_metal_history[newstars.id-1, i] = Z_ind
        endif

;;-- Add the gas SFR during this snapshot + each preceeding snapshot at the correct SSP-index
        ;; the snapshots have delta_age of 2e7 years, here we spread the SFR out as a top hat into SSP bins with t-delta_t/2<t<t+delta_t/2
        ;; this has to be a cumulative for loop, otherwise non-unique ind_ssp's don't get counted
        gas_IMH = ind_metal_history[gas.id -1,*]
        for j=0,i do begin
          ;;SF_length is the time that particles keep forming stars
          ;;For the current snapshot, SF_length is only half of the delta_age
          if j lt i then SF_length = delta_age else SF_length = delta_age/2.0
          ; if SF_length lt delta_age then print,i,j

          for mm=0, N_z-1 do begin

             ; ind_metal = where(ind_metal_history[*,j] eq mm, mcount)
      	     ind_metal = where(gas_IMH[*,j] eq mm, mcount)
             ; help, ind_metal
             ; if j eq 0 then print, mm, n_elements(ind_metal)
             if mcount gt 0 then begin ;;do have gas particles with this metallicity
               ;all the SSPs between this+0.5 and this-0.5 snapshot and all particles attached this metal value
               ind_all = where(age_ssp*1e9 gt (i-j-0.5)*delta_age and age_ssp*1e9 lt (i-j+0.5)*delta_age, count)

               if count le 1 then $ ; only 1 ssp bin matches the age of the gas. Put all mass from this snapshot into the closest bin
                  gassfh[ind_metal,ind_ssp[j],mm] = gassfh[ind_metal,ind_ssp[j], mm]+ sfr[gas[ind_metal].id-1,j]*SF_length $
               else begin           ; many ssp bins lie within the delta_age of the snapshots. Share out mass across SSPs, weighted by SSP age bin size
                  delta_age_ssp_all = total(delta_age_ssp[ind_all])
                  for k=0,count-1 do gassfh[ind_metal,ind_all[k], mm] = gassfh[ind_metal,ind_all[k], mm]+sfr[gas[ind_metal].id-1,j]*SF_length*delta_age_ssp[ind_all[k]]/delta_age_ssp_all
                  ;tmp = 0. & for k=0,count-1 do tmp = tmp+total(sfr[gas.id-1,j]*delta_age*delta_age_ssp[ind_all[k]]/delta_age_ssp_all)
                  ;print, total(sfr[gas.id-1,j]*delta_age), tmp
               endelse
            endif ; if mcount gt 0
          endfor ; for mm
        endfor ; for j

;;--  check for mass conservation
        masscum1 = masscum1+total(gas.sfr)*delta_age
        masscum2 = total(gassfh)

        ;;-- check for young stars
        ind_tmp = where(age_ssp*1e3 le 45)
        print, 'gas mass in stars <45Myr:', total(gassfh[*,ind_tmp])
        ind_tmp = where(age_ssp*1e3 gt 45 and age_ssp*1e3 lt 90)
        print, 'gas mass in stars 45-90Myr:', total(gassfh[*,ind_tmp])
        ind_tmp = where(age_ssp*1e3 gt 90 and age_ssp*1e3 lt 140)
        print, 'gas mass in stars 90-140Myr:', total(gassfh[*,ind_tmp])
        ; ind_tmp = where(age_ssp gt 3 and age_ssp lt 13.7)
        ; print, 'gas mass in stars 3-13.7Gyr:', total(gassfh[*,ind_tmp])
        ; ind_tmp = where(age_ssp gt 13.7)
        ; print, 'gas mass in stars 13.7Gyr+:', total(gassfh[*,ind_tmp])

     endif else gassfh = -1

;;------------------------------------------------------------------
;;-- Assign SFH of star particles
;;------------------------------------------------------------------

;; some star particles are already in place at the start of the
;; simulation. These have an age that can be used to assign a single
;; SSP.

;; new star particles are formed from gas particles, and have the same ID
;; as the progenitor gas particle in order to keep track of their SFH

     if nnewstars gt 0 then begin ; some star particles

        ; newstars=stars[ind_newstars]
        stars_IMH = ind_metal_history[newstars.id -1,*]
;;-- Snapshot 0: these stars formed in first timestep, so they don't have a gas SFH. Assign SFR(now) = median of other gas particles
        if i eq 0 then sfr[stars[ind_newstars].id-1,i] = median(gas.sfr);stars.mass/delta_age
        if i eq 0 then print, "i=0"
        if i eq 0 then continue
;;-- Use the SFH of the progenitor gas particles to assign the SFH to the star particles
;;-- add the gas SFR during all preceeding snapshots at the correct SSP-index
        for j=0,i-1 do begin ;; for new stars, sfr[newstars.id-1, i] = 0, so we don't bother to run j=i

          ;;SF_length is the time that particles keep forming stars
          ;;For the current snapshot, SF_length is only half of the delta_age
          if j lt i then SF_length = delta_age else SF_length = delta_age/2.0
          ; if SF_length lt delta_age then print,i,j

          for mm=0, N_z-1 do begin
      	    ind_metal = where(stars_IMH[*,j] eq mm, mcount)
            ; ind_metal = where(ind_metal_history[*,j] eq mm, mcount)
            if mcount gt 0 then begin ;;do have gas particles with this metallicity
               ; ind_all = where(age_ssp*1e9 gt (j-0.5)*delta_age and age_ssp*1e9 lt (j+0.5)*delta_age,count)
               ind_all = where(age_ssp*1e9 gt (i-j-0.5)*delta_age and age_ssp*1e9 lt (i-j+0.5)*delta_age,count);all the SSPs between this+0.5 and this-0.5 snapshot

               if count le 1 then $ ; only 1 ssp bin matches the age of the gas. Put all mass from this snapshot into the closest bin
                  newstarsfh[ind_metal,ind_ssp[j],mm] = newstarsfh[ind_metal,ind_ssp[j],mm]+ sfr[newstars.id-1,j]*SF_length $
               else begin           ; many ssp bins lie within the delta_age of the snapshots. Share out mass across SSPs, weighted by SSP age bin size
                  delta_age_ssp_all = total(delta_age_ssp[ind_all])
                  for k=0,count-1 do  newstarsfh[ind_metal,ind_all[k],mm]=newstarsfh[ind_metal,ind_all[k],mm]+sfr[newstars[ind_metal].id-1,j]*SF_length*delta_age_ssp[ind_all[k]]/delta_age_ssp_all
               endelse

            endif ;if mcount gt 0
          endfor ; for mm
        endfor  ; for j

;;-- check for mass conservation
        masscum2 = masscum2+total(newstarsfh)

        ;;-- check for young stars
        ind_tmp = where(age_ssp*1e3 le 45)
        print, 'new star mass in stars <45Myr:', total(newstarsfh[*,ind_tmp])
        ind_tmp = where(age_ssp*1e3 gt 45 and age_ssp*1e3 lt 90)
        print, 'new star mass in stars 45-90Myr:', total(newstarsfh[*,ind_tmp])
        ind_tmp = where(age_ssp*1e3 gt 90 and age_ssp*1e3 lt 140)
        print, 'new star mass in stars 90-140Myr:', total(newstarsfh[*,ind_tmp])
        ; ind_tmp = where(age_ssp gt 3 and age_ssp lt 13.7)
        ; print, 'new star mass in stars 3-13.7Gyr:', total(newstarsfh[*,ind_tmp])
        ; ind_tmp = where(age_ssp gt 13.7)
        ; print, 'new star mass in stars 13.7Gyr+:', total(newstarsfh[*,ind_tmp])

     endif else newstarsfh = -1

     ;;-- check for gadget3 young stars
     ind = where(stars[ind_newstars].age*1e3 le 45)
     if ind[0] ne -1 then print, 'Star mass with G3 age<45Myr:', total(stars[ind_newstars[ind]].mass)
     ind = where(stars[ind_newstars].age*1e3 gt 45 and stars[ind_newstars].age*1e3 lt 90)
     if ind[0] ne -1 then print, 'Star mass with G3 age 45-90Myr:', total(stars[ind_newstars[ind]].mass)

;;------------------------------------------------------------------
;;-- output file

     save, gassfh, newstarsfh, ind_metal_history,sfr, file=outfile, /compress

     if not(keyword_set(quiet)) then splog, 'Saved file: ' + outfile

;;------------------------------------------------------------------
;;-- Basic checks on mass conservation and SFR. See Note above

     ;; print, i
     ;; print, i, ' total mass of new stars  ', total(stars[ind_newstars].mass), form='(I3.3, A30, E0.4)'
     ;; print, i, ' total sfr in gas ', masscum1, form='(I3.3, A30, E0.4)'
     ;; print, i, ' ratio sum_sfr_gas/M* ',  masscum1 / total(stars[ind_newstars].mass), form='(I3.3, A30, F0.2)'
     ;; print, i, ' Sum SFH ', total(newstarsfh)+total(gassfh), form='(I3.3, A30, E0.4)'
     ;; print, i, ' ratio sum_SFH/M* ',  (total(newstarsfh)+total(gassfh))/total(stars[ind_newstars].mass),form='(I3.3, A30, F0.2)'
     ;; print, i, ' %iles of sum_SFH*/M* ',vw_percentile(total(newstarsfh,2)/stars[ind_newstars].mass),form='(I3.3, A30, 3(F0.2,1X))'


;;-- check SFR of gas particles vs. that output directly by Gadget vs. SFH

;  NOTE!! BECAUSE OF THE UNEVENLY SPACED SSP BINS IT IS NOT POSSIBLE TO RECOVER THE SFR(t) FROM THE GASSFH+NEWSTARSFH FILE
;  see notes for more details
     sfr_track[i] = total(gas.sfr)
    ;  plot, sfr_global[0,*],sfr_global[1,*]
    ;  oplot, snap_time0+(findgen(nsnap))*delta_age/1d9,sfr_track, color=cgcolor('red'),psym=1,thick=4

  endfor

  ps1,psfile
  plot, sfr_global[0,*],sfr_global[1,*], title='SFH from sfr.txt [line] vs. total(gas.sfr) [red]'
  oplot, snap_time0+(findgen(nsnap))*delta_age/1d9,sfr_track, color=cgcolor('red'),psym=1,thick=4

  plot, sfr_global[0,*],sfr_global[1,*], title='SFH from sfr.txt [line] vs. total(gas.sfr) [red]',/ylog,yrange= [0.01,max(sfr_global[1,*])]
  oplot, snap_time0+(findgen(nsnap))*delta_age/1d9,sfr_track, color=cgcolor('red'),psym=1,thick=4
  ps2

  if not(keyword_set(quiet)) then splog,'total time taken (min): ', (systime(1)-time)/60.


END
