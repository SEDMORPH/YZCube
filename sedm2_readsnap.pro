;+
; NAME:
;	SEDM2_READSNAP
;
; PURPOSE:
;       This procedure reads a Gadget3 hdf5 snapshot file into IDL
;       structures
;
; CALLING SEQUENCE:
;
;	SEDM2_READSNAP, Filename
;
; INPUTS:
;	Filename: name and directory of Gadget file to be read in
;

; KEYWORD PARAMETERS:
;       GETGAS:  read in and return gas structure
;       GETSTARS: read in and return stars structure
;       GETHALO: read in and return DM halo structure
;
;       STAR: structure containing all newly formed star particle information
;       GAS: structure containing all gas particle information
;       HALO: structure containing all dark matter particle information
;       SFR: smoothed SFR of entire simulation: col1 = time, col2=sfr (Msol/yr)
;       SNAP_TIME: physical time of snapshot (Gyr)
;
; Note that structures return as -1 if no particles are present
;
;
; EXAMPLE:
;       SEDM2_READSNAP, filename,  stars=stars, gas=gas, halo=halo, sfr=sfr, snap_time=snap_time
;
; NOTES:
; All input is in code units
; All output is in physical units
;
; positions=kpc/h to convert to kpc pos=pos/hubparam
; velocities=km/s
; masses=10^10 M_sun/hubparam to convert to 10^10 M_sun  mass=mass/hubparam
; time=0.97781307/hubparam Gyr to convert to Gyr time=time*0.97781307/hubparam

; Particles are ordered in the following order:
; Gas,DM,Disk,Bulge,Newly formed Stars, BHs
;
; MODIFICATION HISTORY:
; 	Written by:	Vivienne Wild July 2016
;-


PRO SEDM2_READSNAP, filename,  stars=stars, gas=gas, halo=halo, sfr=sfr, snap_time=snap_time,getstars=getstars, getgas=getgas,gethalo=gethalo

;;-- check files exist
  if file_test(filename) eq 0 then message, 'Snapshot file not found'+ filename

  indir = (strsplit(filename,'snap',/extract,/regex))[0]
  SFR_logfile = indir+'sfr.txt'
  if file_test(sfr_logfile) eq 0 then message, 'SFR file not found'+ SFR_logfile

;;------------------------------------------------------------------
;;-- Read HDF5 file. To see what is there, use h5_browser
;;------------------------------------------------------------------

  @sedm2_codeunits.inc          ; import code units from include file

  smooth_time=100               ; nbins

;;-- read header information
  if strmatch(filename,'~*') eq 1 then filename = file_Search(filename) ;have to do this because h5f_open can't cope with ~
  file_id = H5F_OPEN(filename)

  AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'NumPart_ThisFile')
  NumPart = H5A_READ(AttrID)
  ngas = NumPart[0]
  nhalo = NumPart[1]
  nstars = NumPart[4]

  AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
  Snap_Time = H5A_READ(AttrID)*SimnUnitTime/hubparam


;;-- DM particles

  if NHalo gt 0 and keyword_set(gethalo) then begin
     halo = replicate({x:0.0, y:0.0, z:0.0, mass:0.0, vx:0.0, vy:0.0, vz:0.0, $
                       id:0L,pot:0.0}, Nhalo)
     DMCoordinates = H5D_READ(H5D_OPEN(file_id, 'PartType1/Coordinates'))
     DMParticleIDs =  H5D_READ(H5D_OPEN(file_id, 'PartType1/ParticleIDs'))
     DMVelocity =  H5D_READ(H5D_OPEN(file_id, 'PartType1/Velocities'))
     DMMasses = H5D_READ(H5D_OPEN(file_id, 'PartType1/Masses'))
     DMPotential = H5D_READ(H5D_OPEN(file_id, 'PartType1/Potential'))

     halo.x=reform(DMCoordinates[0,*])/hubparam
     halo.y=reform(DMCoordinates[1,*])/hubparam
     halo.z=reform(DMCoordinates[2,*])/hubparam

     halo.vx=reform(DMVelocity[0,*])
     halo.vy=reform(DMVelocity[1,*])
     halo.vz=reform(DMVelocity[2,*])

     halo.id = DMParticleIDs
     halo.mass = MassUnit*DMMasses/hubparam

     halo.pot = DMpotential
  endif else halo = -1


;;-- gas particles
  if Ngas gt 0 and keyword_set(getgas) then begin
     gas = replicate({x:0.0, y:0.0, z:0.0,vx:0.0, vy:0.0, vz:0.0,mass:0.0, sph_smooth:0.0,sfr:0.0, id:0L, rho:0.0, nh:0.0,metal:fltarr(12),ind_vel:0L, ind_Z:0},Ngas) ;ind_Z added for Metallicity, we have 7 models in current BC03

     GasCoordinates =  H5D_READ(H5D_OPEN(file_id, 'PartType0/Coordinates'))
     GasVelocity =  H5D_READ(H5D_OPEN(file_id, 'PartType0/Velocities'))
     GasIDs =  H5D_READ(H5D_OPEN(file_id, 'PartType0/ParticleIDs'))
     GasZ = H5D_READ(H5D_OPEN(file_id, 'PartType0/Metallicity')) ;12 values per particle
     GasSFR = H5D_READ(H5D_OPEN(file_id, 'PartType0/StarFormationRate'))
     GasSPHsmooth = H5D_READ(H5D_OPEN(file_id, 'PartType0/SmoothingLength'))
     GasDensity = H5D_READ(H5D_OPEN(file_id, 'PartType0/Density'))
     GasMasses = H5D_READ(H5D_OPEN(file_id, 'PartType0/Masses'))
     GasNeutralH =  H5D_READ(H5D_OPEN(file_id, 'PartType0/NeutralHydrogenAbundance'))

     gas.x=reform(GasCoordinates[0,*])/hubparam
     gas.y=reform(GasCoordinates[1,*])/hubparam
     gas.z=reform(GasCoordinates[2,*])/hubparam

     gas.vx=reform(GasVelocity[0,*])
     gas.vy=reform(GasVelocity[1,*])
     gas.vz=reform(GasVelocity[2,*])

     gas.rho = MassUnit*GasDensity*hubparam^2 ; gas density in 10^10 M_sun h^2 / kpc^3 (Ngas)
     gas.nh = GasNeutralH                     ; Neutral gas fraction (Ngas)
     gas.metal = MassUnit*GasZ/hubparam
     gas.sfr = GasSFR
     gas.id = GasIDs
     gas.sph_smooth = GasSphSmooth/hubparam ;kpc/h
     gas.mass = MassUnit*GasMasses/hubparam
  endif else gas = -1

;;-- Star particles
  if Nstars gt 0 and keyword_Set(getstars) then begin
     stars = replicate({x:0.0, y:0.0, z:0.0, vx:0.0, vy:0.0, vz:0.0, mass:0.0,metal:fltarr(12), age:0.0, met:0.0, id:0L,ind_vel:0L, ind_ssp:0L, ind_Z:0},Nstars) ;ind_ssp added for old stars in buildsed,;ind_Z added for Metallicity, we have 7 models in current BC03


     StarCoordinates =  H5D_READ(H5D_OPEN(file_id, 'PartType4/Coordinates'))
     StarVelocity =  H5D_READ(H5D_OPEN(file_id, 'PartType4/Velocities'))
     StarIDs =  H5D_READ(H5D_OPEN(file_id, 'PartType4/ParticleIDs'))
     StarZ = H5D_READ(H5D_OPEN(file_id, 'PartType4/Metallicity')) ;12 values per particle
     ;; note that we use initial mass here, because the SSPs that we use account
     ;; for the decline in stellar mass due to mass loss
     StarMasses = H5D_READ(H5D_OPEN(file_id, 'PartType4/InitialMass')) ;initial mass (also in file is current "masses")
     StarAge =  H5D_READ(H5D_OPEN(file_id, 'PartType4/StellarFormationTime'))

     stars.x=reform(StarCoordinates[0,*])/hubparam
     stars.y=reform(StarCoordinates[1,*])/hubparam
     stars.z=reform(StarCoordinates[2,*])/hubparam

     stars.vx=reform(StarVelocity[0,*])
     stars.vy=reform(StarVelocity[1,*])
     stars.vz=reform(StarVelocity[2,*])

     stars.id = StarIDs
     stars.metal = MassUnit*StarZ/hubparam
     stars.mass = MassUnit*StarMasses/hubparam
     stars.age = snap_time - StarAge*SimnUnitTime/hubparam

  endif else stars = -1

;;--

  H5F_CLOSE, file_id

;;------------------------------------------------------------------
;;-- SFR (not smoothed, just needed to normalise SFHs)
;; note that SFR_tot_Msun/yr is already in physical units
;; this is the sum over the SFR of the gas particles

  readcol, SFR_logfile, sfr_time, SFR_tot_Msun_yr, form='(F,X,X,F,X)',/silent
  sfr_time = sfr_time*SimnUnitTime/hubparam ;Gyr


  sfr = fltarr(2,n_elements(sfr_time))
  sfr[0,*] = sfr_time
  sfr[1,*] = smooth(SFR_tot_Msun_yr,smooth_time) ;boxcar smooth, median doesn't work

END
