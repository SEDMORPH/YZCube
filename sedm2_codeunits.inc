;;-- cosmology
  hubparam = 0.71               
  Omega_m = 0.3                 
  Omega_L = 0.7

;;-- code units
  UnitTime_in_s = 3.08568e+16
  GigaYear_in_s=1.0d9*24.*3600.*365.2425
  Massunit=1d10
  UnitTime_in_Gyr=(UnitTime_in_s/GigaYear_in_s)/hubparam ;; VW 25/01/17 why /h here??
  SimnUnitTime = 0.97781307

;;-- let use a more robust method
  ;;;snap_deltatime = 0.02         ;time between snapshots

;;-- conversions
  Lsol_in_erg = 3.846d33
  pc_in_m = 3.08567758d16
  Mpc_in_cm =3.08567758d24
