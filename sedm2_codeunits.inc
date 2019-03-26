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


;;========= old, original, un-interpolated model=========
;n_Z=7
;Z_models = replicate({keys:'', values:0.0},n_Z)
;Z_models.keys=['22', '32', '42', '52', '62', '72', '82']
;Z_models.values=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05, 0.1]
;;========= old, original, un-interpolated model=========



;;========= interpolated model=========
n_Z=14
Z_models = replicate({keys:'', values:0.0},n_Z)
Z_models.keys=['22', '32', '42', '52', '55', '58', '62', '63', '64', '65', '67', '68', '72', '82']
Z_models.values=[0.000100, 0.000400, 0.004000, 0.008000, 0.012000, 0.016000, 0.020000, 0.023000, 0.026000, 0.030000, 0.035000, 0.040000, 0.050000, 0.100000 ]
;;========= interpolated model=========



