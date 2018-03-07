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


;;-- Metallicity dictionary, use orderedhash instead
;Z_hash=ORDEREDHASH(['m22', 'm32', 'm42', 'm52', 'm62', 'm72', 'm82'],$
;                [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05, 0.1])
;;NOTE, IDL 8.2 does not support hash or orderedhash, use to separate array instead
;Z_keys=['m22', 'm32', 'm42', 'm52', 'm62', 'm72', 'm82']
;Z_values=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05, 0.1]
;n_Z = n_elements(Z_values)

n_Z=7
Z_models = replicate({keys:'', values:0.0},n_Z)
Z_models.keys=['m22', 'm32', 'm42', 'm52', 'm62', 'm72', 'm82']
Z_models.values=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05, 0.1]
