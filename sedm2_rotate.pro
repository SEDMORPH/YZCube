;;** INPUT: filedir
;;
;;** OUTPUT: norien
;;
;;******************************************************************

FUNCTION SEDM2_ROTATE, filedir,norien,seed=seed

  norien = 7

  rot_init = randomu(seed,3)*2d*!DPI ;3 random angles between 0 and 2PI
  
  rotate = fltarr(3,3,norien)+1
  rotate[*,*,0] = VW_FUNC_ROTATE([0.0,0.0,0.0]) ;face-on
  rotate[*,*,1] = VW_FUNC_ROTATE(rot_init)
  rotate[*,*,2] = VW_FUNC_ROTATE(rot_init+[!DPI/2,0,0]) ;rotate around X
  rotate[*,*,3] = VW_FUNC_ROTATE(rot_init+[!DPI,0,0])
  rotate[*,*,4] = VW_FUNC_ROTATE(rot_init+[3*!DPI/2.,0,0])
  rotate[*,*,5] = VW_FUNC_ROTATE(rot_init+[0.,!DPI/2.,0])
  rotate[*,*,6] = VW_FUNC_ROTATE(rot_init+[0.,3*!DPI/2.,0])
  
  return, rotate

END
