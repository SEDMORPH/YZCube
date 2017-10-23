;;-- rotation matrix
;;-- clockwise rotation
FUNCTION VW_FUNC_ROTATE, angles

  t_x = angles[0]               ;phi
  t_y = angles[1]               ;theta
  t_z = angles[2]               ;psi

;;http://en.wikipedia.org/wiki/Rotation_matrix

A = [ [cos(t_y)*cos(t_z), -cos(t_x)*sin(t_z)+sin(t_x)*sin(t_y)*cos(t_z), sin(t_x)*sin(t_z)+cos(t_x)],$
      [cos(t_y)*sin(t_z), cos(t_x)*cos(t_z)+sin(t_x)*sin(t_y)*sin(t_z) , -sin(t_x)*cos(t_z)+cos(t_x)],$
      [-sin(t_y)        , sin(t_x)*cos(t_y)                            , cos(t_x)*cos(t_y)] ]

return, A

END
