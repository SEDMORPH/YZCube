pro ps1c,filename

set_plot,'ps'
device,file=filename,/color,bits_per_pixel=8
device,/portrait,xoffset=2,yoffset=2,ysize=25,xsize=18
!p.charsize=(!x.charsize=(!y.charsize=1.2))


end
