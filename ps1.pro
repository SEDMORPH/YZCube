pro ps1,filename,xsize=xsize,ysize=ysize,xoffset=xoffset, yoffset=yoffset

set_plot,'ps'

;;-- Following instructions here: http://www.astrobetter.com/blog/2009/04/23/making-fonts-better-in-idl-postscript-output/
!p.font = 0

if n_elements(xsize) eq 0 then xsize=20 ;xsize=9  ;; this is ridiculously small
if n_elements(ysize) eq 0 then ysize=xsize/1.5

if n_elements(xoffset) eq 0 then xoffset=0
if n_elements(yoffset) eq 0 then yoffset=0

device,file=filename,/color,bits_per_pixel=8, xsize=xsize, ysize=ysize

!p.thick=(!x.thick=(!y.thick=3))
!p.charthick=3


end
