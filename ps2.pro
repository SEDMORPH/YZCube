pro ps2,noclean=noclean

if NOT KEYWORD_SET(noclean) then cleanplot,/silent
device,/portrait
device,/close
!p.font = -1
set_plot,'x'

return
end
