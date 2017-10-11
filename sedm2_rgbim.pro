PRO SEDM2_RGBIM,image_jansky,ind_filter,rgbim=rgbim

; set parameters
; footnote 7 in Lupton+2004: first set Q (nonlin) -> 0 and choose
; linear stretch (to depend on units of image)
; then adjust Q to bring out brighter features
; 4.5,5.7,7.8
scales= [1,1.3,2]*450000d
nonlinearity= 3

;ind = where(image lt 50)
;image[ind] = 0                  ;remove low SB that we'd never see 

tmp = size(image_jansky,/dimensions)
nx = tmp[0]
ny = tmp[1]

RGBim = fltarr(nx,ny,3)
RGBim[*,*,0] = image_jansky[*,*,ind_filter[0]]
RGBim[*,*,1] = image_jansky[*,*,ind_filter[1]]
RGBim[*,*,2] = image_jansky[*,*,ind_filter[2]]

; rebin
;RGBim= rebin(RGBim,floor(nx*resizefactor),(ny*resizefactor),3)

; scale and set colors
RGBim = nw_scale_rgb(RGBim,scales=scales)

RGBim = nw_arcsinh_fit(RGBim,nonlinearity=nonlinearity)

RGBim = nw_fit_to_box(RGBim,origin=origin)

RGBim = nw_float_to_byte(RGBim)

;RGBim_smooth = SMOOTH(RGBim,3)
cgimage, RGBim,/nointerp,multimargin=1,/keep_aspect_ratio


END
