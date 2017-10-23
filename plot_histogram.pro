;; NORM: Sum_i=1^i=nbins P(n_i<x<m_i) * (m-n)_i = 1
;; PERUNITX: y-axis is nb per unit x axis
;; 
;; NOTE: IT DOESN'T MAKE MUCH SENSE TO HAVE FRACTION PER UNIT X     

PRO PLOT_HISTOGRAM,data,WEIGHTS=WEIGHTS,nbins=nbins,min=min,max=max,normalise=normalise,overplot=overplot,psym=psym,perunitx=perunitx,shade=shade,R=R,loc=x,hist=hist,orientation=orientation,color=color,linestyle=linestyle,spacing=spacing,position=position,xtitle=xtitle,ytitle=ytitle,title=title,xstyle=xstyle,ystyle=ystyle,binsize=binsize,yrange=yrange,xrange=xrange,thick=thick,ylog=ylog,vsurvey=vsurvey, compllim = compllim

;;-- can't remove bins because R won't work
hist = histogram(data,loc=x,reverse=R,nbins=nbins,min=min,max=max,binsize=binsize)
binsize = float(x[1]-x[0])

if n_elements(psym) eq 0  then psym=10

if n_elements(weights) eq 0 then begin
   
   if keyword_set(normalise) then hist = hist/float(total(hist))
   if KEYWORD_SET(perunitx) then  hist = hist/binsize

   if keyword_set(ylog) then begin
      hist = alog10(hist)
      ;ind = where(finite(hist) eq 0)
      ;if ind[0] ne -1 then hist[ind] = -999.
   endif
    
endif else begin
    
    hist2 = fltarr(nbins)
    meanweights = fltarr(nbins)
    for i=0L,nbins-1 do if R[i] ne R[i+1] then begin
       hist2[i]= total(weights[R[R[i]:R[i+1]-1]],/nan)
       meanweights[i] = median(weights[R[R[i]:R[i+1]-1]])
    endif

    if n_elements(vsurvey) ne 0 and n_elements(compllim) ne 0 then hist2[where(1/(meanweights*vsurvey) lt compllim)]=0.0

   if keyword_set(normalise) then hist = hist2/float(total(hist2)) $
   else if KEYWORD_SET(perunitx) then  hist = hist2/binsize else hist=hist2
   if keyword_set(ylog) then begin
      hist = alog10(hist)
      ;ind = where(finite(hist) eq 0)
      ;if ind[0] ne -1 then hist[ind] = -999.
   endif

 endelse

if NOT(KEYWORD_SET(overplot)) then plot,x+binsize/2.,hist,psym=psym,/xs,color=color,linestyle=linestyle,position=position,xtitle=xtitle,ytitle=ytitle,title=title,ystyle=ystyle,yrange=yrange,xrange=xrange,thick=thick $
else oplot, x+binsize/2.,hist,psym=psym,color=color,linestyle=linestyle,thick=thick


    if KEYWORD_SET(shade) then begin
        x = [x,x[nbins-1]+binsize]
        if n_elements(orientation) eq 0 then begin
           for i=0,nbins-1 do polyfill, [x[i],x[i],x[i+1],x[i+1]],[0.0,hist[i],hist[i],0.0],noclip=0,color=color
;           oplot, x+binsize/2.,hist,psym=10
        endif $
        else for i=0,nbins-1 do polyfill, [x[i],x[i],x[i+1],x[i+1]],[0.0,hist[i],hist[i],0.0],/line_fill,orien=orientation,noclip=0,color=color,linestyle=linestyle,spacing=spacing
    endif
    


END
