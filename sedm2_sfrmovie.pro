PRO SEDM2_SFRMOVIE, dir_in, dir_out

  moviefile_sfr = dir_out+'movie_sfr.gif'

  @sedm2_codeunits.inc

  filename = file_search(dir_in+'/*.hdf5',count=nsnap)

;;------------------------------------------------------------------
;;-- make SFR movie
;;------------------------------------------------------------------
  SFR_logfile  = dir_in+'sfr.txt'
  readcol, SFR_logfile, sfr_time, SFR_tot_Msun_yr, form='(F,X,X,F,X)',/silent
  sfr_time = sfr_time*SimnUnitTime/hubparam ;Gyr

  smooth_time=250               ; nbins
  SFR_tot_Msun_yr_smooth = smooth(SFR_tot_Msun_yr,smooth_time)

  window,0,xsize=500,ysize=500
  !x.thick=(!y.thick=2)
  !p.charthick=2
  !x.charsize=(!y.charsize=(!p.charsize=1.2))

  ;;-- get time between snapshots. Silly way to do it, but I
  ;;   don't see another safe way.
    file_id = H5F_OPEN(filename[0])
    AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
    Snap_Time0 = H5A_READ(AttrID)*SimnUnitTime/hubparam

    file_id = H5F_OPEN(filename[150])
    AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
    Snap_Time150 = H5A_READ(AttrID)*SimnUnitTime/hubparam
    snap_deltatime = (Snap_Time150 - Snap_Time0)/150 ;in Gyr

  for index=0,nsnap-1 do begin
     print, index
     time = (index+1)*snap_deltatime+min(sfr_time)
     select=where(sfr_time lt time)

     wset,0
     ind = where(sfr_time gt min(sfr_time)+0.1 and sfr_time lt max(sfr_time)-0.1)

     plot,sfr_time[select], SFR_tot_Msun_yr_smooth[select],/ylog,yr=[0.001,1000],/xs,xtitle='Time / Gyr',position=pos,color=cgcolor('black'),background=cgcolor('white'),xr=minmax(sfr_time)
     xyouts, 0.04,0.25,textoidl('Star formation rate / M_\odot yr^{-1}'),align=0.0, orient=90,color=0,/normal
     image = cgsnapshot(Filename=moviefile_sfr, /GIF,/multiple, repeat_count=0, delay_time=20,/nodialog)

  endfor
  write_gif, moviefile_sfr, /close


END
