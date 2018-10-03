;;*** SEDM2_PCA
;;*
;;* AIM: calculate SDSS-PCs from integrated CSP-spectra, all
;;  orientations (or just 1 orientation in nodust case)
;;*
;;* INPUT: spectra files
;;*
;;******************************************************************

PRO SEDM2_PCA, dir_in, dir_out, tauv, mu_d, dir_pca_data

  @sedm2_codeunits.inc
  ; @sedm2_directories.inc


;;-- set up plotting file
  outstr = '_tauv'+string(tauv,form='(F0.1)')
  outstr = outstr+'_mu'+string(mu_d,form='(F0.1)')

  cell = ''
  ; cell = 'cen_'
  ; cell = 'tracked_cell_'

  style = ''
  ; style = '_eagle'
  ; style = '_eagle_minus'
  ; style='_star_age'
  filename = file_search(dir_out+cell+'spec'+outstr+'_???'+style+'.fits',count=nsnap)

  ;;for code test
  ; nsnap=5

  outfile = dir_out+cell+'pcs'+outstr+style+'.fits'
  psfile = dir_out+cell+'pcs'+outstr+style+'.ps'

;;-- SFR for EQW Halpha
  SFR_logfile  = dir_in+'sfr.txt'
  readcol, SFR_logfile, sfr_time, SFR_tot_Msun_yr, form='(F,X,X,F,X)',/silent
  sfr_time = sfr_time*SimnUnitTime/hubparam ;Gyr
  smooth_time = 100                                ;fairly insensitive to this
  sfr_smooth = smooth(SFR_tot_Msun_yr,smooth_time) ;boxcar smooth, median doesn't work

  tau_young_ha = mu_d*tauv*( (5500./6564.)^0.7) + (1-mu_d)*tauv*( (5500./6564.)^1.3)


  ;;-- get time between snapshots. Silly way to do it, but I
  ;;   don't see another safe way.
    SimnSnaps = file_search(dir_in+'/*.hdf5')
    file_id = H5F_OPEN(SimnSnaps[0])
    AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
    Snap_Time0 = H5A_READ(AttrID)*SimnUnitTime/hubparam

    file_id = H5F_OPEN(SimnSnaps[150])
    AttrID= H5A_OPEN_NAME(H5G_OPEN(file_id, '/Header'), 'Time')
    Snap_Time150 = H5A_READ(AttrID)*SimnUnitTime/hubparam
    snap_deltatime = (Snap_Time150 - Snap_Time0)/150 ;in Gyr
    print, snap_deltatime

;;------------------------------------------------------------------
;;-- calculate SDSS PCs
;;------------------------------------------------------------------

;75. for Stelib/BC03; 26. for IndoUS/BC03; 58. for MILES/BC03
  pcs_tau = fltarr(3,nsnap)
  pcs_notau = fltarr(3,nsnap)
  norm = fltarr(nsnap)

  hd_tau = (hd_notau =  (d4n_tau = (d4n_notau = (eqw_ha = (sfr = fltarr(nsnap))))))

  print, "SEDM2_PCA: Number of snapshots to be processed: ", nsnap

  for i=0,nsnap-1 do begin

     data = mrdfits(filename[i],1,hdr,/silent)

     wave = data.wave
     airtovac,wave

     specstr_tau = {specarr:data.spec_tau,wave:wave,data_disp:58.}
     specstr_notau = {specarr:data.spec_notau,wave:wave,data_disp:58.}

     pcs_notau[*,i] = bc03_projectbc03(specstr_notau,3,25,dir_pca_data, plotspec=0,norm=nn,/usenormgappy,/silent)
     pcs_tau[*,i] = bc03_projectbc03(specstr_tau,3,25,dir_pca_data, plotspec=0,norm=nn,/usenormgappy,/silent)
     norm[i] = nn               ;needs different normalisations!!!

     hd_tau[i] = sedm2_spec_inds(21, wave, data.spec_tau)
     hd_notau[i] = sedm2_spec_inds(21, wave, data.spec_notau)
     d4n_tau[i] = sedm2_spec_inds(101, wave, data.spec_tau)
     d4n_notau[i] = sedm2_spec_inds(101, wave, data.spec_notau)

     cont = data.spec_tau[(where(wave ge 6564))[0]]*3.839*10d^33 ;in erg/s/AA

     time_now = min(sfr_time)+snap_deltatime*i
     ind = where(sfr_time ge time_now-0.001 and sfr_time le time_now) ;past 1d6 years
     sfr[i] =  mean(sfr_smooth[ind])
     eqw_ha[i] = sfr[i]*exp(-tau_young_ha)  * 10d^42 / 7.9 /cont ;erg/s / erg/s/AA -> AA units

  endfor

  ;;------------------------------------------------------------------
  ;;-- PLOTS
  ;;------------------------------------------------------------------

  color=indgen(nsnap)           ;for plotting

  ps1,psfile
  !p.multi=0



  ;;-- PC12 with DR7 galaxies

  restore,dir_pca_data+'DR7/DR7specobj_gal_pcs_noduplicates.sav'
  ind = where(sdss_uniq.pc1err gt 0 and sdss_uniq.pc1err lt 0.25)
  pcs_data = transpose([[sdss_uniq[ind].pc1], [-sdss_uniq[ind].pc2]])
  minx=-7 & miny=-3 & dbinx=0.2 & dbiny=0.1
  hist = hist_nd(pcs_data,[dbinx,dbiny],min=[minx,miny],max=[2,1.5]) ;I think this is doing Min -> Max+binsize
  loghist = alog10(hist)
  ind = where(finite(loghist) eq 0)
  loghist[ind]=-999.

  xx = findgen((size(hist,/dim))[0])*dbinx+minx + dbinx/2.
  yy = findgen((size(hist,/dim))[1])*dbiny+miny + dbiny/2.
  cgloadct,0,/reverse
  cgcontour, loghist,xx,yy,nlevels=20,/fill,missingvalue=-999.,xtitle='PC1',ytitle='PC2',xr=[-7,2],/xs,yr=[-3,2],/ys
  save, xx, yy, loghist, file='DR7_loghist.sav'

  oplot, pcs_notau[0,*],-pcs_notau[1,*],psym=-1
  loadct, 13,ncolors=nsnap
  for i=0,nsnap-1 do plots, pcs_tau[0,i],-pcs_tau[1,i],color=color[i],psym=1

  ;;-- Hdelta-D4000 with DR7 galaxies

  ;data = mrdfits(dir_pca_data+'SDSS_MPA/gal_idxfix_dr7_v5_2.fit',1) ;very slow, need to par down
  ;data2 = replicate({lick_hd_a_sub:0.0,lick_hd_a_sub_err:0.0,D4000_N_SUB:0.0,D4000_N_SUB_ERR:0.0},n_elements(data))
  ;copy_struct, data,data2
  ;data = data2
  ;save, data,file=dir_pca_data+'SDSS_MPA/HdA_D4n_dr7_v5_2.sav'
  restore, dir_pca_data+'SDSS_MPA/HdA_D4n_dr7_v5_2.sav'

  ind = where(data.lick_hd_a_sub_err gt 0 and data.lick_hd_a_sub_err lt 2)
  indices = transpose([[data[ind].D4000_N_SUB],[data[ind].LICK_HD_A_SUB]])
  minx=0.8 & miny = -8 & dbinx = 0.025 & dbiny = 0.3
  hist = hist_nd(indices,[dbinx,dbiny],min=[minx,miny],max=[2.5,12]) ;I think this is doing Min -> Max+binsize
  loghist = alog10(hist)
  ind = where(finite(loghist) eq 0)
  loghist[ind]=-999.

  xx = findgen((size(hist,/dim))[0])*dbinx+minx + dbinx/2.
  yy = findgen((size(hist,/dim))[1])*dbiny+miny + dbiny/2.
  cgloadct, 0,/reverse
  cgcontour, loghist,xx,yy,nlevels=20,/fill,missingvalue=-999.,xtitle='Dn4000',ytitle=textoidl('H\delta_A'),yr=[-8,12],/xs,xr=[0.7,2.5],/ys
  oplot, d4n_notau,hd_notau,psym=-1
  loadct, 13,ncolors=nsnap
  for i=0,nsnap-1 do plots, d4n_tau[i],hd_tau[i],color=color[i],psym=1


;;-- Hdelta Halpha
  plot, sfr_time, sfr_smooth,xtitle='time',ytitle='SFR',/ylog,yr=[0.01,max(sfr_smooth)]
  oplot, findgen(nsnap)*snap_deltatime+min(sfr_time),sfr,color=cgcolor('red')

  plot, Hd_tau, eqw_ha,psym=-3,yr=[0,50],ytitle='EQW_Ha (emission)',xtitle='EQW_Hd (absorption)',title='Dust attenuated'
  for i=0,nsnap-1 do plots, Hd_tau[i],eqw_ha[i],color=color[i],psym=1

  plot, sfr, eqw_ha,xtitle='SFR', ytitle='EQW_Ha (emission)',title='Dust attenuated',psym=-3,/xlog,xr=[0.1,max(sfr_smooth)],/ylog, yr=[1,max(eqw_ha)]
  for i=0,nsnap-1 do plots,sfr[i], eqw_ha[i],psym=1,color=color[i]

  ps2

  ;;------------------------------------------------------------------
  ;;-- OUTPUT
  ;;------------------------------------------------------------------

  struc = {norm:norm, pcs_tau:pcs_tau, pcs_notau:pcs_notau, hd_tau:hd_tau, hd_notau:hd_notau,  d4n_tau:d4n_tau, d4n_notau:d4n_notau}
  mwrfits, struc, outfile,/create

END
