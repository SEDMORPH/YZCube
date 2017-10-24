PRO SEDM2_SDSSMOVIE, dir_in, dir_out, redshift, tauv, mu_d, imagesize,orientation


;;-- set up input and output filenames
  outstr = '_z'+string(redshift,form='(F0.3)')
  outstr = outstr+'_tauv'+string(tauv,form='(F0.1)')
  outstr = outstr+'_mu'+string(mu_d,form='(F0.1)')
  outstr = outstr+'_'+string(imagesize,form='(I0)')+'arcmin'

  moviefile_image = dir_out+'movie_sdssimage'+outstr+'_or'+string(orientation,form='(I0)')+'.gif'
  
;;------------------------------------------------------------------
;;-- make image movie
;;------------------------------------------------------------------
  j = orientation
  
  if not(keyword_Set(noimage)) then begin
     filename = file_search(dir_out+'sdssimage'+outstr+'_???.fits',count=nsnap) 
     window,0,xsize=500,ysize=500  
     
     for i=0,nsnap-1 do begin

;        image_noise = mrdfits(outfile_fits,0)
        image_noise_dust = mrdfits(filename[i],1,/silent)

        norien = (size(image_noise_dust,/dim))[2]
        if j+1 gt norien then MESSAGE, 'GO BACK TO SDSSIMAGE AND COMPUTE MORE ORIENTATIONS (TURN OFF /FACEON KEYWORD)'
        
;        SEDM_RGBIM,  reform(image_noise[*,*,j,*]),[3,2,1],rgbim=rgbim
        SEDM2_RGBIM,  reform(image_noise_dust[*,*,j,*]),[3,2,1],rgbim=rgbim
        image = cgsnapshot(Filename=moviefile_image, /GIF,/multiple, repeat_count=0, delay_time=20,/nodialog)

        ;color_im = color_quan(RGBim[*,*,0], RGBim[*,*,1], RGBim[*,*,2], rmap, gmap, bmap)
        ;write_gif, moviefile_image, color_im, rmap, gmap, bmap, /multiple, repeat_count=0, delay_time=20
        
     endfor
     write_gif, moviefile_image, /close
  endif

END
