; +
; NAME:
; 	SEDM2_CENTERSLIST
;
; PURPOSE:
; 	This procedure reads in all snapshots of a certain fileseq,
;   calcuates the centers of the two halos for all snapshots.
;   The results are printed to a .txt file.
;
; CALLING SEQUENCE:
;
; 	SEDM2_CENTERSLIST, fileseq
;
;
; INPUTS:
; 	Filename: the fileseq of the simulaion
;
;
; KEYWORD PARAMETERS:
; 	INDIR:	Directory where file is stored. Default is current
; 	        directory.
;
; 	OUTDIR:	Directory where output file should be placed. Default
; 	        is input directory.
;
;
; EXAMPLE:
;       SEDM2_HYP_OUTPUT, 'Sc_Scp2_07',                                        $
;                    INDIR='/Volumes/data/gadget_simulations/Sc_Scp2_07/',     $
;                    OUTDIR='/Volumes/data/gadget_simulations/Sc_Scp2_07/'
;
;
; DEPENDENCIES:
;        sedm2_readsnap
;        sedm2_find_centers_merger
;
;
; NOTES:
; MODIFICATION HISTORY:
; 	Written by:	Yirui Zheng, Nov. 2017
; -



PRO SEDM2_CENTERSLIST, fileseq, indir=indir, outdir=outdir, have_BH=have_BH


  if not (keyword_set(indir))  then cd, current=indir
  if not (keyword_set(outdir)) then outdir = indir

;make sure that indir has the '/' on the directory structure
  if (strmid(indir, 0, 1, /reverse_offset) ne '/') then $
     indir = indir+'/'

  filename =  file_search(indir+'*_???.hdf5',count=Nsnap)

  for i=0,Nsnap-1 do $
     filename[i] = (strsplit(filename[i],'/',/extract,count=n))[n-1]
  print, 'SEDM2_HYP_OUTPUT:', string(Nsnap, format='(I03)') +' files found'


  cen_list = fltarr(6, Nsnap)

  for i=0,Nsnap-1 do begin

     if KEYWORD_SET(have_BH) then begin
	 SEDM2_READSNAP, indir+filename[i], BHs=BHs,/getBH, /no_sfr_log
	 nBH = n_elements(BHs.id)
	 if nBH eq 2 then begin ; BH not mergered yet
	     min_BHid = min(BHs.id)
	     min_idx = where(BHs.id eq min_BHid)
	     min_idx = min_idx[0]
	     max_BHid = max(BHs.id)
	     max_idx = where(BHs.id eq max_BHid)
	     max_idx = max_idx[0]
	     cen_list[*,i]= [ BHs[min_idx].x, BHs[min_idx].y, BHs[min_idx].z, BHs[max_idx].x, BHs[max_idx].y, BHs[max_idx].z]
	 endif
	 if nBH eq 1 then begin ;BH mergered
	     cen_list[*,i]= [ BHs[0].x, BHs[0].y, BHs[0].z, BHs[0].x, BHs[0].y, BHs[0].z ]
	 endif
     endif else begin
	 ;;----We only need halo for finding centers
	 SEDM2_READSNAP, indir+filename[i], halo=halo,/gethalo, /no_sfr_log
	 separation=SEDM2_FIND_CENTERS_MERGER(halo=halo, fileseq=fileseq)

	 cen_list[*, i] = separation[0:5]
     endelse


  endfor


  openw, lun, outdir+'centers.txt', /get_lun
  printf, lun, cen_list
  free_lun, lun



END
