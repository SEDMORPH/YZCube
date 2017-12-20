; +
; NAME:
; 	FIND_PARTNUM_RATIO
;
; PURPOSE:
;       This procedure find the particle number ratio of first halo
;       to the second one from the given fileseq
;
; CALLING SEQUENCE:
;
; 	FIND_PARTNUM_RATIO, fileseq
;
; INPUTS:
; 	fileseq: the file sequence of the simulation output
;
;
; OUTPUT:
; 		ratio: the ratio of the particle numbers of the two halos
;
; MODIFICATION HISTORY:
; 	Written by:	Yirui Zheng 17th Nov. 2017
; -


FUNCTION FIND_PARTNUM_RATIO, fileseq

  ratio = 0

;;---split the fileseq to characters
  fileseqlist = []
  for i=0, strlen(fileseq) do $
    fileseqlist =[fileseqlist, strmid(fileseq, i , 1) ]
  ; print, fileseqlist

;;----if the fileseq starts with "2x", then it is the major merger
;;----therefore the ratio = 1
  if ( fileseqlist[0] eq '2') AND (fileseqlist[1] eq 'x') then $
    return, ratio = 1.0

;;Otherwise, this might be a minor merger, the smaller galaxy
;;should be indcated with a "p" follow the Hubble type
;;The mass ratio of the larger galaxy to the smaller one follows "p"
;;The first galaxy comes before the character "_"

    ;;Now search the fileseqlist for p
  ;SmallerFirst = boolean(1) ; a boolean type variable indicates whether the smaller galaxy comes first
	SmallerFirst = 1 ; a boolean type variable indicates whether the smaller galaxy comes first
  for i=0, strlen(fileseq) do begin
    if (fileseqlist[i] eq '_') then SmallerFirst = 0; boolean(0)
    if (fileseqlist[i] eq 'p') then begin
      ratio = float( strmid(fileseq, i+1))
      if ratio eq 0 then $
       print, "Unable to decide the particle number ratio, please check the fileseq!" & $
       print, "Set the ratio to 0. as an error symbol."
      if SmallerFirst then ratio = 1.0/ratio
      return, ratio
    endif
  end
        ; ; Note that the mass ratio may have several digits
        ; ratio_digits_num = 0
        ; while (fileseqlist[i+1] )

    ;;Note that the mass ratio can have two or more digits


;;----If neither of above two situations is met,
;;----it should be a major merger that the two galaxies belong to diffrerent Hubble type.
;;----There should be two '_' and no 'p' in the fileseq
  temp=where(fileseqlist eq '_', count)
  if count eq 2 then return, ratio=1.0

;;----None of above situations is met
;;----return the default 0 as an error symbol
  print, "Unable to decide the particle number ratio, please check the fileseq!"
  print, "Set the ratio to 0. as an error symbol."
  return, ratio
END


; +
; NAME:
; 	FIND_DIVIDING_POINT
;
; PURPOSE:
;       This procedure find the index of the first particle of the second halo
;       in the sorterd haloID list
;
; CALLING SEQUENCE:
;
; 	FIND_DIVIDING_POINT, sorted_idhalo, partnum_ratio=partnum_ratio
;
; INPUTS:
; 	sorted_idhalo: the ascending sorted id list of halo particles
;   partnum_ratio: the particle numbers ratio find the particle number ratio
;                   of first halo to that of the second one in the halo
;
;
; OUTPUT:
; 		Dividing_point: the index of the first particle of the second halo
;                     in the sorterd haloID list
;
; MODIFICATION HISTORY:
; 	Written by:	Yirui Zheng 19th Nov. 2017
; -


FUNCTION FIND_DIVIDING_POINT, sorted_idhalo, partnum_ratio=partnum_ratio

  if NOT(keyword_set(sorted_idhalo)) then $
    print, "please input the sorted_idhalo"

  ;;the default partnum_ratio is 1, meaning the two halo have same particle number
  if NOT(keyword_set(partnum_ratio)) then partnum_ratio = 1.0

  ;;Change the partnum_ratio to float type to ensure the accuracy
  partnum_ratio = float(partnum_ratio)
  ;print, "partnum_ratio in FIND_DIVIDING_POINT:", partnum_ratio

  ;;NOTE: some particles may not be read in the SEDM2_READSNAP, the tot_partnum
  ;;cannot be get by n_elements(sorted_idhalo)
  tot_partnum = ( sorted_idhalo[-1] - sorted_idhalo[0] + 1 )
  ;print, "tot_partnum",tot_partnum
  ;;the ID of real Dividing_point
  DP_id =round(sorted_idhalo[0] + tot_partnum*( partnum_ratio/(1.0+partnum_ratio) ))
  print, "Halo.id of the first particle of the second halo",DP_id


  ;;NOTE: the particle at the real Dividing_point may not be read into the halo
  ;;      in the SEDM2_READSNAP, we have to search the nearest one in that case
  temp = sorted_idhalo[-1]
  while ( DP_id LE temp ) do begin

    Dividing_point = WHERE( sorted_idhalo eq DP_id, count )
    if ( count eq 1 ) then return, Dividing_point

    DP_id++

  endwhile

  ;;If we cannot find any proper Dividing_point, return -1 as an error notice
  print, "Cannot find any proper Dividing_point, return -1 as an error notice"
  return, Dividing_point=-1


END


; +
; NAME:
; 	histcenter
;
; PURPOSE:
;       this function give the approximate center of data that is list 1D position
;
; CALLING SEQUENCE:
;
; 	histcenter, data, binsize= binsize
;
; INPUTS:
;
;
; OUTPUT:
;
; MODIFICATION HISTORY:
; 	Written by:	Yirui Zheng 20th Nov. 2017
; -
FUNCTION histcenter, data, binsize= binsize

  if NOT(keyword_set(binsize)) then binsize = 1.0
  ;;do not consider outskirt particles
  min = cgpercentiles(data, percentiles=0.1)
  max = cgpercentiles(data, percentiles=0.9)
  histcount = cgHistogram(data, binsize=binsize, min=min,max=max)

  count_max = max(histcount)
  center_index = where(histcount eq count_max)
  ;;note: sometimes we have several value of that histcount eq count_max
  ;;use the average to represent them
  center_index = mean(center_index)

  center_pos = min + center_index * binsize + 0.5*binsize

  return, center_pos

END




;+
; NAME:
;	SEDM2_FIND_CENTERS_MERGER
;
; PURPOSE:
;       This procedure finds the central coordinates of two colliding galaxies
;
; CALLING SEQUENCE:
;
;	SEDM2_FIND_CENTERS_MERGER, fname, indir=indir, separation=separation
;
; INPUTS:
;	fname: Filename to be read in
;
; OPTIONAL INPUT:
;		indir: Directory in which the file is located
;       Hubparam: Value of Hubble parameter/100 (default is 0.71)
;
; OUTPUT:
;		separation: a 7 dimensional array:
;		separation[0:2] = coordinates of center_1
;		separation[3:5] = coordinates of center_2
;		separation[6] = separation between centers (sqrt(center_1^2 + center_2^2))

; MODIFICATION HISTORY:
; 	Written by:	Joe Llama 9th April 2014 - adapted from SEDM_FIND_CENTERS_MERGER
;       Adapted by:     Vivienne Wild July 2016 - from hyp_find_centers_merger
;	Re-written by:  Yirui Zheng November 2017 - using another algorithms to divide the halo and decide the center
;-



FUNCTION SEDM2_FIND_CENTERS_MERGER, fname, fileseq=fileseq, indir=indir, halo=halo, silent=silent

;; Set everything up

;; If this procedure is run with halo already defined then we don't need to read the
;; file in again, i.e. if this is being called from hyp_output
if NOT(keyword_set(halo)) then begin
	if NOT(KEYWORD_SET(indir)) then indir='./'
	if file_test(indir+fname) eq 0 then $
			message, 'sedm2_find_centers_merger: Snapshot file not found'+fname
	SEDM2_READSNAP, indir+fname, halo=halo,/gethalo
endif

if NOT(KEYWORD_SET(fileseq)) then begin
  print, "Please set the fileseq for finding centers"
  return, separation=-1.
endif

NHalo  = n_elements(halo.x)

xhalo = halo.x
yhalo = halo.y
zhalo = halo.z

vxhalo = halo.vx
vyhalo = halo.vy
vzhalo = halo.vz

mhalo = halo.mass

idhalo = halo.id

pothalo = halo.pot

sort_id_halo=sort(idhalo)

sorted_idhalo=idhalo(sort_id_halo)
sorted_pothalo=pothalo(sort_id_halo)
sorted_xhalo=xhalo(sort_id_halo)
sorted_yhalo=yhalo(sort_id_halo)
sorted_zhalo=zhalo(sort_id_halo)



ratio = FIND_PARTNUM_RATIO(fileseq)
print, "Particle number ratio of first halo to the second one:", ratio
;the index of the dividing point
DP_index = FIND_DIVIDING_POINT(sorted_idhalo, partnum_ratio=ratio)
idhalo_1=sorted_idhalo(0:DP_index-1L)
pothalo_1=sorted_pothalo(0:DP_index-1L)
xhalo_1=sorted_xhalo(0:DP_index-1L)
yhalo_1=sorted_yhalo(0:DP_index-1L)
zhalo_1=sorted_zhalo(0:DP_index-1L)

idhalo_2=sorted_idhalo(DP_index:Nhalo-1L)
pothalo_2=sorted_pothalo(DP_index:Nhalo-1L)
xhalo_2=sorted_xhalo(DP_index:Nhalo-1L)
yhalo_2=sorted_yhalo(DP_index:Nhalo-1L)
zhalo_2=sorted_zhalo(DP_index:Nhalo-1L)



binsize = 1.0
xc_1 = histcenter(xhalo_1, binsize= binsize)
yc_1 = histcenter(yhalo_1, binsize= binsize)
zc_1 = histcenter(zhalo_1, binsize= binsize)

xc_2 = histcenter(xhalo_2, binsize= binsize)
yc_2 = histcenter(yhalo_2, binsize= binsize)
zc_2 = histcenter(zhalo_2, binsize= binsize)

cen_1 = [xc_1, yc_1, zc_1]
cen_2 = [xc_2, yc_2, zc_2]


if not(keyword_set(silent)) then begin
   print,'First center',cen_1
   print,'Second center',cen_2
endif

rad_separate=sqrt((cen_1(0)-cen_2(0))^2. + $
	              (cen_1(1)-cen_2(1))^2. + $
	              (cen_1(2)-cen_2(2))^2. )

separation = fltarr(7)
separation[0:2] = cen_1
separation[3:5] = cen_2
separation[6]   = rad_separate

return,separation

END
