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
;-



FUNCTION SEDM2_FIND_CENTERS_MERGER, fname, indir=indir, halo=halo, silent=silent

;; Set everything up 

;; If this procedure is run with halo already defined then we don't need to read the 
;; file in again, i.e. if this is being called from hyp_output
if NOT(keyword_set(halo)) then begin
	if NOT(KEYWORD_SET(indir)) then indir='./'
	if file_test(indir+fname) eq 0 then $
			message, 'sedm2_find_centers_merger: Snapshot file not found'+fname
	SEDM2_READSNAP, indir+fname, halo=halo,/gethalo
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

Nhalo_half=Nhalo/2.

idhalo_1=sorted_idhalo(0:Nhalo_half-1L)
pothalo_1=sorted_pothalo(0:Nhalo_half-1L)
xhalo_1=sorted_xhalo(0:Nhalo_half-1L)
yhalo_1=sorted_yhalo(0:Nhalo_half-1L)
zhalo_1=sorted_zhalo(0:Nhalo_half-1L)

idhalo_2=sorted_idhalo(Nhalo_half:Nhalo-1L)
pothalo_2=sorted_pothalo(Nhalo_half:Nhalo-1L)
xhalo_2=sorted_xhalo(Nhalo_half:Nhalo-1L)
yhalo_2=sorted_yhalo(Nhalo_half:Nhalo-1L)
zhalo_2=sorted_zhalo(Nhalo_half:Nhalo-1L)

sort_pot_1=sort(pothalo_1)

potsorted_id_1=idhalo_1(sort_pot_1)
potsorted_pot_1=pothalo_1(sort_pot_1)
potsorted_xhalo_1=xhalo_1(sort_pot_1)
potsorted_yhalo_1=yhalo_1(sort_pot_1)
potsorted_zhalo_1=zhalo_1(sort_pot_1)

sort_pot_2=sort(pothalo_2)

potsorted_id_2=idhalo_2(sort_pot_2)
potsorted_pot_2=pothalo_2(sort_pot_2)
potsorted_xhalo_2=xhalo_2(sort_pot_2)
potsorted_yhalo_2=yhalo_2(sort_pot_2)
potsorted_zhalo_2=zhalo_2(sort_pot_2)

cen_1=[median(potsorted_xhalo_1(0:49)),  $
	   median(potsorted_yhalo_1(0:49)),  $
	   median(potsorted_zhalo_1(0:49))]

cen_2=[median(potsorted_xhalo_2(0:49)),  $
       median(potsorted_yhalo_2(0:49)),  $
       median(potsorted_zhalo_2(0:49))]

if not(keyword_set(silent)) then begin
   print,'First center',median(potsorted_pot_1(0:49)),cen_1
   print,'Second center',median(potsorted_pot_2(0:49)),cen_2
endif

rad_separate=sqrt((cen_1(0)-cen_2(0))^2. + $
	              (cen_1(1)-cen_2(1))^2. + $
	              (cen_1(2)-cen_2(2))^2. )

print,rad_separate
 
separation = fltarr(7)
separation[0:2] = cen_1
separation[3:5] = cen_2
separation[6]   = rad_separate

return,separation

END
