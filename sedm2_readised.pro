FUNCTION SEDM2_READISED, filename, t, lambda
openr, 1, filename, /f77_unformatted
	nseds = 0L
	readu, 1, nseds                 ;nb of SEDS in this file
close, 1
openr, 1, filename, /f77_unformatted
	t = fltarr(nseds)
	readu, 1, nseds, t               ;time steps
	nlambda = 0L
	readu, 1, nlambda               ;nb of wavelength points
close, 1
openr, 1, filename, /f77_unformatted
	readu, 1, nseds, t
	lambda = fltarr(nlambda)
	readu, 1, nlambda, lambda        ;wavelength array
	seds = fltarr(nlambda, nseds)
	n = 0L
	for i= 0, nseds-1 do begin
		flux = fltarr(nlambda)
		readu, 1, nl, flux            ;each sed; units Lsol/AA & Lsol = 3.826 x 10^33 ergs/s
		seds[*, n] = flux
		n = n + 1
	endfor
close, 1
RETURN, SEDS
END
