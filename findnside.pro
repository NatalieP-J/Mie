function findnside, height = height, radius = radius

setdefaultvalue, height, '0.5'
setdefaultvalue, radius, '8.5'
; Read in map information
; Source map directory
ISRFDIR='/cita/d/raid-project/hp/pgmartin/GALPROP/FITS/ISRF/Standard/'

; sky with spectrum integrated over certain filter
; choice of 0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30 for height
;filein = ISRFDIR+'Standard_8.5_0_0.5_Filter_HP.fits.gz'
filein = ISRFDIR+'Standard_'+radius+'_0_'+height+'_Filter_HP.fits.gz'

; read in fits file. filein obviously filename, 1 is fileunit, hfoof
; is the name of the extension to read (this is stored as a keyword)
foofilt = mrdfits(filein, 1, hfoof)
nside = sxpar(hfoof, 'NSIDE')
return, nside
end
