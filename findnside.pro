function findnside
; Read in map information
; Source map directory
ISRFDIR='/cita/d/raid-project/hp/pgmartin/GALPROP/FITS/ISRF/Standard/'

; sky with spectrum integrated over certain filter
; choice of 0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30 for height
;filein = ISRFDIR+'Standard_8.5_0_0.5_Filter_HP.fits.gz'
height = 0
radius = 8.5
sheight = string(height,format = '(I0)')
sradius = string(radius,format = '(F0.1)')
filein = ISRFDIR+'Standard_'+sradius+'_0_'+sheight+'_Filter_HP.fits.gz'

; read in fits file. filein obviously filename, 1 is fileunit, hfoof
; is the name of the extension to read (this is stored as a keyword)
foofilt = mrdfits(filein, 1, hfoof)
nside = sxpar(hfoof, 'NSIDE')
return, nside
end
