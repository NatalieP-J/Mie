; Written May 21 2015 for SURP 2015 project
; this code calculates nside and is not very flexible

function findnside, height = height, radius = radius

; height = height above the midplane (a string)
; radius = radius from the galactic center (a string)

; set default parameter values
setdefaultvalue, height, '0.5'
setdefaultvalue, radius, '8.5'

; radiation field map source directory
ISRFDIR='/cita/d/raid-project/hp/pgmartin/GALPROP/FITS/ISRF/Standard/'

; choice of 0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30 for height
filein = ISRFDIR+'Standard_'+radius+'_0_'+height+'_Filter_HP.fits.gz'

; read in fits file and extract nside keyword
foofilt = mrdfits(filein, 1, hfoof)
nside = sxpar(hfoof, 'NSIDE')
return, nside
end
