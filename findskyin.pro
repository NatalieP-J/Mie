; Written May 21 2015 for SURP 2015 project
; this code imports an all sky radiation field and is not very flexible

function findskyin, height = height, radius = radius, ifilt = ifilt

; height = height above the midplane (a string)
; radius = radius from the galactic center (a string)
; ifilt  = ?

; set default parameter values
setdefaultvalue, height, '0.5'
setdefaultvalue, radius, '8.5'
setdefaultvalue, ifilt, 2

; radiation field map source directory
ISRFDIR='/cita/d/raid-project/hp/pgmartin/GALPROP/FITS/ISRF/Standard/'

; choice of 0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30 for height
filein = ISRFDIR+'Standard_'+radius+'_0_'+height+'_Filter_HP.fits.gz'

; read in fits file
foofilt = mrdfits(filein, 1, hfoof)
; find nside
nside = findnside(height = height, radius = radius)
; calculate the number of pixels
npix = nside2npix(nside)

; create array to hold map
skyfilt = dblarr(npix)
; assign values of the map
skyfilt[*] = (foofilt.spectra)[ifilt, *]

return, skyfilt
end
