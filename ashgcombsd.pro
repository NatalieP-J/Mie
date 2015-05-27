; Written May 26 2015 for SURP 2015 project
; this code produces all sky maps of scattered light from the draco
; nebula using a Henyey-Greenstein phase function for given values of
; the anisotropy parameter g

; heights of frankie illumination maps to use - must be strings for
; reading purposes
heights = ['0.1','0.2','0.5','1'] ; in kpc
; loop over all heights
for h=1, leng(heights)-1 do begin
height = heights[h]

; calculate nside based on height
nside = findnside(height = height)

; direction to Draco nebula from earth
dracoell = 92.24*!dtor
dracobee = 38.43*!dtor

; distance to the Draco nebula from earth
dracodist = [0.5,1.0] ; kpc

; distance to the galactic center from earth
GCdist = 8.5 ; kpc
; coordinates of galactic center in earth's frame
bGC = 0
lGC = 0

; calculate the coordinates of earth in Draco's frame, and the
; scattering angle to earth
scatter = pointingangles(dracodist, GCdist, dracobee, dracoell, bGC, lGC)

learth = scatter['learth']
bearth = scatter['bearth']
scatter = scatter['scatter']
polearth = !dpi/2. - bearth

; convert earth's Galactic coordinates to a pixel location
ang2pix_ring, nside, polearth, learth, earthpix

; calculate the number of pixels
npix = nside2npix(nside)
; calculate domega
domega = 4.*!pi/npix
; create an array of pixel indices
ipix = long(findgen(npix))
; find incident radiation field
skyin = findskyin(height = height)

; find galactic coordinates for each pixel
pix2ang_ring, nside, ipix, pols, ells
; convert polar angle to b angle
bees = !dpi/2. - pols

; create a list of files holding g parameter information
spawn, 'ls *_g.idlsav', filelist
; crop '_g.idlsav' from file names
fnames = filelist.remove(-9)

; cycle over each filename
for i=0, leng(filelist)-1 do begin
fname = fnames[i]
restore, filelist[i]

; create arrays to hold output scattered and isotropic light maps
sumscat = fltarr(npix)
isoscat = fltarr(npix)
timer, /start
; for each pixel of input radiation
for k=0, npix-1 do begin
; calculate the scattering angles from that input radiation pixel to
; each other pixel
uscats = ct(ellillum=ells,beeillum=bees, ellscat=ells[k], beescat=bees[k])
; calculate the phase function at desired values
pf = hg(g = g, u = uscats)
; find the isotropic scattering
iso = hg(g=0, u=uscats)
; add the values for this pixel to total radiation field
sumscat += pf*skyin[k]
isoscat += iso*skyin[k]
endfor
timer, /stop
timer, /print

; choose output directory
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/AllSky/CarbonSilicate/'

; construct a title for the plot and the savefile name
pname = fname+'approx_hg_g'+string(g,format='(F0.2)')+'_h'+$
         height+'kpc'
ptitle = 'Henyey-Greenstein Scattering, g = '+string(g,format='(F0.2)')

; switches for various plot options
; rot - controls how much the output map is rotated by
rot = [0];,180]
; sub - obselete now, originally used to subtract the minimum value
;       from the map. a boolean switch: 0 = do nothing, 1 = do subtraction
sub = [0];,1]
; log - controls whether data is log scaled or not. a boolean switch:
;       0 = do nothing, 1 = do subtraction
log = [0,1]
; iso - controls whether data is normalized by the isotropic value or
;       not. a boolean switch: 0 = do nothing, 1 = do normalization
iso = [0,1]
; iterate over each possible option to produce all possible plots
for r=0, leng(rot)-1 do begin
for s=0, leng(sub)-1 do begin
for l=0, leng(log)-1 do begin
for j=0, leng(iso)-1 do begin
; create a new save file name
npname = pname
; create a new output directory
noutputdir = outputdir
; multiply scattered light by domega
nsumscat = sumscat*domega

; if isotropic normalization requested, modify the array, filename  and output
; directory accordingly
if iso[j] eq 1 then begin
npname += '_iso'
nsumscat /= isoscat*domega 
noutputdir += 'IsoNormal'
endif
; if rotation desired, modify the filename accordingly
if rot[r] ne 0 then npname += '_rot'+string(rot[r],format = '(I0)')
; if min value subtraction requested, modify the array, filename and
; output directory accordingly
if sub[s] eq 1 then begin
npname += '_sub'
nsumscat -= (min(nsumscat))
noutputdir += 'Sub'
endif
; if log scaling requested, modify the array, filename and
; output directory accordingly
if log[l] eq 1 then begin
npname += '_log'
nsumscat = alog10(nsumscat)
noutputdir += 'Log'
endif
; add file extension and output directory
npname += '.png'
npname = noutputdir + '/' + npname
; create and save the map
mollview, nsumscat, grat = [30,30], glsize = 1., rot = rot[r],$
          titleplot = ptitle, subtitle = sptitle, png = npname

; Uncomment if creating a postscript file to convert it to a pdf
;cgfixps, npname
;cgps2pdf, npname
;spawn, 'rm '+npname

endfor
endfor
endfor
endfor

endfor


endfor

end
