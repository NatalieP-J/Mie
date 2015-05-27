; Written May 25 2015 for SURP 2015 project
; this code produces all sky maps of scattered light from the draco
; nebula from phase functions of particle size distributions

; switch to force generation of scattered light arrays
; set to 0 -> do not force generation
; set to 1 -> force generation
forcegen = 0

; heights of frankie illumination maps to use - must be strings for
; reading purposes
heights = ['0.1','0.2','0.5','1'] ; in kpc
; loop over all heights
for h=0, leng(heights)-1 do begin
height = heights[h]
print, 'height '+height + ' kpc'

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
polearth = !dpi/2. - bearth]

; convert earth's Galactic coordinates to a pixel location
ang2pix_ring, nside, polearth, learth, earthpix

; location of source files
sourcedir = '/mnt/raid-project/hp/njones/Mie/miefiles/'

; chosen size distribution
sds = ['compaSil', 'compLamC', 'wdaSil31', 'wdGra31','wdaSil55', 'wdGra55']
; wavelengths for which scattering is calculated
colors = ['red', 'green']
; refraction indices used
particles = ['draine', 'zubko']

; create a list of filenames
fnames = sdgenfilelist(sds = sds, colors = colors, particles = particles)
spawn, 'ls *kpc.idlsav', sumscatfilelist

; specify file extension
exten = '.phsd'
; specify size of header
hsize = 3
; include column titles
colnames = ['angle','m2','m1','s21','d21','pf','pol']

; read in data from each file in fnames
outdata = genfiledata(sourcedir, fnames, exten, hsize, colnames)

; split output into column data and header information
data = outdata['data']
headers = outdata['headers']

; set keys to use in data dictionary to feed plot options
keys = fnames

; specify extension of auxiliary header file and its length
nexten = '.prsd'
nhsize = 16

; read in additional header information from auxiliary file and
; combine with the existing header
for i=0, leng(keys)-1 do begin
addheader = prheader(sourcedir+fnames[i]+nexten, nhsize)
headers[keys[i]] = headers[keys[i]]+addheader
endfor

; choose where to save plots
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/AllSky/'

; initialize index
i=0

while i lt leng(fnames) do begin
fname = fnames[i]

; intialize data dictionary to feed to plotoptions
datadict = dictionary()
; add renormalized phase functions for both colors
datadict[keys[i]] = data[keys[i], 'pf']/(4*!dpi)
; add angle information and scattering angles
datadict['angle'] = data[keys[i], 'angle']*!dtor

; find whether the particles are silicates or carbon
typecheck = strpos(fname,'zubko')
if typecheck ne -1 then ptype = 'carbon'
if typecheck eq -1 then ptype = 'silicate'
; find out which R_V value was used
acheck = strpos(fname, 'wd')
rvtype = ''
if acheck ne -1 then begin
rvcheck = strpos(fname, '31')
if rvcheck ne -1 then rvtype = textoidl('R_{V} = 3.1')
if rvcheck eq -1 then rvtype = textoidl('R_{V} = 5.5')
endif
; find out the colour of the light
ccheck = strpos(fname, 'red')
if ccheck ne -1 then ctype = 'red'
if ccheck eq -1 then ctype = 'green'

; construct a title for the plot and the savefile name
ptitle = 'Mie scattered light'
sptitle = 'for a size distribution of '+ptype+$
         ' particles under '+ctype+' light' 
if rvtype ne '' then sptitle += ' with '+rvtype   

; construct the output plot file name
pname = fname+'_allsky_h'+height+'kpc'

; check whether the required data has already been saved
savecheck = where(sumscatfilelist eq fname+'_h'+height+'kpc.idlsav')

; if force generation is on or the required data is not saved
if ((forcegen ne 0) or (savecheck eq -1)) then begin
; calculate the number of pixels
npix = nside2npix(nside)
; create an array of pixel indices
ipix = long(findgen(npix))
; find incident radiation field
skyin = findskyin(height = height)

; find galactic coordinates for each pixel
pix2ang_ring, nside, ipix, pols, ells
; convert polar angle to b angle
bees = !dpi/2. - pols

; create arrays to hold output scattered and isotropic light maps
sumscat = fltarr(npix)
isoscat = fltarr(npix)
timer, /start
; for each pixel of input radiation
for k=0, npix-1 do begin
; calculate the scattering angles from that input radiation pixel to
; each other pixel
uscats = ct(ellillum=ells[k],beeillum=bees[k], ellscat=ells, beescat=bees)
; interpolate the phase function to find its value at desired
; scattering angles
pf = mie(acos(uscats),angles = datadict['angle'], pf = datadict[fname])
; find the isotropic scattering
iso = hg(g = 0, u = uscats)
; add the values for this pixel to total radiation field
sumscat += pf * skyin[k]
isoscat += iso * skyin[k]
endfor
timer, /stop
timer, /print
; save the arrays for future use
save, sumscat, filename = fname+'_h'+height+'kpc.idlsav'
save, isoscat, filename = fname+'iso_h'+height+'kpc.idlsav'
print, fnames[i]
endif

; if force generation is off and the required file exists
if ((forcegen eq 0) and (savecheck ne -1)) then begin
; restore the relevant arrays
restore, filename = fname+'_h'+height+'kpc.idlsav'
restore, filename = fname+'iso_h'+height+'kpc.idlsav'
endif
; calculate domega
npix = nside2npix(nside)
domega = 4.*!pi/npix

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

i+=1
endwhile
endfor

end
