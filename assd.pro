; direction to Draco nebula from earth
dracoell = 92.24*!dtor
dracobee = 38.43*!dtor

; distance to the Draco nebula from earth
dracodist = [0.5,1.0] ; kpc

; distance to the galactic center from earth
GCdist = 8.5 ; kpc

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
; create arrays to hold numeric (1) and analytic (2) values for g
gpts1 = fltarr(leng(keys))
gpts2 = fltarr(leng(keys))

while i lt leng(fnames) do begin
fname = fnames[i]

; intialize data dictionary to feed to plotoptions
datadict = dictionary()
; add renormalized phase functions for both colors
datadict[keys[i]] = data[keys[i], 'pf']/(4*!dpi)
; add angle information and scattering angles
datadict['angle'] = data[keys[i], 'angle']

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
ccheck = strpos(fname, 'red')
if ccheck ne -1 then ctype = 'red'
if ccheck eq -1 then ctype = 'green'

; construct a title for the plot and the savefile name
ptitle = 'Phase functions for a size distribution of '+ptype+$
         ' particles under '+ctype+' light with '+rvtype   

pname = fname+'_allsky'

nside = findnside(height = '1')
npix = nside2npix(nside)
ipix = long(findgen(npix))

skyin = findskyin(height = '1')

angles = skyin

pixels = long(findgen(npix))
pix2ang_ring, nside, ipix, pols, ells
bees = !dpi/2. - pols

sumscat = skyin
timer, /start
for k=0, npix-1 do begin
uscats = ct(ellillum=ells[k],beeillum=bees[k], ellscat=ells, beescat=bees)
;pf = mie(acos(uscats),angles = datadict['angle'], pf = datadict[keys[i]])
;sumscat +=
endfor
timer, /stop
timer, /print
print, fnames[i]

i+=1
endwhile

end
