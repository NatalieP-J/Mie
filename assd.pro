forcegen = 1

height = '0.5'
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

scatter = pointingangles(dracodist, GCdist, dracobee, dracoell, bGC, lGC)
print, scatter
learth = scatter['learth']
bearth = scatter['bearth']
scatter = scatter['scatter']
polearth = !dpi/2. - bearth
lt0 = where(polearth lt 0)
;polearth[lt0] = 2*!dpi + polearth[lt0]

ang2pix_ring, nside, polearth, learth, earthpix
pix2vec_ring, nside, earthpix, earthvec
skyin = findskyin(height = '1')

l1 = 'earthvecx = ['
for v=0, leng(earthvec[0,*])-1 do begin
if v ne leng(earthvec[0,*])-1 then l1+=string(earthvec[0,v],format='(F0.9)')+','
if v eq leng(earthvec[0,*])-1 then l1+=string(earthvec[0,v],format='(F0.9)')
endfor
l1 += '] &'

l2 = 'earthvecy = ['
for v=0, leng(earthvec[1,*])-1 do begin
if v ne leng(earthvec[1,*])-1 then l2+=string(earthvec[1,v],format='(F0.9)')+','
if v eq leng(earthvec[1,*])-1 then l2+=string(earthvec[1,v],format='(F0.9)')
endfor
l2 += '] &'

l1 = 'earthvecx = [0,1,0] &'
l2 = 'earthvecy = [0,0,1] &'
mycommand = l1+$
            l2+$
            'plots, earthvecx, earthvecy, psym = 6'
mollview, skyin, /log,execute = mycommand
end
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
spawn, 'ls *.idlsav', sumscatfilelist

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
ptitle = 'Mie scattered light'
sptitle = 'for a size distribution of '+ptype+$
         ' particles under '+ctype+' light' 
if rvtype ne '' then sptitle += 'with '+rvtype   

pname = fname+'_allsky_h'+height+'kpc'

savecheck = where(sumscatfilelist eq fname+'_h'+height+'kpc.idlsav')

if ((forcegen ne 0) or (savecheck eq -1)) then begin
npix = nside2npix(nside)
ipix = long(findgen(npix))
skyin = findskyin(height = height)
;mollview, skyin, /log

angles = skyin

pixels = long(findgen(npix))
pix2ang_ring, nside, ipix, pols, ells
bees = !dpi/2. - pols

sumscat = fltarr(npix)
isoscat = fltarr(npix)
timer, /start
for k=0, npix-1 do begin
uscats = ct(ellillum=ells[k],beeillum=bees[k], ellscat=ells, beescat=bees)
pf = mie(acos(uscats),angles = datadict['angle'], pf = datadict[keys[i]])
iso = hg(g = 0, u = uscats)
sumscat += pf * skyin[k]
isoscat += iso * skyin[k]
endfor
timer, /stop
timer, /print
save, sumscat, filename = fname+'_h'+height+'kpc.idlsav'
save, isoscat, filename = fname+'iso_h'+height+'kpc.idlsav'
print, fnames[i]
endif

if ((forcegen eq 0) and (savecheck ne -1)) then begin
restore, filename = fname+'_h'+height+'kpc.idlsav'
restore, filename = fname+'iso_h'+height+'kpc.idlsav'
endif

rot = [0];,180]
sub = [0];,1]
log = [0];,1]
iso = [0];,1]
pname = outputdir+pname
for r=0, leng(rot)-1 do begin
for s=0, leng(sub)-1 do begin
for l=0, leng(log)-1 do begin
for j=0, leng(iso)-1 do begin
npname = pname
nsumscat = sumscat
if iso[j] eq 1 then begin
npname += '_iso'
nsumscat /= isoscat 
endif
if rot[r] ne 0 then npname += '_rot'+string(rot[r],format = '(I0)')
if log[l] eq 1 then begin
npname += '_log'
nsumscat = alog10(nsumscat)
endif
if sub[s] eq 1 then begin
npname += '_sub'
nsumscat -= min(nsumscat)
endif
npname += '.ps'
mollview, nsumscat, grat = [30,30], glsize = 1., rot = rot[r],$
          titleplot = ptitle, subtitle = sptitle, ps = npname;, $
          ;execute = 'plots, '

cgfixps, npname
cgps2pdf, npname
spawn, 'rm '+npname
endfor
endfor
endfor
endfor


i+=1
endwhile

end
