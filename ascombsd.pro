forcegen = 0

height = '0.1'
print, 'height '+height + ' kpc'
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
;skyin = findskyin(height = '1')
;skyin[earthpix] = 1.
;mollview, skyin, /log
;end
; location of source files
sourcedir = '/mnt/raid-project/hp/njones/Mie/miefiles/'

spawn, 'ls *kpc.idlsav', sumscatfilelist
spawn, 'ls *_pf.idlsav', fnamelist

fnames = fnamelist.remove(-10)


; choose where to save plots
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/AllSky/CarbonSilicate/'

; initialize index
i=0
; create arrays to hold numeric (1) and analytic (2) values for g

while i lt leng(fnames) do begin
fname = fnames[i]

restore, fnamelist[i]
restore, fname+'_g.idlsav'

; intialize data dictionary to feed to plotoptions
datadict = dictionary()
; add renormalized phase functions for both colors
datadict[fname] = pf
; add angle information and scattering angles
datadict['angle'] = findgen(181)*!dtor

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
sptitle = 'for a size distribution of carbon and silicate'+$
         ' particles under '+ctype+' light' 
if rvtype ne '' then sptitle += ' with '+rvtype   

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
pf = mie(acos(uscats),angles = datadict['angle'], pf = datadict[fname])
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
npix = nside2npix(nside)
domega = 4.*!pi/npix

rot = [0];,180]
sub = [0];,1]
log = [0,1]
iso = [0,1]
for r=0, leng(rot)-1 do begin
for s=0, leng(sub)-1 do begin
for l=0, leng(log)-1 do begin
for j=0, leng(iso)-1 do begin
npname = pname
noutputdir = outputdir
nsumscat = sumscat*domega
if iso[j] eq 1 then begin
npname += '_iso'
nsumscat /= isoscat*domega 
noutputdir += 'IsoNormal'
endif
if rot[r] ne 0 then npname += '_rot'+string(rot[r],format = '(I0)')
if sub[s] eq 1 then begin
npname += '_sub'
nsumscat -= (min(nsumscat))
noutputdir += 'Sub'
endif
if log[l] eq 1 then begin
npname += '_log'
nsumscat = alog10(nsumscat)
noutputdir += 'Log'
endif
npname += '.ps'
npname = noutputdir + '/' + npname
mollview, nsumscat, grat = [30,30], glsize = 1., rot = rot[r],$
          titleplot = ptitle, subtitle = sptitle, ps = npname
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
