forcegen = 1

height = '0.1'
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

npix = nside2npix(nside)
domega = 4.*!pi/npix
ipix = long(findgen(npix))
skyin = findskyin(height = height)

pixels = long(findgen(npix))
pix2ang_ring, nside, ipix, pols, ells
bees = !dpi/2. - pols

spawn, 'ls *_g.idlsav', filelist
fnames = filelist.remove(-9)

for i=0, leng(filelist)-1 do begin
fname = fnames[i]
restore, filelist[i]

sumscat = fltarr(npix)
isoscat = fltarr(npix)
timer, /start
for k=0, npix-1 do begin
uscats = ct(ellillum=ells,beeillum=bees, ellscat=ells[k], beescat=bees[k])
pf = hg(g = g, u = uscats)
iso = hg(g=0, u=uscats)
sumscat += pf*skyin[k]
isoscat += iso*skyin[k]
endfor
timer, /stop
timer, /print
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/AllSky/CarbonSilicate/'
pname = fname+'approx_hg_g'+string(g,format='(F0.2)')+'_h'+$
         height+'kpc'
ptitle = 'Henyey-Greenstein Scattering, g = '+string(g,format='(F0.2)')

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
noutputdir+='/'
npname += '.ps'
npname = noutputdir+npname
mollview, nsumscat, grat = [30,30], glsize = 1.,$
          ps = npname, $
          titleplot = ptitle

cgfixps, npname
cgps2pdf, npname
spawn, 'rm '+npname
endfor
endfor
endfor
endfor

endfor


end
