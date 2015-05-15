; direction to Draco nebula from earth
dracoell = 90*!dtor
dracobee = 39*!dtor

; distance to the Draco nebula from earth
dracodist = [0.5,1.0] ; kpc

; distance to the galactic center from earth
GCdist = 8.5 ; kpc

; illumination as a beam from the galactic center in radians
ellillum = !dpi-!dpi - atan((dracodist*cos(dracobee))/(GCdist))
beeillum = -atan((dracodist*sin(dracobee))/((GCdist)^2 + $
                 (dracodist*cos(dracobee))^2)^0.5)

; direction to earth from the Draco nebula
earthell = !dpi-dracoell
earthbee = -dracobee

; scattering angle to earth
scatter = acos(ct(ellillum=ellillum,beeillum=beeillum,$
                  ellscat=earthell,beescat=earthbee))

sourcedir = '/mnt/raid-project/hp/njones/Mie/miefiles/'

sds = ['compaSil', 'compLamC', 'wdaSil', 'wdGra']
colors = ['red', 'green']
particles = ['draine', 'zubko']

fnames = sdgenfilelist(sds = sds, colors = colors, particles = particles)

exten = '.phsd'
hsize = 3
colnames = ['angle','m2','m1','s21','d21','pf','pol']

outdata = genfiledata(sourcedir, fnames, exten, hsize, colnames)

data = outdata['data']
headers = outdata['headers']

keys = fnames

outputdir = '/home/njones/Dropbox/Mie/NatalieResults/RGcompare/'

i=0
while i lt leng(fnames) do begin
fname = fnames[i]

datadict = dictionary()
datadict[keys[i]] = data[keys[i], 'pf']/(4*!dpi)
datadict[keys[i+1]] = data[keys[i+1], 'pf']/(4*!dpi)
datadict['angle'] = data[keys[i], 'angle']
datadict['scatter'] = scatter

hgplot = [headers[keys[i],'g'], headers[keys[i+1],'g']]

typecheck = strpos(fname,'zubko')
if typecheck ne -1 then begin
   ptype = 'carbon'
endif
if typecheck eq -1 then begin
   ptype = 'silicon'
endif
ptitle = 'Phase functions for a size distribution of '+ptype+$
         ' particles'   
start = strpos(fname, 'red')
nfname = fname.Remove(start,start+2)
pname = nfname+'_rgphase'

plotoptions, datadict, outputdir, pname, ptitle, $
             hgplot = hgplot, /isoplot, /nplot, /ylog
i+=2
endwhile

end
