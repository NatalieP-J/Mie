; Written May 19 2015 for SURP 2015 project
; This file combines the phase functions for size distributions
; composed of different particle types produced from the same author, 
; and plots each color independently.

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

; produce scattering angles that correspond to a pointing from the
; nebula to earth
scatter = pointingangles(dracodist, GCdist, dracobee, dracoell, bGC, lGC)
scatter = scatter['scatter']

; location of source files
sourcedir = '/mnt/raid-project/hp/njones/Mie/miefiles/'

; chosen size distribution
sds = ['compaSil', 'compLamC', 'wdaSil31', 'wdGra31','wdaSil55', 'wdGra55']
; wavelengths for which scattering is calculated
colors = ['red', 'green']
; refraction indices used
particles = ['draine', 'zubko']

; create a list of file names
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

; set keys to use in data dictionary to feed to plot options
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
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/CarbonSilicate/'

; initialize index
i=0
k=0
; create arrays to hold numeric (1) and analytic (2) values for g
gpts1 = fltarr(leng(keys)/2)
gpts2 = fltarr(leng(keys)/2)

while i lt leng(fnames)-2 do begin
; skip segement so we are only combining data from the same authors
if ((i eq 2)||(i eq 6)) then begin
i+=2
endif

fname = fnames[i]

; initialize data dictionary to feed plotoptions
datadict = dictionary()
; find renormalized phase functions for both particle types
pf1 = data[keys[i], 'pf']/(4*!dpi)
pf2 = data[keys[i+2], 'pf']/(4*!dpi)
; add angle information and scattering angles
datadict['angle'] = data[keys[i], 'angle']
datadict['scatter'] = scatter

; read in the scattering cross section
qs1 = headers[keys[i], 'sca']
qs2 = headers[keys[i+2], 'sca']
; read in g parameter
g1 = headers[keys[i], 'g']
g2 = headers[keys[i+2], 'g']
; read in extinction cross section
ext1 = headers[keys[i], 'ext']
ext2 = headers[keys[i+2], 'ext']
; read in absorption cross section
abs1 = headers[keys[i], 'abs']
abs2 = headers[keys[i+2], 'abs']

; calculate a new phase function, a combination of pf1 and pf2
; weighted by the qs's
pf = addphase(pf1,pf2,qs1,qs2)
g = addphase(g1,g2,qs1,qs2)


; Identify the author
start = strpos(fnames[i],'aSil')
nfname = fname.Remove(start)
if nfname eq 'comp' then begin
pname = 'compiegne' 
endif
if nfname eq 'wd' then begin
; Identify the Rv value for weingartner and draine size distributions
numcheck = strpos(fnames[i], '31')
if numcheck ne -1 then begin
pname = 'wd31'
endif
if numcheck eq -1 then begin
pname = 'wd55'
endif
endif

; Find the color
ccheck = strpos(fnames[i],'green')
if ccheck eq -1 then begin
color = 'red'
pcolor = 'red'
endif
if ccheck ne -1 then begin
color = 'green'
pcolor = 'grn5'
endif

; create an output file name
pname += color

save, pf, filename = pname+'_pf.idlsav'
save, g, filename = pname+'_g.idlsav'

; use output file name as key for new phase function
datadict[pname] = pf
hgplot = dictionary(pname, g)

; construct title for plot
ptitle = 'Phase functions for a size distribution of carbon and silicate particles'

; produce information for legend
nmeanings = ['Mie', 'Henyey-Greenstein', 'Isotropic']
ncolors = [pcolor, pcolor, 'blue']
nlstyles = [0,2,0]
nthicks = [5,5,5]

; call plot options to make desired plot
plotoptions, datadict, outputdir, pname, ptitle, $
             nmeanings, ncolors, nlstyles, nthicks, $
             hgplot = hgplot, /isoplot, /nplot, /ylog,$
             yrange = [0.1,10], xrange = [0,180], /xstyle, $
             xtitle = 'scattering angle', $
             ytitle = 'phase function'

; Calculate g numerically, by integrating the phase function in
; spherical coordinates multiplied by cos(theta)
a = intpftog(datadict, pname)

; add numeric and analytic values of g to arrays for later plotting to
; confirm the values are the same
gpts1[k] = a
gpts2[k] = g

; Integrate the phase function in spherical coordinates to determine
; normalization convention
c = intpf(datadict, pname)
intg = 'Integral of renormalized phase function '+string(c)

newg = g
; Calculate the albedo of the combined phase functions
totalalb = (qs1 + qs2)/(ext1 + ext2)
; Calculate the ratio of scattering to absorption cross sections for
; the combined phase functions
totalqsqa = (qs1 + qs2)/(abs1 + abs2)

; Create strings to write to file
totalalb = 'Total albedo '+string(totalalb)
totalqsqa = 'Ratio of scattering to absorption cross section'+string(totalqsqa)
newg = 'Anisotropy parameter g '+string(newg)

; Write a parameter file with information about the combined phase functions
openw, lun, outputdir+pname+'.param',/get_lun
printf, lun, intg
printf, lun, totalalb
printf, lun, totalqsqa
printf, lun, newg
free_lun, lun

k+=1
i+=1
endwhile

; confirm g is the same whether calculated numerically or analytically
plot, gpts2, gpts1, psym = 6,  xtitle = 'analytic g',$
      ytitle = 'numeric g', yrange = [0.5,0.6], /isotropic

end
