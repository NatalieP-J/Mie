
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

while i lt leng(fnames)-2 do begin
; skip segement so we are only combining data from the same authors
if ((i eq 2)||(i eq 6)) then begin
i+=2
endif

; initialize data dictionary to feed plotoptions
datadict = dictionary()
hgplot = dictionary()

for j=0, 1 do begin
i+=j
fname = fnames[i]

; find renormalized phase functions for both colors
pf1 = data[keys[i], 'pf']/(4*!dpi)
pf2 = data[keys[i+2], 'pf']/(4*!dpi)
; add angle information and scattering angles
datadict['angle'] = data[keys[i], 'angle']
datadict['scatter'] = scatter

; read in the scattering cross section and g parameters
qs1 = headers[keys[i], 'sca']
qs2 = headers[keys[i+2], 'sca']
g1 = headers[keys[i], 'g']
g2 = headers[keys[i+2], 'g']

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
numcheck = strpos(fnames[i],'31')
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
print, pname
; use output file name as key for new phase function
datadict[pname] = pf
hgplot[pname] = g
endfor

; Identify the author
start = strpos(fnames[i],'aSil')
nfname = fname.Remove(start)
if nfname eq 'comp' then begin
pname = 'compiegne_rg' 
endif
if nfname eq 'wd' then begin
; Identify the Rv value for weingartner and draine size distributions
numcheck = strpos(fnames[i],'31')
if numcheck ne -1 then begin
pname = 'wd31_rg'
endif
if numcheck eq -1 then begin
pname = 'wd55_rg'
endif
endif

; construct title for plot
ptitle = 'Phase functions for a size distribution of carbon and silicate particles'

; produce information for legend
nmeanings = ['Mie', 'Henyey-Greenstein','Mie', 'Henyey-Greenstein', 'Isotropic']
ncolors = ['red','red','grn5','grn5', 'blue']
nlstyles = [0,2,0,2,0]
nthicks = [5,5,5,5,5]

; call plot options to make desired plot
plotoptions, datadict, outputdir, pname, ptitle, $
             nmeanings, ncolors, nlstyles, nthicks, $
             hgplot = hgplot, /isoplot, /nplot, /ylog, $
             yrange = [0.1,10], xrange = [0,180], /xstyle, $
             xtitle = 'scattering angle', $
             ytitle = 'phase function'

i+=1
endwhile

end
