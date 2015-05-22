; direction to Draco nebula from earth
dracoell = 92.24*!dtor
dracobee = 38.43*!dtor

; distance to the Draco nebula from earth
dracodist = [0.5,1.0] ; kpc

; distance to the galactic center from earth
GCdist = 8.5 ; kpc
; coordinates of galacit center in earth's frame
bGC = 0
lGC = 0

; produce scattering angles that correspond to a pointing from the
; nebula to earth
scatter = pointingangles(dracodist, GCdist, dracobee, dracoell, bGC, lGC)

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
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/RGcompare/'

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
datadict[keys[i+1]] = data[keys[i+1], 'pf']/(4*!dpi)
; add angle information and scattering angles
datadict['angle'] = data[keys[i], 'angle']
datadict['scatter'] = scatter


; add g values to indicate the desire to overplot the
; Henyey-Greenstein phase function
hgplot = dictionary()
hgplot[keys[i]] = headers[keys[i],'g']
hgplot[keys[i+1]] = headers[keys[i+1],'g']  

; find whether the particles are silicates or carbon
typecheck = strpos(fname,'zubko')
if typecheck ne -1 then begin
   ptype = 'carbon'
endif
if typecheck eq -1 then begin
   ptype = 'silicate'
endif
; find out which R_V value was used
rvcheck = strpos(fname, '31')
if rvcheck ne -1 then begin
   rvtype = textoidl('R_{V} = 3.1')
endif
if rvcheck eq -1 then begin
   rvtype = textoidl('R_{V} = 5.5')
endif

; construct a title for the plot and the savefile name
ptitle = 'Phase functions for a size distribution of '+ptype+$
         ' particles with '+rvtype   
start = strpos(fname, 'red')
nfname = fname.Remove(start,start+2)
pname = nfname+'_rgphase'

; produce information for the legend
nmeanings = ['Mie', 'Henyey-Greenstein','Mie', 'Henyey-Greenstein', 'Isotropic (normalized)']
ncolors = ['red', 'red', 'grn5', 'grn5', 'blue']
nlstyles = [0,2,0,2,0]
nthicks = [5,5,5,5,5]

; call plot options to make the desired plots
plotoptions, datadict, outputdir, pname, ptitle, $
             nmeanings, ncolors, nlstyles, nthicks, $
             hgplot = hgplot, /isoplot, /nplot, /ylog, $
             yrange = [0.1,10], xrange = [0,180], /xstyle, $
             xtitle = 'scattering angle', $
             ytitle = 'phase function';, charsize = 1.5, charthick = 3

; find g numerically
a = intpftog(datadict, keys[i])
b = intpftog(datadict, keys[i+1])

; add numeric and analytic values of g to their respective lists
gpts1[i] = a
gpts1[i+1] = b
gpts2[i] = headers[keys[i],'g']
gpts2[i+1] = headers[keys[i+1],'g']

; confirm that phase function integrates to 1
c = intpf(datadict, keys[i])
d = intpf(datadict, keys[i+1])

i+=2
endwhile

rpsopen, 'sdg.ps',/landscape
plot, gpts2, gpts1, psym = 6, yrange = [0.4,0.7], xtitle = 'analytic g',$
      ytitle = 'numeric g',/isotropic
rpsclose, /high
cgfixps, 'sdg.ps'
cgps2pdf, 'sdg.ps'
spawn, 'rm sdg.ps'
end
