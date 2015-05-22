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

; location of source files
sourcedir = '/mnt/raid-project/hp/njones/Mie/miefiles/'

; particle sizes
sizes = ['p001','p1','p2','p3', 'p4','p5','p6','p7']
;sizes = ['p15772', 'p16828001','p16954','p17422','p18618999','p18838','p19093999','p22543']
; wavelengths for which scattering is calculated
colors = ['red','green']
; refraction indices used
particles = ['rosil','zubko','draine']

; create a list of file names
fnames = genfilelist(sizes = sizes, colors = colors, particles = particles)

; specify file extension
exten = '.pr'
; specify size of header
hsize = 16
; include column titles
colnames = ['angle','m2','m1','s21','d21','pf','pol','xs','xd']

; read in data from each file in fnames
outdata = genfiledata(sourcedir, fnames, exten, hsize, colnames)

; split output into column data and header information
data = outdata['data']
headers = outdata['headers']

; set keys to use in data dictionary to feed to plot options
keys = fnames

; choose where to save plots
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/RGcompare/'

; initialize index
i=0
; create arrays to hold numeric (1) and analytic (2) values for g
gpts1 = fltarr(leng(keys))
gpts2 = fltarr(leng(keys))

while i lt leng(fnames) do begin
fname = fnames[i]

; initialize data dictionary to feed to plotoptions
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

; find the radius of the particles
radius = string(headers[keys[i], 'radius'],format = '(F0.3)')

; find whether the particles are silicates or carbon
typecheck = strpos(fname,'zubko')
if typecheck ne -1 then begin
   ptype = 'carbon'
endif
if typecheck eq -1 then begin
   ptype = 'silicate'
endif

; construct a title for the plot and the savefile name
ptitle = 'Phase functions for a '+ptype+$
         ' particle '+radius+' microns in size'   
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
             yrange = [0.01, 10], xrange = [0,180], /xstyle, $
             xtitle = 'scattering angle', $
             ytitle = 'phase function'

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

rpsopen, 'sinparg.ps',/landscape
plot, gpts2, gpts1, psym = 6,  xtitle = 'analytic g',$
      ytitle = 'numeric g',/isotropic
rpsclose, /high
cgfixps, 'sinparg.ps'
cgps2pdf, 'sinparg.ps'
spawn, 'rm sinparg.ps'
end
