; Written May 19 2015 for SURP 2015 project
; reads .sisd files and plots weighted size distribution vs particle radius

; choose source directory and output directory
sourcedir = '/mnt/raid-project/hp/njones/Mie/miefiles/'
outputdir = '/home/njones/Dropbox/Mie/NatalieResults/SizeDistributions/'

; size distributions
sds = ['compaSil', 'compLamC', 'wdaSil31', 'wdGra31', 'wdaSil55', 'wdGra55']
; colours
colors = ['green']
; particle types
particles = ['draine', 'zubko']
; create list of filenames
fnames = sdgenfilelist(sds = sds, colors = colors, particles = particles)
; file extension
exten = '.sisd'
; number of columns total
ncols = 7
; number of header lines
hsize = 2

; for each of the filenames, perform the following
for i=0, leng(fnames)-1 do begin
fname = fnames[i]
; find author and particle type - remove colour, as size
; distribution doesn't change with colour
start = strpos(fname, 'green')
nfname = fname.Remove(start,start+4)
; read header information
sinfo = prheader(sourcedir+fname+exten, hsize)
; open file for reading columns
openr, lun, sourcedir+fname+exten, /get_lun
; skip past the header
skip_lun, lun, hsize, /lines
; find the number of elements in each array
len = sinfo['numsizes']
; calculate the number of lines associated with each array
nlines = len/ncols
; find whether the each column has the same number of entries
linetest = nlines mod long(nlines)

; if each column has the same number entries
if linetest eq 0 then begin
; prepare array to hold natural log of particle radius
   lnrad = fltarr(ncols, nlines)
; read the natural log of the particle radius
   readf, lun, lnrad
; reshape the array to have only 1 dimension
   lnrad = reform(lnrad, leng(lnrad))
; prepare array to hold particle radius
   rad = fltarr(ncols, nlines)
; skip title line
   skip_lun, lun, 1, /lines
; read particle radius
   readf, lun, rad
; reshape array to have only 1 dimension
   rad = reform(rad, leng(rad))
; prepare array to hold size distribution
   sized = fltarr(ncols, nlines)
; skip title line
   skip_lun, lun, 1, /lines
; read size distribution
   readf, lun, sized
; reshape array to have only 1 dimension
   sized = reform(sized, leng(sized))
; prepare array to hold weighted size distribution
   wsized = fltarr(ncols, nlines)
; skip title line
   skip_lun, lun, 1, /lines
; read weighted size distribution
   readf, lun, wsized
; reshape array to have only 1 dimension
   wsized = reform(wsized, leng(wsized))
; close file
   free_lun, lun
endif

; if columns don't all have the same number of entries
if linetest ne 0 then begin
; number of extra elements in the shorted line
   extras = len mod ncols
; prepare arrays to hold natural log of particle radius
   lnrad1 = fltarr(ncols, nlines)
   lnrad2 = fltarr(extras)
; read natural log of particle radius
   readf, lun, lnrad1
   readf, lun, lnrad2
; combine arrays and reshape to 1 dimension
   lnrad = [reform(lnrad1, leng(lnrad1)), reform(lnrad2, leng(lnrad2))]
; skip title line
   skip_lun, lun, 1, /lines 
; prepare array to hold particle radius
   rad1 = fltarr(ncols, nlines)
   rad2 = fltarr(extras)
; read particle radius
   readf, lun, rad1
   readf, lun, rad2
; combine arrays and reshape to 1 dimension
   rad = [reform(rad1, leng(rad1)), reform(rad2, leng(rad2))]
; skip header line
   skip_lun, lun, 1, /lines 
; prepare arrays to hold size distributions
   sized1 = fltarr(ncols, nlines)
   sized2 = fltarr(extras)
; read size distributions
   readf, lun, sized1
   readf, lun, sized2
; combine arrays and reshape to 1 dimension
   sized = [reform(sized1, leng(sized1)), reform(sized2, leng(sized2))]
; skip header line
   skip_lun, lun, 1, /lines 
; prepare arrays to hold weighted size distributions
   wsized1 = fltarr(ncols, nlines)
   wsized2 = fltarr(extras)
; read weighted size distributions
   readf, lun, wsized1
   readf, lun, wsized2
; combine arrays and reshape to 1 dimension
   wsized = [reform(wsized1, leng(wsized1)), reform(wsized2, leng(wsized2))]
; close file
   free_lun,lun
endif

; create Tex strings
mu = textoidl('\mu')
wsd = textoidl('10^{29} n_{H}^{-1} a^{4} dn/da (cm^3)')

; plot size distribution
rpsopen, outputdir+nfname+'.ps', /color, /landscape
plot, rad, wsized, /xlog, /ylog, xtitle = 'a ('+mu+'m)', ytitle = wsd, title = nfname
rpsclose, /high
cgfixps, outputdir+nfname+'.ps'
cgps2pdf, outputdir+nfname+'.ps'
spawn, 'rm '+outputdir+nfname+'.ps'


endfor

end
