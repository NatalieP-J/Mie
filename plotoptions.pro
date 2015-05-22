pro plotoptions, datadict, outputdir, pname, ptitle,$
                 nmeanings, ncolors, nlstyles, nthicks,$
                 hgplot = hgplot, isoplot = isoplot,$
                 pplot = pplot, nplot=nplot,$
                 mirror = mirror, _extra=extra

; SET DEFAULT VALUES FOR KEYWORD PARAMETERS
setdefaultvalue, hgplot, dictionary('null',0)
setdefaultvalue, isoplot, 0, /boolean
setdefaultvalue, pplot, 0, /boolean
setdefaultvalue, nplot, 0, /boolean
setdefaultvalue, mirror, 0, /boolean


; ROUGH CHECK IF ANGLE IS IN DEGREES - IF SO, CONVERT TO RADIANS
angle = datadict.angle
angle = reform(angle, leng(angle))
if max(abs(angle)) gt !dpi+1 then begin
angle *= !dtor
endif

; EXTRACT KEYS FOR NON-ANGLE ENTRIES IN DATADICT
keys = datadict.keys()
keys = keys[where((keys ne 'angle' and keys ne 'scatter'))]

; IF MIRRORING OF DATA REQUIRED, FLIP DEPENDENT VARIABLE ARRAY AND
; FLIP SIGN AND ORDER OF ANGLE ARRAY
if mirror ne 0 then begin
angle = nmirrordata(angle)
for i=0, leng(keys)-1 do begin
temp = datadict[keys[i]]
temp = reform(temp, leng(temp))
datadict[keys[i]] = mirrordata(temp)
endfor
endif

; PRODUCE TEMPORARY VARIABLES TO CHECK IF HGPLOTTING REQUIRED
a = hgplot.values()
b = a.toarray()

; IF PLOTTING OF THE HENYEY-GREENSTEIN FUNCTION DESIRED, ADD 
; APPROPRIATE ARRAYS TO DATADICT
if total(b) ne 0 then begin
for i=0, leng(hgplot)-1 do begin
hgkey = string(keys[i])+'hg'
datadict[hgkey] = hg(g=hgplot[keys[i]],u=cos(angle))
endfor
endif

; EXTRACT KEYS FOR NON-ANGLE ENTRIES IN DATADICT
keys = datadict.keys()
keys = keys[where((keys ne 'angle' and keys ne 'scatter'))]

; PREPARE ARRAYS TO HOLD LEGEND AND PLOTTING INFORMATION
colors = dictionary()
meanings = dictionary()
lstyles = dictionary()

; CREATE VARIABLES TO HOLD INFORMATION ABOUT THE LARGEST VALUE OF 
; ALL NON-ANGLE ARRAYS IN DATADICT
mkey = keys[0]
mval = 0

for i=0, leng(keys)-1 do begin
pdata = datadict[keys[i]]
; FIND ARRAY CONTAINING THE HIGHEST VALUE OF ALL ARRAYS
if max(pdata) gt mval then begin
mval = max(pdata)
mkey = keys[i]
endif
; FIND THE COLOR OF THE WAVELENGTH (AND THUS COLOR OF THE LINE)
ccheck = strpos(keys[i],'green')
if ccheck eq -1 then begin
colors[keys[i]] = 'red'
endif
if ccheck ne -1 then begin
colors[keys[i]] = 'grn5'
endif
; DETERMINE THE TYPE OF SCATTERING - HG OR MIE - AND USE THE 
; APPROPRIATE LINESTYLE AND LABEL
tcheck = strpos(keys[i], 'hg')
if tcheck eq -1 then begin
meanings[keys[i]] = 'Mie'
lstyles[keys[i]] = 0
endif
if tcheck ne -1 then begin
meanings[keys[i]] = 'Henyey-Greenstein'
lstyles[keys[i]] = 2
endif
endfor

; IF ISOTROPIC FUNCTION TO BE PLOTTED, GENERATE IT AND ADD
; TO DATADICT, ADDING RELEVANT LABELS TO COLORS, MEANINGS AND
; LINESTYLES
if isoplot ne 0 then begin
isovals  = hg(g=0, u = cos(angle))
datadict['iso'] = hg(g=0, u = cos(angle))
colors['iso'] = 'blue'
meanings['iso'] = 'Isotropic'
lstyles['iso'] = 0
endif

; IF PLOT IS TO BE NORMALIZED, COMPUTE ISOTROPIC VALUES TO NORMALIZE TO
if nplot ne 0 then begin
isovals  = hg(g=0, u = cos(angle))
endif

; IF PLOT IS NOT TO BE NORMALIZED, GENERATE ISOVALS AS 1
if nplot eq 0 then begin
isovals = fltarr(leng(angle)) + 1
endif

; EXTRACT KEYS (FOR THE FINAL TIME)
keys = datadict.keys()
keys = keys[where((keys ne 'angle' and keys ne 'scatter'))]

; IF NOT MAKING A POLAR PLOT, USE THE FOLLOWING PLOTTING INSTRUCTIONS
if pplot eq 0 then begin
; CREATE OUTPUT POSTSCRIPT FILE
rpsopen, outputdir+pname+'.ps', /color, /landscape
; DRAW AXES IN BLACK AND SET APPROPRIATE BOUNDS ON X,Y RANGES
plot, angle*!radeg, datadict[mkey]/isovals,/nodata, title = ptitle,$
      _extra = extra
; FOR EACH DATA SET IN THE DICTIONARY, PLOT IT IN THE APPROPRIATE
; COLOR AND LINESTYLE
for i=0, leng(keys)-1 do begin
oplot, angle*!radeg, datadict[keys[i]]/isovals,$ 
       color = cgcolor(colors[keys[i]]), thick = 5,$
       linestyle = lstyles[keys[i]]
endfor
; ADD SCATTERING ANGLES POINTING TOWARD EARTH
vline, datadict['scatter']*!radeg
; ADD LEGEND
legend, nmeanings, linestyle = nlstyles, color = ncolors, thick = nthicks,$
        /right
rpsclose, /high
; ROTATE POSTSCRIPT FILE
cgfixps, outputdir+pname+'.ps'
; CONVERT POSTSCRIPT TO PDF AND REMOVE POSTSCRIPT FILE
cgps2pdf, outputdir+pname+'.ps'
spawn, 'rm '+outputdir+pname+'.ps'
endif

; IF MAKING A POLAR PLOT, USE THE FOLLOWING PLOTTING INSTRUCTIONS
if pplot ne 0 then begin
; CREATE OUTPUT POSTSCIPT FILE
rpsopen, outputdir+pname+'.ps', /color, /landscape
; DRAW AXES IN BLACK AND SET APPROPRIATE BOUNDS ON X,Y RANGES
plot, datadict[mkey]/isovals, angle, /polar,/nodata, title = ptitle,$
      _extra=extra
; FOR EACH DATA SET IN THE DICTIONARY, PLOT IT IN THE APPROPRIATE
; COLOUR AND LINESTYLE
for i=0, leng(keys)-1 do begin
oplot, datadict[keys[i]]/isovals, angle, color = cgcolor(colors[keys[i]]), $
       /polar, linestyle = lstyles[keys[i]], thick = 5
endfor
; ADD SCATTERING ANGLES POINTING TOWARD EARTH
angline, datadict['scatter']
; ADD AXIS LINES
hline, 0
vline, 0
; ADD LEGENDE
legend, nmeanings, linestyle = nlstyles, color = ncolors, thick = nthicks, $
        /right
rpsclose, /high
; ROTATE POSTSCRIPT FILE
cgfixps, outputdir+pname+'.ps'
; CONVERT POSTSCIPT TO PDF AND REMOVE POSTSCRIPT FILE
cgps2pdf, outputdir+pname+'.ps'
spawn, 'rm '+outputdir+pname+'.ps'
endif

end
