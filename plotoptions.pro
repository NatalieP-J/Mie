pro plotoptions, datadict, outputdir, pname, ptitle,$
                 hgplot = hgplot, isoplot = isoplot,$
                 pplot = pplot, nplot=nplot,$
                 mirror = mirror, _extra=extra

setdefaultvalue, hgplot, 0
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

if mirror ne 0 then begin
angle = nmirrordata(angle)
for i=0, leng(keys)-1 do begin
temp = datadict[keys[i]]
temp = reform(temp, leng(temp))
datadict[keys[i]] = mirrordata(temp)
endfor
endif



; IF PLOTTING OF THE HENYEY-GREENSTEIN FUNCTION DESIRED, ADD 
; APPROPRIATE ARRAYS TO DATADICT
if total(hgplot) ne 0 then begin
for i=0, leng(hgplot)-1 do begin
hgkey = string(keys[i])+'hg'
datadict[hgkey] = hg(g=hgplot[i],u=cos(angle))
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
colors[keys[i]] = 'green'
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
rpsopen, outputdir+pname+'.ps', /color, /landscape
plot, angle*!radeg, datadict[mkey]/isovals,/nodata, title = ptitle,$
      _extra = extra
for i=0, leng(keys)-1 do begin
oplot, angle*!radeg, datadict[keys[i]]/isovals, color = cgcolor(colors[keys[i]]), $
       linestyle = lstyles[keys[i]]
endfor
vline, datadict['scatter']*!radeg
meanings = meanings.values()
meanings = meanings.toarray()
lstyles = lstyles.values()
lstyles = lstyles.toarray()
colors = colors.values()
colors = colors.toarray()
legend, meanings, linestyle = lstyles, color = colors
rpsclose, /high
cgfixps, outputdir+pname+'.ps'
cgps2pdf, outputdir+pname+'.ps'
spawn, 'rm '+outputdir+pname+'.ps'
endif

if pplot ne 0 then begin
rpsopen, outputdir+pname+'.ps', /color, /landscape
plot, datadict[mkey]/isovals, angle, /polar,/nodata, title = ptitle,$
      _extra=extra
for i=0, leng(keys)-1 do begin
oplot, datadict[keys[i]]/isovals, angle, color = cgcolor(colors[keys[i]]), $
       /polar, linestyle = lstyles[keys[i]]
endfor
angline, datadict['scatter']
hline, 0
vline, 0
meanings = meanings.values()
meanings = meanings.toarray()
lstyles = lstyles.values()
lstyles = lstyles.toarray()
colors = colors.values()
colors = colors.toarray()
legend, meanings, linestyle = lstyles, color = colors
rpsclose, /high
cgfixps, outputdir+pname+'.ps'
cgps2pdf, outputdir+pname+'.ps'
spawn, 'rm '+outputdir+pname+'.ps'
endif

end
