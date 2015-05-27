; Written May 26 2015 for SURP 2015 project
; this code produces illumination maps for frankie inputs, as well as
; their mirror image about l=0 and the residuals of that mirror image
; with the original

; heights to plot illumination maps for; input as strings to be
; filename readable
heights = ['0.1','0.2','0.5','1'] ; in kpc

; for each height value
for i=0, leng(heights)-1 do begin
; import the illumination map
skyin = findskyin(height = heights[i])
; specify output directory
outdir = '/home/njones/Dropbox/Mie/NatalieResults/AllSky/frankie/'
; specify plot file name
pname = 'illumination_h'+heights[i]+'kpc.png'
; modify names for reversed (about l=0) and residual plots
rname = 'rev'+pname
resname = 'res'+pname
pname = outdir+pname
; Set title and subtitle
ptitle = 'Illumination Map'
stitle = '8.5 kpc from the Galactic Centre and '+heights[i]+' kpc above the midplane'
; produce plot of input radiation
mollview, skyin, grat = [30,30], glsize = 1., png = pname,$
          titleplot = ptitle, subtitle = stitle, /log

; Uncomment to make pdf output, mollview must be saved as postscript
;cgfixps, pname
;cgps2pdf, pname
;spawn, 'rm '+pname

; create mirror image map
reverseskyin = reversehealpix(skyin,heights[i])

; produce mirror image map
mollview, reverseskyin, grat = [30,30],$
          glsize = 1, titleplot = 'Reverse '+ptitle,$
          subtitle = stitle, /log, png = outdir+rname
          
; produce residual map
mollview, abs(reverseskyin-skyin), grat = [30,30], $
          glsize = 1, titleplot = 'Absolute '+ptitle+' Residuals',$
          subtitle = stitle, /log, png = outdir+resname

endfor
end
