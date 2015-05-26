heights = ['0.1','0.2','0.5','1'] ; in kpc

for i=0, leng(heights)-1 do begin
skyin = findskyin(height = heights[i])
outdir = '/home/njones/Dropbox/Mie/NatalieResults/AllSky/frankie/'
pname = 'illumination_h'+heights[i]+'kpc.ps'
pname = outdir+pname
ptitle = 'Illumination Map'
stitle = '8.5 kpc from and '+heights[i]+' kpc above the Galactic Center'
mollview, skyin, grat = [30,30], glsize = 1., ps = pname,$
          titleplot = ptitle, subtitle = stitle, /log

cgfixps, pname
cgps2pdf, pname
spawn, 'rm '+pname

endfor
end
