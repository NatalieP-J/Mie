; Written May 12 2015 for SURP 2015 project

pro rgpcartnorm, fang, fpf1, fhgpf1, fpf2, fhgpf2, scatter, pname, ptitle
; Creates a cartesian plot of three phase functions vs angle
; fang = angle data
; fpf = phase function data
; fhgpf = Henyey - Greenstein phase function data
; scatter = any array of scattering directions
; pname = name of output ps file
; ptitle = plot title

colors = ['red', 'red','green','green', 'blue']
meanings = ['Mie', 'Henyey-Greenstein','Mie','Henyey-Greenstein', 'Isotropic']
linestyles = [0,2,0,2,0]

isovals = hg(g = 0, u = cos(fang))

fang *= !radeg
fpf1 /= isovals
fhgpf1 /= isovals
fpf2 /= isovals
fhgpf2 /= isovals

maxpf1 = (max(fhgpf1) > max(fpf1))
maxpf2 = (max(fhgpf2) > max(fpf2))

maxpf = (maxpf1 > maxpf2)

if maxpf eq maxpf1 then begin
if maxpf1 eq max(fpf1) then begin
plot, fang, fpf1
rpsopen, 'cartnormplot/'+pname+'_rgcartn.ps', /color, /landscape
plot, !x.crange, !y.crange, /nodata, title = ptitle,$
      xtitle = 'Scattering angle', ytitle = 'Normalized phase function'
oplot, fang, fpf1, color = cgcolor(colors[0]), linestyle = linestyles[0]
oplot, fang, fhgpf1, color = cgcolor(colors[1]), linestyle = linestyles[1]
oplot, fang, fpf2, color = cgcolor(colors[2]), linestyle = linestyles[2]
oplot, fang, fhgpf2, color = cgcolor(colors[3]), linestyle = linestyles[3]
endif 
if maxpf1 eq max(fhgpf1) then begin
plot, fang, fhgpf1
rpsopen, 'cartnormplot/'+pname+'_rgcartn.ps', /color, /landscape
plot, !x.crange, !y.crange, /nodata, title = ptitle,$
      xtitle = 'Scattering angle', ytitle = 'Normalized phase function'
oplot, fang, fhgpf1, color = cgcolor(colors[1]), linestyle = linestyles[1]
oplot, fang, fpf1, color = cgcolor(colors[0]), linestyle = linestyles[0]
oplot, fang, fhgpf2, color = cgcolor(colors[3]), linestyle = linestyles[3]
oplot, fang, fpf2, color = cgcolor(colors[2]), linestyle = linestyles[2]
endif
endif

if maxpf eq maxpf2 then begin
if maxpf2 eq max(fpf2) then begin
plot, fang, fpf2
rpsopen, 'cartnormplot/'+pname+'_rgcartn.ps', /color, /landscape
plot, !x.crange, !y.crange, /nodata, title = ptitle,$
      xtitle = 'Scattering angle', ytitle = 'Normalized phase function'
oplot, fang, fpf1, color = cgcolor(colors[0]), linestyle = linestyles[0]
oplot, fang, fhgpf1, color = cgcolor(colors[1]), linestyle = linestyles[1]
oplot, fang, fpf2, color = cgcolor(colors[2]), linestyle = linestyles[2]
oplot, fang, fhgpf2, color = cgcolor(colors[3]), linestyle = linestyles[3]
endif 
if maxpf2 eq max(fhgpf2) then begin
plot, fang, fhgpf2
rpsopen, 'cartnormplot/'+pname+'_rgcartn.ps', /color, /landscape
plot, !x.crange, !y.crange, /nodata, title = ptitle,$
      xtitle = 'Scattering angle', ytitle = 'Normalized phase function'
oplot, fang, fhgpf1, color = cgcolor(colors[1]), linestyle = linestyles[1]
oplot, fang, fpf1, color = cgcolor(colors[0]), linestyle = linestyles[0]
oplot, fang, fhgpf2, color = cgcolor(colors[3]), linestyle = linestyles[3]
oplot, fang, fpf2, color = cgcolor(colors[2]), linestyle = linestyles[2]
endif
endif
oplot, fang,isovals/isovals , color = cgcolor(colors[4])
legend, meanings, linestyle = linestyles, color=colors, /left
vline, scatter*!radeg
rpsclose, /high
cgfixps, 'cartnormplot/'+pname+'_rgcartn.ps'
cgps2pdf, 'cartnormplot/'+pname+'_rgcartn.ps'
spawn, 'rm cartnormplot/'+pname+'_rgcartn.ps'
end
