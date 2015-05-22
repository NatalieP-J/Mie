; Written May 15 2015 for SURP 2015

function sdgenfilelist, sds=sds, colors=colors, particles=particles
; sds = a list of possible size distributions
; colors = a list of colours of light
; particles = a list of particle types as determined by refraction
;             index used
; returns a list of file names

; create array to hold filenames
fnames = strarr(leng(colors), leng(particles), leng(sds))

; construct file names and put them into array
for i=0, leng(sds)-1 do begin &$
   ; ensure only appropriate size distributions and particles
   ; are paired (aSil <-> draine, Gra <-> zubko)
   typecheck = strpos(sds[i],'Sil') &$
   if typecheck ne -1 then begin &$
      k=0 &$
   endif &$
   if typecheck eq -1 then begin &$
      k=1 &$
   endif &$
   for j=0, leng(colors)-1 do begin &$
      fnames[j,k,i] = sds[i]+colors[j]+particles[k] &$
   endfor &$
   endfor

; transform array to be 1 dimensional
fnames = reform(fnames, leng(fnames))
nonempty = where(fnames ne '')
fnames = fnames[nonempty]

return, fnames

end
