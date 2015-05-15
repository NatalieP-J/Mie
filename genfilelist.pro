function genfilelist, sizes=sizes, colors=colors, particles=particles

fnames = strarr(leng(colors), leng(particles), leng(sizes))

; construct file names and put them into array
for i=0, leng(sizes)-1 do begin
   for j=0, leng(particles)-1 do begin
      for k=0, leng(colors)-1 do begin
         fnames[k,j,i] = sizes[i]+colors[k]+particles[j]
      endfor
   endfor
endfor

fnames = reform(fnames, leng(fnames))

return, fnames

end
