function genfiledata, sourcedir, fnames, exten, hsize, colnames

datas = dictionary()
headers = dictionary()

for i=0, leng(fnames)-1 do begin &$
   fname = fnames[i] &$
   data = dictionary() &$
   if leng(hsize) eq 1 then begin &$
      hinfo = prheader(sourcedir+fname+exten, hsize) &$
   endif &$
   if leng(hsize) ne 1 then begin &$
      hinfo = dictionary() &$
      for j=0, leng(hsize)-1 do begin &$
         thinfo = prheader(sourcedir+fname+exten,hsize[j,0],$
                           startline=hsize[j,1]) &$
         hinfo = hinfo + thinfo &$
   endfor &$
   hsize = hsize[0,0] &$
   endif &$
   headers[fname] = hinfo &$
   output = fltarr(leng(colnames), 181) &$
   openr, lun, sourcedir+fname+exten, /get_lun &$
   skip_lun, lun, hsize, /lines &$
   readf, lun, output &$
   free_lun, lun &$
   for k=0, leng(colnames)-1 do begin &$
      data[colnames[k]] = output[k,*] &$
   endfor &$
   datas[fname] = data &$
endfor

outdata = dictionary()
outdata['data'] = datas
outdata['headers'] = headers

return, outdata

end
