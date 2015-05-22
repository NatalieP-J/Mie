; Written May 15 2015 for SURP 2015 project

function genfiledata, sourcedir, fnames, exten, hsize, colnames
; sourcedir = source directory
; fnames = filenames to be read
; exten = file extension
; hsize = size of the header in each file. If header is split, this
;         array contains the number of lines in the first header as
;         the first entry and the startline for the remainder of the
;         header as the second.
; colnames = names of the data columns to be read
; returns a dictionary with two entries, headers and data, which in
; turn contain all header and column data respectively, both keyed to
; the appropriate filename

; intialize the two dictionaries
datas = dictionary()
headers = dictionary()

; cycle through each of the filenames
for i=0, leng(fnames)-1 do begin &$
   fname = fnames[i] &$
   ; create the column data dictionary for that file
   data = dictionary() &$
   if leng(hsize) eq 1 then begin &$
      ; if hsize is a scalar, use prheader to read in
      ; hsize lines as a dictionary
      hinfo = prheader(sourcedir+fname+exten, hsize) &$
   endif &$
   if leng(hsize) ne 1 then begin &$
      ; if hsize is a two by two array, then initialize a total
      ; header dictionary 
      hinfo = dictionary()  &$
      for j=0, leng(hsize)-1 do begin &$
         ; read in a header of approrpriate length starting from the
         ; appropriate line and add it to existing header information
         thinfo = prheader(sourcedir+fname+exten,hsize[j,0],$
                           startline=hsize[j,1]) &$
         hinfo = hinfo + thinfo &$
   endfor &$
   ; reset the header size to tell the header reader 
   ; how many lines to skip at the beginning of the file
   hsize = hsize[0,0] &$
   endif &$
   ; add header information to total header dictionary
   headers[fname] = hinfo &$
   ; create an array to read columns into (currently hardcoded
   ; to a length of 181 but this should be readable from the header
   ; as long as key naming is consistent
   output = fltarr(leng(colnames), 181) &$
   ; open the file, skip hsize lines and read the columns into output
   openr, lun, sourcedir+fname+exten, /get_lun &$
   skip_lun, lun, hsize, /lines &$
   readf, lun, output &$
   free_lun, lun &$
   ; add each column to the dictionary for this file according to the
   ; name given for it
   for k=0, leng(colnames)-1 do begin &$
      data[colnames[k]] = output[k,*] &$
   endfor &$
   ; add the data dictionary for this file to the total data dictionary
   datas[fname] = data &$
endfor

; create output dictionary and pack with headers and data info
outdata = dictionary()
outdata['data'] = datas
outdata['headers'] = headers

return, outdata

end
