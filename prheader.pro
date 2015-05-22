; this function written May 11 2015 for SURP project
function prheader, fname, hsize, startline = startline
; reads in headers
; fname = filename as string
; hsize = number of lines that compose header
; returns a dictionary with keywords corresponding to variable values
; in the header

setdefaultvalue, startline, 0

; open, read in and close the file
openr, lun, fname, /get_lun
header = strarr(hsize)
skip_lun, lun, startline, /lines
readf, lun, header
free_lun, lun

; replace equals signs with white space
header = repstr(header,'=', ' ')
; split each line on whitespace
header = strsplit(header,/extract)

; initialize dictionary to hold header information
hdict = dictionary()

; begin loop over each line in header
for i=0, leng(header)-2 do begin &$
   ; if line is not empty, proceed
   if emptyline(header[i]) ne 1 then begin &$
      ; initialize a string to hold the key
      strkey = '' &$
      ; note how many words the key is
      keylen = 0 &$
      ; begin loop over each element in the line
      for j=0, leng(header[i])-1 do begin &$
         ; define is-a-float to hold result of canfloat on elements
         line = header[i] &$
         isaf = isfloat(line[j]) &$
         ; if canfloat returns zero, element not a float - add it to key
         if isaf eq 0 then begin &$
            ; if the element is the second word in
            ; the key add a leading space
            if keylen gt 0 then begin &$
               strkey += ' '+line[j] &$
               keylen += 1 &$
               endif &$
            if keylen eq 0 then begin &$
               strkey += line[j] &$
               keylen += 1 &$
               endif &$
         endif &$
         ; if canfloat returns a nonzero float, the element is a float
         ; set the float in the dictionary with appropriate key
         if isaf ne 0 then begin &$
            ; convert all nonalphanumeric characters to '_'
            strkey = idl_validname(strkey,/convert_all) &$
            hdict[strkey] = isaf &$
            ; reset the key and its length 
            strkey = '' &$
            keylen = 0 &$
         endif &$
      endfor &$
   endif &$
endfor

return, hdict

end
