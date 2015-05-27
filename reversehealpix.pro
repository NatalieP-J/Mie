; Written May 27 2015 for SURP 2015 project
; this function creates a mirror image of a healpix array around l=0

function reversehealpix, inarray, height
; inarray = input array to be mirrored
; height = height parameter of input array

; create an array of pixel indices (being sure to convert to integers)
ipix = long(findgen(leng(inarray)))
nside = findnside(height=height)

; find the Galactic coordinates of each pixel
pix2ang_ring, nside, ipix, pols, ells

; intialize output array of appropriate shape
outarray=inarray
; for each pixel
for i=0, leng(inarray)-1 do begin
; find the coordinates we wish to match
polmatch = pols[i]
ellmatch = -ells[i]
; calculate the pixel index for those coordinates
ang2pix_ring, nside, polmatch, ellmatch, matchloc
; set output array pixel equal to original value at that pixel index
outarray[i] = inarray[matchloc]
endfor

return, outarray

end
