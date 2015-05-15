; Written May 12 2015 for SURP 2015 project

function ctsky, ellillum = ellillum, beeillum = beeillum, nside = nside
; computes scattering angle at each pixel for an array of pixels
; ellillum, beeillum = galactic coordinates of incoming beam
; ellscat, beescat = galactic coordinates of scattering direction
;
; returns an array containing the cosine of the scattering angle at
; each pixel

; calculate the number of pixels in output map
npix = nside2npix(nside)
uscats = fltarr(npix)

for i = 0, npix-1 do begin
   pix2ang_ring, nside, i, polscat, ellscat
   beescat = !dpi/2. - polscat
   uscat = ct(ellillum=ellillum,beeillum=beeillum,$
              ellscat=ellscat,beescat=beescat)
   uscats[i] = uscat
endfor
return, uscats
end
