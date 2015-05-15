; Written May 13 2015 for SURP 2015 project

fname = zubko.refr

; read in reference file of refraction indices
readcol, fname, wavelength, n, k

interpindices, fname, $
               '1 ! amorphous carbon, interpolated from Zubko 1996',$
               wavelength, n, k, [0.477,0.62]
end
