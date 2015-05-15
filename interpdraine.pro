; Written May 13 2015 for SURP 2015 project

fname = 'draine.refr'

; read in reference file of refraction indices
readcol, fname, wavelength, eps1, eps2, n, k

interpindices, fname,$
               '0 ! silicon, interpolated from draine.refr',$
               wavelength, n, k, [0.477,0.62]
end
