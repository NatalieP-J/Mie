; Written May 13 2015 for SURP 2015 project

fname = 'zubko.refr'

; read in reference file of refraction indices
readcol, fname, wavelength, n, k

reqwv = [0.477,0.62]

; 1 in l1 means the file contains n, 0 means it is n-1
interpindices, fname, $
               '1 ! amorphous carbon, interpolated from Zubko 1996',$
               wavelength, n, k, reqwv
end
