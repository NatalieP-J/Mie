; Written May 13 2015 for SURP 2015 project

pro interpindices, fname, l1, wavelength, n, k, desiredwv
; interpolate refraction indices for desired wavelengths
; fname = fname of reference file
; l1 = first line to print in file
; wavelength = wavelengths for which indices were intially generated
; n = the real part of the refraction index at wavelength
; k = the imaginary part of the refraction index at wavelength
; desiredwv = array of wavelengths to find indices for (in microns)

; interpolate to find real and imaginary parts of the index of
; refraction at desired wavelengths
newn = interpol(n, wavelength, desiredwv ,/quadratic)
newk = interpol(k, wavelength, desiredwv,/quadratic)


; uncomment to see the success of the interpolation
;plot,[0.2,1.4],[0.0,3.0], color = 0, /nodata 
;oplot, wavelength, n, color = cgcolor('red')
;oplot, wavelength, k, color = cgcolor('blue')
;plots, desiredwv, newn, psym = 6, color = cgcolor('red')
;plots, desiredwv, newk, psym = 6, color = cgcolor('blue')

; concatenate arrays for the purposes of writing to file
wavelength = [wavelength,desiredwv]
n = [n,newn]
k = [k,newk]

; create output array
outarr = fltarr(3, leng(n))
outarr[0,*] = wavelength
outarr[1,*] = n
outarr[2,*] = k

; specify output directiory
outputdir = 'miefiles/'

; create a file contaning all data
openw, lun, outputdir+fname, /get_lun
; write header line
printf, lun, l1
printf, lun, outarr, format='(F0.7,1X,F0.7,1X,F0.7)'
close, lun

; create a file containing data about the green wavelength
openw, lun, outputdir+'green'+fname, /get_lun
; write header line
printf, lun, l1
printf, lun, outarr[*,-2], format='(F0.7,1X,F0.7,1X,F0.7)'
close, lun

; create a file containing data about the red wavelength
openw, lun, outputdir+'red'+fname, /get_lun
; write header line
printf, lun, l1
printf, lun, outarr[*,-1], format='(F0.7,1X,F0.7,1X,F0.7)'
close, lun 

end
