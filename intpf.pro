; Written May 20 2015 for SURP 2015 project

function intpf, datadict, key
; datadict = a dictionary containing information about the phase
; function and corresponding scattering angle
; key = the key to access phase function in datadict
; returns the integral of the phase function in spherical coordinates

; access and reshape the phase function data
temp = datadict[key]
temp = reform(temp, leng(temp))
fpf = temp
; access and reshape the angle data - convert to radians
temp = datadict['angle']
temp = reform(temp, leng(temp))
fang = temp*!dtor

; produce array of scattering angles at which to calculate phase functions
thetas = make_array(18001,/index, start = 0, increment = 0.01)
thetas *= !dtor

; interpolate values of the phase function for a more accurate integral
nintg = interpol(fpf,fang,thetas,/spline)

; integrate the phase function in spherical coordinates
return, int_tabulated(thetas,2*!dpi*nintg*sin(thetas))
end
