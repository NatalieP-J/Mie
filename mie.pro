function mie, theta, angles=angles, pf=pf
; theta = angle at which the phase function is desired
; angles = angles at which the phase function is to be computed
; pf = phase function values at angles
; theta and angles may be degrees or radians but must have the same
; units

return, interpol(pf, angles, abs(theta), /spline)
end
