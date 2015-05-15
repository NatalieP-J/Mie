
function hg, g=g, u=u
; g is the asymmetry parameter
; u is the cosine of the scattering angle
; normalized such that 1/2 int hg du = 1.
; if used as 2 dimensional, a further 2 pi (then it is per steradian)
return, ( (1. - g^2)/(1. + g^2 - 2.*g*u)^1.5 )/(4*!pi)
end
