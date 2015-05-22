; Written May 21 2015 for SURP 2015 project

function addphase, pf1, pf2, qs1, qs2
; pf1 = array of phase function values
; pf2 = array of phase functions values
; qs1 = weight for pf1
; qs2 = weight for pf2
; returns the weighted sum of pf1 and pf2

; calculate numerator
numer = qs1*pf1 + qs2*pf2
; calculate denominator
denom = qs1+qs2

return, numer/denom

end
