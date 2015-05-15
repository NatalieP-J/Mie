function ct, ellillum=ellillum, beeillum=beeillum, ellscat=ellscat, beescat=beescat
; evaluate the scattering angle
; between the original direction and the scattered beam
; illumination from ellillum, beeillum is going in the 
; direction to ellillum+!pi,-beeillum
; direction of scattered beam is ellscat, beescat
ellto = ellillum+!dpi
beeto = -beeillum
cbm= cos(beeto - beescat)
cbp= cos(beeto + beescat)
return, cos(ellto - ellscat)*(cbm+cbp)/2. + (cbm-cbp)/2.
end
