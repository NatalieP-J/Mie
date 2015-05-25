function pointingangles, dcloud, dpixel, bcloud, lcloud, bpixel, lpixel
; dcloud = distance from earth to cloud
; dpixel = distance from earth to pixel of interest
; bcloud = galactic latitude from earth to cloud 
; lcloud = galactic longitude from earth to cloud (0 longitude is
;    direction to galactic center)
; bpixel = galactic latitude from earth to illuminating pixel
; lpixel = galactic longitude from earth to illuminating pixel
; returns the scattering angle toward earth corresponding to the above parameter

; CURRENTLY ASSUMES billum = lillum = 0

; z height of cloud from earth
zh = dcloud*sin(bcloud)
; distance to cloud from earth in the plane
dp = dcloud*cos(bcloud)

;print, 'z height', zh, ' distance in plane ', dp

; find critical length (case of right angle triangle in the plane)
dcrit = sqrt(dpixel^2 - dp^2)

; calculate the distance from cloud to pixel in the plane
if lcloud lt !dpi then begin
dillum = sqrt(dpixel^2 + dp^2 - 2*dpixel*dp*cos(lcloud))
endif

if lcloud gt !dpi then begin
dillum = sqrt(dpixel^2 + dp^2 + 2*dpixel*dp*cos(lcloud))
endif

;print, 'distance of cloud to Galactic centre ', dillum

; calculates the angle from the cloud to the illuminating pixel
lillum = 0
billum = -atan(zh/dillum)

; calculate the angle from the cloud to the earth
learth = -asin(dpixel*sin(lcloud)/dillum)
crits = where(dillum lt dcrit)
if total(crits) ne -1 then learth[crits] = !dpi - learth[crits]
lt0 = where(learth lt 0)
if total(lt0) ne -1 then learth[lt0] += 2*!dpi
;print, learth*!radeg
bearth = fltarr(leng(learth))-bcloud

;print, 'longitude cloud to sun ', learth*!radeg

;print, 'elevation illumination angle ', billum*!radeg

; find the scattering angle
uscat = ct(ellillum = lillum, beeillum = billum, ellscat = learth, beescat = bearth)
scatter = acos(uscat)
;print, 'scattering angle ', scatter*!radeg
outdict = dictionary()
outdict['scatter'] = scatter
outdict['learth'] = learth
outdict['bearth'] = bearth

return, outdict
end
