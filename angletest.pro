;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; scattering angle to Sun from perspective of cloud

; Draco
print, 'example for Draco nebula'
ellcloud = 92.24 
beecloud = 38.43
; distance from Sun
ds = 500.


; Spider
print, 'example for Spider'
;ellcloud = 135. 
;beecloud = 40.
; distance from Sun
;ds = 100.


; LH2
print, 'example for LH2'
;ellcloud = 152. 
;beecloud = 53.
; distance from Sun
;ds = 100.

; AG
print, 'example for AG'
;ellcloud = 165. 
;beecloud = 65.
; distance from Sun
;ds = 100.

; Bootes
print, 'example for Bootes'
;ellcloud = 58. 
;beecloud = 68.
; distance from Sun
;ds = 100.

; N1
print, 'example for N1'
;ellcloud = 85. 
;beecloud = 44.
; distance from Sun
;ds = 1000.


ellcloud = ellcloud * !dpi/180.
beecloud = beecloud  * !dpi/180.
print, 'cloud coordinates ', ellcloud * 180./!dpi, beecloud * 180./!dpi

; z height
zh = sin(beecloud)*ds
dp = cos(beecloud)*ds
print, 'z height ', zh, ' distance in plane ', dp

; distance to Galactic Centre
dgc = 8500.

; law of cosines
rcgc = sqrt(dgc*dgc + dp*dp - 2.*dgc*dp*cos(ellcloud))
print, 'distance of cloud to Galactic centre ', rcgc

; case of right angle triangle in plane
rcrit = sqrt(dgc*dgc - dp*dp)
ellcrit = asin(rcrit/dgc)
print, ' critical angle ', ellcrit*180./!dpi
; law of sines as a check
ellcloudscatcrit = !dpi + asin(dgc*sin(ellcrit)/rcrit)
print, 'critcal longitude cloud to Sun', ellcloudscatcrit * 180./!dpi


; law of sines
ellcloudscat = -asin(dgc*sin(ellcloud)/rcgc)
if ellcloudscat lt 0. then ellcloudscat = 2.*!dpi + ellcloudscat
;print, 'longitude cloud to Sun', ellcloudscat * 180./!dpi

if rcgc lt rcrit then ellcloudscat =  !dpi + asin(dgc*sin(ellcloud)/rcgc) 
print, 'longitude cloud to Sun', ellcloudscat * 180./!dpi

beecloudscat = -beecloud

ellscat = ellcloudscat
beescat = beecloudscat


; illumination from this direction in radians
ellillum = 0.
beeillum = - atan(zh/rcgc)
print, 'elevation illumination angle ', beeillum*180./!dpi


uscat = ct(ellillum=ellillum,beeillum=beeillum,ellscat=ellscat,beescat=beescat)
theta = acos(uscat)
print, 'scattering angle', theta*180./!dpi


bpixel = 0
lpixel = 0
ds = [500.,1000.]
a = pointingangles(ds, dgc, beecloud, ellcloud, bpixel, lpixel)

end
