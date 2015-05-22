ISRFDIR='/cita/d/raid-project/hp/pgmartin/GALPROP/FITS/ISRF/Standard/'


; sky with spectrum integrated over certain filter
; choice of 0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30
;filein = ISRFDIR+'Standard_8.5_0_0.5_Filter_HP.fits.gz'
filein = ISRFDIR+'Standard_8.5_0_1_Filter_HP.fits.gz'

foofilt = mrdfits(filein, 1, hfoof)
print, hfoof
; 28 energy bins
; Emin = 0.37
; log bin size 0.1957
; 28 energy bins
; actually, this is INCORRECT, the energies are in microns, and they are
; for astronomical filters, not logarithmically spaced
nside = sxpar(hfoof, 'NSIDE')
order = sxpar(hfoof,'ORDERING')
print, nside, order

foofilte = mrdfits(filein, 2, hfoofe)
print, hfoofe
; this is standard set of filters in microns
print, foofilte.energy

; whole sky in filter at 0.55 microns
ifilt = 2
skyfilt = dblarr(12288)
skyfilt[*] = (foofilt.spectra)[ifilt,*]

print, skyfilt
print, max(skyfilt), 'max skyfilt'
print, min(skyfilt), 'min skyfilt'

; display the full sky plot
; smart enough not to take the log of zero
mollview, skyfilt, /log, grat=[30.,30.], glsize= 1.

; now create a copies of the sky array, to be manipulated

; isoptropic unit intensity
skyiso = skyfilt
help, skyiso
skyiso[*] = 1.

; use the ISRF as illuminating input
skyinput = skyfilt
; for testing purposes use the isotropic as input
;skyinput = skyiso

; make an array to constain the scattered light on the sky
skyoutput = skyfilt

; number of pixels on sky at NSIDE = 32
npix = 12288
npix = nside2npix(32)

; in the all sky integrals, need to know the angular size of the
; illuminating pixel
domega = 4.*!pi/npix
;print, domega

mollview, skyinput, /log

sumsky = 0.
for i = 0,npix-1 do begin
   sumsky = sumsky + skyinput[i]
endfor

; mean and integrated
meansky = sumsky/npix
intsky = sumsky*domega

print,' mean over sky', meansky
; integral of isotropic unity input over whole sky would be 4pi
print, ' integral over sky', intsky, 4.*!pi*meansky

; halt here if just looking at average input properties
;return
;end


; estimate of the asymmetry parameter for this filter
gfilt = 0.5
gfilt = 0.65
uinput = 1.

; phase function tests

pf = hg(g=gfilt, u=uinput)

print, 'g, u, pf',  gfilt, uinput, pf

; get Galactic ell and bee coordinates correponding to the pixel positions
ipix = 0
ipix = 12288/2
ipix = 12288/2 - 1
pix2ang_ring, nside, ipix, pol, ell
; get Galactic b from the complementary polar angle
bee = !dpi/2. - pol
print, 'test l, b in degrees', ell*180./!dpi, bee*180./!dpi

; illumination from this direction in radians
ellillum = 0.
beeillum = 0.
beeillum = -10.*!dpi/180.
; scattering to this direction in radians
; back scattered relative to 0, 0
ellscat = 0.
beescat = 0.
;forward scattered relative to 0, 0
ellscat = !dpi
beescat = 0.
; sideways scattering relative to 0, 0
beescat =!dpi/2.
print, 'test angles in radians', ellillum, beeillum, ellscat, beescat

; u = cosine of scattering angle (ct for cos theta)
uscat = ct(ellillum=ellillum,beeillum=beeillum,ellscat=ellscat,beescat=beescat)
print, 'test u', uscat

pf = hg(g=gfilt, u=uscat)

print, 'test pf for g, u', pf, gfilt, uscat

; display the phase function on the sky for illumination from 0, 0
; also check the integral
sumpf = 0.

for i = 0, npix-1 do begin
; find Galactic angles for this pixel
pix2ang_ring, nside, i, polscat, ellscat
beescat = !dpi/2. - polscat
uscat = ct(ellillum=ellillum,beeillum=beeillum,ellscat=ellscat,beescat=beescat)
pf = hg(g=gfilt, u=uscat)
sumpf = sumpf + pf
skyoutput[i] = pf
endfor

;mollview, skyoutput, /log
;mollview, skyfilt, /log, grat=30., glsize= 1.
;mollview, skyfilt, /log, grat=[45.,15.], glsize= 1.
mollview, skyoutput, /log, grat=[30.,30.], glsize= 1.

; this integral over the whole sky should be 1
intpf = sumpf*domega
print, 'integral of the phase function over sky', intpf

; view from the pole, north on the left, Galactic centre at bottom
;orthview, skyfilt, rot=[0.,90.], /log 
; view from the pole, north on the left, Galactic centre at contact point
;orthview, skyfilt, rot=[0.,90.,90.], /log 
;orthview, skyfilt, rot=[0.,90.,90.], /log, grat=[45.,15.], glsize= 1.
orthview, skyoutput, rot=[0.,90.,90.], /log, grat=[45.,15.], glsize= 1.

print, 'beginning the double integral'

; integral over an illuminating sky defined by skyinput
; use the ISRF as illuminating input
skyinput = skyfilt
; for testing purposes use the isotropic as input
;skyinput = skyiso


; zero the array to constain the scattered light on the sky
skyoutput[*] = 0.

; look for isotropic output = <input> for isoptropic scattering -- yes
;gfilt = 0

for i = 0, npix-1 do begin
; find Galactic angles for scattering directed from this pixel
; test the inner loop
;i = 12288/2
pix2ang_ring, nside, i, polscat, ellscat
beescat = !dpi/2. - polscat

; sum over the illumination
sumscat = 0.
   for j = 0, npix-1 do begin
; use the ISRF as illuminating input
; find the Galactic angles for the illumination
   pix2ang_ring, nside, j, polillum, ellillum
   beeillum = !dpi/2. - polillum
   uscat = ct(ellillum=ellillum,beeillum=beeillum,ellscat=ellscat,beescat=beescat)
   pf = hg(g=gfilt, u=uscat)
   scattered = pf*skyinput[j]
   sumscat = sumscat + scattered
   endfor
skyoutput[i] = sumscat*domega
; print to watch inner loop complete
skyoutputtest = skyoutput[i]
print, 'the integrated scattered light at pixel i', i, skyoutputtest
endfor

; check on conservation over whole sky
suminput = 0.
sumoutput = 0.
for i = 0, npix-1 do begin
suminput = suminput + skyinput[i]
sumoutput = sumoutput + skyoutput[i]
endfor

meaninput = suminput/npix
intinput = suminput*domega
print, 'mean and integral of input over sky',meaninput, intinput

meanoutput = sumoutput/npix
intoutput = sumoutput*domega
;print, 'mean and integral of output over sky', meanoutput, intoutput, skyoutputtest/npix, skyoutputtest*domega
print, 'mean and integral of output over sky', meanoutput, intoutput
; note that if mean agrees, then so will the out integral, since they
; come from the same sum

mollview, skyoutput, /log, grat=[30.,30.], glsize= 1.

; save it
writefits, 'skyfiltVz0.2g0.65.fits', skyoutput, hfoof

; read it in again
footest = readfits('skyfiltVz1g0.65.fits', hfootest)
mollview, footest, /log, grat=[30.,30.], glsize= 1.
;print, footest

; get the value for the Draco cloud
; cloud at ellcloud, beecloud
; radiation coming from the cloud to the earth is going *approximately* in the opposite direction 
; ellcloudscat, beecloudstat
; example for Draco
print, 'example for Draco nebula'
ellcloud = 92.24 * !dpi/180.
beecloud = 38.43  * !dpi/180.
print, 'cloud coordinates', ellcloud * 180./!dpi, beecloud * 180./!dpi

ellcloudscat = ellcloud + !dpi
beecloudscat = -beecloud
print, 'direction cloud to Sun', ellcloudscat * 180./!dpi, beecloudscat * 180./!dpi

; colatitude
polcloudscat = !dpi/2 - beecloudscat
print, 'colatitude', polcloudscat * 180./!dpi
; find the pixel corresponding to this
icloudscatpix = 0
ang2pix_ring, nside, polcloudscat, ellcloudscat, icloudscatpix
; check: this will not be perfect because we found the containing pixel,
; which will not have quite the same centre coordinates as the chosen direction
pix2ang_ring, nside, icloudscatpix, theta, phi
print, ' colatitude', theta * 180./!dpi
print, 'direction cloud to Sun', phi * 180./!dpi - 180., 90. - theta * 180./!dpi

; recover the value from the map
;outputcloudscat = skyoutput[icloudscatpix]
outputcloudscat = footest[icloudscatpix]
;print, 'log of value scattered', alog10(outputcloudscat), 'ratio to mean', outputcloudscat/meanoutput
print, 'log of value scattered', alog10(outputcloudscat), 'ratio to mean', outputcloudscat/5.853e+8

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
ellcloud = 135. 
beecloud = 40.
; distance from Sun
ds = 100.


; LH2
print, 'example for LH2'
ellcloud = 152. 
beecloud = 53.
; distance from Sun
ds = 100.

; AG
print, 'example for AG'
ellcloud = 165. 
beecloud = 65.
; distance from Sun
ds = 100.

; Bootes
print, 'example for Bootes'
ellcloud = 58. 
beecloud = 68.
; distance from Sun
ds = 100.

; N1
print, 'example for N1'
ellcloud = 85. 
beecloud = 44.
; distance from Sun
ds = 100.


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

ellscat = ellcloudscat
beescat = beecloudscat


; illumination from this direction in radians
ellillum = 0.
beeillum = - atan(zh/rcgc)
print, 'elevation illumination angle ', beeillum*180./!dpi


uscat = ct(ellillum=ellillum,beeillum=beeillum,ellscat=ellscat,beescat=beescat)
theta = acos(uscat)
print, 'scattering angle', theta*180./!dpi


end
