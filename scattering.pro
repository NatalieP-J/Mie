; testing the ISRF Healpix files
; note the bintable environments

; Coordinate of interest:
GLONcoord = 92.24 ; degrees
GLATcoord = 38.43 ; degrees

D = 0. ; kpc
Rsun = 8.5 ; kpc
DtoRratio = D/Rsun

; Anticoordinate (exact solution)
;  Could also assume DtoRratio is 0
GLONanti = 2.*!pi - asin(sqrt( sin(GLONcoord/!radeg)^2 / (1 + DtoRratio^2 - 2*DtoRratio*cos(GLONcoord/!radeg))))
if GLONanti lt 0 then GLONanti = GLONanti + 2*!pi
if GLONanti gt 2.*!pi then GLONanti = GLONanti - 2.*!pi

GLATanti = !pi/2. - (-GLATcoord)/!radeg ; as theta

isotropic  = 0.

Rdistfoo = ['0','2','4','6','8.5','10','12','16','20','30']
Rdist = '8.5' ; distance from galactic centre

Zdistfoo = ['0','0.1','0.2','0.5','1','2','5','10','20','30']
Zdist = '1' ; distance above plane in kpc
Zdist = '0.5'
Zdist = '1'

g = 0.65 ; Henyey-Greenstein scattering parameter
g = 0.5
g=0

gfoo = [0,0.2,0.5,0.65]

scatlight = dblarr(n_elements(zdistfoo), n_elements(gfoo))
scatfrac = scatlight

for Zind = 0, n_elements(Zdistfoo)-1 do begin

   Zdist = Zdistfoo[Zind]

   FILT = 0.55d

   ISRFDIR='/cita/d/raid-project/hp/pgmartin/GALPROP/FITS/ISRF/Standard/'

   selection = 'Direct'
   selection = 'Scattered'
   selection = 'Thermal'
   selection = 'Total'
   selection = 'Filter'

; Healpix files
   filein = ISRFDIR+'Standard_'+Rdist+'_0_'+Zdist+'_'+selection+'_HP.fits.gz'

   foo1 = mrdfits(filein, 1, hfoo1)
   ordering = sxpar(hfoo1, 'ORDERING')
   nside = sxpar(hfoo1, 'NSIDE')
   foo2 = mrdfits(filein, 2, hfoo2)

   wavel = foo2.energy          ; in microns

   indfilt = where(wavel eq FILT)

   filtmap = foo1[*].spectra[indfilt[0]]

mollview, filtmap, max=1e9

   pixvec = indgen(n_elements(filtmap))

   if strtrim(ordering,2) eq 'NEST' then begin
      pix2ang_nest, nside, pixvec, theta,phi
      ang2pix_nest, nside, GLATanti, GLONanti, antipix
   endif
   if strtrim(ordering,2) eq 'RING' then begin
      pix2ang_ring, nside, pixvec, theta,phi
      ang2pix_ring, nside, GLATanti, GLONanti, antipix
   endif

   glat = !pi/2.-theta
   glon = phi
   ind = where(phi gt !pi and phi lt 2.*!pi)
   glon[ind] = glon[ind]-2.*!pi

   for gind = 0, n_elements(gfoo)-1 do begin

      g = gfoo[gind]

; Calculate scattered light for each pixel
      scatmap = dblarr(n_elements(glon))
      for i=0, n_elements(glon)-1 do begin
; Incident radiation
         radin = filtmap[i]
         if keyword_set(isotropic) then radin = 1.0
; Location of incident radiation
         glatin = glat[i]
         glonin = glon[i]
; theta is measured relative to the anti-direction, so distances are
; adjusted accordingly
         theta = !pi - sphdist(glonin, glatin, glon, glat)
         hg =  (1. - g^2)/(1.0+g^2 - 2.0*g*cos(theta))^1.5
; Normalize
         hg = hg / total(hg)
         scatmap = scatmap + hg*radin
      endfor

mollview, scatmap, /graticule, /glsize,/log

      print, GLONcoord, GLATcoord, scatmap[antipix], scatmap[antipix]/mean(scatmap)

      scatlight[zind, gind] = scatmap[antipix]
      scatfrac[zind, gind] = scatmap[antipix]/mean(scatmap)
   endfor
endfor


set_plot, 'ps'
device, file = 'draco_scatlightmodel.eps', /encapsulated, /color
plot, zdistfoo, scatlight[*,0], /ylog, /xlog, xrange=[0.1,30], xst=1, $
      yrange=[5e6,8e8], yst=1, xtitle = 'z (kpc)', ytitle = 'Scattered light'
oplot, zdistfoo, scatlight[*,2], thick=4
oplot, zdistfoo, scatlight[*,1], thick=2
oplot, zdistfoo, scatlight[*,3], thick=6
xyouts, 10, 5e8, 'R = 8.5 kpc'
xyouts, 1, 10.^7.0, 'g = 0.0'
xyouts, 1, 10.^7.07, 'g = 0.2'
xyouts, 1, 10.^7.14, 'g = 0.5'
xyouts, 1, 10.^7.21, 'g = 0.65'
device,/close
set_plot, 'x'


set_plot, 'ps'
device, file = 'draco_scatlightmodel2.eps', /encapsulated, /color
plot, zdistfoo, scatfrac[*,0], /ylog, /xlog, xrange=[0.1,30], xst=1, $
      yrange=[0.1,1.0], yst=1, xtitle = 'z (kpc)', ytitle = 'Scattered light/Mean'
oplot, zdistfoo, scatfrac[*,2], thick=4
oplot, zdistfoo, scatfrac[*,1], thick=2
oplot, zdistfoo, scatfrac[*,3], thick=6
xyouts, 10, 5e8, 'R = 8.5 kpc'
xyouts, 0.2, 10.^(-0.4), 'g = 0.0'
xyouts, 0.2, 10.^(-0.43), 'g = 0.2'
xyouts, 0.2, 10.^(-0.46), 'g = 0.5'
xyouts, 0.2, 10.^(-0.49), 'g = 0.65'
device,/close
set_plot, 'x'

end
