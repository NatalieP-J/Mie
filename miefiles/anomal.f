c utility for use with mieuvx
      subroutine anomal(x,qext,qabs,qphase,xistar,delta,beta)
c
c in anomalous diffraction: qext and qabs calculated - qscatt by subtraction
c in rayleigh-gans:         qscatt and qabs calculated
c in mie:                   qext and qscatt calculated
c
      implicit real*8 (a-h,o-z)
      complex*16 cw,cbigk,ci
c anomalous diffraction: x>>1 and |m-1|<<1, any xi,xii
c original approach see Martin 1970. MN 149, 221
      xi=2.d0*x*delta
      xii=2.d0*x*beta
c xistar small is the basis for rayleigh-gans, any x, m-1
      xistar=dsqrt(xi*xi+xii*xii)
c alternative approach see martin 1978 p 23
      ci=(0.d0,1.d0)
      cw=-dcmplx(xi,xii)*ci
      call bigk(cw,cbigk)
      qext=4.d0*dreal(cbigk)
      qphase=4.d0*dimag(cbigk)
      cw=dcmplx(2.d0*xii,0.d0)
      call bigk(cw,cbigk)
      qabs=2.d0*dreal(cbigk)
c ?? put g in here - analytic version not known
      return
      end
c
      subroutine bigk(cw,cbigk)
c see martin 1978 p 23
      implicit real*8 (a-h,o-z)
      complex*16 cw,cbigk,cappr
c non-vax; use generic function
      d=abs(cw)
      if(d.lt.1.d-2)  go to 1
      cbigk=0.5d0+(cdexp(-cw)*(1.d0+cw)-1.d0)/(cw*cw)
      return
 1    continue
c avoid severe loss of precision for small cw; expand exponential
c coefficients are 1/n! - 1/(n+1)! = 1/(n+1)(n-1)! ;n=2,3,4,5,6,7
      cappr=cw*((1.d0/3.d0)
     1    - cw*((1.d0/8.d0)
     2    - cw*((1.d0/30.d0)
     3    - cw*((1.d0/144.d0)
     4    - cw*((1.d0/840.d0)
     5    - cw* (1.d0/5760.d0) ))))) 
c accurate to  (1+ order cw**6)
c      write(3,*) d,cbigk,cappr
      cbigk=cappr
      return
      end
c
