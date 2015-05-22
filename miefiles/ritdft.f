c utility for use with mieuvx/sd
      subroutine ritodf(nr,ni,eps1,eps2)
      implicit real*8 (a-h,o-z)
      real*8 ni,nr
c refractive index to dielectric function
      eps1 = nr*nr - ni*ni
      eps2 = 2.d0*nr*ni
      return
      end
c
      subroutine dftori(nr,ni,eps1,eps2)
      implicit real*8 (a-h,o-z)
      real*8 ni,nr
c dielectric function to refractive index 
      eps = dsqrt(eps2*eps2 + eps1*eps1)
      nr = dsqrt((eps + eps1)/2.d0)
      ni = dsqrt((eps - eps1)/2.d0)
      return
      end
c
