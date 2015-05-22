c size distributions for use with mieuvxsd
c
c
c check this code if anything other than power law is used
c
c?? gamma distribution temporarily removed july 99
c?? to avoid dependency on NAG library
c
      subroutine sizeds
      implicit real*8(a-h,o-z)
c the rlog(i) entered is assumed to be logarithmic
      dimension r(288)
c number of size distributions limited to 2 here
      common/sizds/ rm(2 ),rn(2 ),rp(2 ),
     z rin( 288),w( 288),sizdis(2 , 288),nsd(2 ),ngt,numsd,irlog
c
c irlog=1: the rin(i) entered is assumed to be logarithmic
      do 3002 i=1,ngt
         r(i)=rin(i)
 3002 continue
      if(irlog.eq.0) go to 3001
      do 3000 i=1,ngt
         r(i)=dexp(rin(i))
 3000 continue
 3001 continue
c
      do 1000 n=1,numsd
      nsdn=nsd(n)
   50 go to(10,11,12,13,14,15,16,17),nsdn
c****************** the following (nsd=1) is diermendjians cloud model *
   10 do 1 i=1,ngt
    1  sizdis(n,i)=((r(i)/rm(n))**6)*dexp(-6.d0*r(i)/rm(n))
      go to 1000
c****************** the following (nsd=2) is haze model m **************
   11 do 2 i=1,ngt
    2 sizdis(n,i)=5.33d0*1.d+4*r(i)*dexp(-8.944d0*dsqrt(r(i)))
      go to 1000
c****************** the following (nsd=3) is the gamma distribution ****
   12 continue
c??
      if(nsdn.eq.3) then
         write(6,*) ' abort: gamma distr not available'
         stop
      endif
      aa=1.d0/(rm(n)*rn(n)   )
      gm=(1.d0-2.d0*rn(n))/rn(n)
      do 3 i=1,ngt
c??      sizdis(n,i)=dexp(gm*dlog(aa)+(gm-1)*dlog(r(i))-aa*r(i)
c??     1	-s14abf(gm,ifail))
c nag library log of gamma function
3	continue
c??	if(ifail.gt.0) write(3,2000) ifail
2000	format(1x,'ifail=',i4,'in gamma distribution')
      go to 1000
c****************** the following (nsd=4) is the flat distribution *****
   13 do 4 i=1,ngt
      sizdis(n,i)=1.d0
    4 continue
      go to 1000
c**************the following (nsd=5) is the log-normal distribution ****
   14 sqrt2p=dsqrt(3.141592653589793d0)
      temp=sqrt2p*rn(n)
c????? r2 not defined      rp(n)=r2
      a2=dlog10(rm(n))
      a3=2.d0*rn(n)*rn(n)
      do 100 i=1,ngt
      a1=dlog10(r(i)/(1.d0-r(i)/rp(n)))
  100 sizdis(n,i)=1.d0/temp     *dexp(-(a1-a2)**2/a3)
      go to 1000
c********the following (nsd=6) is the generalized cloud model***********
15    continue
      b=0.5d0
      do 5 i=1,ngt
      aint=dlog(r(i)/rm(n))
5     sizdis(n,i)=dexp(rn(n)*aint)*  dexp(-b*dexp(rp(n)*aint))
      go to 1000
c ************** the following (nsd=7) is a power law *********************
   16 continue
      do 200 i=1,ngt
         aint=dlog(r(i)/rp(n))
         sizdis(n,i)=rm(n)*dexp(-rn(n)*aint)
 200  continue
      go to 1000
c **************************************************************
 17   continue
 1000 continue
      return
      end
