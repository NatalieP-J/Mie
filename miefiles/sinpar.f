c sinpar.f
c consistency checks updated july 1999
c t1 updated mildly 19 oct 1992
c utility for mieuvx.f and mieuvxsd.f
      subroutine sinpar(x,qext,qphase,qscat,ctbrqs,qback,iflag)
      implicit real*8 (a-h,o-z)
      real*8 ni,nr
      complex*16 wn,wn1,wn2,sman,sman1,smbn,smbn1,nc,a
      complex*16 rrf,rrfx,tc1,tc2,cdum1,cdum2
      complex*16 nc2,nc212,ci
      complex*16 dcmplx
c 91 is .ge. jx 
c 181 is .ge.jxf
      parameter(n91=91,n181=181)
c must match main program; can be 1,2 in the cloudy application
      dimension taun(n91),taun1(n91),taun2(n91)
      parameter(nmxlim=16000)
      dimension pix1(n91),pix2(n91),pix3(n91),a(nmxlim)
      common/cmn/ p1f(n91),p2f(n91),p3f(n91),p4f(n91),p1b(n91),
     zp2b(n91 ),p3b(n91),p4b(n91),cstht(n91),si2tht(n91),thetd(n181),
     znr,ni,jx,nosini
      common/sinprc/wnr,wni,wn1r,wn1i,smanr,smani,sman1r,sman1i,smbnr,
     asmbni,smbn1r,smbn1i
c just assigns real part of a complex number
      equivalence (wn,wnr),(wn1,wn1r),(sman,smanr),(sman1,sman1r),
     a (smbn,smbnr),(smbn1,smbn1r)
      iflag=0
   20 nc=dcmplx(nr,-ni)
      rrf=1.d0/nc
      rx = 1.0d0/x
      rrfx = rrf * rx
c
c  t1 is the number of terms nmx2 that will be needed to obtain convergence
c  try to minimize this, because the a(n) downwards recursion has to
c  start at nmx1 larger than this
c
c major loop series is summed to nmx2, or less when converged
c nmx1 is used for a(n) only, n up to nmx2.
c must start evaluation sufficiently above nmx2 that a(nmx1)=(0.,0.)
c is a good approximation
c
c
corig with slight modification for extreme UV and X-ray, n near 1., large x
corig      t1=x*dmax1( 1.1d0,dsqrt(nr*nr+ni*ni) )*
corig     1(1.d0+0.02d0*dmax1(dexp(-x/100.d0)*x/10.d0,dlog10(x)))
c
c rules like those of wiscombe 1980 are slightly more efficient
      xrd=dexp(dlog(x)/3.d0)
c the final number in t1 was 1., 2. for large x, and 3. is needed sometimes
c see also idnint use below
      t1=x+4.d0*xrd+3.d0
      if(x.le.8.d0 .or. x.ge.4200.d0)  go to 9999 
c      t1=t1+0.05d0*xrd
c was 0., then 1., then 2., now 3. for intermediate x
c 19 oct 1992
      t1=t1+0.05d0*xrd+3.d0
 9999 continue
      t1=t1*1.01d0
c
c the original rule of dave for starting the downwards recursion was
c to start at 1.1*|mx| + 1, i.e. with the original version of t1
corig      nmx1=1.10d0*t1
c 
c try a simpler, less costly one, as in bohren and huffman, p 478
c this is the form for use with wiscombe rules for t1
c tests: it produces the same results as the more costly version
c
      nmx1=idnint( dmax1(t1,x*dsqrt(nr*nr+ni*ni)) ) +15
c
      if ( nmx1 .le. nmxlim-1 ) go to 21
      write(3,8) nmx1
 8    format(' failed: increase nmxlim, nmx1=',i10)
      iflag=2
      return
c normally call exit; here handle error flag externally
c
   21 nmx2=idnint( t1 )
corig      if ( nmx1  .gt. 150 ) go to 22
corig      nmx1 = 150
corig      nmx2 = 135
c
c try a more efficient scheme
      if(nmx2.gt.4) go to 22
      nmx2=4
      nmx1=idnint( dmax1(4.d0,x*dsqrt(nr*nr+ni*ni)) ) +15
c downwards recursion for logarithmic derivative
   22 a(nmx1+1)=(0.d0,0.d0)
c
c note that with the method of lentz 1976 (appl opt 15, 668), it would be
c possible to find a(nmx2) directly, and start the downwards recursion there
c however, there is not much in it with above form for nmx1 which uses just
      do 23 n=1,nmx1
         nn=nmx1-n+1
   23 a(nn)=(nn+1)*rrfx-1.d0/((nn+1)*rrfx+a(nn+1))
c
c initialize loop with calculations for n=1
c angle loop
      do 30 j=1,jx
         pix1(j)=0.d0
         pix2(j)=1.d0
         taun2(j)=0.d0
   30 taun1(j)=cstht(j)
c
      t1=dcos(x)
      t2=dsin(x)
      wn2=dcmplx(t1,-t2)
      wn1 = dcmplx(t2,t1)
      wn=rx*wn1-wn2
      tc1=a(1)*rrf+rx
      tc2=a(1)*nc+rx
      sman=(tc1*wnr-wn1r)/(tc1*wn-wn1)
      smbn=(tc2*wnr-wn1r)/(tc2*wn-wn1)
c
c small x; above calculations subject to rounding errors
c see bohren and huffman p 131 
c wiscombe 1980 appl opt 19, 1505 gives alternative formulation
      xcut=3.d-04
      if(x.ge.xcut) go to 1001
      ci=(0.d0,1.d0)
      nc2=nc*nc
      nc212=(nc2-1.d0)/(nc2+2.d0)
      x3=x*x*x
      x5=x3*x*x
c note change sign convention for m = n - ik here
      sman=ci*2.d0*x3*nc212*(1.d0/3.d0+x*x*0.2d0*(nc2-2.d0)/(nc2+2.d0))
     1   +  4.d0*x5*x*nc212*nc212/9.d0     
      smbn=ci*x5*(nc2-1.d0)/45.d0
c
 1001 sman1r=smanr
      sman1i=smani
      smbn1r=smbnr
      smbn1i=smbni
      t1=1.5d0
      smanr=t1 * smanr
      smani=t1 * smani
      smbnr=t1 * smbnr
      smbni=t1 * smbni
c angle loop
      do 60 j = 1,jx
         p=pix2(j)
         t=taun1(j)
         smanrp=smanr*p
         smanip=smani*p
         smbnrp=smbnr*p
         smbnip=smbni*p
         smbnrt=smbnr*t
         smbnit=smbni*t
         smanrt=smanr*t
         smanit=smani*t
         p1f(j)=smanrp+smbnrt
         p2f(j)=smanip+smbnit
         p3f(j)=smbnrp+smanrt
         p4f(j)=smbnip+smanit
         p1b(j)=smanrp-smbnrt
         p2b(j)=smanip-smbnit
         p3b(j)=smbnrp-smanrt
   60 p4b(j)=smbnip-smanit
c case n=1; note previous multiplication of sman and smbn by t1=1.5
      qext=2.d0*( smanr+smbnr)
      qphase=2.d0*( smani+smbni)
      nsqbk=-1
      qbackr=-2.d0*( smanr-smbnr)
      qbacki=-2.d0*( smani-smbni)
      qscat=(smanr*smanr+smani*smani+smbnr*smbnr+smbni*smbni)/.75d0
      ctbrqs = 0.0d0
      n=2
c************************ Major loop begins here ***********************
   65 t1=2*n-1
      t3=t1+2
      wn2=wn1
      wn1=wn
      wn=t1*rx*wn1-wn2
      cdum1=a(n)
      cdum2=n*rx
      tc1=cdum1*rrf+cdum2
      tc2=cdum1*nc+cdum2
      sman=(tc1*wnr-wn1r)/(tc1*wn-wn1)
      smbn=(tc2*wnr-wn1r)/(tc2*wn-wn1)
c
c small x, n=2
c see bohren and huffman p 131
      if(x.ge.xcut) go to 1002
      if(n.ne.2) go to 1002
c note change sign convention for m = n - ik here
      sman=ci*x5*(nc2-1.d0)/(15.d0*(2.d0*nc2+3.d0))
      smbn=(0.d0,0.d0)
c
 1002 eqext=t3*(smanr+smbnr)
      qext=qext+eqext
      eqpha=t3*(smani+smbni)
      qphase=qphase+eqpha
      nsqbk=-nsqbk
      eqbr=t3*(smanr-smbnr)*nsqbk
      qbackr=qbackr+eqbr 
      eqbi=t3*(smani-smbni)*nsqbk
      qbacki=qbacki+eqbi
      tx=smanr*smanr+smani*smani+smbnr*smbnr+smbni*smbni
      eqscat=t3*tx
      qscat=qscat+eqscat
      t2=n-1
c angle loop
      do 70 j = 1,jx
         t1pix2=t1*pix2(j)
         csthtj=cstht(j)
         pix1j=pix1(j)
         pix3(j)=(t1pix2*csthtj-n*pix1j)/t2
   70 taun(j)=csthtj*(pix3(j)-pix1j)-t1pix2*si2tht(j)+taun2(j)
c
      t5=n
      t4=t1/(t5*t2)
      t2=(t2*(t5+1.d0))/t5
      ectb=t2*(sman1r*smanr+sman1i*smani+smbn1r*smbnr+smbn1i*
     1   smbni)+t4*(sman1r*smbn1r+sman1i*smbn1i)
      ctbrqs=ctbrqs+ectb 
c angle loop
      t2=n*(n+1)
      t1=t3/t2
      k=(n/2)*2
      do 80 j = 1,jx
         p=t1*pix3(j)
         t=t1*taun(j)
         smanrp=smanr*p
         smanip=smani*p
         smbnrp=smbnr*p
         smbnip=smbni*p
         smanrt=smanr*t
         smanit=smani*t
         smbnrt=smbnr*t
         smbnit=smbni*t
c forward hemisphere
         p1f(j)=p1f(j)+smanrp+smbnrt
         p2f(j)=p2f(j)+smanip+smbnit
         p3f(j)=p3f(j)+smbnrp+smanrt
         p4f(j)=p4f(j)+smbnip+smanit
c backwards hemisphere
         if(k.eq.n)go to 75
c odd n
         p1b(j)=p1b(j)+smanrp-smbnrt
         p2b(j)=p2b(j)+smanip-smbnit
         p3b(j)=p3b(j)+smbnrp-smanrt
         p4b(j)=p4b(j)+smbnip-smanit
         go to 80
c even n
 75      p1b(j)=p1b(j)-smanrp+smbnrt
         p2b(j)=p2b(j)-smanip+smbnit
         p3b(j)=p3b(j)-smbnrp+smanrt
         p4b(j)=p4b(j)-smbnip+smanit
 80   continue
c
c check convergence
c could decrease for large x and small m-1 in UV - X-ray; probably negligible
      if(tx.ge.1.d-14) go to 101
c looks good but check relative convergence
      eqext=dabs(eqext/qext)
      eqpha=dabs(eqpha/qphase)
      eqscat=dabs(eqscat/qscat)
      eqbr=dabs(eqbr/qbackr)
      eqbi=dabs(eqbi/qbacki)
      ectb=dabs(ectb/ctbrqs)
      if(n.eq.2) ectb=0.d0
c leave out eqbr/i, which are sometimes least well converged
      error=dmax1(eqext,eqpha,eqscat,ectb)
c put a milder constraint on eqbr/i
      error1=dmax1(eqbr,eqbi)
      if(error.lt.1.d-07 .and. error1.lt.1.d-04) go to 100
c
c not sufficiently converged
c     
c      write(3,*) n
c      write(3,889) eqext,eqpha,eqscat,ectb
c      write(3,889) eqbr,eqbi,tx
c
c cut out after n=2 for small x, since approximation is being used
      if(x.ge.xcut) go to 101
      write(3,987) n
 987  format(' small x failed to converge for n =',i10)
      ifail=1
      go to 100
 101  continue
c
c shift for next use of recurrence relations
c
c angle loop
      do 90 j = 1,jx
         pix1(j)=pix2(j)
         pix2(j)=pix3(j)
         taun2(j)=taun1(j)
   90 taun1(j)=taun(j)
c
      smbn1=smbn
      sman1=sman
   66 n=n+1
      if(n.le.nmx2)go to 65
c********************** major loop ends here *******************
      write(3,888) nmx2
  888 format(' failed to converge for nmx2=',i10)
      iflag=1
c don't exit
c
c *******converged***************8
 100  continue
      if(nosini.eq.1) write(3,889) error,error1
 889  format(1x,' relative errors ',7(1x,1pd12.4))
c
c *** consistency checks for items calculated from scattering amplitudes****
c twice forward scattering amplitude, to check on qext, qphase
      if(dabs(si2tht(1)).gt.1.d-10) go to 678
      fsr=2.d0*p1f(1)
      fsi=2.d0*p2f(1)
      sumre=dabs(fsr/qext-1.d0)
      sumim=dabs(fsi/qphase-1.d0)
      sumrei=dmax1(sumre,sumim)
      if(sumrei.gt.1.d-14) 
     1    write(3,679) sumre,sumim,fsr,qext,fsi,qphase
 679  format(1x,'forwards consistency failure ',6(1x,1pd11.4))
c twice backward amplitude
      bsr=2.d0*p3b(1)
      bsi=2.d0*p4b(1)
      sumre=dabs(bsr/qbackr-1.d0)
      sumim=dabs(bsi/qbacki-1.d0)
      sumrei=dmax1(sumre,sumim)
      if(sumrei.gt.1.d-14) 
     1    write(3,680) sumre,sumim,bsr,qbackr,bsi,qbacki
 680  format(1x,'backwards consistency failure ',6(1x,1pd11.4))
 678  continue
c
c angle loop
c collect combinations of real and imaginary for scattering matrix elements
c note convention for 1 and 2
      do 120 j = 1, jx
         t1=p1f(j)
         t2=p2f(j)
         t3=p3f(j)
         t4=p4f(j)
         p1f(j)=t3*t3+t4*t4
         p2f(j)=t1*t1+t2*t2
         p3f(j)=t1*t3+t2*t4
         p4f(j)=t2*t3-t4*t1
c backward
         t1=p1b(j)
         t2=p2b(j)
         t3=p3b(j)
         t4=p4b(j)
         p1b(j)=t3*t3+t4*t4
         p2b(j)=t1*t1+t2*t2
         p3b(j)=t1*t3+t2*t4
         p4b(j)=t2*t3-t4*t1
 120  continue
c     
c      check sums: see bohren and huffman  p. 478
      do 121 j=1,jx
         t1=p1f(j)
         t2=p2f(j)
         t3=p3f(j)
         t4=p4f(j)
         t11=(t1+t2)/2.d0
         t12=(t1-t2)/2.d0
         sum=dabs((t12*t12+t3*t3+t4*t4)/(t11*t11)-1.d0)
         sum1=0.d0
         sum2=0.d0
         if(j.eq.1) sum1=dabs(t12/t11)
         if(j.eq.1) sum2=dabs(t4/t11)
         sumt=dmax1(sum,sum1,sum2)
         if(sumt.gt.1.d-12)
     z write(3,681) j,sum,sum1,sum2
 681     format(1x,'forwards scattering failure ',i3,3(1x,1pd12.4))
c backward
         t1=p1b(j)
         t2=p2b(j)
         t3=p3b(j)
         t4=p4b(j)
         t11=(t1+t2)/2.d0
         t12=(t1-t2)/2.d0
         sum=dabs((t12*t12+t3*t3+t4*t4)/(t11*t11)-1.d0)
         sum1=0.d0
         sum2=0.d0
         if(j.eq.1) sum1=dabs(t12/t11)
         if(j.eq.1) sum2=dabs(t4/t11)
         sumt=dmax1(sum,sum1,sum2)
c slacken this from 1.d-11
         if(sumt.gt.1.d-07)
     z    write(3,682) j,sum,sum1,sum2
 682     format(1x,'backwards scattering failure ',i3,3(1x,1pd12.4))
 121  continue
c
c renormalize
 131  t1=2.d0*rx*rx
      qext=qext*t1
      qphase=qphase*t1
      qback=(qbackr*qbackr+qbacki*qbacki)*rx*rx
      qscat=qscat*t1
      if(dabs(ni).gt.1.d-14) go to 132
c real refractive index check
      sum=dabs(qscat/qext-1.d0)
      if(sum.gt.1.d-13)
     z write(3,133) qext,qscat,sum
 133  format(1x,'real refractive index failure ',3(1x,1pd15.7))
 132  ctbrqs=2.d0*ctbrqs*t1
      if(nosini.eq.1) write(3,988) n,nmx2,nmx1
 988  format(1x,' n, nmx2, nmx1', 3i6)
      return
      end
