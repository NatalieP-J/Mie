c mieuvxsd - like mieuvx, with integration over size distribution
c i/o updated for linux 2014
c max1 changed to amax1 july 99
c can use logarithmic increments in size, so as to cover many decades
c oct 1988 for uv - x-ray extinction, 
c excluding anomalous diffraction check that was in mieuvx
c     this version reads in real and imaginary parts of the refractive
c     index, with imaginary part positive (nridf= 1) or nr-1 (nridf=0) or
c     real and imaginary parts of the dielectric function (nridf = -1)
c dec 1988: added qback (removed in ...sd)
c approximation for small x in sinpar; 
c qphase, better convergence checking
c
      implicit real*8 (a-h,o-z)
      real*8 ni,nr,nr1,numpar
      logical lopen
c 96 is .ge. ng (even)
c 288 is .ge. ngt
      dimension xx(96),a(96),rr(96),ww(96),qscat(288)
c 91 is .ge. jx 
c 181 is .ge.jxf
      common/cmn/ p1f(91 ),p2f(91 ),p3f(91 ),p4f(91 ),p1b(91 ),p2b(91 ),
     zp3b(91 ),p4b(91 ),cstht(91 ),si2tht(91 ),thetd(181 )
     z,nr,ni,jx,nosini
c
c number of size distributions limited to 2 here
c 2 is .ge. numsd
      common/sizds/ rm(2 ),rn(2 ),rp(2 ),
     z r( 288),w( 288),sizdis(2 , 288),nsd(2 ),ngt,numsd,irlog
c
      common/sizint/ numpar(2 ),sigsca(2 )
     z,sigext(2 ),pizero(2 ),cosbar(2 ),rbarex(2 ),rbarsc(2 ),rbar10(2 )
     z,rbar32(2 ),rbar76(2 ),rscmom(2 ),rthmom(2 ),rthmm(2 ),rsqeff(2 ),
     z rbarte(2 ),rbartf(2 ),sigpha(2 ),sigabs(2 ),sigpr(2 ),vol(2 )
     z,cover(2 ),p1(2 ,91),p2(2 ,91),p3(2 ,91),p4(2 ,91)
c
c simple association of i/o logical units with global names, in fortran
c link with -li77 explicitly first
c see ioinit(3f); see fortran library routines, page 173
c sample shell usage to attach to unit 1:
c setenv fort01 iofilename in c-shell
c or fort01=iofilename; export fort01 in bourne shell
c
c iris
c
c change shell script to get around this f77 construction
c      call ioinit(.true.,.false.,.false.,'fort',.false.)
c
c May 2015, standard i/o to go with new shell script mieuvx.sh
c temporary file names: "t" for temporary
c INPUT
c 2 input wavelength and refractive index
      open(unit=2,file='tsd2.refr',status='old')
c 4 input mie, rho, size distribution
      open(unit=4,file='tsd4.uvxsd',status='old')
c OUTPUT
c 3 printed output
      open(unit=3,file='tsd3.prsd',status='old')
c 1 single particle output for verification
      open(unit=1,file='tsd1.sd',status='old')
c
c
c*****punch is positive for punching
      zz=0.d0
      punch=-1.
 700  format(1x)
      pi=3.141592653589793d0
c	refractive index of surrounding medium (real)
      refmed = 1.d0
      twopi=2.d0*pi
      read(4,*) nosina,nosini,nopol,nofin,norpr
      write(3,8320) nosina,nosini,nopol,nofin,norpr
 8319 format(i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x)
 8320 format(' nosina=',i2,' nosini=',i2,' nopol=',i2,
     1    ' nofin=',i2,' norpr=',i2)
c*********************compute scattering angles*************************
c jx=0, gives jx=1, jxf=2, angles 0, 180 
c jx=1, gives jx=2, jxf=3, angles 0,90, 90,180, with 90 actually computed twice
c in sinpar (set up for general case which might not hit 90)
      read(4,*) jx
      dth=0.d0
      if(jx.gt.0) dth=90.d0/dfloat(jx)
      jx=jx+1
c      jxf = max1(2 , 2*jx-1)
      jxf = max(2 , 2*jx-1)
c??dec      jxf = imax(2 , 2*jx-1)
      do 216 j=1,jx
      thetd(j)=dfloat(j-1)*dth
      cstht(j)=dcos(thetd(j)*pi/180.d0)
      c=cstht(j)
      si2tht(j)=1.d0-c*c
      l=jxf-(j-1)
 216  thetd(l)=180.d0-thetd(j)
      write(3,698) jx,jxf
 698  format(' number of scattering angles = ',i3,1x,i3)
c
      read(4,*) rho
      write(3,2003) rho
 2003 format(' density =', f12.6)
c
c set up the radii
c read the lower and upper sizes in microns, same for each size distribution
      read(4,*) irlog,r1,r2
c irlog=1, use logarithmic increments in size, so as to cover many decades
      if(irlog.eq.0) go to 6123
      r1=dlog(r1)
      r2=dlog(r2)
 6123 continue
c read in the basic gaussian quadrature order, and the number of partitions
      read(4,*) ng,ipart
      write(3,6122)ng,ipart
 6122 format(' ng=',i2,' ipart=',i3)
      eps=1.d-16
      call  gausle(ng,eps,csa,tsa,xx,a)
      ngt=ng*ipart
      rnge=r2-r1
      div=rnge/ipart
      do 2220 i=1,ipart
      xbot=(i-1)*div+r1
      xtop=xbot+div
      lower=(i-1)*ng
      call  gaussi(ng,xbot,xtop,xx,a,rr,ww)
      do 2220 kk=1,ng
      kkk=lower+kk
      r(kkk)=rr(kk)
 2220 w(kkk)=ww(kk)
      if(irlog.eq.0) go to 2332
c add the factor r in the integrand, from converting to logarithmic integral
      do 2331 kk=1,ngt
      w(kk)=w(kk)*dexp(r(kk))
 2331 continue
 2332 continue
      if(irlog.ne.0) write(3,*) ' natural log radii'
      if(irlog.eq.0) write(3,*) '  radii'
      write(3,4444) (r(k),k=1,ngt)
 4444 format(1x,7(1x,1pd10.3))
c
c*******choose the number of different size distributions and their para
c******the wavelength and upper and lower sizes will be the same for eac
      read(4,*) numsd
      write(3,4726) numsd
 4726 format(' numsd='i4)
      do 830 n=1,numsd
      read(4,*) nsd(n),rm(n),rn(n),rp(n)
c      write(3,4729) nsd(n),rm(n),rn(n),rp(n)
c 4729 format(' nsd=',i3,' rm=',1pd12.5,' rn=',1pd12.5,' rp=',1pd12.5)
 830  continue
c
      call sizeds
c integrate the factors which do not depend on the cross-sections
c initialize
      do 8833 n=1,numsd
         numpar(n)=0.d0
         rbar10(n)=0.d0
         rbar32(n)=0.d0
         rbar76(n)=0.d0
         rbarte(n)=0.d0
         rbartf(n)=0.d0
 8833 continue
c integrate over size
      do 8830 jco=1,ngt
         rrr=r(jco)
         if(irlog.ne.0) rrr=dexp(rrr)
         rrr3=rrr*rrr*rrr
         do 8829 n=1,numsd
            sw=sizdis(n,jco)*w(jco)
c number of particles
            numpar(n)=numpar(n)+sw
c for size
            rbar10(n)=rbar10(n)+rrr*sw
c for volume
            rbar32(n)=rbar32(n)+rrr3*sw
c for area
            rbarte(n)=rbarte(n)+rrr*rrr*sw
c r**6
            rbartf(n)=rbartf(n)+rrr3*rrr3*sw
c r**7
            rbar76(n)=rbar76(n)+rrr3*rrr3*rrr*sw
 8829    continue
 8830 continue
c renormalize
      do 8831 n=1,numsd
         vol(n)=4.d0*pi*rbar32(n)/3.d0
c could get area here
c area(n)=4.d0*pi*rbarte(n)
         rbar10(n)=rbar10(n)/numpar(n)
         rbar32(n)=rbar32(n)/rbarte(n)
         rbar76(n)=rbar76(n)/rbartf(n)
         write(3,4729) nsd(n),rm(n),rn(n),rp(n)
 4729    format(' nsd=',i3,' rm=',1pd12.5,' rn=',1pd12.5,
     1       ' rp=',1pd12.5)
         write(3,2118)numpar(n),vol(n),rbar10(n),rbar32(n),rbar76(n)
 2118    format(1x,' num=',1pd10.4,' vol=',1pd10.4,' r10=',
     1     1pd10.4,' r32=',1pd10.4,' r76=',1pd10.4)
 8831 continue
c
c for plotting and hazy, and qphase for kramers-kronig checks
c ibase = 12 used consistently below
c maximum 2 size distributions; open both even if 2nd not to be used
      open(unit=13,file='tsd13.posd',status='old')
      open(unit=14,file='tsd14.posd',status='old')
c first line in each plot file
      ibase=12
      do 8820 n=1,numsd
         iter=ibase+n
c iris
c         inquire(unit=iter,opened=lopen)
c         if(lopen .eqv. .false.) go to 9001
c         write(3,*) ' writing unit',iter,lopen
         vminus=-vol(n)
c needed as alternative for size in plotting programme,
c and r32 average for block data in cloudy
         write(iter,2007) vminus,rbar32(n)
 8820 continue
c
c this is for choosing ref index or diel function input
      read(2,*) nridf 
      write(3,700)
c****loop here on wavelength
 90   continue
      write(3,700)
c     read in the wavelength (or x or wavenumber)
c     and the complex refractive index or dielectric function of sphere
      read(2,*,end=777) wavlen,nr,ni
c this allows greater accuracy reading in, when nr is close to 1.
c
c checked accuracy - 1986 fundamental ocnstants
c R (infinity scale)
      ryd=9.11267053d-02/wavlen
      ev=ryd*13.6056981d0
      wavnum=1.d0/wavlen
      coefin=wavlen*wavlen/pi
      write(3,1001) wavlen,wavnum,ryd,ev
 1001 format(' vac length=',1pe10.4,' numb=',1pe10.4,
     1    ' ryd=',1pe10.4,' ev=',1pe10.4)
      if(nridf .lt. 0) go to 2000
      if(nridf .eq. 0) nr = nr+1.d0
      call ritodf(nr,ni,eps1,eps2)
      go to 2001
 2000 eps1 = nr
      eps2 = ni
      call dftori(nr,ni,eps1,eps2)
 2001 nr1=nr-1.d0
      write(3,2002) nr,nr1,ni,eps1,eps2
 2002 format(' n=',1pe11.4,' n-1=',1pe11.4,' k=',1pe10.4,
     1  ' ep1=',1pe11.4,' ep2=',1pe10.4)
c     allow for refractive index of medium
      nr = nr/refmed
      ni = ni/refmed
c
c ****************** loop over size initializes  here *********************
c
c zero the items to be integrated
      do 833 n=1,numsd
      rbarex(n)=0.d0
      rbarsc(n)=0.d0
      cosbar(n)=0.d0
      sigext(n)=0.d0
      sigsca(n)=0.d0
      rscmom(n)=0.d0
      rthmom(n)=0.d0
      rsqeff(n)=0.d0
      sigpha(n)=0.d0
      do 833 i=1,jxf
      p1(n,i)=0.d0
      p2(n,i)=0.d0
      p3(n,i)=0.d0
  833 p4(n,i)=0.d0
      jco=0
c     
c************* jco loop (new particle) begins here *********************
 2004 jco=jco+1
      rjco=r(jco)
      if(irlog.ne.0) rjco=dexp(rjco)
      size=rjco
      x=rjco*twopi*wavnum*refmed
      call sinpar(x,qext,qphase,qscatt,ctbrqs,qback,iflag)
c
c 0 normal exit, 1 failure to converge, 2 not even tried
c for exit 1 and 2, stop and figure out why
      if(iflag.gt.0) write(3,*) ' failure iflag',iflag
      if(iflag.gt.0) stop
c
      qscat(jco)=qscatt
c kept in an array for use in avarage sizes below
  707 qabs=qext-qscat(jco)
      cosb=ctbrqs/qscat(jco)
c**************print out single particle results, if desired************
      if(nosini.eq.7) go to 206
      qpr = qext - ctbrqs
      albedo=qscat(jco)/qext
c assumes size in microns, so this is is cm-1
      qovera = qext/(size*1.d-4)
c this is cross-section per unit volume
      coverv=qovera *0.75d0
c cm2 gm-1
      rkappa =coverv/rho 
      write(3,45) qext,qphase,qscatt,qabs
      write(3,2010)cosb,albedo
 45   format(1x,'qext=',1pe12.4,' qphase=',1pe12.4,
     1    ' qscat=',1pe12.4,' qabs=',1pe12.4)
 2010 format(' g =',1pe12.4,' albedo =',1pe12.4)
      write(3,2006) qovera,coverv,rkappa
 2006 format(' qext/a =',1pe12.4,' c/v =',1pe12.4,' kappa =',1pe12.4)
c
c for diagnostic purposes
      write(1,2007) wavnum,qext,qscatt,qabs,
     1    cosb,albedo,qpr,qphase,wavlen,ryd,ev
 206  continue
      if(nosina.eq.7) go to 205
c
c *************************** single size angular functions**************
 3000 format(6(d13.6))
      fjxf=dfloat(jxf)
      x = x / refmed
      if(punch.gt.0.) write(7,3000) x,qext,qscatt,cosb,fjxf
      write(3,367)
  367 format(1x,' scat. angle   ', '      m2         ',    '   m1
     a  ',    '   s21        ',  '     d21        ',    'phase func.   '
     b,   ' polarization  ','     xs','             xd')
      do 774 jj=1,jx
         xi=.5d0*(p1f(jj)+p2f(jj))
c*****q is m2-m1, so q=p1-p2  (note interchange of indices)
c*****p is -q, so p would be p2-p1 (here xp is for q really)
         xp=.5d0*(p1f(jj)-p2f(jj)) /xi
         xs=p3f(jj)/xi
         xd=p4f(jj)/xi
c*****for convention of book, sign of d is changed, so that f43=-d21
         xd=-xd
c*****phase function normalized so that integral over all angles is 4 pi
         xi=4.d0*xi/(x*x*qscatt)
         if(punch.gt.0.) write(7,3000)thetd(jj),xi,xp,xs,xd
      write(3,40) thetd(jj),p1f(jj),p2f(jj),p3f(jj),p4f(jj)
     z ,xi,xp,xs,xd
c      write(3,40) thetd(jj),xp
 774  continue
 40   format(1x,8(1pd12.5,3x),1pd12.5)
c      jmx = max1(1,jx - 1)
      jmx = max(1,jx - 1)
      do 210 j = 1,jmx
c         jj = max1(1,jx - j)
         jj = max(1,jx - j)
         jp=j+jx
         xi=.5d0*(p1b(jj)+p2b(jj))
c*****q is m2-m1, so q=p1-p2  (note interchange of indices)
c*****p is -q, so p would be p2-p1 (here xp is for q really)
         xp=.5d0*(p1b(jj)-p2b(jj)) /xi
         xs=p3b(jj)/xi
         xd=p4b(jj)/xi
c*****for convention of book, sign of d is changed, so that f43=-d21
         xd=-xd
c*****phase function normalized so that integral over all angles is 4 pi
         xi=4.d0*xi/(x*x*qscatt)
         if(punch.gt.0.) write(7,3000)thetd(jp),xi,xp,xs,xd
         write(3,40) thetd(jp),p1b(jj),p2b(jj),p3b(jj),p4b(jj)
     z ,xi,xp,xs,xd
c       write(3,40) thetd(jp),xp
 210  continue
 205  continue
c
c *********************** integrate over size  **************************
      do 829 n=1,numsd
      sw=sizdis(n,jco)*w(jco)
      do 2221 i=1,jx
      p1(n,i)=p1f(i)*sw+p1(n,i)
      p2(n,i)=p2f(i)*sw+p2(n,i)
      p3(n,i)=p3f(i)*sw+p3(n,i)
 2221 p4(n,i)=p4f(i)*sw+p4(n,i)
      jxx=jx+1
c already defined, including case jxoriginal =0     jxf=2*jx-1
      kk=0
      do 2222 i=jxx,jxf
      kk=kk+1
c      k=max1(1,i-2*kk)
      k=max(1,i-2*kk)
      p1(n,i)=p1b(k)*sw+p1(n,i)
      p2(n,i)=p2b(k)*sw+p2(n,i)
      p3(n,i)=p3b(k)*sw+p3(n,i)
 2222 p4(n,i)=p4b(k)*sw+p4(n,i)
      rrr=r(jco)
      if(irlog.ne.0) rrr=dexp(rrr)
      sigext(n)=sigext(n)+rrr*rrr*qext*sw
      sigpha(n)=sigpha(n)+rrr*rrr*qphase*sw
      temp6=rrr*rrr*qscat(jco)*sw
      sigsca(n)=sigsca(n)+temp6
      cosbar(n)=cosbar(n)+cosb*temp6
c cross-section averaged radii
      rbarex(n)=rbarex(n)+rrr*rrr*rrr*qext*sw
      rbarsc(n)=rbarsc(n)+rrr*temp6
  829 continue
      if(jco.ne.ngt)go to 2004
c********************* jco loop (new particle) ends here ***************
c renormalize
      do 828 n=1,numsd
         sigext(n)=pi*sigext(n)
         cover(n)=sigext(n)/vol(n)
         sigpha(n)=pi*sigpha(n)
         sigsca(n)=pi*sigsca(n)
         sigabs(n)=sigext(n) - sigsca(n)
         pizero(n)=sigsca(n)/sigext(n)
         cosbar(n)=pi*cosbar(n)/sigsca(n)
         sigpr(n)=sigext(n) - cosbar(n)*sigsca(n)
         rbarex(n)=rbarex(n)*pi/sigext(n)
         rbarsc(n)=rbarsc(n)*pi/sigsca(n)
         do 2050 i=1,jxf
            p1(n,i)=p1(n,i)*coefin/sigsca(n)
            p2(n,i)=p2(n,i)*coefin/sigsca(n)
            p3(n,i)=p3(n,i)*coefin/sigsca(n)
 2050    p4(n,i)=p4(n,i)*coefin/sigsca(n)
c      write(3,2121) nsd(n),rm(n),rn(n),rp(n)
c 2121 format(1x,' following for nsd=',i3,' rm=',1pd12.5,' rn=',1pd12.5,
c     z     ' rp=',1pd12.5)
c ************ print out polarization matrix if desired ****************
      if(nopol.eq.7) go to 2555
      write(3,3367)
 3367 format(1x,' scat. angle   ','       m2         ','       m1
     a  ','       s21        ','       d21        ','    phase func.   '
     b,'    polarization  ')
      do 2005 i=1,jxf
         xi=(p1(n,i)+p2(n,i))*.5d0
c not dimensioned, not used later        phase(i)=xi
c note sign change, and xs and xd additions would be made here, as for single
         xp=(p2(n,i)-p1(n,i))/xi*.5d0
 2005 write(3,2226)thetd(i),p1(n,i),p2(n,i),p3(n,i),p4(n,i),xi,xp
 2226 format(1x,7(1pd13.4,3x))
c
c **********
 2555 continue
      if(nofin.eq.7) go to 2556
c      write(3,2008)rbarex(n),rbarsc(n),rbar10(n),rbar32(n),rbar76(n)
c     1     1pd12.4,' rbar32=',1pd12.4,' rbar76=',1pd12.4)
c      write(3,2227)sigsca(n),sigext(n),cosbar(n),pizero(n),numpar(n)
      write(3,2227)n,sigsca(n),sigext(n),cosbar(n),pizero(n)
      write(3,2228)sigabs(n),sigpha(n),sigpr(n)
      write(3,2008)rbarex(n),rbarsc(n)
  828 continue
c 2008 format(1x,'rbarex=',1pd12.4,' rbarsc=',1pd12.4,' rbar10=',
 2008 format(1x,'rex=',1pd10.4,' rsc=',1pd10.4)
 2227 format(1x,'nsd=',i2,
     1' sca=',1pd12.5,' ext=',1pd12.5,'  g=',
     1  1pd12.5,' albedo=',1pd12.5)
c 2227 format(1x,'sigsca=',1pd12.5,' sigext=',1pd12.5,' cosbar=',
c     1  1pd12.5,' pizero=',1pd12.5,' numpar=',1pd12.5)
 2228 format(1x,'       abs=',1pd12.5,' pha=',1pd12.5,
     1    ' pr=',1pd12.5)
c 2228 format(1x,'sigabs=',1pd12.5,' sigpha=',1pd12.5,
c     1    ' sigpr=',1pd12.5)
c
 2556 continue
c ********* properties of the size distribution, effective radii etc ********
      if(norpr.eq.7) go to 502 
      do 500 n=1,numsd
      top=0.d0
      do 501 j=1,ngt
      sw=sizdis(n,j)*w(j)
      rj=r(j)
      if(irlog.ne.0) rj=dexp(rj)
  501 top=top+(rj-rbar10(n))**2*sw
  500 rthmm(n)=top/numpar(n)/rbar10(n)**2
      do 7621 n=1,numsd
      do 7622 jco=1,ngt
         rjco=r(jco)
         if(irlog.ne.0) rjco=dexp(rjco)
      t0=rjco*rjco*sizdis(n,jco)*w(jco)
      t1=t0*qscat(jco)*pi
      rscmom(n)=rscmom(n)+((rjco-rbarsc(n))**2)*t1
      rsqeff(n)=rsqeff(n)+((rjco-rbar32(n))**2)*t0
 7622 rthmom(n)=rthmom(n)+((rjco-rbarsc(n))**3)*t1
      rscmom(n)=rscmom(n)/sigsca(n)/rbarsc(n)**2
      rsqeff(n)=rsqeff(n)/rbarte(n)/rbar32(n)**2
 7621 rthmom(n)=rthmom(n)/sigsca(n)/rbarsc(n)**3
      write(3,7631)(rm(n),n=1,numsd)
      write(3,7632)(rn(n),n=1,numsd)
      write(3,7633)(rp(n),n=1,numsd)
      write(3,7634)(rbarsc(n),n=1,numsd)
      write(3,7635)(rscmom(n),n=1,numsd)
      write(3,7636)(rthmom(n),n=1,numsd)
      write(3,7637)(rbar32(n),n=1,numsd)
      write(3,7638)(rsqeff(n),n=1,numsd)
 7631 format('     rm='10(f10.6,1x))
 7632 format('     rn='10(f10.6,1x))
 7633 format('     rp='10(f10.6,1x))
 7634 format(' rbarsc='10(f10.6,1x))
 7635 format(' rscmom='10(f10.6,1x))
 7636 format(' rthmom='10(f10.6,1x))
 7637 format(' rbar32='10(f10.6,1x))
 7638 format(' rsqeff='10(f10.6,1x))
 502  continue
c********************* end of size loop ********************************8
c
c*** input for plotting and hazy, and qphase for kramers-kronig checks****
      zz=0.d0
      ibase=12
      do 8800 n=1,numsd
         iter=ibase+n
         write(iter,2007) wavnum,sigext(n),sigsca(n),
     z        sigabs(n),cosbar(n),pizero(n),sigpr(n),
     z        sigpha(n),wavlen,ryd,ev
 8800 continue
 2007 format(11(1x,1pe13.6))
c* * * * * *
c* * * * * * * * end of wavelength loop * * * * * * * * * * * * * * * * 
      go to 90
c
 777  continue
c put expected endings on plot datasets
      do 8810 n=1,numsd
      write(ibase+n,2007)zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz
 8810 continue
      go to 9002
 9001 write(3,*) ' not enough files specified to match nsd ', iter
 9002 stop
      end
c
