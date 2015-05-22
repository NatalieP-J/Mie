c mieuvx
c Oct 1988 for UV - X-ray extinction, including anomalous diffraction check
c     this version reads in real and imaginary parts of the refractive
c     index, with imaginary part positive (nridf= 1) or nr-1 (nridf=0) or
c     real and imaginary parts of the dielectric function (nridf = -1)
c Dec 1988: added qback; approximation for small x; 
c qphase, better convergence checking
c
c in anomalous diffraction: qext and qabs calculated - qscatt by subtraction
c in rayleigh-gans:         qscatt and qabs calculated
c in mie:                   qext and qscatt calculated
c
      implicit real*8 (a-h,o-z)
      real*8 ni,nr,nr1
      common/cmn/ p1f(91 ),p2f(91 ),p3f(91 ),p4f(91 ),p1b(91 ),p2b(91 ),
     zp3b(91 ),p4b(91 ),cstht(91 ),si2tht(91 ),thetd(181 ),
     znr,ni,jx,nosini
c
c simple association of i/o logical units with global names, in fortran
c link with -lI77 explicitly first
c see ioinit(3F); see fortran library routines, page 173
c sample shell usage to attach to unit 1:
c setenv fort01 iofilename in C-shell
c or fort01=iofilename; export fort01 in Bourne shell
c
c change shell script to get around this f77 construction
c      call ioinit(.true.,.false.,.false.,'fort',.false.)
c
c Oct 2014, standard i/o to go with new shell script mieuvx.sh
c temporary file names
      open(unit=1,file='t1.po',status='old')
c 1 plotting output for plotq and  for hazy
      open(unit=2,file='t2.refr',status='old')
c 2 input wavelength and refractive index
      open(unit=3,file='t3.pr',status='old')
c 3 printed output
      open(unit=4,file='t4.uvx',status='old')
c 4 input mie, rho, sizes
c
c*****punch is positive for punching
      zz=0.d0
      punch=-1.
 700  format(1x)
      pi=3.141592653589793d0
c	refractive index of surrounding medium (real)
      refmed = 1.d0
      twopi=2.d0*pi
      read(4,*) nosina,nosini
      write(3,8320) nosina,nosini
 8319 format(i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x,i3,7x)
 8320 format(' nosina=',i2,' nosini=',i2)
c*********************compute scattering angles*************************
c jx=0, gives jx=1, jxf=2, angles 0, 180 
c jx=1, gives jx=2, jxf=3, angles 0,90, 90,180, with 90 actually computed twice
c in sinpar (set up for general case which might not hit 90)
      jx=0
      jx=90
      dth=0.d0
      if(jx.gt.0) dth=90.d0/dfloat(jx)
      jx=jx+1
c      jxf=max1(2,2*jx-1)
      jxf=max(2,2*jx-1)
      do 216 j=1,jx
      thetd(j)=dfloat(j-1)*dth
      cstht(j)=dcos(thetd(j)*pi/180.d0)
      c=cstht(j)
      si2tht(j)=1.d0-c*c
      l=jxf-(j-1)
 216  thetd(l)=180.d0-thetd(j)
      write(3,698) jxf
 698  format(' number of scattering angles = ',i3)
c
      read(4,*) rho
      write(3,2003) rho
 2003 format(' density =', f12.6)
c*****loop to here for new size, when given negative wavelength*****
 776  continue
      write(3,700)
c     read in the radius or set it in same units as wavelength
      read(4,*,end=777) size
c*****size = 1.d0 would change the wavenumber into the size parameter x
      write(3,775) size
 775  format(' radius = ',1pe12.6)
      write(1,2007) size,size
c this is for choosing ref index or diel function input
      read(2,*) nridf 
      write(3,700)
c****loop here on wavelength
 90   continue
      write(3,700)
c     read in the wavelength (or x or wavenumber)
c     and the complex refractive index or dielectric function of sphere
      read(2,*,end=877) wavlen,nr,ni
      go to 878
c put spacer for different size in output file
 877  write(1,2007)zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz
      rewind 2
      go to 776
c this allows greater accuracy reading in, when nr is close to 1.
 878  continue
c
c checked accuracy - 1986 fundamental ocnstants
c R (infinity scale)
      ryd=9.11267053d-02/wavlen
      ev=ryd*13.6056981d0
      wavnum=1.d0/wavlen
      x= size *wavnum*twopi
      write(3,1001) wavlen,wavnum,x
 1001 format(' vacuum wavelength =',1pe12.4,' wavenumber =',1pe12.4,
     1 ' x =',1pe12.4)
      write(3,1002) ryd,ev
 1002 format(' rydberg =',1pe12.4,' ev =',1pe12.4)
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     allow for refractive index of medium
      nr = nr/refmed
      ni = ni/refmed
      x = x * refmed
      call sinpar(x,qext,qphase,qscatt,ctbrqs,qback,iflag)
c 0 normal exit, 1 failure to converge, 2 not even tried
c for exit 2, see if anomalous diffraction is available
      if(iflag.gt.0) write(3,*) ' failure iflag',iflag
      if(iflag.eq.2) go to 1707
c pass exit 1 through here, but do not put result in plot file (write(1))
 707  qabs=qext-qscatt
      cosb=ctbrqs/qscatt
      qpr = qext - ctbrqs
      albedo=qscatt/qext
c assumes size in microns, so this is is cm-1
      qovera = qext/(size*1.d-4)
c this is cross-section per unit volume
      coverv=qovera *0.75d0
c cm2 gm-1
      rkappa =coverv/rho 
c**************print out single particle results, if desired************
      if(nosini.eq.7) go to 206
      write(3,45) qext,qphase,qscatt,qabs
      og=1.d0-cosb
      write(3,2010)cosb,albedo,og
 45   format(1x,'qext=',1pe12.4,' qphase=',1pe12.4,
     1    ' qscat=',1pe12.4,' qabs=',1pe12.4)
 2010 format(' g =',1pe12.4,' albedo =',1pe12.4,' 1-g =',1pe12.4)
      write(3,2006) qovera,coverv,rkappa
 2006 format(' qext/a =',1pe12.4,' c/v =',1pe12.4,' kappa =',1pe12.4)
c
      if(x.lt.100.d0) go to 134
      qblack=x*ni
      if(qblack.lt.0.1d0) go to 134
c
c a check, using backscattering approximation for large x,
c see bohren and huffman p. 123 and  p. 478
      qbackl=((nr-1.d0)*(nr-1.d0)+ni*ni)/((nr+1.d0)*(nr+1.d0)+ni*ni)
      qbackl=qback/qbackl -1.d0
      write(3,2016) qback,qbackl,x,qblack
 2016 format(' qback =',1pe12.4,' rel -1 =',1pe12.4,2(1x,1pd12.4))
 134  continue
c
c*** input for plotting and hazy, and qphase for kramers-kronig checks****
      if(iflag.eq.0) write(1,2007) wavnum,qext,qscatt,qabs,
     1    cosb,albedo,qpr,qphase,wavlen,ryd,ev
 2007 format(11(1x,1pe13.6))
c* * * * * *
 206  continue
      if(nosina.eq.7) go to 205
c
c *************************** angular functions**************
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
         write(3,40) thetd(jj),p1f(jj),p2f(jj),p3f(jj),p4f(jj),
     z xi,xp,xs,xd
c         write(3,40) thetd(jj),xp
 774     continue
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
         write(3,40) thetd(jp),p1b(jj),p2b(jj),p3b(jj),p4b(jj),
     z xi,xp,xs,xd
c         write(3,40) thetd(jp),xp
 210  continue
c
 205  continue
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * * *
c
c anomalous diffraction -- x >> 1 and |m-1| << 1 but any phase shift
c
 1707 rmodm1 = dsqrt(nr1*nr1+ni*ni)
      if(x.lt.100.d0 .or. rmodm1.gt.0.001d0) go to 2060
c*******print out single particle results, if desired*****
      if(nosini.eq.7) go to 2060
      write(3,*) ' anomalous diffraction '
      delta= -nr1
      beta=ni
      call anomal(x,aqext,aqabs,aqphas,xistar,delta,beta)
      aqscat=aqext-aqabs
c      cosb=ctbrqs/qscatt
c      qpr = qext - ctbrqs
      cosb=zz
      qpr=zz
c??
c g is almost unity for large x and therefore 1-g is very small
      cosb=1.
      qpr=1.e-36
c
      albedo=aqscat/aqext
c assumes size in microns, so this is is cm-1
c      qovera = aqext/(size*1.d-4)
c      rkappa = qovera *0.75d0 / rho
      write(3,45) aqext,aqphas,aqscat,aqabs
      write(3,3010)xistar,albedo
 3010 format(' xistar =',1pe12.4,' albedo =',1pe12.4)
c      write(3,2006) qovera,rkappa
      if(iflag.eq.0) go to 1708
      write(3,*) iflag,
     1    ' iflag failure: plot anomalous diffraction'
c don't do this substitution for a size integration
      write(1,2007) wavnum,aqext,aqscat,aqabs,cosb,albedo,qpr,
     1    aqphas,wavlen,ryd,ev
      if(iflag.eq.2) go to 2060
c  do comparison with mie, whether iflag=0 or 1
 1708 cqe=dabs(aqext/qext -1.d0)
      cqs=dabs(aqscat/qscatt -1.d0)
      cqa=0.d0
      if(qabs.gt.1.d-15) cqa=dabs(aqabs/qabs -1.d0)
      cqph=dabs(aqphas/qphase -1.d0)
      cmax=dmax1(cqe,cqs,cqa,cqph)
c
c crude 0.1 percent check on worst case
      if(cmax.lt.1.d-03) go to 2060
      write(3,*)  iflag,
     1    ' failure: poor anomalous diffraction -- mie agreement '
      write(3,3011) cqe,cqs,cqa,cqph
 3011 format(' cqe,cqs,cqa,cqph=  ',4(1x,1pe12.4))
c **************end of anomalous diffraction**************
c
 2060 continue
      go to 90
  777 stop
      end
c
