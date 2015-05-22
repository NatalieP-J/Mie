c gauss.f  
c utility routines for mieuvxsd.f
      subroutine gaussi(ng,xbot,xtop,xx,a,rr,ww)
      implicit real*8(a-h,o-z)
      dimension  xx(ng),a(ng),rr(ng),ww(ng)
      bpa=(xtop+xbot)/2.d0
      bma=(xtop-xbot)/2.d0
      do 1 i=1,ng
      rr(i)=bpa+bma*xx(ng+1-i)
      ww(i)=bma*a(i)
1     continue
      return
      end
      subroutine gausle(nn,eps,csa,tsa,x,a)
      implicit real*8(a-h,o-z)
      dimension x(nn),a(nn),c(96)
      fn=dfloat(nn)
      nn2=nn/2
      csa=0.d0
      cc=2.d0
      c(1)=0.d0
      do 1 j=2,nn
      fj=dfloat(j)
      c(j)=(fj-1.d0)*(fj-1.d0)/(4.d0*(fj-0.5d0)*(fj-1.5d0))
1     cc=cc*c(j)
      do 12 i=1,nn2
      if (i-1) 12,2,3
c   largest  zero
2     xt=1.d0-2.78d0/(4.d0+fn*fn)
      go to 11
3     if (i-2) 12,4,5
c   second  zero
4     xt=xt-4.1d0*(1.d0+0.06d0*(1.d0-8.d0/fn))*(1.d0-xt)
      go to 11
5     if (i-3) 12,6,7
c   third  zero
6     xt=xt-1.67d0*(1.d0+0.22d0*(1.d0-8.d0/fn))*(x(1)-xt)
      go to 11
c     middle zeros
7     xt=3.d0*(x(i-1)-x(i-2))+x(i-3)
11    iter=0
21    iter=iter+1
24    pn1=1.d0
      pn=xt
      dp1=0.d0
      dpn=1.d0
      do 26 j=2,nn
      q=xt*pn-c(j)*pn1
      dq=xt*dpn-c(j)*dp1  +pn
      pn1=pn
      pn=q
      dp1=dpn
      dpn=dq
26    continue
25    d=pn/dpn
      xt=xt-d
      if (dabs(d)-eps) 23,23,22
22    if (iter-20) 21,23,23
23    x(i)=xt
      x(nn+1-i)=-xt
      a(i)=cc/(dpn*pn1)
      a(nn+1-i)=a(i)
      csa=csa+a(i)
12    continue
      tsa=1.d0-csa
      if(dabs(tsa).gt.1.d-16) write(3,30) csa,tsa
30    format(1x,'gausle error',2d26.16)
      return
      end
