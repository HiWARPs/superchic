      subroutine lightlightpol(p,mu,u,t,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp
      integer i,p
      complex*16 b0fqsww,b0ftsww,b0fusww
      complex*16 c0ccqswww,c0cctswww,c0ccuswww
      complex*16 d0ccccqstswwww,d0ccccqsuswwww,d0ccccustswwww
      complex*16 b0fqsff,b0ftsff,b0fusff
      complex*16 c0ccqsfff,c0cctsfff,c0ccusfff
      complex*16 d0ccccqstsffff,d0ccccqsusffff,d0ccccustsffff
      complex*16 wm2
      complex*16 ztest1,ztest2,ztest3,ztest4
      complex*16 zmprot,zfs,zft,zfu
      complex*16 znpp,znpm,zrh,zppt,zpmt
      double precision qfarr(12)
      integer lmin,lmax
      
      include 'lbylamps.f'
      include 'lbylfac.f'
      include 'lep.f'
      include 'pvfun.f'
      
      include 'pi.f'
      include 'zi.f'
      include 'vars.f'
      include 'ewpars.f'
      include 'norm.f'
      include 'mpip.f'
      include 'lbylloop.f'
     
      
      mz=91.1876d0
      stw=0.48076d0
      ctw=dsqrt(1d0-stw**2)
      mw=mz*ctw
      rwm2=mw**2
      wm2 = dcmplx(rwm2,-1d-30)

      if(mu.gt.2d0*mb)then
         qf=35d0/81d0
      else
         qf=34d0/81d0
      endif

      qf=qf*3d0
      ql=3d0
      
      norm=(8d0*alpha**2)**2
      norm=norm/16d0/pi/mx**2
      norm=norm/4d0
      norm=norm*conv
      
      norm=dsqrt(norm)

      sh=mx**2
      
ccccccccccccccc
      
      qfarr( 1) = -1d0
      qfarr( 2) = -1d0
      qfarr( 3) = -1d0
      qfarr( 4) = -1d0/3d0
      qfarr( 5) =  2d0/3d0
      qfarr( 6) = -1d0/3d0
      qfarr( 7) =  2d0/3d0
      qfarr( 8) = -1d0/3d0
      qfarr( 9) =  2d0/3d0
 
      do i=1,9
         mars(i)=rmf1(i)**2
         marsi(i)=dcmplx(rmf1(i)**2,-1d-30)
         cqfa(i)=qfarr(i)**4
         if(i.gt.3)cqfa(i)=cqfa(i)*3d0
      enddo

      if(p.eq.1)then

      qedamp(1,1,1,1) = 0d0
      qedamp(1,1,1,2) = 0d0
      qedamp(1,1,2,2) = 0d0
      qedamp(1,2,1,2) = 0d0
      qedamp(1,2,2,1) = 0d0

      zppt=0d0
      zpmt=0d0

      if(loop.eq.'tot')then
         lmin=1
         lmax=9
      elseif(loop.eq.'lep')then
         lmin=1
         lmax=3
      elseif(loop.eq.'quark'.or.loop.eq.'sum')then
         lmin=4
c     lmax=9
         lmin=4
         lmax=9
      elseif(loop.eq.'pion')then
         lmin=10
         lmax=10
      elseif(loop.eq.'meson')then
         lmin=11
         lmax=11
      elseif(loop.eq.'w')then
         goto 666
      elseif(loop.eq.'sum')then
         call rhint(1,mx,out1)
         call rhint(2,mx,out2)
         znsum=out1+zi*out2
         lmin=4
         lmax=9
      endif
         
      do i=lmin,lmax                  ! number of fermion loop particles
         
      mgen2 = mars(i)
      mgen2i = marsi(i)
      cqfactor = cqfa(i)


      if(i.eq.10)then
         
         mgen2=mpip**2
         mgen2i=dcmplx(mpip**2,-1d-30)

c$$$         mgen2=mars(1)
c$$$         mgen2i=marsi(1)

         cqfactor=1d0

         b0fqsff = b0f2m(sh,mgen2i)
         b0ftsff = b0f2m(t,mgen2i)
         b0fusff = b0f2m(u,mgen2i)
         
         
         c0ccqsfff = -c01(0d0,0d0,-sh,mgen2i,mgen2i,mgen2i)
         c0cctsfff = -c01(0d0,0d0,-t,mgen2i,mgen2i,mgen2i)
         c0ccusfff = -c01(0d0,0d0,-u,mgen2i,mgen2i,mgen2i)
         d0ccccqstsffff = d0404M(sh,t,mgen2i)
         d0ccccqsusffff = d0404M(sh,u,mgen2i)
         d0ccccustsffff = d0404M(u,t,mgen2i)

         qedamp(1,1,2,2) = qedamp(1,1,2,2) - 0.5d0*(
     &        +1d0-2d0*mgen2**2*(d0ccccqstsffff+d0ccccqsusffff
     &        +d0ccccustsffff))

      qedamp(1,1,1,2) = qedamp(1,1,1,2) - 0.5d0*(
     & +1d0-mgen2*(sh**2+t**2+u**2)*(c0ccqsfff/u/t
     &                              +c0cctsfff/sh/u+c0ccusfff/sh/t)
     & -mgen2*((2d0*mgen2+sh*t/u)*d0ccccqstsffff
     &         +(2d0*mgen2+sh*u/t)*d0ccccqsusffff
     &         +(2d0*mgen2+t*u/sh)*d0ccccustsffff))

      qedamp(1,1,1,1) = qedamp(1,1,1,1) + 0.5d0*(
     &     1d0-(u-t)/sh*(b0fusff-b0ftsff)
     &     +2d0*mgen2*u*t/sh*d0ccccustsffff
     &     -2d0*mgen2/sh*(1d0+u*t/2d0/mgen2/sh)*(2d0*t*c0cctsfff
     &     +2d0*u*c0ccusfff-u*t*d0ccccustsffff)
     &     +2d0*mgen2**2*(d0ccccqstsffff+d0ccccqsusffff+d0ccccustsffff))
      
!      s <-> u
      
      qedamp(1,2,2,1) = qedamp(1,2,2,1) + 0.5d0*(
     &     1d0-(sh-t)/u*(b0fqsff-b0ftsff)
     &     +2d0*mgen2*sh*t/u*d0ccccqstsffff
     &     -2d0*mgen2/u*(1d0+sh*t/2d0/mgen2/u)*(2d0*t*c0cctsfff
     &     +2d0*sh*c0ccqsfff-sh*t*d0ccccqstsffff)
     &     +2d0*mgen2**2*(d0ccccqstsffff+d0ccccqsusffff+d0ccccustsffff))
      
!      s <-> t
      
      qedamp(1,2,1,2) = qedamp(1,2,1,2) + 0.5d0*(
     &     1d0-(u-sh)/t*(b0fusff-b0fqsff)
     &     +2d0*mgen2*u*sh/t*d0ccccqsusffff
     &     -2d0*mgen2/t*(1d0+u*sh/2d0/mgen2/t)*(2d0*sh*c0ccqsfff
     &     +2d0*u*c0ccusfff-u*sh*d0ccccqsusffff)
     &     +2d0*mgen2**2*(d0ccccqstsffff+d0ccccqsusffff+d0ccccustsffff))

      elseif(i.eq.11)then

         cqfactor=1d0

         mprot=0.938d0
         zmprot=dcmplx(mprot**2,-1d-30)

         fpi0=0.092d0
         feta=0.093d0
         fetap=0.074d0

         gampi0=7.8d-9
         gameta=1.31d-6
         gametap=0.197d-3

         mpi0=134.9766d-3
         meta=547.862d-3
         metap=957.78d-3

         zfs=1d0
         zft=1d0
         zfu=1d0

         zfs=2d0*mprot**2*c01(0d0,0d0,-sh,zmprot,zmprot,zmprot)
         zfs=zfs**2
         zft=2d0*mprot**2*c01(0d0,0d0,-t,zmprot,zmprot,zmprot)
         zft=zft**2
         zfu=2d0*mprot**2*c01(0d0,0d0,-u,zmprot,zmprot,zmprot)
         zfu=zfu**2
         
         qedamp(1,1,1,1)=1d0/(sh-mpi0**2+zi*mpi0*gampi0)/fpi0**2
         qedamp(1,1,1,1)=qedamp(1,1,1,1)
     &        +1d0/(sh-meta**2+zi*meta*gameta)/feta**2
         qedamp(1,1,1,1)=qedamp(1,1,1,1)
     &        +1d0/(sh-metap**2+zi*metap*gametap)/fetap**2
         qedamp(1,1,1,1)=qedamp(1,1,1,1)*(-sh**2/4d0/pi**2)/8d0*zfs

         qedamp(1,1,1,2)=0d0

         qedamp(1,1,2,2)=(sh**2/(sh-mpi0**2+zi*mpi0*gampi0)*zfs
     &        +t**2/(t-mpi0**2+zi*mpi0*gampi0)*zft
     &        +u**2/(u-mpi0**2+zi*mpi0*gampi0)*zfu)/fpi0**2
         qedamp(1,1,2,2)=qedamp(1,1,2,2)
     &        +(sh**2/(sh-meta**2+zi*mpi0*gameta)*zfs
     &        +t**2/(t-meta**2+zi*mpi0*gameta)*zft
     &        +u**2/(u-meta**2+zi*mpi0*gameta)*zfu)/feta**2
         qedamp(1,1,2,2)=qedamp(1,1,2,2)
     &        +(sh**2/(sh-metap**2+zi*mpi0*gametap)*zfs
     &        +t**2/(t-metap**2+zi*mpi0*gametap)*zft
     &        +u**2/(u-metap**2+zi*mpi0*gametap)*zfu)/fetap**2
         qedamp(1,1,2,2)=qedamp(1,1,2,2)/4d0/pi**2/8d0

         qedamp(1,2,2,1)=1d0/(t-mpi0**2+zi*mpi0*gampi0)/fpi0**2
         qedamp(1,2,2,1)=qedamp(1,2,2,1)
     &        +1d0/(t-meta**2+zi*meta*gameta)/feta**2
         qedamp(1,2,2,1)=qedamp(1,2,2,1)
     &        +1d0/(t-metap**2+zi*metap*gametap)/fetap**2
         qedamp(1,2,2,1)=qedamp(1,2,2,1)*(-t**2/4d0/pi**2)/8d0*zft

         qedamp(1,2,1,2)=1d0/(u-mpi0**2+zi*mpi0*gampi0)/fpi0**2
         qedamp(1,2,1,2)=qedamp(1,2,1,2)
     &        +1d0/(u-meta**2+zi*meta*gameta)/feta**2
         qedamp(1,2,1,2)=qedamp(1,2,1,2)
     &        +1d0/(u-metap**2+zi*metap*gametap)/fetap**2
         qedamp(1,2,1,2)=qedamp(1,2,1,2)*(-u**2/4d0/pi**2)/8d0*zfu

         
      elseif (mgen2.eq.0) then
         
      qedamp(1,1,2,2) = qedamp(1,1,2,2) + cqfactor  

      qedamp(1,1,1,2) = qedamp(1,1,1,2) + cqfactor

      qedamp(1,1,1,1) = qedamp(1,1,1,1) + cqfactor*(  ! ++++
     & -1d0 + (t-u)/sh*(LOG(-u/sh) - LOG(-t/sh))
     & -(1d0/2 - u*t/sh**2)*((LOG(-u/sh)-LOG(-t/sh))**2 + pi**2))
      
      qedamp(1,2,2,1) = qedamp(1,2,2,1) + cqfactor*(  ! +-+-
     & -1d0 - zi*pi*(t-sh)/u
     & -((1d0 + zi*pi)*(t-sh)/u + 2d0*zi*pi*(t/u)**2)*LOG(-t/sh)
     & -(1d0/2 - sh*t/u**2)*LOG(-t/sh)**2)

      qedamp(1,2,1,2) = qedamp(1,2,1,2) + cqfactor*(   ! +--+
     & -1d0 - zi*pi*(u-sh)/t
     & -((1d0 + zi*pi)*(u-sh)/t + 2d0*zi*pi*(u/t)**2)*LOG(-u/sh)
     & -(1d0/2 - sh*u/t**2)*LOG(-u/sh)**2) 

      else

c$$$         if(loop.eq.'sum')then
c$$$
c$$$            thp=-sh*1d-20
c$$$            uhp=-thp-sh
c$$$
c$$$            b0fqsff = b0f2m(sh,mgen2i)
c$$$            b0ftsff = b0f2m(thp,mgen2i)
c$$$            b0fusff = b0f2m(uhp,mgen2i)
c$$$            
c$$$      
c$$$            c0ccqsfff = -c01(0d0,0d0,-sh,mgen2i,mgen2i,mgen2i)
c$$$            c0cctsfff = -c01(0d0,0d0,-thp,mgen2i,mgen2i,mgen2i)
c$$$            c0ccusfff = -c01(0d0,0d0,-uhp,mgen2i,mgen2i,mgen2i)
c$$$            d0ccccqstsffff = d0404M(sh,thp,mgen2i)
c$$$            d0ccccqsusffff = d0404M(sh,uhp,mgen2i)
c$$$            d0ccccustsffff = d0404M(uhp,thp,mgen2i)
c$$$
c$$$            znpp=(
c$$$     &           -1d0+(uhp-thp)/sh*(b0fusff-b0ftsff)
c$$$     &           +(4d0*mgen2/sh+2d0*(thp*uhp/sh**2-1d0/2d0))*
c$$$     &           (uhp*c0ccusfff+thp*c0cctsfff)
c$$$     &           -2d0*mgen2*sh*(mgen2/sh-1d0/2d0)*(d0ccccqstsffff
c$$$     &                               +d0ccccqsusffff+d0ccccustsffff)
c$$$     &           -thp*uhp*(4d0*mgen2/sh+thp*uhp/sh**2-1d0/2d0)*
c$$$     &           d0ccccustsffff)
c$$$
c$$$            znpm=(
c$$$     &           -1d0+(sh-thp)/uhp*(b0fqsff-b0ftsff)
c$$$     &           +(4d0*mgen2/uhp+2d0*(thp*sh/uhp**2-1d0/2d0))*
c$$$     &           (sh*c0ccqsfff+thp*c0cctsfff)
c$$$     &           -2d0*mgen2*uhp*(mgen2/uhp-1d0/2d0)*(d0ccccqstsffff
c$$$     &           +d0ccccqsusffff+d0ccccustsffff)
c$$$     &           -thp*sh*(4d0*mgen2/uhp+thp*sh/uhp**2-1d0/2d0)*
c$$$     &           d0ccccqstsffff)
c$$$
c$$$c            zppt=zppt+znpp*cqfactor*8d0
c$$$c           zpmt=zpmt+znpm*cqfactor*8d0
c$$$
c$$$
c$$$               call mint(1,1,i,dsqrt(sh),out1)
c$$$               call mint(1,2,i,dsqrt(sh),out2)
c$$$               zppt=zppt+out1+zi*out2
c$$$               call mint(2,1,i,dsqrt(sh),out1)
c$$$               call mint(2,2,i,dsqrt(sh),out2)
c$$$               zpmt=zpmt+out1+zi*out2
c$$$               
c$$$c               znpp=zppt
c$$$c               znpm=zpmt
c$$$            
c$$$c     if(mx.gt.5d0)then
c$$$c$$$            if(i.le.8)then
c$$$c$$$
c$$$c$$$               call mint(1,1,i,dsqrt(sh),out1)
c$$$c$$$               call mint(1,2,i,dsqrt(sh),out2)
c$$$c$$$               zppt=out1+zi*out2
c$$$c$$$               call mint(2,1,i,dsqrt(sh),out1)
c$$$c$$$               call mint(2,2,i,dsqrt(sh),out2)
c$$$c$$$               zpmt=out1+zi*out2
c$$$c$$$               
c$$$c$$$c               print*,dsqrt(mgen2),cqfactor,(1d0/3d0)**4*3d0
c$$$c$$$c               print*,mx,znpp*8d0*cqfactor,znpm*8d0*cqfactor
c$$$c$$$c               print*,zppt,zpmt
c$$$c$$$               print*,i,mx,znpp*8d0*cqfactor/zppt,znpm*8d0*cqfactor/zpmt
c$$$c$$$               print*,''
c$$$c$$$            endif
c$$$c     endif
c$$$c               print*,zrh,znpp,znpm
c$$$
c$$$
c$$$c$$$
c$$$c$$$            if(cdabs(znpm).gt.0d0.or.cdabs(znpm).eq.0d0)then
c$$$c$$$            else
c$$$c$$$               print*,i,dsqrt(sh),znpm,zpmt
c$$$c$$$            endif
c$$$            
c$$$
c$$$c$$$            print*,i,dsqrt(sh),znpp,znpm
c$$$c$$$            print*,zppt,zpmt,zrh
c$$$c$$$            print*,''
c$$$c            print*,'tt',znpp/cqfactor,znpm/cqfactor
c$$$            
c$$$         endif
         
      b0fqsff = b0f2m(sh,mgen2i)
      b0ftsff = b0f2m(t,mgen2i)
      b0fusff = b0f2m(u,mgen2i)

      
      c0ccqsfff = -c01(0d0,0d0,-sh,mgen2i,mgen2i,mgen2i)
      c0cctsfff = -c01(0d0,0d0,-t,mgen2i,mgen2i,mgen2i)
      c0ccusfff = -c01(0d0,0d0,-u,mgen2i,mgen2i,mgen2i)
      d0ccccqstsffff = d0404M(sh,t,mgen2i)
      d0ccccqsusffff = d0404M(sh,u,mgen2i)
      d0ccccustsffff = d0404M(u,t,mgen2i)
 
      qedamp(1,1,2,2) = qedamp(1,1,2,2) + cqfactor*(
     & +1d0-2d0*mgen2**2*(d0ccccqstsffff+d0ccccqsusffff+d0ccccustsffff))

      qedamp(1,1,1,2) = qedamp(1,1,1,2) + cqfactor*(
     & +1d0-mgen2*(sh**2+t**2+u**2)*(c0ccqsfff/u/t
     &                              +c0cctsfff/sh/u+c0ccusfff/sh/t)
     & -mgen2*((2d0*mgen2+sh*t/u)*d0ccccqstsffff
     &         +(2d0*mgen2+sh*u/t)*d0ccccqsusffff
     &         +(2d0*mgen2+t*u/sh)*d0ccccustsffff))

c$$$      if(loop.eq.'sum'.and.i.le.8)then
c$$$
c$$$         qedamp(1,1,1,1) = qedamp(1,1,1,1) + cqfactor*(
c$$$     & -1d0+(u-t)/sh*(b0fusff-b0ftsff)
c$$$     & +(4d0*mgen2/sh+2d0*(t*u/sh**2-1d0/2d0))*(u*c0ccusfff+t*c0cctsfff)
c$$$     & -2d0*mgen2*sh*(mgen2/sh-1d0/2d0)*(d0ccccqstsffff
c$$$     &                               +d0ccccqsusffff+d0ccccustsffff)
c$$$     & -t*u*(4d0*mgen2/sh+t*u/sh**2-1d0/2d0)*d0ccccustsffff)*znpp
c$$$
c$$$      
c$$$      qedamp(1,2,2,1) = qedamp(1,2,2,1) + cqfactor*(
c$$$     & -1d0+(sh-t)/u*(b0fqsff-b0ftsff)
c$$$     & +(4d0*mgen2/u+2d0*(t*sh/u**2-1d0/2d0))*(sh*c0ccqsfff+t*c0cctsfff)
c$$$     & -2d0*mgen2*u*(mgen2/u-1d0/2d0)*(d0ccccqstsffff
c$$$     &                        +d0ccccqsusffff+d0ccccustsffff)
c$$$     &     -t*sh*(4d0*mgen2/u+t*sh/u**2-1d0/2d0)*d0ccccqstsffff)*znpm


c$$$      
c$$$      else
         
      qedamp(1,1,1,1) = qedamp(1,1,1,1) + cqfactor*(
     & -1d0+(u-t)/sh*(b0fusff-b0ftsff)
     & +(4d0*mgen2/sh+2d0*(t*u/sh**2-1d0/2d0))*(u*c0ccusfff+t*c0cctsfff)
     & -2d0*mgen2*sh*(mgen2/sh-1d0/2d0)*(d0ccccqstsffff
     &                               +d0ccccqsusffff+d0ccccustsffff)
     & -t*u*(4d0*mgen2/sh+t*u/sh**2-1d0/2d0)*d0ccccustsffff)

      
      qedamp(1,2,2,1) = qedamp(1,2,2,1) + cqfactor*(
     & -1d0+(sh-t)/u*(b0fqsff-b0ftsff)
     & +(4d0*mgen2/u+2d0*(t*sh/u**2-1d0/2d0))*(sh*c0ccqsfff+t*c0cctsfff)
     & -2d0*mgen2*u*(mgen2/u-1d0/2d0)*(d0ccccqstsffff
     &                        +d0ccccqsusffff+d0ccccustsffff)
     &     -t*sh*(4d0*mgen2/u+t*sh/u**2-1d0/2d0)*d0ccccqstsffff)

c$$$      endif

      qedamp(1,2,1,2) = qedamp(1,2,1,2) + cqfactor*(
     & -1d0+(u-sh)/t*(b0fusff-b0fqsff)
     & +(4d0*mgen2/t+2d0*(sh*u/t**2-1d0/2d0))*(u*c0ccusfff+sh*c0ccqsfff)
     & -2d0*mgen2*t*(mgen2/t-1d0/2d0)*(d0ccccqstsffff
     &                        +d0ccccqsusffff+d0ccccustsffff)
     & -sh*u*(4d0*mgen2/t+sh*u/t**2-1d0/2d0)*d0ccccqsusffff)

c$$$      write(6,*)i,
c$$$     & -1d0+(u-t)/sh*(b0fusff-b0ftsff)
c$$$     & +(4d0*mgen2/sh+2d0*(t*u/sh**2-1d0/2d0))*(u*c0ccusfff+t*c0cctsfff)
c$$$     & -2d0*mgen2*sh*(mgen2/sh-1d0/2d0)*(d0ccccqstsffff
c$$$     &                               +d0ccccqsusffff+d0ccccustsffff)
c$$$     &     -t*u*(4d0*mgen2/sh+t*u/sh**2-1d0/2d0)*d0ccccustsffff,
c$$$     &-1d0+(sh-t)/u*(b0fqsff-b0ftsff)
c$$$     & +(4d0*mgen2/u+2d0*(t*sh/u**2-1d0/2d0))*(sh*c0ccqsfff+t*c0cctsfff)
c$$$     & -2d0*mgen2*u*(mgen2/u-1d0/2d0)*(d0ccccqstsffff
c$$$     &                        +d0ccccqsusffff+d0ccccustsffff)
c$$$     & -t*sh*(4d0*mgen2/u+t*sh/u**2-1d0/2d0)*d0ccccqstsffff


c$$$      write(6,*)i,  ! ++++
c$$$     & -1d0 + (t-u)/sh*(LOG(-u/sh) - LOG(-t/sh))
c$$$     &     -(1d0/2 - u*t/sh**2)*((LOG(-u/sh)-LOG(-t/sh))**2 + pi**2),
c$$$     & -1d0 - zi*pi*(u-sh)/t
c$$$     & -((1d0 + zi*pi)*(u-sh)/t + 2d0*zi*pi*(u/t)**2)*LOG(-u/sh)
c$$$     &     -(1d0/2d0 - sh*u/t**2)*LOG(-u/sh)**2

c$$$      write(6,*)-1d0 - zi*pi*(u-sh)/t
c$$$     &     -((1d0 + zi*pi)*(u-sh)/t + 2d0*zi*pi*(u/t)**2)*LOG(-u/sh)
c$$$     &     -(1d0/2d0 - sh*u/t**2)*LOG(-u/sh)**2
c$$$      write(6,*) -1d0 - zi*pi*(u-sh)/t
c$$$     &     ,-((1d0 + zi*pi)*(u-sh)/t + 2d0*zi*pi*(u/t)**2)*LOG(-u/sh),
c$$$     &     -(1d0/2d0 - sh*u/t**2)*LOG(-u/sh)**2
c$$$
c$$$      write(6,*)-1d0+2d0*pi*zi*sh/t+zi*pi
c$$$      write(6,*)LOG(-u/sh)*(1d0 + zi*pi)*(u-sh)/t
c$$$     &     ,2d0*zi*pi*(u/t)**2*LOG(-u/sh)
      
      endif
      enddo

ccccccccc  Normalize

      if(loop.eq.'sum')then

c         shp=2d0
         
         kappa_fermi=0.01d0
         eta_fermi=1d0/7d0

c         thp=0d0
         
c         print*,shp
c     &        ,fermi(-thp,kappa_fermi*shp,eta_fermi*kappa_fermi*shp)

         
         call rhint(1,dsqrt(sh),out1)
         call rhint(2,dsqrt(sh),out2)
         zrh=out1+zi*out2
         zrh=zrh/alpha**2

         znpp=zrh/8d0
c         /zppt
         znpm=zrh/8d0
c         /zpmt

c     print*,sh,t,znpp

         znpp=(znpp-1d0)*(fermi(-t,kappa_fermi*sh
     &        ,eta_fermi*kappa_fermi*sh)+fermi(sh+t,kappa_fermi*sh
     &        ,eta_fermi*kappa_fermi*sh))+1d0
         znpm=(znpm-1d0)*(fermi(-t,kappa_fermi*sh
     &        ,eta_fermi*kappa_fermi*sh)+fermi(sh+t,kappa_fermi*sh
     &        ,eta_fermi*kappa_fermi*sh))+1d0

c$$$         znpp=(znpp-1d0)*(fermi(-t,kappa_fermi*sh
c$$$     &        ,eta_fermi*kappa_fermi*sh)+fermi(sh+t,kappa_fermi*sh
c$$$     &        ,eta_fermi*kappa_fermi*sh))+1d0
c$$$         znpm=(znpm-1d0)*(fermi(-t,kappa_fermi*sh
c$$$     &        ,eta_fermi*kappa_fermi*sh)+fermi(sh+t,kappa_fermi*sh
c$$$     &        ,eta_fermi*kappa_fermi*sh))+1d0

c         print*,sh,t,znpp
c         print*,''
         
         qedamp(1,1,1,1)=qedamp(1,1,1,1)*znpp
         qedamp(1,2,1,2)=qedamp(1,2,1,2)*znpm

      endif
         
ccccccccc      
      
      
      qedamp(1,1,2,1) = qedamp(1,1,1,2)
      qedamp(1,2,1,1) = qedamp(1,1,1,2)
      qedamp(1,2,2,2) = qedamp(1,1,1,2)
      qedamp(2,1,1,1) = qedamp(1,2,2,2)
      qedamp(2,1,1,2) = qedamp(1,2,2,1)
      qedamp(2,1,2,1) = qedamp(1,2,1,2)
      qedamp(2,1,2,2) = qedamp(1,2,1,1)
      qedamp(2,2,1,1) = qedamp(1,1,2,2)
      qedamp(2,2,1,2) = qedamp(1,1,2,1)
      qedamp(2,2,2,1) = qedamp(1,1,1,2)
      qedamp(2,2,2,2) = qedamp(1,1,1,1)


cccccccc

      
      
 666  ewamp(1,1,1,1) = 0d0
      ewamp(1,1,1,2) = 0d0
      ewamp(1,1,2,2) = 0d0
      ewamp(1,2,1,2) = 0d0
      ewamp(1,2,2,1) = 0d0

      b0fqsww = b0f2m(sh,wm2)
      b0ftsww = b0f2m(t,wm2)
      b0fusww = b0f2m(u,wm2)
      c0ccqswww = -c01(0d0,0d0,-sh,wm2,wm2,wm2)
      c0cctswww = -c01(0d0,0d0,-t,wm2,wm2,wm2)
      c0ccuswww = -c01(0d0,0d0,-u,wm2,wm2,wm2)
      d0ccccqstswwww = d0404M(sh,t,wm2)
      d0ccccqsuswwww = d0404M(sh,u,wm2)
      d0ccccustswwww = d0404M(u,t,wm2)
 
      ewamp(1,1,2,2) = 
     & +1d0-2d0*rwm2**2*(d0ccccqstswwww+d0ccccqsuswwww+d0ccccustswwww)

      ewamp(1,1,1,2) = 
     & +1d0-rwm2*(sh**2+t**2+u**2)*( 1d0/u/t*c0ccqswww
     &                             +1d0/sh/u*c0cctswww
     &                             +1d0/sh/t*c0ccuswww)
     & -rwm2*((2d0*rwm2+sh*t/u)*d0ccccqstswwww
     &       +(2d0*rwm2+sh*u/t)*d0ccccqsuswwww
     &       +(2d0*rwm2+t*u/sh)*d0ccccustswwww)

      ewamp(1,1,1,1) = 
     & -1d0+(u-t)/sh*(b0fusww-b0ftsww)
     & +(4d0*rwm2/sh+2d0*(t*u/sh**2-4d0/3d0))*(u*c0ccuswww+t*c0cctswww)
     & -(2d0*rwm2*sh*(rwm2/sh-4d0/3d0)+2d0/3d0*sh**2)*
     &                        ( d0ccccqstswwww
     &                         +d0ccccqsuswwww
     &                         +d0ccccustswwww )
     & -t*u*(4d0*rwm2/sh+t*u/sh**2-4d0/3d0)*d0ccccustswwww

      ewamp(1,2,2,1) = 
     & -1d0+(sh-t)/u*(b0fqsww-b0ftsww)
     & +(4d0*rwm2/u+2d0*(t*sh/u**2-4d0/3d0))*(sh*c0ccqswww+t*c0cctswww)
     & -(2d0*rwm2*u*(rwm2/u-4d0/3d0)+2d0/3d0*u**2)*
     &                        ( d0ccccqstswwww
     &                         +d0ccccqsuswwww
     &                         +d0ccccustswwww )
     & -t*sh*(4d0*rwm2/u+t*sh/u**2-4d0/3d0)*d0ccccqstswwww

      ewamp(1,2,1,2) = 
     & -1d0+(u-sh)/t*(b0fusww-b0fqsww)
     & +(4d0*rwm2/t+2d0*(sh*u/t**2-4d0/3d0))*(u*c0ccuswww+sh*c0ccqswww)
     & -(2d0*rwm2*t*(rwm2/t-4d0/3d0)+2d0/3d0*t**2)*
     &                        ( d0ccccqstswwww
     &                         +d0ccccqsuswwww
     &                         +d0ccccustswwww )
     & -sh*u*(4d0*rwm2/t+sh*u/t**2-4d0/3d0)*d0ccccqsuswwww

      ewamp(1,1,2,1) = ewamp(1,1,1,2)
      ewamp(1,2,1,1) = ewamp(1,1,1,2)
      ewamp(1,2,2,2) = ewamp(1,1,1,2)
      ewamp(2,1,1,1) = ewamp(1,2,2,2)
      ewamp(2,1,1,2) = ewamp(1,2,2,1)
      ewamp(2,1,2,1) = ewamp(1,2,1,2)
      ewamp(2,1,2,2) = ewamp(1,2,1,1)
      ewamp(2,2,1,1) = ewamp(1,1,2,2)
      ewamp(2,2,1,2) = ewamp(1,1,2,1)
      ewamp(2,2,2,1) = ewamp(1,1,1,2)
      ewamp(2,2,2,2) = ewamp(1,1,1,1)
      
      endif


           qedamp(1,1,2,1) = qedamp(1,1,1,2)
      qedamp(1,2,1,1) = qedamp(1,1,1,2)
      qedamp(1,2,2,2) = qedamp(1,1,1,2)
      qedamp(2,1,1,1) = qedamp(1,2,2,2)
      qedamp(2,1,1,2) = qedamp(1,2,2,1)
      qedamp(2,1,2,1) = qedamp(1,2,1,2)
      qedamp(2,1,2,2) = qedamp(1,2,1,1)
      qedamp(2,2,1,1) = qedamp(1,1,2,2)
      qedamp(2,2,1,2) = qedamp(1,1,2,1)
      qedamp(2,2,2,1) = qedamp(1,1,1,2)
      qedamp(2,2,2,2) = qedamp(1,1,1,1)
      
      if(p.eq.1)then    
         pp=qedamp(1,1,1,1)
         mm=qedamp(1,1,2,2)
         pm=qedamp(1,1,1,2)
         mp=qedamp(1,1,2,1)


c         print*,1,cdabs(pp),cdabs(mm),cdabs(pm),cdabs(mp)
         
      elseif(p.eq.2)then
         pp=qedamp(2,2,1,1)
         mm=qedamp(2,2,2,2)
         pm=qedamp(2,2,1,2)
         mp=qedamp(2,2,2,1)

c         print*,2,cdabs(pp),cdabs(mm),cdabs(pm),cdabs(mp)
         
      elseif(p.eq.3)then         
         pp=qedamp(1,2,1,1)
         mm=qedamp(1,2,2,2)
         pm=qedamp(1,2,1,2)
         mp=qedamp(1,2,2,1)

c         print*,3,cdabs(pp),cdabs(mm),cdabs(pm),cdabs(mp)
         
      elseif(p.eq.4)then      
         pp=qedamp(2,1,1,1)
         mm=qedamp(2,1,2,2)
         pm=qedamp(2,1,1,2)
         mp=qedamp(2,1,2,1)

c         print*,4,cdabs(pp),cdabs(mm),cdabs(pm),cdabs(mp)

c         stop
         
      endif

      
      efac=1.5d0

      if(loop.eq.'tot'.or.loop.eq.'w')then
         
      if(p.eq.1)then
         pp=pp-ewamp(1,1,1,1)*efac
         mm=mm-ewamp(1,1,2,2)*efac
         pm=pm-ewamp(1,1,1,2)*efac
         mp=mp-ewamp(1,1,2,1)*efac         
      elseif(p.eq.2)then
         pp=pp-ewamp(2,2,1,1)*efac
         mm=mm-ewamp(2,2,2,2)*efac
         pm=pm-ewamp(2,2,1,2)*efac
         mp=mp-ewamp(2,2,2,1)*efac
      elseif(p.eq.3)then         
         pp=pp-ewamp(1,2,1,1)*efac
         mm=mm-ewamp(1,2,2,2)*efac
         pm=pm-ewamp(1,2,1,2)*efac
         mp=mp-ewamp(1,2,2,1)*efac
      elseif(p.eq.4)then         
         pp=pp-ewamp(2,1,1,1)*efac
         mm=mm-ewamp(2,1,2,2)*efac
         pm=pm-ewamp(2,1,1,2)*efac
         mp=mp-ewamp(2,1,2,1)*efac
      endif

      endif

      pp=pp*norm
      mm=mm*norm
      pm=pm*norm
      mp=mp*norm
      
 888  return
      end
