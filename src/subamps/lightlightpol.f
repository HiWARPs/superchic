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
      complex*16 zmes(5),zlep(5)
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
      include 'interpolate.f'

      include 'bottomonium.f' ! include bottomonium resonances
      
      
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
      
      qfarr( 1) = -1d0     !
      qfarr( 2) = -1d0     ! lepton charges
      qfarr( 3) = -1d0     !
      qfarr( 4) = -1d0/3d0    !
      qfarr( 5) =  2d0/3d0    !
      qfarr( 6) = -1d0/3d0    ! quark charges
      qfarr( 7) =  2d0/3d0    !
      qfarr( 8) = -1d0/3d0    !
      qfarr( 9) =  2d0/3d0    !
 
      do i=1,9
         mars(i)=rmf1(i)**2
         marsi(i)=dcmplx(rmf1(i)**2,-1d-30)
         cqfa(i)=qfarr(i)**4                 ! charge factors
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

      if(loop.eq.'tot_quark')then
         lmin=1
         lmax=9
      elseif(loop.eq.'tot_sum')then
         lmin=1
         lmax=9
      elseif(loop.eq.'lep')then
         lmin=1
         lmax=3
      elseif(loop.eq.'quark')then
         lmin=4
         lmax=9
      elseif(loop.eq.'pion')then
         lmin=10
         lmax=10
      elseif(loop.eq.'tot_meson')then
         lmin=1
         lmax=11
      elseif(loop.eq.'meson')then
c         lmin=10
         lmin=4
         lmax=11
      elseif(loop.eq.'w')then
         goto 666
      elseif(loop.eq.'sum')then
         lmin=4
         lmax=9
      endif

         
      do i=lmin,lmax                  ! number of fermion loop particles
         
      mgen2 = mars(i)
      mgen2i = marsi(i)
      cqfactor = cqfa(i)


      if(i.eq.10)then

         
         if(interpolate)then
            zmes(1)=qedamp(1,1,2,2)
            zmes(2)=qedamp(1,1,1,2)
            zmes(3)=qedamp(1,1,1,1)
            zmes(4)=qedamp(1,2,2,1)
            zmes(5)=qedamp(1,2,1,2)
         endif
         
         qedamp(1,1,2,2)=0d0
         qedamp(1,1,1,2)=0d0
         qedamp(1,1,1,1)=0d0
         qedamp(1,2,2,1)=0d0
         qedamp(1,2,1,2)=0d0
         
         mgen2=mpip**2
         mgen2i=dcmplx(mpip**2,-1d-30)

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

         lamexp=2d0
         
         zfs=1d0
         zft=1d0
         zfu=1d0

cccc    Have set form factors to 1 as interpolating via Fermi distribution!        

c$$$         zfs=2d0*mprot**2*c01(0d0,0d0,-sh,zmprot,zmprot,zmprot)
c$$$c         zfs=dexp(-(sh-mpi0**2)**2/lamexp**4) ! Szczurek et al.
c$$$         zfs=zfs**2
c$$$         zft=2d0*mprot**2*c01(0d0,0d0,-t,zmprot,zmprot,zmprot)
c$$$c         zft=dexp((t-mpi0**2)/lamexp**2) ! Szczurek et al.
c$$$         zft=zft**2
c$$$         zfu=2d0*mprot**2*c01(0d0,0d0,-u,zmprot,zmprot,zmprot)
c$$$c         zfu=dexp((u-mpi0**2)/lamexp**2)  ! Szczurek et al.
c$$$         zfu=zfu**2
         
         qedamp(1,1,1,1)=qedamp(1,1,1,1)
     &        +1d0/(sh-mpi0**2+zi*gampi0*sh**2/mpi0**3)/fpi0**2
         qedamp(1,1,1,1)=qedamp(1,1,1,1)
     &        +1d0/(sh-meta**2+zi*gameta*sh**2/meta**3)/feta**2
         qedamp(1,1,1,1)=qedamp(1,1,1,1)
     &        +1d0/(sh-metap**2+zi*gametap*sh**2/metap**3)/fetap**2
         qedamp(1,1,1,1)=qedamp(1,1,1,1)*(-sh**2/4d0/pi**2)/8d0*zfs

c         qedamp(1,1,1,2)=0d0

         qedamp(1,1,2,2)=qedamp(1,1,2,2)
     &        +(sh**2/(sh-mpi0**2+zi*gampi0*sh**2/mpi0**3)*zfs
     &        +t**2/(t-mpi0**2)*zft
     &        +u**2/(u-mpi0**2)*zfu)/fpi0**2
         qedamp(1,1,2,2)=qedamp(1,1,2,2)
     &        +(sh**2/(sh-meta**2+zi*gameta*sh**2/meta**3)*zfs
     &        +t**2/(t-meta**2)*zft
     &        +u**2/(u-meta**2)*zfu)/feta**2
         qedamp(1,1,2,2)=qedamp(1,1,2,2)
     &        +(sh**2/(sh-metap**2+zi*gametap*sh**2/metap**3)*zfs
     &        +t**2/(t-metap**2)*zft
     &        +u**2/(u-metap**2)*zfu)/fetap**2
         qedamp(1,1,2,2)=qedamp(1,1,2,2)/4d0/pi**2/8d0

         qedamp(1,2,2,1)=qedamp(1,2,2,1)
     &        +1d0/(t-mpi0**2)/fpi0**2
         qedamp(1,2,2,1)=qedamp(1,2,2,1)
     &        +1d0/(t-meta**2)/feta**2
         qedamp(1,2,2,1)=qedamp(1,2,2,1)
     &        +1d0/(t-metap**2)/fetap**2
         qedamp(1,2,2,1)=qedamp(1,2,2,1)*(-t**2/4d0/pi**2)/8d0*zft

         qedamp(1,2,1,2)=qedamp(1,2,1,2)
     &        +1d0/(u-mpi0**2)/fpi0**2
         qedamp(1,2,1,2)=qedamp(1,2,1,2)
     &        +1d0/(u-meta**2)/feta**2
         qedamp(1,2,1,2)=qedamp(1,2,1,2)
     &        +1d0/(u-metap**2)/fetap**2
         qedamp(1,2,1,2)=qedamp(1,2,1,2)*(-u**2/4d0/pi**2)/8d0*zfu

         if(interpolate)then
            
            xi_fermi=sh*t*u
            eta_fermi=-(sh*t+t*u+sh*u)
            
            x0_fermi=kappa
            w_fermi=1d0/5d0*x0_fermi
            
            if(fermivar.eq.'xi')then
               fermipp=fermi(xi_fermi,x0_fermi,w_fermi)
            elseif(fermivar.eq.'xi/eta')then
               fermipp=fermi(xi_fermi/eta_fermi,x0_fermi,w_fermi)
            else
               print*,'Incorrect fermi var.'
               stop
            endif


            qedamp(1,1,2,2)=zmes(1)*(1d0-fermipp)
     &           +qedamp(1,1,2,2)*fermipp
            qedamp(1,1,1,2)=zmes(2)*(1d0-fermipp)
     &           +qedamp(1,1,1,2)*fermipp
            qedamp(1,1,1,1)=zmes(3)*(1d0-fermipp)
     &           +qedamp(1,1,1,1)*fermipp
            qedamp(1,2,2,1)=zmes(4)*(1d0-fermipp)
     &           +qedamp(1,2,2,1)*fermipp
            qedamp(1,2,1,2)=zmes(5)*(1d0-fermipp)
     &           +qedamp(1,2,1,2)*fermipp

         endif
         
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
c         Fermion boxes
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

      qedamp(1,2,1,2) = qedamp(1,2,1,2) + cqfactor*(
     & -1d0+(u-sh)/t*(b0fusff-b0fqsff)
     & +(4d0*mgen2/t+2d0*(sh*u/t**2-1d0/2d0))*(u*c0ccusfff+sh*c0ccqsfff)
     & -2d0*mgen2*t*(mgen2/t-1d0/2d0)*(d0ccccqstsffff
     &                        +d0ccccqsusffff+d0ccccustsffff)
     & -sh*u*(4d0*mgen2/t+sh*u/t**2-1d0/2d0)*d0ccccqsusffff)

      
      endif
c     Dissection of Lepton boxes
      if(i.eq.3)then
         zlep(1)=qedamp(1,1,2,2)
         zlep(2)=qedamp(1,1,1,2)
         zlep(3)=qedamp(1,1,1,1)
         zlep(4)=qedamp(1,2,2,1)
         zlep(5)=qedamp(1,2,1,2)
         qedamp(1,1,2,2)=0d0
         qedamp(1,1,1,2)=0d0
         qedamp(1,1,1,1)=0d0
         qedamp(1,2,2,1)=0d0
         qedamp(1,2,1,2)=0d0
      endif
      
      enddo

ccccccccc  Normalize

      if(sumr)then
         
         call rhint(1,dsqrt(sh),out1)
         call rhint(2,dsqrt(sh),out2)
         zrh=out1+zi*out2
         zrh=zrh/alpha**2

         znpp=zrh/8d0  ! match to normalization of qedamp
         znpm=zrh/8d0

         xi_fermi=sh*t*u
         eta_fermi=-(sh*t+t*u+sh*u)
         
         x0_fermi=kappa
         w_fermi=x0_fermi/5d0
         
         if(fermivar.eq.'xi')then
            fermipp=fermi(xi_fermi,x0_fermi,w_fermi)
         elseif(fermivar.eq.'xi/eta')then
            fermipp=fermi(xi_fermi/eta_fermi,x0_fermi,w_fermi)
         else
            print*,'Incorrect fermi var.'
            stop
         endif
         fermipm=fermipp

         if(interpolate)then
            qedamp(1,1,1,1)=znpp*fermipp+qedamp(1,1,1,1)*(1d0-fermipp)
            qedamp(1,2,1,2)=znpm*fermipm+qedamp(1,2,1,2)*(1d0-fermipm)
         else
            qedamp(1,1,1,1)=znpp
            qedamp(1,2,1,2)=znpm
         endif

      endif

ccccccccc Bottomonium states
     
      if(bottomonium)then
c     pseudoscalars eta_b_1S eta_b_2S
         fretab1S=5.87d-5 
         fretab2S=5.86d-5
         
         fretab1S=fretab1S/alpha**2 ! match to normalization of qedamp
         fretab2S=fretab2S/alpha**2 ! match to normalization of qedamp

         gametab1S=17.9d-3
         gametab2S=8.34d-3
      
         metab1S=9.399d0
         metab2S=9.999d0

         etab1S=gametab1S/(metab1S**3)
         etab2S=gametab2S/(metab2S**3)
c     scalars chi_b0_1P chi_b0_2P
         frchib01P=5.87d-5
         frchib02P=5.41d-5

         frchib01P=frchib01P/alpha**2  ! match to normalization of qedamp
         frchib02P=frchib02P/alpha**2  ! match to normalization of qedamp

         gamchib01P=3.39d-3
         gamchib02P=3.54d-3

         mchib01P=9.85944d0
         mchib02P=10.2325d0

         chib01P=gamchib01P/(mchib01P**3)
         chib02P=gamchib02P/(mchib02P**3)

         qedamp(1,1,1,1)=qedamp(1,1,1,1)
     &        +(etab1S*fretab1S/(sh-metab1S**2
     &        +zi*etab1S*sh**2)
     &        +etab2S*fretab2S/(sh-metab2S**2
     &        +zi*etab2S*sh**2)
     &        +chib01P*frchib01P/(sh-mchib01P**2
     &        +zi*chib01P*sh**2)
     &        +chib02P*frchib02P/(sh-mchib02P**2
     &        +zi*chib02P*sh**2))*(-sh**2*16d0*pi)/8d0

c     qedamp(1,1,1,2)=0d0

         qedamp(1,1,2,2)=qedamp(1,1,2,2)
     &        +(-etab1S*fretab1S*(sh**2/(sh-metab1S**2
     &        +zi*etab1S*sh**2)
     &        +t**2/(t-metab1S**2)
     &        +u**2/(u-metab1S**2))
     &        -etab2S*fretab2S*(sh**2/(sh-metab2S**2
     &        +zi*etab2S*sh**2)
     &        +t**2/(t-metab2S**2)
     &        +u**2/(u-metab2S**2))
     &        +chib01P*frchib01P*(sh**2/(sh-mchib01P**2
     &        +zi*chib01P*sh**2)
     &        +t**2/(t-mchib01P**2)
     &        +u**2/(u-mchib01P**2))
     &        +chib02P*frchib02P*(sh**2/(sh-mchib02P**2
     &        +zi*chib02P*sh**2)
     &        +t**2/(t-mchib02P**2)
     &        +u**2/(u-mchib02P**2)))*(-16d0*pi)/8d0

         qedamp(1,2,2,1)=qedamp(1,2,2,1)
     &        +(etab1S*fretab1S/(t-metab1S**2)
     &        +etab2S*fretab2S/(t-metab2S**2)
     &        +chib01P*frchib01P/(t-mchib01P**2)
     &        +chib02P*frchib02P/(t-mchib02P**2))*(-t**2*16d0*pi)/8d0

         qedamp(1,2,1,2)=qedamp(1,2,1,2)
     &        +(etab1S*fretab1S/(u-metab1S**2)
     &        +etab2S*fretab2S/(u-metab2S**2)
     &        +chib01P*frchib01P/(u-mchib01P**2)
     &        +chib02P*frchib02P/(u-mchib02P**2))*(-u**2*16d0*pi)/8d0
      endif

ccccccccc

      if(loop.eq.'tot_quark'.or.loop.eq.'tot_sum'.or.loop.eq.'tot_meson'
     &     .or.loop.eq.'lep')then
         
         qedamp(1,1,2,2)=qedamp(1,1,2,2)+zlep(1)
         qedamp(1,1,1,2)=qedamp(1,1,1,2)+zlep(2)
         qedamp(1,1,1,1)=qedamp(1,1,1,1)+zlep(3)
         qedamp(1,2,2,1)=qedamp(1,2,2,1)+zlep(4)
         qedamp(1,2,1,2)=qedamp(1,2,1,2)+zlep(5)

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

      if(loop.eq.'tot_quark'.or.loop.eq.'tot_sum'.or.loop.eq.'w')then
         
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
