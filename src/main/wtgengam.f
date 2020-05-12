ccc   calls subprocess amplitude
      subroutine wtgengam
      implicit double precision (a-z)
      integer p,i,itot
      complex*16 pp,mm,pm,mp

      include 'polarization.f'
      include 'proc.f'
      include 'ppamp.f'
      include 'mandelstam.f'
      include 'vars.f'
      include 'lbylloop.f'
      
      if(proc.eq.54.or.proc.eq.55)then
         do p=1,pol                       
            call wwpol(p,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp
            pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2
         enddo
      elseif(proc.eq.56.or.proc.eq.57.or.proc.eq.58.or.proc.eq.61
     &   .or.proc.eq.73.or.proc.eq.74.or.proc.eq.75.or.
     &        proc.eq.76)then
         do p=1,pol                       
            call llpol(p,mx,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp    
            pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2   
         enddo
      elseif(proc.eq.59)then
cc         call lightlightpol(1,mx,uh,th,pp,mm,pm,mp)
c         thp=-1d-2
c         call lightlightpol(1,mx,-mx**2-thp,thp,pp,mm,pm,mp)
c         thp=-1d-5*mx**2
c         call lightlightpol(1,mx,-mx**2-thp,thp,pp,mm,pm,mp)
c         thp=-1d-10
c         call lightlightpol(1,mx,uh,th,pp,mm,pm,mp)
c     stop

c$$$         itot=100
c$$$         do i=1,itot
c$$$            mxxmin=5d0
c$$$            mxxmax=15d0
c$$$            mx=mxxmin+(mxxmax-mxxmin)*dble(i)/dble(itot)
c$$$            thp=-1d-2*mx**2
c$$$            uhp=-mx**2-thp
c$$$            call lightlightpol(1,mx,uhp,thp,ppp,mm,pm,mp)
c$$$         enddo
c$$$
c$$$         stop
         
         do p=1,pol                       
            call lightlightpol(p,mx,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp
            pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2
         enddo
      elseif(proc.eq.60)then
          call higgsgam(mx,pp,mm,pm,mp)
            ppa(1)=pp
            mma(1)=mm
            pma(1)=pm
            mpa(1)=mp
            pincarr(1)=cdabs(pp)**2+cdabs(mm)**2           
     &           +cdabs(pm)**2+cdabs(mp)**2
      elseif(proc.eq.68)then
          call alp(mx,pp,mm,pm,mp)
            ppa(1)=pp
            mma(1)=mm
            pma(1)=pm
            mpa(1)=mp
            pincarr(1)=cdabs(pp)**2+cdabs(mm)**2           
     &           +cdabs(pm)**2+cdabs(mp)**2
         elseif(proc.eq.69.or.proc.eq.70)then
            call monop(mx,pp,mm,pm,mp)
            ppa(1)=pp
            mma(1)=mm
            pma(1)=pm
            mpa(1)=mp
            pincarr(1)=cdabs(pp)**2+cdabs(mm)**2           
     &           +cdabs(pm)**2+cdabs(mp)**2            
         elseif(proc.eq.71.or.proc.eq.72)then
            do p=1,pol                       
               call mmpol(p,mx,uh,th,pp,mm,pm,mp)
               ppa(p)=pp
               mma(p)=mm
               pma(p)=pm
               mpa(p)=mp
               pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &              +cdabs(pm)**2+cdabs(mp)**2
            enddo
         endif

      return
      end
