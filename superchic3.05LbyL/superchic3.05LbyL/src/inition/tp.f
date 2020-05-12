      function tpz(qt)
      implicit double precision(a-y)
      integer n,ntot

      include 'pi.f'
      include 'ion.f'
      
      sum=0d0

      btmax=rzg*10d0
      
      ntot=100000
      hb=btmax/dble(ntot)
      
      do n=1,ntot

         bt=(dble(n)-0.5d0)*hb

         wt=rhoxyint(1,bt)
         wt=wt*bt*besj0(bt*qt)
         wt=wt*2d0*pi*hb

         sum=sum+wt
         
      enddo

      tpz=sum

      return
      end

      function tpn(qt)
      implicit double precision(a-y)
      integer n,ntot

      include 'pi.f'
      include 'ion.f'
  
      
      sum=0d0

      btmax=rzg*10d0
      
      ntot=100000
      hb=btmax/dble(ntot)
      
      do n=1,ntot

         bt=(dble(n)-0.5d0)*hb

         wt=rhoxyint(2,bt)
         wt=wt*bt*besj0(bt*qt)
         wt=wt*2d0*pi*hb

         sum=sum+wt
         
      enddo
      
      tpn=sum

      return
      end
