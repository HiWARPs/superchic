      function screen(qt)
      implicit double precision(a-y)
      integer iphi,ibt,nphi,nbt

      include 'pi.f'
      include 'ion.f'
      include 'beam.f'

      btmax=10d0*rzg
      btmax=100d0
      
      nbt=10000
      hbt=btmax/dble(nbt)

      sum=0d0
      
      do ibt=1,nbt

         bt=(dble(ibt)-0.5d0)*hbt
         
         if(beam.eq.'ion')wt=1d0-opacpbint(bt)
         if(beam.eq.'ionp')wt=1d0-opacpbpint(bt)
         
c$$$         if(bt.gt.2d0*rzg)then
c$$$            wt=0d0
c$$$         else
c$$$            wt=1d0
c$$$         endif
         
         wt=wt*besj0(qt*bt)
         wt=-wt*bt*hbt/2d0/pi
         
         sum=sum+wt
         
      enddo

      screen=sum
      
      return
      end
