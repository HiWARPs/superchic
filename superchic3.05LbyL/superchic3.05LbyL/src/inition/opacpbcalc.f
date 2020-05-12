      subroutine opacpbcalc
      implicit double precision(a-y)
      integer i
      
      include 'opacpbpars.f'
      include 'ion.f'
      include 'pi.f'
      include 'ionqcd.f'
      include 'qcd.f'

      ioppb=200

      btmax=5d0*rzg
      btmin=1d0*rzg
      
      hb=btmax/dble(ioppb)

      if(ionqcd.eq.'incoh')sigin=6.5d0 ! fm^2
      if(ionqcd.eq.'coh')sigin=6.74d0 ! fm^2
      
      do i=1,ioppb+1

         bt=(dble(i)-1d0)*hb
         
         opacpbarr(1,i)=bt
         opac=opacpb(bt)
         if(qcd)opac=opac*(1d0-1d0/nshell)**2
         if(opac.gt.100d0)then
            opacpbarr(2,i)=0d0
         else
            opacpbarr(2,i)=dexp(-opac/2d0)
         endif
         
      enddo
      
      return
      end
