ccc   re-initialise incoming nucleon momenta in nucleon-nucleon cms
      subroutine AAinit
      implicit double precision(a-y)
      
      include 'mom.f'
      include 'vars.f'
      include 'varsi.f'
      include 'mp.f'
      include 'mion.f'
      include 'ion.f'
      include 'sAA.f'

      beta=dsqrt(1d0-4d0*mp**2/s)
      saa=an**2*s
      rtsaa=dsqrt(saa)
      
      q(1,1)=0d0
      q(2,1)=0d0
      q(4,1)=rtsaa/2d0
      q(3,1)=rtsaa/2d0*beta
      
      q(1,2)=0d0
      q(2,2)=0d0
      q(4,2)=rtsaa/2d0
      q(3,2)=-rtsaa/2d0*beta
      
      return
      end
 
