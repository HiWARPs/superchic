ccc   gamgam --> ALP amplitude
      subroutine alp(mx,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp

      include 'pi.f'
      include 'ewpars.f'
      include 'norm.f'
      include 'gax.f'
      
      norma=gax/2d0*mx**2
      norma=norma*dsqrt(conv)
      norma=norma/2d0   ! to get correct M^2/4
      
      pp=norma
      mm=-norma
      
      if(alpt.eq.'ps')then
         mm=-norma
      elseif(alpt.eq.'sc')then
         mm=norma
      endif
      pm=0d0
      mp=0d0
      
      return
      end







