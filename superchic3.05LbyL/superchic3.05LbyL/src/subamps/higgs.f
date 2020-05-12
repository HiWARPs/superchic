ccc   gg --> SM higgs subprocess amplitude
      subroutine higgs(mx,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp

      include 'pi.f'
      include 'ewpars.f'
      include 'normh.f'

      pp=normh
      mm=normh
      mp=0d0
      pm=0d0
      
      return
      end







