      subroutine ioninit
      implicit double precision(a-y)
      integer isurv
      character*100 dum

      include 'pi.f'
      include 'vars.f'
      include 'intag.f'
      include 'pdfinf.f'
      include 'ionqcd.f'
      include 'beam.f'
      include 'qcd.f'

      pi=dacos(-1d0)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*,'Initialising ion screening effects...'
      
      call ionpars
      call opacpars(rtsnn)
      call rhonorm
      call rhoxycalc
      call tpcalc
      call opacpcalc
      if(beam.eq.'ion')then
         call opacpbcalc
         if(ionqcd.eq.'incoh')then
            call tpqcdcalc
            call s2qcdion
         endif
      elseif(beam.eq.'ionp')then
         call opacpbpcalc
         if(ionqcd.eq.'incoh')then
            call tpqcdcalc
            call s2qcdionp
         endif
      endif
      if(qcd)then
         if(ionqcd.eq.'coh')call screencalc
      else
         call screencalc
      endif
      
      print*,'...finished!'

      return
      end
