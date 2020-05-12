      function betaionex(t)
      implicit double precision(a-y)

      include 'onechannel.f'

      expon=(b*(a-t))**c-(a*b)**c
      
      betaionex=dexp(-expon)
      
      return
      end
