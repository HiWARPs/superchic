ccccc EPA form factors (proton)
      subroutine formfacgam(io,t1,t2,out)
      implicit double precision(a-y)
      integer i1,i2,io

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      
      qsq=(x1**2*mp**2+t1)/(1d0-x1)
      qsqp=(x2**2*mp**2+t2)/(1d0-x2)      

ccccccccc

      qsqmin=mp**2*x1**2/(1d0-x1)

      a0e=0.98462d0
      a1e=0.68414d0
      a2e=0.01933d0
      a0m=0.28231d0
      a1m=1.34919d0
      a2m=0.55473d0
      ge=(a0e/(1d0+qsq/a1e)**2+(1d0-a0e)/(1d0+qsq/a2e)**2)**2
      gm=(a0m/(1d0+qsq/a1m)**2+(1d0-a0m)/(1d0+qsq/a2m)**2)**2
      gm=gm*7.78d0
      
      fe=(4d0*mp**2*ge+qsq*gm)/(4d0*mp**2+qsq)
      fm=gm

      ww1p=((1d0-x1)*(1d0-qsqmin/qsq)*fe)/qsq
      ww1pa=x1**2/2d0*fm/qsq

      ww1=((1d0-x1)*fe)/qsq**2/(1d0-x1)
      ww1=ww1/pi/(1d0-x1)/137d0

      ww1p=ww1p/pi/(1d0-x1)/137d0
      ww1pa=ww1pa/pi/(1d0-x1)/137d0

cccccccccc

      qsqmin=mp**2*x2**2/(1d0-x2)
      
      ge=(a0e/(1d0+qsqp/a1e)**2+(1d0-a0e)/(1d0+qsqp/a2e)**2)**2
      gm=(a0m/(1d0+qsqp/a1m)**2+(1d0-a0m)/(1d0+qsqp/a2m)**2)**2
      gm=gm*7.78d0
      
      fe=(4d0*mp**2*ge+qsqp*gm)/(4d0*mp**2+qsqp)
      fm=gm

      ww2p=((1d0-x2)*(1d0-qsqmin/qsqp)*fe)
     &/qsqp
      ww2pa=x2**2/2d0*fm/qsqp

      ww2=((1d0-x2)*fe)/qsqp**2/(1d0-x2)
      ww2=ww2/pi/(1d0-x2)/137d0

      ww2p=ww2p/pi/(1d0-x2)/137d0
      ww2pa=ww2pa/pi/(1d0-x2)/137d0

      if(io.eq.1)then
         out=dsqrt(ww1*ww2)
      elseif(io.eq.2)then
         out=dsqrt(ww1p*ww2pa+ww2p*ww1pa+ww1pa*ww2pa)
      endif

      return
      end
