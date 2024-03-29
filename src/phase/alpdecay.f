ccc   generates rho meson invariant mass distribution according to modified BW 
      subroutine bwalp(mout,jalp)
      implicit double precision(a-y)

      include 'mres.f'
      include 'mpip.f'
      include 'gax.f'
      include 'pi.f'

ccccc calculate width

      width=gax**2*malp**3/64d0/pi

      width=0.5d0
      
ccccc      
      
      mminr=mres-10d0*width
      mmaxr=mres+10d0*width
      
      almin=datan(-(mres**2-mminr**2)/width/mres)
      almax=datan(-(mres**2-mmaxr**2)/width/mres)

      norm=almax-almin
      
      rm1=ran2()
      al1=almin+(almax-almin)*rm1
      mout=dsqrt(dtan(al1)*mres*width+mres**2)

      jalp=almax-almin
      jalp=jalp*width*mres
      jalp=jalp*(1d0+dtan(al1)**2)
 
cccccc up to this point just change of variables
ccccc
ccccc  now BW

      jalp=jalp*mout*width/norm
      jalp=jalp/((mres**2-mout**2)**2+mout**2*
     &     width**2)
      
      return
      end
