ccccc EPA form factors (proton)
      subroutine formfacgamion(t1,t2,out)
      implicit double precision(a-y)
      integer i1,i2,io

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'mom.f'
      include 'ion.f'
      include 'mn.f'
      include 'rho0.f'
      include 'mion.f'

      q0=0.71d0

c      mion=mp*az+(an-az)*mn
      
      x1i=x1
c      /an
      x2i=x2
c      /an
      
      qsq1=(x1i**2*mion**2+t1)/(1d0-x1i)
      qsq2=(x2i**2*mion**2+t2)/(1d0-x2i)    
      
ccccccccc

      f1=1d0/(1d0+qsq1/q0)**2
      f1=f1/(t1+x1i**2*mion**2)
      f1=f1*tpint(1,dsqrt(qsq1))
      f1=f1*dsqrt((1d0-x1i)/137d0/pi)
      
      f2=1d0/(1d0+qsq2/q0)**2
      f2=f2/(t2+x2i**2*mion**2)
      f2=f2*tpint(1,dsqrt(qsq2))      
      f2=f2*dsqrt((1d0-x2i)/137d0/pi)
      
      out=f1*f2
      
      return
      end
