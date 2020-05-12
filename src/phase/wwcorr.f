ccc   spin correlations in W+W- production (leptonic decays)
      subroutine wwcorr(wt,wtt)
      implicit double precision(a-y)
      double precision rhow(9),rhowc(9)
      integer mm,jj,nn
      complex*16 wt(10),rhoww(9)

      include 'ewpars.f'
      include 'pi.f'
      include 'partonmom2.f'
      include 'partonmom4.f'

      do mm=1,9
         rhoww(mm)=cdabs(wt(mm))**2
      enddo

      cost1=(p1(1)*paa(1)+p1(2)*paa(2)+p1(3)*paa(3))
     &/dsqrt((paa(1)**2+paa(2)**2+paa(3)**2)*
     &(p1(1)**2+p1(2)**2+p1(3)**2))
      sint1=dsqrt(1d0-cost1**2)

      cost2=(p2(1)*pbb(1)+p2(2)*pbb(2)+p2(3)*pbb(3))
     &/dsqrt((pbb(1)**2+pbb(2)**2+pbb(3)**2)*
     &(p2(1)**2+p2(2)**2+p2(3)**2))
      sint2=dsqrt(1d0-cost2**2)

      rhowp0=-sint1*dsqrt(3d0/2d0)
      rhowpm=(1d0+cost1)*dsqrt(3d0/4d0)
      rhowpp=(1d0-cost1)*dsqrt(3d0/4d0)
      rhowm0=-sint2*dsqrt(3d0/2d0)
      rhowmp=(1d0+cost2)*dsqrt(3d0/4d0)
      rhowmm=(1d0-cost2)*dsqrt(3d0/4d0)
      
      rhow(1)=rhowpp*rhowmp
      rhow(2)=rhowpp*rhowmm
      rhow(3)=rhowpm*rhowmp
      rhow(4)=rhowpm*rhowmm
      rhow(5)=rhowp0*rhowmp  
      rhow(6)=rhowp0*rhowmm      
      rhow(7)=rhowpp*rhowm0
      rhow(8)=rhowpm*rhowm0
      rhow(9)=rhowp0*rhowm0  

      wtt=0d0

      do mm=1,9  
            wtt=wtt+rhoww(mm)*rhow(mm)**2
      enddo

      return
      end
