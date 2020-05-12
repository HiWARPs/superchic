ccc   integrates bare + screened amplitude over k_t 
ccc   (QCD induced processes)
      subroutine schimcion(qtmax,p1x,p1y,p2x,p2y,out)
      implicit double precision(a-y)
      integer jx,jy,i1,i2,p,nphi,nqt,jqt,jphi,nnphi,nk
      complex*16 out(10),x0(10),x00(10),x0p(10)
      complex*16 screen(2,2)
      complex*16 outt(0:4000,10)

      include 'nchan.f'
      include 'mt.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'nsurv.f'
      include 'pi.f'
      include 'zarr.f'
      include 'beam.f'
      include 'rho0.f'
      include 'ionqcd.f'
      
      
      call setmu(mu)

      do p=1,pol
         out(p)=0d0
      enddo
      
      nphi=20
      
      hphi=2d0*pi/dble(nphi)

      nk=200
      hqt=qtmax/dble(nk)
      
      do i1=0,nk
         do p=1,10
            outt(i1,p)=0d0
         enddo
      enddo

      call bare(mu,p1x,p1y,p2x,p2y,x0p)
      
      if(sfac)then
         
        do jqt=0,nk
            
            qt=hqt*dble(jqt)
            sc=screeningionint(qt)
            
           do jphi=1,nphi
             
               phiq=(dble(jphi)-0.5d0)*hphi
               
               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)

                wt=hphi*qt*hqt
               
                p1xp=p1x-tpx
                p1yp=p1y-tpy
                t12=p1xp**2+p1yp**2
                p2xp=tpx+p2x
                p2yp=tpy+p2y
                t22=p2xp**2+p2yp**2

                x00p=betaionex(-t12)*betaionex(-t22)

            if(beam.eq.'ion')then
               tp1=tpint(1,dsqrt(t12))+tpint(2,dsqrt(t12))
               tp2=tpint(1,dsqrt(t22))+tpint(2,dsqrt(t22))
               x00p=x00p*tp1*tp2
            elseif(beam.eq.'ionp')then
               tp1=tpint(1,dsqrt(t12))+tpint(2,dsqrt(t12))
               x00p=x00p*tp1
            endif
               
           do p=1,pol
              
              x0(p)=x0p(p)
              x0(p)=x0(p)*x00p
   
              outt(jqt,p)=outt(jqt,p)+x0(p)*wt*sc
              
           enddo
           
        enddo
      enddo
      
      do jqt=0,nk-2,2
         do p=1,pol
            out(p)=out(p)+(outt(jqt,p)+outt(jqt+1,p)*4d0
     &           +outt(jqt+2,p))/3d0
         enddo
      enddo
      
      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      x00p=betaionex(-t11)*betaionex(-t22)

      if(beam.eq.'ion')then
         tp1=tpint(1,dsqrt(t11))+tpint(2,dsqrt(t11))
         tp2=tpint(1,dsqrt(t22))+tpint(2,dsqrt(t22))
         if(ionqcd.eq.'coh')x00p=x00p*tp1*tp2
      elseif(beam.eq.'ionp')then
         tp1=tpint(1,dsqrt(t11))+tpint(2,dsqrt(t11))
         if(ionqcd.eq.'coh')x00p=x00p*tp1
      endif
      
      do p=1,pol
         out(p)=out(p)+x0p(p)*x00p
      enddo

      
      return
      end
