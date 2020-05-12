ccc   integrates bare + screened amplitude over k_t 
ccc   (photoproduction processes)
      subroutine schimcphotionp(p1x,p1y,p2x,p2y,out)
      implicit double precision(a-y)
      integer jx,jy,i1,i2,p,nphi,nqt,jqt,jphi
      complex*16 out(10),x0(10),x00(10),out1(10),out2(10),x01(10)
      complex*16 screen(2,2)

      include 'nchan.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'pi.f'
      include 'nsurv.f'
      include 'beam.f'

      do p=1,pol
         out(p)=0d0
         out1(p)=0d0
         out2(p)=0d0
      enddo

      nphi=s2int
      nqt=s2int*4

      hphi=2d0*pi/dble(nphi)
      hqt=2d0/dble(nqt)

      if(sfac)then

         do jqt=1,nqt
            
            qt=(dble(jqt)-0.5d0)*hqt
            
            tp2=qt**2
            
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
           
           do p=1,pol
              
              call formfacphotionp(1,t12,t22,x00p)
              call formfacphotionp(2,t12,t22,x00pp)

              x0(p)=x00p
              x01(p)=x00pp
              
              out(p)=out(p)+x0(p)*wt*sc*p2xp
              out1(p)=out1(p)+x0(p)*wt*sc*p2yp
              out2(p)=out2(p)+x01(p)*wt*sc

           enddo

      enddo
      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2
      
      call formfacphotionp(1,t11,t22,x00p)
      call formfacphotionp(2,t11,t22,x00pp)

      do p=1,pol
         out(p)=out(p)+x00p*p2x
         out1(p)=out1(p)+x00p*p2y
         out2(p)=out2(p)+x00pp
      enddo
      
      do p=1,pol
         out(p)=dsqrt(cdabs(out(p))**2+cdabs(out1(p))**2
     &        +cdabs(out2(p))**2)
      enddo

      return
      end
