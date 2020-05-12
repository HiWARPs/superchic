ccc   integrates bare + screened amplitude over k_t 
ccc   (two-photon induced processes)
      subroutine schimcgamion(p1x,p1y,p2x,p2y,outg)
      implicit double precision(a-y)
      implicit complex*16(z)
      integer jx,jy,i1,i2,p,i,nphi,nqt,jqt,jphi
      complex*16 out(4,10),x0(10),x00(10),outg(10)
      complex*16 screen(2,2)
      integer nni,nk
      complex*16 outt(0:200,0:10)

      include 'ppamp.f'
      include 'nchan.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'proc.f'
      include 'zi.f'
      include 'mandelstam.f'
      include 'pi.f'
      include 'nsurv.f'
      include 'inparticle.f'
      include 'zarr.f'
      include 'beam.f'
      
      nphi=s2int
      nqt=s2int*4

      qtmax=0.5d0
      xggmin=1d-2
      ypmax=dlog(xggmin**2+qtmax**2)
      ypmin=dlog(xggmin**2)

      nphi=20
      
      hphi=2d0*pi/dble(nphi)
      hqt=2d0/dble(nqt)

      nk=100
      hy=(ypmax-ypmin)/dble(nk)
      hqt=qtmax/dble(nk)
      
      zoutg=0d0
      
      do i1=0,nk
         do p=0,10
            outt(i1,p)=0d0
         enddo
      enddo

      do p=1,pol
         outg(p)=0d0
         do i=1,4
            out(i,p)=0d0
         enddo
      enddo

      call wtgengam
 
      if(sfac)then

         do jqt=0,nk

            qt=(dble(jqt))*hqt

            yp=ypmin+hy*dble(jqt)
            qt=dexp(yp)-xggmin**2
            tp2=qt
            qt=dsqrt(dabs(qt))
           
            sc=screeningionint(qt)
        
            do jphi=1,nphi
               
               phiq=(dble(jphi)-0.5d0)*hphi
               
               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)
               wt=hphi*hy*(tp2+xggmin**2)/2d0
               
               p1xp=p1x-tpx
               p1yp=p1y-tpy
               t12=p1xp**2+p1yp**2
               p2xp=tpx+p2x
               p2yp=tpy+p2y
               t22=p2xp**2+p2yp**2

               if(beam.eq.'ion')then
                  call formfacgamion(t12,t22,x00p)
               elseif(beam.eq.'ionp')then
                  call formfacgamionp(1,t12,t22,x00p)
                  call formfacgamionp(2,t12,t22,x00p2)
               endif
               
           do p=1,pol
                    
              x0(p)=x00p
              
              zout=-0.5d0*(ppa(p)+mma(p))*(p1xp*p2xp+p1yp*p2yp)
     &             -0.5d0*zi*(ppa(p)-mma(p))*(p1xp*p2yp-p2xp*p1yp)
     &             +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &             +zi*(p1xp*p2yp+p1yp*p2xp))*mpa(p)
     &             +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &             -zi*(p1xp*p2yp+p1yp*p2xp))*pma(p)
                     
              zout=zout*2d0
              
              outt(jqt,p)=outt(jqt,p)+x0(p)*wt
     &             *sc*zout

              if(p.eq.1)then
                 outt(jqt,0)=outt(jqt,0)+x00p2*wt*sc
              endif
              
           enddo
        enddo
     
      enddo

      do jqt=0,nk-2,2
         do p=1,pol
            outg(p)=outg(p)+(outt(jqt,p)+outt(jqt+1,p)*4d0
     &           +outt(jqt+2,p))/3d0
         enddo
         zoutg=zoutg+(outt(jqt,0)+outt(jqt+1,0)*4d0
     &        +outt(jqt+2,0))/3d0
      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      if(beam.eq.'ion')then
         call formfacgamion(t11,t22,x00p)
         x00p2=0d0
      elseif(beam.eq.'ionp')then
         call formfacgamionp(1,t11,t22,x00p)
         call formfacgamionp(2,t11,t22,x00p2)
      endif
         
      do p=1,pol
         zout=-0.5d0*(ppa(p)+mma(p))*(p1x*p2x+p1y*p2y)
     &        -0.5d0*zi*(ppa(p)-mma(p))*(p1x*p2y-p2x*p1y)
     &        +0.5d0*(p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x))*mpa(p)
     &        +0.5d0*(p1x*p2x-p1y*p2y-zi*(p1x*p2y+p1y*p2x))*pma(p)

         outg(p)=outg(p)+zout*x00p*2d0

         dbl=cdabs(outg(p))
         dbl=dsqrt(dbl**2+cdabs(zoutg+x00p2)**2*pincarr(p))
         outg(p)=dbl
         
      enddo

      return
      end
