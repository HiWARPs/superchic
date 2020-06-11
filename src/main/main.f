ccc   calculates CEP cross section
      function cs(rarr,wgt)
      implicit double precision (a-z)
      complex*16 wt(10),wtn(10),wtd(10)
      double precision rarr(10),ran(5)
      integer i,p,icut
        
      include 'polvecs.f'
      include 'gencuts.f'
      include 'pi.f'
      include 'unweighted.f'
      include 'bin.f'
      include 'x.f'
      include 'zi.f'
      include 'vars.f'
      include 'survpars.f'
      include 'mom.f'
      include 'mandelstam.f'
      include 'range.f'
      include 'polarization.f'
      include 'proc.f'
      include 'mq.f'
      include 'norm.f'
      include 'meta.f'
      include 'mres.f'
      include 'decay.f'
      include 'mp.f'
      include 'quarkonia.f'
      include 'scorr.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'mpip.f'
      include 'mkp.f'
      include 'mphi.f'
      include 'gamma.f'
      include 'ewpars.f'
      include 'eff.f'
      include 'brs.f'
      include 'wmax.f'
      include 'unwsurv.f'
      include 'survin.f'
      include 'widths.f'
      include 'surv.f'
      include 'beam.f'
      include 'mion.f'
      include 'ionqcd.f'
      include 's2qcd.f'
      include 'ion.f'
      include 'rho0.f'
      include 'mcharg.f'

      wtt=0d0

      if(beam.eq.'ionp')call pAinit

      if(dps.eq.1)then
         mx=mres
         if(proc.eq.48)call bwrho(mx,jrho)
         if(fwidth)then
            if(proc.gt.20.and.proc.lt.38)then
               call bwchi(mx,jchi)
            endif
            if(proc.eq.68)call bwalp(mx,jalp)
            if(proc.eq.69.or.proc.eq.70)then
               rm=rarr(5)
               mx=mmin**2+(mmax**2-mmin**2)*rm
               mx=dsqrt(mx)
               call bwmono(mx,jmono)
               jmono=jmono*(mmax**2-mmin**2)
            endif
         endif
      else
         if(photo.or.gamma)then
         else
            if(mmin.lt.2d0)then
               print*,'WARNING: need mmin > 2 GeV for QCD processes'
               print*,'Please edit in input card'
               stop
            endif
         endif
         if(mmax.gt.rts)mmax=rts
         if(mmin.lt.2d0*mq)mmin=2d0*mq
         rm=rarr(5)
         mx=mmin+(mmax-mmin)*rm
         if(mmin.lt.1d-3)mmin=1d-3
         msub=1d0/1d0*(1d0/mmax**1d0+
     &           (1d0/mmin**1d0-1d0/mmax**1d0)*rm)
         mx=(1d0/msub/1d0)**(1d0/1d0)
      endif


      if(beam.eq.'prot')then
         mpp1=mp
         mpp2=mp
      elseif(beam.eq.'ion')then
         mpp1=mion
         mpp2=mion
      elseif(beam.eq.'ionp')then
         mpp1=mp
         mpp2=mion
      elseif(beam.eq.'el')then
         mpp1=me
         mpp2=me
      endif

      if(photo)then
         
         r1=rarr(2)
         r3=rarr(3)
         r4=rarr(4)
         r5=rarr(5)
         
         r2=ran2()
         
         ptmax=dsqrt(5d0)
         ptmin=0d0
         
         pt2sq=(ptmax-ptmin)*r1+ptmin
         pt2sq=pt2sq**2   
            
         phi1=2d0*pi*r2
         phi2=2d0*pi*r3+phi1
         
         rmx=dsqrt(pt2sq+mx**2)
         wtpt=2d0*dsqrt(pt2sq)*(ptmax-ptmin)
         
         xgmin=(mx/rts)**2
      
         if(r5.gt.0.5d0)then
            ypmax=dlog(xgmin**2*mpp1**2+ptmax**2)
            ypmin=dlog(xgmin**2*mpp1**2)
            yp=(ypmax-ypmin)*r4+ypmin
            pt1sq=dexp(yp)-xgmin**2*mpp1**2
            wt1=xgmin**2*mpp1**2+pt1sq
         else
            ypmax=dlog(xgmin**2*mpp2**2+ptmax**2)
            ypmin=dlog(xgmin**2*mpp2**2)
            yp=(ypmax-ypmin)*r4+ypmin
            pt1sq=pt2sq         
            pt2sq=dexp(yp)-xgmin**2*mpp2**2       
            wt1=xgmin**2*mpp2**2+pt2sq
         endif
         
         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)
         
         ptxx=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(ptxx+mx**2)
            
         if(beam.eq.'ionp')then
            ymax=dlog(rts/rmx)
            ymin=-ymax
         endif
            
         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry
            
         wty=ymax-ymin

         if(beam.eq.'ionp')r5=0d0
         
         if(r5.gt.0.5d0)then    ! photon emitted from q(1,k) 
            xgam=rmx*dexp(yx)/rts ! photon mom. fraction 
            wpsi=dsqrt(xgam*s)  ! proton-photon cms energy
            xglu=(rmx/wpsi)**2  ! gluon mom. fraction
            x1=xgam
            x2=xglu
            prot=1
         else                   ! photon emitted from q(2,k)
            xgam=rmx*dexp(-yx)/rts
            wpsi=dsqrt(xgam*s)  
            xglu=(rmx/wpsi)**2
            x2=xgam
            x1=xglu
            prot=2
         endif
         
         bpsi=bpsi0+4d0*alphapb*dlog(wpsi/w0b)

         if(x1.gt.1d0)goto 777
         if(x2.gt.1d0)goto 777
         
      elseif(gamma)then

         r2=rarr(2)
         r3=rarr(3)
         r4=rarr(4)
           
         r1=ran2()
         
         phi1=2d0*pi*r1
         phi2=2d0*pi*r2+phi1
         
         if(beam.eq.'prot'.or.beam.eq.'ion'.or.beam.eq.'ionp')then
            ptmax=dsqrt(3d0)       
         elseif(beam.eq.'el')then
            ptmax=dsqrt(50d0)       
         endif

         ptmin=0d0
            
         xgmin=(mx/rts)**2

         ypmax1=dlog(xgmin**2*mpp1**2+ptmax**2)
         ypmin1=dlog(xgmin**2*mpp1**2)
         
         ypmax2=dlog(xgmin**2*mpp2**2+ptmax**2)
         ypmin2=dlog(xgmin**2*mpp2**2)

         yp=(ypmax1-ypmin1)*r3+ypmin1
         ypp=(ypmax2-ypmin2)*r4+ypmin2
         
         pt1sq=dexp(yp)-xgmin**2*mpp1**2
         pt2sq=dexp(ypp)-xgmin**2*mpp2**2

         if(pt1sq.lt.0d0)then
            pt1sq=0d0
         endif
         
         if(pt2sq.lt.0d0)then
            pt2sq=0d0
         endif
         
         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)
         
         ptxx=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(ptxx+mx**2)
         

         if(beam.eq.'ionp')then
            ymax=dlog(rts/rmx)
            ymin=-ymax
         endif
            
         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry
         
         wty=ymax-ymin
         
         x1=rmx*dexp(yx)/rts    ! photon 1 mom. fraction 
         x2=rmx*dexp(-yx)/rts   ! photon 2 mom. fraction 

         if(x1.gt.1d0)goto 777
         if(x2.gt.1d0)goto 777

         qsq=(x1**2*mpp1**2+pt1sq)/(1d0-x1)
         qsqp=(x2**2*mpp2**2+pt2sq)/(1d0-x2)

         wpsi=dsqrt(x1*x2*s)
         
      else
         
         r1=rarr(2)
         r2=rarr(3)
         r4=rarr(4)
        
         r3=ran2()
         
         ptmax=dsqrt(2d0)
         ptmin=0d0
         
         pt1sq=r1*ptmax
         pt2sq=r2*ptmax

         pt1sq=ptmin+(ptmax-ptmin)*r1
         pt2sq=ptmin+(ptmax-ptmin)*r2
         
         pt1sq=pt1sq**2
         pt2sq=pt2sq**2
         
         phi1=2d0*pi*r3
         phi2=2d0*pi*r4+phi1
            
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)
         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         
         ptxsq=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(mx**2+ptxsq)

cccccccccccccccccccccc

         if(beam.eq.'ionp')then
            ymax=dlog(rts/rmx)
            ymin=-ymax
         endif

         wty=ymax-ymin
         
         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry
         
         x1=rmx/rts*dexp(yx)
         x2=rmx/rts*dexp(-yx)

         if(x1.gt.1d0)goto 777
         if(x2.gt.1d0)goto 777
 
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      aa1=(1d0-x1)*rts/dsqrt(2d0)
      aa2=(1d0-x2)*rts/dsqrt(2d0)
      cc1=0.5d0*(pt2sq+mpp2**2)
      cc2=0.5d0*(pt1sq+mpp1**2)

c     impose massive on-shell condition by solving
c                   p1+ + cc1/p2- = aa1
c                   p2- + cc2/p1+ = aa2 

      root1sq=(cc1-cc2-aa1*aa2)**2-4d0*cc2*aa1*aa2
      root2sq=(cc2-cc1-aa1*aa2)**2-4d0*cc1*aa1*aa2
      if(root1sq.le.0d0.or.root2sq.le.0d0)then
         weight=0d0
         goto 777
      endif
      p1p=(cc2-cc1+aa1*aa2+dsqrt(root1sq))/(2d0*aa2)
      p2m=(cc1-cc2+aa1*aa2+dsqrt(root2sq))/(2d0*aa1)
      p1m=(pt1sq+mpp1**2)/(2d0*p1p)
      p2p=(pt2sq+mpp2**2)/(2d0*p2m)

      q(1,3)=pt1x
      q(2,3)=pt1y
      q(3,3)=(p1p-p1m)/dsqrt(2d0)
      q(4,3)=(p1p+p1m)/dsqrt(2d0)

      q(1,4)=pt2x
      q(2,4)=pt2y
      q(3,4)=(p2p-p2m)/dsqrt(2d0)
      q(4,4)=(p2p+p2m)/dsqrt(2d0)

      do i=1,4
         q(i,5)=q(i,1)+q(i,2)-q(i,3)-q(i,4)
      enddo
      
ccccccccccccccccccccccccccccccccc

         if(dps.eq.2)then
            call twojetps(mx,mq,rarr(6),ps,uh,th)
         elseif(dps.eq.3)then
            call threejetps(mx,mq,ps)
         elseif(dps.eq.12)then
            if(proc.eq.19)then
               call twojetpsm(mx,mpsi,mpsip,ps,uh,th)
            else
               call twojetpsm(mx,meta,metap,ps,uh,th)
            endif
         endif

cccccccccccccccccc

cccc decays

ccccccccccccccccccc

         wt2=1d0
         wt3=1d0
         wt4=1d0
         wt6=1d0

         if(proc.eq.1.or.proc.eq.60.or.proc.eq.67)then
            call twobody(1,5,6,7,mb,mb,wt2)
         elseif(proc.eq.68.or.proc.eq.69.or.proc.eq.70)then
            call twobody(1,5,6,7,0d0,0d0,wt2)
         elseif(proc.eq.18.or.proc.eq.19.or.proc.eq.20)then
            call twobody(1,6,8,9,mmu,mmu,wt2a)
            call twobody(1,7,10,11,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.eq.21.or.proc.eq.22.or.proc.eq.23)then
            call twobody(1,5,6,7,0d0,mpsi,wt2a)
            call twobody(2,7,8,9,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.gt.23.and.proc.lt.29)then
            call twobody(1,5,6,7,m2b,m2b,wt2)
         elseif(proc.eq.29.or.proc.eq.30.or.proc.eq.31)then
            call fourbody(mpip,mpip,wt4)
         elseif(proc.eq.32.or.proc.eq.33.or.proc.eq.34)then
            call fourbody(mpip,mkp,wt4)
         elseif(proc.eq.35.or.proc.eq.36.or.proc.eq.37)then
            call sixbody(mpip,wt6)
         elseif(proc.eq.39.or.proc.eq.40.or.proc.eq.41)then
            call twobody(1,5,6,7,0d0,mups,wt2a)
            call twobody(2,7,8,9,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.gt.41.and.proc.lt.47)then
            call twobody(1,5,6,7,m2b,m2b,wt2)
         elseif(proc.eq.48)then
            call twobody(1,5,6,7,mpip,mpip,wt2)
         elseif(proc.eq.49)then
            call twobody(1,5,6,7,mkp,mkp,wt2)
         elseif(proc.eq.50.or.proc.eq.51.or.proc.eq.52)then
            call twobody(1,5,6,7,mmu,mmu,wt2)
         elseif(proc.eq.53)then
            call threebody(1,5,6,7,8,mpsi,mpip,mpip,wt3)
            call twobody(1,6,9,10,mmu,mmu,wt2a)
            wt3=wt3*wt2a
         elseif(proc.eq.54.or.proc.eq.61)then
            call twobodyw(6,8,9,0d0,mmu)
            call twobodyw(7,10,11,0d0,mmu)
         elseif(proc.eq.55.or.proc.eq.62)then
            call twobodyw(6,8,9,0d0,me)
            call twobodyw(7,10,11,0d0,me)
         elseif(proc.eq.73)then
            call threebody(1,6,8,9,10,mneut,mmu,0d0,wt3a)
            call threebody(1,7,11,12,13,mneut,mmu,0d0,wt3b)
            wt3=wt3a*wt3b
         elseif(proc.eq.74)then
            call threebody(1,6,8,9,10,mneut,0d0,0d0,wt3a)
            call threebody(1,7,11,12,13,mneut,0d0,0d0,wt3b)
            wt3=wt3a*wt3b
         elseif(proc.eq.75)then
            call threebody(1,6,8,9,10,mneut,mmu,0d0,wt3a)
            call threebody(2,7,11,12,13,mneut,0d0,0d0,wt3b)
            wt3=wt3a*wt3b
         elseif(proc.eq.76)then
            call twobody(1,6,8,9,mneut,mmu,wt2a)
            call twobody(1,7,10,11,mneut,mmu,wt2b)
            wt2=wt2a*wt2b
         endif

         if(beam.eq.'ionp')call pAboost
         
ccccccccccccccccccc  cuts ccccccccccccccccc

         
         neff0=neff0+1

         if(gencuts)then
            call cut(icut)
            if(icut.eq.0)goto 777
         endif

         neff=neff+1

ccccccccccccccccccccccccccccccccccccccccccc

          if(proc.eq.22.or.proc.eq.25.or.proc.eq.27.or.proc.eq.30
     &        .or.proc.eq.33.or.proc.eq.36)then
             call genpol1(5,echi1)
          elseif(proc.eq.23.or.proc.eq.26.or.proc.eq.28.or.proc.eq.
     &            31.or.proc.eq.34.or.proc.eq.37)then
             call genpol2
          endif

          if(proc.eq.40.or.proc.eq.43.or.proc.eq.45)then
             call genpol1(5,echi1)
          elseif(proc.eq.41.or.proc.eq.44.or.proc.eq.46)then
             call genpol2
          endif
          
ccccccccc

    

          
 456      if(photo)then
             if(beam.eq.'prot')then
                call schimcphot(pt1x,pt1y,pt2x,pt2y,wt)
             elseif(beam.eq.'ionp')then
                call schimcphotionp(pt1x,pt1y,pt2x,pt2y,wt)
             elseif(beam.eq.'ion')then
                print*,'Photoproduction not currently available for AA'
                stop
             endif
          elseif(gamma)then
             if(beam.eq.'prot'.or.beam.eq.'el')then
                call schimcgam(pt1x,pt1y,pt2x,pt2y,wt)
             elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
                call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt)
             endif
          else
             call wtgen
             if(beam.eq.'prot')then               
                call schimc(pt1x,pt1y,pt2x,pt2y,wt)
             elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
                if(ionqcd.eq.'incoh')then
                   if(sfac)then
                      call schimc(pt1x,pt1y,pt2x,pt2y,wt)
                      do p=1,pol
                         if(beam.eq.'ion')wt(p)=wt(p)*dsqrt(s2qcd)*an
                         if(beam.eq.'ionp')wt(p)=wt(p)*dsqrt(s2qcd*an)
                      enddo
                   else
                      ktcut=1.5d0
                      call schimcion(ktcut,pt1x,pt1y,pt2x,pt2y,wt)
                      do p=1,pol
                         if(beam.eq.'ion')wt(p)=wt(p)*dsqrt(s2qcd)*an
                         if(beam.eq.'ionp')wt(p)=wt(p)*dsqrt(s2qcd*an)
                      enddo
                   endif      
                elseif(ionqcd.eq.'coh')then
                   ptdif=dsqrt((pt1x-pt2x)**2+(pt1y-pt2y)**2)
                   ktcut=ptdif
                   if(ktcut.lt.1d0)ktcut=1.5d0
                   call schimcion(ktcut,pt1x,pt1y,pt2x,pt2y,wt)
                   if(sfac)then
                      sfac=.false.
                      call schimc(pt1x,pt1y,pt2x,pt2y,wtd)
                      sfac=.true.
                      call schimc(pt1x,pt1y,pt2x,pt2y,wtn)
                      do p=1,pol
                         wt(p)=wt(p)*cdabs(wtn(p))/cdabs(wtd(p))
                      enddo
                   endif
                endif
             endif
          endif

         wtt=0d0
         do p=1,pol
            wtt=wtt+cdabs(wt(p))**2
         enddo

         wtpol=1d0

         if(scorr)then
            if(proc.eq.18.or.proc.eq.19.or.proc.eq.20)then
               call jpsidecay(wt,wtt)
            endif
            if(proc.eq.21.or.proc.eq.39)then
               call chic0decay3(wtc0)
               wtt=wtt*wtc0
            endif
            if(proc.eq.22.or.proc.eq.40)then
               do i=4,6
                  wt(i)=conjg(wt(i-3))
               enddo
               call chic1decay3(wt,wtt)
            endif
           if(proc.eq.23.or.proc.eq.41)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay3(wt,wtt)
           endif
           if(proc.eq.25.or.proc.eq.43)then
              do i=4,6
                 wt(i)=conjg(wt(i-3))
              enddo
              call chic1decay2s(wt,wtt)
           endif
           if(proc.eq.26.or.proc.eq.44)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay2s(wt,wtt)
           endif
           if(proc.eq.27.or.proc.eq.45)then
              do i=4,6
                 wt(i)=conjg(wt(i-3))
              enddo
              call chic1decay2f(wt,wtt)
           endif
           if(proc.eq.28.or.proc.eq.46)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay2f(wt,wtt)
           endif
           if(proc.eq.50.or.proc.eq.51.or.proc.eq.52)then
              call jpsidecayphot(wtpol)
           endif
           if(proc.eq.54.or.proc.eq.55.or.proc.eq.61.or.proc.eq.62)then
              call wwcorr(wt,wtt)
           endif
         endif
    
         wtt=wtt*wt2*wt3*wt4*wt6

         if(decays)then
            do i=1,nbr
               wtt=wtt*br(i)
            enddo
         endif
         
         if(photo)then
            wtt=wtt*wty
            wtt=wtt*wtpt
            wtt=wtt*(ypmax-ypmin)*wt1
            wtt=wtt*(wpsi/w0)**delta*normp*bpsi
            if(scorr)wtt=wtt*wtpol
            if(proc.eq.48)wtt=wtt*jrho
         elseif(gamma)then
            wtt=wtt*(ypmax1-ypmin1)*(ypmax2-ypmin2)
            wtt=wtt*(xgmin**2*mpp1**2+pt1sq)*(xgmin**2*mpp2**2+pt2sq)
            wtt=wtt*wty
            wtt=wtt*2d0/mx
            if(dps.eq.1)wtt=wtt*pi/2d0/mx**3
            if(dps.eq.2)wtt=wtt*mx**2*(1d0/mmin**1d0-1d0/mmax**1d0)/1d0
            if(scorr)wtt=wtt*wtpol
            if(fwidth)then
               if(proc.eq.69.or.proc.eq.70)wtt=wtt*jmono
            endif
         else
            wtt=wtt*wty
            wtt=wtt*4d0*ptmax**2*dsqrt(pt1sq*pt2sq)*pi**2
            if(fwidth)then
               if(proc.gt.20.and.proc.lt.38)then
                  wtt=wtt*jchi
               elseif(proc.eq.68)then
                  wtt=wtt*jalp
               endif
            endif
         endif
         
         wtt=wtt/sym

         if(photo)goto 888
         if(gamma)goto 888

ccccccccccccc 1 body phase space

         if(dps.eq.1)then
            wtt=wtt/(16d0**2*pi**5)
         endif

         if(dps.gt.1)then
            wtt=wtt/(64d0*pi**2)*ps
            wtt=wtt*2d0*mx
            wtt=wtt/(16d0**2*pi**6)
            wtt=wtt*mx**2*(1d0/mmin**1d0-1d0/mmax**1d0)/1d0
         endif

         wtt=wtt*conv*surv


         
 888     val=wtt*wgt
         if(bin)then
         call binit(val)
         endif

         if(calcmax)then
            if(wmax.lt.wtt*wgt*ren)then
               wmax=wtt*wgt*ren
               iw=iw+1
            endif
         endif

         if(unw)then
            runw=ran2()
            call unweight(wtt*wgt*ren,runw)
         endif

 777     cs=wtt

         
      return
      end

   
