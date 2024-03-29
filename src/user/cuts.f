      subroutine cut(icut)
      implicit double precision(a-y)
      double precision p1(4),p2(4),p3(4)
      integer icut,jflag,i
      logical accut

      include 'gencuts.f'
      include 'vars.f'
      include 'mom.f'
      include 'proc.f'
      include 'jetalg.f'
      include 'decay.f'
      include 'pi.f'
      include 'ptXcuts.f'
      include 'range.f'
      include 'interpolate.f'
      include 'tetraquark.f'

      icut=0

ccccccccccccccccccccccccccccccccccccccccccccc
ccc
cc    Place user-defined cuts here if needed
cc
cc    if(..)return - return for failed cut
cc
cccccccccccccccccccccccccccccccccccccccccccccc

      if(yx.gt.ymax)return
      if(yx.lt.ymin)return
      if(dsqrt(q(1,5)**2+q(2,5)**2).gt.ptxmax)return

      if(run3)then
         if(q(4,6).lt.1.5d0)return
         if(q(4,7).lt.1.5d0)return
      endif
         
      if(decay4)then

         if(proc.eq.53)then
            
            et1=dsqrt(q(1,7)**2+q(2,7)**2)
            et2=dsqrt(q(1,8)**2+q(2,8)**2)
            et3=dsqrt(q(1,9)**2+q(2,9)**2)
            et4=dsqrt(q(1,10)**2+q(2,10)**2)
            
            pmod1=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
            pmod2=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
            pmod3=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod4=dsqrt(q(1,10)**2+q(2,10)**2+q(3,10)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,7))/(pmod1-q(3,7)))
            eta2=0.5d0*dlog((pmod2+q(3,8))/(pmod2-q(3,8)))
            eta3=0.5d0*dlog((pmod3+q(3,9))/(pmod3-q(3,9)))
            eta4=0.5d0*dlog((pmod4+q(3,10))/(pmod4-q(3,10)))

         elseif(proc.eq.18.or.proc.eq.19.or.proc.eq.20)then

            et1=dsqrt(q(1,8)**2+q(2,8)**2)
            et2=dsqrt(q(1,9)**2+q(2,9)**2)
            et3=dsqrt(q(1,10)**2+q(2,10)**2)
            et4=dsqrt(q(1,11)**2+q(2,11)**2)
            
            pmod1=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
            pmod2=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod3=dsqrt(q(1,10)**2+q(2,10)**2+q(3,10)**2)
            pmod4=dsqrt(q(1,11)**2+q(2,11)**2+q(3,11)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,8))/(pmod1-q(3,8)))
            eta2=0.5d0*dlog((pmod2+q(3,9))/(pmod2-q(3,9)))
            eta3=0.5d0*dlog((pmod3+q(3,10))/(pmod3-q(3,10)))
            eta4=0.5d0*dlog((pmod4+q(3,11))/(pmod4-q(3,11)))

            elseif(proc.gt.28.and.proc.lt.35)then

            et1=dsqrt(q(1,6)**2+q(2,6)**2)
            et2=dsqrt(q(1,7)**2+q(2,7)**2)
            et3=dsqrt(q(1,8)**2+q(2,8)**2)
            et4=dsqrt(q(1,9)**2+q(2,9)**2)
            
            pmod1=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
            pmod2=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
            pmod3=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
            pmod4=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,6))/(pmod1-q(3,6)))
            eta2=0.5d0*dlog((pmod2+q(3,7))/(pmod2-q(3,7)))
            eta3=0.5d0*dlog((pmod3+q(3,8))/(pmod3-q(3,8)))
            eta4=0.5d0*dlog((pmod4+q(3,9))/(pmod4-q(3,9)))

         elseif(proc.eq.74)then

            et1=dsqrt(q(1,9)**2+q(2,9)**2)
            et2=dsqrt(q(1,10)**2+q(2,10)**2)
            et3=dsqrt(q(1,12)**2+q(2,12)**2)
            et4=dsqrt(q(1,13)**2+q(2,13)**2)
            
            pmod1=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod2=dsqrt(q(1,10)**2+q(2,10)**2+q(3,10)**2)
            pmod3=dsqrt(q(1,12)**2+q(2,12)**2+q(3,12)**2)
            pmod4=dsqrt(q(1,13)**2+q(2,13)**2+q(3,13)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,9))/(pmod1-q(3,9)))
            eta2=0.5d0*dlog((pmod2+q(3,10))/(pmod2-q(3,10)))
            eta3=0.5d0*dlog((pmod3+q(3,12))/(pmod3-q(3,12)))
            eta4=0.5d0*dlog((pmod4+q(3,13))/(pmod4-q(3,13)))
            
         endif
         
         if(et1.lt.ptamin4)return
         if(et2.lt.ptbmin4)return
         if(et3.lt.ptcmin4)return
         if(et4.lt.ptdmin4)return
         if(eta1.gt.etaamax4)return
         if(eta2.gt.etabmax4)return
         if(eta3.gt.etacmax4)return
         if(eta4.gt.etadmax4)return
         if(eta1.lt.etaamin4)return
         if(eta2.lt.etabmin4)return
         if(eta3.lt.etacmin4)return
         if(eta4.lt.etadmin4)return

         
      elseif(decay6)then
         
         et1=dsqrt(q(1,6)**2+q(2,6)**2)
         et2=dsqrt(q(1,7)**2+q(2,7)**2)
         et3=dsqrt(q(1,8)**2+q(2,8)**2)
         et4=dsqrt(q(1,9)**2+q(2,9)**2)
         et5=dsqrt(q(1,10)**2+q(2,10)**2)
         et6=dsqrt(q(1,11)**2+q(2,11)**2)
         
         pmod1=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
         pmod2=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
         pmod3=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
         pmod4=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
         pmod5=dsqrt(q(1,10)**2+q(2,10)**2+q(3,10)**2)
         pmod6=dsqrt(q(1,11)**2+q(2,11)**2+q(3,11)**2)
            
         eta1=0.5d0*dlog((pmod1+q(3,6))/(pmod1-q(3,6)))
         eta2=0.5d0*dlog((pmod2+q(3,7))/(pmod2-q(3,7)))
         eta3=0.5d0*dlog((pmod3+q(3,8))/(pmod3-q(3,8)))
         eta4=0.5d0*dlog((pmod4+q(3,9))/(pmod4-q(3,9)))
         eta5=0.5d0*dlog((pmod5+q(3,10))/(pmod5-q(3,10)))
         eta6=0.5d0*dlog((pmod6+q(3,11))/(pmod6-q(3,11)))

         if(et1.lt.ptamin6)return
         if(et2.lt.ptbmin6)return
         if(et3.lt.ptcmin6)return
         if(et4.lt.ptdmin6)return
         if(et5.lt.ptemin6)return
         if(et6.lt.ptfmin6)return
         if(eta1.gt.etaamax6)return
         if(eta2.gt.etabmax6)return
         if(eta3.gt.etacmax6)return
         if(eta4.gt.etadmax6)return
         if(eta5.gt.etaemax6)return
         if(eta6.gt.etafmax6)return
         if(eta1.lt.etaamin6)return
         if(eta2.lt.etabmin6)return
         if(eta3.lt.etacmin6)return
         if(eta4.lt.etadmin6)return
         if(eta5.lt.etaemin6)return
         if(eta6.lt.etafmin6)return

      elseif(decay3)then

         if(proc.eq.75)then

            et1=dsqrt(q(1,9)**2+q(2,9)**2)
            et2=dsqrt(q(1,12)**2+q(2,12)**2)
            et3=dsqrt(q(1,13)**2+q(2,13)**2)
            
            pmod1=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod2=dsqrt(q(1,12)**2+q(2,12)**2+q(3,12)**2)
            pmod3=dsqrt(q(1,13)**2+q(2,13)**2+q(3,13)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,9))/(pmod1-q(3,9)))
            eta2=0.5d0*dlog((pmod2+q(3,12))/(pmod2-q(3,12)))
            eta3=0.5d0*dlog((pmod3+q(3,13))/(pmod3-q(3,13)))

            
         else
            
            et1=dsqrt(q(1,6)**2+q(2,6)**2)
            et2=dsqrt(q(1,8)**2+q(2,8)**2)
            et3=dsqrt(q(1,9)**2+q(2,9)**2)
            
            pmod1=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
            pmod2=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
            pmod3=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,6))/(pmod1-q(3,6)))
            eta2=0.5d0*dlog((pmod2+q(3,8))/(pmod2-q(3,8)))
            eta3=0.5d0*dlog((pmod3+q(3,9))/(pmod3-q(3,9)))

      endif
      
      if(et1.lt.ptamin3)return
      if(et2.lt.ptbmin3)return
      if(et3.lt.ptcmin3)return
      if(eta1.gt.etaamax3)return
      if(eta2.gt.etabmax3)return
      if(eta3.gt.etacmax3)return
      if(eta1.lt.etaamin3)return
      if(eta2.lt.etabmin3)return
      if(eta3.lt.etacmin3)return
           
         
      elseif(dps.eq.2.or.dps.eq.12.or.decay2)then

      delphi=(q(1,6)*q(1,7)+q(2,6)*q(2,7))/dsqrt((q(1,6)**2+q(2,6)**2))
     &     /dsqrt((q(1,7)**2+q(2,7)**2))
      if(dabs(delphi).gt.1d0)delphi=delphi/dabs(delphi)
      delphi=dacos(delphi)
      acoab=1d0-delphi/pi
      if(acoab.gt.acoabmax)return 
      
         
         if(proc.eq.54.or.proc.eq.55.or.proc.eq.76)then
            pmod6=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod7=dsqrt(q(1,11)**2+q(2,11)**2+q(3,11)**2)
            eta6=0.5d0*dlog((pmod6+q(3,9))/(pmod6-q(3,9)))
            eta7=0.5d0*dlog((pmod7+q(3,11))/(pmod7-q(3,11)))
            y6=0.5d0*dlog((q(4,9)+q(3,9))/(q(4,9)-q(3,9)))
            y7=0.5d0*dlog((q(4,11)+q(3,11))/(q(4,11)-q(3,11)))
            pt6=dsqrt(q(1,9)**2+q(2,9)**2)
            pt7=dsqrt(q(1,11)**2+q(2,11)**2)
         elseif(proc.eq.73)then
            pmod6=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod7=dsqrt(q(1,12)**2+q(2,12)**2+q(3,12)**2)
            eta6=0.5d0*dlog((pmod6+q(3,9))/(pmod6-q(3,9)))
            eta7=0.5d0*dlog((pmod7+q(3,12))/(pmod7-q(3,12)))
            y6=0.5d0*dlog((q(4,9)+q(3,9))/(q(4,9)-q(3,9)))
            y7=0.5d0*dlog((q(4,12)+q(3,12))/(q(4,12)-q(3,12)))
            pt6=dsqrt(q(1,9)**2+q(2,9)**2)
            pt7=dsqrt(q(1,12)**2+q(2,12)**2)
         else
            pmod6=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
            pmod7=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
            eta6=0.5d0*dlog((pmod6+q(3,6))/(pmod6-q(3,6)))
            eta7=0.5d0*dlog((pmod7+q(3,7))/(pmod7-q(3,7)))
            y6=0.5d0*dlog((q(4,6)+q(3,6))/(q(4,6)-q(3,6)))
            y7=0.5d0*dlog((q(4,7)+q(3,7))/(q(4,7)-q(3,7)))
            pt6=dsqrt(q(1,6)**2+q(2,6)**2)
            pt7=dsqrt(q(1,7)**2+q(2,7)**2)
         endif

      if(pt6.lt.ptamin)return
      if(pt7.lt.ptbmin)return
      if(eta6.lt.etaamin)return
      if(eta7.lt.etabmin)return
      if(eta6.gt.etaamax)return
      if(eta7.gt.etabmax)return

      elseif(dps.eq.3)then

      et1=dsqrt(q(1,6)**2+q(2,6)**2)
      et2=dsqrt(q(1,7)**2+q(2,7)**2)
      et3=dsqrt(q(1,8)**2+q(2,8)**2)

      pmod1=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
      pmod2=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
      pmod3=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)

      eta1=0.5d0*dlog((pmod1+q(3,6))/(pmod1-q(3,6)))
      eta2=0.5d0*dlog((pmod2+q(3,7))/(pmod2-q(3,7)))
      eta3=0.5d0*dlog((pmod3+q(3,8))/(pmod3-q(3,8)))

      if(et1.lt.ptamin3)return
      if(et2.lt.ptbmin3)return
      if(et3.lt.ptcmin3)return
      if(eta1.gt.etaamax3)return
      if(eta2.gt.etabmax3)return
      if(eta3.gt.etacmax3)return
      if(eta1.lt.etaamin3)return
      if(eta2.lt.etabmin3)return
      if(eta3.lt.etacmin3)return
     
cccccccccc kt/anti-kt/durham alg

      if(jalg.eq.'kt')then

      d1=q(1,6)**2+q(2,6)**2
      dmin1=d1
      ijet=1
      d2=q(1,7)**2+q(2,7)**2
      if(d2.lt.dmin1)ijet=2
      dmin1=min(d2,dmin1)
      d3=q(1,8)**2+q(2,8)**2
      if(d3.lt.dmin1)ijet=3
      dmin1=min(d3,dmin1)

      elseif(jalg.eq.'antikt')then
         
      d1=1d0/(q(1,6)**2+q(2,6)**2)
      ijet=1   
      dmin1=d1
      d2=1d0/(q(1,7)**2+q(2,7)**2)
      if(d2.lt.dmin1)ijet=2
      dmin1=min(d2,dmin1)
      d3=1d0/(q(1,8)**2+q(2,8)**2)
      if(d3.lt.dmin1)ijet=3
      dmin1=min(d3,dmin1)

      endif

ccccccccc  

      d12=min(d1,d2)*((eta1-eta2)**2+dphi(6,7)**2)/rjet**2
      dmin2=d12
      d13=min(d1,d3)*((eta1-eta3)**2+dphi(6,8)**2)/rjet**2
      dmin2=min(d13,dmin2)
      d23=min(d2,d3)*((eta2-eta3)**2+dphi(7,8)**2)/rjet**2
      dmin2=min(d23,dmin2)

      if(dmin1.gt.dmin2)then
         jflag=2
      else
         if(ijet.eq.1)dmin1p=min(d2,d3)
         if(ijet.eq.2)dmin1p=min(d1,d3)
         if(ijet.eq.3)dmin1p=min(d1,d2)
         if(ijet.eq.1)dmin2p=d23
         if(ijet.eq.2)dmin2p=d13
         if(ijet.eq.3)dmin2p=d12

         jflag=1
         
         if(dmin1p.gt.dmin2p)jflag=2
      endif

      if(jflag.eq.2)return
           
      endif

       icut=1

      return
      end
 
