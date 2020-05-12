      function opacpbint(r)
      implicit double precision(a-y)
      integer it
      
      include 'opacpbpars.f'
      
      rmax=opacpbarr(1,ioppb)
      rmin=opacpbarr(1,1)

      if(r.gt.rmax)then
         opacpbint=1d0
      elseif(r.lt.rmin)then
         opacpbint=0d0
      else
         rbin=opacpbarr(1,2)-opacpbarr(1,1)
         it=nint((r-rmin)/rbin)
         if(dble(it).gt.((r-rmin)/rbin))then
            it=it-1
         endif
         m=(opacpbarr(2,it+2)-opacpbarr(2,it+1))
         m=m/(opacpbarr(1,it+2)-opacpbarr(1,it+1))
         del=r-opacpbarr(1,it+1)
         opacpbint=m*del+opacpbarr(2,it+1)
      endif
         
      return
      end
