      function screeningionint(ri)
      implicit double precision(a-y)
      integer it

      include 'scionpars.f'
      
      rmin=scionarr(1,1)
      rmax=scionarr(1,itot)

      if(ri.lt.dexp(rmin))then
         screeningionint=scionarr(2,1)
         return
      endif
      
      r=dlog(ri)
      r=r-rmin
      
      if(ri.gt.dexp(rmax))then
         screeningionint=0d0
      else
         rbin=scionarr(1,2)-scionarr(1,1)
         it=nint(r/rbin)
         if(dble(it).gt.(r/rbin))then
            it=it-1
         endif
         m=(scionarr(2,it+2)-scionarr(2,it+1))
         m=m/(scionarr(1,it+2)-scionarr(1,it+1))
         del=r-scionarr(1,it+1)+rmin
         screeningionint=m*del+scionarr(2,it+1)

      endif

      return
      end
