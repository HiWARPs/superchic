      subroutine screencalc
      implicit double precision(a-y)
      integer i

      include 'scionpars.f'
      
      qmin=1d-3
      qmax=2d0

      lqmin=dlog(qmin)
      lqmax=dlog(qmax)

      itot=1900
    
      do i=1,itot+1

         lq=lqmin+(lqmax-lqmin)*dble(i-1)/dble(itot)
         qt=dexp(lq)

         scionarr(1,i)=lq
         scionarr(2,i)=screen(qt)
         
      enddo
      
      return
      end
