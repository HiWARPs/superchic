ccc   interpolates opacity from array
      subroutine rhint(i,mxx,out)
      implicit double precision(a-y)
      integer it,i

      include 'RHS4.f'

      mx=mxx-0.3d0
      
      inc=5d-2
      it=nint(mx/inc)
      if(dble(it).gt.(mx/inc))then
         it=it-1
      endif

      m=RHarr(i+1,it+2)-RHarr(i+1,it+1)
      del=mxx-RHarr(1,it+1)
      out=m*del+RHarr(i+1,it+1)

c      write(6,*)it,m,mx,del,RHarr(1,it+1)
      
      return
      end

      subroutine RHread
      implicit double precision(a-y)
      integer i

      include 'RHS4.f'

c      open(10,file='../src/subamps/RHS4.dat')
      open(10,file='../src/subamps/Msrule.dat')

c     do i=1,991
      do i=1,1995
         read(10,*)RHarr(1,i),RHarr(2,i),RHarr(3,i)
c         RHarr(1,i)=dsqrt(RHarr(1,i))
      enddo

c$$$      call rhint(1,0.55d0,out)
c$$$      print*,out,RHarr(2,2)
      close(10)

      open(10,file='../amps/MQdFrw.dat')
      
      do i=1,1000
         read(10,*)Mpparr(1,4,i),Mpparr(2,4,i),Mpparr(3,4,i)
     &        ,Mpmarr(2,4,i),Mpmarr(3,4,i)
         Mpmarr(1,4,i)=Mpparr(1,4,i)
      enddo

      close(10)

      open(10,file='../amps/MQuFrw.dat')
      
      do i=1,1000
         read(10,*)Mpparr(1,5,i),Mpparr(2,5,i),Mpparr(3,5,i)
     &        ,Mpmarr(2,5,i),Mpmarr(3,5,i)
         Mpmarr(1,5,i)=Mpparr(1,5,i)
      enddo

      close(10)

      open(10,file='../amps/MQsFrw.dat')
      
      do i=1,1000
         read(10,*)Mpparr(1,6,i),Mpparr(2,6,i),Mpparr(3,6,i)
     &        ,Mpmarr(2,6,i),Mpmarr(3,6,i)
         Mpmarr(1,6,i)=Mpparr(1,6,i)
      enddo

      close(10)

      open(10,file='../amps/MQcFrw.dat')
      
      do i=1,1000
         read(10,*)Mpparr(1,7,i),Mpparr(2,7,i),Mpparr(3,7,i)
     &        ,Mpmarr(2,7,i),Mpmarr(3,7,i)
         Mpmarr(1,7,i)=Mpparr(1,7,i)
      enddo

      close(10)
      
      open(10,file='../amps/MQbFrw.dat')
      
      do i=1,1000
         read(10,*)Mpparr(1,8,i),Mpparr(2,8,i),Mpparr(3,8,i)
     &        ,Mpmarr(2,8,i),Mpmarr(3,8,i)
         Mpmarr(1,8,i)=Mpparr(1,8,i)
      enddo

      close(10)

      open(10,file='../amps/MQtFrw.dat')
      
      do i=1,1000
         read(10,*)Mpparr(1,9,i),Mpparr(2,9,i),Mpparr(3,9,i)
     &        ,Mpmarr(2,9,i),Mpmarr(3,9,i)
         Mpmarr(1,9,i)=Mpparr(1,9,i)
      enddo

      close(10)

      return
      end
      
ccc   interpolates opacity from array
      subroutine mint(i,k,j,mx,out)
      implicit double precision(a-y)
      integer it,i,j,k

      include 'RHS4.f'

      inc=1d-2
      it=nint(mx/inc)
      if(dble(it).gt.(mx/inc))then
         it=it-1
      endif

      if(i.eq.1)then
         m=Mpparr(k+1,j,it+2)-Mpparr(k+1,j,it+1)
         del=mx-Mpparr(1,j,it+1)
         out=m*del+Mpparr(k+1,j,it+1)
      elseif(i.eq.2)then
         m=Mpmarr(k+1,j,it+2)-Mpmarr(k+1,j,it+1)
         del=mx-Mpmarr(1,k+1,it+1)
         out=m*del+Mpmarr(k+1,j,it+1)
      endif
         
c      write(6,*)it,m,mx,del,Mpmarr(1,it+1)
      
      return
      end


      function fermi(x,x0,w)
      implicit double precision(a-y)

      fermi=1d0/(1d0+dexp((x-x0)/w))

      return
      end
