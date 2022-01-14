! Name: Anantha Rao
! Reg no: 20181044
! This program provides the Euler Method to solve the differential equation y' = f(x,y)

program eulermethod
      implicit none
      real*8:: x,y, xwant, h, fxy
      integer:: niter, i

      x = 0.0; y=0.0; h=0.001; xwant=1.57 
      niter = int((xwant - x)/h)
      open(unit=21, file="eulermethod.dat")

      do i=1,niter
      fxy = 1+y**2
      x = x+h
      y = y + h*fxy

      write(21,*) dfloat(i)*h, y
      end do

end program eulermethod

!function subroutine
real function f(x,y)
          real*8:: x,y
          f = 1 + y**2
end function
