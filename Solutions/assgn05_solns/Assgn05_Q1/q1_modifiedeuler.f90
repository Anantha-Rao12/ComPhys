! Name: Anantha Rao
! Reg no: 20181044
! This program provides the Improved Euler Method to solve the differential equation y' = f(x,y)

program modifiedeulermethod
      implicit none
      real*8:: x,y,y1, xwant, h, fxy, fxy1
      integer:: niter, i

      x= 0.0; y=0.0; h=0.001; xwant=1.57
      niter= int((xwant-x)/h)   
      open(unit=21, file="modified_euler_method.dat")

      do i=1,niter
      fxy = 1+y**2
      y1 = y + (h/2)*fxy
      fxy1 = 1 + y1**2
      y = y + h*fxy1

      write(21,*) dfloat(i)*h, y
      end do
end program modifiedeulermethod
