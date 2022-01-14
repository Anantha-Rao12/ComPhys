program improvedeulermethod
    implicit none
      real*8:: x, y, y1, xwant, h, fxy, fxy1, fe
      integer:: niter, i

      x= 0.0; y=0.0; h=0.001; xwant=1.57
      niter = int((xwant -x)/h)
      open(unit=21, file="improved_euler_method.dat")

      do i=1,niter
      fxy = 1 + y**2
      y1 = y + h*fxy
      fe = 1 + y1**2
      fxy1 = (fxy+fe)/2
      y = y + h*fxy1

      write(21,*) dfloat(i)*h, y
      end do

end program improvedeulermethod
