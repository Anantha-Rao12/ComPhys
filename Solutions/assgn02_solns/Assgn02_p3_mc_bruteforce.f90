program multi_dimensional_bruteforce_mc

      !Author : Anantha Rao
      ! Reg no : 20181044
      ! 6D (x1,x2,x3,x4,x5,x6) Integral
      implicit none
      integer:: n,i,j
      real*8 :: x(6), y, func, int_mc, var, sigma, p(6)
      ! p : Dummy argument for random number generator (6 values)
      real*8 :: length, volume

      length=5.0d0 ! use one side of each dimension (0,5)
      volume=(2.0d0*length)**6  ! volume of 6D space

      open(unit=1, file="20181044_multi_brute_mc.dat")
      n=1

      ! Initiate variables


print *,"Welcome to this program.& 
& This program solves mutidimensional integral using Monte-Carlo rule"

    7 int_mc=0.0d0
      var=0.0d0
      sigma=0.0d0
      write(*,10) n
      10 format("Computig for n=",i10)

      do i=1,n
        call random_number(p)
        do j=1,6
            x(j) = -length + 2.0d0*length*p(j)
        end do
        int_mc=int_mc+func(x)
        sigma=sigma+func(x)*func(x)
      end do

      int_mc = int_mc/real(n)
      sigma=sigma/real(n)
      var=sigma-int_mc*int_mc

      int_mc = volume*int_mc
      sigma=volume*sqrt(var/real(n))

      write(1,*) n, "", int_mc, "", sigma

      ! begin automation
      n=n*10
      if (n .lt. 1000000000 ) goto 7
      print*,"Done! Output stored in file 20181044_multi_brute_mc.dat"
      ! end automation


end program


real*8 function func(x)
      implicit none
      real*8::x(6), xx, yy, xy
      real*8::a,b
      a=1.0d0
      b=0.5d0

      xx= x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
      yy= x(4)*x(4) + x(5)*x(5) + x(6)*x(6)
      xy= (x(1)-x(4))**2 + (x(2)-x(5))**2 + (x(3)-x(6))**2
      func=exp(-a*xx - a*yy - b*xy)

end function
