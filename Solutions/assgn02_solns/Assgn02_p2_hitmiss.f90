program hitmiss_mc_int
      ! Program to evaluate I= 0^3 exp(3)dx

      implicit none
      integer:: counter, nbins, npts_in_curve
      ! counter : counter variable
      ! nbins : number of bins
      ! npts_in_curve : number of points that lie within the region of
      ! interest

      real*8:: x, p, y, integral, actual_value
      ! x : varialbe storing random number of x-axis
      ! y : varialbe storing random number of y-axis
      ! p : dummy variable for the random number generator
      ! integral : value of the computed integral
      ! actual_value : Analytic value of the integral

     print *,"Welcome to this program.&
         & This program does integration on I=0^3 exp(x)dx using Monte-Carlo Hit & Miss method on different number of bins (n) " 
      
      nbins=1
      actual_value=exp(3.0d0)-1
      7 npts_in_curve=0
      open(unit=1, file="20181044_exp_accep_rej.dat")
      write(*,10) nbins
      10 format("Computing for nbins=", i10)
 
      do counter=1,nbins
          ! find random value of "x" between 0 and 3
          call random_number(p)
          x=3.0d0*p

          ! find value of function between 0 and exp(3)
          call random_number(p)
          y=exp(3.0)*p

          ! If y at exp(x) is below the curve we accept the number
          if (y .lt. exp(x)) then
              npts_in_curve=npts_in_curve+1
          end if
      end do

      !multiply with area of rectangle and divide by the no. of cycles
      integral=3.0d0*exp(3.0)*(real(npts_in_curve)/real(nbins))
      write(1,*) nbins,"",integral,"",abs(integral-actual_value)

      ! automated
      nbins=nbins*10
      if (nbins .le. 100000000 ) goto 7
      print*, "Done! Output stored in 20181044_exp_accep_rej.dat"
      end program hitmiss_mc_int
