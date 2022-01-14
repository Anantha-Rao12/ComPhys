program trapezoidal_with_different_n
    ! Author : Anantha Rao (REg no: 20181044)
    ! Program to compute the Integral f(x) with different number of binsizes in powers of 5 
      implicit none
      integer:: n,i
      ! n : number of bins
      ! i : counter

      real*8::a,b,h,func, fa,fb,trap_summ,log_error, pi
      ! a : Lower limit of the integral
      ! b : Upper limit of the integral
      ! h : Bin size
      ! trap_summ : Value of the integral using trapezoidal technique


      open(unit=1, file='20181044_trap_pi.dat')
      n = 1
      a = 0
      b = 1
      pi = 2.0d0*ASIN(1.0d0)
      ! write(*,*) "Give the value of n"
      ! read(*,*) n

      ! write(*,*) "Give the lower limit of the integral"
      ! read(*,*) a
      ! write(*,*) "Give the upper limit of the integral"
      ! read(*,*) b


      print *,"Welcome to this program.&
          & This program does integration on I=0^1 4.0/(1+x^2) using trapezoidal rule different number of bins"

      do
        n=n*5
        h=(b-a)/real(n)
        write(*,10) n
        10 format("Computing for n=",i10)

        fa = func(a)/2.0d0
        fb = func(b)/2.0d0

        trap_summ=0.0d0

        do i=1, n-1
            trap_summ = trap_summ+func(a+h*i)
        end do
      
      trap_summ=4.0d0*(trap_summ+fa+fb)*h
      log_error=LOG(ABS(pi - trap_summ))
      write(1,*) LOG(real(n)),"","",trap_summ,"",log_error
      if (n .ge. 1000000000) exit
      end do
      print*, "Done!, Output stored in trap_pi.dat"

END PROGRAM

real*8 function func(x)
      !evaluates f(x)=1/(1+x^2)
      implicit none
      real*8::x
      func=1.0d0/(1.0d0+x*x)
end function

