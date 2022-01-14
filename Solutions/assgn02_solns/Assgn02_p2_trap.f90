program trapezoidal_with_different_n
    ! Author : Anantha Rao (Reg no: 20181044)
    ! Program to compute the Integral f(x)=exp(x) with different number of bin sizes in powers of 10
      implicit none
      integer:: n,i
      ! n : number of bins
      ! i : counter

      real*8::a,b,h,func, fa,fb,trap_summ, error, actual_value
      ! a : Lower limit of the integral
      ! b : Upper limit of the integral
      ! h : Bin size
      ! trap_summ : Value of the integral using trapezoidal technique
      ! error : Value of the error between computed and analytical value
      ! actual_value : Analytical value of the integral


      open(unit=7, file='20181044_trap_e3.dat')
      n = 1
      a = 0
      b = 3
      actual_value = exp(3.0d0)-1
      ! write(*,*) "Give the value of n"
      ! read(*,*) n
      ! write(*,*) "Give the lower limit of the integral"
      ! read(*,*) a
      ! write(*,*) "Give the upper limit of the integral"
      ! read(*,*) b

      print *,"Welcome to this program.&
          & This program does integration on I=0^3 exp(x)dx using trapezoidal rule on different number of bins (n)"

      do
        h=(b-a)/real(n)
        write(*,10) n
        10 format("Computing for n=",i10)

        fa = func(a)/2.0d0
        fb = func(b)/2.0d0
        trap_summ=0.0d0

        do i=1, n-1
            trap_summ = trap_summ+func(a+h*i)
        end do
      
      trap_summ=(trap_summ+fa+fb)*h
      error=ABS(actual_value - trap_summ)
      write(7,*) n,"","",trap_summ,"",error
      if (n .ge. 1000000000) exit
      n=n*10
      end do
      print*, "Done!, Output stored in 20181044_trap_e^3.dat"

END PROGRAM

real*8 function func(x)
      !evaluates f(x)=e^x
      implicit none
      real*8::x
      func=exp(x)
end function

