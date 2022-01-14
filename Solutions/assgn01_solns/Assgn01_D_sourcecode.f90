subroutine trapezoidal (b, trap_sum)
    ! Author : Anantha Rao (Reg no : 20181044)
    ! Program to compute the value of erf(x=b)
    ! Erf(x) = 2/sqrt(pi)*Int_0^x e^{-x^2}
    implicit none
    integer :: i,n
    ! n : number of bins
    ! i : counter

    REAL*8:: a,b,h, fa, fb, trap_sum, func
    ! a : Lower limit of Error function
    ! b : Upper limit of Error function
    ! fa : erf(x=a)
    ! fb : erf(x=b)
    ! trap_sum : value of the integral
    ! func : external function

    real*8, parameter :: pi = 2.0d0*asin(1.0d0)
    n = 10000   
    a = 0.0d0
    h = (b-a)/real(n)
    fa = func(a)/2.0d0
    fb = func(b)/2.0d0
    trap_sum = 0.0d0

    do i=1, n-1
      trap_sum=trap_sum+func(a+h*i)
    end do
    trap_sum=(trap_sum+fa+fb)*h
end subroutine

real*8 function func(x)
    ! Evaluates exp(-x^2)
    implicit none
    real*8::x
    real*8, parameter :: pi = 2.0d0*asin(1.0d0)    
    func=exp(-(x*x))
end function

program main
    implicit none
    real*8 :: b, result
    ! b : upper limit of the integral
    real*8, parameter :: pi = 2.0d0*asin(1.0d0) 
    print*,"Computing Error function values..."
    b=-3.0d0
    5   b=b + 0.005d0
        print*, "Computing Erf(x) for x=",b
    call trapezoidal(b, result)
    open(1,file = '20181044_D_erfvalues.dat')
    write(1,*) b, 2.0d0*result/(sqrt(pi))
    if (b .lt. 3.0d0) goto 5
    print*,"Output written to D_erfvalues.dat"
end program main

