program composite_trapezoidal
    ! Author : Anantha S Rao (Reg no : 20181044)
    ! Program to evaluate the integral I = Int_0^pi sin(x)
    implicit none
    real*8,external::func
    ! An external function f in double precision

    real*8::a,b,h,trapezoidal_integral
    ! a : Lower limit of the integral
    ! b : Upper limit of the integral
    ! h : bin size
    ! trapezoidal_integral : Value of the integral using trapezoidal technique

    integer::n,i
    ! n : number of bins 
    ! i : counter
    
    n = 1000
    a = 0
    b = 2.0d0*ASIN(1.0d0)
    print 10, n
    10 format("Welcome to this program. &
         & This program does integration on I= Int_0^pi sin(x) using trapezoidal rule with n=",i5," bins" )

    h=(b-a)/real(n)   ! size of each bin
    trapezoidal_integral=(h/2)*(func(a)+func(b))

    do i=1,n-1
        trapezoidal_integral=trapezoidal_integral+h*func(a+i*h)
    end do
    print *, "The integrated value using trapezoidal rule is:",trapezoidal_integral
    !5 format("The integrated value using trapezoidal rule is:",
end program composite_trapezoidal

real*8 function func(x)
    implicit none
    real*8::x
    func=SIN(x)
end function
