program composite_trapezoidal
    ! Author : Anantha S Rao (Reg no : 20181044)
    ! Program to evaluate the integral I = Int_a^b Gaus(0,1)    
    implicit none
    real*8,external::func
    ! An external function f in double precision

    real*8::a,b,h,pi,trapezoidal_integral
    ! a : Lower limit of the integral
    ! b : Upper limit of the integral
    ! h : bin size
    ! trapezoidal_integral : Value of the integral using trapezoidal technique

    integer::n,i
    ! n : number of bins 
    ! i : counter
    
    pi= 2.0d0*ASIN(1.0d0)
    n = 1000
    print 10, n
    10 format("Welcome to this program. &
         & This program does integration on I= Int_a^b Gaus(0,1) using trapezoidal rule with n=",i5," bins" )

    a = -3
    b = 3
    h=(b-a)/real(n)   ! size of each bin
    trapezoidal_integral=(h/2)*(func(a)+func(b))
    do i=1,n-1
        trapezoidal_integral=trapezoidal_integral+h*func(a+i*h)
    end do
    print *, "The integrated value using trapezoidal rule in (-3,3) is:",trapezoidal_integral

    a = -5
    b = 5
    h=(b-a)/real(n)   ! size of each bin
    trapezoidal_integral=(h/2)*(func(a)+func(b))
    do i=1,n-1
        trapezoidal_integral=trapezoidal_integral+h*func(a+i*h)
    end do
    print *, "The integrated value using trapezoidal rule in (-5,5) is:",trapezoidal_integral

end program composite_trapezoidal

real*8 function func(x)
    implicit none
    real*8::x, pi
    pi= 2.0d0*ASIN(1.0d0)
    func= (1.0d0/SQRT(2.0d0*pi))*EXP(-(x*x)/2.0d0)
end function
