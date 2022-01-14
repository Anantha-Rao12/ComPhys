! Name: Anantha Rao
! Reg no: 20181044
! This program provides the Runge-Kutta4 Method to solve the differential equation y' = f(x,y)

program rungekutta
    implicit none
    real:: x,y,xwant,h,m1,m2,m3,m4,f
    integer:: niter,i

    x=0.0; y=0.0; xwant=1.57; h=0.01; 
    niter=int((xwant-x)/h)
    open(unit=21, file="q1_rk4.dat")

    do i=1,niter
        m1=h*f(x,y)
        m2=h*f(x+0.5*h,y+0.5*m1)
        m3=h*f(x+0.5*h,y+0.5*m2)
        m4=h*f(x+h,y+m3)
        x=x+h
        y=y+(m1+2.0*m2+2.0*m3+m4)*(1/6.0)
        write(21,*) dfloat(i)*h,y
    end do
end program

!define your function here
real function f(x,y)
real:: x,y
f=1+y**2

end function
