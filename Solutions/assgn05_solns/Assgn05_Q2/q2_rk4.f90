!Name : Anantha Rao
!Reg no: 20181044

! Runge Kutta 4th order to solve set of first order differential equations

PROGRAM rk4
    IMPLICIT none

    ! Declarations
    ! t : time (initial value given)
    ! x : position 
    ! v : velocity (=dx/dt)
    ! h : step-size (=0.01)
    ! f : dx/dt 
    ! g : dv/dt
    ! E : energy

    real*8:: t,x,v,twant, h, x1,x2,x3, x4, v1, v2, v3, v4, f, g, E
    integer:: niter, i

    t=0.0; x=0.1; v=1.0; h=0.01; twant=50; niter=int((twant-t)/h)

    open(unit=50, file="q2_rk4.dat")
    E = (v**2)/2 - cos(x) 
    do i=1,niter
        write(50,*) dfloat(i)*h, x, v, E

        x1=h*f(t,x,v)
        v1=h*g(t,x,v)

        x2=h*f(t+0.5*h,x+0.5*x1, v)
        v2=h*g(t+0.5*h, x, v+0.5*v1)

        x3=h*f(t+0.5*h,x+0.5*x2, v)
        v3=h*g(t+0.5*h, x, v+0.5*v2)
        
        x4=h*f(t+h,x+x3, v)
        v4=h*g(t+h, x, v+v3)
        
        x=x+(x1+2.0*x2+2.0*x4+x4)*(1/6.0)
        v=v+(v1+ 2.0*v2+ 2.0*v3+ v4)*(1/6.0)
        ! E = (x**2+v**2)/2
        E = (v**2)/2 - cos(x)

    end do
end program



!define your function here
    real*8 function f(t,x,v)
    real*8:: t,x,v
    f=v
end function

!define your function here
    real*8 function g(t,x,v)
    real*8:: t,x,v, k,m
    k=1; m=1
    g=-sin(x)
end function

