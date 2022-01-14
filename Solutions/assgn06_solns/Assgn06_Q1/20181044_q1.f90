! Program to simulate a string of 50 beads/particles with given initial conditions and Periodic boundry conditions
! Author: Anantha Rao
! Reg no:20181044

program coupled_DE
    implicit none
    integer,parameter:: n=50, niter=2000
    integer:: i, j, a(n), b(n)
    real*8:: t=0.0d0, dt=0.02d0, y_old(n), y_new(n), v_old(n), v_new(n)
    real*8:: f, f0(n), f1(n), f2(n), f3(n), g0(n), g1(n), g2(n), g3(n)

    open(42, file="coupled_DE_positions.csv", action="write")
    y_old=0.0d0; v_old=0.0d0
    do i=1,n
        if(i==1 .or. i==26)  y_old(i)=0.8d0
        write(42,*)t,",",i,",",y_old(i)
    end do
        write(42,*)"--------------------------------------------------------------------"

    do j=1,niter
            t=t+dt
        do i=1,n !RK4 step 1
            a(i)=i-1; b(i)=i+1
            if(i==1) a(i)=50; if(i==50) b(i)=1  !PBC
            f0(i)=v_old(i)
            g0(i)=f(y_old(a(i)),y_old(i),y_old(b(i)))
            y_new(i)= y_old(i)+f0(i)*dt/2.0d0
            v_new(i)=v_old(i)+g0(i)*dt/2.0d0
        end do

        do i=1,n !RK4 step 2
            f1(i)=v_new(i)
            g1(i)=f(y_new(a(i)),y_new(i),y_new(b(i)))
            y_new(i)= y_old(i)+f1(i)*dt/2.0d0
            v_new(i)=v_old(i)+g1(i)*dt/2.0d0
        end do
        
        do i=1,n !RK4 step 3
            f2(i)=v_new(i)
            g2(i)=f(y_new(a(i)),y_new(i),y_new(b(i)))
            y_new(i)= y_old(i)+f2(i)*dt
            v_new(i)=v_old(i)+g2(i)*dt
        end do
        
        do i=1,n !RK4 step 4
            f3(i)=v_new(i)
            g3(i)=f(y_new(a(i)),y_new(i),y_new(b(i)))
            y_new(i)=y_old(i)+ dt/6.0d0 *(f0(i)+2.0d0*f1(i)+2.0d0*f2(i)+f3(i))
            v_new(i)=v_old(i)+ dt/6.0d0 *(g0(i)+2.0d0*g1(i)+2.0d0*g2(i)+g3(i))
            write(42,*)t,",",i,",",y_new(i)
        end do
            write(42,*)"---------------------------------------------------------------------"
         
         y_old=y_new
         v_old=v_new
    end do
    close(89)
end program coupled_DE

real*8 function f(a,b,c) result(reslt)
    implicit none
    real*8::a,b,c
    reslt=a+c-2.0d0*b
end function
