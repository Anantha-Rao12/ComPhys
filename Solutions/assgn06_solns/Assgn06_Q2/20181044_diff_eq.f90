! Program to solve a second order DE using Gauss Seidel method
! Author: Anantha Rao
! Reg no: 20181044

program diff_eq
implicit none

real*8, parameter:: f_x=0.0d0, end_x=1.0d0, dx=0.01 ! INTERVAL
real*8 :: f_y, end_y
integer, parameter:: nop=int((end_x-f_x)/dx) ! NO OF GRID POINTS IN THE INTERVAL
integer :: i,j,k, ll,cond

! Pre factors for calculation 
real*8 :: x(nop), y(nop), y_old(nop), limit
real*8, parameter :: dum1=1.0/(2.0 - 10.0d0*dx*dx)
real*8, parameter :: dum2 = (1.0d0-2.5*dx)
real*8, parameter :: dum3 = (1.0d0+2.5*dx)
      real*8, parameter :: dum4 = -10.0d0*dx*dx

! Solving differential equation: y'' -5y' + 10y = 10x
! y'(i) = (y(i+1) - y(i-1))/(2*dx)
! y'' = (y(i+1) + y(i-1) - 2*y(i))/h^2

write(*,*) 'no of grid points', nop
! INITIALIZATION    
     limit = 0.0001 ! Convergence condition
     ll = 0 ! Specifies how many iterations you need for convergence 
     cond = 0! Variable whose value decides exit from the do loop

    open(080, file='b_value_dx_01_limit_0001.dat', status='unknown')
    x(1)=f_x; x(nop)=end_x; y(1)=0.0d0; y(nop)=40.0d0 ! SPECIFYING BOUNDARY CONDITIONS

    ! SETTING UP THE GRID
    do i=2,nop-1
        x(i) = x(i-1) +dx 
        ! Setting y-values at each point
        ! Use straight line with slope (end_y-f_y)/(end_x-f_x)
        y(i) = (end_y - f_y)*x(i)/(end_x - f_x) ! Initial condition 
    end do

    !write(*,*) (x(j),y(j), j=1,nop)

    do 
        ll=ll+1
        if(cond==1) exit

        ! write(*,*) 'll',ll

        y_old=y
        do i=2,nop-1 ! NOTE: Boundary points are not updated
        ! JACOBI METHOD
        ! y(i) = dum1*(dum2*y_old(i+1) + dum3*y_old(i-1) + dum4*x(i)) ! UPDATE STEP

        ! GAUSS SEIDEL METHOD
        y(i) = dum1*(dum2*y_old(i+1) + dum3*y(i-1) + dum4*x(i)) ! UPDATE STEP
        end do

        cond=1
        do i=2,nop-1
            if(abs(y_old(i) - y(i)).ge.limit) cond=0 ! CONDITION TO CHECK CONV
        end do

    end do

    write(*,*) 'No of iterations required to achieve convergence', ll

    do i=1,nop
    write(80,*) x(i), y(i)
    end do
    close(80)

end program diff_eq
