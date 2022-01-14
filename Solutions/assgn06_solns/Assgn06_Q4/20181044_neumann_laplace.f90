! Program that solves the Laplace equation using Neumann's Boundary conditions
! Author: Anantha Rao
! Reg no: 20181044

program laplace_neumann
implicit none

      integer, parameter :: lx=34, ly=34 ! BOUNDARIES at x,y=0 and 34
      real*8 :: old_temp(1:lx, 1:ly), temp(1:lx, 1:ly)
      integer:: i,j,ii,jj,kk
      real*8 :: dx,dy, prefactor, shift_temp
      real*8 :: A(ly), B(ly), C(lx), D(lx)
      character(len=30) :: filename
      integer :: iunit, ci, cii, test, counter

      A = -70.0d0; B=-40.0d0; C=20.0d0; D=-10.0d0
      old_temp=0.0d0 ! at ZEROTH iteration

      dx=1.0d0; dy=dx
      test=0; counter=0; prefactor=(0.5d0*dx*dx*dy*dy)/(dx*dx + dy*dy)
      
      do  !LOOP OVER ITERATIONS
        counter=counter+1
        test=0

        ! UPDATE BOUNDARIES: NEUMANN CONDITIONS: LEAVING OUT CORNERS
        do j=2, ly-1
          temp(1,j) = 0.25d0*(2.0d0*old_temp(2,j) - 2.0d0*dx*A(j) + old_temp(1,j+1) + old_temp(1,j-1))
          temp(lx,j) = 0.25d0*(2.0d0*old_temp(lx-1,j) + 2.0d0*dx*B(j) + old_temp(lx,j+1) + old_temp(lx,j-1))
        end do

        do i=2,lx-1
          temp(i,1) =0.25d0*( old_temp(i+1,1)+ old_temp(i-1,1)+ 2.0d0*old_temp(i,2)- 2.0d0*dx*C(i))
          temp(i,ly) = 0.25d0*( old_temp(i+1,ly)+ old_temp(i-1,ly) + 2.0d0*old_temp(i,ly-1) + 2.0d0*dx*D(i))
        end do

        ! UPDATE THE VALUES AT 4 CORNERS
        temp(1,1)   = 0.5d0*(old_temp(1,2)     - dx*C(1)  + old_temp(2,1) - dx*A(1))
        temp(1,ly)  = 0.5d0*(old_temp(1,ly-1)  + dx*D(1)  + old_temp(2,ly)- dx*A(ly))
        temp(lx,1)  = 0.5d0*(old_temp(lx-1,1)  + dx*B(1)  +old_temp(lx,2)-dx*C(lx))
        temp(lx,ly) = 0.5d0*(old_temp(lx-1,ly) + dx*B(ly) +old_temp(lx,ly-1) + dx*D(lx))

        write(*,*) 'before update step', counter
        ! UPDATE STEP: THE INTERIOR OF THE LATTICE
        do jj=2, ly-1
        do ii=2, lx-1
            temp(ii,jj) = 0.25*( old_temp(ii-1,jj) + old_temp(ii+1,jj) + old_temp(ii,jj-1) + old_temp(ii,jj+1))
        end do
        end do

        write(*,*) 'before convergence check'
        ! CHECK FOR CONVERGENCE AT EACH LATTICE SITE
        do jj=1,ly
        do ii=1,lx
          if ((abs(temp(ii,jj) - old_temp(ii,jj))).gt.0.00001d0) test=1
        end do
        end do

        if(test.eq.0) exit ! EXIT CONDITION
        shift_temp =2000 - temp(1,1) ; temp=temp+shift_temp
        old_temp = temp  !AFTER TEST CONDITION

      end do
      write(*,*) 'counter',counter


      ! WRITE THE CONVERGED RESULT: SOLUTION TO THE EQN WITH BC
      iunit=71; ci=10000
      ! write(filename,'("initialize_",i0,".dat")') ci

      open(iunit, file='data.dat', status='new')
      do ii=1,lx
        do jj=1,ly
            write(iunit,*) ii,jj,temp(ii,jj)
        end do
      end do
      close(iunit)

!     open(newunit=iunit, file='data.dat')
!      do ii=1,lx
!        do jj=1,ly
!            write(iunit,*) ii,jj,temp(ii,jj)
!        end do
!      end do
!      close(iunit)


end program laplace_neumann
