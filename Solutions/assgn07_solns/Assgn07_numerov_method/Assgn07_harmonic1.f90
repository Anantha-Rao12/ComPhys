!---------------------------------------------------------------
program harmonic1
!---------------------------------------------------------------
  !
  !       Solution of the quantum harmonic oscillator.
  !       Forward and backward integration with Numerov method.
  !       Solution matching at a classical turning point.
  !       Eigenvalue search using the shooting method.
  !       Adimensional units: x = (mK/hbar^2)^(1/4) X
  !                           e = E/(hbar omega)
  !
  implicit none
  ! In next line, "dp is the kind having at least 14 significan digits,
  ! range between 10^-200 and 10^200: in practice, double precision
  integer, parameter :: dp = selected_real_kind(14,200)
  !
  integer :: mesh, i, icl
  integer :: nodes, hnodes, ncross, kkk, n_iter
  real(dp) :: xmax, dx, ddx12, norm, arg, djump, fac
  real(dp) :: eup, elw, e
  real(dp), allocatable :: x(:), y(:), p(:), vpot(:), f(:)
  character (len=80) :: fileout
  !
  ! read input data
  !
  write(*,"('Max value for x (typical value: 10) > ')", advance='no') 
  read (*,*) xmax
  write(*,"('Number of grid points (typically a few hundreds) > ')", advance='no')
  read (*,*) mesh
  !
  ! allocate arrays, initialize grid
  !
  allocate ( x(0:mesh), y(0:mesh), p(0:mesh), vpot(0:mesh), f(0:mesh) )
  dx =  xmax/mesh 
  ddx12=dx*dx/12.0_dp
  !
  ! set up the potential (must be even w.r.t. x=0)
  !
  do i = 0, mesh
     x(i) = float(i) * dx
     vpot(i) = 0.5_dp * x(i)*x(i)
  end do
  !
  write(*,"('output file name > ')", advance='no') 
  read (*,'(a)') fileout
  ! output to fort.7 file by default if fileout not provided
  if ( fileout /= ' ' ) &
      open (7,file=fileout,status='unknown',form='formatted')

  ! this is the entry point for a new eigenvalue search
  search_loop: do 
     !
     ! read number of nodes (stop if < 0)
     !         
     write(*,"('nodes (type -1 to stop) > ')", advance='no') 
     read (*,*) nodes
     if (nodes < 0) then
        close(7)
        deallocate ( f, vpot, p, y, x )
        stop 
     end if
     !
     ! set initial lower and upper bounds to the eigenvalue
     !
     eup=maxval (vpot(:))
     elw=minval (vpot(:))
     !
     ! Set trial energy
     !
     write(*,"('Trial energy (0=search with bisection) > ')", advance='no') 
     read (*,*) e
     if ( e == 0.0_dp ) then
        ! search eigenvalues with bisection (max 1000 iterations)
        e = 0.5_dp * (elw + eup)
        n_iter = 1000
     else
        ! test a single energy value (no bisection)
        n_iter = 1
     endif
     
     iterate: do kkk = 1, n_iter
        !
        ! set up the f-function used by the Numerov algorithm
        ! and determine the position of its last crossing, i.e. change of sign
        ! f < 0 means classically allowed   region
        ! f > 0 means classically forbidden region
        !
        f(0)=ddx12*(2.0_dp*(vpot(0)-e))
        icl=-1
        do i=1,mesh
           f(i)=ddx12*2.0_dp*(vpot(i)-e)
           ! beware: if f(i) is exactly zero the change of sign is not observed
           ! the following line is a trick to prevent missing a change of sign 
           ! in this unlikely but not impossible case:
           if ( f(i) == 0.0_dp) f(i)=1.d-20
           ! store the index 'icl' where the last change of sign has been found
           if ( f(i) /= sign(f(i),f(i-1)) ) icl=i
        end do
        
        if (icl >= mesh-2) then
           deallocate ( f, vpot, p, y, x )
           print *, 'Error: last change of sign too far'
           stop 1
        else if (icl < 1) then
           deallocate ( f, vpot, p, y, x )
           print *, 'Error: no classical turning point'
           stop 1
        end if
        !
        ! f(x) as required by the Numerov algorithm (note f90 array syntax)
        !
        f = 1.0_dp - f
        y = 0.0_dp
        !
        ! determination of the wave-function in the first two points 
        !
        hnodes = nodes/2
        !
        ! beware the integer division: 1/2 = 0 !
        ! if nodes is even, there are 2*hnodes nodes
        ! if nodes is odd,  there are 2*hnodes+1 nodes (one is in x=0)
        ! hnodes is thus the number of nodes in the x>0 semi-axis (x=0 excepted)
        !
        if (2*hnodes == nodes) then
           ! even number of nodes: wavefunction is even
           y(0) = 1.0_dp
           ! assume f(-1) = f(1)
           y(1) = 0.5_dp*(12.0_dp-10.0_dp*f(0))*y(0)/f(1)
        else
           ! odd number of nodes: wavefunction is odd
           y(0) = 0.0_dp
           y(1) = dx
        end if
        !
        ! outward integration and count number of crossings
        !
        ncross=0
        do i =1,icl-1
           y(i+1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
           if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
        end do
        fac = y(icl)
        !
        if (2*hnodes == nodes) then
           ! even number of nodes: no node in x=0
           ncross = 2*ncross
        else
           ! odd number of nodes: node in x=0
           ncross = 2*ncross+1
        end if
        !
        ! check number of crossings
        !
        if ( n_iter > 1 ) then
           if (ncross /= nodes) then
              ! Incorrect number of crossings: adjust energy
              if ( kkk == 1) &
                   print '("Bisection         Energy       Nodes  Discontinuity")'
              print '(i5,f25.15,i5)', kkk, e, ncross
              if (ncross > nodes) then
                 ! Too many crossings: current energy is too high,
                 ! lower the upper bound
                 eup = e
              else 
                 ! Too few crossings: current energy is too low,
                 ! raise the lower bound
                 elw = e
              end if
              ! New trial value:
              e = 0.5_dp * (eup+elw)
              ! go to beginning of do loop, don't perform inward integration
              cycle
           end if
        else
           print *, e, ncross, nodes
        end if
        !
        ! Correct number of crossings: proceed to inward integration
        !
        ! Determination of the wave-function in the last two points
        ! assuming y(mesh+1) = 0
        !
        y(mesh) = dx
        y(mesh-1) = (12.0_dp-10.0_dp*f(mesh))*y(mesh)/f(mesh-1)
        norm = 1.0d100
        do i = mesh-1,icl+1,-1
           y(i-1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i+1)*y(i+1))/f(i-1)
           ! the following lines prevent overflows if starting from too far
           if ( abs(y(i-1)) > norm ) then
              y(i-1:mesh) = y(i-1:mesh) / norm
           endif
        end do
        !
        ! rescale function to match at the classical turning point (icl)
        !
        fac = fac/y(icl)
        y(icl:) = y(icl:)*fac
        !
        ! normalize on the [-xmax,xmax] segment
        ! the x=0 point must be counted once
        !
        norm = (2.0_dp*dot_product (y, y) - y(0)*y(0)) * dx 
        y = y / sqrt(norm)
        !
        if ( n_iter > 1 ) then
           !  calculate the discontinuity in the first derivative
           !  y'(i;RIGHT) - y'(i;LEFT)
           djump = ( y(icl+1) + y(icl-1) - (14.0_dp-12.0_dp*f(icl))*y(icl) ) / dx
           print '(i5,f25.15,i5,f14.8)', kkk, e, nodes, djump
           if (djump*y(icl) > 0.0_dp) then
              ! Energy is too high --> choose lower energy range
              eup = e
           else
              ! Energy is too low  --> choose upper energy range
              elw = e
           endif
           e = 0.5_dp * (eup+elw)
           ! ---- convergence test
           if ( eup-elw < 1.d-10) exit iterate
        endif
        !
     end do iterate
     !
     ! ---- convergence has been achieved (or it wasn't required) -----
     !
     ! Calculation of the classical probability density for energy e:
     !
     norm = 0.0_dp
     p(icl:) = 0.0_dp
     do i=0,icl
        arg = (e - x(i)**2/2.0_dp)
        if ( arg > 0.0_dp) then
           p(i) = 1.0_dp/sqrt(arg)
        else
           p(i) = 0.0_dp
        end if
        norm = norm + 2.0_dp*dx*p(i)
     enddo
     ! The point at x=0 must be counted once:
     norm = norm - dx*p(0)
     ! Normalize p(x) so that  Int p(x)dx = 1
     p(:icl-1) = p(:icl-1)/norm
     ! lines starting with # ignored by gnuplot
     write (7,'("#   x         y(x)            y(x)^2      classical p(x)    V")')
     ! x<0 region:
     do i=mesh,1,-1
        ! if exponent is > 99, the format X.Y-100 is misinterpreted by gnuplot
        if ( abs(y(i)) < 1.0D-50 ) y(i) = 0.0_dp
        write (7,'(f8.3,3e16.8,f12.6)') & 
             -x(i), (-1)**nodes*y(i), y(i)*y(i), p(i), vpot(i)
     enddo
     ! x>0 region:
     do i=0,mesh
        write (7,'(f8.3,3e16.8,f12.6)') &
             x(i), y(i), y(i)*y(i), p(i), vpot(i)
     enddo
     ! two blank lines separating blocks of data, useful for gnuplot plotting
     write (7,'(/)')
     
  end do search_loop
  
end program harmonic1


