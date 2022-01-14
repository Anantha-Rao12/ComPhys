!---------------------------------------------------------------
program potentialwell
  !---------------------------------------------------------------
  !
  !     Lowest energy levels of the finite (symmetric) potential well
  !     via expansion on a plane-wave basis set and diagonalization
  !     Units: hbar^2/2m = 1
  !     Requires lapack dsyev

  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi=3.14159265358979_dp
  integer :: n, npw
  real(dp) :: v0, a, b
  real(dp), allocatable :: kn(:), e(:), h(:,:), work (:)
  real(dp) :: x, dx, norm, prob
  complex(dp) :: f
  integer :: i, j, nr, lwork, info
  !
  !       Input data
  !
  !     Potential well: V(x)=-V_0 for |x|<b/2, V(x)=0 for |x|>b/2
  !
  write (*,"('Parameters for potential well: V_0, b > ')", advance="no")
  read (*,*) v0, b
  if ( v0 <= 0.0_dp .or. b <= 0.0_dp) stop ' wrong input parameters, stopping '
  write (*,"('   V_0, b =',2f10.4)") v0, b
  !
  !     Plane waves between -a/2 < x < a/2, k_i=+-2*pi*i/a, i=0,1,...,n
  !
  do
     write (*,"('Parameters for plane waves: a, n > ')", advance="no")
     ! To exit the loop, type 0 0 or something similar
     read (*,*) a, n
     if ( a <= b .or. n <= 0) stop ' wrong input parameters '
     write (*,"('a, n=',f8.4,i6)") a, n
     npw = 2*n+1
     allocate (kn(npw), e(npw), work(3*npw), h(npw,npw) )
     !
     !       Assign values of k_n: n=0,+1,-1,+2,-2, etc
     !
     kn(1) = 0.0_dp
     do i=2,npw-1,2
        kn(i  ) = (i/2)*2.0_dp*pi/a
        kn(i+1) =-(i/2)*2.0_dp*pi/a
     end do
     !       cleanup
     h(:,:) = 0.0_dp
     !
     !       Assign values of the matrix elements of the hamiltonian 
     !       on the plane wave basis
     !
     do i=1,npw
        do j=1,npw
           if ( i ==j ) then
              h(i,j) = kn(i)**2 - v0/a*b
           else
              h(i,j) = -v0/a * sin( (kn(j)-kn(i))*b/2.0_dp ) / (kn(j)-kn(i))*2.0_dp
           end if
           !print  '(2i4,f12.6)', i,j, h(i,j)
        end do
     end do
     !
     !       Solution [expansion coefficients are stored into h(j,i)
     !                 j=basis function index, i= eigenvalue index]
     !
     lwork = 3*npw
     !       The "leading dimension of array" h is its first dimension
     !       For dynamically allocated matrices, see the "allocate" command
     call dsyev ( 'V', 'U', npw, h, npw, e, work, lwork, info )
     if (info /= 0) stop 'H-matrix diagonalization failed '
     !
     write (*,"('   Lowest eigenvalues:',3f12.6)") (e(i),i=1,3)
     !
     !       Write to output file the lowest-energy state:
     !
     open (7,file='gs-wfc.out',status='unknown',form='formatted')
     write(7,'("        x      |psi(x)|^2  Re[psi(x)]   Im[psi(x)]")')
     dx = 0.01_dp
     nr = nint(a/2.0_dp/dx)
     norm = 0.d0
     do i=-nr, nr
        x = dx*i
        f = 0.d0
        do j=1,npw
           f = f + h(j,1)*exp((0.0,1.0)*kn(j)*x)/sqrt(a)
        end do
        prob = f*conjg(f)
        norm = norm + prob*dx
        write(7,'(f12.6,3f12.6)') x, prob, f
     end do
     !       verify normalization (if desired):
     write (*,"('   norm: ',f12.6)") norm
     close(7)
     deallocate ( h, work, e, kn)
  end do
  !
end program potentialwell


