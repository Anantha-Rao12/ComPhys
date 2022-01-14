program cubic_lattice
    implicit none 
    integer :: i, j, k, nx=20, ny=10, nz=10
    real*8 :: x, y, z
    open(23, file='cubic.xyz', action='write')
    do i=1,nx
        do j=1,ny
            do k=1,10
                x=dfloat(i)-1.0d0
                y=2.0d0*dfloat(j)-1.0d0
                z=2.0d0*dfloat(k)-1.0d0
                write(23,*) x, y, z
            end do
        end do
    end do
    close(23)
end program cubic_lattice

