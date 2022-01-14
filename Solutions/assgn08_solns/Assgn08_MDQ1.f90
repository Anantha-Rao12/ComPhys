! Fortran code to perform moelcular dynamics simulations
! Author: Anantha Rao
! Reg no:20181044

module param
    implicit none 
    integer, parameter :: n_part=1000
    real*8, parameter:: T = 1.0d0, m=1.0d0, kbT=1.0d0
    real*8, parameter:: sigma = 1.0d0, eps = 4.0d0, rc = 2.5d0
    real*8, parameter:: tf=16.0d0, dt=0.005d0
    real*8, parameter:: lx=20.0d0, ly=20.0d0, lz=20.0d0
    real*8, parameter:: sigma6=sigma**6, sigma12=sigma**12, sigma2=sigma*sigma
    real*8, parameter:: fc=eps*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
    real*8, parameter:: ufc = fc*rc+eps*(((sigma/rc)**12)-(sigma/rc)**6)
    real*8, parameter:: lxby2=lx/2.0d0, lyby2=ly/2.0d0, lzby2=lz/2.0d0
    real*8 :: vx(n_part), vy(n_part), vz(n_part)
    real*8 :: x(n_part), y(n_part), z(n_part)
    real*8 :: xp(n_part), yp(n_part), zp(n_part)
    real*8 :: fx(n_part), fy(n_part), fz(n_part)
    real*8 :: xnew(n_part), ynew(n_part), znew(n_part)
    real*8 :: xx(n_part), yy(n_part), zz(n_part)
    real*8 :: dx, dy, dz, dr
    real*8 :: avg_vx,avg_vy, avg_vz, sumv2, vxbar, vybar, vzbar
    real*8 :: randm, fs, force, ke, pe, te, lj, ljforces, theoryke, scalef
    real*8 :: x1, x2, y1, y2, z1, z2, tt
    integer :: i,j,k 
end module param

program md_verlet
    use param 
    implicit none
    open (12, file='cubic.txt', action='read')
    open(13, file='data.dat', action='write')
    open(14, file='momt.dat', action='write')
    do i = 1, n_part
        read(12,*) x(i), y(i), z(i)
    end do
    write(*,*) n_part, "Particle Positions Initialized"
    call vel_init
    call calcforces
    write(*,*) n_part, "Particle Velocities Initialized"
    write(*,*) "Potential Energy", pe/dfloat(n_part)
    write(*,*) "Kinetic Energy", ke/dfloat(n_part)
    
    do i = 1, n_part
        xp(i) = x(i)-vx(i)*dt
        yp(i) = y(i)-vy(i)*dt
        zp(i) = z(i)-vz(i)*dt
    end do
    
    tt=0.0d0
    do while (tt<= tf)
        call calcforces 
        call integrate
        write(13,*)tt, ke/dfloat(n_part), pe/dfloat(n_part), (ke+pe)/dfloat(n_part)
        write(14,*) tt, avg_vx, avg_vy, avg_vz
        write(*,*) tt, avg_vx, avg_vy, avg_vz
        tt=tt+dt
    end do 
    close(12)
    close(13)
    close(14)
end program md_verlet
subroutine vel_init
    use param 
    implicit none
    avg_vx = 0.0d0; avg_vy=0.0d0; avg_vz=0.0d0
    sumv2=0.0d0
    ke=0.0d0
    do i = 1, n_part
        call random_number(randm); vx(i) = (randm-0.5d0)
        call random_number(randm); vy(i) = (randm-0.5d0)
        call random_number(randm); vz(i) = (randm-0.5d0)
    end do
    
    do i = 1, n_part
        avg_vx=avg_vx+vx(i)
        avg_vy=avg_vy+vy(i)
        avg_vz=avg_vz+vz(i)
    end do
    
    avg_vx=avg_vx/dfloat(n_part)
    avg_vy=avg_vy/dfloat(n_part)
    avg_vz=avg_vz/dfloat(n_part)
    
    vxbar = 0.0d0
    vybar = 0.0d0
    vzbar = 0.0d0
    
    do i = 1, n_part
        vx(i) = (avg_vx-vx(i))
        vy(i) = (avg_vy-vy(i))
        vz(i) = (avg_vz-vz(i))
    end do 
    
    do i = 1, n_part
        vxbar=vxbar+vx(i)
        vybar=vybar+vy(i)
        vzbar=vzbar+vz(i)
    end do
    
    do i = 1, n_part
        sumv2=sumv2+(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
    end do
    
    vxbar=vxbar/dfloat(n_part)
    vybar=vybar/dfloat(n_part)
    vzbar=vzbar/dfloat(n_part)
    fs = dsqrt((3.0d0*T*(n_part-1.0d0))/sumv2)
    do i = 1, n_part
        vx(i) = vx(i)*fs
        vy(i) = vy(i)*fs
        vz(i) = vz(i)*fs
        ke=ke+(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
    end do
    ke = ke*0.5d0

end subroutine vel_init

subroutine calcforces 
    use param 
    implicit none
    pe = 0.0d0
    do i = 1, n_part
        fx(i) = 0.0d0; fy(i) = 0.0d0; fz(i) = 0.0d0
    end do
    
    do i = 1, n_part-1
        x1=x(i); y1=y(i); z1=z(i)
        do j = i+1, n_part
            x2=x(j); y2=y(j); z2=z(j)
            dx = x1-x2; dy = y1-y2; dz=z1-z2
            if(abs(dx).ge.lxby2) dx=(lx-abs(dx))*((-1.0d0*dx/abs(dx)))
            if(abs(dy).ge.lyby2) dy=(ly-abs(dy))*((-1.0d0*dy/abs(dy)))
            if(abs(dz).ge.lzby2) dz=(lz-abs(dz))*((-1.0d0*dz/abs(dz)))
            
            dr = dsqrt(dx*dx+dy*dy+dz*dz)
            if (dr <= rc) then 
                lj = eps*((sigma/dr)**12-(sigma/dr)**6)-ufc+fc*dr
                pe = pe + lj
                ljforces = eps*((12.0d0*((sigma12)/(dr)**13))-(6.0d0*((sigma6/(dr)**7))))-fc
                fx(i) = fx(i) + ljforces*(dx/dr)
                fy(i) = fy(i) + ljforces*(dy/dr)
                fz(i) = fz(i) + ljforces*(dz/dr)
                fx(j) = fx(j) - ljforces*(dx/dr)
                fy(j) = fy(j) - ljforces*(dy/dr)
                fz(j) = fz(j) - ljforces*(dz/dr)
                
            end if
        end do 
    end do
end subroutine calcforces

subroutine integrate
    use param 
    implicit none
    real*8::dt2
    dt2=dt**2
    avg_vx=0.0d0; avg_vy=0.0d0; avg_vz=0.0d0; ke=0.0d0
    do i = 1, n_part
        xnew(i) = 2.0d0*x(i)-xp(i)+fx(i)*dt2
        ynew(i) = 2.0d0*y(i)-yp(i)+fy(i)*dt2
        znew(i) = 2.0d0*z(i)-zp(i)+fz(i)*dt2
        vx(i) = (xnew(i)-xp(i))/(2.0d0*dt)
        vy(i) = (ynew(i)-yp(i))/(2.0d0*dt)
        vz(i) = (znew(i)-zp(i))/(2.0d0*dt)
        ke = ke + 0.5d0*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
        avg_vx=avg_vx+(vx(i)/dfloat(n_part))
        avg_vy=avg_vy+(vy(i)/dfloat(n_part))
        avg_vz=avg_vz+(vz(i)/dfloat(n_part))
        
    end do
    
    !if (mod(int(tt/0.005d0),100)==0) then
    !    theoryke=1.50d0*dfloat(n_part)*kbT
    !    scalef=dsqrt(theoryke/ke)
    !    vx=vx*scalef
    !    vy=vy*scalef
    !    vz=vz*scalef
    !    ke=0.0d0
    !    do i=1,n_part
    !        ke=ke+0.50d0*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
    !    end do 
    !endif
    
    do i = 1, n_part
        xnew(i) = xp(i) + 2.0d0*vx(i)*dt
        ynew(i) = yp(i) + 2.0d0*vy(i)*dt
        znew(i) = zp(i) + 2.0d0*vz(i)*dt
        xp(i)=x(i)
        yp(i)=y(i)
        zp(i)=z(i)
        
        x(i) = xnew(i)
        y(i) = ynew(i)
        z(i) = znew(i)
        
        if (x(i)>lx) then 
            x(i) = x(i)-lx
            xp(i) = xp(i)-lx
        elseif(x(i)<0.0d0) then 
            x(i) = x(i)+lx  
            xp(i) = xp(i)+lx
        end if
        if (y(i)>ly) then 
            y(i) = y(i)-ly
            yp(i) = yp(i)-ly
        elseif(z(i)<0.0d0) then 
            z(i) = z(i)+lz 
            zp(i) = zp(i)+lz
        end if
        if (z(i)>lz) then 
            z(i) = z(i)-lz
            zp(i) = zp(i)-lz
        elseif(x(i)<0.0d0) then 
            z(i) = z(i)+lz  
            zp(i) = zp(i)+lz
        end if
    end do
end subroutine integrate
