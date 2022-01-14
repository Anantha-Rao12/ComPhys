! Program to generate a lattice of spin-1/2 particles and study its
! Ferromagnetic properties 

program ising_3d
    implicit none

    integer:: i,j,k,L,p,a,b,c,d,f,g,niter,time,mm,nn,oo,N
    real*8:: r,E,M,mag,Ei,Ef,dE,u,h

    real*8:: T, J_ising=1.0 !assigning value to relevant parameters: k_B, J_ising

    integer, dimension(:,:,:), allocatable :: spin
    integer:: seed
    ! character(len=30):: charac_a, charac_b ! charac_b stores the name dump_pos

    seed=44859
    ! charac_b = "store_config"

    print*, "Welcome to this program! This program simulates a 3D ising model with all spins iniitally in the random state"
    print*, "Enter the value of KbT in units of J_ising"
    read*, T
    print*, "Enter the number of lattice points in one dimension for (LxLxL)"
    read*,L
    print*, "Enter the number of iterations"
    read*, niter

    allocate(spin(L,L,L))
    E=0.0d0 ! Instantaneous Energy of the lattice
    M=0.0d0 ! Instantaneous magnetization of the lattice 
    N = L*L*L ! Total number of spins in lattice

    ! Initialize your lattice
    open(unit=71, file="initial_ising_3D.dat")
    p=0
    do i=1,L
        do j=1,L
            do k=1,L
             call RANDOM_NUMBER(r)
                ! spin(k,j,i) = -1
                if (r < 0.5) then 
                    spin(i,j,k) = -1
                else
                    spin(i,j,k)=1
                end if
            ! Write down this configuration 
            !write(71, fmt="(5g10.8)") float(i), float(j), float(p), float(spin(j,i))
            write(71,*) float(i), float(j), float(k), float(p), float(spin(i,j,k))
            end do
        end do
    end do

    close(71)

    ! Calculate Initial magnetization and energy 
    do i=1,L
        do j=1,L  
            do k=1,L
            a=i+1; b=i-1; c=j+1; d=j-1; f=k+1; g=k-1 !identifying the 4 neighbors of spin(j,i) 

            ! Setting PBCs (Periodic Boundary conditions) 
            if (i==L)a=1
            if (i==1)b=L
            if (j==L)c=1
            if (j==1)d=L
            if (k==L)f=1
            if (k==1)g=L

            E=E-J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,f)+spin(i,j,g)))
            M=M+spin(i,j,k)
            end do
        end do
    end do

    ! Calculate instantaneous magnetization per spin 
    mag = M/float(N)
    E=E*0.5d0
    print*, "Initial energy E, E per spin is: ",E, E/float(N)
    print*, "Initial magentization M, M per spin is: ", M,M/float(N)

    ! INITIALIZATION COMPLETE
    ! -----------------------------
    ! EVOLVE TO REACH EQUILIBRIUM   

    open(unit=10, file="ising_3d_T_L20_random.dat")
    do time=i,niter ! loop over number of MCS
        do mm=1,L
            do nn=1,L
                do oo=1,L
                call random_number(r); i=int(r*float(L))+1 !choosing a lattice site
                call random_number(r); j=int(r*float(L))+1
                call random_number(r); k=int(r*float(L))+1

                a=i+1;b=i-1;c=j+1;d=j-1;f=k+1;g=k-1 !identify the neighbors of spin(j,i)   

                if(i==L)a=1; if(i==1)b=L; if(j==L)c=1; if(j==1)d=L; if(k==L)f=1; if(k==1)g=L !PBC
                ! BEFORE TRAIL FLIP
                Ei = -J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,f)+spin(i,j,g)))

                spin(i,j,k) =-spin(i,j,k) !TRIAL FLIP

                ! AFTER TRAIL FLIP
                Ef = -J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,f)+spin(i,j,g)))
                dE = Ef-Ei

                if(dE <=0.0) then
                    E=E+dE
                    M=M+(2.0*float(spin(i,j,k)))
                else
                    u=exp(-dE/(T))
                    call random_number(h)
                    if (h<u) then
                        E=E+dE
                        M=M+(2.0*float(spin(i,j,k))) ! Instantaneous mag. Of entire lattice
                    else
                        spin(i,j,k) = -spin(i,j,k) ! Trial flip not accepted; E & M not updated
                    end if
                end if
                end do
            end do
        end do
    write(10,*) time, M/float(N), E/float(N) ! Writing down mag and E/float(N) with no of iterations
    end do
    close(10)


end program ising_3d

