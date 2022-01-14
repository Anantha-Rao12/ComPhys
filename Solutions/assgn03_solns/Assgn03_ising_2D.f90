! Program to generate a lattice of spin-1/2 particles and study its
! Ferromagnetic properties 

program ising_2d
    implicit none

    integer:: i,j,L,p,a,b,c,d,niter,time,mm,nn,N
    real:: r,q,E,M,mag,Ei,Ef,dE,u,h

    real:: T, J_ising=1.0 !assigning value to relevant parameters: k_B, J_ising

    integer, dimension(:,:), allocatable :: spin
    integer:: seed
    character(len=30):: charac_a, charac_b ! charac_b stores the name dump_pos

    seed=44859
    charac_b = "store_config"

    print*, "Enter the value of KbT is units of J_ising"
    read*, T
    print*, "Enter the number of lattice points in one dimension"
    read*,L
    print*, "Enter the number of iterations"
    read*, niter

    allocate(spin(L,L))
    E=0.0d0 ! Instantaneous Energy of the lattice
    M=0.0d0 ! Instantaneous magnetization of the lattice 
    N = L*L ! Total numer of spins in lattice

    ! Initialize your lattice
    open(unit=71, file="initial_ising.dat")
    p=0
    do i=1,L
        do j=1,L
            ! call RANDOM_NUMBER(r)
                spin(j,i) = 1
                ! if (r < 0.5) then 
                !     spin(j,i) = -1
                ! else
                !     spin(j,i)=1
                ! End if
            ! Write down this configuration 
            !write(71, fmt="(4g10.8)") float(i), float(j), float(p), float(spin(j,i))
            write(71,*) float(i), float(j), float(p), float(spin(j,i))
        end do
    end do
    close(71)

    ! Calculate Initial magnetization and energy 
    do i=1,L
        do j=1,L  
            a=i+1; b=i-1; c=j+1; d=j-1 !identifying the 4 neighbors of spin(j,i) 

            ! Setting PBCs (Periodic Boundary conditions) 
            if (i==L)a=1
            if (i==1)b=L
            if (j==1)d=L
            if (j==L)c=1

            E=E-J_ising*float((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))
            M=M+spin(i,j)
        end do
    end do

    ! Calculate instantaneous magnetization per spin 
    mag = M/float(N)
    E=E*0.5d0
    print*, "Initial energy E, E per spin is: ",E, E/float(N)
    print*, "Initial magentization M, M per spin is: ", M,mag

    ! INITIALIZATION COMPLETE
    ! -----------------------------
    ! EVOLVE TO REACH EQUILIBRIUM   

    open(unit=10, file="ising_T_Linit_random.dat")
    do time=i,niter ! loop over number of MCS
        do mm=1,L
            do nn=1,L
                call random_number(r); i=int(r*float(L))+1 !choosing a lattice site
                call random_number(r); j=int(r*float(L))+1

                a=i+1;b=i-1;c=j+1;d=j-1 !identify the neighbors of spin(j,i)   

                if(i==L)a=1; if(i==1)b=L; if(j==1)d=L; if(j==L)c=1 !PBC
                ! BEFORE TRAIL FLIP
                Ei = -J_ising*float((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))

                spin(i,j) =-spin(i,j) !TRIAL FLIP

                ! AFTER TRAIL FLIP
                Ef = -J_ising*float((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))
                dE = Ef-Ei

                if(dE <=0.0) then
                    E=E+dE
                    M=M+(2.0*float(spin(i,j)))
                else
                    u=exp(-dE/(T))
                    call random_number(h)
                    if (h<u) then
                        E=E+dE
                        M=M+(2.0*float(spin(i,j))) ! Instantaneous mag. Of entire lattice
                    else
                        spin(i,j) = -spin(i,j) ! Trial flip not accepted; E & M not updated
                    end if
                end if
            end do
        end do
    write(10,*) time, M/float(N), E/float(N) ! Writing down E and M with no of iterations
    end do
    close(10)


end program ising_2d

