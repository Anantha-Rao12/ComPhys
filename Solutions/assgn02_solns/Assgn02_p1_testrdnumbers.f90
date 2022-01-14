program test_randomnumbers
      implicit none
      integer:: nbins, k, i, j
      real*8:: p, xavg, xsqr_avg,sigma, x_ik
      real*8, dimension(:), allocatable::index_array,rdn_array,ck_array,ck

      nbins = 10000
      allocate(rdn_array(nbins))
      allocate(index_array(nbins))
      allocate(ck_array(nbins-1))
      allocate(ck(nbins-1))

      ! index_array stores the indices 1:n
      do i=1,nbins-1
        index_array(i)=i
      end do

      ! rdn_array stores random numbers from the distribution
      ! 2exp(-2x)dx. We use -0.5log(x) for the same
      do i=1,nbins
        call random_number(p)
        rdn_array(i) = -0.50d0*log(p)
      end do
      
      ! We write the array of random numbers in a file
      open(unit=5, file="20181044_p1_randomnumbers.dat")
      do i=1,nbins
        write(5,*) rdn_array(i)
      end do
      print*,"Random numbers written to 20181044_p1_randomnumbers.dat "

      ! xavg stores the averages value of the random variable
      ! xsqr_avg stores the average value of the squares
      xavg=0
      xsqr_avg=0
      do i=1,nbins
        xavg = xavg+rdn_array(i)
        xsqr_avg = xsqr_avg+(rdn_array(i)*rdn_array(i))
      end do
      xavg = xavg/real(nbins)
      xsqr_avg = xsqr_avg/real(nbins)

      sigma=sqrt(xsqr_avg - (xavg*xavg))

      write(*,*) "The mean is =", xavg
      write(*,*) "The std dev is =", sigma

      ! Now we compute the Correlation values
      ! ck_arry stores the correlation values
      
      print*, "Computing the autocorrelation of the random numbers ..."
      x_ik=0
      do i=1, nbins-1
        k = i
        do j=1,nbins-k
            x_ik = x_ik + rdn_array(j)*rdn_array(j+k)
        end do
        x_ik = x_ik/real(nbins-k)
        ck_array(i) = x_ik        
      end do

      ! Now we compute the Correlation values and store that in ck
      do i=1, nbins-1
        ck(i) = (ck_array(i) - (xavg*xavg))/(sigma*sigma)
      end do

      open(unit=7, file="20181044_p1_ckvalues.dat")
      do i=1,nbins-1
        write(7,*) index_array(i), ck(i)
      end do
      print*,"Done! Autocorrelation data stored in 20181044_p1_ckvalues.dat"

       ! Print scatter plot
       print*, "Computing scatterplot"
       open(unit=9, file="20181044_p1_scattervalues.dat")
       do i=1, nbins-1
           write(9,*) rdn_array(i+1), rdn_array(i)
       end do
       print*,"Scatter plot data avaiable at 20181044_p1_scattervalues.dat"

end program test_randomnumbers
