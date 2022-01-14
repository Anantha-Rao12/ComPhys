program oned_mc_impsampling

      !Author : Anantha Rao
      ! Reg no : 20181044
      integer:: n,i
      real*8:: p, crude_mc, x, func, sigma, var, length, actual

      length=3.0d0
      actual=exp(3.0d0) -1

      n = 1
      open(unit=1, file="20181044_p2_MC_sampling.dat")

      7 crude_mc = 0.0d0
      sigma=0.0d0

      do i=1,n
        call random_number(p)
        x = 3.0d0*p**(1.0/3.0)
        crude_mc = crude_mc+func(x)
        sigma=sigma+func(x)*func(x)
      end do

      crude_mc=crude_mc/real(n)
      sigma=sigma/real(n)
      var = sigma - crude_mc*crude_mc

      crude_mc = length*crude_mc
      var=length*sqrt(var/real(n))
      write(*,*) n,"",crude_mc," ", var
      write(1,*) n,"",crude_mc,"",abs(crude_mc-actual)

      !automation
      n=n*10
      if (n .le. 100000000) goto 7
      print*, "Output written to 20181044_p2_MC_impsampling.dat"
      ! end automation

      end program oned_mc_impsampling

real*8 function func(x)
    implicit none
    real*8 ::x
    func=3.0d0*x**(-2.0)*exp(x)
end function



