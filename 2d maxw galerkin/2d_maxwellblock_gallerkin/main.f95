program maxblock_gallerkin

    !18/03/2016
    use constants
    use funcs
    use binom

    implicit none

    real(8) :: intime
    real(8) :: timeinit=0.d0
    real(8), allocatable, dimension(:) :: time
    !real(8), dimension(9) :: yinit
    integer :: numstep
    integer :: indexkeep

    real(8), allocatable, dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), allocatable, dimension(:) :: intensity_x2, intensity_y2, intensity

    real(8), allocatable, dimension (:,:) :: loss, pump, Atest
    integer :: rho_dim=200
    integer :: phi_dim=30
    real(8), allocatable, dimension(:) :: rho, phi
    real(8) :: rn,phin
    real(8) :: rho_0=0.1
    integer :: numkeep

    integer :: p,m,i
    real(8) :: ln, Alag
    integer :: factorial,factorial_ros, choose
    real(8) :: binomial

    external factorial_ros, choose
!    real(8) :: ln

    integer :: z,c,j
    real(8) :: tempvar
    integer, parameter :: debug=1

    !************************************************

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    !1-Program Scheme. Generate the tables needed for the integrations.
    !2-From the tables: the absissas for rho(from the generalized laguerre rule), and the absissas for phi(from the chebychev rule), and the weights.
    !3-Generate a table for the nonlinear term
    !4-Proyect my functions with the laguerre-chevy interpolation.
    !5-evolve all the coefficients,
    !6-for each time, keep the coeffs in a table. And keep the funcition for a fixed phi=0.(how good is this??)


    !1) -To do: Create the table, calling the programs. And find a way to read the values.
        !Problems: a) The binomial coefficient doens-t work well if there is a factorial > 16
                  !b)
    !2) -To do: Make two array with the absissas and Rearange them.
    !3) -To do: Use the cuadrature to generate the nonlinear table.
    !4) -To do: Use the quadratures, and make a table for the coefficients(n,l,i).
    !5) -To do: Loop runge kutta for the coeffs.
    !6) -To do: keep everything nice and clean.

    allocate(rho(rho_dim))
    call linspace(rho,0d0,20d0,rho_dim,1)
    allocate(phi(phi_dim))
    call linspace(phi,0d0,2*pi,phi_dim,1)

    allocate(loss(phi_dim,rho_dim))!phi en las filas, rho en las colmnas
    do c=1,rho_dim
        do j=1,phi_dim
             loss(j,c)=5d0+4d0*tanh(5d0*(rho(c)-rho_0))
        end do
    end do


    allocate(pump(phi_dim,rho_dim))!phi en las filas, rho en las colmnas
    do c=1,rho_dim
        do j=1,phi_dim
             pump(j,c)=5d0-tanh(5d0*(rho(c)-rho_0))/2d0
        end do
    end do

!    call gen_laguerre_test() !generalized lagurre rule(quadrature)
!    call cheby1()           !order n chebychev rule

    print*,shape(loss)
    print*, 'factorial', factorial(19)
    print*, 'factorial_ros', factorial_ros(19)
    print*, 'binomial', binomial(26,8)
    print*, 'choose', choose(26,8)
    print*, 'ln', ln(1,1,1d0)
    print*, 'Alag', Alag(1,1,1,0d0,1d0)
  !  print*, 'binomial from stack', combo(8,15)
    !test to plot some A
    allocate(Atest(phi_dim,rho_dim))
    do c=1,rho_dim
        do j=1, phi_dim
            Atest(j,c)=Alag(0,0,0,rho(c),phi(j))
        end do
    end do

    print*, shape(Atest)
    open(1,file='atest.in', form='unformatted')
        write(1) Atest
    close(1)


end program

integer function factorial(num) !breaks after 16
implicit none
    integer :: c,j
    integer, intent(in) :: num
    c = 1
    DO j = 2, num
        c = c*j
    END DO
    factorial=c
endfunction

double precision function binomial(n,k)
implicit none
    integer, intent(in) :: n,k
    integer :: factorial
    if (n.ge.k) then
        binomial=dble(factorial(n)/(factorial(n-k)*factorial(k)))
    else
    end if
end function


 function factorial_ros (n) result (res)

    implicit none
    integer, intent (in) :: n
    integer :: res
    integer :: i

    res = product ((/(i, i = 1, n)/))

  end function factorial_ros

  function choose (n, k) result (res)
 !rosetta code
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    integer :: res
    integer :: factorial_ros
    external factorial_ros
    res = factorial_ros (n) / (factorial_ros (k) * factorial_ros (n - k))

  end function choose


  !end function
real*8 function ln(p,m,x)
! The recursion relation is
!(n + 1) Ln+l(x) = (2n + 1 - x) Ln(x) - nLn(x)   (3.9.7)
!where L0(x) = 1
!L1(x) = 1 - x

    real(8) :: lntemp
    real(8),intent(in) :: x
    integer,intent(in) :: m,p
    integer :: z
    real(8) :: ln0,ln1,ln2

    if (p.eq.0) then
        ln0=1
        ln=ln0
    elseif (p.eq.1) then
        ln1=1+m+x
        ln=ln1
    elseif (p.gt.1) then
        ln0=1
        ln1=1+m+x
        do z=2,p
            ln2=((2*(z-1)+1+m-x)*ln1-ln0*(z-1+m))/(z)
            ln1=ln2
            ln0=ln1
        enddo
        ln=ln2
    end if
endfunction

real*8 function Alag(p,m,i,rn,phin)
real(8) :: ln
real(8), intent(in) :: rn,phin
integer, intent(in) :: p,m,i
integer :: factorial
    Alag=0
    if (m.eq.0) then
    Alag=2*exp(-rn**2)*ln(p,0,2*rn**2)*(1/sqrt(2*pi))
    elseif (i.eq.0) then
    Alag=2*((2*rn**2)**(m/2))*sqrt(factorial(p)*(1./factorial(p+m)))*exp(-rn**2)*ln(p,m,2*(rn**2))*(1/sqrt(pi))*sin(m*phin)
    elseif (i.eq.1) then
    Alag=2*((2*rn**2)**(m/2))*sqrt(factorial(p)*(1./factorial(p+m)))*exp(-rn**2)*ln(p,m,2*(rn**2))*(1/sqrt(pi))*cos(m*phin)
    end if
endfunction
