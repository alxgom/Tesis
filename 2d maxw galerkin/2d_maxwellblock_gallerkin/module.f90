 !*********************************************************************************
  MODULE constants
    implicit none
    real(8), parameter :: a=2d0
    real(8), parameter :: gperp=10**8d0 !#gamma perpendicular, loss rate
    real(8), parameter :: tempscale=1*(10.**6)/gperp !#scale to micro seconds
    real(8), parameter :: wscale=1000*gperp/(10.**6) !#scale frequency to khz
    real(8), parameter :: mu=0.25d0/10**4, Dphi0=0.0d0
    real(8), parameter :: k=0.9*10**7d0/gperp, g=2.5*10.**4/gperp, D0=a*k/mu!, w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale, w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale, atest=D0*mu/k,
    real(8) :: d=1.0d0
    !real(8) :: m=0.02d0
    real(8) :: wf=0.0045d0
    real(8) :: w_res, w!, atest
    integer, parameter :: savefile=1


! TOL : Tolerance for functions and Romberg integration
! NRAD: Number of points used in radial Gauss quadratures
! NANG: Number of points used in angular Gauss quadratures
!
      DOUBLE COMPLEX, PARAMETER   :: IM = (0.,1.)
      real(8), PARAMETER :: pi = 3.141592653589793d0
      DOUBLE PRECISION, PARAMETER :: TOL = 1.d-12
      DOUBLE PRECISION, PARAMETER :: TINY = 1.d-12
      INTEGER, PARAMETER          :: NRAD = 87
      INTEGER, PARAMETER          :: NANG = 87
      SAVE

    contains

    subroutine comparams()
    !'''parameters to compare with the results'''
        w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale !#resonance frequency
        !atest=D0*mu/k
        w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale !#Relaxation oscilations frequency
    end subroutine

    subroutine saveparams()
        if (savefile.eq.1) then
            open (1,file="scales.in",form='unformatted')
                write(1)  wf*wscale, Dphi0, w_res , k, mu, d, g, D0, a, wf, wscale, tempscale
            close (1)
        endif
    end subroutine

  END MODULE constants
 !************************************************************************************
  MODULE rungekutta
!
! ord: order of the Runge-Kutta time integration
!
      INTEGER, PARAMETER :: ord = 4
      SAVE

  END MODULE
  !**************************************************************************************

 module funcs
    implicit none
    public :: linspace

    contains

    subroutine linspace(x,x_start, x_end, x_len,dir)
    !******************************************************************************
    !linearly spaced array named x, from x_start value to x_end, with #x_len values.
    !dir=1 ---> from min to max.    dir=2 ---> from max to min.
    !******************************************************************************
        real(8), dimension(:), intent(inout) :: x
        real(8) :: x_start, x_end
        integer :: x_len, i
        real(8) :: dx
        integer :: dir

        dx=(x_end - x_start)/(x_len-1)
        if (dir.eq.1) then
            do i=1,x_len
                x(i)=x_start+(i-1)*dx
            end do
        end if
        if (dir.eq.2) then
            do i=1,x_len
                x(i)=x_end-(i-1)*dx
            end do
        end if
        !x(1:x_len)=[(x_start+(i-1)*dx)),i=1,x_len]
    end subroutine

subroutine Laguerre(p,m,rho,L)

    real(8),dimension(:), intent(in) :: rho
    real(8),allocatable, dimension(:,:),intent(out) :: L
    real(8),dimension(size(rho)) :: L0,L1,L2
    integer :: i
    integer, intent(in) :: p,m
    integer :: dimtot

    dimtot=0
    do i=0,p
        dimtot=2*i+1
    end do

    allocate(L(size(rho),dimtot))


end subroutine

subroutine sln(xout,p,m,xin)
! The recursion relation is
!(n + 1) Ln+l(x) = (2n + 1 - x +m) Ln(x) - (n+m)Ln(x)   (3.9.7)
!where L0(x) = 1
!L1(x) = 1 +m+x
    real(8) :: lntemp
    real(8),intent(in) :: xin
    integer,intent(in) :: m,p
    real(8),intent(out) :: xout
    integer :: z
    real(8) :: ln0,ln1,ln2

    if (p.eq.0) then
        ln0=1
        xout=ln0
    elseif (p.eq.1) then
        ln1=1+m+xin
        xout=ln1
    elseif (p.gt.1) then
        ln0=1
        ln1=1+m+xin
        do z=2,p
            ln2=((2*(z-1)+1+m-xin)*ln1-ln0*(z-1+m))/(z)
            ln1=ln2
            ln0=ln1
        enddo
        xout=ln2
    end if
endsubroutine


!*****************************************************************
       SUBROUTINE zerobes(n,x1,x2,root)
!-----------------------------------------------------------------
!!
!! Computes the roots of the spherical Besselj functions between
!! x1 and x2. Uses a combination of bisection and Newton-Raphson
!! to find the zeros. The error in the root is controled by the
!! variable TOL in the module constants.
!!
!! Parameters
!!    n   : order of the Bessel function (>=0)
!!    x1  : beginning of the interval
!!    x2  : end of the interval
!!    root: at the output contains the root
!!
!      USE constants
!      IMPLICIT NONE
!
      DOUBLE PRECISION, intent(IN)  :: x1,x2
      DOUBLE PRECISION, intent(OUT) :: root
!      DOUBLE PRECISION              :: df,dx,dxold,f,fh
!      DOUBLE PRECISION              :: fl,temp,xh,xl,order
      INTEGER, INTENT(IN)           :: n
!      INTEGER, PARAMETER            :: MAXIT = 100
!      INTEGER                       :: j
!
!      order = n+0.5d0
!!      CALL sln(x1,order,fl,df)
! !     CALL sln(x2,order,fh,df)
!      IF ((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) THEN
!         OPEN(1,file='error.log',position='append')
!         WRITE(1,*) 'zerobes: interval contains no change of sign'
!         CLOSE(1)
!      ENDIF
!      IF (ABS(fl).eq.0.) THEN
!         root = x1
!         RETURN
!      ELSE IF (ABS(fh).eq.0.) THEN
!         root = x2
!         RETURN
!      ELSE IF(fl.lt.0.) THEN
!         xl = x1
!         xh = x2
!      ELSE
!         xh = x1
!         xl = x2
!      ENDIF
!      root = .5*(x1+x2)
!      dxold = ABS(x2-x1)
!      dx = dxold
!  !    CALL sln(root,order,f,df)
! Z1 : DO j = 1,MAXIT
!         IF (((root-xh)*df-f)*((root-xl)*df-f).gt.0. &
!            .or. abs(2.*f).gt.abs(dxold*df)) THEN
!            dxold = dx
!            dx = 0.5*(xh-xl)
!            root = xl+dx
!            IF (xl.eq.root) RETURN
!         ELSE
!            dxold = dx
!            dx = f/df
!            temp = root
!            root = root-dx
!            IF (temp.eq.root) RETURN
!         ENDIF
!         IF (abs(dx).lt.TOL) RETURN
!         CALL sln(root,order,f,df)
!         IF (f.lt.0.) THEN
!            xl = root
!         ELSE
!            xh = root
!         ENDIF
!      ENDDO Z1
!      OPEN(1,file='error.log',position='append')
!      WRITE(1,*) 'zerobes: root exceeding maximum iterations'
!      CLOSE(1)
!
!      RETURN
      END SUBROUTINE zerobes
!
!
end module
