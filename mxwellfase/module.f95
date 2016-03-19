!=================================================================
! MODULES for SPHERE code
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu
!=================================================================

!=================================================================

  MODULE resolution
!
! q: maximum number of zeros in the radial direction
! l: maximum number of l in the spherical harmonics
!
      INTEGER, PARAMETER :: q = 5
      INTEGER, PARAMETER :: l = 5
      SAVE

  END MODULE

!=================================================================

  MODULE rungekutta
!
! ord: order of the Runge-Kutta time integration
!
      INTEGER, PARAMETER :: ord = 4
      SAVE

  END MODULE

!=================================================================

  MODULE constants
    implicit none
    real(8), parameter :: a=2d0
    real(8), parameter :: gperp=10**8d0 !#gamma perpendicular, loss rate
    real(8), parameter :: tempscale=1*(10.**6)/gperp !#scale to micro seconds
    real(8), parameter :: wscale=1000*gperp/(10.**6) !#scale frequency to khz
    real(8), parameter :: mu=0.25d0/10**4, Dphi0=0.0d0, d=1.0d0
    real(8), parameter :: k=0.9*10**7d0/gperp, g=2.5*10.**4/gperp, D0=a*k/mu, m=0.02d0,wf=0.00420d0!, w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale, w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale, atest=D0*mu/k,


    SAVE

  END MODULE constants
!=================================================================


  MODULE random
      CONTAINS
       DOUBLE PRECISION FUNCTION randu(idum)
!
! Uniform distributed random numbers between -1 and
! 1. The seed idum must be between 0 and the value
! of mask

       INTEGER, PARAMETER  :: iq=127773,ir=2836,mask=123459876
       INTEGER, PARAMETER  :: ia=16807,im=2147483647
       INTEGER             :: idum
       INTEGER             :: k
       DOUBLE PRECISION, PARAMETER :: am=1.d0/im

       idum = ieor(idum,mask)
       k = idum/iq
       idum = ia*(idum-k*iq)-ir*k
       IF (idum.lt.0) idum = idum+im
       randu = am*idum
       randu = (randu-.5d0)*2.d0
       idum = ieor(idum,mask)

       END FUNCTION randu

 END MODULE random


!=================================================================

!
!  MODULE functions
!!
!! Defines functions to be integrated using Romberg integration
!! or Gauss quadratures. The data type 'arguments' allows
!! optional arguments to be passed to these functions.
!!
!      TYPE ARGUMENTS
!         INTEGER          :: l1,l2,l3
!         INTEGER          :: m1,m2,m3
!         DOUBLE PRECISION :: lam1,lam2,lam3
!      END TYPE ARGUMENTS
!!
!   CONTAINS
!
!      DOUBLE PRECISION FUNCTION TRIPLER(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2
!         DOUBLE PRECISION             :: sj3,sjp
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp)
!         CALL sphbes(arg%l3,arg%lam3*x,sj3,sjp)
!         TRIPLER = sj1*sj2*sj3*x
!      RETURN
!      END FUNCTION TRIPLER
!
!      DOUBLE PRECISION FUNCTION TRIPLERD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2
!         DOUBLE PRECISION             :: sj3,sjp
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp)
!         CALL sphbes(arg%l3,arg%lam3*x,sj3,sjp)
!         sjp = sjp*arg%lam3*x+sj3
!         TRIPLERD = sj1*sj2*sjp
!      RETURN
!      END FUNCTION TRIPLERD
!
!      DOUBLE PRECISION FUNCTION TRIPLERDD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2,sj3
!         DOUBLE PRECISION             :: sjp2,sjp3
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp2)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp2)
!         CALL sphbes(arg%l3,arg%lam3*x,sj3,sjp3)
!         sjp2 = sjp2*arg%lam2*x+sj2
!         sjp3 = sjp3*arg%lam3*x+sj3
!         TRIPLERDD = sj1*sjp2*sjp3/x
!      RETURN
!      END FUNCTION TRIPLERDD
!
!      DOUBLE PRECISION FUNCTION TRIPLEX(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2,plg3
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL plgndr(arg%l2,arg%m2,x,plg2,2)
!         CALL plgndr(arg%l3,arg%m3,x,plg3,2)
!         TRIPLEX = plg1*plg2*plg3/(1-x**2)
!      RETURN
!      END FUNCTION TRIPLEX
!
!      DOUBLE PRECISION FUNCTION TRIPLEXD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2,plg3
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL plgndr(arg%l2,arg%m2,x,plg2,2)
!         CALL legder(arg%l3,arg%m3,x,plg3,2)
!         TRIPLEXD = plg1*plg2*plg3/(1-x**2)
!      RETURN
!      END FUNCTION TRIPLEXD
!
!      DOUBLE PRECISION FUNCTION TRIPLEXDD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2,plg3
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL legder(arg%l2,arg%m2,x,plg2,2)
!         CALL legder(arg%l3,arg%m3,x,plg3,2)
!         TRIPLEXDD = plg1*plg2*plg3/(1-x**2)
!      RETURN
!      END FUNCTION TRIPLEXDD
!
!      DOUBLE PRECISION FUNCTION DOUBLER(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2
!         DOUBLE PRECISION             :: sjp
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp)
!         DOUBLER = sj1*sj2*x**2
!      RETURN
!      END FUNCTION DOUBLER
!
!      DOUBLE PRECISION FUNCTION DOUBLERD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2
!         DOUBLE PRECISION             :: sjp
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp)
!         sjp = sjp*arg%lam2*x+sj2
!         DOUBLERD = sj1*sjp*x
!      RETURN
!      END FUNCTION DOUBLERD
!
!      DOUBLE PRECISION FUNCTION DOUBLERDD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2
!         DOUBLE PRECISION             :: sjp1,sjp2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp1)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp2)
!         sjp1 = sjp1*arg%lam1*x+sj1
!         sjp2 = sjp2*arg%lam2*x+sj2
!         DOUBLERDD = sjp1*sjp2
!      RETURN
!      END FUNCTION DOUBLERDD
!
!      DOUBLE PRECISION FUNCTION DOUBLEQ(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2
!         DOUBLE PRECISION             :: sjp
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp)
!         DOUBLEQ = sj1*sj2*x
!      RETURN
!      END FUNCTION DOUBLEQ
!
!      DOUBLE PRECISION FUNCTION DOUBLEQD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: sj1,sj2
!         DOUBLE PRECISION             :: sjp
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL sphbes(arg%l1,arg%lam1*x,sj1,sjp)
!         CALL sphbes(arg%l2,arg%lam2*x,sj2,sjp)
!         sjp = sjp*arg%lam2*x+sj2
!         DOUBLEQD = sj1*sjp
!      RETURN
!      END FUNCTION DOUBLEQD
!
!      DOUBLE PRECISION FUNCTION DOUBLEW(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL plgndr(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEW = plg1*plg2/SQRT(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEW
!
!      DOUBLE PRECISION FUNCTION DOUBLEWD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL legder(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEWD = plg1*plg2/SQRT(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEWD
!
!      DOUBLE PRECISION FUNCTION DOUBLEWDD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL legder(arg%l1,arg%m1,x,plg1,2)
!         CALL legder(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEWDD = plg1*plg2/SQRT(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEWDD
!
!      DOUBLE PRECISION FUNCTION DOUBLEX(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL plgndr(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEX = x*plg1*plg2/(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEX
!
!      DOUBLE PRECISION FUNCTION DOUBLEXD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL legder(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEXD = x*plg1*plg2/(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEXD
!
!      DOUBLE PRECISION FUNCTION DOUBLEXDD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL legder(arg%l1,arg%m1,x,plg1,2)
!         CALL legder(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEXDD = x*plg1*plg2/(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEXDD
!
!      DOUBLE PRECISION FUNCTION DOUBLEY(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL plgndr(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEY = x*plg1*plg2/SQRT(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEY
!
!
!      DOUBLE PRECISION FUNCTION DOUBLEYD(x,arg)
!         DOUBLE PRECISION, INTENT(IN) :: x
!         DOUBLE PRECISION             :: plg1,plg2
!         TYPE(ARGUMENTS), INTENT(IN)  :: arg
!         CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!         CALL legder(arg%l2,arg%m2,x,plg2,2)
!         DOUBLEYD = x*plg1*plg2/SQRT(1-x**2)
!      RETURN
!      END FUNCTION DOUBLEYD
!
!   !   DOUBLE PRECISION FUNCTION DOUBLEZ(x,arg)
!   !      DOUBLE PRECISION, INTENT(IN) :: x
!   !!      DOUBLE PRECISION             :: plg1,plg2
!    !!     TYPE(ARGUMENTS), INTENT(IN)  :: arg
!  !       CALL plgndr(arg%l2,arg%m2,x,plg2,2)
!     !    CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!     !    DOUBLEZ = plg1*plg2
!     ! RETURN
!     ! END FUNCTION DOUBLEZ
!
!     DOUBLE PRECISION FUNCTION DOUBLEZD(x,arg)
!        DOUBLE PRECISION, INTENT(IN) :: x
!       DOUBLE PRECISION             :: plg1,plg2
!       TYPE(ARGUMENTS), INTENT(IN)  :: arg
!      CALL plgndr(arg%l1,arg%m1,x,plg1,2)
!     CALL legder(arg%l2,arg%m2,x,plg2,2)
!      DOUBLEZD = plg1*plg2
!      RETURN
!      END FUNCTION DOUBLEZD
!
!  END MODULE functions
!=================================================================
!!
