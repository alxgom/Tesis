program main
    USE constants
    USE rungekutta
    use my_lib
    use funcs

    IMPLICIT NONE
    !************************************************
    !set variables for main program

    real(8) :: intime
    real(8) :: timeinit=0.d0
    real(8), allocatable, dimension(:) :: time
    real(8), dimension(9) :: yinit
    integer :: numstep
    integer :: indexkeep

    real(8), allocatable, dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), allocatable, dimension(:) :: intensity_x2, intensity_y2, intensity

    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1


    INTEGER NMAX,NSTPMX
    PARAMETER (NMAX=9,NSTPMX=20000)! Maximum number of functions and
    !maximum number of values to
    !be stored.
    REAL*8 xx(NSTPMX),y(NMAX,NSTPMX)
    COMMON /path/ xx,y
    !************************************************

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)

    intime=500.*25*gperp/10**6                   !normalized integration time
    call initial(intime,numstep, timeinit , time)!sets numstep, and time array
    numkeep=numstep*10/25                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1

    if (debug.eq.1) then
        print*, indexkeep,'<', numstep
        print*, 'shape time: ', shape(time), '=', 'number of time steps: ', numstep
    end if

    !***********************************************************************************************************
    !newton integration
    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
    call newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)

    allocate(intensity_x2(numkeep),intensity_y2(numkeep),intensity(numkeep))

    print*, size(exr(indexkeep:numstep)), '=', numkeep
    print*, indexkeep+numkeep-1,'=',numstep
    print*, 'indexkeep',indexkeep

    do z=indexkeep,numstep,1
        intensity_x2(z-indexkeep+1)=exr(z)**2+exi(z)**2
        intensity_y2(z-indexkeep+1)=eyr(z)**2+eyi(z)**2
        intensity(z-indexkeep+1)=sqrt(intensity_x2(z-indexkeep+1)+intensity_y2(z-indexkeep+1))
    enddo

    if (debug.eq.1) then
        print*,
        print*, 'size intensity: ', size(intensity), '=', 'numkeep:', numkeep
    end if

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) time(indexkeep: numstep)
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) exr(indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='pop_euler.in',form='unformatted')
        WRITE(1) pop(indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensity
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensity_x2
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensity_y2
    CLOSE(1)
    !**************************************************************************************************************

    call local_max(intensity,numkeep,3,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    if (debug.eq.1) then
        print*,
        print*, 'integration time=', intime*tempscale, 'microseconds'
        print*, 'integration steps:  numstep=', numstep
    end if

   ! call neararray(yinit,numstep,time,dt)



   !********************************
!   function f(x,t)
!        real*8,intent(in) :: x
!        real*8,intent(in) :: t
!        real*8 :: ka
!
!        f=x*t !test function
!    endfunction

    !**************************************


   ! call rk4mio(yinit, time, numstep,f, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
    call rkdumb(yinit,9,0d0,intime,20000,derivs)
    print*,xx
end


    !'''parameters for normalization'''
!    a=2d0
!    gperp=10**8d0 !#gamma perpendicular, loss rate
!    tempscale=1*(10.**6)/gperp !#scale to micro seconds
!    wscale=1000*gperp/(10.**6) !#scale frequency to khz
!
!    !'''parameters for the equation'''
!    k=0.9*10**7d0/gperp !#normalized loss rate
!    mu=0.25d0/10**4 !#g
!    Dphi0=0.0d0 !#phase shift [-pi,pi]
!    d=1.0d0 !#detuning
!    g=2.5*10.**4/gperp !#*((2*pi)**2) #sigma parallel, normalized loss rate
!    D0=a*k/mu !#Pump
!    m=0.02d0 !#modulation amplitud [0,1]
!    wf=0.00420d0
