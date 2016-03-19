!program swipe_w
!
!    USE constants
!    USE rungekutta
!    use my_lib
!    use funcs
!
!    IMPLICIT NONE
!
!    real(8) :: intime
!    real(8) :: timeinit=0.d0
!    real(8), allocatable, dimension(:) :: time
!    real(8), dimension(9) :: yinit
!    integer :: numstep
!    integer :: numkeep
!    integer :: indexkeep
!    real(8), allocatable, dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
!    real(8), allocatable, dimension(:) :: intensity_x2, intensity_y2, intensity
!
!    real(8), parameter :: wfmin=0.00350, wfmax=0.00410
!    integer, parameter :: len_wfn=3
!    real(8), dimension (1:len_wfn) :: wfn
!    integer(4), allocatable, dimension(:) :: indexdata
!    real(8), allocatable, dimension(:) :: peakdata
!    real(8), allocatable, dimension(:) :: wdata
!    integer :: i
!    integer :: z
!    integer, parameter :: debug=1
!
!    call comparams()                             !parameters to compare with the expected solutions
!    call saveparams()                            !saves the used parameters to a bin file, to be read by python.
!
!    yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
!
!    intime=500.*25*gperp/10**6                   !normalized integration time
!    call initial(intime,numstep, timeinit , time)!sets numstep, and time array
!    numkeep=numstep*10/25                        !number of steps to keep on file(transitory)
!    indexkeep=numstep-numkeep+1
!
!    call linspace(wfn,wfmin,wfmax,len_wfn,2)    !i make wfn, a linearly spaced freq vector of dim (1,3) for a swipe
!
!    if (debug.eq.1) then
!    print*, 'frequency swipe steps:  shape wfn=', shape(wfn)
!    endif
!
!
!    yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
!
!    intime=500.*25*gperp/10**6                   !normalized integration time
!    call initial(intime,numstep, timeinit , time)!sets numstep, and time array
!    numkeep=numstep*10/25                        !number of steps to keep on file(transitory)
!    indexkeep=numstep-numkeep+1
!
!
!    do i=1,len_wfn
!        wf=wfn(i)
!        if (i.eq.0) then
!
!            yinit=yinit
!            timeinit=0d0
!
!    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
!            call newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
!
!            allocate(intensity_x2(numkeep),intensity_y2(numkeep),intensity(numkeep))
!            do z=indexkeep,numstep,1   !set the first intensity
!                intensity_x2(z-indexkeep+1)=exr(z)**2+exi(z)**2
!                intensity_y2(z-indexkeep+1)=eyr(z)**2+eyi(z)**2
!                intensity(z-indexkeep+1)=sqrt(intensity_x2(z-indexkeep+1)+intensity_y2(z-indexkeep+1))
!            enddo
!
!        else
!    !        yinit=(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
!    !        timeinit= time(numstep)
!    !        call initial(intime,numstep, timeinit , time)                !sets numstep, and time array
!!    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
!!    call newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
!!
!!    allocate(intensity_x2(numkeep),intensity_y2(numkeep),intensity(numkeep))
!!
!!    print*, size(exr(indexkeep:numstep)), '=', numkeep
!!    print*, indexkeep+numkeep-1,'=',numstep
!!    print*, 'indexkeep',indexkeep
!!
!!    do z=indexkeep,numstep,1
!!        !print*,'z=',z
!!        intensity_x2(z-indexkeep+1)=exr(z)**2+exi(z)**2
!!        intensity_y2(z-indexkeep+1)=eyr(z)**2+eyi(z)**2
!!        intensity(z-indexkeep+1)=sqrt(intensity_x2(z-indexkeep+1)+intensity_y2(z-indexkeep+1))
!     !enddo
!    endif
!    enddo
!    call local_max_bif(intensity,numkeep,3,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)
!
!
!!        end if
!!    end do
!!    allocate(indexdata)
!!    OPEN(1, file='peak_index.in', ACCESS='STREAM', FORM='UNFORMATTED')
!!        READ (1) indexdata
!!    CLOSE(1)
!!    OPEN(10, file='peak_index_fort.in', FORM='UNFORMATTED')
!!        write(10) indexdata
!!    CLOSE(10)
!!    allocate(peakdata)
!!    OPEN(1, file='peak.in', ACCESS='STREAM', FORM='UNFORMATTED')
!!        READ (1) peakdata
!!    CLOSE(1)
!!    OPEN(10, file='peak_fort.in', FORM='UNFORMATTED')
!!        write(10) peakdata
!!    CLOSE(10)
!!
!end program
!!****************************************************************************************
!
!
!
!
!
!
!
!
!
!
!
