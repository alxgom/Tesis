module funcs
    implicit none

    contains

!****************************************************************************************

subroutine local_max_bif_v1(valdata,size_data,compare_size,count_peak,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    implicit none
    integer(4), intent(in) :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer, intent(in) :: debug
    integer, intent(in) :: size_data
    real(8), dimension(size_data), intent(in) :: valdata
    integer(4) :: j
    real(8), dimension(compare_size) :: localdata
    real(8) :: peak
    integer(4), intent(inout) :: count_peak
    integer(4) :: peak_index
    integer(4),dimension(1) :: maxind
    integer(4), allocatable, dimension(:) :: data
    integer(4) :: midpoint

    if (debug.eq.1) then
        print*,
        print*, '*************************************'
        print*, 'debug info for local_max_bif:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2
    !count_peak=0
    OPEN(1,file='peak_index.in',form='unformatted',status='replace',access='stream')
    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            peak=localdata(midpoint)
            count_peak=count_peak+1
            peak_index=j+midpoint-1
            WRITE(1) peak_index
            WRITE(2) peak
        end if
    end do
    CLOSE(1)

    if (debug.eq.1) then
        print*, 'last peak index: ', peak_index
        print*, 'amount of peaks', count_peak
        print*, '*************************************'
        print*,
    endif

    !-need  to add a way to test if the new peak is repeated or not. then if not, add no the file the new one.
    !-should i save the peaks? or only the index?.. if i have the index the i still have to gou though all the array.
    !if i save the peaks i dont have to spend memory on the arrays
    !- A cheaper way would be to keep only 5 steps at a time of the array, and look there for maxima.
    !So i don't have to use so much memory on the fields.

end subroutine

!****************************************************************************************

subroutine local_max_bif_v2(valdata,size_data,compare_size,count_peak,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    use constants !changing parameters for the bifurcation
    implicit none

    integer(4), intent(in) :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer(4), intent(in) :: size_data
    real(8), dimension(size_data), intent(in) :: valdata
    integer(4), intent(inout) :: count_peak
    real(8), dimension(compare_size) :: localdata
    real(8) :: peak
    integer(4) :: peak_index
    integer(4) :: midpoint
    integer(4) :: j
    integer(4),dimension(1) :: maxind
    integer(4),  intent(in) :: debug

    if (debug.eq.1) then
        print*,
        print*, '*************************************'
        print*, 'debug info for local_max_bif:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2

    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            peak=localdata(midpoint)
            count_peak=count_peak+1
            peak_index=j+midpoint-1
            WRITE(2) peak
            WRITE(3) wf
            WRITE(4) m
        end if
    end do

    if (debug.eq.1) then
        print*, 'last peak index: ', peak_index
        print*, 'amount of peaks', count_peak
        print*, '*************************************'
        print*,
    endif

    return
    !-need  to add a way to test if the new peak is repeated or not. then if not, add no the file the new one.
    !-should i save the peaks? or only the index?.. if i have the index the i still have to gou though all the array.
    !if i save the peaks i dont have to spend memory on the arrays
    !- A cheaper way would be to keep only 5 steps at a time of the array, and look there for maxima.
    !So i don't have to use so much memory on the fields.

end subroutine

!****************************************************************************************

subroutine local_max_bif_cos(valdata,valdata2,size_data,compare_size,count_peak,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    use constants !changing parameters for the bifurcation
    implicit none

    integer(4), intent(in) :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer(4), intent(in) :: size_data
    real(8), dimension(size_data), intent(in) :: valdata,valdata2
    integer(4), intent(inout) :: count_peak
    real(8), dimension(compare_size) :: localdata
    real(8) :: peak
    integer(4) :: peak_index
    integer(4) :: midpoint
    integer(4) :: j
    integer(4),dimension(1) :: maxind
    integer(4),  intent(in) :: debug

    if (debug.eq.1) then
        print*,
        print*, '*************************************'
        print*, 'debug info for local_max_bif_coseno:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2

    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            !peak=localdata(midpoint)
            count_peak=count_peak+1
            peak_index=j+midpoint-1
            peak=valdata2(peak_index)
            WRITE(2) peak
        end if
    end do

    if (debug.eq.1) then
        print*, 'last peak index: ', peak_index
        print*, 'amount of peaks', count_peak
        print*, '*************************************'
        print*,
    endif

    return
    !-need  to add a way to test if the new peak is repeated or not. then if not, add no the file the new one.
    !-should i save the peaks? or only the index?.. if i have the index the i still have to gou though all the array.
    !if i save the peaks i dont have to spend memory on the arrays
    !- A cheaper way would be to keep only 5 steps at a time of the array, and look there for maxima.
    !So i don't have to use so much memory on the fields.

end subroutine

!****************************************************************************************

subroutine derivs(x,y,dydx)
    use constants
    implicit none
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

    dydx(1)=(-k*y(1)+mu*y(5))
    dydx(2)=(-k*y(2)+mu*y(6))
    dydx(3)=(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m*cos(wf*x))
    dydx(4)=(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m*cos(wf*x))
    dydx(5)=(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine

!****************************************************************************************

subroutine derivsn(x,n,y,dydx)
    use constants
    implicit none
    integer*4,intent(inout) :: n
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

    dydx(1)=(-k*y(1)+mu*y(5))
    dydx(2)=(-k*y(2)+mu*y(6))
    dydx(3)=(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m*cos(wf*x))
    dydx(4)=(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m*cos(wf*x))
    dydx(5)=(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine

!****************************************************************************************

subroutine derivsn_rev(x,n,y,dydx)
    use constants
    implicit none
    integer*4,intent(inout) :: n
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

    dydx(1)=-(-k*y(1)+mu*y(5))
    dydx(2)=-(-k*y(2)+mu*y(6))
    dydx(3)=-(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m*cos(wf*x))
    dydx(4)=-(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m*cos(wf*x))
    dydx(5)=-(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=-(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=-(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=-(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=-(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine


end module


module my_lib
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
end module

!****************************************************************************************

module not_in_use
    implicit none

    contains

    subroutine derivs_var(x,y,dydx)
    use constants
    implicit none
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

    dydx(1)=(-k*y(1)+mu*y(5))
    dydx(2)=(-k*y(2)+mu*y(6))
    dydx(3)=(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m*cos(wf*x))
    dydx(4)=(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m*cos(wf*x))
    dydx(5)=(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine

    subroutine neararray(yinit,numstep,time)
!
    use constants
    implicit none
    real(8), allocatable, dimension(:) :: exr_near, exi_near,eyr_near,eyi_near,rxr_near,rxi_near,ryr_near,ryi_near,pop_near  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), dimension(9) ,intent(in) :: yinit
    integer, intent(in) :: numstep
    real(8), dimension(numstep), intent(in) :: time
    real(8), allocatable, dimension(:) :: intensity_x2_near, intensity_y2_near, intensity_near
    integer :: i
!    real(8),intent(in) :: dt
    real(8) :: temp

    exr_near(1)=yinit(1)+10*epsilon(yinit(1))
    exi_near(1)=yinit(2)
    eyr_near(1)=yinit(3)
    eyi_near(1)=yinit(4)
    rxr_near(1)=yinit(5)
    ryi_near(1)=yinit(6)
    ryr_near(1)=yinit(7)
    ryi_near(1)=yinit(8)
    pop_near(1)=yinit(9)

    do i=1, numstep
        exr_near(i+1)=exr_near(i)+dt*(-k*exr_near(i)+mu*rxr_near(i))
        exi_near(i+1)=exi_near(i)+dt*(-k*exi_near(i)+mu*rxi_near(i))
        eyr_near(i+1)=eyr_near(i)+dt*(-k*eyr_near(i)+mu*ryr_near(i))-eyi_near(i)*(Dphi0+m*cos(wf*time(i)))
        eyi_near(i+1)=eyi_near(i)+dt*(-k*eyi_near(i)+mu*ryi_near(i))+eyr_near(i)*(Dphi0+m*cos(wf*time(i)))
        rxr_near(i+1)=rxr_near(i)+dt*(-(rxr_near(i)-d*rxi_near(i))+exr_near(i)*pop_near(i))
        rxi_near(i+1)=rxi_near(i)+dt*(-(rxi_near(i)+d*rxr_near(i))+exi_near(i)*pop_near(i))
        ryr_near(i+1)=ryr_near(i)+dt*(-(ryr_near(i)-d*ryi_near(i))+eyr_near(i)*pop_near(i))
        ryi_near(i+1)=ryi_near(i)+dt*(-(ryi_near(i)+d*ryr_near(i))+eyi_near(i)*pop_near(i))
        temp=(-g*(pop_near(i)-D0+(exr_near(i)*rxr_near(i)+exi_near(i)*rxi_near(i)+eyr_near(i)*ryr_near(i)+eyi_near(i)*ryi_near(i))))
        pop_near(i+1)=pop_near(i)+dt*temp
    enddo

    allocate(intensity_x2_near(numstep),intensity_y2_near(numstep),intensity_near(numstep))
    do i=1,numstep
        intensity_x2_near(i)=exr_near(i)**2+exi_near(i)**2
        intensity_y2_near(i)=eyr_near(i)**2+eyi_near(i)**2
        intensity_near(i)=sqrt(intensity_x2_near(i)+intensity_y2_near(i))
    enddo

    OPEN(1,file='intensity_near.in',form='unformatted')
        WRITE(1) intensity_near
    CLOSE(1)

end subroutine

    subroutine newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
    !Newton method integration.

    use constants
    implicit none

    real(8), dimension(numstep), intent(inout) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), dimension(9), intent(in) :: yinit
    real(8) :: dtn
    integer, intent(in) :: numstep
    real(8), dimension(numstep), intent(in) :: time
    integer :: i

    integer, optional :: debug

    dtn=(time(numstep)-time(1))/(numstep-1)

    if (debug.eq.1) then
        print*,
        print*, '******************************************************************'
        print*, 'subroutine debug info for newtonint:'
        print*,
        print*, 'time step:  dt=', dt*tempscale
        print*, 'time steps:', size(time)
        print*, 'integration from ', time(1), 'to', time(size(time))
        !print*, 'Initial values:', yinit
    end if

    exr(1)=yinit(1)
    exi(1)=yinit(2)
    eyr(1)=yinit(3)
    eyi(1)=yinit(4)
    rxr(1)=yinit(5)
    ryi(1)=yinit(6)
    ryr(1)=yinit(7)
    ryi(1)=yinit(8)
    pop(1)=yinit(9)

    do i=1,numstep
        exr(i+1)=exr(i)+dtn*(-k*exr(i)+mu*rxr(i))
        exi(i+1)=exi(i)+dtn*(-k*exi(i)+mu*rxi(i))
        eyr(i+1)=eyr(i)+dtn*(-k*eyr(i)+mu*ryr(i))-eyi(i)*(Dphi0+m*cos(wf*time(i)))
        eyi(i+1)=eyi(i)+dtn*(-k*eyi(i)+mu*ryi(i))+eyr(i)*(Dphi0+m*cos(wf*time(i)))
        rxr(i+1)=rxr(i)+dtn*(-(rxr(i)-d*rxi(i))+exr(i)*pop(i))
        rxi(i+1)=rxi(i)+dtn*(-(rxi(i)+d*rxr(i))+exi(i)*pop(i))
        ryr(i+1)=ryr(i)+dtn*(-(ryr(i)-d*ryi(i))+eyr(i)*pop(i))
        ryi(i+1)=ryi(i)+dtn*(-(ryi(i)+d*ryr(i))+eyi(i)*pop(i))
        pop(i+1)=pop(i)+dtn*(-g*(pop(i)-D0+(exr(i)*rxr(i)+exi(i)*rxi(i)+eyr(i)*ryr(i)+eyi(i)*ryi(i))))
    enddo

    if (debug.eq.1) then
        !print*, 'variables size ', exr(1:20)
        print*, '******************************************************************'
    end if
end subroutine

    subroutine initial(numstep, timeinit, time)
!**********************************************************************************************
!sets the amount of time steps and the time array to use, from the total time to integrate (intime), and the initial time (timeinit)

    use constants
    use my_lib

    implicit none

!    real(8), intent(in) :: intime
    integer, intent(out) :: numstep
    real(8) :: jump
    real(8), intent(in) :: timeinit
    real(8), allocatable, dimension(:), intent(inout) :: time

    jump=.5d0
    numstep=int(intime/jump)
    allocate(time(numstep))
    call linspace(time,timeinit, timeinit+intime, numstep,1)
end subroutine

    subroutine local_max(valdata,size_data,compare_size,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    implicit none
    integer(4), intent(in) :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer, intent(in) :: debug
    integer, intent(in) :: size_data
    real(8), dimension(size_data), intent(in) :: valdata
    integer(4) :: j
    real(8), dimension(compare_size) :: localdata
    real(8) :: peak
    integer(4) :: count_peak, peak_index
    integer(4), dimension(1) :: maxind
    integer(4), allocatable, dimension(:) :: data
    integer(4) :: midpoint

    if (debug.eq.1) then
        print*,
        print*, '*************************************'
        print*, 'debug info for local_max:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2
    count_peak=0
    OPEN(1,file='peak_index.in',form='unformatted',status='unknown',access='stream')
    OPEN(2,file='peak.in',form='unformatted',status='unknown',access='stream')
    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            peak=localdata(midpoint)
            count_peak=count_peak+1
            peak_index=j+midpoint-1
                WRITE(1) peak_index
                WRITE(2) peak
        end if
    end do
    CLOSE(1)
    CLOSE(2)


    if (debug.eq.1) then
        print*, 'last peak index: ', peak_index
        print*, 'amount of peaks', count_peak
        print*, '*************************************'
        print*,
    endif


    allocate(data(count_peak))
    OPEN(1, file='peak_index.in', ACCESS='STREAM', FORM='UNFORMATTED')
    READ (1) data
    CLOSE(1)
    OPEN(1, file='peak_index_fort.in', FORM='UNFORMATTED')
    write(1) data
    CLOSE(1)
end subroutine

!    subroutine w_swipe_1()! for rk45
!    !bug: no agarra bien el segundo tiempo.
!    use constants
!    implicit none
!
!    real*8 :: t_start, t_stop,t_out,mvar,wvar
!    integer :: i
!    real*8 :: w_start, w_stop, h
!
!    t_start = 0.0D+00
!    t_stop = 500.*40*gperp/10**6
!    !asi como esta ni hace falta iterar en el espacio dos veces, pero si poner el nuevo tiempo despues de cada paso
!    w_start=0.0038
!    w_stop=0.0045
!    h=(w_stop-w_start)/2000d0
!    do i=1,20
!        !w swipe
!        mvar=0.02
!        wvar=w_start+h
!        call maxintrk45f_var(t_start,t_stop,t_out,mvar,wvar)
!        print*, 'i step t_out:', i, '-->', t_out
!        print*, 'i step t_stop:', i, '-->', t_stop
!        t_start=t_out
!        t_stop=t_out+t_stop
!        print*,i
!    end do
!
!end subroutine

    subroutine maxintegeuler()!euler integration
    USE constants
    USE rungekutta
    use my_lib
    use funcs

    IMPLICIT NONE
    !************************************************
    !set variables for main program

 !   real(8) :: intime
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
    !********************************************************************************
    print*, ' -------------------------'
    print*, '|-Euler method program:   |'
    print*, ' -------------------------'

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.
    yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
  !  intime=500.*10*gperp/10**6                   !normalized integration time
    call initial(numstep, timeinit , time)!sets numstep, and time array
    numkeep=numstep*10/25                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1

    if (debug.eq.1) then
        print*, indexkeep,'<', numstep
        print*, 'shape time: ', shape(time), '=', 'number of time steps: ', numstep
        print*, 'intime:', intime
        print*, 'yinit:', yinit
    end if

    !***********************************************************************************************************
    !newton integration
    allocate(exr(numstep),exi(numstep),eyr(numstep),eyi(numstep),rxr(numstep),rxi(numstep),ryr(numstep),ryi(numstep),pop(numstep))
    call newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
    allocate(intensity_x2(numkeep),intensity_y2(numkeep),intensity(numkeep))
!
    print*, size(exr(indexkeep:numstep)), '=', numkeep
    print*, indexkeep+numkeep-1,'=',numstep
    print*, 'indexkeep',indexkeep
!
    do z=indexkeep,numstep,1
        intensity_x2(z-indexkeep+1)=exr(z)**2+exi(z)**2
        intensity_y2(z-indexkeep+1)=eyr(z)**2+eyi(z)**2
        intensity(z-indexkeep+1)=sqrt(intensity_x2(z-indexkeep+1)+intensity_y2(z-indexkeep+1))
    enddo
!
    if (debug.eq.1) then
        print*,
        print*, 'size intensity: ', size(intensity), '=', 'numkeep:', numkeep
    end if
!
    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time_euler.in',form='unformatted')
        WRITE(2) time(indexkeep: numstep)
    CLOSE(2)

    OPEN(1,file='exr_euler.in',form='unformatted')
        WRITE(1) exr(indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='pop_euler.in',form='unformatted')
        WRITE(1) pop(indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='intensity_euler.in',form='unformatted')
        WRITE(1) intensity
    CLOSE(1)

    OPEN(1,file='intensity_x2_euler.in',form='unformatted')
        WRITE(1) intensity_x2
    CLOSE(1)

    OPEN(1,file='intensity_y2_euler.in',form='unformatted')
        WRITE(1) intensity_y2
    CLOSE(1)
    !**************************************************************************************************************
!
    call local_max(intensity,numkeep,1,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    if (debug.eq.1) then
        print*,
        print*, 'integration time=', intime*tempscale, 'microseconds'
        print*, 'integration steps:  numstep=', numstep
    end if

   ! call neararray(yinit,numstep,time,dt)


end subroutine

    subroutine maxintegrk45f()
    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE


    integer ( kind = 4 ), parameter :: neqn = 9 !2

    real ( kind = 8 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    external r8_f2
    real ( kind = 8 ) relerr
    real ( kind = 8 ) t
    real ( kind = 8 ) t_out
    real ( kind = 8 ) t_start
    real ( kind = 8 ) t_stop
    real ( kind = 8 ) y(neqn)
    real ( kind = 8 ) yp(neqn)

    real*8, allocatable, dimension(:,:) :: intensitys, yy
    REAL*8, allocatable, dimension(:) :: xx
    integer :: numstep
    integer :: indexkeep
    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Runge kutta 4/5 felhbelg method:'
    write ( *, '(a)' ) '  Solve a vector equation using R8_RKF45:'
    write ( *, '(a)' ) ' '

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    t_start = 0.0D+00
    t_stop = 500.*40*gperp/10**6

    n_step = int(t_stop)
    allocate(xx(n_step),yy(neqn,n_step)) !mine
    t = 0.0D+00
    t_out = 0.0D+00

    y=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)  !valor inicial de y.

    call derivs ( t, y, yp )

    !  write ( *, '(a)' ) ' '
    !  write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
    !  write ( *, '(a)' ) ' '
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    do i_step = 1, n_step
    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
        + real (          i_step - 1, kind = 8 ) * t_stop ) &
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

    call r8_rkf45 ( derivs, neqn, y, yp, t, t_out, relerr, 10d0**(-4), flag )

    xx(i_step)=t_out !mine
    yy(:,i_step)=y  !mine
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do

    numstep=n_step
    numkeep=numstep*5/40                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1
    allocate(intensitys(3,numkeep))
    do z=indexkeep,numstep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z-indexkeep+1)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z-indexkeep+1)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z-indexkeep+1)=sqrt(intensitys(1,z-indexkeep+1)+intensitys(2,z-indexkeep+1)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', n_step
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) xx(indexkeep: numstep)
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) yy(1,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='pop.in',form='unformatted')
        WRITE(1) yy(9,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensitys(3,:)
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensitys(1,:)
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensitys(2,:)
    CLOSE(1)

    return

end subroutine

    subroutine maxintrk45f_var(t_start,t_stop,t_out)
    USE constants
    use my_lib
    use funcs
    use varparams

    IMPLICIT NONE


    integer ( kind = 4 ), parameter :: neqn = 9 !2

    real ( kind = 8 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    external r8_f2
    real ( kind = 8 ) relerr
    real ( kind = 8 ) t
    real ( kind = 8 ), intent(out) :: t_out
    real ( kind = 8 ), intent(inout) :: t_start
    real ( kind = 8 ), intent(in) :: t_stop
    real ( kind = 8 ) y(neqn)
    real ( kind = 8 ) yp(neqn)

    real*8, allocatable, dimension(:,:) :: intensitys, yy
    REAL*8, allocatable, dimension(:) :: xx
    integer :: numstep
    integer :: indexkeep
    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Runge kutta 4/5 felhbelg method:'
    write ( *, '(a)' ) '  Solve a vector equation using R8_RKF45:'
    write ( *, '(a)' ) ' '

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    n_step = int(t_stop)
    allocate(xx(n_step),yy(neqn,n_step)) !mine
    t = 0.0D+00
    t_out = 0.0D+00

    y=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)  !valor inicial de y.

    call derivs_var ( t, y, yp )

    !  write ( *, '(a)' ) ' '
    !  write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
    !  write ( *, '(a)' ) ' '
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    do i_step = 1, n_step
    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
        + real (          i_step - 1, kind = 8 ) * t_stop ) &
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

    call r8_rkf45 ( derivs, neqn, y, yp, t, t_out, relerr, 10d0**(-4), flag )

    xx(i_step)=t_out !mine
    yy(:,i_step)=y  !mine
    !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do

    numstep=n_step
    numkeep=numstep*5/40                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1
    allocate(intensitys(3,numkeep))
    do z=indexkeep,numstep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z-indexkeep+1)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z-indexkeep+1)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z-indexkeep+1)=sqrt(intensitys(1,z-indexkeep+1)+intensitys(2,z-indexkeep+1)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', n_step
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    return

end subroutine

    subroutine maxrk4f90()
     !17/3/2016
    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE

    integer ( kind = 4 ), parameter :: n = 9
   ! real ( kind = 8 ), parameter :: dt = 1d0
    !external derivsn
    real ( kind = 8 ) t0
    real ( kind = 8 ) t1
    real ( kind = 8 ), parameter :: tmax = 500.*20*gperp/10**6
    real ( kind = 8 ) u1(n)
    real ( kind = 8 ) u0(n)

    real*8, allocatable, dimension(:,:) :: intensitys, yy
    REAL*8, allocatable, dimension(:) :: xx
    integer :: numstep
    integer :: indexkeep
    integer :: numkeep
    integer :: z
    integer, parameter :: debug=1
    integer i_step

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RK4_PRB'
    write ( *, '(a)' ) '  FORTRAN90 version.'
    write ( *, '(a)' ) '  Test the RK4 library.'

    t0 = 0.0D+00
    u0=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
    numstep=int(tmax/dt)
    allocate(xx(numstep),yy(n,numstep))
    i_step=1
    xx(1)=t0 !mine0
    yy(:,1)=u0  !mine


    do
        !
        !  Print (T0,U0).
        !  Stop if we've exceeded TMAX.
        if ( tmax <= t0 ) then
          exit
        end if
        !
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        i_step=i_step+1
        xx(i_step)=t1 !mine
        yy(:,i_step)=u1  !mine
        !  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do

    print*,''
    PRINT*, 'check: dim(xx)=numstep', shape(xx),'=', numstep
    print*,''

    numkeep=numstep*5/20                        !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1
    allocate(intensitys(3,numkeep))
    do z=indexkeep,numstep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z-indexkeep+1)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z-indexkeep+1)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z-indexkeep+1)=sqrt(intensitys(1,z-indexkeep+1)+intensitys(2,z-indexkeep+1)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', numstep
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) xx(indexkeep: numstep)
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) yy(1,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='pop.in',form='unformatted')
        WRITE(1) yy(9,indexkeep: numstep)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensitys(3,:)
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensitys(1,:)
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensitys(2,:)
    CLOSE(1)

    deallocate(xx,yy,intensitys)
    return

end subroutine

    subroutine maxrk4f90_straight()

    !17/3/2016
    !Base script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.
    ! This subroutine is the basis for the other routines used to study the equations.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the number of steps i will keep before the last.)

    !Input <-- dt(time step size), t0, u0(initial field values), tmax
    !Input <-- The parameters used for the normalizations are in Constants module.
    !Input <-- Numkeep(the number of steps i keep for display)

    !Output --> xx(numkeep) the times used in the integration. yy(9,numkeep) the field values.
    !Output --> Intensitys(3,numkeep) The intensity of the Ex field, the Ey field, and the total intensity.

    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE

    integer ( kind = 4 ), parameter :: n = 9 !(number of fields)
   ! real ( kind = 8 ), parameter :: dt = 1d0 !time step
    real ( kind = 8 ) t0 !initial time.
    real ( kind = 8 ) t1
    real ( kind = 8 ), parameter :: tmax = 500.*20*gperp/10**6 ! End time
    real ( kind = 8 ) u1(n)
    real ( kind = 8 ) u0(n) !Initial conditions

    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    REAL*8, allocatable, dimension(:) :: xx                 !Times
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z
    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.
    integer i_step, k_step           !integration step and keep step.

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RK4_Phasemod Maxwell'

    t0 = 0.0D+00
    u0=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
    numstep=int(tmax/dt)
    numkeep=numstep*4/20                         !number of steps to keep on file(transitory)
    allocate(xx(numkeep),yy(n,numkeep))
    indexkeep=numstep-numkeep+1       !set indexkeep
    i_step=1
    k_step=1
    do
        !  Stop if we've exceeded TMAX.
        if ( tmax <= t0 ) then
          exit
        end if
        !
        !  Otherwise, advance to time T1, and have RK4 estimate
        !  the solution U1 there.
        t1 = t0 + dt
        call rk4vec ( t0, n, u0, dt, derivsn, u1 )
        !  Shift the data to prepare for another step.

        t0 = t1
        u0(1:n) = u1(1:n)

        if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
            xx(k_step)=t1    !new time step
            yy(:,k_step)=u1  !new fields
            k_step=k_step+1
        end if
        i_step=i_step+1
    end do

    allocate(intensitys(3,numkeep))
    do z=1,numkeep,1  ! calculate the intensity for the last 'numkeep' values of y
        intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
        intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
        intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
    enddo
    print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', numkeep
    print*, 'shape intensity:', shape(intensitys(3,:))
    call local_max(intensitys(3,:),numkeep,5,debug) !searches for local maximum of the intensity, for the biffurcatoin diagrams.(prints to bynary file)

    !**************************************************************************************
    !write the variables to a binary file, to be read by python for plotting and analysis.

    OPEN(2,file='time.in',form='unformatted')
        WRITE(2) xx(:)
    CLOSE(2)

    OPEN(1,file='exr.in',form='unformatted')
        WRITE(1) yy(1,:)
    CLOSE(1)

    OPEN(1,file='pop.in',form='unformatted')
        WRITE(1) yy(9,:)
    CLOSE(1)

    OPEN(1,file='intensity.in',form='unformatted')
        WRITE(1) intensitys(3,:)
    CLOSE(1)

    OPEN(1,file='intensity_x2.in',form='unformatted')
        WRITE(1) intensitys(1,:)
    CLOSE(1)

    OPEN(1,file='intensity_y2.in',form='unformatted')
        WRITE(1) intensitys(2,:)
    CLOSE(1)

    deallocate(xx,yy,intensitys)
    return

end subroutine
    subroutine maxrk4f90_swipe()


    !18/3/2016

    !Routine used to swipe different values of w or m  in order to map the intensity dinamics
    !script for a runge kutta 4 (fixed step) integration of the Maxwell bloch equations with phase modulation.

    !Integrates the MX-BL phase mod equations with RK4 and output the fields in a matrix yy, only for the times
    !explicited in numkeep(the number of steps i will keep before the last.)

    !Input <-- dt(time step size), t0, u0(initial field values), tmax
    !Input <-- The parameters used for the normalizations are in Constants module.
    !Input <-- Numkeep(the number of steps i keep for display)

    !Output --> xx(numkeep) the times used in the integration. yy(9,numkeep) the field values.
    !Output --> Intensitys(3,numkeep) The intensity of the Ex field, the Ey field, and the total intensity.

    USE constants
    use my_lib
    use funcs

    IMPLICIT NONE

    !Rk4
    integer ( kind = 4 ), parameter :: n = 9 !(number of fields)
   ! real ( kind = 8 ), parameter :: dt = 1d0 !time step
    real ( kind = 8 ) t0 !initial time.
    real ( kind = 8 ) t1
    real ( kind = 8 ) tmax ! End time
    real ( kind = 8 ) u1(n)
    real ( kind = 8 ) u0(n) !Initial conditions

    !fields to keep
   ! real*8, parameter :: intime = 500.*20*gperp/10**6

    real*8, allocatable, dimension(:,:) :: intensitys, yy   !Fields
    real*8, allocatable, dimension(:) :: xx                 !Times
    integer :: numstep              !number of integration steps
    integer :: indexkeep            !index of the first step i keep
    integer :: numkeep              !number of steps i keep
    integer :: z
    integer, parameter :: debug=1    ! 1 to display debug info, 0 to display a cleaner output.
    integer i_step, k_step           !integration step and keep step.

    !SWIPE
    integer(4) :: count_peak
    real*8 :: w_stop, h
    integer :: i
    real*8, allocatable, dimension(:) :: data1,data2

    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RK4_Phasemod Maxwell'

    count_peak=0
    tmax = intime
    t0 = 0.0D+00
    u0=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
    numstep=int(tmax/dt)
    numkeep=numstep*3/20              !number of steps to keep on file(transitory)
    indexkeep=numstep-numkeep+1       !set indexkeep

    w_stop=0.0047
    wf=0.0035
    h=(w_stop-w)/20d0

    m=0.02d0
    do i=1,4
        allocate(xx(numkeep),yy(n,numkeep))


        print*, 't0=', t0, '', 'tmax=', tmax
        print*, 'w=', w
        print*, i

        i_step=1
        k_step=1

        do
            !  Stop if we've exceeded TMAX.
            if ( tmax <= t0 ) then
              exit
            end if
            !
            !  Otherwise, advance to time T1, and have RK4 estimate
            !  the solution U1 there.
            t1 = t0 + dt
            call rk4vec ( t0, n, u0, dt, derivsn, u1 )
            !  Shift the data to prepare for another step.

            t0 = t1
            u0(1:n) = u1(1:n)

            if ( i_step >= indexkeep ) then !if step > indexkeep, keep the values.
                xx(k_step)=t1    !new time step
                yy(:,k_step)=u1  !new fields
                k_step=k_step+1
            end if
            i_step=i_step+1
        end do

        print*, 'RK4 exit ok'

        allocate(intensitys(3,numkeep))
        do z=1,numkeep,1  ! calculate the intensity for the last 'numkeep' values of y
            intensitys(1,z)=yy(1,z)**2+yy(2,z)**2 !|E_x|^2
            intensitys(2,z)=yy(3,z)**2+yy(4,z)**2 !|E_y|^2
            intensitys(3,z)=sqrt(intensitys(1,z)+intensitys(2,z)) !|E|
        enddo

        !Until here is the integration

        print*, 'intensitys exit ok'
        print*, 'shape exr:', shape(yy(1,:)) ,' ', '=', numkeep
        print*, 'shape intensity:', shape(intensitys(3,:))
        print*, 'indexkeep=', indexkeep

        !**************************************************************************************
        !write the variables to a binary file, to be read by python for plotting and analysis.
        call local_max_bif_v2(intensitys(3,:),numkeep,5,count_peak,debug)

        !i dont need the index for the peaks.
        deallocate(xx,yy,intensitys)
        wf=wf+h
        tmax=t0+intime
    end do
!
    allocate(data1(count_peak))
    OPEN(2, file='peak.in', ACCESS='STREAM', FORM='UNFORMATTED')
    READ(2) data1
    CLOSE(2)
    print*, data1
    OPEN(3, file='peak_fort.in', FORM='UNFORMATTED' )
    write(3) data1
    CLOSE(3)

!    allocate(data2(count_peak))
!    OPEN(4, file='varparam_m0.in', ACCESS='STREAM', FORM='UNFORMATTED')
!    READ(4) data2
!    CLOSE(4)
!    OPEN(5, file='varparam_m0_fort.in', FORM='UNFORMATTED' )
!    write(5) data2
!    CLOSE(5)

    OPEN(3, file='varw0.in', ACCESS='STREAM', FORM='UNFORMATTED')
    READ(3) data1
    CLOSE(3)
    OPEN(3, file='varparam_w0_fort.in', FORM='UNFORMATTED' )
    write(3) data1
    CLOSE(3)

    deallocate(data1,data2)
   return

end subroutine

end module


!SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
!    INTEGER n,NMAX
!    REAL*8 h,x,dydx(n),y(n),yout(n)
!    EXTERNAL derivs
!    PARAMETER (NMAX=50)! Set to the maximum number of functions.
!    !Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
!    !the fourth-order Runge-Kutta method to advance the solution over an interval h and return
!    !the incremented variables as yout(1:n), which need not be a distinct array from y. The
!    !user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
!    INTEGER i
!    REAL*8 h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
!    hh=h*0.5
!    h6=h/6.
!    xh=x+hh
!    !16.1 Runge-Kutta Method 707
!    !visit website http://www.nr.com or call 1-800-872-7423 (North America only),
!    !or send email to trade@cup.cam.ac.uk (outside North America).
!    !readable files (including this one) to any server
!    !computer, is strictly prohibited. To order Numerical Recipes books,
!    !diskettes, or CDROMs
!    !Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machineCopyright
!    !(C) 1986-1992 by Cambridge University Press.
!    !Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
!    !Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
!    do i=1,n !First step.
!        yt(i)=y(i)+hh*dydx(i)
!    enddo
!    call derivs(xh,yt,dyt) !Second step.
!    do i=1,n
!        yt(i)=y(i)+hh*dyt(i)
!    enddo
!    call derivs(xh,yt,dym) !Third step.
!    do  i=1,n
!        yt(i)=y(i)+h*dym(i)
!        dym(i)=dyt(i)+dym(i)
!    enddo
!    call derivs(x+h,yt,dyt) !Fourth step.
!    do  i=1,n !Accumulate increments with proper weights.
!        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
!    enddo
!    return
!END
!
!
!SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs)
!    INTEGER nstep,nvar,NMAX,NSTPMX
!    PARAMETER (NMAX=9,NSTPMX=200000)! Maximum number of functions and
!    !maximum number of values to
!    !be stored.
!    REAL*8 x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX)
!    EXTERNAL derivs
!    COMMON /path/ xx,y !Storage of results.
!    !C USES rk4
!    !Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
!    !advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
!    !evaluates derivatives. Results are stored in the common block path. Be sure to dimension
!    !the common block appropriately.
!    INTEGER i,k
!    REAL*8 h,x,dv(NMAX),v(NMAX)
!    do  i=1,nvar !Load starting values.
!        v(i)=vstart(i)
!        y(i,1)=v(i)
!    enddo
!    xx(1)=x1
!    x=x1
!    h=(x2-x1)/nstep
!    do  k=1,nstep !Take nstep steps.
!        call derivs(x,v,dv)
!        call rk4(v,dv,nvar,x,h,v,derivs)
!        if(x+h.eq.x) then! pause! ’stepsize not significant in rkdumb’
!        end if
!        x=x+h
!        xx(k+1)=x !Store intermediate steps.
!        do  i=1,nvar
!            y(i,k+1)=v(i)
!    enddo
!    enddo
!    !print*,xx
!    return
!        !708 Chapter 16. Integration of Ordinary Differential Equations
!        !visit website http://www.nr.com or call 1-800-872-7423 (North America only),
!        !or send email to trade@cup.cam.ac.uk (outside North America).
!        !readable files (including this one) to any server
!        !computer, is strictly prohibited. To order Numerical Recipes books,
!        !diskettes, or CDROMs
!        !Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machineCopyright
!        !(C) 1986-1992 by Cambridge University Press.
!        !Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
!        !Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
!
!endsubroutine


!
!SUBROUTINE rkdumbmod(vstart,nvar,x1,x2,nstep,derivs,nmax,xx,y)
!    implicit none
!    INTEGER, intent(in) ::  nstep,nvar,NMAX
!    ! Maximum number of functions and
!    !maximum number of values to
!    !be stored.
!    REAL*8 :: x1,x2
!    real*8, dimension(nvar) :: vstart
!    real*8, allocatable, dimension(:) , intent(out) :: xx
!    real*8, allocatable, dimension(:,:), intent(out) :: y
!    EXTERNAL derivs
!    !COMMON /path/ xx,y !Storage of results.
!    !C USES rk4
!    !Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
!    !advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
!    !evaluates derivatives. Results are stored in the common block path. Be sure to dimension
!    !the common block appropriately.
!    INTEGER i,k
!    REAL*8 h,x,dv(NMAX),v(NMAX)
!    allocate(y(nmax,nstep),xx(nstep))
!    do  i=1,nvar !Load starting values.
!        v(i)=vstart(i)
!        y(i,1)=v(i)
!    enddo
!    xx(1)=x1
!    x=x1
!    h=(x2-x1)/nstep
!    do  k=1,nstep !Take nstep steps.
!        call derivs(x,v,dv)
!        call rk4(v,dv,nvar,x,h,v,derivs)
!        if(x+h.eq.x) then! pause! ’stepsize not significant in rkdumb’
!        end if
!        x=x+h
!        xx(k+1)=x !Store intermediate steps.
!        do  i=1,nvar
!            y(i,k+1)=v(i)
!    enddo
!    enddo
!    return
!endsubroutine


!***************************************************************************************
!subroutine rk4mio(yinit, time, numstep,f, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
!    !Newton method integration.
!
!    use constants
!
!    real(8), dimension(numstep), intent(inout) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
!    real(8), dimension(numstep) :: y
!    real(8), dimension(9), intent(in) :: yinit
!    real(8) :: dt
!    integer, intent(in) :: numstep
!    real(8), dimension(numstep), intent(in) :: time
!    integer :: i
!
!    integer :: ord=4
!    external :: f
!    real*8 :: ftemp
!
!    integer, optional :: debug
!
!    dt=(time(numstep)-time(1))/(numstep-1)
!
!    y(1)=yinit(1)!initial condition
!
!
!    !************** integrates ******************
!    do c=1,numstep-1!number of time steps
!        y(c+1)=y(c)
!        do i=1,ord
!            y(c+1)=y(c+1)+(1/6d0)*h*ftemp
!            if (c.eq.1) then
!                ftemp=f(time(c),y(c))!k1
!                elseif (c.eq.2) then
!                ftemp=2d0*f(time(c)+0.5d0*h, y(c)+0.5d0*h*ftemp)!k2
!                elseif (c.eq.3) then
!                    ftemp=2d0*f(time(c)+0.5d0*h, y(c)+0.5d0*h*ftemp)!k3
!                    elseif (c.eq.4) then
!                        ftemp=f(time(c)+h, y(c)+h*ftemp)!k4
!            end if
!        end do
!    end do
!
!end subroutine

