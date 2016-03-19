module funcs

    contains

subroutine local_max(valdata,size_data,compare_size,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    integer(4), intent(in) :: compare_size !even number, greater than 3 (3 or 5 give fast, good results on soft data.)
    integer, intent(in) :: debug
    integer, intent(in) :: size_data
    real(8), dimension(size_data), intent(in) :: valdata
    integer(4) :: j
    real(8), dimension(compare_size) :: localdata
    real(8) :: peak
    integer(4) :: count_peak, peak_index
    integer(4),dimension(1) :: maxind
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

!****************************************************************************************

subroutine local_max_bif_v1(valdata,size_data,compare_size,count_peak,debug)
      !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
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
    OPEN(1, file='peak_index.in', ACCESS='STREAM', FORM='UNFORMATTED', status='old')
    READ (1) data
    CLOSE(1)
    OPEN(1, file='peak_index_fort.in', FORM='UNFORMATTED' )
    write(1) data
    CLOSE(1)
    !-need  to add a way to test if the new peak is repeated or not. then if not, add no the file the new one.
    !-should i save the peaks? or only the index?.. if i have the index the i still have to gou though all the array.
    !if i save the peaks i dont have to spend memory on the arrays
    !- A cheaper way would be to keep only 5 steps at a time of the array, and look there for maxima.
    !So i don't have to use so much memory on the fields.

end subroutine

subroutine local_max_bif_v2(valdata,size_data,compare_size,count_peak,debug)
    !************************************************
    ! makes an dim(compare_size) array from data. if the maximum value is in the middle of the array (localdata(2)) saves the value to a binary file
    !************************************************
    use varparams !changing parameters for the bifurcation

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
    integer(4) :: midpoint
!    real*8, allocatable, dimension(:) :: data1,data2,data3
!    real*8 :: w0, m0
!    common /varparam/ w0,m0

    if (debug.eq.1) then
        print*,
        print*, '*************************************'
        print*, 'debug info for local_max_bif:  '
        print*,
        !print*, 'local_max_biff goes from : ', size_data-evaluate_lenght, 'to ', size_data-compare_size
    end if

    midpoint=(1+compare_size)/2
    !count_peak=0
!    OPEN(1,file='peak_index.in',form='unformatted',status='unknown',access='stream')
    OPEN(2,file='peak.in',form='unformatted',status='unknown',access='stream')
    OPEN(3,file='varw0.in',form='unformatted',status='unknown',access='stream')
    OPEN(4,file='varm0.in',form='unformatted',status='unknown',access='stream')

    do j=1,size_data-compare_size
        localdata=valdata(j:j+compare_size-1)
        maxind=maxloc(localdata)
        if (maxind(1).eq.(midpoint)) then
            peak=localdata(midpoint)
            count_peak=count_peak+1
            peak_index=j+midpoint-1
!                WRITE(1) peak_index
                WRITE(2) peak
                WRITE(3) w0
!                print*,'peak=', peak,'|','w=', w0
                WRITE(4) m0
        end if
    end do
!    CLOSE(1)
    CLOSE(2)
    CLOSE(3)
    CLOSE(4)

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
    !cant i close the files outside??

end subroutine

!****************************************************************************************

subroutine initial(intime,numstep, timeinit, time)
!**********************************************************************************************
!sets the amount of time steps and the time array to use, from the total time to integrate (intime), and the initial time (timeinit)

    use constants
    use my_lib

    implicit none

    real(8), intent(in) :: intime
    integer, intent(out) :: numstep
    real(8) :: jump
    real(8), intent(in) :: timeinit
    real(8), allocatable, dimension(:), intent(inout) :: time

    jump=.5d0
    numstep=int(intime/jump)
    allocate(time(numstep))
    call linspace(time,timeinit, timeinit+intime, numstep,1)
end subroutine

!****************************************************************************************

subroutine newtonint(yinit, time, numstep, exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop,debug)
    !Newton method integration.

    use constants

    real(8), dimension(numstep), intent(inout) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), dimension(9), intent(in) :: yinit
    real(8) :: dt
    integer, intent(in) :: numstep
    real(8), dimension(numstep), intent(in) :: time
    integer :: i

    integer, optional :: debug

    dt=(time(numstep)-time(1))/(numstep-1)

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
        exr(i+1)=exr(i)+dt*(-k*exr(i)+mu*rxr(i))
        exi(i+1)=exi(i)+dt*(-k*exi(i)+mu*rxi(i))
        eyr(i+1)=eyr(i)+dt*(-k*eyr(i)+mu*ryr(i))-eyi(i)*(Dphi0+m*cos(wf*time(i)))
        eyi(i+1)=eyi(i)+dt*(-k*eyi(i)+mu*ryi(i))+eyr(i)*(Dphi0+m*cos(wf*time(i)))
        rxr(i+1)=rxr(i)+dt*(-(rxr(i)-d*rxi(i))+exr(i)*pop(i))
        rxi(i+1)=rxi(i)+dt*(-(rxi(i)+d*rxr(i))+exi(i)*pop(i))
        ryr(i+1)=ryr(i)+dt*(-(ryr(i)-d*ryi(i))+eyr(i)*pop(i))
        ryi(i+1)=ryi(i)+dt*(-(ryi(i)+d*ryr(i))+eyi(i)*pop(i))
        pop(i+1)=pop(i)+dt*(-g*(pop(i)-D0+(exr(i)*rxr(i)+exi(i)*rxi(i)+eyr(i)*ryr(i)+eyi(i)*ryi(i))))
    enddo

    if (debug.eq.1) then
        !print*, 'variables size ', exr(1:20)
        print*, '******************************************************************'
    end if
end subroutine

!********************************************************************

subroutine neararray(yinit,numstep,time,dt)
!
    use constants
    implicit none
    real(8), allocatable, dimension(:) :: exr_near, exi_near,eyr_near,eyi_near,rxr_near,rxi_near,ryr_near,ryi_near,pop_near  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
    real(8), dimension(9) ,intent(in) :: yinit
    integer, intent(in) :: numstep
    real(8), dimension(numstep), intent(in) :: time
    real(8), allocatable, dimension(:) :: intensity_x2_near, intensity_y2_near, intensity_near
    integer :: i
    real(8),intent(in) :: dt
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


subroutine derivs_var(x,y,dydx)
    use constants
    implicit none
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

    dydx(1)=(-k*y(1)+mu*y(5))
    dydx(2)=(-k*y(2)+mu*y(6))
    dydx(3)=(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m*cos(w*x))
    dydx(4)=(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m*cos(w*x))
    dydx(5)=(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine


subroutine derivsn_var(x,n,y,dydx)
    use constants
    use varparams
    implicit none
    integer*4,intent(inout) :: n
    real*8, intent(in), dimension(9) :: y
    real*8, intent(out), dimension(9) :: dydx
    real*8, intent(in) :: x

!    real*8 :: w0, m0
!    common /varparam/ w0,m0

    dydx(1)=(-k*y(1)+mu*y(5))
    dydx(2)=(-k*y(2)+mu*y(6))
    dydx(3)=(-k*y(3)+mu*y(7))-y(4)*(Dphi0+m0*cos(w0*x))
    dydx(4)=(-k*y(4)+mu*y(8))+y(3)*(Dphi0+m0*cos(w0*x))
    dydx(5)=(-(y(5)-d*y(6))+y(1)*y(9))
    dydx(6)=(-(y(6)+d*y(5))+y(2)*y(9))
    dydx(7)=(-(y(7)-d*y(8))+y(3)*y(9))
    dydx(8)=(-(y(8)+d*y(7))+y(4)*y(9))
    dydx(9)=(-g*(y(9)-D0+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))))
end subroutine

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

end module




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
