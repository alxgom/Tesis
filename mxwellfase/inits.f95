subroutine initial(init,time,y)
!    use my_lib
!
!    implicit none
!    use globalpar
!    !real(8),intent(in) :: gperp !esto lo podria cambiar
!    real(8) :: intime
!    real(8), dimension(9), intent*(out) :: yinit
!    real(8), allocatable, intent(inout), dimension(1) :: y
!    real(8), allocatable, intent(inout), dimension(1) :: time
!    real(8), allocatable, intent(inout), dimension(1) :: init
!    real(8), allocatable, dimension(:), intent(out) :: timeinit
!    character :: new
!    character :: last
!    integer,intent(in) :: init
!
!    call globalparams
!
!    intime=500.*15*gperp/10**6 !#integration time FOR TRANDITORY
!    jump=1.
!    numstep=int(intime/jump)
!    print*, 'Type 1 for new initial conditions, or type 2 to continue with the last value'
!    !read*, init
!    if (init.eq.1) then
!        !'''User defined initial condition'''
!        linspace(timeinit,0d0, intime, numstep,1)
!        yinit=(/1., 1., 1.,-1., 1.,1., 1.,-1.9, 6.65973518*10**3/)
!    endif
!    if (init.eq.2) then
!        !'''initial condition from last simulation'''
!        linspace(timeinit,time(size(time)) ,intime*5/15+time(size(time)) , numstep, 1)
!        !yinit=y[-1]
!        yinit=y(size(y))
!    endif
!    return yinit, timeinit

end initial
