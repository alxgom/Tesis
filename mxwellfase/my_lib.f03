module my_lib
    implicit none
    public :: linspace

    contains

    subroutine linspace(x,x_start, x_end, x_len,dir)
    !******************************************************************************
    !
    !dir=1 ---> from min to max.    dir=2 ---> from max to min.
    !******************************************************************************
        real(8), dimension(:), intent(inout) :: x
        real(8) :: x_start, x_end
        integer :: x_len, i
        real(8) :: dx
        integer :: dir

        dx=(x_end - x_start)/(x_len-1)
        print*, 'dx=', dx
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

