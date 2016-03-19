 module binom
 implicit none

 integer, parameter :: iknd = selected_real_kind(31)
 real(iknd), allocatable ::  mm(:,:)

 contains

 recursive function combo(n,k) result(cmb)
 real (kind=iknd) :: cmb
 integer, intent(in) :: n,k
    if (k == n) then
    cmb = real(1,16)
 else if (k == 1) then
    cmb = real(n,16)
 else if (mm(n,k) /=0)  then
    cmb = mm(n,k)
 else if ((k /= 1) .and. (k /= n)) then
    cmb = combo(n-1,k-1) + combo(n-1,k)
    mm(n,k) = cmb
 end if
 end function

 end module
