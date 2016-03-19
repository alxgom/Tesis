program main
    implicit none
    real :: co
    integer :: j=0, k=0

    call lm_polynomial_coefficients( j, k, co )
    print *, co


end program

subroutine lm_polynomial_coefficients ( n, m, c )

!*****************************************************************************80
!
!! LM_POLYNOMIAL_COEFFICIENTS: coefficients of Laguerre polynomial Lm(n,m,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the
!    Laguerre polynomials of degree 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,0) = real ( m + 1, kind = 8 )
  c(1,1) = -1.0D+00

  do i = 2, n

    c(i,0:i) = ( real (   m + 2 * i - 1, kind = 8 ) * c(i-1,0:i)   &
               + real ( - m     - i + 1, kind = 8 ) * c(i-2,0:i) ) &
               / real (           i,     kind = 8 )

    c(i,1:i) = c(i,1:i) - c(i-1,0:i-1) / real ( i, kind = 8 )

  end do

  return
end

