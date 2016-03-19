!****************************************************
!*   Program to demonstrate Laguerre coefficients   *
!* ------------------------------------------------ *
!* Reference: BASIC Scientific Subroutines, Vol. II *
!* by F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981.    *
!*                                                  *
!*            F90 Version by J.-P. Moreau, Paris.   *
!*                    (www.jpmoreau.fr)             *
!* ------------------------------------------------ *
!* SAMPLE RUN:                                      *
!*                                                  *
!* Laguerre polynomial coefficients for order 2     *
!*                                                  *
!*  A( 0) =         2                               *
!*  A( 1) =        -4                               *
!*  A( 2) =         1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 3     *
!*                                                  *
!*  A( 0) =         6                               *
!*  A( 1) =       -18                               *
!*  A( 2) =         9                               *
!*  A( 3) =        -1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 4     *
!*                                                  *
!*  A( 0) =        24                               *
!*  A( 1) =       -96                               *
!*  A( 2) =        72                               *
!*  A( 3) =       -16                               *
!*  A( 4) =         1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 5     *
!*                                                  *
!*  A( 0) =       120                               *
!*  A( 1) =      -600                               *
!*  A( 2) =       600                               *
!*  A( 3) =      -200                               *
!*  A( 4) =        25                               *
!*  A( 5) =        -1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 6     *
!*                                                  *
!*  A( 0) =       720                               *
!*  A( 1) =     -4320                               *
!*  A( 2) =      5400                               *
!*  A( 3) =     -2400                               *
!*  A( 4) =       450                               *
!*  A( 5) =       -36                               *
!*  A( 6) =         1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 7     *
!*                                                  *
!*  A( 0) =      5040                               *
!*  A( 1) =    -35280                               *
!*  A( 2) =     52920                               *
!*  A( 3) =    -29400                               *
!*  A( 4) =      7350                               *
!*  A( 5) =      -882                               *
!*  A( 6) =        49                               *
!*  A( 7) =        -1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 8     *
!*                                                  *
!*  A( 0) =     40320                               *
!*  A( 1) =   -322560                               *
!*  A( 2) =    564480                               *
!*  A( 3) =   -376320                               *
!*  A( 4) =    117600                               *
!*  A( 5) =    -18816                               *
!*  A( 6) =      1568                               *
!*  A( 7) =       -64                               *
!*  A( 8) =         1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 9     *
!*                                                  *
!*  A( 0) =    362880                               *
!*  A( 1) =  -3265920                               *
!*  A( 2) =   6531840                               *
!*  A( 3) =  -5080320                               *
!*  A( 4) =   1905120                               *
!*  A( 5) =   -381024                               *
!*  A( 6) =     42336                               *
!*  A( 7) =     -2592                               *
!*  A( 8) =        81                               *
!*  A( 9) =        -1                               *
!*                                                  *
!* Laguerre polynomial coefficients for order 10    *
!*                                                  *
!*  A( 0) =   3628800                               *
!*  A( 1) = -36288000                               *
!*  A( 2) =  81648000                               *
!*  A( 3) = -72576000                               *
!*  A( 4) =  31752000                               *
!*  A( 5) =  -7620480                               *
!*  A( 6) =   1058400                               *
!*  A( 7) =    -86400                               *
!*  A( 8) =      4050                               *
!*  A( 9) =      -100                               *
!*  A(10) =         1                               *
!*                                                  *
!****************************************************
! Laguerre Polynomials
! --------------------
!
! Laguerre polynomials are defined over the interval 0 <= x <= inf. The weighting
! function is the exponential w(x) = exp(-x):
!
!          inf.
!          Sum exp(-x) Ln(x) Lm(x) dx = 0   for n <> m     (3.9.6)
!           0
!		                      = 1   for n = m
!
! The recursion relation is
!
!          (n + 1) Ln+l(x) = (2n + 1 - x) Ln(x) - nLn(x)   (3.9.7)
!
!           where L0(x) = 1
!		  L1(x) = 1 - x
!
! Simply by changing a few lines in the LEGNDRE subroutine, a program for evaluating
! Laguerre polynomial coefficients can be generated [see program Laguerre).
!
! Note that the Laguerre polynomia1 coefficients are aIl integers, alternate in sign,
! and have a leading coefficient  an = (-1)^n. .
!------------------------------------------------------------------------------------
subroutine Laguerre()

real*8  A(0:10)
real*8  B(0:10,0:10)
integer n,k

   print *,' '
   do n = 2, 10
     write(*,25) n

     call Laguerre_Coeff(n,A,B)

     do k = 0, n
       write(*,50) k, A(k)
     end do
     print *,' '
     if (n<10) pause ' <Enter> to continue...'
   end do
   print *,' '

   stop
25 format(' Laguerre polynomial coefficients for order',i2/)
50 format('   A(',i2,') = ',f10.0)

end subroutine


!***************************************************
!* Laguerre polynomial coefficients evaluation by  *
!*'* means of recursion relation. The order of the *
!* polynomial is n. The coefficients are returned  *
!*'* in A(i).                                      *
!***************************************************
Subroutine Laguerre_Coeff(n,A,B)
  integer i,j,n
  real*8  A(0:10), B(0:10,0:10)
  !Establish l0 and l1 coefficients
  B(0,0)=1.d0 ; B(1,0)=1.d0 ; B(1,1)=-1.d0
  !Return if order is less than two
  if (n>1) then
    do i = 2, n
      B(i,0)=(2*i-1)*B(i-1,0)-(i-1)*(i-1)*B(i-2,0)
      do j = 1, i
        !Basic recursion relation
        B(i,j)=(2*i-1)*B(i-1,j)-B(i-1,j-1)-(i-1)*(i-1)*B(i-2,j)
      end do
    end do
    do i = 0, n
      A(i)=B(n,i)
    end do
  end if

  return
end


! End of file laguerre.f90
