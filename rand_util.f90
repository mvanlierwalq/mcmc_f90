module rand_util

use shr_kind_mod, only: r8 => shr_kind_r8
use shr_RandNum_mod, only: ShrRandGen

implicit none
private

public :: random_real
public :: kissinit
public :: latin_random

class(ShrRandGen), allocatable :: rand_gen

contains

! Returns a random real number uniformly distributed over the range [0, 1).
real(r8) function random_real()
  real(r8) :: out_arr(1,1)

  call rand_gen%random(out_arr)
  random_real = out_arr(1,1)

end function random_real


subroutine kissinit(iinit)

  use shr_RandNum_mod, only: ShrKissRandGen

  integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
  integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

  integer(i4b) idum,ia,im,iq,ir,iinit
  integer(i4b) k,c1
  integer(i4b) :: seed(1,4)
  real(r8b)    rdum
  parameter (ia=16807,im=2147483647,iq=127773,ir=2836)

!!! Test integer representation !!!
  c1=-8
  c1=ishftc(c1,-3)
  if (c1.ne.536870911) then
     print *,'Nonstandard integer representation. Stoped.'
     stop 'ERROR IN RKISS'
  endif

  idum=iinit
  idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
  if (idum.eq.0) idum=1
  if (idum.ge.IM) idum=IM-1

  k=(idum)/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum = idum + IM
  if (idum.lt.1) then
     seed(1,1) = idum+1
  else 
     seed(1,1) = idum
  endif
  k=(idum)/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum = idum + IM
  if (idum.lt.1) then 
     seed(1,2) = idum+1
  else 
     seed(1,2) = idum
  endif
  k=(idum)/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum = idum + IM
  if (idum.lt.1) then
     seed(1,3) = idum+1
  else 
     seed(1,3) = idum
  endif
  k=(idum)/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum = idum + IM
  if (idum.lt.1) then
     seed(1,4) = idum+1
  else 
     seed(1,4) = idum
  endif

  if (allocated(rand_gen)) then
     ! Because RandNum doesn't use F2003 destructors, have to explicitly
     ! call this before replacing the value of rand_gen.
     call rand_gen%finalize()
     deallocate(rand_gen)
  end if

  allocate(rand_gen, source=ShrKissRandGen(seed))

  ! This is presumably in place to iterate past the LCG-generated seed?
  rdum = random_real()

end subroutine kissinit

subroutine perm_uniform(n, p)
!*****************************************************************************80
!
!! PERM_UNIFORM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n)

  do i = 1, n
    p(i) = i
  end do

  do i = 1, n - 1
    j = i4_uniform_ab ( i, n )
    k    = p(i)
    p(i) = p(j)
    p(j) = k
  end do

end subroutine perm_uniform

subroutine latin_random(dim_num, point_num, x)
!*****************************************************************************80
!
!! LATIN_RANDOM returns points in a Latin Random square.
!
!  Discussion:
!
!    In each spatial dimension, there will be exactly one
!    point whose coordinate value lies between consecutive
!    values in the list:
!
!      ( 0, 1, 2, ..., point_num ) / point_num
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) perm(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
!
!  Pick DIM_NUM * POINT_NUM random numbers between 0 and 1.
!
  !call r8mat_uniform_01 ( dim_num, point_num, seed, x )
  DO i=1,dim_num
    DO j=1,point_num
      x(i,j) = random_real()
    ENDDO
  ENDDO
!
!  For spatial dimension I,
!    pick a random permutation of 1 to POINT_NUM,
!    force the corresponding I-th components of X to lie in the
!    interval ( PERM(J)-1, PERM(J) ) / POINT_NUM.
!
  do i = 1, dim_num

    call perm_uniform ( point_num, perm )

    do j = 1, point_num
      x(i,j) = ( real ( perm(j) - 1, kind = 8 ) + x(i,j) ) &
               / real ( point_num, kind = 8 )
    end do

  end do

end subroutine latin_random

function i4_uniform_ab ( a, b )
!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform_ab
  real ( kind = 4 ) r
  real(8)           randy    !..double prec uniform(0,1)
  integer ( kind = 4 ) value

  !if ( seed == 0 ) then
  !  write ( *, '(a)' ) ' '
  !  write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
  !  write ( *, '(a)' ) '  Input value of SEED = 0.'
  !  stop 1
  !end if

  !k = seed / 127773

  !seed = 16807 * ( seed - k * 127773 ) - k * 2836

  !if ( seed < 0 ) then
  !  seed = seed + i4_huge
  !end if

  !r = real ( seed, kind = 4 ) * 4.656612875E-10
  randy = random_real()
  r  = real ( randy, kind = 4 )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00_4 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00_4 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00_4 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

end function i4_uniform_ab

end module rand_util
