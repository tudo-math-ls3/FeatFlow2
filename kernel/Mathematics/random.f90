!##############################################################################
!# ****************************************************************************
!# <name> random </name>
!# ****************************************************************************
!# <purpose>
!# This module implements a pseudo-random number generator (rng).
!# The intention of this module is to provide an rng which generates a
!# (warning: paradoxon ahead) reproducible random number stream independent of
!# the platform and compiler in use.
!#
!# The algorithm used in this implementation is a 128-bit Xor-Shift agorithm,
!# as described in: George Marsaglia: Xorshift RNGs,
!# http://www.jstatsoft.org/v08/i14/paper
!#
!# The following routines can be found here:
!#
!# 1.) rng_init
!#     -> Initialises the rng data structure with a user-defined seed value.
!#
!# 2.) rng_initByClock
!#     -> Initialises the rng data structure using the current clock count
!#        as seed and returns the seed if desired.
!#
!# 3.) rng_get
!#     -> Returns the next random number in the stream.
!#
!# </purpose>
!##############################################################################

module random

use fsystem

implicit none

private

public :: t_random
public :: rng_init
public :: rng_initByClock
public :: rng_get
public :: rng_get_int32
public :: rng_get_int32_range
public :: rng_get_double
public :: rng_get_double_range
public :: rng_get_single
public :: rng_get_single_range

!<constants>

!<constantblock description="Default values for the random number generator">

  integer(I32), parameter, public :: RNG_DEFAULT_S = 428147976_I32
  integer(I32), parameter, public :: RNG_DEFAULT_X = 362436069_I32
  integer(I32), parameter, public :: RNG_DEFAULT_Y = 521288629_I32
  integer(I32), parameter, public :: RNG_DEFAULT_Z =  88675123_I32

!</constantblock>

!</constants>

!<types>

!<typeblock description="Random number generator data structure">

  type t_random

    integer(I32) :: s = RNG_DEFAULT_S
    integer(I32) :: x = RNG_DEFAULT_X
    integer(I32) :: y = RNG_DEFAULT_Y
    integer(I32) :: z = RNG_DEFAULT_Z

  end type

!</typeblock>

!</types>

interface rng_get
  module procedure rng_get_int32
  module procedure rng_get_int32_range
  module procedure rng_get_double
  module procedure rng_get_double_range
  module procedure rng_get_single
  module procedure rng_get_single_range
end interface

contains

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_init(rrng, iseed)

!<description>
  ! Initialises the random number generator data structure.
!</description>

!<output>
  ! The random number generator data structure to be initialised.
  type(t_random), intent(out) :: rrng
!</output>

!<input>
  ! OPTIONAL: The seed for the rng. If not given, the default seed is used.
  integer(I32), optional, intent(in) :: iseed
!</input>

!</subroutine>

    ! rrng is automatically initialised with default values, since it is
    ! intent(out), so we only need to set the seed - if it is given.
    if(present(iseed)) &
      rrng%s = iseed

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine rng_initByClock(rrng, iseed)

!<description>
  ! Initialises the random number generator data structure using the current
  ! processor clock count as a seed. The seed that was used can optionally be
  ! returned by this routine to assure reproducibility of the random number
  ! stream.
!</description>

!<output>
  ! The random number generator data structure to be initialised.
  type(t_random), intent(out) :: rrng
!</output>

!<output>
  ! OPTIONAL: Recieves the seed which was used for the random number generator.
  integer(I32), optional, intent(out) :: iseed
!</output>

!</subroutine>

  integer(I32) :: t

    ! get current processor clock
    call system_clock(t)

    ! return seed if desired
    if(present(iseed)) &
      iseed = t

    ! initialise the rng
    call rng_init(rrng, t)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_advance(rrng)

!<description>
  ! Advances the rng by one step. This routines does not need to be called by
  ! the user, as it is internally called by the rng_get_* routines.
!</description>

!<inputoutput>
  ! The rng data structure to be advanced.
  type(t_random), intent(inout) :: rrng
!</inputoutput>

!</subroutine>

  integer(I32) :: t

    t = ieor(rrng%s, ishft(rrng%s, 11))
    rrng%s = rrng%x
    rrng%x = rrng%y
    rrng%y = rrng%z
    rrng%z = ieor(ieor(rrng%z, ishft(rrng%z,-19)), ieor(t, ishft(t,-8)))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_get_int32(rrng, ix)

!<description>
  ! Returns the next non-negative 32-bit integer from the random number stream.
  ! The returned number is within the range {0,...,2^31-1}.
!</description>

!<inputoutput>
  ! The rng data structure.
  type(t_random), intent(inout) :: rrng
!</inputoutput>

!<output>
  ! Recieves the next random number within the range {0,...,2^31-1}.
  integer(I32), intent(out) :: ix
!</output>

!</subroutine>

  ! Sign bit mask: This mask has all bits set to 1 except for the sign bit.
  integer(I32), parameter :: cmask = huge(0_I32)

    ! advance the rng
    call rng_advance(rrng)

    ! mask out the sign bit to obtain a value in the range {0, ..., 2^31-1}
    ix = iand(rrng%z, cmask)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_get_int32_range(rrng, ix, imin, imax)

!<description>
  ! Returns the next 32-bit integer within a user-defined range from the
  ! random number stream.
  ! Notes:
  ! 1. For technical reasons, the range |imin-imax| must not exceed 2^31.
  ! 2. If imax < imin, then imin and imax are swapped.
  ! 3. If imin = imax, then the result ix is always set to imin.
!</description>

!<inputoutput>
  ! The rng data structure.
  type(t_random), intent(inout) :: rrng
!</inputoutput>

!<output>
  ! Recieves the next random number within the range {imin,...,imax}.
  integer(I32), intent(out) :: ix
!</output>

!<input>
  ! The range in which the random number ix should be contained. See notes in
  ! this routines description for more details.
  integer(I32), intent(in) :: imin, imax
!</input>

!</subroutine>

  integer(I32) :: t

    ! Note: It is intentional that we first get an int from the rng - and thus
    !       advance the rng - and check for imin = imax later.

    ! calculate a non-negative integer
    call rng_get_int32(rrng, t)

    ! transform into desired range
    if(imin .lt. imax) then
      ix = imin + mod(t, imax-imin+1)
    else if(imax .lt. imin) then
      ix = imax + mod(t, imin-imax+1)
    else
      ix = imin
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_get_double(rrng, dx)

!<description>
  ! Returns the next double precision value in interval [0,1] (boundaries
  ! included) from the random number stream.
!</description>

!<inputoutput>
  ! The rng data structure.
  type(t_random), intent(inout) :: rrng
!</inputoutput>

!<output>
  ! Recieves the next random number within the interval [0,1].
  real(DP), intent(out) :: dx
!</output>

!</subroutine>

  integer(I32) :: t

  ! calculate scaling factor: [0,2^31-1] -> [0,1]
  real(DP), parameter :: scl = 1.0_DP / real(huge(0_I32),DP)

    ! get a 32bit int
    call rng_get_int32(rrng, t)

    ! and scale it
    dx = scl * real(t,DP) 

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_get_double_range(rrng, dx, dmin, dmax)

!<description>
  ! Returns the next double precision value in a user-defined range from the
  ! random number stream.
!</description>

!<inputoutput>
  ! The rng data structure.
  type(t_random), intent(inout) :: rrng
!</inputoutput>

!<output>
  ! Recieves the next random number within the interval [dmin,dmax].
  real(DP), intent(out) :: dx
!</output>

!<input>
  ! The range in which the random number dx should be contained.
  real(DP), intent(in) :: dmin, dmax
!</input>

!</subroutine>

    call rng_get_double(rrng, dx)

    ! transform to range: [0,1] -> [dmin,dmax]
    dx = dmin + dx*(dmax-dmin)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_get_single(rrng, fx)

!<description>
  ! Returns the next single precision value in interval [0,1] (boundaries
  ! included) from the random number stream.
!</description>

!<inputoutput>
  ! The rng data structure.
  type(t_random), intent(inout) :: rrng
!</inputoutput>

!<output>
  ! Recieves the next random number within the interval [0,1].
  real(SP), intent(out) :: fx
!</output>

!</subroutine>

  real(DP) :: dx

    ! call double-precision version
    call rng_get_double(rrng, dx)

    ! and transform to single-precision
    fx = real(dx, SP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine rng_get_single_range(rrng, fx, fmin, fmax)

!<description>
  ! Returns the next single precision value in a user-defined range from the
  ! random number stream.
!</description>

!<inputoutput>
  ! The rng data structure.
  type(t_random), intent(inout) :: rrng
!</inputoutput>

!<output>
  ! Recieves the next random number within the interval [fmin,fmax].
  real(SP), intent(out) :: fx
!</output>

!<input>
  ! The range in which the random number dx should be contained.
  real(SP), intent(in) :: fmin, fmax
!</input>

!</subroutine>

  real(DP) :: dx

    ! call double-precision version
    call rng_get_double_range(rrng, dx, real(fmin,DP), real(fmax,DP))

    ! and transform to single-precision
    fx = real(dx, SP)

  end subroutine

end module