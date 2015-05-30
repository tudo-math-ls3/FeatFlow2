!##############################################################################
!# ****************************************************************************
!# <name> derivatives </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some definition to define derivative types.
!# Basically there are three different types of derivative quantifiers:
!# for 1D, 2D and 3D discretisations. The definitions are used e.g. to
!# specify which types of derivatives a Finite Element function should
!# calculate; see element.f90 for a deeper look how to work with these
!# constants!
!#
!# The following routines can be found in this module:
!#
!# 1.) der_igetID
!#     -> Interpret a string as derivative identifier. This is typically
!#        used for parsing INI files.
!# </purpose>
!##############################################################################

module derivatives

!$ use omp_lib
  use fsystem
  use genoutput

  implicit none

  private

  public :: der_igetID

!<constants>

!<constantblock description="Descriptors to identify derivative types in 1D">

  ! function value in term
  integer, parameter, public :: DER_FUNC1D       = 1

  ! 1st derivative in term
  integer, parameter, public :: DER_DERIV1D_X    = 2

  ! 2nd derivative in term
  integer, parameter, public :: DER_DERIV1D_XX   = 3

  ! 3rd derivative in term
  integer, parameter, public :: DER_DERIV1D_XXX  = 4

!</constantblock>

!<constantblock description="Descriptors to identify derivative types in 2D">

  ! function value in term
  integer, parameter, public :: DER_FUNC         = 1
  integer, parameter, public :: DER_FUNC2D       = 1

  ! x derivative in term
  integer, parameter, public :: DER_DERIV_X      = 2
  integer, parameter, public :: DER_DERIV2D_X    = 2

  ! y derivative in term
  integer, parameter, public :: DER_DERIV_Y      = 3
  integer, parameter, public :: DER_DERIV2D_Y    = 3

  ! 2nd x derivative in term
  integer, parameter, public :: DER_DERIV_XX     = 4
  integer, parameter, public :: DER_DERIV2D_XX   = 4

  ! xy derivative in term
  integer, parameter, public :: DER_DERIV_XY     = 5
  integer, parameter, public :: DER_DERIV2D_XY   = 5

  ! 2nd y derivative in term
  integer, parameter, public :: DER_DERIV_YY     = 6
  integer, parameter, public :: DER_DERIV2D_YY   = 6

  ! 3rd x derivative in term
  integer, parameter, public :: DER_DERIV_XXX    = 7
  integer, parameter, public :: DER_DERIV2D_XXX  = 7

  ! xxy derivative in term
  integer, parameter, public :: DER_DERIV_XXY    = 8
  integer, parameter, public :: DER_DERIV2D_XXY  = 8

  ! xyy derivative in term
  integer, parameter, public :: DER_DERIV_XYY    = 9
  integer, parameter, public :: DER_DERIV2D_XYY  = 9

  ! 3rd y derivative in term
  integer, parameter, public :: DER_DERIV_YYY    = 10
  integer, parameter, public :: DER_DERIV2D_YYY  = 10

!</constantblock>

!<constantblock description="Descriptors to identify derivative types in 3D">

  ! function value in term
  integer, parameter, public :: DER_FUNC3D      = 1

  ! x derivative in term
  integer, parameter, public :: DER_DERIV3D_X   = 2

  ! y derivative in term
  integer, parameter, public :: DER_DERIV3D_Y   = 3

  ! z derivative in term
  integer, parameter, public :: DER_DERIV3D_Z   = 4

  ! 2nd x derivative in term
  integer, parameter, public :: DER_DERIV3D_XX  = 5

  ! xy derivative in term
  integer, parameter, public :: DER_DERIV3D_XY  = 6

  ! 2nd y derivative in term
  integer, parameter, public :: DER_DERIV3D_YY  = 7

  ! xz derivative in term
  integer, parameter, public :: DER_DERIV3D_XZ  = 8

  ! yz derivative in term
  integer, parameter, public :: DER_DERIV3D_YZ  = 9

  ! zz derivative in term
  integer, parameter, public :: DER_DERIV3D_ZZ  = 10

  ! 3rd x derivative in term
  integer, parameter, public :: DER_DERIV3D_XXX = 11

  ! xxy derivative in term
  integer, parameter, public :: DER_DERIV3D_XXY = 12

  ! xxz derivative in term
  integer, parameter, public :: DER_DERIV3D_XXZ = 13

  ! xyy derivative in term
  integer, parameter, public :: DER_DERIV3D_XYY = 14
  
  ! xyz derivative in term
  integer, parameter, public :: DER_DERIV3D_XYZ = 15

  ! xzz derivative in term
  integer, parameter, public :: DER_DERIV3D_XZZ = 16

  ! 3rd y derivative in term
  integer, parameter, public :: DER_DERIV3D_YYY = 17

  ! yyz derivative in term
  integer, parameter, public :: DER_DERIV3D_YYZ = 18

  ! zzy derivative in term
  integer, parameter, public :: DER_DERIV3D_ZZY = 19

  ! 3rd z derivative in term
  integer, parameter, public :: DER_DERIV3D_ZZZ = 20
  
!</constantblock>

!<constantblock description="General constants.">

  ! Number of derivative identifiers.
  integer, parameter, public :: DER_MAXNDER  = 20

!</constantblock>

!</constants>

contains

  !****************************************************************************

!<function>

  integer function der_igetID(sderName, bcheck)

!<description>
  ! This routine returns the id of the derivative type to a given
  ! derivative type name. It is case-insensitive.
!</description>

!<result>
  ! id of the derivative type
!</result>

!<input>

  ! derivative type name - one of the DER_xxxx constants.
  character (LEN=*) :: sderName

  ! Check the derivative type. If set to TRUE and the derivative type
  ! name is invalid, the program is not stopped, but 0 is returned.
  logical, intent(in), optional :: bcheck

!</input>

!</function>

  character(len=len(sderName)+1) :: sder
  logical :: bchk

  ! SELECT CASE is not allowed for strings (although supported by a majority
  ! of compilers), therefore we have to use a couple of IF-commands :(
  ! select case(trim(sys_upcase(sderName)))

  sder = trim(sys_upcase(sderName))

  ! 1D cas 
  if (sder .eq. "DER_FUNC1D") then
    der_igetID = DER_FUNC1D
  elseif (sder .eq. "DER_DERIV1D_X") then
    der_igetID = DER_DERIV1D_X
  elseif (sder .eq. "DER_DERIV1D_XX") then
    der_igetID = DER_DERIV1D_XX
  elseif (sder .eq. "DER_DERIV1D_XXX") then
    der_igetID = DER_DERIV1D_XXX

    ! 2D case
  elseif (sder .eq. "DER_FUNC") then
    der_igetID = DER_FUNC
  elseif (sder .eq. "DER_FUNC2D") then
    der_igetID = DER_FUNC2D
  elseif (sder .eq. "DER_DERIV_X") then
    der_igetID = DER_DERIV_X
  elseif (sder .eq. "DER_DERIV2D_X") then
    der_igetID = DER_DERIV2D_X
  elseif (sder .eq. "DER_DERIV_Y") then
    der_igetID = DER_DERIV_Y
  elseif (sder .eq. "DER_DERIV2D_Y") then
    der_igetID = DER_DERIV2D_Y
  elseif (sder .eq. "DER_DERIV_XX") then
    der_igetID = DER_DERIV_XX
  elseif (sder .eq. "DER_DERIV2D_XX") then
    der_igetID = DER_DERIV2D_XX
  elseif (sder .eq. "DER_DERIV_XY") then
    der_igetID = DER_DERIV_XY
  elseif (sder .eq. "DER_DERIV2D_XY") then
    der_igetID = DER_DERIV2D_XY
  elseif (sder .eq. "DER_DERIV_YY") then
    der_igetID = DER_DERIV_YY
  elseif (sder .eq. "DER_DERIV2D_YY") then
    der_igetID = DER_DERIV2D_YY
  elseif (sder .eq. "DER_DERIV_XXX") then
    der_igetID = DER_DERIV_XXX
  elseif (sder .eq. "DER_DERIV2D_XXX") then
    der_igetID = DER_DERIV2D_XXX
  elseif (sder .eq. "DER_DERIV_XXY") then
    der_igetID = DER_DERIV_XXY
  elseif (sder .eq. "DER_DERIV2D_XXY") then
    der_igetID = DER_DERIV2D_XXY
  elseif (sder .eq. "DER_DERIV_XYY") then
    der_igetID = DER_DERIV_XYY
  elseif (sder .eq. "DER_DERIV2D_XYY") then
    der_igetID = DER_DERIV2D_XYY
  elseif (sder .eq. "DER_DERIV_YYY") then
    der_igetID = DER_DERIV_YYY
  elseif (sder .eq. "DER_DERIV2D_YYY") then
    der_igetID = DER_DERIV2D_YYY

    ! 3D case
  elseif (sder .eq. "DER_FUNC3D") then
    der_igetID = DER_FUNC3D
  elseif (sder .eq. "DER_DERIV3D_X") then
    der_igetID = DER_DERIV3D_X
  elseif (sder .eq. "DER_DERIV3D_Y") then
    der_igetID = DER_DERIV3D_Y
  elseif (sder .eq. "DER_DERIV3D_Z") then
    der_igetID = DER_DERIV3D_Z
  elseif (sder .eq. "DER_DERIV3D_XX") then
    der_igetID = DER_DERIV3D_XX
  elseif (sder .eq. "DER_DERIV3D_XY") then
    der_igetID = DER_DERIV3D_XY
  elseif (sder .eq. "DER_DERIV3D_XZ") then
    der_igetID = DER_DERIV3D_XZ
  elseif (sder .eq. "DER_DERIV3D_YY") then
    der_igetID = DER_DERIV3D_YY
  elseif (sder .eq. "DER_DERIV3D_YZ") then
    der_igetID = DER_DERIV3D_YZ
  elseif (sder .eq. "DER_DERIV3D_ZZ") then
    der_igetID = DER_DERIV3D_ZZ
  elseif (sder .eq. "DER_DERIV3D_XXX") then
    der_igetID = DER_DERIV3D_XXX
  elseif (sder .eq. "DER_DERIV3D_XXY") then
    der_igetID = DER_DERIV3D_XXY
  elseif (sder .eq. "DER_DERIV3D_XXZ") then
    der_igetID = DER_DERIV3D_XXZ
  elseif (sder .eq. "DER_DERIV3D_XYY") then
    der_igetID = DER_DERIV3D_XYY
  elseif (sder .eq. "DER_DERIV3D_XYY") then
    der_igetID = DER_DERIV3D_XYY
  elseif (sder .eq. "DER_DERIV3D_XYZ") then
    der_igetID = DER_DERIV3D_XYZ
  elseif (sder .eq. "DER_DERIV3D_XZZ") then
    der_igetID = DER_DERIV3D_XZZ
  elseif (sder .eq. "DER_DERIV3D_YYY") then
    der_igetID = DER_DERIV3D_YYY
  elseif (sder .eq. "DER_DERIV3D_YYZ") then
    der_igetID = DER_DERIV3D_YYZ
  elseif (sder .eq. "DER_DERIV3D_ZZZ") then
    der_igetID = DER_DERIV3D_ZZZ

  else
    bchk = .false.
    if (present(bcheck)) bchk = bcheck
    if (.not. bchk) then
      call output_line('Unknown derivative type: '//sderName,&
                      OU_CLASS_ERROR,OU_MODE_STD,'der_igetID')
      call sys_halt()
    else
      der_igetID = 0
    end if
  end if

  end function der_igetID
  
end module
