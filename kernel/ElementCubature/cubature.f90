!##############################################################################
!# ****************************************************************************
!# <name> cubature </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a collection of cubature formulas. There are cubature
!# formulas for 1D, 2D and 3D included.
!#
!# A cubature formula is identified by its cubature ID CUB_xxxx.
!# By calling cub_getCubPoints with one of the ID`s, one gets the coordinates
!# and weights of that specific cubature formula on the reference element
!# (either <tex>$[-1,1]^2$</tex> or the triangle <tex>$(0,0), (0,1), (1,0)$</tex>
!# (depending on whether it is a cubature formulas for triangles or quadrilaterals.
!#
!# Furthermore, the routine cub_igetID allows to translate a string into
!# a cubature formula ID.
!#
!# The module supports standard cubature formulas and summed cubature formulas
!# which stem from a regular refinement of the reference element. The cubature
!# points of an underlying standard cubature formula are transferred to all
!# subelements of the refined reference element, thus realising a higher
!# order cubature also for nonsmooth functions. For a summed formula, the
!# ID contains additional information about the element refinement.
!#
!# The following routines can be found in this module:
!#
!#  1.) cub_igetID
!#      -> Interpret a string as cubature type identifier. This is typically
!#         used for parsing INI files.
!#
!#  2.) cub_igetNumPts
!#      -> Returns the number of points for a given cubature type identifier.
!#
!#  3.) cub_igetShape
!#      -> Returns the element shape identifier on which a cubature formula
!#         is defined.
!#
!#  4.) cub_igetCoordDim
!#      -> Returns the dimension of the coordinates for a given cubature
!#         type identifier.
!#
!#  5.) cub_getCubature
!#      -> Returns the cubature points and weights for a given cubature type
!#         identifier.
!#
!#  6.) cub_getCubPoints
!#      => DEPRECATED: See note below.
!#      -> Get information about cubature points for a specific cubature
!#         formula: Coordinates, number and weights.
!#
!#  7.) cub_getStdCubType
!#      -> For summed cubature formulas: Determine the underlying standard
!#         cubature formula
!#
!#  8.) cub_getRefLevels
!#      -> Determine the number of refinements of the reference element
!#         for summed cubature formulas.
!#
!#  9.) cub_getSummedCubType
!#      -> Get the ID of a summed cubature formula.
!#
!# 10.) cub_getName
!#      -> Returns the name of the cubature formula as a string.
!#
!# 11.) cub_resolveGenericCubType
!#      -> Resolves a generic cubature formula ID to the actual
!#         cubature formula ID.
!#
!# 12.) cub_isExtended
!#      -> Determine if a cubature formula ID refers to a generic cubature formula.
!#
!#
!# A note on cub_getCubPoints and cub_getCubature
!# ----------------------------------------------
!# The old routine 'cub_getCubPoints' and the parameter 'CUB_MAXCUBP' have
!# been marked as deprecated for several reasons:
!# 1. The CUB_MAXCUBP constant leads to a kind of 'sloppy' programming style,
!#    the dimension of array variables on the stack is set to this value which
!#    may lead to stack overflow once cubature rules with *many* points are
!#    added in future.
!#
!# 2. The cub_getCubPoints routine returned the arrays in Feat 1.x style, so
!#    that the array had to be transposed before it could be used by other
!#    Feat 2 kernel routines, which had to be done by the caller.
!#
!# 3. All cubature points and weights are 'hard-coded' in the cub_getCubPoints
!#    routines with a precision of maximal 16 digits, so the rules may be
!#    inaccurate when running Feat 2 on a platform with higher precision.
!#
!# Therefore a new routine 'cub_getCubature' has been implemented to replace
!# the old routine 'cub_getCubPoints'. To see the differences between the call
!# of the old 'cub_getCubPoints' and the new 'cub_getCubature' routine, see
!# the two code examples below.
!#
!# Remarks:
!# The new 'cub_getCubature' routine supports all cubature formulas that the
!# old 'cub_getCubPoints' supports, as it acts as a wrapper for the old routine
!# in the case that no new implementation for the given cubature formula exists.
!# Also, the new routine supports even more cubature formulas, e.g. Gauss-rules
!# up to degree 12 for 1D lines, 2D quadrilaterals and 3D hexahedra, which are
!# not supported by the old routine.
!#
!# Previously, one used a code like the one below to get the cubature formula:
!# <code>
!#
!#   ! Cubature points
!#   real(DP), dimension(CUB_MAXCUBP,4) :: Dxi
!#   real(DP), dimension(4,CUB_MAXCUBP) :: Dpoints
!#
!#   ! Cubature weights
!#   real(DP), dimension(CUB_MAXCUBP) :: Domega
!#
!#   ! Number of points
!#   integer :: ncubp
!#
!#     ! Get trapezoidal rule for quads
!#     call cub_getCubPoints(CUB_TRZ, ncubp, Dxi, Domega)
!#
!#     ! Transpose points
!#     Dpoints = transpose(Dxi)
!#
!# </code>
!#
!# The new recommended style is:
!# <code>
!#
!#   ! Cubature points
!#   real(DP), dimension(:,:), allocatable :: Dpoints
!#
!#   ! Cubature weights
!#   real(DP), dimension(:), allocatable :: Domega
!#
!#   ! Number of points and coordinate dimension
!#   integer :: ncubp, ncdim
!#
!#     ! Get number of points and coordinate dimension
!#     ncubp = cub_igetNumPts(CUB_TRZ)
!#     ncdim = cub_igetCoordDim(CUB_TRZ)
!#
!#     ! Allocate arrays
!#     allocate(Dpoints(ncdim,ncubp))
!#     allocate(Domega(ncubp))
!#
!#     ! Get trapezoidal rule for quads
!#     call cub_getCubature(CUB_TRZ, Dpoints, Domega)
!#
!#     ! Do not forget to deallocate Dpoints and Domega after you used them!
!#
!# </code>
!#
!# About the cubature ID
!# ---------------------
!# The cubature ID is interpreted as a bitfield that codes the necessary
!# information about the cubature rule. The following bits are used:
!#
!#  Bit 28-30: Type of the cubature formula.
!#             =0: Standard and summed cubature formula.
!#             =4: Identifies special constants like automatic cubature rules etc.
!#  Bit 16-27: For standard and summed cubature formulas:
!#             Refinement level of the reference element.
!#  Bit 0-15:  For standard and summed cubature formulas:
!#             Identifier for the underlying cubature formula.
!#
!# Standard and summed cubature formulas differ only in Bit 16-29 which
!# encodes a refinement level of the reference element; if this is set to a
!# value > 0, the cubature formula is summed, otherwise it is without summing.
!#
!# About generic cubature formulas
!# -------------------------------
!# There are a couple of generic cubature formulas available.
!# They are named "CUB_GEN_xxxx". These formulas apply for all types
!# of elements (1D, 2D, 3D, quad, triangle, hexa, ...) and define a
!# standard way for numerical integration. For the actual integration,
!# these formulas have to be "resolved" to the actual, dimension- and
!# element-shape dependent cubature formula. this can be done with the
!# routine "cub_resolveGenericCubType" which calculates the cubature
!# formula specifier according to the element shape.
!#
!# Note however, that the "CUB_GEN_AUTO" and "CUB_GEN_AUTO_LUMPMASS"
!# specifiers cannot be resolved that way. These specifiers are reserved
!# for being used with special finite elements and are therefore element
!# dependent. Due to the fact that this module does not know about
!# finite elements, appropriate finite element modules must implement
!# the behaviour if these constants are specified -- with a replacement
!# for "cub_resolveGenericCubType".
!#
!# </purpose>
!##############################################################################

module cubature

  use fsystem
  use basicgeometry
  use genoutput

  implicit none
  
  private
  
!<constants>

!<constantblock variable="ccubType" description="General values">

  ! Undefined cubature formula
  integer(I32), parameter, public :: CUB_UNDEFINED = 0
  
!</constantblock>

!<constantblock description="Special types of cubature formulas.">

  ! internal mask for special cubature formula types.
  integer(I32), parameter         :: CUB_TP_MASK = ishft(7_I32,28)

  ! Standard or summed cubature formula.
  integer(I32), parameter, public :: CUB_TP_STD = 0
  
  ! Automatic cubature formula.
  integer(I32), parameter, public :: CUB_TP_AUTO = ishft(1_I32,28)

  ! Deprecated cubature rule.
  integer(I32), parameter, public :: CUB_TP_DEPR = ishft(2_I32,28)

!</constantblock>

!<constantblock variable="ccubType" description="Generic cubature formulas">

  ! Automatic cubature formula
  integer(I32), parameter, public :: CUB_GEN_AUTO = CUB_TP_AUTO + 1

  ! Automatic cubature formula, lumping the mass matrix (if possible)
  integer(I32), parameter, public :: CUB_GEN_AUTO_LUMPMASS = CUB_TP_AUTO + 2

  ! Automatic cubature formula, 1-point Gauss formula
  integer(I32), parameter, public :: CUB_GEN_AUTO_G1 = CUB_TP_AUTO + 101

  ! Automatic cubature formula, 2-point Gauss formula
  integer(I32), parameter, public :: CUB_GEN_AUTO_G2 = CUB_TP_AUTO + 102

  ! Automatic cubature formula, 3-point Gauss formula
  integer(I32), parameter, public :: CUB_GEN_AUTO_G3 = CUB_TP_AUTO + 103

  ! Automatic cubature formula, 4-point Gauss formula
  integer(I32), parameter, public :: CUB_GEN_AUTO_G4 = CUB_TP_AUTO + 104

  ! Automatic cubature formula, 5-point Gauss formula
  integer(I32), parameter, public :: CUB_GEN_AUTO_G5 = CUB_TP_AUTO + 105

  ! Automatic cubature formula, 6-point Gauss formula
  integer(I32), parameter, public :: CUB_GEN_AUTO_G6 = CUB_TP_AUTO + 106

  ! Automatic cubature formula, Trapezoidal rule
  integer(I32), parameter, public :: CUB_GEN_AUTO_TRZ = CUB_TP_AUTO + 107

  ! Automatic cubature formula, Midpoint rule
  integer(I32), parameter, public :: CUB_GEN_AUTO_MID = CUB_TP_AUTO + 108

  ! Automatic cubature formula, Simpson rule
  integer(I32), parameter, public :: CUB_GEN_AUTO_SIMPSON = CUB_TP_AUTO + 109
  
!</constantblock>

!<constantblock variable="ccubType" description="1D formulas">

  ! 1-point Gauss formula, 1D, degree = 2, ncubp = 1
  integer(I32), parameter, public :: CUB_G1_1D = 101
 
  ! trapezoidal rule, 1D, degree = 2, ncubp = 2
  integer(I32), parameter, public :: CUB_TRZ_1D = 102

  ! 2-point Gauss formula, 1D, degree = 4, ncubp = 2
  integer(I32), parameter, public :: CUB_G2_1D = 103

  ! 3-point Gauss formula, 1D, degree = 6, ncubp = 3
  integer(I32), parameter, public :: CUB_G3_1D = 104

  ! 4-point Gauss formula, 1D, degree = 8, ncubp = 4
  integer(I32), parameter, public :: CUB_G4_1D = 105

  ! 5-point Gauss formula, 1D, degree = 10, ncubp = 5
  integer(I32), parameter, public :: CUB_G5_1D = 106

  ! Simpson-rule, 1D, degree = 4, ncubp = 3
  integer(I32), parameter, public :: CUB_SIMPSON_1D = 107

  ! 6-point Gauss formula, 1D, degree = 12, ncubp = 6
  integer(I32), parameter, public :: CUB_G6_1D = 108

  ! Pulcherima, 1D, degree = 4, ncubp = 4
  integer(I32), parameter, public :: CUB_PULCHERIMA_1D = 109

  ! Milne rule, 1D, degree = 5, ncubp = 5
  integer(I32), parameter, public :: CUB_MILNE_1D = 110

  ! 6-point rule, 1D, degree = 6, ncubp = 6
  integer(I32), parameter, public :: CUB_6POINT_1D = 111

  ! Weddle rule, 1D, degree = 7, ncubp = 7
  integer(I32), parameter, public :: CUB_WEDDLE_1D = 112

!</constantblock>

!<constantblock variable="ccubType" description="2D formulas, quad">

  ! 1x1 Gauss formula, degree = 2, ncubp = 1
  integer(I32), parameter, public :: CUB_G1X1 = 201
  integer(I32), parameter, public :: CUB_G1_2D = CUB_G1X1

  ! trapezoidal rule, degree = 2, ncubp = 4
  integer(I32), parameter, public :: CUB_TRZ = 202
  integer(I32), parameter, public :: CUB_TRZ_2D = 202

  ! midpoint rule, degree = 2, ncubp = 4
  integer(I32), parameter, public :: CUB_MID = 203
  integer(I32), parameter, public :: CUB_MID_2D = 203

  ! 2x2 Gauss formula, degree = 4, ncubp = 4
  integer(I32), parameter, public :: CUB_G2X2 = 204
  integer(I32), parameter, public :: CUB_G2_2D = CUB_G2X2

  ! Newton formula 1, degree = 4, ncubp = 4
  integer(I32), parameter, public :: CUB_NS1 = 205

  ! Newton formula 2, degree = 5, ncubp = 6
  integer(I32), parameter, public :: CUB_NS2 = 206

  ! Newton formula 3, degree = 6, ncubp = 7
  integer(I32), parameter, public :: CUB_NS3 = 207

  ! 3x3 Gauss formula, degree = 6, ncubp = 9
  integer(I32), parameter, public :: CUB_G3X3 = 208
  integer(I32), parameter, public :: CUB_G3_2D = CUB_G3X3

  ! Gauss formula, degree = 7, ncubp = 12
  integer(I32), parameter, public :: CUB_G = 209

  ! 4x4 Gauss formula, degree = 8, ncubp = 16
  integer(I32), parameter, public :: CUB_G4X4 = 210
  integer(I32), parameter, public :: CUB_G4_2D = CUB_G4X4

  ! 5x5 Gauss formula, degree = 10, ncubp = 25
  integer(I32), parameter, public :: CUB_G5X5 = 211
  integer(I32), parameter, public :: CUB_G5_2D = CUB_G5X5

  ! piecewise 1x1 Gauss formula, degree = 2, ncubp = 4
  integer(I32), parameter, public :: CUB_PG1X1 = 212

  ! piecewise trapezoidal rule, degree = 2, ncubp = 9
  integer(I32), parameter, public :: CUB_PTRZ = 213

  ! piecewise 2x2 Gauss formula, degree = 4, ncubp = 16
  integer(I32), parameter, public :: CUB_PG2X2 = 214

  ! piecewise 3x3 Gauss formula, degree = 6, ncubp = 36
  integer(I32), parameter, public :: CUB_PG3X3 = 215
  
  ! Simpson rule (corners and element midpoint), degree = 3, ncubp = 9
  integer(I32), parameter, public :: CUB_SIMPSON = 216
  integer(I32), parameter, public :: CUB_SIMPSON_2D = 216

  ! Simpson 3/8 rule (corners and 1/3 + 2/3), degree = 3, ncubp = 16
  integer(I32), parameter, public :: CUB_3_8 = 217

  ! 6x6 Gauss formula, degree = 12, ncubp = 36
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer(I32), parameter, public :: CUB_G6_2D = 218

!</constantblock>

!<constantblock variable="ccubType" description="2D formulas, tri">

  ! 1-point Gauss formula, triangle, degree = 2, ncubp = 1
  integer(I32), parameter, public :: CUB_G1_T = 250

  ! trapezoidal rule, triangle, degree = 2, ncubp = 3
  integer(I32), parameter, public :: CUB_TRZ_T = 251

  ! 3-point Gauss formula, triangle, degree = 3, ncubp = 3
  integer(I32), parameter, public :: CUB_G3_T = 252

  ! Collatz formula (Gauss formula, edge midpoints), degree = 3, ncubp = 3
  integer(I32), parameter, public :: CUB_Collatz = 253
  integer(I32), parameter, public :: CUB_Collatz_T = 253
  integer(I32), parameter, public :: CUB_G3MP_T = 253

  ! Vertices, midpoints, center, degree = 4, ncubp = 7
  integer(I32), parameter, public :: CUB_VMC = 254
!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, hexa">

  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer(I32), parameter, public :: CUB_G1_3D = 301

  ! midpoints of areas, 3D, degree = 2, ncubp = 6
  integer(I32), parameter, public :: CUB_MIDAREA_3D = 302

  ! trapezoidal rule, 3D, degree = 2, ncubp = 8
  integer(I32), parameter, public :: CUB_TRZ_3D = 303

  ! 2-point Gauss formula, 3D, degree = 4, ncubp = 8
  integer(I32), parameter, public :: CUB_G2_3D = 304
  
  ! 3-point Gauss formula, 3D, degree = 6, ncubp = 27
  integer(I32), parameter, public :: CUB_G3_3D = 305
  
  ! 4-point Gauss formula, 3D, degree = 8, ncubp = 64
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer(I32), parameter, public :: CUB_G4_3D = 306
  
  ! 5-point Gauss formula, 3D, degree = 10, ncubp = 125
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer(I32), parameter, public :: CUB_G5_3D = 307

  ! 6-point Gauss formula, degree = 12, ncubp = 216
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer(I32), parameter, public :: CUB_G6_3D = 308

!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, tetra">
  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer(I32), parameter, public :: CUB_G1_3D_T = 351
  
  ! trapezoidal rule, 3D, degree = 2, ncubp = 4
  integer(I32), parameter, public :: CUB_TRZ_3D_T = 352
  
  ! 4-point Stroud rule, 3D, degree = 2, ncubp = 4
  integer(I32), parameter, public :: CUB_S2_3D_T = 353
  
  ! 10-point Stroud rule, 3D, degree = 3, ncubp = 10
  integer(I32), parameter, public :: CUB_S3_3D_T = 354
  
  ! 15-point Stroud rule, 3D, degree = 5, ncubp = 15
  integer(I32), parameter, public :: CUB_S5_3D_T = 355
  
!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, pyramid">
  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer(I32), parameter, public :: CUB_G1_3D_Y = 401
  
  ! trapezoidal rule, 3D, degree = 2, ncubp = 5
  integer(I32), parameter, public :: CUB_TRZ_3D_Y = 402
  
!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, prism">
  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer(I32), parameter, public :: CUB_G1_3D_R = 451
  
  ! trapezoidal rule, 3D, degree = 2, ncubp = 6
  integer(I32), parameter, public :: CUB_TRZ_3D_R = 452
  
  ! 3x2-point Gauss formula, 3D, degree = 3 (maybe even 4?), ncubp = 6
  ! This formula is the 'cross-product' of the 2D CUB_G3_T formula
  ! for triangles and the 1D CUB_G2_1D formula.
  integer(I32), parameter, public :: CUB_G2_3D_R = 453
  
!</constantblock>

  ! DEPRECATED: maximal size of cubature node field
  integer, parameter, public :: CUB_MAXCUBP = 36

  ! DEPRECATED: maximal number of cubature points in 1D
  integer, parameter, public :: CUB_MAXCUBP_1D = 6
  
  ! DEPRECATED: maximal number of cubature points in 2D
  integer, parameter, public :: CUB_MAXCUBP_2D = 36
  
  ! DEPRECATED: maximal number of cubature points in 3D
  integer, parameter, public :: CUB_MAXCUBP_3D = 27

!</constants>

  public :: cub_igetID
  public :: cub_igetNumPts
  public :: cub_igetShape
  public :: cub_igetCoordDim
  public :: cub_getCubature
  public :: cub_getCubPoints
  public :: cub_getStdCubType
  public :: cub_getRefLevels
  public :: cub_getSummedCubType
  public :: cub_getName
  public :: cub_resolveGenericCubType
  public :: cub_isExtended
  
contains

  !****************************************************************************

!<function>

  integer(I32) function cub_igetID(scubName, bcheck)
  
!<description>
  ! This routine returns the cubature id to a given cubature formula name. It is
  ! case-insensitive.
!</description>

!<result>
  ! id of the cubature formula
!</result>

!<input>

  !cubature formula name - one of the CUB_xxxx constants.
  character (LEN=*) :: scubName

  ! Check the cubature type. If set to TRUE and the cubature
  ! name is invalid, the program is not stopped, but 0 is returned.
  logical, intent(in), optional :: bcheck

!</input>
  
!</function>

  character(len=len(scubName)+1) :: scub
  logical :: bchk
  
  ! SELECT CASE is not allowed for strings (although supported by a majority
  ! of compilers), therefore we have to use a couple of IF-commands :(
  ! select case(trim(sys_upcase(scubName)))

  scub = trim(sys_upcase(scubName))

  ! Generic formulas
  if (scub .eq. "AUTO") then
    cub_igetID=CUB_GEN_AUTO
  else if (scub .eq. "AUTO_LUMPMASS") then
    cub_igetID=  CUB_GEN_AUTO_LUMPMASS
  else if (scub .eq. "AUTO_G1") then
    cub_igetID=CUB_GEN_AUTO_G1
  else if (scub .eq. "AUTO_G2") then
    cub_igetID=CUB_GEN_AUTO_G2
  else if (scub .eq. "AUTO_G3") then
    cub_igetID=CUB_GEN_AUTO_G3
  else if (scub .eq. "AUTO_G4") then
    cub_igetID=CUB_GEN_AUTO_G4
  else if (scub .eq. "AUTO_G5") then
    cub_igetID=CUB_GEN_AUTO_G5
  else if (scub .eq. "AUTO_G6") then
    cub_igetID=CUB_GEN_AUTO_G6
  else if (scub .eq. "AUTO_TRZ") then
    cub_igetID=CUB_GEN_AUTO_TRZ
  else if (scub .eq. "AUTO_MID") then
    cub_igetID=CUB_GEN_AUTO_MID
  else if (scub .eq. "AUTO_SIMPSON") then
    cub_igetID=CUB_GEN_AUTO_SIMPSON

  ! 1D-formulas
  else if (scub .eq. "G1_1D") then
    cub_igetID=CUB_G1_1D
  else if (scub .eq. "TRZ_1D") then
    cub_igetID =CUB_TRZ_1D
  else if (scub .eq. "G2_1D") then
    cub_igetID=CUB_G2_1D
  else if (scub .eq. "G3_1D") then
    cub_igetID=CUB_G3_1D
  else if (scub .eq. "G4_1D") then
    cub_igetID=CUB_G4_1D
  else if (scub .eq. "G5_1D") then
    cub_igetID=CUB_G5_1D
  else if (scub .eq. "SIMPSON_1D") then
    cub_igetID=CUB_SIMPSON_1D
  else if (scub .eq. "G6_1D") then
    cub_igetID=CUB_G6_1D
  else if (scub .eq. "PULCHERIMA_1D") then
    cub_igetID=CUB_PULCHERIMA_1D
  else if (scub .eq. "MILNE_1D") then
    cub_igetID=CUB_MILNE_1D
  else if (scub .eq. "6POINT_1D") then
    cub_igetID=CUB_6POINT_1D
  else if (scub .eq. "WEDDLE_1D") then
    cub_igetID=CUB_WEDDLE_1D

  ! 2D-fomulas, quadrilateral
  else if (scub .eq. "G1X1" .or. scub .eq. "G1_2D") then
    cub_igetID=CUB_G1_2D
  else if (scub .eq. "TRZ" .or. scub .eq. "TRZ_2D") then
    cub_igetID=CUB_TRZ_2D
  else if (scub .eq. "MID" .or. scub .eq. "MID_2D") then
    cub_igetID=CUB_MID_2D
  else if (scub .eq. "SIMPSON" .or. scub .eq. "SIMPSON_2D") then
    cub_igetID=CUB_SIMPSON_2D
  else if (scub .eq. "G2X2" .or. scub .eq. "G2_2D") then
    cub_igetID=CUB_G2_2D
  else if (scub .eq. "NS1") then
    cub_igetID=CUB_NS1
  else if (scub .eq. "NS2") then
    cub_igetID=CUB_NS2
  else if (scub .eq. "NS3") then
    cub_igetID=CUB_NS3
  else if (scub .eq. "G3X3" .or. scub .eq. "G3_2D") then
    cub_igetID=CUB_G3_2D
  else if (scub .eq. "G") then
    cub_igetID=CUB_G
  else if (scub .eq. "G4X4" .or. scub .eq. "G4_2D") then
    cub_igetID=CUB_G4_2D
  else if (scub .eq. "G5X5" .or. scub .eq. "G5_2D") then
    cub_igetID=CUB_G5_2D
  else if (scub .eq. "PG1X1") then
    cub_igetID=CUB_PG1X1
  else if (scub .eq. "PTRZ") then
    cub_igetID=CUB_PTRZ
  else if (scub .eq. "PG2X2") then
    cub_igetID=CUB_PG2X2
  else if (scub .eq. "PG3X3") then
    cub_igetID=CUB_PG3X3
  else if (scub .eq. "G6_2D") then
    cub_igetID=CUB_G6_2D
    
  ! 2D-formulas, triangle
  else if (scub .eq. "G1_T") then
    cub_igetID=CUB_G1_T
  else if (scub .eq. "TRZ_T") then
    cub_igetID=CUB_TRZ_T
  else if (scub .eq. "G3_T") then
    cub_igetID=CUB_G3_T
  else if ((scub .eq. "G3MP_T") .or. (scub .eq. "COLLATZ")) then
    cub_igetID=CUB_G3MP_T
  else if (scub .eq. "VMC") then
    cub_igetID=CUB_VMC

  ! 3D-formulas, hexahedron
  else if (scub .eq. "G1_3D") then
    cub_igetID=CUB_G1_3D
  else if (scub .eq. "MIDAREA_3D") then
    cub_igetID =CUB_MIDAREA_3D
  else if (scub .eq. "TRZ_3D") then
    cub_igetID=CUB_TRZ_3D
  else if (scub .eq. "G2_3D") then
    cub_igetID=CUB_G2_3D
  else if (scub .eq. "G3_3D") then
    cub_igetID=CUB_G3_3D
  else if (scub .eq. "G4_3D") then
    cub_igetID=CUB_G4_3D
  else if (scub .eq. "G5_3D") then
    cub_igetID=CUB_G5_3D
  else if (scub .eq. "G6_3D") then
    cub_igetID=CUB_G6_3D
  
  ! 3D-formulas, tetrahedron
  else if (scub .eq. "G1_3D_T") then
    cub_igetID=CUB_G1_3D_T
  else if (scub .eq. "TRZ_3D_T") then
    cub_igetID=CUB_TRZ_3D_T
  else if (scub .eq. "S2_3D_T") then
    cub_igetID=CUB_S2_3D_T
  else if (scub .eq. "S3_3D_T") then
    cub_igetID=CUB_S3_3D_T
  else if (scub .eq. "S5_3D_T") then
    cub_igetID=CUB_S5_3D_T
  
  ! 3D-formulas, pyramid
  else if (scub .eq. "G1_3D_Y") then
    cub_igetID=CUB_G1_3D_Y
  else if (scub .eq. "TRZ_3D_Y") then
    cub_igetID=CUB_TRZ_3D_Y
  
  ! 3D-formulas, prism
  else if (scub .eq. "G1_3D_R") then
    cub_igetID=CUB_G1_3D_R
  else if (scub .eq. "TRZ_3D_R") then
    cub_igetID=CUB_TRZ_3D_R
  else if (scub .eq. "G2_3D_R") then
    cub_igetID=CUB_G2_3D_R

  else
    bchk = .false.
    if (present(bcheck)) bchk = bcheck
    if (.not. bchk) then
      call output_line('Unknown cubature formula: '//scubname,&
                      OU_CLASS_ERROR,OU_MODE_STD,'cub_igetID')
      call sys_halt()
    else
      cub_igetID = CUB_UNDEFINED
    end if
  end if
    
  end function cub_igetID

  ! ***************************************************************************

!<function>

  character(len=32) function cub_getName(ccubature) result(sname)

!<description>
  ! This function returns a string which represents the cubature formula name.
  !
  ! This function is constructed in a way such that for any valid cubature
  ! rule identifier ccubature, the equation
  ! <verb>
  !   cub_igetID(cub_getName(ccubature)) = ccubature
  ! </verb>
  ! holds.
  !
  ! However, please note that for a string scubature representing a valid
  ! cubature rule name, the equation
  ! <verb>
  !   cub_getName(cub_igetID(scubature)) = scubature
  ! </verb>
  ! does *not* hold in the general case as some cubature rules have multiple
  ! valid names, e.g. G2_2D and G2X2.
!</description>

!<input>
  ! The cubature formula whose name is to be returned.
  integer(I32), intent(in) :: ccubature
!</input>

!<result>
  ! A string representing the name of the element.
!</result>

!</function>

    select case(ccubature)
    ! General cubature formulas
    case (CUB_GEN_AUTO)
      sname ="AUTO"
    case (CUB_GEN_AUTO_LUMPMASS)
      sname ="AUTO_LUMPMASS"
    case (CUB_GEN_AUTO_G1)
      sname ="AUTO_G1"
    case (CUB_GEN_AUTO_G2)
      sname ="AUTO_G2"
    case (CUB_GEN_AUTO_G3)
      sname ="AUTO_G3"
    case (CUB_GEN_AUTO_G4)
      sname ="AUTO_G4"
    case (CUB_GEN_AUTO_G5)
      sname ="AUTO_G5"
    case (CUB_GEN_AUTO_G6)
      sname ="AUTO_G6"
    case (CUB_GEN_AUTO_TRZ)
      sname ="AUTO_TRZ"
    case (CUB_GEN_AUTO_MID)
      sname ="AUTO_MID"
    case (CUB_GEN_AUTO_SIMPSON)
      sname ="AUTO_SIMPSON"

    ! 1D formulas, line
    case (CUB_G1_1D)
      sname = 'G1_1D'
    case (CUB_G2_1D)
      sname = 'G2_1D'
    case (CUB_G3_1D)
      sname = 'G3_1D'
    case (CUB_G4_1D)
      sname = 'G4_1D'
    case (CUB_G5_1D)
      sname = 'G5_1D'
    case (CUB_G6_1D)
      sname = 'G6_1D'
    case (CUB_TRZ_1D)
      sname = 'TRZ_1D'
    case (CUB_SIMPSON_1D)
      sname = 'SIMPSON_1D'
    case (CUB_PULCHERIMA_1D)
      sname = 'PULCHERIMA_1D'
    case (CUB_MILNE_1D)
      sname = 'MILNE_1D'
    case (CUB_6POINT_1D)
      sname = '6POINT_1D'
    case (CUB_WEDDLE_1D)
      sname = 'WEDDLE_1D'

    ! 2D formulas, quadrilateral
    case (CUB_G1_2D)      ! alias: CUB_G1X1
      sname = 'G1_2D'
    case (CUB_G2_2D)      ! alias: CUB_G2X2
      sname = 'G2_2D'
    case (CUB_G3_2D)      ! alias: CUB_G3X3
      sname = 'G3_2D'
    case (CUB_G4_2D)      ! alias: CUB_G4X4
      sname = 'G4_2D'
    case (CUB_G5_2D)      ! alias: CUB_G5X5
      sname = 'G5_2D'
    case (CUB_G6_2D)
      sname = 'G6_2D'
    case (CUB_TRZ_2D)     ! alias: CUB_TRZ
      sname = 'TRZ_2D'
    case (CUB_MID_2D)     ! alias: CUB_MID
      sname = 'MID_2D'
    case (CUB_SIMPSON_2D) ! alias: CUB_SIMPSON
      sname = 'SIMPSON_2D'
    case (CUB_NS1)
      sname = 'NS1'
    case (CUB_NS2)
      sname = 'NS2'
    case (CUB_NS3)
      sname = 'NS3'
    case (CUB_G)
      sname = 'G'
    case (CUB_PG1X1)
      sname = 'PG1X1'
    case (CUB_PG2X2)
      sname = 'PG2X2'
    case (CUB_PG3X3)
      sname = 'PG3X3'
    case (CUB_PTRZ)
      sname = 'PTRZ'

    ! 2D formulas, triangle
    case (CUB_G1_T)
      sname = 'G1_T'
    case (CUB_G3_T)
      sname = 'G3_T'
    case (CUB_TRZ_T)
      sname = 'TRZ_T'
    case (CUB_G3MP_T)
      sname = 'CUB_G3MP_T'
    case (CUB_VMC)
      sname = 'VMC'

    ! 3D formulas, hexahedron
    case (CUB_G1_3D)
      sname = 'G1_3D'
    case (CUB_G2_3D)
      sname = 'G2_3D'
    case (CUB_G3_3D)
      sname = 'G3_3D'
    case (CUB_G4_3D)
      sname = 'G4_3D'
    case (CUB_G5_3D)
      sname = 'G5_3D'
    case (CUB_G6_3D)
      sname = 'G6_3D'
    case (CUB_MIDAREA_3D)
      sname = 'MIDAREA_3D'
    case (CUB_TRZ_3D)
      sname = 'TRZ_3D'

    ! 3D formulas, tetrahedron
    case (CUB_G1_3D_T)
      sname = 'G1_3D_T'
    case (CUB_TRZ_3D_T)
      sname = 'TRZ_3D_T'
    case (CUB_S2_3D_T)
      sname = 'S2_3D_T'
    case (CUB_S3_3D_T)
      sname = 'S3_3D_T'
    case (CUB_S5_3D_T)
      sname = 'S5_3D_T'

    ! 3D formulas, pyramid
    case (CUB_G1_3D_Y)
      sname = 'G1_3D_Y'
    case (CUB_TRZ_3D_Y)
      sname = 'TRZ_3D_Y'

    ! 3D formulas, prism
    case (CUB_G1_3D_R)
      sname = 'G1_3D_R'
    case (CUB_G2_3D_R)
      sname = 'G2_3D_R'
    case (CUB_TRZ_3D_R)
      sname = 'TRZ_3D_R'

    case default
      sname = '-unknown-'

    end select

  end function

  !****************************************************************************

!<function>

  elemental integer(I32) function cub_getStdCubType(ccubType) result(n)
  
!<description>
  ! Returns the underlying standard cubature formula if ccubType
  ! identifies a summed cubature formula.
!</description>

!<result>
  ! id of the underlying cubature formula.
!</result>

!<input>
  ! Cubature type identifier of ad adaptive cubature formula.
  integer(I32), intent(in) :: ccubType
!</input>
  
!</function>
  
    ! Blend out the highest 16 bit.
    if (iand(ccubType,CUB_TP_MASK) .eq. CUB_TP_STD) then
      n = iand(ccubType,int(2**16-1,I32))
    else
      n = ccubType
    end if
  
  end function

  !****************************************************************************

!<function>

  elemental integer(I32) function cub_getSummedCubType(ccubType,nlevels) result(n)
  
!<description>
  ! Creates the identifier of a summed cubature rule based on a
  ! standard cubature formula and a number of refinement levels of the
  ! reference element.
!</description>

!<result>
  ! id of the underlying cubature formula.
!</result>

!<input>
  ! Cubature type identifier of ad adaptive cubature formula.
  integer(I32), intent(in) :: ccubType
  
  ! Number of refinements of the reference element.
  integer, intent(in) :: nlevels
!</input>
  
!</function>
  
    ! Include the number of levels into ccubType.
    if (iand(ccubType,CUB_TP_MASK) .eq. CUB_TP_STD) then
      n = ior(cub_getStdCubType(ccubType),ishft(int(max(0,nlevels),I32),16))
    else
      n = 0
    end if
  
  end function

  !****************************************************************************

!<function>

  elemental integer function cub_getRefLevels(ccubType) result(n)
  
!<description>
  ! Returns the number of refinement levels of the reference element.
  ! This is used for summed cubature formulas to determine the refinement
  ! level.
!</description>

!<result>
  ! Number of refinement levels
!</result>

!<input>
  ! Cubature type identifier of ad adaptive cubature formula.
  integer(I32), intent(in) :: ccubType
!</input>
  
!</function>
    if ((iand(ccubType,CUB_TP_MASK) .eq. CUB_TP_STD) .or. &
        (iand(ccubType,CUB_TP_MASK) .eq. CUB_TP_AUTO)) then
      ! Get the level from the highest 16 bit.
      n = iand(ccubType,NOT(CUB_TP_MASK))
      n = ishft(n,-16)
    else
      n = 0
    end if
  
  end function

  !****************************************************************************

!<function>

  integer function cub_getRefElements(ccubType) result(n)
  
!<description>
  ! Auxiliary. Returns for a given summed cubature formula the number of
  ! elements into which the reference element is subdivided in every refinement
  ! step.
!</description>

!<result>
  ! Number of refinement levels
!</result>

!<input>
  ! Cubature type identifier of ad adaptive cubature formula.
  integer(I32), intent(in) :: ccubType
!</input>
  
!</function>

    integer(I32) :: ccubStd

    ! Get the underlying cubature formula
    ccubStd = cub_getStdCubType(ccubType)

    ! Get the shape of the element
    select case(cub_igetShape(ccubStd))
    case (BGEOM_SHAPE_LINE)
      ! One element is refined into two subelements
      n = 2

    case (BGEOM_SHAPE_QUAD,BGEOM_SHAPE_TRIA)
      ! One element is refined into four subelements
      n = 4
      
    case (BGEOM_SHAPE_HEXA)
      ! One element is refined into eight subelements
      n = 8
      
    case default
      call output_line ('Unsupported element.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'cub_getRefElements')
      call sys_halt()
    end select

!!$    integer :: nreflevels
!!$
!!$    ! Get the underlying cubature formula an number of refinement levels.
!!$    nreflevels = cub_getRefLevels(ccubType)
!!$
!!$    ! Based on this information, calculate the number of cubature points.
!!$    n = 0
!!$    select case (cub_getStdCubType(ccubType))
!!$      case (CUB_G1_1D,CUB_G2_1D,CUB_G3_1D,CUB_G4_1D,CUB_G5_1D,CUB_TRZ_1D)
!!$        ! One element is refined into two subelements
!!$        n = 2
!!$      case (CUB_G1_2D,CUB_G2_2D,CUB_G3_2D,CUB_G4_2D,CUB_G5_2D,CUB_TRZ_2D)
!!$        ! One element is refined into four subelements
!!$        n = 4
!!$      case default
!!$        call output_line ('Unsupported element.', &
!!$                          OU_CLASS_ERROR,OU_MODE_STD,'cub_getRefElements')
!!$        call sys_halt()
!!$    end select
  
  end function

  !****************************************************************************

!<function>

  integer function cub_igetNumPts(ccubType) result(n)
  
!<description>
  ! This routine returns the number of cubature points for a given cubature id.
!</description>

!<result>
  ! the number of points for the formula.
!</result>

!<input>
  ! Cubature type identifier
  integer(I32), intent(in) :: ccubType
!</input>
  
!</function>

    integer(I32) :: cstdCubType
    integer :: nreflevels,k
    
    ! Get the underlying cubature formula an dnumber of refinement levels.
    cstdCubType = cub_getStdCubType(ccubType)
    nreflevels = cub_getRefLevels(ccubType)

    select case(cstdCubType)
    ! -= 1D Line Formulas =-
    case (CUB_G1_1D)
      n = 1
    case (CUB_G2_1D,CUB_TRZ_1D)
      n = 2
    case (CUB_G3_1D,CUB_SIMPSON_1D)
      n = 3
    case (CUB_G4_1D,CUB_PULCHERIMA_1D)
      n = 4
    case (CUB_G5_1D,CUB_MILNE_1D)
      n = 5
    case (CUB_G6_1D,CUB_6POINT_1D)
      n = 6
    case (CUB_WEDDLE_1D)
      n = 7
    
    ! -= 2D Triangle Formulas =-
    case (CUB_G1_T)
      n = 1
    case (CUB_G3_T,CUB_TRZ_T,CUB_G3MP_T)
      n = 3
    case (CUB_VMC)
      n = 7
    
    ! -= 2D Quadrilateral Formulas =-
    case (CUB_G1_2D)
      n = 1
    case (CUB_G2_2D,CUB_TRZ,CUB_MID,CUB_NS1,CUB_PG1X1)
      n = 4
    case (CUB_NS2)
      n = 6
    case (CUB_NS3)
      n = 7
    case (CUB_G3_2D,CUB_PTRZ,CUB_SIMPSON)
      n = 9
    case (CUB_G)
      n = 12
    case (CUB_G4_2D,CUB_PG2X2,CUB_3_8)
      n = 16
    case (CUB_G5_2D)
      n = 25
    case (CUB_G6_2D,CUB_PG3X3)
      n = 36
      
    ! -= 3D Tetrahedron Formulas =-
    case (CUB_G1_3D_T)
      n = 1
    case (CUB_S2_3D_T,CUB_TRZ_3D_T)
      n = 4
    case (CUB_S3_3D_T)
      n = 10
    case (CUB_S5_3D_T)
      n = 15
      
    ! -= 3D Hexahedron Formulas =-
    case (CUB_G1_3D)
      n = 1
    case (CUB_MIDAREA_3D)
      n = 6
    case (CUB_G2_3D,CUB_TRZ_3D)
      n = 8
    case (CUB_G3_3D)
      n = 27
    case (CUB_G4_3D)
      n = 64
    case (CUB_G5_3D)
      n = 125
    case (CUB_G6_3D)
      n = 216
    
    ! -= 3D Pyramid Formulas =-
    case (CUB_G1_3D_Y)
      n = 1
    case (CUB_TRZ_3D_Y)
      n = 5
    
    ! -= 3D Prism Formulas =-
    case (CUB_G1_3D_R)
      n = 1
    case (CUB_TRZ_3D_R,CUB_G2_3D_R)
      n = 6
    
    case default
      call output_line('Unknown cubature type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'cub_igetNumPts')
      call sys_halt()
    end select
    
    ! Based on this information, if this is a summed cubature formula,
    ! calculate the number of cubature points.
    if (nreflevels .gt. 0) then
      select case (cstdCubType)
        case (CUB_G1_1D,CUB_G2_1D,CUB_G3_1D,CUB_G4_1D,CUB_G5_1D)
          ! All points are inner points. Total number =
          ! number per element * #elements in the reference element.
          ! This is an element in 1D, so every refinement brings 2 new
          ! elements.
          n = n * (cub_getRefElements(ccubType))**nreflevels

        case (CUB_G1_2D,CUB_G2_2D,CUB_G3_2D,CUB_G4_2D,CUB_G5_2D)
          ! All points are inner points. Total number =
          ! number per element * #elements in the reference element.
          ! This is a quad element in 2D, so every refinement brings 4 new
          ! elements.
          n = n * (cub_getRefElements(ccubType))**nreflevels
          
        case (CUB_G1_3D,CUB_G2_3D,CUB_G3_3D,CUB_G4_3D,CUB_G5_3D)
          ! All points are inner points. Total number =
          ! number per element * #elements in the reference element.
          ! This is a hexa element in 3D, so every refinement brings 8 new
          ! elements.
          n = n * (cub_getRefElements(ccubType))**nreflevels

        case (CUB_TRZ_1D)
          ! Points are in the corners of the elements.
          n = 2**nreflevels + 1
          
        case (CUB_TRZ_2D)
          ! Points are in the corners of the elements.
          n = 4**nreflevels + 2*2**nreflevels + 1

        case (CUB_TRZ_3D)
          ! Points are in the corners of the elements.
          n = (2**nreflevels + 1) * (4**nreflevels + 2*2**nreflevels + 1)

        case (CUB_G3_T)
          ! Points are in the midpoints of the elements. The computation of
          ! the total number of cubature points is more complicated for
          ! triangles. The number of midpoints follows from the Euler formula:
          !   "vertices - edges + (elements+1) = 2"
          ! The number of centers is simply the number of elements.
          ! So we first have to compute the number of vertices at the
          ! corners which can be computed from the following summation
          !    $sum_{i=0}^{k} ( sum_{j=0}^{k-i} 1 )$
          ! where k=2^nreflevels. This formula follows from considering
          ! all possible combinations of barycentric coordinates.
          ! An explicit expression of this formula is as follows:
          ! k*(k+1) - (1/2)*(k+1)^2 + (3/2)*k + 3/2
          k = 2**nreflevels
          n = k*(k+1) - ((k+1)**2)/2 + 3*k/2 + 3/2
          ! Now compute the number of edges from Euler formula
          n = n+4**nreflevels-1

        case (CUB_G1_T,CUB_G3MP_T)
          ! All points are inner points. Total number =
          ! number per element * #elements in the reference element.
          ! This is a triangle in 2D, so every refinement brings 4 new
          ! elements.
          n = n * (cub_getRefElements(ccubType))**nreflevels
          
        case (CUB_TRZ_T)
          ! Points are in the corners of the elements. The computation of
          ! the total number of cubature points is more complicated for
          ! triangles. It can be computed from the following summation
          !    $sum_{i=0}^{k} ( sum_{j=0}^{k-i} 1 )$
          ! where k=2^nreflevels. This formula follows from considering
          ! all possible combinations of barycentric coordinates.
          ! An explicit expression of this formula is as follows:
          ! k*(k+1) - (1/2)*(k+1)^2 + (3/2)*k + 3/2
          k = 2**nreflevels
          n = k*(k+1) - ((k+1)**2)/2 + 3*k/2 + 3/2

        case (CUB_VMC)
          ! Points are in the corners of the elements, in the midpoints
          ! and in the center. The computation of the total number of
          ! cubature points is more complicated for triangles.
          ! It can be computed from the following summation
          ! $sum_{i=0}^{k} ( sum_{j=0}^{k-i} 1 )$
          ! where k=2^nreflevels. This formula follows from considering
          ! all possible combinations of barycentric coordinates.
          ! An explicit expression of this formula is as follows:
          ! k*(k+1) - (1/2)*(k+1)^2 + (3/2)*k + 3/2
          k = 2**nreflevels
          n = k*(k+1) - ((k+1)**2)/2 + 3*k/2 + 3/2
          
          ! We still have to add the number of midpoints and centers.
          ! The number of midpoints follows from the Euler formula:
          !   "vertices - edges + (elements+1) = 2"
          ! The number of centers is simply the number of elements.
          n = 2*(n+4**nreflevels)-1

        case default
          
          call output_line ('Unsupported summed cubature formula.', &
                            OU_CLASS_ERROR,OU_MODE_STD,'cub_igetNumPts')
          call sys_halt()
      end select
    end if

  end function cub_igetNumPts


  !****************************************************************************

!<function>

  integer(I32) function cub_igetShape(ccubType) result(ishp)
  
!<description>
  ! This function returns the element shape identifier for a given cubature id.
  ! The shape identifier is one of the BGEOM_SHAPE_XXXX constants defined in
  ! basicgeometry.f90.
!</description>

!<result>
  ! One of the BGEOM_SHAPE_XXXX shape identifiers.
!</result>

!<input>
  ! Cubature type identifier
  integer(I32), intent(in) :: ccubType
!</input>
  
!</function>

    integer(I32) :: cstdCubType
    
    ! Get the underlying cubature formula an dnumber of refinement levels.
    cstdCubType = cub_getStdCubType(ccubType)

    if(cstdCubType .le. 100) then
      ! Invalid identifier
      ishp = BGEOM_SHAPE_UNKNOWN
    else if(cstdCubType .le. 200) then
      ! 1D Line
      ishp = BGEOM_SHAPE_LINE
    else if(cstdCubType .lt. 250) then
      ! 2D Quadrilateral
      ishp = BGEOM_SHAPE_QUAD
    else if(cstdCubType .le. 300) then
      ! 2D Triangle
      ishp = BGEOM_SHAPE_TRIA
    else if(cstdCubType .le. 350) then
      ! 3D Hexahedron
      ishp = BGEOM_SHAPE_HEXA
    else if(ccubType .le. 400) then
      ! 3D Tetrahedron
      ishp = BGEOM_SHAPE_TETRA
    else if(cstdCubType .le. 450) then
      ! 3D Pyramid
      ishp = BGEOM_SHAPE_PYRA
    else if(cstdCubType .le. 500) then
      ! 3D Prism
      ishp = BGEOM_SHAPE_PRISM
    else
      ! Invalid identifier
      ishp = BGEOM_SHAPE_UNKNOWN
    end if

  end function

  !****************************************************************************

!<function>

  integer(I32) function cub_resolveGenericCubType(cgeoShape,ccubType) result(ccubTypeAct)
  
!<description>
  ! Maps an automatic cubature formula to the actual cubature
  ! formula. If not possible, CUB_UNDEFINED is returned.
!</description>

!<result>
  ! The actual cubature formula or CUB_UNDEFINED, if the cubature
  ! formula could not be determined.
!</result>

!<input>
  ! Geometrical shape identifier which defines the shape of the
  ! underlying element primitive. One of the BGEOM_SHAPE_xxxx constants.
  integer(I32), intent(in) :: cgeoShape

  ! Cubature type identifier; one of the CUB_GEN_xxxx constants.
  integer(I32), intent(in) :: ccubType
!</input>

!<remark>
  ! The automatic cubature formulas CUB_GEN_AUTO and CUB_GEN_AUTO_LUMPMASS
  ! cannot be resolved with this routine! These depend on the discretisation
  ! and thus, must be determined with the corresponding cubature
  ! determination routine which involves the underlying finite element!
!</remark>
  
!</function>

    ccubTypeAct = CUB_UNDEFINED

    select case (cgeoShape)
    case (BGEOM_SHAPE_LINE)
      select case (ccubType)
      case (CUB_GEN_AUTO_G1)
        ccubTypeAct = CUB_G1_1D
        
      case (CUB_GEN_AUTO_G2)
        ccubTypeAct = CUB_G2_1D
        
      case (CUB_GEN_AUTO_G3)
        ccubTypeAct = CUB_G3_1D
        
      case (CUB_GEN_AUTO_G4)
        ccubTypeAct = CUB_G4_1D
        
      case (CUB_GEN_AUTO_G5)
        ccubTypeAct = CUB_G5_1D
        
      case (CUB_GEN_AUTO_G6)
        ccubTypeAct = CUB_G6_1D
        
      case (CUB_GEN_AUTO_TRZ)
        ccubTypeAct = CUB_TRZ_1D
        
      case (CUB_GEN_AUTO_MID)
        ccubTypeAct = CUB_G1_1D
        
      case (CUB_GEN_AUTO_SIMPSON)
        ccubTypeAct = CUB_SIMPSON_1D
        
      end select

    case (BGEOM_SHAPE_QUAD)
      select case (ccubType)
      case (CUB_GEN_AUTO_G1)
        ccubTypeAct = CUB_G1_2D
        
      case (CUB_GEN_AUTO_G2)
        ccubTypeAct = CUB_G2_2D
        
      case (CUB_GEN_AUTO_G3)
        ccubTypeAct = CUB_G3_2D
        
      case (CUB_GEN_AUTO_G4)
        ccubTypeAct = CUB_G4_2D
        
      case (CUB_GEN_AUTO_G5)
        ccubTypeAct = CUB_G5_2D
        
      case (CUB_GEN_AUTO_G6)
        ccubTypeAct = CUB_G6_2D
        
      case (CUB_GEN_AUTO_TRZ)
        ccubTypeAct = CUB_TRZ_2D
        
      case (CUB_GEN_AUTO_MID)
        ccubTypeAct = CUB_MID_2D
        
      case (CUB_GEN_AUTO_SIMPSON)
        ccubTypeAct = CUB_SIMPSON_2D
        
      end select

    case (BGEOM_SHAPE_TRIA)
      select case (ccubType)
      case (CUB_GEN_AUTO_G1)
        ccubTypeAct = CUB_G1_T

      case (CUB_GEN_AUTO_G2)
        ! Not implemented. Take G3.
        ccubTypeAct = CUB_G3_T
        
      case (CUB_GEN_AUTO_G3)
        ccubTypeAct = CUB_G3_T
        
      case (CUB_GEN_AUTO_TRZ)
        ccubTypeAct = CUB_TRZ_T
        
      case (CUB_GEN_AUTO_MID)
        ccubTypeAct = CUB_G3MP_T
        
      case (CUB_GEN_AUTO_SIMPSON)
        ccubTypeAct = CUB_VMC
        
      end select

    case (BGEOM_SHAPE_HEXA)
      select case (ccubType)
      case (CUB_GEN_AUTO_G1)
        ccubTypeAct = CUB_G1_3D
        
      case (CUB_GEN_AUTO_G2)
        ccubTypeAct = CUB_G2_3D
        
      case (CUB_GEN_AUTO_G3)
        ccubTypeAct = CUB_G3_3D
        
      case (CUB_GEN_AUTO_G4)
        ccubTypeAct = CUB_G4_3D
        
      case (CUB_GEN_AUTO_G5)
        ccubTypeAct = CUB_G5_3D
        
      case (CUB_GEN_AUTO_G6)
        ccubTypeAct = CUB_G6_3D
        
      case (CUB_GEN_AUTO_TRZ)
        ccubTypeAct = CUB_TRZ_3D
        
      case (CUB_GEN_AUTO_MID)
        ccubTypeAct = CUB_MIDAREA_3D
        
      end select

    case (BGEOM_SHAPE_TETRA)
      select case (ccubType)
      case (CUB_GEN_AUTO_G1)
        ccubTypeAct = CUB_G1_3D_T
        
      case (CUB_GEN_AUTO_G2)
        ccubTypeAct = CUB_S2_3D_T
        
      case (CUB_GEN_AUTO_G3)
        ccubTypeAct = CUB_S3_3D_T
        
      case (CUB_GEN_AUTO_G5)
        ccubTypeAct = CUB_S5_3D_T
        
      case (CUB_GEN_AUTO_TRZ)
        ccubTypeAct = CUB_TRZ_3D_T
        
      end select

    case (BGEOM_SHAPE_PYRA)
      select case (ccubType)
      case (CUB_GEN_AUTO_G1)
        ccubTypeAct = CUB_G1_3D_Y
        
      case (CUB_GEN_AUTO_G2)
        ccubTypeAct = CUB_TRZ_3D_Y
        
      end select

    case (BGEOM_SHAPE_PRISM)
      select case (ccubType)
      case (CUB_GEN_AUTO_G1)
        ccubTypeAct = CUB_G1_3D_R
        
      case (CUB_GEN_AUTO_TRZ)
        ccubTypeAct = CUB_TRZ_3D_R
        
      end select

    end select

  end function

  !****************************************************************************

!<function>

  integer function cub_isExtended(ccubType) result(igeneric)
  
!<description>
  ! Determins if ccubType refers to an extended cubature rule.
!</description>

!<result>
  ! =0, if ccubType directly specifies a cubature rule.
  ! =1, if ccubType refers to a generic (CUB_GEN_xxxx) cubature formula
  ! =2, if ccubType refers to a generic cubature formula which cannot be
  !     resolved due to implementation specific details, e.g. FEM information.
  !     I.e., ccubType is =CUB_GEN_AUTO or =CUB_GEN_AUTO_LUMPMASS.
  ! =3, if ccubType refers to a deprecated cubature rule.
!</result>

!<input>
  ! Cubature type identifier; one of the CUB_GEN_xxxx constants.
  integer(I32), intent(in) :: ccubType
!</input>

!</function>

    select case (iand(ccubType,CUB_TP_MASK))
    case (CUB_TP_AUTO)
      igeneric = 1
      if ((ccubType .eq. CUB_GEN_AUTO) .or.&
          (ccubType .eq. CUB_GEN_AUTO_LUMPMASS)) then
        igeneric = 2
      end if
    case (CUB_TP_DEPR)
      igeneric = 3
    case default
      igeneric = 0
    end select

  end function

  !****************************************************************************

!<function>

  integer function cub_igetCoordDim(ccubType) result(n)
  
!<description>
  ! This routine returns the dimension of the coordinates for a given
  ! cubature type identifier.
!</description>

!<result>
  ! The coordinate dimension of the formula.
!</result>

!<input>
  ! Cubature type identifier
  integer(I32), intent(in) :: ccubType
!</input>
  
!</function>

    ! Get the shape of the cubature id
    select case(cub_igetShape(cub_getStdCubType(ccubType)))
    case(BGEOM_SHAPE_LINE)
      ! 1D line formula -> reference coordinates
      n = 1
    case(BGEOM_SHAPE_TRIA)
      ! 2D triangle formula -> barycentric coordinates
      n = 3
    case(BGEOM_SHAPE_QUAD)
      ! 2D quadrilateral formula -> reference coordinates
      n = 2
    case(BGEOM_SHAPE_TETRA)
      ! 3D tetrahedron formula -> barycentric coordinates
      n = 4
    case(BGEOM_SHAPE_HEXA)
      ! 3D hexahedron formula -> reference coordinates
      n = 3
    case(BGEOM_SHAPE_PYRA)
      ! 3D pyramid formula -> reference coordinates
      n = 3
    case(BGEOM_SHAPE_PRISM)
      ! 3D prism formula -> reference coordinates
      n = 3
    case default
      ! unknown formula
      n = 0
    end select

  end function cub_igetCoordDim

  !****************************************************************************

!<subroutine>

  recursive subroutine cub_getCubature(ccubType, Dpoints, Domega)

!<description>
  ! This routine returns the cubature point coordinates and the corresponding
  ! weights for a given cubature type identifier.
!</description>

!<input>
  ! id of the cubature formula to be set
  integer(I32), intent(in) :: ccubType
!</input>

!<output>
  ! Coordinates of the cubature points.
  ! The dimension is assumed to be at least Dpoints(1:NDIM,1:NPTS), where:
  ! NDIM is the dimension of a point, which is returned by cub_igetCoordDim.
  ! NPTS is the number of cubatue points, which is returned by cub_igetNumPts.
  real(DP), dimension(:,:), intent(out) :: Dpoints
  
  ! For every cubature point the corresponding cubature weight.
  ! The dimension is assumed to be at least Domega(1:NPTS) (see Dpoints).
  real(DP), dimension(:), intent(out) :: Domega
!</output>

!</subroutine>

  ! local variables for wrapper
  real(DP), dimension(CUB_MAXCUBP,4) :: Dxi
  integer :: i,j,k,l,ncubp
  
  ! Auxiliary arrays for Gauss-Legendre rules
  real(DP), dimension(6) :: Dv, Dw
  integer :: ndim, npts

  ! Variables needed for summed cubature formulas.
  integer(I32) :: cstdCubType
  integer :: nreflevels,icubp,ncubpts,nrefcubpts,icubpStart
  integer :: isubelement, nsubelements, isubx, isuby, isubz
  real(DP) :: dweight,dedgelen
  real(DP), dimension(:,:), allocatable :: DpointsLocal,DpointsRef
  real(DP), dimension(:), allocatable :: DomegaLocal
  
    ! Get the underlying cubature formula an dnumber of refinement levels.
    cstdCubType = cub_getStdCubType(ccubType)
    nreflevels = cub_getRefLevels(ccubType)
  
    if (nreflevels .le. 0) then
    
      ! Standard cubature formula.

      ndim = 0
      npts = 0
      
      ! Okay, let us see what cubature rule we have here...
      select case(cstdCubType)
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 1D LINE CUBATURE RULES
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      case (CUB_G1_1D)
        ndim = 1
        npts = 1
      
      case (CUB_G2_1D)
        ndim = 1
        npts = 2
      
      case (CUB_G3_1D)
        ndim = 1
        npts = 3
      
      case (CUB_G4_1D)
        ndim = 1
        npts = 4
      
      case (CUB_G5_1D)
        ndim = 1
        npts = 5
      
      case (CUB_G6_1D)
        ndim = 1
        npts = 6

      case(CUB_TRZ_1D) ! trapezoidal rule
        Dpoints(1,1) = -1.0_DP
        Dpoints(1,2) =  1.0_DP
        Domega(1) = 1.0_DP
        Domega(2) = 1.0_DP

      case(CUB_SIMPSON_1D) ! Simpson rule
        Dpoints(1,1) = -1.0_DP
        Dpoints(1,2) =  0.0_DP
        Dpoints(1,3) =  1.0_DP
        Domega(1) = 1.0_DP / 3.0_DP
        Domega(2) = 4.0_DP / 3.0_DP
        Domega(3) = 1.0_DP / 3.0_DP

      case(CUB_PULCHERIMA_1D) ! Pulcherima (3/7-Simpson rule)
        Dpoints(1,1) = -1.0_DP
        Dpoints(1,2) = -1.0_DP / 3.0_DP
        Dpoints(1,3) =  1.0_DP / 3.0_DP
        Dpoints(1,4) =  1.0_DP
        Domega(1) = 1.0_DP / 4.0_DP
        Domega(2) = 3.0_DP / 4.0_DP
        Domega(3) = 3.0_DP / 4.0_DP
        Domega(4) = 1.0_DP / 4.0_DP

      case(CUB_MILNE_1D) ! Milne rule
        Dpoints(1,1) = -1.0_DP
        Dpoints(1,2) = -0.5_DP
        Dpoints(1,3) =  0.0_DP
        Dpoints(1,4) =  0.5_DP
        Dpoints(1,5) =  1.0_DP
        Domega(1) =  7.0_DP / 45.0_DP
        Domega(2) = 32.0_DP / 45.0_DP
        Domega(3) = 12.0_DP / 45.0_DP
        Domega(4) = 32.0_DP / 45.0_DP
        Domega(5) =  7.0_DP / 45.0_DP
        
      case(CUB_6POINT_1D) ! 6-point rule
        Dpoints(1,1) = -1.0_DP
        Dpoints(1,2) = -0.6_DP
        Dpoints(1,3) = -0.2_DP
        Dpoints(1,4) =  0.2_DP
        Dpoints(1,5) =  0.6_DP
        Dpoints(1,6) =  1.0_DP
        Domega(1) = 19.0_DP / 144.0_DP
        Domega(2) = 75.0_DP / 144.0_DP
        Domega(3) = 50.0_DP / 144.0_DP
        Domega(4) = 50.0_DP / 144.0_DP
        Domega(5) = 75.0_DP / 144.0_DP
        Domega(6) = 19.0_DP / 144.0_DP

      case(CUB_WEDDLE_1D) ! Weddle rule
        Dpoints(1,1) = -1.0_DP
        Dpoints(1,2) = -2.0_DP / 3.0_DP
        Dpoints(1,3) = -1.0_DP / 3.0_DP
        Dpoints(1,4) =  0.0_DP
        Dpoints(1,5) =  1.0_DP / 3.0_DP
        Dpoints(1,6) =  2.0_DP / 3.0_DP
        Dpoints(1,7) =  1.0_DP
        Domega(1) =  41.0_DP / 420.0_DP
        Domega(2) = 216.0_DP / 420.0_DP
        Domega(3) =  27.0_DP / 420.0_DP
        Domega(4) = 272.0_DP / 420.0_DP
        Domega(5) =  27.0_DP / 420.0_DP
        Domega(6) = 216.0_DP / 420.0_DP
        Domega(7) =  41.0_DP / 420.0_DP

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 2D TRIANGLE CUBATURE RULES
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 2D QUADRILATERAL CUBATURE RULES
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      case (CUB_G1_2D)
        ndim = 2
        npts = 1
      
      case (CUB_G2_2D)
        ndim = 2
        npts = 2
      
      case (CUB_G3_2D)
        ndim = 2
        npts = 3
      
      case (CUB_G4_2D)
        ndim = 2
        npts = 4
      
      case (CUB_G5_2D)
        ndim = 2
        npts = 5
      
      case (CUB_G6_2D)
        ndim = 2
        npts = 6

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 3D TETRAHEDRON CUBATURE RULES
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 3D HEXAHEDRON CUBATURE RULES
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      case (CUB_G1_3D)
        ndim = 3
        npts = 1
      
      case (CUB_G2_3D)
        ndim = 3
        npts = 2
      
      case (CUB_G3_3D)
        ndim = 3
        npts = 3
      
      case (CUB_G4_3D)
        ndim = 3
        npts = 4
      
      case (CUB_G5_3D)
        ndim = 3
        npts = 5
      
      case (CUB_G6_3D)
        ndim = 3
        npts = 6

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 3D PYRAMID CUBATURE RULES
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 3D PRISM CUBATURE RULES
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! OLD INTERFACE WRAPPER
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      case default
    
        ! call old routine
        call cub_getCubPoints(ccubType, ncubp, Dxi, Domega)
        
        ! transpose the point array
        do i = 1, min(4,size(Dpoints,1))
          do j = 1, min(ncubp,size(Dpoints,2))
            Dpoints(i,j) = Dxi(j,i)
          end do
        end do
        
        ! And get out of here
        return
        
      end select
      
      ! Gauss-Legendre rule ?
      if((ndim .gt. 0) .and. (npts .gt. 0)) then
      
        ! Okay, get the corresponding Gauss-Legendre rule
        call cub_auxGaussLegendre(npts,Dv,Dw)
        
        ! What dimension should the rule be for?
        select case(ndim)
        case (NDIM1D)
          ! 1D Gauss rule
          do i = 1, npts
            Dpoints(1,i) = Dv(i)
            Domega(i) = Dw(i)
          end do
        
        case (NDIM2D)
          ! 2D Gauss rule for quadrilaterals
          l = 1
          do i = 1, npts
            do j = 1, npts
              Dpoints(1,l) = Dv(i)
              Dpoints(2,l) = Dv(j)
              Domega(l) = Dw(i)*Dw(j)
              l = l+1
            end do
          end do
        
        case(NDIM3D)
          ! 3D Gauss rule for hexahedra
          l = 1
          do i = 1, npts
            do j = 1, npts
              do k = 1, npts
                Dpoints(1,l) = Dv(i)
                Dpoints(2,l) = Dv(j)
                Dpoints(3,l) = Dv(k)
                Domega(l) = Dw(i)*Dw(j)*Dw(k)
                l = l+1
              end do
            end do
          end do
          
        end select
      
      end if
  
    else
    
      ! This is a summed cubature formula, so determining the cubature points
      ! is slightly more complicated.
      !
      ! In a first step, calculate the points on the reference element
      ! of the standard cubature formula.
      ncubpts = cub_igetNumPts(cstdCubType)
      
      allocate(DpointsLocal(cub_igetCoordDim(cstdCubType),ncubpts))
      allocate(DomegaLocal(ncubpts))
      
      call cub_getCubature(cstdCubType, DpointsLocal, DomegaLocal)
      
      select case (cstdCubType)
        case (CUB_G1_1D,CUB_G2_1D,CUB_G3_1D,CUB_G4_1D,CUB_G5_1D)
          ! All points are innter points. Total number =
          ! number per element * #elements in the reference element.
          ! This is an element in 1D, so every refinement brings 2 new
          ! elements.
          nsubelements = cub_getRefElements(ccubType) ** nreflevels

          ! Transfer the weights.
          dweight = 1.0_DP/real(nsubelements,dp)
          do isubelement = 1,nsubelements
            do icubp = 1,ncubpts
              Domega((isubelement-1)*ncubpts+icubp) = DomegaLocal(icubp)*dweight
            end do
          end do
          
          ! Transfer the point coordinates.
          ! We have 2^nlevels elements in every dimension
          nsubelements = 2**nreflevels
          
          ! Length of the edge of every subelement.
          ! Note that the reference element is [-1,1]^2 !
          dedgelen = 2.0_DP/real(nsubelements,dp)
          
          do isubx = 0,nsubelements-1
            do icubp = 1,ncubpts
              Dpoints(1,isubx*ncubpts+icubp) = &
                  DpointsLocal(1,icubp)*0.5_DP*dedgelen - 1.0_DP + 0.5_DP*dedgelen &
                  + real(isubx,dp)*dedgelen
            end do
          end do
          
        case (CUB_TRZ_1D)
          ! Points are in the corners of the (sub-)elements.
          !
          ! How many subelements do we have?
          ! We have 2^nlevels elements in every dimension
          nsubelements = 2**nreflevels
          
          ! Manually initialise the coordinates and weights.
          dweight = 1.0_DP/real(cub_getRefElements(ccubType) ** nreflevels,dp)
          do isubx = 0,nsubelements
            Dpoints(1,isubx+1) = &
                real(isubx,dp)*2.0_DP/real(nsubelements,dp)-1.0_DP
          end do

          do isubx = 1,nsubelements-1
            Domega(isubx+1) = dweight*2.0_DP
          end do
          
          Domega(1)              = dweight
          Domega(nsubelements+1) = dweight

        case (CUB_G1_2D,CUB_G2_2D,CUB_G3_2D,CUB_G4_2D,CUB_G5_2D)
          ! All points are inner points. Total number =
          ! number per element * #elements in the reference element.
          ! This is a quad element in 2D, so every refinement brings 4 new
          ! elements.
          nsubelements = cub_getRefElements(ccubType) ** nreflevels

          ! Transfer the weights.
          dweight = 1.0_DP/real(nsubelements,dp)
          do isubelement = 1,nsubelements
            do icubp = 1,ncubpts
              Domega((isubelement-1)*ncubpts+icubp) = DomegaLocal(icubp)*dweight
            end do
          end do

          ! Transfer the point coordinates.
          ! We have 2^nlevels elements in every dimension
          nsubelements = 2**nreflevels
          
          ! Length of the edge of every subelement.
          ! Note that the reference element is [-1,1]^2 !
          dedgelen = 2.0_DP/real(nsubelements,dp)
          
          do isuby = 0,nsubelements-1
            do isubx = 0,nsubelements-1
              do icubp = 1,ncubpts
                Dpoints(1,(isuby*nsubelements+isubx)*ncubpts+icubp) = &
                    DpointsLocal(1,icubp)*0.5_DP*dedgelen - 1.0_DP + 0.5_DP*dedgelen &
                    + real(isubx,dp)*dedgelen

                Dpoints(2,(isuby*nsubelements+isubx)*ncubpts+icubp) = &
                    DpointsLocal(2,icubp)*0.5_DP*dedgelen - 1.0_DP + 0.5_DP*dedgelen &
                    + real(isuby,dp)*dedgelen
              end do
            end do
          end do
          
        case (CUB_TRZ_2D)
          ! Points are in the corners of the (sub-)elements.
          !
          ! How many subelements do we have?
          ! We have 2^nlevels elements in every dimension
          nsubelements = 2**nreflevels
          
          ! Manually initialise the coordinates and weights.
          dweight = 1.0_DP/real(cub_getRefElements(ccubType) ** nreflevels,dp)
          do isuby = 0,nsubelements
            do isubx = 0,nsubelements
              Dpoints(1,isuby*(nsubelements+1)+isubx+1) = &
                  real(isubx,dp)*2.0_DP/real(nsubelements,dp)-1.0_DP
              Dpoints(2,isuby*(nsubelements+1)+isubx+1) = &
                  real(isuby,dp)*2.0_DP/real(nsubelements,dp)-1.0_DP
            end do
          end do

          do isuby = 1,nsubelements-1
            do isubx = 1,nsubelements-1
              Domega(isuby*(nsubelements+1)+isubx+1) = dweight*4.0_DP
            end do
          end do
          
          do isuby = 1,nsubelements-1
            Domega(isuby+1)                               = dweight*2.0_DP
            Domega(nsubelements*(nsubelements+1)+isuby+1) = dweight*2.0_DP
            Domega(isuby*(nsubelements+1)+1)              = dweight*2.0_DP
            Domega(isuby*(nsubelements+1)+nsubelements+1) = dweight*2.0_DP
          end do

          Domega(1)                                            = dweight
          Domega(nsubelements+1)                               = dweight
          Domega(nsubelements*(nsubelements+1)+1)              = dweight
          Domega(nsubelements*(nsubelements+1)+nsubelements+1) = dweight

        case (CUB_G1_3D,CUB_G2_3D,CUB_G3_3D,CUB_G4_3D,CUB_G5_3D)
          ! All points are inner points. Total number =
          ! number per element * #elements in the reference element.
          ! This is a hexa element in 3D, so every refinement brings 8 new
          ! elements.
          nsubelements = cub_getRefElements(ccubType) ** nreflevels

          ! Transfer the weights.
          dweight = 1.0_DP/real(nsubelements,dp)
          do isubelement = 1,nsubelements
            do icubp = 1,ncubpts
              Domega((isubelement-1)*ncubpts+icubp) = DomegaLocal(icubp)*dweight
            end do
          end do
          
          ! Transfer the point coordinates.
          ! We have 2^nlevels elements in every dimension
          nsubelements = 2**nreflevels
          
          ! Length of the edge of every subelement.
          ! Note that the reference element is [-1,1]^2 !
          dedgelen = 2.0_DP/real(nsubelements,dp)
          
          do isubz = 0,nsubelements-1
            do isuby = 0,nsubelements-1
              do isubx = 0,nsubelements-1
                do icubp = 1,ncubpts
                  Dpoints(1,(isubz*nsubelements**2+&
                             isuby*nsubelements+isubx)*ncubpts+icubp) = &
                      DpointsLocal(1,icubp)*0.5_DP*dedgelen - 1.0_DP + 0.5_DP*dedgelen &
                      + real(isubx,dp)*dedgelen
                  
                  Dpoints(2,(isubz*nsubelements**2+&
                             isuby*nsubelements+isubx)*ncubpts+icubp) = &
                      DpointsLocal(2,icubp)*0.5_DP*dedgelen - 1.0_DP + 0.5_DP*dedgelen &
                      + real(isuby,dp)*dedgelen

                  Dpoints(3,(isubz*nsubelements**2+&
                             isuby*nsubelements+isubx)*ncubpts+icubp) = &
                      DpointsLocal(3,icubp)*0.5_DP*dedgelen - 1.0_DP + 0.5_DP*dedgelen &
                      + real(isubz,dp)*dedgelen
                end do
              end do
            end do
          end do

        case (CUB_TRZ_3D)
          ! Points are in the corners of the (sub-)elements.
          !
          ! How many subelements do we have?
          ! We have 2^nlevels elements in every dimension
          nsubelements = 2**nreflevels

          ! Manually initialise the coordinates and weights.
          dweight = 1.0_DP/real(cub_getRefElements(ccubType) ** nreflevels,dp)
          do isubz = 0,nsubelements
            do isuby = 0,nsubelements
              do isubx = 0,nsubelements
                Dpoints(1,isubz*(nsubelements+1)**2+&
                          isuby*(nsubelements+1)+isubx+1) = &
                    real(isubx,dp)*2.0_DP/real(nsubelements,dp)-1.0_DP
                Dpoints(2,isubz*(nsubelements+1)**2+&
                          isuby*(nsubelements+1)+isubx+1) = &
                    real(isuby,dp)*2.0_DP/real(nsubelements,dp)-1.0_DP
                Dpoints(3,isubz*(nsubelements+1)**2+&
                          isuby*(nsubelements+1)+isubx+1) = &
                    real(isubz,dp)*2.0_DP/real(nsubelements,dp)-1.0_DP
              end do
            end do
          end do
          
          ! all weights
          do isubz = 1,nsubelements-1
            do isuby = 1,nsubelements-1
              do isubx = 1,nsubelements-1
                Domega(isubz*(nsubelements+1)**2+&
                       isuby*(nsubelements+1)+isubx+1) = dweight*8.0_DP
              end do
            end do
          end do

          ! interior points at front face
          do isuby = 1,nsubelements-1
            do isubx = 1,nsubelements-1
              Domega(isuby*(nsubelements+1)+isubx+1) = dweight*4.0_DP
            end do
          end do

          ! interior points at backward face
          do isuby = 1,nsubelements-1
            do isubx = 1,nsubelements-1
              Domega(nsubelements*(nsubelements+1)**2+&
                     isuby*(nsubelements+1)+isubx+1) = dweight*4.0_DP
            end do
          end do

          ! interior points at bottom face
          do isubz = 1,nsubelements-1
            do isubx = 1,nsubelements-1
              Domega(isubz*(nsubelements+1)**2+isubx+1) = dweight*4.0_DP
            end do
          end do

          ! interior points at top face
          do isubz = 1,nsubelements-1
            do isubx = 1,nsubelements-1
              Domega(isubz*(nsubelements+1)**2+&
                     nsubelements*(nsubelements+1)+isubx+1) = dweight*4.0_DP
            end do
          end do

          ! interior points at left face
          do isubz = 1,nsubelements-1
            do isuby = 1,nsubelements-1
              Domega(isubz*(nsubelements+1)**2+&
                     isuby*(nsubelements+1)+1) = dweight*4.0_DP
            end do
          end do

          ! interior points at right face
          do isubz = 1,nsubelements-1
            do isuby = 1,nsubelements-1
              Domega(isubz*(nsubelements+1)**2+&
                     isuby*(nsubelements+1)+nsubelements+1) = dweight*4.0_DP
            end do
          end do

          ! boundary points at front face
          do isuby = 1,nsubelements-1
            Domega(isuby+1)                               = dweight*2.0_DP
            Domega(nsubelements*(nsubelements+1)+isuby+1) = dweight*2.0_DP
            Domega(isuby*(nsubelements+1)+1)              = dweight*2.0_DP
            Domega(isuby*(nsubelements+1)+nsubelements+1) = dweight*2.0_DP
          end do

          ! boundary points at backward face
          do isuby = 1,nsubelements-1
            Domega(nsubelements*(nsubelements+1)**2+&
                   isuby+1)                               = dweight*2.0_DP
            Domega(nsubelements*(nsubelements+1)**2+&
                   nsubelements*(nsubelements+1)+isuby+1) = dweight*2.0_DP
            Domega(nsubelements*(nsubelements+1)**2+&
                   isuby*(nsubelements+1)+1)              = dweight*2.0_DP
            Domega(nsubelements*(nsubelements+1)**2+&
                   isuby*(nsubelements+1)+nsubelements+1) = dweight*2.0_DP
          end do
          
          ! boundary points at missing face in z-direction
          do isubz = 1,nsubelements-1
            Domega(isubz*(nsubelements+1)**2+1)              = dweight*2.0_DP
            Domega(isubz*(nsubelements+1)**2+&
                          nsubelements*(nsubelements+1)+1)   = dweight*2.0_DP
            Domega(isubz*(nsubelements+1)**2+nsubelements+1) = dweight*2.0_DP
            Domega(isubz*(nsubelements+1)**2+&
                          nsubelements*(nsubelements+1)+&
                          nsubelements+1)                    = dweight*2.0_DP
          end do
          
          ! corner points
          Domega(1)                                            = dweight
          Domega(nsubelements+1)                               = dweight
          Domega(nsubelements*(nsubelements+1)+1)              = dweight
          Domega(nsubelements*(nsubelements+1)+nsubelements+1) = dweight
          Domega(nsubelements*(nsubelements+1)**2+1)           = dweight
          Domega(nsubelements*(nsubelements+1)**2+&
                 nsubelements+1)                               = dweight
          Domega(nsubelements*(nsubelements+1)**2+&
                 nsubelements*(nsubelements+1)+1)              = dweight
          Domega(nsubelements*(nsubelements+1)**2+&
                 nsubelements*(nsubelements+1)+1)              = dweight
          Domega(nsubelements*(nsubelements+1)**2+&
                 nsubelements*(nsubelements+1)+nsubelements+1) = dweight

      case (CUB_G3_T)
        ! Points are in the midpoints of the (sub-)elements.
        !
        ! How many subelements do we have?
        ! We have 2^nlevels elements in every barycentric dimension
        nsubelements = 2**nreflevels

        ! Manually initialise the coordinates.
        icubp = 1
        do isubz = 0, 2*nsubelements
          if (mod(isubz,2).eq.0) then
            do isuby = 1, 2*nsubelements-isubz, 2
              Dpoints(1,icubp) = 1.0_DP-real(isuby,dp)/real(2*nsubelements,dp)-&
                                        real(isubz,dp)/real(2*nsubelements,dp)
              Dpoints(2,icubp) = real(isuby,dp)/real(2*nsubelements,dp)
              Dpoints(3,icubp) = real(isubz,dp)/real(2*nsubelements,dp)
              icubp = icubp+1
            end do
          else
            do isuby = 0, 2*nsubelements-isubz
              Dpoints(1,icubp) = 1.0_DP-real(isuby,dp)/real(2*nsubelements,dp)-&
                                        real(isubz,dp)/real(2*nsubelements,dp)
              Dpoints(2,icubp) = real(isuby,dp)/real(2*nsubelements,dp)
              Dpoints(3,icubp) = real(isubz,dp)/real(2*nsubelements,dp)
              icubp = icubp+1
            end do
          end if
        end do
        
        ! Store number of cubature points on refined element
        nrefcubpts = icubp-1

        ! Manually initialise the weights.
        dweight = 1.0_DP/(6.0_DP*real(nsubelements**2,dp))

        ! all weights
        do icubp = 1, nrefcubpts
          Domega(icubp) = dweight*2.0_DP
        end do

        ! weights at the boundary
        do isubelement = 1, nsubelements
          Domega(isubelement)                         = dweight
          Domega(3*nsubelements*(isubelement+1)-&
                 (3*(isubelement+1)**2)/2+&
                 (9*isubelement+3)/2-3*nsubelements)  = dweight
          Domega(3*nsubelements*(isubelement+1)-&
                 (3*(isubelement+1)**2)/2+&
                 (13*isubelement+1)/2-5*nsubelements) = dweight
        end do

      case (CUB_G1_T,CUB_G3MP_T)
        ! All points are inner points. Total number =
        ! number per element * #elements in the reference element.
        ! This is a triangle element in 2D, so every refinement brings 4 new
        ! elements.
        nsubelements = cub_getRefElements(ccubType) ** nreflevels
        dweight = 1.0_DP/real(nsubelements,dp)

        ! Initialise the number of cubature points
        nrefcubpts = 0

        ! Initialise the coordinates of the reference triangle
        allocate(DpointsRef(3,3))
        DpointsRef(:,1) = (/1.0_DP, 0.0_DP, 0.0_DP/)
        DpointsRef(:,2) = (/0.0_DP, 1.0_DP, 0.0_DP/)
        DpointsRef(:,3) = (/0.0_DP, 0.0_DP, 1.0_DP/)

        ! Recursively initialise the coordinates and weights
        call cub_auxTriangle(DpointsRef, DpointsLocal, DomegaLocal,&
            dweight, ncubpts, nreflevels, Dpoints, Domega, nrefcubpts)
        
        ! Deallocate temporal memory
        deallocate(DpointsRef)

      case (CUB_TRZ_T)
        ! Points are in the corners of the (sub-)elements.
        !
        ! How many subelements do we have?
        ! We have 2^nlevels elements in every barycentric dimension
        nsubelements = 2**nreflevels

        ! Manually initialise the coordinates.
        icubp = 1
        do isubz = 0, nsubelements
          do isuby = 0, nsubelements-isubz
            Dpoints(1,icubp) = 1.0_DP-real(isuby,dp)/real(nsubelements,dp)-&
                                      real(isubz,dp)/real(nsubelements,dp)
            Dpoints(2,icubp) = real(isuby,dp)/real(nsubelements,dp)
            Dpoints(3,icubp) = real(isubz,dp)/real(nsubelements,dp)
            icubp = icubp+1
          end do
        end do

        ! Store number of cubature points on refined element
        nrefcubpts = icubp-1

        ! Manually initialise the weights.
        dweight = 1.0_DP/(6.0_DP*real(nsubelements**2,dp))

        ! all weights
        do icubp = 1, nrefcubpts
          Domega(icubp) = dweight*6.0_DP
        end do

        ! weights at the boundary
        do isubelement = 1, nsubelements
          Domega(isubelement)                                       = dweight*3.0_DP
          Domega(nsubelements*(isubelement+1)-((isubelement+1)**2)/2+&
                 (5*isubelement+1)/2-nsubelements)                  = dweight*3.0_DP
          Domega(nsubelements*(isubelement+1)-((isubelement+1)**2)/2+&
                 (5*isubelement+1)/2-nsubelements+1)                = dweight*3.0_DP
        end do

        ! weights at the corners
        Domega(1)              = dweight
        Domega(nsubelements+1) = dweight
        Domega(nrefcubpts)     = dweight

      case (CUB_VMC)
        ! Points are in the corners, midpoints and center of the (sub-)elements.
        !
        ! How many subelements do we have?
        ! We have 2^nlevels elements in every barycentric dimension
        nsubelements = 2**nreflevels
        
        ! Manually initialise the coordinates of the corners.
        icubp = 1
        do isubz = 0, nsubelements
          do isuby = 0, nsubelements-isubz
            Dpoints(1,icubp) = 1.0_DP-real(isuby,dp)/real(nsubelements,dp)-&
                                      real(isubz,dp)/real(nsubelements,dp)
            Dpoints(2,icubp) = real(isuby,dp)/real(nsubelements,dp)
            Dpoints(3,icubp) = real(isubz,dp)/real(nsubelements,dp)
            icubp = icubp+1
          end do
        end do

        ! Store number of cubature points on refined element
        nrefcubpts = icubp-1

        ! Manually initialise the weights of the corners.
        dweight = 0.025_DP/real(nsubelements**2,dp)

        ! all weights
        do icubp = 1, nrefcubpts
          Domega(icubp) = dweight*6.0_DP
        end do

        ! weights at the boundary
        do isubelement = 1, nsubelements
          Domega(isubelement)                                       = dweight*3.0_DP
          Domega(nsubelements*(isubelement+1)-((isubelement+1)**2)/2+&
                 (5*isubelement+1)/2-nsubelements)                  = dweight*3.0_DP
          Domega(nsubelements*(isubelement+1)-((isubelement+1)**2)/2+&
                 (5*isubelement+1)/2-nsubelements+1)                = dweight*3.0_DP
        end do

        ! weights at the corners
        Domega(1)              = dweight
        Domega(nsubelements+1) = dweight
        Domega(nrefcubpts)     = dweight

        ! Manually initialise the coordinates of the midpoints.
        do isubz = 0, 2*nsubelements
          if (mod(isubz,2).eq.0) then
            do isuby = 1, 2*nsubelements-isubz, 2
              Dpoints(1,icubp) = 1.0_DP-real(isuby,dp)/real(2*nsubelements,dp)-&
                                        real(isubz,dp)/real(2*nsubelements,dp)
              Dpoints(2,icubp) = real(isuby,dp)/real(2*nsubelements,dp)
              Dpoints(3,icubp) = real(isubz,dp)/real(2*nsubelements,dp)
              icubp = icubp+1
            end do
          else
            do isuby = 0, 2*nsubelements-isubz
              Dpoints(1,icubp) = 1.0_DP-real(isuby,dp)/real(2*nsubelements,dp)-&
                                        real(isubz,dp)/real(2*nsubelements,dp)
              Dpoints(2,icubp) = real(isuby,dp)/real(2*nsubelements,dp)
              Dpoints(3,icubp) = real(isubz,dp)/real(2*nsubelements,dp)
              icubp = icubp+1
            end do
          end if
        end do
        
        ! Store number of cubature points on refined element
        icubpStart = nrefcubpts;   nrefcubpts = icubp-1

        ! Manually initialise the weights.
        dweight = 1.0_DP/(15.0_DP*real(nsubelements**2,dp))

        ! all weights
        do icubp = icubpStart+1, nrefcubpts
          Domega(icubp) = dweight*2.0_DP
        end do

        ! weights at the boundary
        do isubelement = 1, nsubelements
          Domega(icubpStart+isubelement)              = dweight
          Domega(icubpStart+3*nsubelements*(isubelement+1)-&
                 (3*(isubelement+1)**2)/2+&
                 (9*isubelement+3)/2-3*nsubelements)  = dweight
          Domega(icubpStart+3*nsubelements*(isubelement+1)-&
                 (3*(isubelement+1)**2)/2+&
                 (13*isubelement+1)/2-5*nsubelements) = dweight
        end do

        ! Initialise the inner cubature points
        nsubelements = cub_getRefElements(ccubType) ** nreflevels
        dweight = 1.0_DP/real(nsubelements,dp)

        ! Initialise the coordinates of the reference triangle
        allocate(DpointsRef(3,3))
        DpointsRef(:,1) = (/1.0_DP, 0.0_DP, 0.0_DP/)
        DpointsRef(:,2) = (/0.0_DP, 1.0_DP, 0.0_DP/)
        DpointsRef(:,3) = (/0.0_DP, 0.0_DP, 1.0_DP/)
        
        ! Initialise the local cubature point
        DpointsLocal(1,1) = 1.0_DP/3.0_DP
        DpointsLocal(2,1) = 1.0_DP/3.0_DP
        DpointsLocal(3,1) = 1.0_DP/3.0_DP

        ! Initialise the local weight
        DomegaLocal(1) = 0.225_DP

        ! Recursively initialise the coordinates and weights at the centroids
        call cub_auxTriangle(DpointsRef, DpointsLocal, DomegaLocal,&
            dweight, 1, nreflevels, Dpoints, Domega, nrefcubpts)
      
        ! Deallocate temporal memory
        deallocate(DpointsRef)

      case default
        call output_line ('Unsupported element.', &
                          OU_CLASS_ERROR,OU_MODE_STD,'cub_getCubature')
        call sys_halt()
      end select
        
      ! Release memory, that is it.
      
      deallocate(DomegaLocal)
      deallocate(DpointsLocal)
    
    end if
    
    ! That is it
    
  contains
  
    ! -------------------------------------------------------
  
    pure subroutine cub_auxGaussLegendre(n,Dv,Dw)
    integer, intent(in) :: n
    real(DP), dimension(:), intent(out) :: Dv, Dw
    
    real(DP) :: daux
    
      ! How many points for Gauss-formula?
      select case(n)
      case(1)
        Dv(1) =  0.0_DP
        Dw(1) =  2.0_DP
      
      case(2)
        Dv(1) = -sqrt(1.0_DP / 3.0_DP)
        Dv(2) =  sqrt(1.0_DP / 3.0_DP)
        Dw(1) =  1.0_DP
        Dw(2) =  1.0_DP
      
      case(3)
        Dv(1) = -sqrt(0.6_DP)
        Dv(2) =  0.0_DP
        Dv(3) =  sqrt(0.6_DP)
        Dw(1) =  5.0_DP / 9.0_DP
        Dw(2) =  8.0_DP / 9.0_DP
        Dw(3) =  5.0_DP / 9.0_DP
      
      case(4)
        daux = sqrt(4.8_DP)
        Dv(1) = -sqrt((3.0_DP + daux) / 7.0_DP)
        Dv(2) = -sqrt((3.0_DP - daux) / 7.0_DP)
        Dv(3) =  sqrt((3.0_DP - daux) / 7.0_DP)
        Dv(4) =  sqrt((3.0_DP + daux) / 7.0_DP)
        daux = sqrt(30.0_DP)
        Dw(1) = (18.0_DP - daux) / 36.0_DP
        Dw(2) = (18.0_DP + daux) / 36.0_DP
        Dw(3) = (18.0_DP + daux) / 36.0_DP
        Dw(4) = (18.0_DP - daux) / 36.0_DP
      
      case(5)
        daux = 2.0_DP * sqrt(10.0_DP / 7.0_DP)
        Dv(1) = -sqrt(5.0_DP + daux) / 3.0_DP
        Dv(2) = -sqrt(5.0_DP - daux) / 3.0_DP
        Dv(3) =  0.0_DP
        Dv(4) =  sqrt(5.0_DP - daux) / 3.0_DP
        Dv(5) =  sqrt(5.0_DP + daux) / 3.0_DP
        daux = 13.0_DP * sqrt(70.0_DP)
        Dw(1) = (322.0_DP - daux) / 900.0_DP
        Dw(2) = (322.0_DP + daux) / 900.0_DP
        Dw(3) =  128.0_DP / 225.0_DP
        Dw(4) = (322.0_DP + daux) / 900.0_DP
        Dw(5) = (322.0_DP - daux) / 900.0_DP

      case(6)
        Dv(1) = -0.932469514203152_DP
        Dv(2) = -0.661209386466265_DP
        Dv(3) = -0.238619186083197_DP
        Dv(4) =  0.238619186083197_DP
        Dv(5) =  0.661209386466265_DP
        Dv(6) =  0.932469514203152_DP
        Dw(1) =  0.171324492379170_DP
        Dw(2) =  0.360761573048139_DP
        Dw(3) =  0.467913934572691_DP
        Dw(4) =  0.467913934572691_DP
        Dw(5) =  0.360761573048139_DP
        Dw(6) =  0.171324492379170_DP
    
      end select
      
    end subroutine cub_auxGaussLegendre

    ! -------------------------------------------------------
    
    recursive subroutine cub_auxTriangle(DpointsRef, DpointsLocal,&
        DomegaLocal, dweight, ncubpts, ireflevel, Dpoints, Domega, nrefcubpts)
      real(DP), dimension(:,:), intent(in) :: DpointsRef
      real(DP), dimension(:,:), intent(in) :: DpointsLocal
      real(DP), dimension(:), intent(in)   :: DomegaLocal
      real(DP), intent(in) :: dweight
      integer, intent(in)  :: ncubpts, ireflevel
      real(DP), dimension(:,:), intent(inout) :: Dpoints
      real(DP), dimension(:), intent(inout)   :: Domega
      integer, intent(inout) :: nrefcubpts

      real(DP), dimension(3,3) :: DpointsAux

      if (ireflevel .eq. 0) then

        ! Manually initialise the coordinates and weights for the (sub-)element
        do icubp = 1, ncubpts
          Dpoints(:,nrefcubpts+icubp) = DpointsLocal(1,icubp)*DpointsRef(:,1)+&
                                        DpointsLocal(2,icubp)*DpointsRef(:,2)+&
                                        DpointsLocal(3,icubp)*DpointsRef(:,3)
          Domega(nrefcubpts+icubp)    = DomegaLocal(icubp)*dweight
        end do
        nrefcubpts = nrefcubpts+ncubpts

      else
        
        ! Recursively process the subelement in the center
        DpointsAux(:,1) = (DpointsRef(:,1)+DpointsRef(:,2))/2.0_DP
        DpointsAux(:,2) = (DpointsRef(:,1)+DpointsRef(:,3))/2.0_DP
        DpointsAux(:,3) = (DpointsRef(:,2)+DpointsRef(:,3))/2.0_DP
        call cub_auxTriangle(DpointsAux, DpointsLocal, DomegaLocal,&
            dweight, ncubpts, ireflevel-1, Dpoints, Domega, nrefcubpts)

        ! Recursively process the first outer subelement
        DpointsAux(:,1) =  DpointsRef(:,1)
        DpointsAux(:,2) = (DpointsRef(:,1)+DpointsRef(:,2))/2.0_DP
        DpointsAux(:,3) = (DpointsRef(:,1)+DpointsRef(:,3))/2.0_DP
        call cub_auxTriangle(DpointsAux, DpointsLocal, DomegaLocal,&
            dweight, ncubpts, ireflevel-1, Dpoints, Domega, nrefcubpts)

        ! Recursively process the second outer subelement
        DpointsAux(:,1) =  DpointsRef(:,2)
        DpointsAux(:,2) = (DpointsRef(:,2)+DpointsRef(:,1))/2.0_DP
        DpointsAux(:,3) = (DpointsRef(:,2)+DpointsRef(:,3))/2.0_DP
        call cub_auxTriangle(DpointsAux, DpointsLocal, DomegaLocal,&
            dweight, ncubpts, ireflevel-1, Dpoints, Domega, nrefcubpts)
        
        ! Recursively process the third outer subelement
        DpointsAux(:,1) =  DpointsRef(:,3)
        DpointsAux(:,2) = (DpointsRef(:,3)+DpointsRef(:,1))/2.0_DP
        DpointsAux(:,3) = (DpointsRef(:,3)+DpointsRef(:,2))/2.0_DP
        call cub_auxTriangle(DpointsAux, DpointsLocal, DomegaLocal,&
            dweight, ncubpts, ireflevel-1, Dpoints, Domega, nrefcubpts)
        
      end if

    end subroutine cub_auxTriangle

  end subroutine cub_getCubature

  !****************************************************************************

!<subroutine>

  subroutine cub_getCubPoints(ccubType, ncubp, Dxi, Domega)

!<description>
  ! This routine initialises the coordinates and weight fields according
  ! to the selected cubature formula. The integration domain is <tex>$[-1,1]^n$</tex>.
  ! In the case of one-dimensional integration, only cub_dxi(i,1) is used.
  ! The coordinates of the cubature points on triangles are given in
  ! barycentric coordinates.
!</description>

!<input>
  ! id of the cubature formula to be set
  integer(I32), intent(in) :: ccubType
!</input>
  
!<output>
  ! number of cubature points; =0: error, unknown cubature formula
  integer , intent(out) :: ncubp
  
  ! Coordinates of the cubature points.
  ! 1D: Dxi(1..ncubp,1)=coordinates,
  ! 2D: Quadrilaterals:
  !        Dxi(1..ncubp,1)=x-coord,
  !        Dxi(1..ncubp,2)=y-coord
  !     Triangles:
  !        Dxi(1..ncubp,1)=1st barycentric coordinate,
  !        Dxi(1..ncubp,2)=2nd barycentric coordinate,
  !        Dxi(1..ncubp,3)=3rd barycentric coordinate
  ! 3D: Hexahedra:
  !       Dxi(1..ncubp,1)=x-coord,
  !       Dxi(1..ncubp,2)=y-coord,
  !       Dxi(1..ncubp,3)=z-coord
  !     Tetrahedra:
  !        Dxi(1..ncubp,1)=1st barycentric coordinate,
  !        Dxi(1..ncubp,2)=2nd barycentric coordinate,
  !        Dxi(1..ncubp,3)=3rd barycentric coordinate,
  !        Dxi(1..ncubp,4)=4th barycentric coordinate
  real(DP), dimension(:,:), intent(out) :: Dxi
  
  ! For every cubature point the corresponding cubature weight
  real(DP), dimension(:), intent(out) :: Domega
  
!</output>

!</subroutine>

  ! local variables
  integer :: i
  
  select case (ccubType)

  !1D cubature formulas
  case (CUB_G1_1D)
    Dxi(1,1)  = 0.0_DP
     
    Domega(1) = 2.0_DP
     
    ncubp     = 1
     
  case(CUB_TRZ_1D)
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  =  1.0_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    
    ncubp     =  2
     
  case(CUB_G2_1D)
    Dxi(1,1)  = -0.577350269189626_DP
    Dxi(2,1)  =  0.577350269189626_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    
    ncubp     =  2
     
  case(CUB_G3_1D)
    Dxi(1,1)  = -0.774596669241483_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(3,1)  =  0.774596669241483_DP
    
    Domega(1) =  0.555555555555556_DP
    Domega(2) =  0.888888888888889_DP
    Domega(3) =  0.555555555555556_DP
    
    ncubp     =  3
    
  case(CUB_G4_1D)
    Dxi(1,1)  = -0.861136311594053_DP
    Dxi(2,1)  = -0.339981043584856_DP
    Dxi(3,1)  =  0.339981043584856_DP
    Dxi(4,1)  =  0.861136311594053_DP
    
    Domega(1) =  0.347854845137454_DP
    Domega(2) =  0.652145154862546_DP
    Domega(3) =  0.652145154862546_DP
    Domega(4) =  0.347854845137454_DP
    
    ncubp     =  4
    
  case(CUB_G5_1D)
    Dxi(1,1)  = -0.906179845938664_DP
    Dxi(2,1)  = -0.538469310105683_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(4,1)  =  0.538469310105683_DP
    Dxi(5,1)  =  0.906179845938664_DP
    
    Domega(1) =  0.236926885056189_DP
    Domega(2) =  0.478628670499366_DP
    Domega(3) =  0.568888888888889_DP
    Domega(4) =  0.478628670499366_DP
    Domega(5) =  0.236926885056189_DP
    
    ncubp = 5
    
  case(CUB_SIMPSON_1D)
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(3,1)  =  1.0_DP
    
    Domega(1) =  0.333333333333333_DP
    Domega(2) =  1.333333333333333_DP
    Domega(3) =  0.333333333333333_DP
    
    ncubp     =  3
    
  case(CUB_PULCHERIMA_1D)
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  = -0.333333333333333_DP
    Dxi(3,1)  =  0.333333333333333_DP
    Dxi(4,1)  =  1.0_DP
    
    Domega(1) =  0.25_DP
    Domega(2) =  0.75_DP
    Domega(3) =  0.75_DP
    Domega(4) =  0.25_DP
    
    ncubp     =  4

  case(CUB_MILNE_1D)
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  = -0.5_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(4,1)  =  0.5_DP
    Dxi(5,1)  =  1.0_DP
    
    Domega(1) =  0.155555555555555_DP
    Domega(2) =  0.711111111111111_DP
    Domega(3) =  0.266666666666666_DP
    Domega(4) =  0.711111111111111_DP
    Domega(5) =  0.155555555555555_DP
    
    ncubp     =  5

  case(CUB_6POINT_1D)
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  = -0.6_DP
    Dxi(3,1)  = -0.2_DP
    Dxi(4,1)  =  0.2_DP
    Dxi(5,1)  =  0.6_DP
    Dxi(6,1)  =  1.0_DP
    
    Domega(1) =  0.131944444444444_DP
    Domega(2) =  0.520833333333333_DP
    Domega(3) =  0.347222222222222_DP
    Domega(4) =  0.347222222222222_DP
    Domega(5) =  0.520833333333333_DP
    Domega(6) =  0.131944444444444_DP
    
    ncubp     =  6

  case(CUB_WEDDLE_1D)
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  = -0.666666666666666_DP
    Dxi(3,1)  = -0.333333333333333_DP
    Dxi(4,1)  =  0.0_DP
    Dxi(5,1)  =  0.333333333333333_DP
    Dxi(6,1)  =  0.666666666666666_DP
    Dxi(7,1)  =  1.0_DP
    
    Domega(1) =  0.0976190476190476_DP
    Domega(2) =  0.514285714285714_DP
    Domega(3) =  0.0642857142857143_DP
    Domega(4) =  0.647619047619048_DP
    Domega(5) =  0.0642857142857143_DP
    Domega(6) =  0.514285714285714_DP
    Domega(7) =  0.0976190476190476_DP
    
    ncubp     =  7

  case(CUB_G6_1D)
    Dxi(1,1)  = -0.932469514203152_DP
    Dxi(2,1)  = -0.661209386466265_DP
    Dxi(3,1)  = -0.238619186083197_DP
    Dxi(4,1)  =  0.238619186083197_DP
    Dxi(5,1)  =  0.661209386466265_DP
    Dxi(6,1)  =  0.932469514203152_DP
    
    Domega(1) =  0.171324492379170_DP
    Domega(2) =  0.360761573048139_DP
    Domega(3) =  0.467913934572691_DP
    Domega(4) =  0.467913934572691_DP
    Domega(5) =  0.360761573048139_DP
    Domega(6) =  0.171324492379170_DP
    
    ncubp     =  6
    
  !2D cubature formulas
  case (CUB_G1X1)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.0_DP
    
    Domega(1) =  4.0_DP
    
    ncubp     =  1
    
  case (CUB_TRZ)
    Dxi(1,1)  = -1.0_DP
    Dxi(1,2)  = -1.0_DP
    Dxi(2,1)  =  1.0_DP
    Dxi(2,2)  = -1.0_DP
    Dxi(3,1)  =  1.0_DP
    Dxi(3,2)  =  1.0_DP
    Dxi(4,1)  = -1.0_DP
    Dxi(4,2)  =  1.0_DP
    
    Domega(1) = 1.0_DP
    Domega(2) = 1.0_DP
    Domega(3) = 1.0_DP
    Domega(4) = 1.0_DP
    
    ncubp     = 4
    
  case (CUB_MID)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  = -1.0_DP
    Dxi(2,1)  =  1.0_DP
    Dxi(2,2)  =  0.0_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  =  1.0_DP
    Dxi(4,1)  = -1.0_DP
    Dxi(4,2)  =  0.0_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    Domega(3) =  1.0_DP
    Domega(4) =  1.0_DP
    
    ncubp     =  4
     
  case (CUB_G2X2)
    Dxi(1,1)  =  0.577350269189626_DP
    Dxi(1,2)  =  0.577350269189626_DP
    Dxi(2,1)  = -0.577350269189626_DP
    Dxi(2,2)  =  0.577350269189626_DP
    Dxi(3,1)  = -0.577350269189626_DP
    Dxi(3,2)  = -0.577350269189626_DP
    Dxi(4,1)  =  0.577350269189626_DP
    Dxi(4,2)  = -0.577350269189626_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    Domega(3) =  1.0_DP
    Domega(4) =  1.0_DP
    
    ncubp     =  4
     
  case (CUB_NS1)
    Dxi(1,1)  =  0.816496580927726_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(2,1)  = -0.816496580927726_DP
    Dxi(2,2)  =  0.0_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  = -0.816496580927726_DP
    Dxi(4,1)  =  0.0_DP
    Dxi(4,2)  =  0.816496580927726_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    Domega(3) =  1.0_DP
    Domega(4) =  1.0_DP
    
    ncubp     =  4
    
  case (CUB_NS2)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.934172358962716_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(2,2)  = -0.356822089773090_DP
    Dxi(3,1)  =  0.774596669241483_DP
    Dxi(3,2)  =  0.390885162530070_DP
    Dxi(4,1)  = -0.774596669241483_DP
    Dxi(4,2)  =  0.390885162530070_DP
    Dxi(5,1)  =  0.774596669241483_DP
    Dxi(5,2)  = -0.852765377881771_DP
    Dxi(6,1)  = -0.774596669241483_DP
    Dxi(6,2)  = -0.852765377881771_DP
    
    Domega(1) =  0.491365692888926_DP
    Domega(2) =  1.28641208488885_DP
    Domega(3) =  0.761883709085613_DP
    Domega(4) =  0.761883709085613_DP
    Domega(5) =  0.349227402025498_DP
    Domega(6) =  0.349227402025498_DP
    
    ncubp     =  6

  case (CUB_NS3)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(2,2)  =  0.9660917830792960_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  = -0.966091783079296_DP
    Dxi(4,1)  =  0.774596669241483_DP
    Dxi(4,2)  =  0.774596669241483_DP
    Dxi(5,1)  = -0.774596669241483_DP
    Dxi(5,2)  =  0.774596669241483_DP
    Dxi(6,1)  =  0.774596669241483_DP
    Dxi(6,2)  = -0.774596669241483_DP
    Dxi(7,1)  = -0.774596669241483_DP
    Dxi(7,2)  = -0.774596669241483_DP
    
    Domega(1) =  1.142857142857143_DP
    Domega(2) =  0.317460317460317_DP
    Domega(3) =  0.317460317460317_DP
    Domega(4) =  0.555555555555556_DP
    Domega(5) =  0.555555555555556_DP
    Domega(6) =  0.555555555555556_DP
    Domega(7) =  0.555555555555556_DP
    
    ncubp     =  7

  case (CUB_G3X3)
    Dxi(1,1)  =  0.774596669241483_DP
    Dxi(1,2)  =  0.774596669241483_DP
    Dxi(2,1)  = -0.774596669241483_DP
    Dxi(2,2)  =  0.774596669241483_DP
    Dxi(3,1)  =  0.774596669241483_DP
    Dxi(3,2)  = -0.774596669241483_DP
    Dxi(4,1)  = -0.774596669241483_DP
    Dxi(4,2)  = -0.774596669241483_DP
    Dxi(5,1)  =  0.774596669241483_DP
    Dxi(5,2)  =  0.0_DP
    Dxi(6,1)  = -0.774596669241483_DP
    Dxi(6,2)  =  0.0_DP
    Dxi(7,1)  =  0.0_DP
    Dxi(7,2)  =  0.774596669241483_DP
    Dxi(8,1)  =  0.0_DP
    Dxi(8,2)  = -0.774596669241483_DP
    Dxi(9,1)  =  0.0_DP
    Dxi(9,2)  =  0.0_DP
    
    Domega(1) =  0.308641975308642_DP
    Domega(2) =  0.308641975308642_DP
    Domega(3) =  0.308641975308642_DP
    Domega(4) =  0.308641975308642_DP
    Domega(5) =  0.493827160493827_DP
    Domega(6) =  0.493827160493827_DP
    Domega(7) =  0.493827160493827_DP
    Domega(8) =  0.493827160493827_DP
    Domega(9) =  0.790123456790123_DP
    
    ncubp     =  9
    
  case (CUB_G)
    Dxi(1,1)  =  0.92582009977255146_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(2,1)  = -0.92582009977255146_DP
    Dxi(2,2)  =  0.0_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  =  0.92582009977255146_DP
    Dxi(4,1)  =  0.0_DP
    Dxi(4,2)  = -0.92582009977255146_DP
    
    Dxi(5,1)  =  0.38055443320831566_DP
    Dxi(5,2)  =  0.38055443320831566_DP
    Dxi(6,1)  =  0.38055443320831566_DP
    Dxi(6,2)  = -0.38055443320831566_DP
    Dxi(7,1)  = -0.38055443320831566_DP
    Dxi(7,2)  =  0.38055443320831566_DP
    Dxi(8,1)  = -0.38055443320831566_DP
    Dxi(8,2)  = -0.38055443320831566_DP
    
    Dxi(9,1)  =  0.80597978291859874_DP
    Dxi(9,2)  =  0.80597978291859874_DP
    Dxi(10,1) =  0.80597978291859874_DP
    Dxi(10,2) = -0.80597978291859874_DP
    Dxi(11,1) = -0.80597978291859874_DP
    Dxi(11,2) =  0.80597978291859874_DP
    Dxi(12,1) = -0.80597978291859874_DP
    Dxi(12,2) = -0.80597978291859874_DP
    
    Domega(1) =  0.24197530864197531_DP
    Domega(2) =  0.24197530864197531_DP
    Domega(3) =  0.24197530864197531_DP
    Domega(4) =  0.24197530864197531_DP
    
    Domega(5) =  0.52059291666739446_DP
    Domega(6) =  0.52059291666739446_DP
    Domega(7) =  0.52059291666739446_DP
    Domega(8) =  0.52059291666739446_DP
    
    Domega(9) =  0.23743177469063023_DP
    Domega(10)=  0.23743177469063023_DP
    Domega(11)=  0.23743177469063023_DP
    Domega(12)=  0.23743177469063023_DP
    
    ncubp     =  12
    
  case (CUB_G4X4)
    Dxi(1,1)  =  0.861136311594053003_DP
    Dxi(1,2)  =  0.339981043584855994_DP
    Dxi(2,1)  = -0.861136311594053003_DP
    Dxi(2,2)  = -0.339981043584855994_DP
    Dxi(3,1)  = -0.861136311594053003_DP
    Dxi(3,2)  =  0.339981043584855994_DP
    Dxi(4,1)  =  0.861136311594053003_DP
    Dxi(4,2)  = -0.339981043584855994_DP
    Dxi(5,1)  = -0.339981043584855994_DP
    Dxi(5,2)  = -0.861136311594053003_DP
    Dxi(6,1)  =  0.339981043584855994_DP
    Dxi(6,2)  =  0.861136311594053003_DP
    Dxi(7,1)  =  0.339981043584855994_DP
    Dxi(7,2)  = -0.861136311594053003_DP
    Dxi(8,1)  = -0.339981043584855994_DP
    Dxi(8,2)  =  0.861136311594053003_DP
    
    Dxi(9,1)  = -0.339981043584855994_DP
    Dxi(9,2)  =  0.339981043584855994_DP
    Dxi(10,1) =  0.339981043584855994_DP
    Dxi(10,2) = -0.339981043584855994_DP
    Dxi(11,1) =  0.339981043584855994_DP
    Dxi(11,2) =  0.339981043584855994_DP
    Dxi(12,1) = -0.339981043584855994_DP
    Dxi(12,2) = -0.339981043584855994_DP
    
    Dxi(13,1) =  0.861136311594053003_DP
    Dxi(13,2) = -0.861136311594053003_DP
    Dxi(14,1) = -0.861136311594053003_DP
    Dxi(14,2) =  0.861136311594053003_DP
    Dxi(15,1) = -0.861136311594053003_DP
    Dxi(15,2) = -0.861136311594053003_DP
    Dxi(16,1) =  0.861136311594053003_DP
    Dxi(16,2) =  0.861136311594053003_DP
    
    Domega(1) =  0.226851851851851888_DP
    Domega(2) =  0.226851851851851888_DP
    Domega(3) =  0.226851851851851888_DP
    Domega(4) =  0.226851851851851888_DP
    Domega(5) =  0.226851851851851888_DP
    Domega(6) =  0.226851851851851888_DP
    Domega(7) =  0.226851851851851888_DP
    Domega(8) =  0.226851851851851888_DP
    
    Domega(9) =  0.425293303010694096_DP
    Domega(10)=  0.425293303010694096_DP
    Domega(11)=  0.425293303010694096_DP
    Domega(12)=  0.425293303010694096_DP
    Domega(13)=  0.121002993285602101_DP
    Domega(14)=  0.121002993285602101_DP
    Domega(15)=  0.121002993285602101_DP
    Domega(16)=  0.121002993285602101_DP
    
    ncubp     =  16

  case (CUB_G5X5)
    Dxi(1,1)  =  0.906179845938664005_DP
    Dxi(1,2)  =  0.906179845938664005_DP
    Dxi(2,1)  =  0.906179845938664005_DP
    Dxi(2,2)  = -0.906179845938664005_DP
    Dxi(3,1)  = -0.906179845938664005_DP
    Dxi(3,2)  =  0.906179845938664005_DP
    Dxi(4,1)  = -0.906179845938664005_DP
    Dxi(4,2)  = -0.906179845938664005_DP
    
    Dxi(5,1)  = -0.906179845938664005_DP
    Dxi(5,2)  = -0.538469310105682997_DP
    Dxi(6,1)  = -0.538469310105682997_DP
    Dxi(6,2)  = -0.906179845938664005_DP
    Dxi(7,1)  =  0.906179845938664005_DP
    Dxi(7,2)  =  0.538469310105682997_DP
    Dxi(8,1)  = -0.906179845938664005_DP
    Dxi(8,2)  =  0.538469310105682997_DP
    Dxi(9,1)  =  0.538469310105682997_DP
    Dxi(9,2)  = -0.906179845938664005_DP
    Dxi(10,1) = -0.538469310105682997_DP
    Dxi(10,2) =  0.906179845938664005_DP
    Dxi(11,1) =  0.538469310105682997_DP
    Dxi(11,2) =  0.906179845938664005_DP
    Dxi(12,1) =  0.906179845938664005_DP
    Dxi(12,2) = -0.538469310105682997_DP
    
    Dxi(13,1) =  0.538469310105682997_DP
    Dxi(13,2) =  0.0_DP
    Dxi(14,1) =  0.0_DP
    Dxi(14,2) =  0.538469310105682997_DP
    Dxi(15,1) = -0.538469310105682997_DP
    Dxi(15,2) =  0.0_DP
    Dxi(16,1) =  0.0_DP
    Dxi(16,2) = -0.538469310105682997_DP
    
    Dxi(17,1) =  0.538469310105682997_DP
    Dxi(17,2) =  0.538469310105682997_DP
    Dxi(18,1) =  0.538469310105682997_DP
    Dxi(18,2) = -0.538469310105682997_DP
    Dxi(19,1) = -0.538469310105682997_DP
    Dxi(19,2) =  0.538469310105682997_DP
    Dxi(20,1) = -0.538469310105682997_DP
    Dxi(20,2) = -0.538469310105682997_DP
    
    Dxi(21,1) =  0.906179845938664005_DP
    Dxi(21,2) =  0.0_DP
    Dxi(22,1) =  0.0_DP
    Dxi(22,2) =  0.906179845938664005_DP
    Dxi(23,1) = -0.906179845938664005_DP
    Dxi(23,2) =  0.0_DP
    Dxi(24,1) =  0.0_DP
    Dxi(24,2) = -0.906179845938664005_DP
    
    Dxi(25,1) =  0.0_DP
    Dxi(25,2) =  0.0_DP
    
    Domega(1) =  0.561343488624285944E-1_DP
    Domega(2) =  0.561343488624285944E-1_DP
    Domega(3) =  0.561343488624285944E-1_DP
    Domega(4) =  0.561343488624285944E-1_DP
    
    Domega(5) =  0.113399999999999834_DP
    Domega(6) =  0.113399999999999834_DP
    Domega(7) =  0.113399999999999834_DP
    Domega(8) =  0.113399999999999834_DP
    
    Domega(9) =  0.113399999999999834_DP
    Domega(10)=  0.113399999999999834_DP
    Domega(11)=  0.113399999999999834_DP
    Domega(12)=  0.113399999999999834_DP
    
    Domega(13)=  0.272286532550750485_DP
    Domega(14)=  0.272286532550750485_DP
    Domega(15)=  0.272286532550750485_DP
    Domega(16)=  0.272286532550750485_DP
    
    Domega(17)=  0.229085404223990666_DP
    Domega(18)=  0.229085404223990666_DP
    Domega(19)=  0.229085404223990666_DP
    Domega(20)=  0.229085404223990666_DP
    
    Domega(21)=  0.134785072387520868_DP
    Domega(22)=  0.134785072387520868_DP
    Domega(23)=  0.134785072387520868_DP
    Domega(24)=  0.134785072387520868_DP
    
    Domega(25)=  0.323634567901234682_DP
    
    ncubp     =  25
    
  case (CUB_PG1X1)
    Dxi(1,1)  =  0.5_DP
    Dxi(1,2)  =  0.5_DP
    Dxi(2,1)  =  0.5_DP
    Dxi(2,2)  = -0.5_DP
    Dxi(3,1)  = -0.5_DP
    Dxi(3,2)  =  0.5_DP
    Dxi(4,1)  = -0.5_DP
    Dxi(4,2)  = -0.5_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    Domega(3) =  1.0_DP
    Domega(4) =  1.0_DP
    
    ncubp     =  4

  case (CUB_PTRZ)
    ! vertices
    Dxi (1,1) = -1.0_DP
    Dxi (1,2) = -1.0_DP
    Dxi (2,1) =  1.0_DP
    Dxi (2,2) = -1.0_DP
    Dxi (3,1) =  1.0_DP
    Dxi (3,2) =  1.0_DP
    Dxi (4,1) = -1.0_DP
    Dxi (4,2) =  1.0_DP
    ! edge midpoints
    Dxi (5,1) =  0.0_DP
    Dxi (5,2) = -1.0_DP
    Dxi (6,1) =  1.0_DP
    Dxi (6,2) =  0.0_DP
    Dxi (7,1) =  0.0_DP
    Dxi (7,2) =  1.0_DP
    Dxi (8,1) = -1.0_DP
    Dxi (8,2) =  0.0_DP
    ! quad mipoint
    Dxi (9,1) =  0.0_DP
    Dxi (9,2) =  0.0_DP

    Domega(1) =  0.25_DP
    Domega(2) =  0.25_DP
    Domega(3) =  0.25_DP
    Domega(4) =  0.25_DP
    Domega(5) =  0.5_DP
    Domega(6) =  0.5_DP
    Domega(7) =  0.5_DP
    Domega(8) =  0.5_DP
    Domega(9) =  1.0_DP

    ncubp     =  9

  case (CUB_PG2X2)
    Dxi(1,1)  =  0.7886751345948130_DP
    Dxi(1,2)  =  0.7886751345948130_DP
    Dxi(2,1)  =  0.2113248654051870_DP
    Dxi(2,2)  =  0.7886751345948130_DP
    Dxi(3,1)  =  0.2113248654051870_DP
    Dxi(3,2)  =  0.2113248654051870_DP
    Dxi(4,1)  =  0.7886751345948130_DP
    Dxi(4,2)  =  0.2113248654051870_DP
    
    Dxi(5,1)  = -0.7886751345948130_DP
    Dxi(5,2)  =  0.7886751345948130_DP
    Dxi(6,1)  = -0.2113248654051870_DP
    Dxi(6,2)  =  0.7886751345948130_DP
    Dxi(7,1)  = -0.2113248654051870_DP
    Dxi(7,2)  =  0.2113248654051870_DP
    Dxi(8,1)  = -0.7886751345948130_DP
    Dxi(8,2)  =  0.2113248654051870_DP
    
    Dxi(9,1)  =  0.7886751345948130_DP
    Dxi(9,2)  = -0.7886751345948130_DP
    Dxi(10,1) =  0.2113248654051870_DP
    Dxi(10,2) = -0.7886751345948130_DP
    Dxi(11,1) =  0.2113248654051870_DP
    Dxi(11,2) = -0.2113248654051870_DP
    Dxi(12,1) =  0.7886751345948130_DP
    Dxi(12,2) = -0.2113248654051870_DP
    
    Dxi(13,1) = -0.7886751345948130_DP
    Dxi(13,2) = -0.7886751345948130_DP
    Dxi(14,1) = -0.2113248654051870_DP
    Dxi(14,2) = -0.7886751345948130_DP
    Dxi(15,1) = -0.2113248654051870_DP
    Dxi(15,2) = -0.2113248654051870_DP
    Dxi(16,1) = -0.7886751345948130_DP
    Dxi(16,2) = -0.2113248654051870_DP
    
    Domega(1) =  0.25_DP
    Domega(2) =  0.25_DP
    Domega(3) =  0.25_DP
    Domega(4) =  0.25_DP
    Domega(5) =  0.25_DP
    Domega(6) =  0.25_DP
    Domega(7) =  0.25_DP
    Domega(8) =  0.25_DP
    Domega(9) =  0.25_DP
    Domega(10)=  0.25_DP
    Domega(11)=  0.25_DP
    Domega(12)=  0.25_DP
    Domega(13)=  0.25_DP
    Domega(14)=  0.25_DP
    Domega(15)=  0.25_DP
    Domega(16)=  0.25_DP
    
    ncubp     =  16
    
  case (CUB_PG3X3)
    Dxi(1,1)  =  0.8872983346207415_DP
    Dxi(1,2)  =  0.8872983346207415_DP
    Dxi(2,1)  =  0.1127016653792585_DP
    Dxi(2,2)  =  0.8872983346207415_DP
    Dxi(3,1)  =  0.8872983346207415_DP
    Dxi(3,2)  =  0.1127016653792585_DP
    Dxi(4,1)  =  0.1127016653792585_DP
    Dxi(4,2)  =  0.1127016653792585_DP
    Dxi(5,1)  =  0.8872983346207415_DP
    Dxi(5,2)  =  0.5_DP
    Dxi(6,1)  =  0.1127016653792585_DP
    Dxi(6,2)  =  0.5_DP
    Dxi(7,1)  =  0.5_DP
    Dxi(7,2)  =  0.8872983346207415_DP
    Dxi(8,1)  =  0.5_DP
    Dxi(8,2)  =  0.1127016653792585_DP
    Dxi(9,1)  =  0.5_DP
    Dxi(9,2)  =  0.5_DP
    
    Dxi(10,1) = -0.8872983346207415_DP
    Dxi(10,2) =  0.8872983346207415_DP
    Dxi(11,1) = -0.1127016653792585_DP
    Dxi(11,2) =  0.8872983346207415_DP
    Dxi(12,1) = -0.8872983346207415_DP
    Dxi(12,2) =  0.1127016653792585_DP
    Dxi(13,1) = -0.1127016653792585_DP
    Dxi(13,2) =  0.1127016653792585_DP
    Dxi(14,1) = -0.8872983346207415_DP
    Dxi(14,2) =  0.5_DP
    Dxi(15,1) = -0.1127016653792585_DP
    Dxi(15,2) =  0.5_DP
    Dxi(16,1) = -0.5_DP
    Dxi(16,2) =  0.8872983346207415_DP
    Dxi(17,1) = -0.5_DP
    Dxi(17,2) =  0.1127016653792585_DP
    Dxi(18,1) = -0.5_DP
    Dxi(18,2) =  0.5_DP
    
    Dxi(19,1) =  0.8872983346207415_DP
    Dxi(19,2) = -0.8872983346207415_DP
    Dxi(20,1) =  0.1127016653792585_DP
    Dxi(20,2) = -0.8872983346207415_DP
    Dxi(21,1) =  0.8872983346207415_DP
    Dxi(21,2) = -0.1127016653792585_DP
    Dxi(22,1) =  0.1127016653792585_DP
    Dxi(22,2) = -0.1127016653792585_DP
    Dxi(23,1) =  0.8872983346207415_DP
    Dxi(23,2) = -0.5_DP
    Dxi(24,1) =  0.1127016653792585_DP
    Dxi(24,2) = -0.5_DP
    Dxi(25,1) =  0.5_DP
    Dxi(25,2) = -0.8872983346207415_DP
    Dxi(26,1) =  0.5_DP
    Dxi(26,2) = -0.1127016653792585_DP
    Dxi(27,1) =  0.5_DP
    Dxi(27,2) = -0.5_DP
    
    Dxi(28,1) = -0.8872983346207415_DP
    Dxi(28,2) = -0.8872983346207415_DP
    Dxi(29,1) = -0.1127016653792585_DP
    Dxi(29,2) = -0.8872983346207415_DP
    Dxi(30,1) = -0.8872983346207415_DP
    Dxi(30,2) = -0.1127016653792585_DP
    Dxi(31,1) = -0.1127016653792585_DP
    Dxi(31,2) = -0.1127016653792585_DP
    Dxi(32,1) = -0.8872983346207415_DP
    Dxi(32,2) = -0.5_DP
    Dxi(33,1) = -0.1127016653792585_DP
    Dxi(33,2) = -0.5_DP
    Dxi(34,1) = -0.5_DP
    Dxi(34,2) = -0.8872983346207415_DP
    Dxi(35,1) = -0.5_DP
    Dxi(35,2) = -0.1127016653792585_DP
    Dxi(36,1) = -0.5_DP
    Dxi(36,2) = -0.5_DP
    
    do i=0,27,9
      Domega(1+i) = 0.7716049382716050E-1_DP
      Domega(2+i) = 0.7716049382716050E-1_DP
      Domega(3+i) = 0.7716049382716050E-1_DP
      Domega(4+i) = 0.7716049382716050E-1_DP
      
      Domega(5+i) = 0.1234567901234570_DP
      Domega(6+i) = 0.1234567901234570_DP
      Domega(7+i) = 0.1234567901234570_DP
      Domega(8+i) = 0.1234567901234570_DP
    enddo
    
    Domega(9) =  0.1975308641975310_DP
    Domega(18)=  0.1975308641975310_DP
    Domega(27)=  0.1975308641975310_DP
    Domega(36)=  0.1975308641975310_DP
    
    ncubp     =  36
    
  case(CUB_SIMPSON)
    Dxi(1,1)  = -1.0_DP
    Dxi(1,2)  = -1.0_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(2,2)  = -1.0_DP
    Dxi(3,1)  =  1.0_DP
    Dxi(3,2)  = -1.0_DP
    Dxi(4,1)  = -1.0_DP
    Dxi(4,2)  =  0.0_DP
    Dxi(5,1)  =  0.0_DP
    Dxi(5,2)  =  0.0_DP
    Dxi(6,1)  =  1.0_DP
    Dxi(6,2)  =  0.0_DP
    Dxi(7,1)  = -1.0_DP
    Dxi(7,2)  =  1.0_DP
    Dxi(8,1)  =  0.0_DP
    Dxi(8,2)  =  1.0_DP
    Dxi(9,1)  =  1.0_DP
    Dxi(9,2)  =  1.0_DP

    Domega(1) = 0.111111111111111_DP
    Domega(2) = 0.444444444444444_DP
    Domega(3) = 0.111111111111111_DP
    Domega(4) = 0.444444444444444_DP
    Domega(5) = 1.777777777777778_DP
    Domega(6) = 0.444444444444444_DP
    Domega(7) = 0.111111111111111_DP
    Domega(8) = 0.444444444444444_DP
    Domega(9) = 0.111111111111111_DP
    ncubp     = 9

  case(CUB_3_8)
    Dxi(1,1)   = -1.0_DP
    Dxi(1,2)   = -1.0_DP
    Dxi(2,1)   = -0.33333333333033_DP
    Dxi(2,2)   = -1.0_DP
    Dxi(3,1)   =  0.33333333333333_DP
    Dxi(3,2)   = -1.0_DP
    Dxi(4,1)   =  1.0_DP
    Dxi(4,2)   = -1.0_DP
    Dxi(5,1)   = -1.0_DP
    Dxi(5,2)   = -0.33333333333033_DP
    Dxi(6,1)   = -0.33333333333033_DP
    Dxi(6,2)   = -0.33333333333033_DP
    Dxi(7,1)   =  0.33333333333333_DP
    Dxi(7,2)   = -0.33333333333033_DP
    Dxi(8,1)   =  1.0_DP
    Dxi(8,2)   = -0.33333333333033_DP
    Dxi(9,1)   = -1.0_DP
    Dxi(9,2)   =  0.33333333333033_DP
    Dxi(10,1)  = -0.33333333333033_DP
    Dxi(10,2)  =  0.33333333333333_DP
    Dxi(11,1)  =  0.33333333333033_DP
    Dxi(11,2)  =  0.33333333333333_DP
    Dxi(12,1)  =  1.0_DP
    Dxi(12,2)  =  0.33333333333333_DP
    Dxi(13,1)  = -1.0_DP
    Dxi(13,2)  =  1.0_DP
    Dxi(14,1)  = -0.33333333333333_DP
    Dxi(14,2)  =  1.0_DP
    Dxi(15,1)  =  0.33333333333333_DP
    Dxi(15,2)  =  1.0_DP
    Dxi(16,1)  =  1.0_DP
    Dxi(16,2)  =  1.0_DP


    Domega(1)  = 0.0625_DP
    Domega(2)  = 0.1875_DP
    Domega(3)  = 0.1875_DP
    Domega(4)  = 0.0625_DP
    Domega(5)  = 0.1875_DP
    Domega(6)  = 0.5625_DP
    Domega(7)  = 0.5625_DP
    Domega(8)  = 0.1875_DP
    Domega(9)  = 0.1875_DP
    Domega(10) = 0.5625_DP
    Domega(11) = 0.5625_DP
    Domega(12) = 0.1875_DP
    Domega(13) = 0.0625_DP
    Domega(14) = 0.1875_DP
    Domega(15) = 0.1875_DP
    Domega(16) = 0.0625_DP
    ncubp     = 16
    
  !triangle cubature formulas
  case(CUB_G1_T)
    Dxi(1,1)  =  0.3333333333333333_DP
    Dxi(1,2)  =  0.3333333333333333_DP
    Dxi(1,3)  =  0.3333333333333333_DP
    
    Domega(1) =  0.5_DP
    
    ncubp     =  1
    
  case(CUB_TRZ_T)
    Dxi(1,1)  =  1.0_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(1,3)  =  0.0_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(2,2)  =  1.0_DP
    Dxi(2,3)  =  0.0_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  =  0.0_DP
    Dxi(3,3)  =  1.0_DP
    
    Domega(1) =  0.1666666666666667_DP
    Domega(2) =  0.1666666666666667_DP
    Domega(3) =  0.1666666666666667_DP
    
    ncubp     =  3
    
  case(CUB_G3MP_T)
    Dxi(1,1)  =  0.5_DP
    Dxi(1,2)  =  0.5_DP
    Dxi(1,3)  =  0.0_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(2,2)  =  0.5_DP
    Dxi(2,3)  =  0.5_DP
    Dxi(3,1)  =  0.5_DP
    Dxi(3,2)  =  0.0_DP
    Dxi(3,3)  =  0.5_DP
    
    Domega(1) =  0.1666666666666667_DP
    Domega(2) =  0.1666666666666667_DP
    Domega(3) =  0.1666666666666667_DP
    
    ncubp     =  3
    
  case(CUB_G3_T)
    Dxi(1,1)  =  0.6666666666666667_DP
    Dxi(1,2)  =  0.1666666666666667_DP
    Dxi(1,3)  =  0.1666666666666667_DP
    Dxi(2,1)  =  0.1666666666666667_DP
    Dxi(2,2)  =  0.6666666666666667_DP
    Dxi(2,3)  =  0.1666666666666667_DP
    Dxi(3,1)  =  0.1666666666666667_DP
    Dxi(3,2)  =  0.1666666666666667_DP
    Dxi(3,3)  =  0.6666666666666667_DP
    
    Domega(1) =  0.1666666666666667_DP
    Domega(2) =  0.1666666666666667_DP
    Domega(3) =  0.1666666666666667_DP
    
    ncubp     =  3
    
  case(CUB_VMC)
    !center
    Dxi(1,1)  =  0.3333333333333333_DP
    Dxi(1,2)  =  0.3333333333333333_DP
    Dxi(1,3)  =  0.3333333333333333_DP
    !midpoints
    Dxi(2,1)  =  0.5_DP
    Dxi(2,2)  =  0.5_DP
    Dxi(2,3)  =  0.0_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  =  0.5_DP
    Dxi(3,3)  =  0.5_DP
    Dxi(4,1)  =  0.5_DP
    Dxi(4,2)  =  0.0_DP
    Dxi(4,3)  =  0.5_DP
    !vertices
    Dxi(5,1)  =  1.0_DP
    Dxi(5,2)  =  0.0_DP
    Dxi(5,3)  =  0.0_DP
    Dxi(6,1)  =  0.0_DP
    Dxi(6,2)  =  1.0_DP
    Dxi(6,3)  =  0.0_DP
    Dxi(7,1)  =  0.0_DP
    Dxi(7,2)  =  0.0_DP
    Dxi(7,3)  =  1.0_DP

    Domega(1) =  0.225_DP
    Domega(2) =  0.0666666666666667_DP
    Domega(3) =  0.0666666666666667_DP
    Domega(4) =  0.0666666666666667_DP
    Domega(5) =  0.025_DP
    Domega(6) =  0.025_DP
    Domega(7) =  0.025_DP

    ncubp     =  7

  ! 3D cubature formulas

  case(CUB_G1_3D)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(1,3)  =  0.0_DP
    
    Domega(1) =  8.0_DP
    
    ncubp      =  1
    
  case(CUB_MIDAREA_3D)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(1,3)  = -1.0_DP
    
    Dxi(2,1)  =  0.0_DP
    Dxi(2,2)  =  -1.0_DP
    Dxi(2,3)  =  0.0_DP

    Dxi(3,1)  =  1.0_DP
    Dxi(3,2)  =  0.0_DP
    Dxi(3,3)  =  0.0_DP

    Dxi(4,1)  =  0.0_DP
    Dxi(4,2)  =  1.0_DP
    Dxi(4,3)  =  0.0_DP

    Dxi(5,1)  =  -1.0_DP
    Dxi(5,2)  =  0.0_DP
    Dxi(5,3)  =  0.0_DP

    Dxi(6,1)  =  0.0_DP
    Dxi(6,2)  =  0.0_DP
    Dxi(6,3)  =  1.0_DP

    Domega(1) =  1.333333333333333_DP
    Domega(2) =  1.333333333333333_DP
    Domega(3) =  1.333333333333333_DP
    Domega(4) =  1.333333333333333_DP
    Domega(5) =  1.333333333333333_DP
    Domega(6) =  1.333333333333333_DP

    ncubp     =  6

  case(CUB_TRZ_3D)
    Dxi(1,1)  = -1.0_DP
    Dxi(1,2)  = -1.0_DP
    Dxi(1,3)  = -1.0_DP

    Dxi(2,1)  =  1.0_DP
    Dxi(2,2)  = -1.0_DP
    Dxi(2,3)  = -1.0_DP

    Dxi(3,1)  =  1.0_DP
    Dxi(3,2)  =  1.0_DP
    Dxi(3,3)  = -1.0_DP

    Dxi(4,1)  = -1.0_DP
    Dxi(4,2)  =  1.0_DP
    Dxi(4,3)  = -1.0_DP

    Dxi(5,1)  = -1.0_DP
    Dxi(5,2)  = -1.0_DP
    Dxi(5,3)  =  1.0_DP

    Dxi(6,1)  =  1.0_DP
    Dxi(6,2)  = -1.0_DP
    Dxi(6,3)  =  1.0_DP

    Dxi(7,1)  =  1.0_DP
    Dxi(7,2)  =  1.0_DP
    Dxi(7,3)  =  1.0_DP

    Dxi(8,1)  = -1.0_DP
    Dxi(8,2)  =  1.0_DP
    Dxi(8,3)  =  1.0_DP

    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    Domega(3) =  1.0_DP
    Domega(4) =  1.0_DP
    Domega(5) =  1.0_DP
    Domega(6) =  1.0_DP
    Domega(7) =  1.0_DP
    Domega(8) =  1.0_DP

    ncubp     =  8
    
  case(CUB_G2_3D)

    Dxi(1,1)  =  0.577350269189626_DP
    Dxi(1,2)  =  0.577350269189626_DP
    Dxi(1,3)  = -0.577350269189626_DP

    Dxi(2,1)  = -0.577350269189626_DP
    Dxi(2,2)  =  0.577350269189626_DP
    Dxi(2,3)  = -0.577350269189626_DP

    Dxi(3,1)  = -0.577350269189626_DP
    Dxi(3,2)  = -0.577350269189626_DP
    Dxi(3,3)  = -0.577350269189626_DP

    Dxi(4,1)  =  0.577350269189626_DP
    Dxi(4,2)  = -0.577350269189626_DP
    Dxi(4,3)  = -0.577350269189626_DP

    Dxi(5,1)  =  0.577350269189626_DP
    Dxi(5,2)  =  0.577350269189626_DP
    Dxi(5,3)  =  0.577350269189626_DP

    Dxi(6,1)  = -0.577350269189626_DP
    Dxi(6,2)  =  0.577350269189626_DP
    Dxi(6,3)  =  0.577350269189626_DP

    Dxi(7,1)  = -0.577350269189626_DP
    Dxi(7,2)  = -0.577350269189626_DP
    Dxi(7,3)  =  0.577350269189626_DP

    Dxi(8,1)  =  0.577350269189626_DP
    Dxi(8,2)  = -0.577350269189626_DP
    Dxi(8,3)  =  0.577350269189626_DP

    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    Domega(3) =  1.0_DP
    Domega(4) =  1.0_DP
    Domega(5) =  1.0_DP
    Domega(6) =  1.0_DP
    Domega(7) =  1.0_DP
    Domega(8) =  1.0_DP

    ncubp     =  8

  case(CUB_G3_3D)
    
    Dxi(1,1)  =  0.774596669241483_DP
    Dxi(1,2)  =  0.774596669241483_DP
    Dxi(1,3)  = -0.774596669241483_DP
    
    Dxi(2,1)  = -0.774596669241483_DP
    Dxi(2,2)  =  0.774596669241483_DP
    Dxi(2,3)  = -0.774596669241483_DP
    
    Dxi(3,1)  =  0.774596669241483_DP
    Dxi(3,2)  = -0.774596669241483_DP
    Dxi(3,3)  = -0.774596669241483_DP
    
    Dxi(4,1)  = -0.774596669241483_DP
    Dxi(4,2)  = -0.774596669241483_DP
    Dxi(4,3)  = -0.774596669241483_DP
    
    Dxi(5,1)  =  0.774596669241483_DP
    Dxi(5,2)  =  0.0_DP
    Dxi(5,3)  = -0.774596669241483_DP
    
    Dxi(6,1)  = -0.774596669241483_DP
    Dxi(6,2)  =  0.0_DP
    Dxi(6,3)  = -0.774596669241483_DP
    
    Dxi(7,1)  =  0.0_DP
    Dxi(7,2)  =  0.774596669241483_DP
    Dxi(7,3)  = -0.774596669241483_DP
    
    Dxi(8,1)  =  0.0_DP
    Dxi(8,2)  = -0.774596669241483_DP
    Dxi(8,3)  = -0.774596669241483_DP
    
    Dxi(9,1)  =  0.0_DP
    Dxi(9,2)  =  0.0_DP
    Dxi(9,3)  = -0.774596669241483_DP
    
    Dxi(10,1) =  0.774596669241483_DP
    Dxi(10,2) =  0.774596669241483_DP
    Dxi(10,3) =  0.0_DP
    
    
    Dxi(11,1) = -0.774596669241483_DP
    Dxi(11,2) =  0.774596669241483_DP
    Dxi(11,3) =  0.0_DP
    
    Dxi(12,1) =  0.774596669241483_DP
    Dxi(12,2) = -0.774596669241483_DP
    Dxi(12,3) =  0.0_DP
    
    Dxi(13,1) = -0.774596669241483_DP
    Dxi(13,2) = -0.774596669241483_DP
    Dxi(13,3) =  0.0_DP
    
    Dxi(14,1) =  0.774596669241483_DP
    Dxi(14,2) =  0.0_DP
    Dxi(14,3) =  0.0_DP
    
    Dxi(15,1) = -0.774596669241483_DP
    Dxi(15,2) =  0.0_DP
    Dxi(15,3) =  0.0_DP
    
    Dxi(16,1) =  0.0_DP
    Dxi(16,2) =  0.774596669241483_DP
    Dxi(16,3) =  0.0_DP
    
    Dxi(17,1) =  0.0_DP
    Dxi(17,2) = -0.774596669241483_DP
    Dxi(17,3) =  0.0_DP
    
    Dxi(18,1) =  0.0_DP
    Dxi(18,2) =  0.0_DP
    Dxi(18,3) =  0.0_DP
    
    Dxi(19,1) =  0.774596669241483_DP
    Dxi(19,2) =  0.774596669241483_DP
    Dxi(19,3) =  0.774596669241483_DP
    
    Dxi(20,1) = -0.774596669241483_DP
    Dxi(20,2) =  0.774596669241483_DP
    Dxi(20,3) =  0.774596669241483_DP
    
    Dxi(21,1) =  0.774596669241483_DP
    Dxi(21,2) = -0.774596669241483_DP
    Dxi(21,3) =  0.774596669241483_DP
    
    Dxi(22,1) = -0.774596669241483_DP
    Dxi(22,2) = -0.774596669241483_DP
    Dxi(22,3) =  0.774596669241483_DP
    
    Dxi(23,1) =  0.774596669241483_DP
    Dxi(23,2) =  0.0_DP
    Dxi(23,3) =  0.774596669241483_DP
    
    Dxi(24,1) = -0.774596669241483_DP
    Dxi(24,2) =  0.0_DP
    Dxi(24,3) =  0.774596669241483_DP
    
    Dxi(25,1) =  0.0_DP
    Dxi(25,2) =  0.774596669241483_DP
    Dxi(25,3) =  0.774596669241483_DP
    
    Dxi(26,1) =  0.0_DP
    Dxi(26,2) = -0.774596669241483_DP
    Dxi(26,3) =  0.774596669241483_DP
    
    Dxi(27,1) =  0.0_DP
    Dxi(27,2) =  0.0_DP
    Dxi(27,3) =  0.774596669241483_DP

    Domega(1) =  0.17146776406036_DP
    Domega(2) =  0.17146776406036_DP
    Domega(3) =  0.17146776406036_DP
    Domega(4) =  0.17146776406036_DP
    Domega(5) =  0.27434842249657_DP
    Domega(6) =  0.27434842249657_DP
    Domega(7) =  0.27434842249657_DP
    Domega(8) =  0.27434842249657_DP
    Domega(9) =  0.43895747599451_DP
    Domega(10)=  0.27434842249657_DP
    Domega(11)=  0.27434842249657_DP
    Domega(12)=  0.27434842249657_DP
    Domega(13)=  0.27434842249657_DP
    Domega(14)=  0.43895747599451_DP
    Domega(15)=  0.43895747599451_DP
    Domega(16)=  0.43895747599451_DP
    Domega(17)=  0.43895747599451_DP
    Domega(18)=  0.70233196159122_DP
    Domega(19)=  0.17146776406036_DP
    Domega(20)=  0.17146776406036_DP
    Domega(21)=  0.17146776406036_DP
    Domega(22)=  0.17146776406036_DP
    Domega(23)=  0.27434842249657_DP
    Domega(24)=  0.27434842249657_DP
    Domega(25)=  0.27434842249657_DP
    Domega(26)=  0.27434842249657_DP
    Domega(27)=  0.43895747599451_DP

    ncubp     =  27

  ! tetrahedra
  case(CUB_G1_3D_T)
    Dxi(1,1)  =  0.25_DP
    Dxi(1,2)  =  0.25_DP
    Dxi(1,3)  =  0.25_DP
    Dxi(1,4)  =  0.25_DP
    
    Domega(1) = 0.1666666666666667_DP
    
    ncubp = 1
    
  case(CUB_TRZ_3D_T)
    Dxi(1,1)  =  1.0_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(1,3)  =  0.0_DP
    Dxi(1,4)  =  0.0_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(2,2)  =  1.0_DP
    Dxi(2,3)  =  0.0_DP
    Dxi(2,4)  =  0.0_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  =  0.0_DP
    Dxi(3,3)  =  1.0_DP
    Dxi(3,4)  =  0.0_DP
    Dxi(4,1)  =  0.0_DP
    Dxi(4,2)  =  0.0_DP
    Dxi(4,3)  =  0.0_DP
    Dxi(4,4)  =  1.0_DP

    Domega(1) =  0.0416666666666667_DP
    Domega(2) =  0.0416666666666667_DP
    Domega(3) =  0.0416666666666667_DP
    Domega(4) =  0.0416666666666667_DP
    
    ncubp     =  4

  case(CUB_S2_3D_T)
    Dxi(1,1)  = 0.5854101966249685_DP
    Dxi(1,2)  = 0.1381966011250105_DP
    Dxi(1,3)  = 0.1381966011250105_DP
    Dxi(1,4)  = 0.1381966011250105_DP
    Dxi(2,1)  = 0.1381966011250105_DP
    Dxi(2,2)  = 0.5854101966249685_DP
    Dxi(2,3)  = 0.1381966011250105_DP
    Dxi(2,4)  = 0.1381966011250105_DP
    Dxi(3,1)  = 0.1381966011250105_DP
    Dxi(3,2)  = 0.1381966011250105_DP
    Dxi(3,3)  = 0.5854101966249685_DP
    Dxi(3,4)  = 0.1381966011250105_DP
    Dxi(4,1)  = 0.1381966011250105_DP
    Dxi(4,2)  = 0.1381966011250105_DP
    Dxi(4,3)  = 0.1381966011250105_DP
    Dxi(4,4)  = 0.5854101966249685_DP
    
    Domega(1) = 0.0416666666666667_DP
    Domega(2) = 0.0416666666666667_DP
    Domega(3) = 0.0416666666666667_DP
    Domega(4) = 0.0416666666666667_DP
    
    ncubp = 4
  
  case(CUB_S3_3D_T)
    Dxi( 1,1) = 0.1438564719343849_DP
    Dxi( 1,2) = 0.1438564719343849_DP
    Dxi( 1,3) = 0.1438564719343849_DP
    Dxi( 1,4) = 0.5684305841968454_DP
    Dxi( 2,1) = 0.1438564719343849_DP
    Dxi( 2,2) = 0.1438564719343849_DP
    Dxi( 2,3) = 0.5684305841968454_DP
    Dxi( 2,4) = 0.1438564719343849_DP
    Dxi( 3,1) = 0.1438564719343849_DP
    Dxi( 3,2) = 0.5684305841968454_DP
    Dxi( 3,3) = 0.1438564719343849_DP
    Dxi( 3,4) = 0.1438564719343849_DP
    Dxi( 4,1) = 0.5684305841968454_DP
    Dxi( 4,2) = 0.1438564719343849_DP
    Dxi( 4,3) = 0.1438564719343849_DP
    Dxi( 4,4) = 0.1438564719343849_DP
    Dxi( 5,1) = 0.5_DP
    Dxi( 5,2) = 0.5_DP
    Dxi( 5,3) = 0.0_DP
    Dxi( 5,4) = 0.0_DP
    Dxi( 6,1) = 0.5_DP
    Dxi( 6,2) = 0.0_DP
    Dxi( 6,3) = 0.5_DP
    Dxi( 6,4) = 0.0_DP
    Dxi( 7,1) = 0.5_DP
    Dxi( 7,2) = 0.0_DP
    Dxi( 7,3) = 0.0_DP
    Dxi( 7,4) = 0.5_DP
    Dxi( 8,1) = 0.0_DP
    Dxi( 8,2) = 0.5_DP
    Dxi( 8,3) = 0.5_DP
    Dxi( 8,4) = 0.0_DP
    Dxi( 9,1) = 0.0_DP
    Dxi( 9,2) = 0.5_DP
    Dxi( 9,3) = 0.0_DP
    Dxi( 9,4) = 0.5_DP
    Dxi(10,1) = 0.0_DP
    Dxi(10,2) = 0.0_DP
    Dxi(10,3) = 0.5_DP
    Dxi(10,4) = 0.5_DP

    Domega( 1) = 0.0362941783134010_DP
    Domega( 2) = 0.0362941783134010_DP
    Domega( 3) = 0.0362941783134010_DP
    Domega( 4) = 0.0362941783134010_DP
    Domega( 5) = 0.0035816589021771_DP
    Domega( 6) = 0.0035816589021771_DP
    Domega( 7) = 0.0035816589021771_DP
    Domega( 8) = 0.0035816589021771_DP
    Domega( 9) = 0.0035816589021771_DP
    Domega(10) = 0.0035816589021771_DP
    
    ncubp = 10
  
  case(CUB_S5_3D_T)
    Dxi( 1,1)  = 0.25_DP
    Dxi( 1,2)  = 0.25_DP
    Dxi( 1,3)  = 0.25_DP
    Dxi( 1,4)  = 0.25_DP
    Dxi( 2,1)  = 0.7240867658418309_DP
    Dxi( 2,2)  = 0.0919710780527230_DP
    Dxi( 2,3)  = 0.0919710780527230_DP
    Dxi( 2,4)  = 0.0919710780527230_DP
    Dxi( 3,1)  = 0.0919710780527230_DP
    Dxi( 3,2)  = 0.7240867658418309_DP
    Dxi( 3,3)  = 0.0919710780527230_DP
    Dxi( 3,4)  = 0.0919710780527230_DP
    Dxi( 4,1)  = 0.0919710780527230_DP
    Dxi( 4,2)  = 0.0919710780527230_DP
    Dxi( 4,3)  = 0.7240867658418309_DP
    Dxi( 4,4)  = 0.0919710780527230_DP
    Dxi( 5,1)  = 0.0919710780527230_DP
    Dxi( 5,2)  = 0.0919710780527230_DP
    Dxi( 5,3)  = 0.0919710780527230_DP
    Dxi( 5,4)  = 0.7240867658418309_DP
    Dxi( 6,1)  = 0.0406191165111103_DP
    Dxi( 6,2)  = 0.3197936278296299_DP
    Dxi( 6,3)  = 0.3197936278296299_DP
    Dxi( 6,4)  = 0.3197936278296299_DP
    Dxi( 7,1)  = 0.3197936278296299_DP
    Dxi( 7,2)  = 0.0406191165111103_DP
    Dxi( 7,3)  = 0.3197936278296299_DP
    Dxi( 7,4)  = 0.3197936278296299_DP
    Dxi( 8,1)  = 0.3197936278296299_DP
    Dxi( 8,2)  = 0.3197936278296299_DP
    Dxi( 8,3)  = 0.0406191165111103_DP
    Dxi( 8,4)  = 0.3197936278296299_DP
    Dxi( 9,1)  = 0.3197936278296299_DP
    Dxi( 9,2)  = 0.3197936278296299_DP
    Dxi( 9,3)  = 0.3197936278296299_DP
    Dxi( 9,4)  = 0.0406191165111103_DP
    Dxi(10,1)  = 0.4436491673103708_DP
    Dxi(10,2)  = 0.0563508326896292_DP
    Dxi(10,3)  = 0.0563508326896292_DP
    Dxi(10,4)  = 0.4436491673103708_DP
    Dxi(11,1)  = 0.4436491673103708_DP
    Dxi(11,2)  = 0.0563508326896292_DP
    Dxi(11,3)  = 0.4436491673103708_DP
    Dxi(11,4)  = 0.0563508326896292_DP
    Dxi(12,1)  = 0.4436491673103708_DP
    Dxi(12,2)  = 0.4436491673103708_DP
    Dxi(12,3)  = 0.0563508326896292_DP
    Dxi(12,4)  = 0.0563508326896292_DP
    Dxi(13,1)  = 0.0563508326896292_DP
    Dxi(13,2)  = 0.0563508326896292_DP
    Dxi(13,3)  = 0.4436491673103708_DP
    Dxi(13,4)  = 0.4436491673103708_DP
    Dxi(14,1)  = 0.0563508326896292_DP
    Dxi(14,2)  = 0.4436491673103708_DP
    Dxi(14,3)  = 0.0563508326896292_DP
    Dxi(14,4)  = 0.4436491673103708_DP
    Dxi(15,1)  = 0.0563508326896292_DP
    Dxi(15,2)  = 0.4436491673103708_DP
    Dxi(15,3)  = 0.4436491673103708_DP
    Dxi(15,4)  = 0.0563508326896292_DP
    
    Domega( 1) = 0.0197530864197531_DP
    Domega( 2) = 0.0119895139631698_DP
    Domega( 3) = 0.0119895139631698_DP
    Domega( 4) = 0.0119895139631698_DP
    Domega( 5) = 0.0119895139631698_DP
    Domega( 6) = 0.0115113678710454_DP
    Domega( 7) = 0.0115113678710454_DP
    Domega( 8) = 0.0115113678710454_DP
    Domega( 9) = 0.0115113678710454_DP
    Domega(10) = 0.0088183421516755_DP
    Domega(11) = 0.0088183421516755_DP
    Domega(12) = 0.0088183421516755_DP
    Domega(13) = 0.0088183421516755_DP
    Domega(14) = 0.0088183421516755_DP
    Domega(15) = 0.0088183421516755_DP
    
    ncubp = 15

  ! pyramid
  case(CUB_G1_3D_Y)
    Dxi(1,1) = 0.0_DP
    Dxi(1,2) = 0.0_DP
    Dxi(1,3) = 0.5_DP
    
    Domega(1) = 1.3333333333333333_DP
    
    ncubp = 1
  
  case(CUB_TRZ_3D_Y)
    Dxi(1,1) = -1.0_DP
    Dxi(1,2) = -1.0_DP
    Dxi(1,3) =  0.0_DP
    Dxi(2,1) =  1.0_DP
    Dxi(2,2) = -1.0_DP
    Dxi(2,3) =  0.0_DP
    Dxi(3,1) =  1.0_DP
    Dxi(3,2) =  1.0_DP
    Dxi(3,3) =  0.0_DP
    Dxi(4,1) = -1.0_DP
    Dxi(4,2) =  1.0_DP
    Dxi(4,3) =  0.0_DP
    Dxi(5,1) =  0.0_DP
    Dxi(5,2) =  0.0_DP
    Dxi(5,3) =  1.0_DP
    
    Domega(1) = 0.26666666666666666_DP
    Domega(2) = 0.26666666666666666_DP
    Domega(3) = 0.26666666666666666_DP
    Domega(4) = 0.26666666666666666_DP
    Domega(5) = 0.26666666666666666_DP
    
    ncubp = 5
  
  ! prism
  case(CUB_G1_3D_R)
    Dxi(1,1) = 0.3333333333333333_DP
    Dxi(1,2) = 0.3333333333333333_DP
    Dxi(1,3) = 0.0_DP
    
    Domega(1) = 1.0_DP
    
    ncubp = 1
  
  case(CUB_TRZ_3D_R)
    Dxi(1,1) =  0.0_DP
    Dxi(1,2) =  0.0_DP
    Dxi(1,3) = -1.0_DP
    Dxi(2,1) =  1.0_DP
    Dxi(2,2) =  0.0_DP
    Dxi(2,3) = -1.0_DP
    Dxi(3,1) =  0.0_DP
    Dxi(3,2) =  1.0_DP
    Dxi(3,3) = -1.0_DP
    Dxi(4,1) =  0.0_DP
    Dxi(4,2) =  0.0_DP
    Dxi(4,3) =  1.0_DP
    Dxi(5,1) =  1.0_DP
    Dxi(5,2) =  0.0_DP
    Dxi(5,3) =  1.0_DP
    Dxi(6,1) =  0.0_DP
    Dxi(6,2) =  1.0_DP
    Dxi(6,3) =  1.0_DP
    
    Domega(1) = 0.16666666666666666_DP
    Domega(2) = 0.16666666666666666_DP
    Domega(3) = 0.16666666666666666_DP
    Domega(4) = 0.16666666666666666_DP
    Domega(5) = 0.16666666666666666_DP
    Domega(6) = 0.16666666666666666_DP
    
    ncubp = 6
  
  case(CUB_G2_3D_R)
    Dxi(1,1)  =  0.5_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(1,3)  = -0.577350269189626_DP
    Dxi(2,1)  =  0.5_DP
    Dxi(2,2)  =  0.5_DP
    Dxi(2,3)  = -0.577350269189626_DP
    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  =  0.5_DP
    Dxi(3,3)  = -0.577350269189626_DP
    Dxi(4,1)  =  0.5_DP
    Dxi(4,2)  =  0.0_DP
    Dxi(4,3)  =  0.577350269189626_DP
    Dxi(5,1)  =  0.5_DP
    Dxi(5,2)  =  0.5_DP
    Dxi(5,3)  =  0.577350269189626_DP
    Dxi(6,1)  =  0.0_DP
    Dxi(6,2)  =  0.5_DP
    Dxi(6,3)  =  0.577350269189626_DP
    
    Domega(1) = 0.16666666666666666_DP
    Domega(2) = 0.16666666666666666_DP
    Domega(3) = 0.16666666666666666_DP
    Domega(4) = 0.16666666666666666_DP
    Domega(5) = 0.16666666666666666_DP
    Domega(6) = 0.16666666666666666_DP
    
    ncubp = 6

  case default
    call output_line('Unknown cubature type!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'cub_getCubPoints')
    call sys_halt()
  end select
   
  end subroutine cub_getCubPoints

end module
