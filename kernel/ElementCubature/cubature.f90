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
!# By calling cub_getCubPoints with one of the ID's, one gets the coordinates
!# and weights of that specific cubature formula on the reference element
!# (either $[-1,1]^2$ or the triangle $(0,0), (0,1), (1,0)$ (depending
!# on whether it's a cubature formulas for triangles or quadrilaterals.
!#
!# Furthermore, the routine cub_igetID allows to translate a string into
!# a cubature formula ID.
!#
!# The following routines can be found in this module:
!#
!# 1.) cub_igetID
!#     -> Interpret a string as cubature type identifier. This is typically
!#        used for parsing INI files.
!#
!# 2.) cub_igetNumPts
!#     -> Returns the number of points for a given cubature type identifier.
!#
!# 3.) cub_igetShape
!#     -> Returns the element shape identifier on which a cubature formula
!#        is defined.
!#
!# 4.) cub_igetCoordDim
!#     -> Returns the dimension of the coordinates for a given cubature
!#        type identifier.
!#
!# 5.) cub_getCubature
!#     -> Returns the cubature points and weights for a given cubature type
!#        identifier.
!#
!# 6.) cub_getCubPoints
!#     => DEPRECATED: See note below.
!#     -> Get information about cubature points for a specific cubature
!#        formula: Coordinates, number and weights.
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
!#     ! Don't forget to deallocate Dpoints and Domega after you used them!
!#
!# </code>
!#
!# </purpose>
!##############################################################################

module cubature

  use fsystem
  use basicgeometry

  implicit none
  
!<constants>

!<constantblock variable="ccubType" description="1D formulas">  

  ! 1-point Gauss formula, 1D, degree = 2, ncubp = 1
  integer, parameter :: CUB_G1_1D = 101
 
  ! trapezoidal rule, 1D, degree = 2, ncubp = 2
  integer, parameter :: CUB_TRZ_1D = 102

  ! 2-point Gauss formula, 1D, degree = 4, ncubp = 2
  integer, parameter :: CUB_G2_1D = 103

  ! 3-point Gauss formula, 1D, degree = 6, ncubp = 3
  integer, parameter :: CUB_G3_1D = 104

  ! 4-point Gauss formula, 1D, degree = 8, ncubp = 4
  integer, parameter :: CUB_G4_1D = 105

  ! 5-point Gauss formula, 1D, degree = 10, ncubp = 5
  integer, parameter :: CUB_G5_1D = 106

  ! Simpson-rule, 1D, degree = 4, ncubp = 3
  integer, parameter :: CUB_SIMPSON_1D = 107

  ! 6-point Gauss formula, 1D, degree = 12, ncubp = 6
  integer, parameter :: CUB_G6_1D = 108

!</constantblock>

!<constantblock variable="ccubType" description="2D formulas, quad">  

  ! 1x1 Gauss formula, degree = 2, ncubp = 1
  integer, parameter :: CUB_G1X1 = 201
  integer, parameter :: CUB_G1_2D = CUB_G1X1

  ! trapezoidal rule, degree = 2, ncubp = 4
  integer, parameter :: CUB_TRZ = 202

  ! midpoint rule, degree = 2, ncubp = 4
  integer, parameter :: CUB_MID = 203

  ! 2x2 Gauss formula, degree = 4, ncubp = 4
  integer, parameter :: CUB_G2X2 = 204
  integer, parameter :: CUB_G2_2D = CUB_G2X2

  ! Newton formula 1, degree = 4, ncubp = 4 
  integer, parameter :: CUB_NS1 = 205

  ! Newton formula 2, degree = 5, ncubp = 6
  integer, parameter :: CUB_NS2 = 206

  ! Newton formula 3, degree = 6, ncubp = 7
  integer, parameter :: CUB_NS3 = 207

  ! 3x3 Gauss formula, degree = 6, ncubp = 9
  integer, parameter :: CUB_G3X3 = 208
  integer, parameter :: CUB_G3_2D = CUB_G3X3

  ! Gauss formula, degree = 7, ncubp = 12
  integer, parameter :: CUB_G = 209

  ! 4x4 Gauss formula, degree = 8, ncubp = 16
  integer, parameter :: CUB_G4X4 = 210
  integer, parameter :: CUB_G4_2D = CUB_G4X4

  ! 5x5 Gauss formula, degree = 10, ncubp = 25
  integer, parameter :: CUB_G5X5 = 211
  integer, parameter :: CUB_G5_2D = CUB_G5X5

  ! piecewise 1x1 Gauss formula, degree = 2, ncubp = 4
  integer, parameter :: CUB_PG1X1 = 212

  ! piecewise trapezoidal rule, degree = 2, ncubp = 9
  integer, parameter :: CUB_PTRZ = 213

  ! piecewise 2x2 Gauss formula, degree = 4, ncubp = 16 
  integer, parameter :: CUB_PG2X2 = 214

  ! piecewise 3x3 Gauss formula, degree = 6, ncubp = 36
  integer, parameter :: CUB_PG3X3 = 215
  
  ! Simpson rule (corners and element midpoint), degree = 3, ncubp = 9
  integer, parameter :: CUB_SIMPSON = 216

  ! Simpson 3/8 rule (corners and 1/3 + 2/3), degree = 3, ncubp = 16
  integer, parameter :: CUB_3_8 = 217

  ! 6x6 Gauss formula, degree = 12, ncubp = 36
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer, parameter :: CUB_G6_2D = 218

!</constantblock>

!<constantblock variable="ccubType" description="2D formulas, tri">

  ! 1-point Gauss formula, triangle, degree = 2, ncubp = 1
  integer, parameter :: CUB_G1_T = 250

  ! trapezoidal rule, triangle, degree = 2, ncubp = 3 
  integer, parameter :: CUB_TRZ_T = 251

  ! 3-point Gauss formula, triangle, degree = 3, ncubp = 3
  integer, parameter :: CUB_G3_T = 252

  ! Collatz formula, degree = 3, ncubp = 3
  integer, parameter :: CUB_Collatz = 253

  ! vertices, midpoints, center, degree = 4, ncubp = 7 
  integer, parameter :: CUB_VMC = 254
!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, hexa">

  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer, parameter :: CUB_G1_3D = 301

  ! midpoints of areas, 3D, degree = 2, ncubp = 6
  integer, parameter :: CUB_MIDAREA_3D = 302

  ! trapezoidal rule, 3D, degree = 2, ncubp = 8
  integer, parameter :: CUB_TRZ_3D = 303

  ! 2-point Gauss formula, 3D, degree = 4, ncubp = 8
  integer, parameter :: CUB_G2_3D = 304
  
  ! 3-point Gauss formula, 3D, degree = 6, ncubp = 27
  integer, parameter :: CUB_G3_3D = 305
  
  ! 4-point Gauss formula, 3D, degree = 8, ncubp = 64
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer, parameter :: CUB_G4_3D = 306
  
  ! 5-point Gauss formula, 3D, degree = 10, ncubp = 125
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer, parameter :: CUB_G5_3D = 307

  ! 6-point Gauss formula, degree = 12, ncubp = 216
  ! NOT SUPPORTED BY 'cub_getCubPoints' !
  integer, parameter :: CUB_G6_3D = 308

!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, tetra">
  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer, parameter :: CUB_G1_3D_T = 351
  
  ! trapezoidal rule, 3D, degree = 2, ncubp = 4
  integer, parameter :: CUB_TRZ_3D_T = 352
  
  ! 4-point Stroud rule, 3D, degree = 2, ncubp = 4
  integer, parameter :: CUB_S2_3D_T = 353
  
  ! 10-point Stroud rule, 3D, degree = 3, ncubp = 10
  integer, parameter :: CUB_S3_3D_T = 354
  
  ! 15-point Stroud rule, 3D, degree = 5, ncubp = 15
  integer, parameter :: CUB_S5_3D_T = 355
  
!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, pyramid">
  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer, parameter :: CUB_G1_3D_Y = 401
  
  ! trapezoidal rule, 3D, degree = 2, ncubp = 5
  integer, parameter :: CUB_TRZ_3D_Y = 402
  
!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, prism">
  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  integer, parameter :: CUB_G1_3D_R = 451
  
  ! trapezoidal rule, 3D, degree = 2, ncubp = 6
  integer, parameter :: CUB_TRZ_3D_R = 452
  
  ! 3x2-point Gauss formula, 3D, degree = 3 (maybe even 4?), ncubp = 6
  ! This formula is the 'cross-product' of the 2D CUB_G3_T formula
  ! for triangles and the 1D CUB_G2_1D formula.
  integer, parameter :: CUB_G2_3D_R = 453
  
!</constantblock>

  ! DEPRECATED: maximal size of cubature node field
  integer, parameter :: CUB_MAXCUBP = 36

  ! DEPRECATED: maximal number of cubature points in 1D
  integer, parameter :: CUB_MAXCUBP_1D = 6
  
  ! DEPRECATED: maximal number of cubature points in 2D
  integer, parameter :: CUB_MAXCUBP_2D = 36
  
  ! DEPRECATED: maximal number of cubature points in 3D
  integer, parameter :: CUB_MAXCUBP_3D = 27

!</constants>

contains

  !****************************************************************************

!<function>  
  integer function cub_igetID(scubName)
  
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

  !</input>
  
!</function>

  select case(trim(sys_upcase(scubName)))

  ! 1D-formulas
  case("G1_1D")
    cub_igetID=CUB_G1_1D
  case("TRZ_1D")
    cub_igetID =CUB_TRZ_1D
  case("G2_1D")
    cub_igetID=CUB_G2_1D
  case("G3_1D")
    cub_igetID=CUB_G3_1D
  case("G4_1D")
    cub_igetID=CUB_G4_1D
  case("G5_1D")
    cub_igetID=CUB_G5_1D
  case("SIMPSON_1D")
    cub_igetID=CUB_SIMPSON_1D
  case("G6_1D")
    cub_igetID=CUB_G6_1D

  ! 2D-fomulas, quadrilateral
  case ("G1X1","G1_2D")
    cub_igetID=CUB_G1_2D
  case ("TRZ")
    cub_igetID=CUB_TRZ
  case ("MID")
    cub_igetID=CUB_MID
  case ("G2X2","G2_2D")
    cub_igetID=CUB_G2_2D
  case ("NS1")
    cub_igetID=CUB_NS1
  case ("NS2")
    cub_igetID=CUB_NS2
  case ("NS3")
    cub_igetID=CUB_NS3
  case ("G3X3","G3_2D")
    cub_igetID=CUB_G3_2D
  case ("G")
    cub_igetID=CUB_G
  case ("G4X4","G4_2D")
    cub_igetID=CUB_G4_2D
  case ("G5X5","G5_2D")
    cub_igetID=CUB_G5_2D
  case ("PG1X1")
    cub_igetID=CUB_PG1X1
  case ("PTRZ")
    cub_igetID=CUB_PTRZ
  case ("PG2X2")
    cub_igetID=CUB_PG2X2
  case ("PG3X3")
    cub_igetID=CUB_PG3X3
  case ("G6_2D")
    cub_igetID=CUB_G6_2D
    
  ! 2D-formulas, triangle
  case ("G1_T")
    cub_igetID=CUB_G1_T
  case ("TRZ_T")
    cub_igetID=CUB_TRZ_T
  case ("G3_T")
    cub_igetID=CUB_G3_T
  case ("COLLATZ")
    cub_igetID=CUB_COLLATZ
  case ("VMC")
    cub_igetID=CUB_VMC

  ! 3D-formulas, hexahedron
  case("G1_3D")
    cub_igetID=CUB_G1_3D
  case("MIDAREA_3D")
    cub_igetID =CUB_MIDAREA_3D
  case("TRZ_3D")
    cub_igetID=CUB_TRZ_3D
  case("G2_3D")
    cub_igetID=CUB_G2_3D
  case("G3_3D")
    cub_igetID=CUB_G3_3D
  case("G4_3D")
    cub_igetID=CUB_G4_3D
  case("G5_3D")
    cub_igetID=CUB_G5_3D
  case("G6_3D")
    cub_igetID=CUB_G6_3D
  
  ! 3D-formulas, tetrahedron
  case("G1_3D_T")
    cub_igetID=CUB_G1_3D_T
  case("TRZ_3D_T")
    cub_igetID=CUB_TRZ_3D_T
  case("S2_3D_T")
    cub_igetID=CUB_S2_3D_T
  case("S3_3D_T")
    cub_igetID=CUB_S3_3D_T
  case("S5_3D_T")
    cub_igetID=CUB_S5_3D_T
  
  ! 3D-formulas, pyramid
  case("G1_3D_Y")
    cub_igetID=CUB_G1_3D_Y
  case("TRZ_3D_Y")
    cub_igetID=CUB_TRZ_3D_Y
  
  ! 3D-formulas, prism
  case("G1_3D_R")
    cub_igetID=CUB_G1_3D_R
  case("TRZ_3D_R")
    cub_igetID=CUB_TRZ_3D_R
  case("G2_3D_R")
    cub_igetID=CUB_G2_3D_R

  case default
    print *,'Error: Unknown cubature formula: ',scubname
    call sys_halt()
  end select
    
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
  integer, intent(IN) :: ccubType
!</input>
  
!</function>

    select case(ccubType)
    ! -= 1D Line Formulas =-
    case (CUB_G1_1D)
      n = 1
    case (CUB_G2_1D,CUB_TRZ_1D)
      n = 2
    case (CUB_G3_1D,CUB_SIMPSON_1D)
      n = 3
    case (CUB_G4_1D)
      n = 4
    case (CUB_G5_1D)
      n = 5
    case (CUB_G6_1D)
      n = 6
    
    ! -= 2D Triangle Formulas =-
    case (CUB_G1_T)
      n = 1
    case (CUB_G3_T,CUB_TRZ_T,CUB_Collatz)
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
      print *, 'Error: Unknown cubature formula'
      call sys_halt()
    end select

  end function


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
  integer, intent(IN) :: ccubType
!</input>
  
!</function>

    if(ccubType .le. 100) then
      ! Invalid identifier
      ishp = BGEOM_SHAPE_UNKNOWN
    else if(ccubType .le. 200) then
      ! 1D Line
      ishp = BGEOM_SHAPE_LINE
    else if(ccubType .lt. 250) then
      ! 2D Quadrilateral
      ishp = BGEOM_SHAPE_QUAD
    else if(ccubType .le. 300) then
      ! 2D Triangle
      ishp = BGEOM_SHAPE_TRIA
    else if(ccubType .le. 350) then
      ! 3D Hexahedron
      ishp = BGEOM_SHAPE_HEXA
    else if(ccubType .le. 400) then
      ! 3D Tetrahedron
      ishp = BGEOM_SHAPE_TETRA
    else if(ccubType .le. 450) then
      ! 3D Pyramid
      ishp = BGEOM_SHAPE_PYRA
    else if(ccubType .le. 500) then
      ! 3D Prism
      ishp = BGEOM_SHAPE_PRISM
    else
      ! Invalid identifier
      ishp = BGEOM_SHAPE_UNKNOWN
    end if


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
  integer, intent(IN) :: ccubType
!</input>
  
!</function>

    ! Get the shape of the cubature id
    select case(cub_igetShape(ccubType))
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

  end function

  !****************************************************************************

!<subroutine>

  subroutine cub_getCubature(ccubType, Dpoints, Domega)

!<descrption>
  ! This routine returns the cubature point coordinates and the corresponding
  ! weights for a given cubature type identifier.
!</descrption>

!<input>
  ! id of the cubature formula to be set
  integer, intent(IN) :: ccubType
!</input>

!<output>
  ! Coordinates of the cubature points.
  ! The dimension is assumed to be at least Dpoints(1:NDIM,1:NPTS), where:
  ! NDIM is the dimension of a point, which is returned by cub_igetCoordDim.
  ! NPTS is the number of cubatue points, which is returned by cub_igetNumPts.
  real(DP), dimension(:,:), intent(OUT) :: Dpoints
  
  ! For every cubature point the corresponding cubature weight.
  ! The dimension is assumed to be at least Domega(1:NPTS) (see Dpoints).
  real(DP), dimension(:), intent(OUT) :: Domega
!</output>

!</subroutine>

  ! local variables for wrapper
  real(DP), dimension(CUB_MAXCUBP,4) :: Dxi
  integer :: i,j,k,l,ncubp
  
  ! Auxiliary arrays for Gauss-Legendre rules
  real(DP), dimension(5) :: Dv, Dw
  integer :: ndim, npts
  
    ndim = 0
    npts = 0
    
    ! Okay, let's see what cubature rule we have here...
    select case(ccubType)
    
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
          Dpoints(i,j) = Dxi(i,j)
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
      case (1)
        ! 1D Gauss rule
        do i = 1, npts
          Dpoints(1,i) = Dv(i)
          Domega(i) = Dw(i)
        end do
      
      case (2)
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
      
      case(3)
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
    
    ! That's it
    
  contains
  
  ! -------------------------------------------------------
  
    pure subroutine cub_auxGaussLegendre(n,Dv,Dw)
    integer, intent(IN) :: n
    real(DP), dimension(:), intent(OUT) :: Dv, Dw
    
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
        Dv(1) = -sqrt(3.0_DP + daux) / 7.0_DP
        Dv(2) = -sqrt(3.0_DP - daux) / 7.0_DP
        Dv(3) =  sqrt(3.0_DP - daux) / 7.0_DP
        Dv(4) =  sqrt(3.0_DP + daux) / 7.0_DP
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
        Dw(2) = (322.0_DP - daux) / 900.0_DP
        Dw(3) =  128.0_DP / 225.0_DP
        Dw(4) = (322.0_DP - daux) / 900.0_DP
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

  end subroutine cub_getCubature

  !****************************************************************************

!<subroutine>

  subroutine cub_getCubPoints(ccubType, ncubp, Dxi, Domega)

!<description>
  ! This routine initializes the coordinates and weight fields according 
  ! to the selected cubature formula. The integration domain is $[-1,1]^n$.
  ! In the case of one-dimensional integration, only cub_dxi(i,1) is used. 
  ! The coordinates of the cubature points on triangles are given in 
  ! barycentric coordinates.
!</description>

!<input>
  ! id of the cubature formula to be set
  integer, intent(IN) :: ccubType
!</input>
  
!<output>
  ! number of cubature points; =0: error, unknown cubature formula
  integer , intent(OUT) :: ncubp
  
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
  real(DP), dimension(:,:), intent(OUT) :: Dxi
  
  ! For every cubature point the corresponding cubature weight
  real(DP), dimension(:), intent(OUT) :: Domega
  
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
    Dxi (1,1) = -1.0_DP
    Dxi (1,2) = -1.0_DP
    Dxi (2,1) =  0.0_DP
    Dxi (2,2) = -1.0_DP
    Dxi (3,1) =  0.0_DP
    Dxi (3,2) =  0.0_DP
    Dxi (4,1) = -1.0_DP
    Dxi (4,2) =  0.0_DP
    
    Dxi (5,1) =  1.0_DP
    Dxi (5,2) = -1.0_DP
    Dxi (6,1) =  1.0_DP
    Dxi (6,2) =  0.0_DP
    Dxi (7,1) =  0.0_DP
    Dxi (7,2) =  0.0_DP
    Dxi (8,1) =  0.0_DP
    Dxi (8,2) = -1.0_DP
    
    Dxi (9,1) =  1.0_DP
    Dxi (9,2) =  1.0_DP
    Dxi(10,1) =  0.0_DP
    Dxi(10,2) =  1.0_DP
    Dxi(11,1) =  0.0_DP
    Dxi(11,2) =  0.0_DP
    Dxi(12,1) =  1.0_DP
    Dxi(12,2) =  0.0_DP
    
    Dxi(13,1) = -1.0_DP
    Dxi(13,2) =  1.0_DP
    Dxi(14,1) = -1.0_DP
    Dxi(14,2) =  0.0_DP
    Dxi(15,1) =  0.0_DP
    Dxi(15,2) =  0.0_DP
    Dxi(16,1) =  0.0_DP
    Dxi(16,2) =  1.0_DP
    
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
    
  case(CUB_G3_T)
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
    
  case(CUB_Collatz)
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
    print *,'Error: unknown cubature formula: ',ccubType
    call sys_halt()
  end select
   
  end subroutine cub_getCubPoints

end module
