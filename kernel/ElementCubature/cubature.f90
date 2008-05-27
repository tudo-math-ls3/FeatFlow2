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
!# 2.) cub_getCubPoints
!#     -> Get information about cubature points for a specific cubature
!#        formula: Coordinates, number and weights.
!# </purpose>
!##############################################################################

MODULE cubature

  USE fsystem
  USE basicgeometry

  IMPLICIT NONE
  
!<constants>

!<constantblock variable="ccubType" description="1D formulas">  

  ! 1-point Gauss formula, 1D, degree = 2, ncubp = 1
  INTEGER, PARAMETER :: CUB_G1_1D = 101
 
  ! trapezoidal rule, 1D, degree = 2, ncubp = 2
  INTEGER, PARAMETER :: CUB_TRZ_1D = 102

  ! 2-point Gauss formula, 1D, degree = 4, ncubp = 2
  INTEGER, PARAMETER :: CUB_G2_1D = 103

  ! 3-point Gauss formula, 1D, degree = 6, ncubp = 3
  INTEGER, PARAMETER :: CUB_G3_1D = 104

  ! 4-point Gauss formula, 1D, degree = 8, ncubp = 4
  INTEGER, PARAMETER :: CUB_G4_1D = 105

  ! 5-point Gauss formula, 1D, degree = 10, ncubp = 5
  INTEGER, PARAMETER :: CUB_G5_1D = 106

  ! Simpson-rule, 1D, degree = 4, ncubp = 3
  INTEGER, PARAMETER :: CUB_SIMPSON_1D = 107

  ! 6-point Gauss formula, 1D, degree = 12, ncubp = 6
  INTEGER, PARAMETER :: CUB_G6_1D = 108

!</constantblock>

!<constantblock variable="ccubType" description="2D formulas, quad">  

  ! 1x1 Gauss formula, degree = 2, ncubp = 1
  INTEGER, PARAMETER :: CUB_G1X1 = 201

  ! trapezoidal rule, degree = 2, ncubp = 4
  INTEGER, PARAMETER :: CUB_TRZ = 202

  ! midpoint rule, degree = 2, ncubp = 4
  INTEGER, PARAMETER :: CUB_MID = 203

  ! 2x2 Gauss formula, degree = 4, ncubp = 4
  INTEGER, PARAMETER :: CUB_G2X2 = 204

  ! Newton formula 1, degree = 4, ncubp = 4 
  INTEGER, PARAMETER :: CUB_NS1 = 205

  ! Newton formula 2, degree = 5, ncubp = 6
  INTEGER, PARAMETER :: CUB_NS2 = 206

  ! Newton formula 3, degree = 6, ncubp = 7
  INTEGER, PARAMETER :: CUB_NS3 = 207

  ! 3x3 Gauss formula, degree = 6, ncubp = 9
  INTEGER, PARAMETER :: CUB_G3X3 = 208

  ! Gauss formula, degree = 7, ncubp = 12
  INTEGER, PARAMETER :: CUB_G = 209

  ! 4x4 Gauss formula, degree = 8, ncubp = 16
  INTEGER, PARAMETER :: CUB_G4X4 = 210

  ! 5x5 Gauss formula, degree = 10, ncubp = 25
  INTEGER, PARAMETER :: CUB_G5X5 = 211

  ! piecewise 1x1 Gauss formula, degree = 2, ncubp = 4
  INTEGER, PARAMETER :: CUB_PG1X1 = 212

  ! piecewise trapezoidal rule, degree = 2, ncubp = 9
  INTEGER, PARAMETER :: CUB_PTRZ = 213

  ! piecewise 2x2 Gauss formula, degree = 4, ncubp = 16 
  INTEGER, PARAMETER :: CUB_PG2X2 = 214

  ! piecewise 3x3 Gauss formula, degree = 6, ncubp = 36
  INTEGER, PARAMETER :: CUB_PG3X3 = 215
  
  ! Simpson rule (corners and element midpoint), degree = 3, ncubp = 9
  INTEGER, PARAMETER :: CUB_SIMPSON = 216

  ! Simpson 3/8 rule (corners and 1/3 + 2/3), degree = 3, ncubp = 16
  INTEGER, PARAMETER :: CUB_3_8 = 217

!</constantblock>

!<constantblock variable="ccubType" description="2D formulas, tri">

  ! 1-point Gauss formula, triangle, degree = 2, ncubp = 1
  INTEGER, PARAMETER :: CUB_G1_T = 250

  ! trapezoidal rule, triangle, degree = 2, ncubp = 3 
  INTEGER, PARAMETER :: CUB_TRZ_T = 251

  ! 3-point Gauss formula, triangle, degree = 3, ncubp = 3
  INTEGER, PARAMETER :: CUB_G3_T = 252

  ! Collatz formula, degree = 3, ncubp = 3
  INTEGER, PARAMETER :: CUB_Collatz = 253

  ! vertices, midpoints, center, degree = 4, ncubp = 7 
  INTEGER, PARAMETER :: CUB_VMC = 254
!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, hexa">

  ! 1-point Gauss formula, 3D, degree = 2, ncubp = 1
  INTEGER, PARAMETER :: CUB_G1_3D = 301

  ! midpoints of areas, 3D, degree = 2, ncubp = 6
  INTEGER, PARAMETER :: CUB_MIDAREA_3D = 302

  ! trapezoidal rule, 3D, degree = 2, ncubp = 8
  INTEGER, PARAMETER :: CUB_TRZ_3D = 303

  ! 2-point Gauss formula, 3D, degree = 4, ncubp = 8
  INTEGER, PARAMETER :: CUB_G2_3D = 304
  
  ! 3-point Gauss formula, 3D, degree = 6, ncubp = 27
  INTEGER, PARAMETER :: CUB_G3_3D = 305

!</constantblock>

!<constantblock variable="ccubType" description="3D formulas, tetra">
  ! 1-point Gauss formula, 3D, degree = 1, ncubp = 1
  INTEGER, PARAMETER :: CUB_G1_3D_T = 350
  
  ! trapezoidal rule, 3D, degree = 2, ncubp = 4
  INTEGER, PARAMETER :: CUB_TRZ_3D_T = 351
  
  ! 4-point Stroud rule, 3D, degree = 2, ncubp = 4
  INTEGER, PARAMETER :: CUB_S2_3D_T = 352
  
  ! 10-point Stroud rule, 3D, degree = 3, ncubp = 10
  INTEGER, PARAMETER :: CUB_S3_3D_T = 353
  
  ! 15-point Stroud rule, 3D, degree = 5, ncubp = 15
  INTEGER, PARAMETER :: CUB_S5_3D_T = 354
  
!</constantblock>

  ! maximal size of cubature node field
  INTEGER, PARAMETER :: CUB_MAXCUBP = 36

  ! maximal number of cubature points in 1D
  INTEGER, PARAMETER :: CUB_MAXCUBP_1D = 6
  
  ! maximal number of cubature points in 2D
  INTEGER, PARAMETER :: CUB_MAXCUBP_2D = 36
  
  ! maximal number of cubature points in 3D
  INTEGER, PARAMETER :: CUB_MAXCUBP_3D = 15

!</constants>

CONTAINS

  !****************************************************************************

!<function>  
  INTEGER FUNCTION cub_igetID(scubName)
  
!<description>
  ! This routine returns the cubature id to a given cubature formula name. It is 
  ! case-insensitive. 
!</description>

!<result>
  ! id of the cubature formula
!</result>

  !<input>

  !cubature formula name - one of the CUB_xxxx constants.
  CHARACTER (LEN=*) :: scubName

  !</input>
  
!</function>

  SELECT CASE(TRIM(sys_upcase(scubName)))

  ! 1D-formulas
  CASE("G1_1D")
    cub_igetID=CUB_G1_1D
  CASE("TRZ_1D")
    cub_igetID =CUB_TRZ_1D
  CASE("G2_1D")
    cub_igetID=CUB_G2_1D
  CASE("G3_1D")
    cub_igetID=CUB_G3_1D
  CASE("G4_1D")
    cub_igetID=CUB_G4_1D
  CASE("G5_1D")
    cub_igetID=CUB_G5_1D
  CASE("SIMPSON_1D")
    cub_igetID=CUB_SIMPSON_1D
  CASE("G6_1D")
    cub_igetID=CUB_G6_1D

  ! 2D-fomulas, quadrilateral
  CASE ("G1X1")
    cub_igetID=CUB_G1X1
  CASE ("TRZ")
    cub_igetID=CUB_TRZ
  CASE ("MID")
    cub_igetID=CUB_MID
  CASE ("G2X2")
    cub_igetID=CUB_G2X2
  CASE ("NS1")
    cub_igetID=CUB_NS1
  CASE ("NS2")
    cub_igetID=CUB_NS2
  CASE ("NS3")
    cub_igetID=CUB_NS3
  CASE ("G3X3")
    cub_igetID=CUB_G3X3
  CASE ("G")
    cub_igetID=CUB_G
  CASE ("G4X4")
    cub_igetID=CUB_G4X4
  CASE ("G5X5")
    cub_igetID=CUB_G5X5
  CASE ("PG1X1")
    cub_igetID=CUB_PG1X1
  CASE ("PTRZ")
    cub_igetID=CUB_PTRZ
  CASE ("PG2X2")
    cub_igetID=CUB_PG2X2
  CASE ("PG3X3")
    cub_igetID=CUB_PG3X3
    
  ! 2D-formulas, triangle
  CASE ("G1_T")
    cub_igetID=CUB_G1_T
  CASE ("TRZ_T")
    cub_igetID=CUB_TRZ_T
  CASE ("G3_T")
    cub_igetID=CUB_G3_T
  CASE ("COLLATZ")
    cub_igetID=CUB_COLLATZ
  CASE ("VMC")
    cub_igetID=CUB_VMC

  ! 3D-formulas, hexahedron
  CASE("G1_3D")
    cub_igetID=CUB_G1_3D
  CASE("MIDAREA_3D")
    cub_igetID =CUB_MIDAREA_3D
  CASE("TRZ_3D")
    cub_igetID=CUB_TRZ_3D
  CASE("G2_3D")
    cub_igetID=CUB_G2_3D
  CASE("G3_3D")
    cub_igetID=CUB_G3_3D
  
  ! 3D-formulas, tetrahedron
  CASE("G1_3D_T")
    cub_igetID=CUB_G1_3D_T
  CASE("TRZ_3D_T")
    cub_igetID=CUB_TRZ_3D_T
  CASE("S2_3D_T")
    cub_igetID=CUB_S2_3D_T
  CASE("S3_3D_T")
    cub_igetID=CUB_S3_3D_T
  CASE("S5_3D_T")
    cub_igetID=CUB_S5_3D_T

  CASE DEFAULT
    PRINT *,'Error: Unknown cubature formula: ',scubname
    CALL sys_halt()
  END SELECT
    
  END FUNCTION 

  !****************************************************************************

!<subroutine>

  SUBROUTINE cub_getCubPoints(ccubType, ncubp, Dxi, Domega)

!<description>
  ! This routine initializes the coordinates and weight fields according 
  ! to the selected cubature formula. The integration domain is $[-1,1]^n$.
  ! In the case of one-dimensional integration, only cub_dxi(i,1) is used. 
  ! The coordinates of the cubature points on triangles are given in 
  ! barycentric coordinates.
!</description>

!<input>
  ! id of the cubature formula to be set
  INTEGER, INTENT(IN) :: ccubType
!</input>
  
!<output>
  ! number of cubature points; =0: error, unknown cubature formula
  INTEGER , INTENT(OUT) :: ncubp
  
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
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dxi
  
  ! For every cubature point the corresponding cubature weight
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Domega
  
!</output>

!</subroutine>    

  ! local variables
  INTEGER :: i 
  
  SELECT CASE (ccubType)

  !1D cubature formulas
  CASE (CUB_G1_1D)
    Dxi(1,1)  = 0.0_DP
     
    Domega(1) = 2.0_DP
     
    ncubp     = 1
     
  CASE(CUB_TRZ_1D)
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  =  1.0_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    
    ncubp     =  2
     
  CASE(CUB_G2_1D)
    Dxi(1,1)  = -0.577350269189626_DP
    Dxi(2,1)  =  0.577350269189626_DP
    
    Domega(1) =  1.0_DP
    Domega(2) =  1.0_DP
    
    ncubp     =  2
     
  CASE(CUB_G3_1D)
    Dxi(1,1)  = -0.774596669241483_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(3,1)  =  0.774596669241483_DP
    
    Domega(1) =  0.555555555555556_DP
    Domega(2) =  0.888888888888889_DP
    Domega(3) =  0.555555555555556_DP
    
    ncubp     =  3
    
  CASE(CUB_G4_1D)
    Dxi(1,1)  = -0.861136311594053_DP
    Dxi(2,1)  = -0.339981043584856_DP
    Dxi(3,1)  =  0.339981043584856_DP
    Dxi(4,1)  =  0.861136311594053_DP
    
    Domega(1) =  0.347854845137454_DP
    Domega(2) =  0.652145154862546_DP
    Domega(3) =  0.652145154862546_DP
    Domega(4) =  0.347854845137454_DP
    
    ncubp     =  4
    
  CASE(CUB_G5_1D)
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
    
  CASE(CUB_SIMPSON_1D) 
    Dxi(1,1)  = -1.0_DP
    Dxi(2,1)  =  0.0_DP
    Dxi(3,1)  =  1.0_DP
    
    Domega(1) =  0.333333333333333_DP
    Domega(2) =  1.333333333333333_DP
    Domega(3) =  0.333333333333333_DP
    
    ncubp     =  3
    
  CASE(CUB_G6_1D)
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
  CASE (CUB_G1X1)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.0_DP
    
    Domega(1) =  4.0_DP
    
    ncubp     =  1
    
  CASE (CUB_TRZ)
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
    
  CASE (CUB_MID)
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
     
  CASE (CUB_G2X2)
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
     
  CASE (CUB_NS1)
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
    
  CASE (CUB_NS2)
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

  CASE (CUB_NS3)
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

  CASE (CUB_G3X3)
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
    
  CASE (CUB_G)
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
    
  CASE (CUB_G4X4)
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

  CASE (CUB_G5X5)
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
    
  CASE (CUB_PG1X1)
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

  CASE (CUB_PTRZ)
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

  CASE (CUB_PG2X2)
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
    
  CASE (CUB_PG3X3)
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
    
  CASE(CUB_SIMPSON)
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

  CASE(CUB_3_8)
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
  CASE(CUB_G1_T)
    Dxi(1,1)  =  0.3333333333333333_DP
    Dxi(1,2)  =  0.3333333333333333_DP
    Dxi(1,3)  =  0.3333333333333333_DP
    
    Domega(1) =  0.5_DP
    
    ncubp     =  1
    
  CASE(CUB_TRZ_T)
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
    
  CASE(CUB_G3_T)
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
    
  CASE(CUB_Collatz)
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
    
  CASE(CUB_VMC)
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

  CASE(CUB_G1_3D)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(1,3)  =  0.0_DP
    
    Domega(1) =  8.0_DP
    
    ncubp      =  1
    
  CASE(CUB_MIDAREA_3D)
    Dxi(1,1)  =  0.0_DP
    Dxi(1,2)  =  0.0_DP
    Dxi(1,3)  = -1.0_DP
    
    Dxi(2,1)  = -1.0_DP
    Dxi(2,2)  =  0.0_DP
    Dxi(2,3)  =  0.0_DP

    Dxi(3,1)  =  0.0_DP
    Dxi(3,2)  = -1.0_DP
    Dxi(3,3)  =  0.0_DP

    Dxi(4,1)  =  1.0_DP
    Dxi(4,2)  =  0.0_DP
    Dxi(4,3)  =  0.0_DP

    Dxi(5,1)  =  0.0_DP
    Dxi(5,2)  =  1.0_DP
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

  CASE(CUB_TRZ_3D)
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
    
  CASE(CUB_G2_3D)

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

  CASE(CUB_G3_3D)
    
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
  CASE(CUB_G1_3D_T)
    Dxi(1,1)  =  0.25_DP
    Dxi(1,2)  =  0.25_DP
    Dxi(1,3)  =  0.25_DP
    Dxi(1,4)  =  0.25_DP
    
    Domega(1) = 0.1666666666666667_DP
    
    ncubp = 1
    
  CASE(CUB_TRZ_3D_T)
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

  CASE(CUB_S2_3D_T)
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
  
  CASE(CUB_S3_3D_T)
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
  
  CASE(CUB_S5_3D_T)
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

  CASE DEFAULT 
    PRINT *,'Error: unknown cubature formula: ',ccubType
    CALL sys_halt()
  END SELECT
   
  END SUBROUTINE cub_getCubPoints

END MODULE
