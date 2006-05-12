!##############################################################################
!# ****************************************************************************
!# <name> BasicGeometry </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic entity of the whole FEAT world: a point,
!# a coordinate system and a set of routines to maintain that.
!# </purpose>
!##############################################################################

MODULE BasicGeometry

  USE fsystem

  IMPLICIT NONE


  ! One could ask: Why do we need the following two constants...?
  ! The answer is simple: Then we can more easily GREP for dimsnion-dependent
  ! quantities :-)
  
!<constants>
!<constantblock description="Dimension constants">
  
  ! Dimension constant for 2D triangulations.
  INTEGER, PARAMETER :: NDIM2D = 2

  ! Dimension constant for 3D triangulations.
  INTEGER, PARAMETER :: NDIM3D = 3

!</constantblock>
!<!constants>
  
  !<types>

  !<typeblock>
  
  ! The point structure for 2D points.
  TYPE t_point2D
  
    ! X-coordinate
    REAL(DP) :: X
    
    ! Y-coordinate
    REAL(DP) :: Y
  
  END TYPE

  !<typeblock>
  
  !</typeblock>
  
  ! The point structure for 3D points.
  TYPE t_point3D
  
    ! X-coordinate
    REAL(DP) :: X
    
    ! Y-coordinate
    REAL(DP) :: Y

    ! Z-coordinate
    REAL(DP) :: Z
    
  END TYPE

  !</typeblock>

  !<typeblock>
  
  ! The point structure for regular 2D coordinate systems, consisting of
  ! a rotated X- and Y-axes
  TYPE t_coordinateSystem3D
  
    ! Coordinates of the origin
    TYPE(t_point2D) :: origin = t_point2D(0.0_DP,0.0_DP)
    
    ! Rotation angle; 0..2*PI
    REAL(DP) :: rotation = 0.0_DP
    
    ! precalculated value: sin(rotation); for quicker calculations
    REAL(DP) :: sin_rotation = 0.0_DP

    ! precalculated value: cos(rotation); for quicker calculations
    REAL(DP) :: cos_rotation = 1.0_DP
    
    ! scaling factor of the coordinate system; usually = 1.0
    REAL(DP) :: scalingFactor = 1.0_DP
    
  END TYPE

  !</typeblock>
  !</types>

END MODULE