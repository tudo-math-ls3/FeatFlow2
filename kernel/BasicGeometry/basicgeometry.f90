!##############################################################################
!# ****************************************************************************
!# <name> basicgeometry </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic entity of the whole FEAT world: a point,
!# a coordinate system and a set of routines to maintain that.
!#
!# The following routines can be found in this module:
!#
!# 1.) bgeom_transformPoint2D
!#     -> Transforms a 2D point with coordinates relative to a given coordinate
!#        system to world coordinates.
!#
!# 2.) bgeom_transformBackPoint2D
!#     -> Transforms a 2D point with world coordinates to coordinates relative
!#        to a given coordinate system.
!# </purpose>
!##############################################################################

MODULE basicgeometry

  USE fsystem

  IMPLICIT NONE


  ! One could ask: Why do we need the following two constants...?
  ! The answer is simple: Then we can more easily GREP for dimension-dependent
  ! quantities :-)
  
!<constants>
!<constantblock description="Dimension constants">
  
  ! Dimension constant for 1D triangulations.
  INTEGER, PARAMETER :: NDIM1D = 1

  ! Dimension constant for 2D triangulations.
  INTEGER, PARAMETER :: NDIM2D = 2

  ! Dimension constant for 3D triangulations.
  INTEGER, PARAMETER :: NDIM3D = 3

!</constantblock>
!</constants>
  
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
  TYPE t_coordinateSystem2D
  
    ! Coordinates of the origin
    REAL(DP), DIMENSION(2) :: Dorigin = (/0.0_DP,0.0_DP/)
    
    ! Rotation angle; 0..2*PI
    REAL(DP) :: drotation = 0.0_DP
    
    ! precalculated value: sin(rotation); for quicker calculations
    REAL(DP) :: dsin_rotation = 0.0_DP

    ! precalculated value: cos(rotation); for quicker calculations
    REAL(DP) :: dcos_rotation = 1.0_DP
    
    ! scaling factor of the coordinate system; usually = 1.0
    REAL(DP) :: dscalingFactor = 1.0_DP
    
  END TYPE

  !</typeblock>
  !</types>

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE bgeom_transformPoint2D(rcoordSys, DpointIn, DpointOut)

!<description>
  ! Transforms a given 2D point with coordinates relative to the given
  ! coordinate system to world coordinates.
!</description>

!<input>
  ! The 2D coordinate system that is to be applied to the given point.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys

  ! The 2D point that is to be transformed.
  REAL(DP), DIMENSION(:),      INTENT(IN)  :: DpointIn
!</input>

!<output>
  ! The transformed 2D point.
  REAL(DP), DIMENSION(:),      INTENT(OUT) :: DpointOut
!</output>

!</subroutine>

    ! Check if the rotation of the coordinate system is non-zero.
    IF (rcoordSys%drotation .NE. 0.0_DP) THEN

      ! Apply rotation to the input coordinates.
      DpointOut(1) = (rcoordSys%dcos_rotation * DpointIn(1)) - &
                     (rcoordSys%dsin_rotation * DpointIn(2))
      DpointOut(2) = (rcoordSys%dsin_rotation * DpointIn(1)) + &
                     (rcoordSys%dcos_rotation * DpointIn(2))

    ELSE

      ! No rotation in the coordinate system, so simply copy the input coords.
      DpointOut(1) = DpointIn(1)
      DpointOut(2) = DpointIn(2)

    END IF

    ! Now scale the given coordinates by the scaling factor and translate them
    ! by the origin of our coordinate system.
    DpointOut(1) = (rcoordSys%dscalingFactor * DpointOut(1)) + &
                    rcoordSys%Dorigin(1)
    DpointOut(2) = (rcoordSys%dscalingFactor * DpointOut(2)) + &
                    rcoordSys%Dorigin(2)

    ! That's it!

  END SUBROUTINE


  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE bgeom_transformBackPoint2D(rcoordSys, DpointIn, DpointOut)

!<description>
  ! Transform a 2D point given in world coordinates to coordinates relative
  ! to a given 2D coordinate system.
  ! This subroutine is the inverse routine of bgeom_transformPoint2D.
!</description>

!<input>
  ! The 2D coordinate system of our input point.
  TYPE(t_coordinateSystem2D), INTENT(IN)  :: rcoordSys

  ! The 2D point that is to be transformed, relative to the given coordinate
  ! system.
  REAL(DP), DIMENSION(:),     INTENT(IN)  :: DpointIn
!</input>

!<output>
  ! The transformed 2D point, in world coordinates.
  REAL(DP), DIMENSION(:),     INTENT(OUT) :: DpointOut
!</output>

!</subroutine>

    ! local variables
    REAL(DP) :: X,Y

    ! Translate the given point by the negatives of our coordinate system
    ! origin.
    X = DpointIn(1) - rcoordSys%Dorigin(1)
    Y = DpointIn(2) - rcoordSys%Dorigin(2)

    ! Now scale the point by the inverse of our coordinate system sclaing
    ! factor, of course only if it is non-zero ...
    IF (rcoordSys%dscalingFactor .NE. 0.0_DP) THEN

      X = X / rcoordSys%dscalingFactor
      Y = Y / rcoordSys%dscalingFactor

    ELSE

      ! ... or set the result to zero and exit subroutine otherwise.
      DpointOut(1) = 0.0_DP
      DpointOut(2) = 0.0_DP
      RETURN

    ENDIF

    ! Check if the rotation of the coordinate system is non-zero.
    IF (rcoordSys%drotation .NE. 0.0_DP) THEN

      ! Apply rotation to the input coordinates.
      DpointOut(1) = ( rcoordSys%dcos_rotation * X) + (rcoordSys%dsin_rotation * Y)
      DpointOut(2) = (-rcoordSys%dsin_rotation * X) + (rcoordSys%dcos_rotation * Y)

    ELSE

      ! No rotation in the coordinate system, so simply copy the input coords.
      DpointOut(1) = X
      DpointOut(2) = Y

    END IF

    ! That's it

  END SUBROUTINE

END MODULE
