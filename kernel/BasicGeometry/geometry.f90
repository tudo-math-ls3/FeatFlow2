!##############################################################################
!# ****************************************************************************
!# <name> geometry </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises some basic 2D/3D geometry objects...
!# Foobar...
!#
!# The following routines can be found in this module:
!#
!#  1.) geom_init_circle
!#      -> Creates a 2D circle object.
!#
!#  2.) geom_init_square
!#      -> Creates a 2D square object.
!#
!#  3.) geom_done
!#      -> Destroys a geometric object.
!#
!# </purpose>
!##############################################################################
 
MODULE geometry

  USE basicgeometry

  IMPLICIT NONE

!<constants>

!<constantblock description="Geometry type identifiers">

  ! Composed geometry object with sub-objects
  INTEGER, PARAMETER :: GEOM_COMPOSED  = -1

  ! No geometry object
  INTEGER, PARAMETER :: GEOM_NONE      = 0

  ! Circle object
  INTEGER, PARAMETER :: GEOM_CIRCLE    = 1

  ! Square object
  INTEGER, PARAMETER :: GEOM_SQUARE    = 2

  ! Ellipse object
  INTEGER, PARAMETER :: GEOM_ELLIPSE   = 3

  ! Rectangle object
  INTEGER, PARAMETER :: GEOM_RECTANGLE = 4

!</constantblock>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<types>
!<typeblock>

  ! This structure realises the subnode for the composed geometry object.
  
  TYPE t_geometryComposed
    
    ! Number of sub-objects in the composed object.
    INTEGER :: nsubobjects
    
    ! An array of sub-objects with are composed.
    TYPE(t_geometryObject), DIMENSION(:), POINTER :: rsubobjects
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the circle object.
  
  TYPE t_geometryCircle
    
    ! Radius of the circle. Must be positive.
    REAL(DP) :: dradius
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>

  ! This structure realises the subnode for the square object
  TYPE t_geometrySquare
  
    ! Length of each edge of the square. Must be positive.
    REAL(DP) :: dlength
    
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! Geometry object structure for 2D and 3D geometry objects.
  TYPE t_geometryObject

    ! Type identifier of the geometry object.
    ! One of the GEOM_**** constants defined above.
    INTEGER                    :: ctype = GEOM_NONE

    ! Dimension of the coordinate system.
    ! May be NDIM2D or NDIM3D (as defined in basicgeometry.f90).
    INTEGER                    :: ndimension = NDIM2D
    
    ! The 2D coordinate system for this geometry object.
    ! Used for 2D geometry objects - undefined for 3D geometry objects.
    TYPE(t_coordinateSystem2D) :: rcoord2D
    
    ! The 3D coordinate system for this geometry object.
    ! Used for 3D geometry objects - undefined for 2D geometry objects.
    !TYPE(t_coordinateSystem3D) :: rcoord3D
    
    ! A boolean which tells us whether the geometry object is inverted or not.
    LOGICAL                    :: binverted
    
    ! Structure for the composed geometry object
    TYPE(t_geometryComposed)   :: rcomposed
    
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! = 2D object subnodes -
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! 2DStructure for the circle object (2D)
    TYPE(t_geometryCircle)     :: rcircle
    
    ! Structure for the square object (2D)
    TYPE(t_geometrySquare)     :: rsquare
    
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! = 3D object subnodes -
    ! -=-=-=-=-=-=-=-=-=-=-=
    ! Todo...
    
  END TYPE
  
!</typeblock>
!</types>

  INTERFACE geom_init_circle
    MODULE PROCEDURE geom_init_circle
    MODULE PROCEDURE geom_init_circle_2
  END INTERFACE
  
  INTERFACE geom_init_square
    MODULE PROCEDURE geom_init_square
    MODULE PROCEDURE geom_init_square_2
  END INTERFACE

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE geom_init_circle(rgeomObject, rcoordSys, dradius, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D circle.
!</description>

!<input>
  ! A 2D coordinate system for the circle.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! A radius for the circle.
  REAL(DP),                    INTENT(IN)  :: dradius
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! Destroy the previous geometry object.
    !geom_done(rgeomObject)

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a circle.
    rgeomObject%ctype = GEOM_CIRCLE
    
    ! Store the coordinate system.
    rgeomObject%rcoord2D = rcoordSys
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the radius of the circle
    rgeomObject%rcircle%dradius = dradius
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************

  SUBROUTINE geom_init_circle_2(rgeomObject, dradius, Dorigin, drotation, &
                                dscalingFactor, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D circle.
!</description>

!<input>
  ! A radius for the circle.
  REAL(DP),                          INTENT(IN)  :: dradius
  
  ! OPTIONAL: The origin of the circle.
  ! Is set to (/ 0.0_DP, 0.0_DP /) if not given.
  REAL(DP), DIMENSION(:), OPTIONAL,  INTENT(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation of the circle.
  ! Is set to 0.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: drotation
  
  ! OPTIONAL: The scaling factor of the circle.
  ! Is set to 1.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: dscalingFactor
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,                 INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),            INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! Destroy the previous geometry object.
    !geom_done(rgeomObject)

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a circle.
    rgeomObject%ctype = GEOM_CIRCLE
    
    ! Now we need to create the coordinate system.
    CALL bgeom_initCoordSys2D (rgeomObject%rcoord2D, Dorigin, drotation, &
                                 dscalingFactor)
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the radius of the circle
    rgeomObject%rcircle%dradius = dradius
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
    
!<subroutine>

  SUBROUTINE geom_init_square(rgeomObject, rcoordSys, dlength, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D square.
!</description>

!<input>
  ! A 2D coordinate system for the square.
  TYPE(t_coordinateSystem2D),  INTENT(IN)  :: rcoordSys
  
  ! The edge length for the square.
  REAL(DP),                    INTENT(IN)  :: dlength
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,           INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),      INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! Destroy the previous geometry object.
    !geom_done(rgeomObject)

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a square.
    rgeomObject%ctype = GEOM_SQUARE
    
    ! Store the coordinate system.
    rgeomObject%rcoord2D = rcoordSys
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the edge length of the square.
    rgeomObject%rsquare%dlength = dlength
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
  SUBROUTINE geom_init_square_2(rgeomObject, dlength, Dorigin, drotation, &
                                dscalingFactor, binverted)

!<description>
  ! Creates a t_geometryObject representing a 2D square.
!</description>

!<input>
  ! The edge length for the square.
  REAL(DP),                          INTENT(IN)  :: dlength
  
  ! OPTIONAL: The origin of the square.
  ! Is set to (/ 0.0_DP, 0.0_DP /) if not given.
  REAL(DP), DIMENSION(:), OPTIONAL,  INTENT(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation of the square.
  ! Is set to 0.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: drotation
  
  ! OPTIONAL: The scaling factor of the square.
  ! Is set to 1.0_DP if not given.
  REAL(DP), OPTIONAL,                INTENT(IN)  :: dscalingFactor
  
  ! OPTIONAL: A boolean telling us whether the object is inverted.
  ! Is set to .FALSE. if not given.
  LOGICAL, OPTIONAL,                 INTENT(IN)  :: binverted

!</input>

!<output>
  ! A t_geometryObject structure to be written.
  TYPE(t_geometryObject),            INTENT(OUT) :: rgeomObject

!</output>

!</subroutine>

    ! Destroy the previous geometry object.
    !geom_done(rgeomObject)

    ! The dimension is 2D.
    rgeomObject%ndimension = NDIM2D
    
    ! We want a square.
    rgeomObject%ctype = GEOM_SQUARE
    
    ! Now we need to create the coordinate system.
    CALL bgeom_initCoordSys2D (rgeomObject%rcoord2D, Dorigin, drotation, &
                                 dscalingFactor)
    
    ! Is our object inverted?
    IF (PRESENT(binverted)) THEN
      rgeomObject%binverted = binverted
    ELSE
      rgeomObject%binverted = .FALSE.
    END IF
    
    ! Store the edge length of the square.
    rgeomObject%rsquare%dlength = dlength
    
    ! That's it!
  
  END SUBROUTINE

  ! ***************************************************************************
      
!<subroutine>

  RECURSIVE SUBROUTINE geom_done(rgeomObject)
  
!<description>
  ! This subroutine destroys a geometry object.
  ! Although this routine is redundant for simple geometry objects, it must
  ! must be called for composed geometry objects.
!</description>

!<inputoutput>
  ! The geometry object structure that is to be destroyed.
  TYPE(t_geometryObject), INTENT(INOUT) :: rgeomObject

!</inputoutput>

!</subroutine>


  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ! This routine is just a dummy for now.  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  END SUBROUTINE

END MODULE