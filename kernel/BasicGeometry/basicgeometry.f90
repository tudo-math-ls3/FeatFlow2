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
!#  1.) bgeom_initCoordSys2D
!#      -> Creates a 2D coordinate system from an optionally given origin,
!#         rotation and scaling factor.
!#
!#  2.) bgeom_transformPoint2D
!#      -> Transforms a 2D point with coordinates relative to a given 
!#         coordinate system to world coordinates.
!#
!#  3.) bgeom_transformBackPoint2D
!#      -> Transforms a 2D point with world coordinates to coordinates relative
!#         to a given coordinate system.
!# </purpose>
!##############################################################################

module basicgeometry

  use fsystem

  implicit none


  ! One could ask: Why do we need the following two constants...?
  ! The answer is simple: Then we can more easily GREP for dimension-dependent
  ! quantities :-)
  
!<constants>
!<constantblock description="Dimension constants">
  
  ! Dimension constant for 1D triangulations.
  integer, parameter :: NDIM1D = 1

  ! Dimension constant for 2D triangulations.
  integer, parameter :: NDIM2D = 2

  ! Dimension constant for 3D triangulations.
  integer, parameter :: NDIM3D = 3

!</constantblock>
!</constants>
  
  !<types>

  !<typeblock>
  
  ! The point structure for 2D points.
  type t_point2D
  
    ! X-coordinate
    real(DP) :: X
    
    ! Y-coordinate
    real(DP) :: Y
  
  end type

  !<typeblock>
  
  !</typeblock>
  
  ! The point structure for 3D points.
  type t_point3D
  
    ! X-coordinate
    real(DP) :: X
    
    ! Y-coordinate
    real(DP) :: Y

    ! Z-coordinate
    real(DP) :: Z
    
  end type

  !</typeblock>

  !<typeblock>
  
  ! The point structure for regular 2D coordinate systems, consisting of
  ! a rotated X- and Y-axes
  type t_coordinateSystem2D
  
    ! Coordinates of the origin
    real(DP), dimension(2) :: Dorigin = (/0.0_DP,0.0_DP/)
    
    ! Rotation angle; 0..2*PI
    real(DP) :: drotation = 0.0_DP
    
    ! precalculated value: sin(rotation); for quicker calculations
    real(DP) :: dsin_rotation = 0.0_DP

    ! precalculated value: cos(rotation); for quicker calculations
    real(DP) :: dcos_rotation = 1.0_DP
    
    ! scaling factor of the coordinate system; usually = 1.0
    real(DP) :: dscalingFactor = 1.0_DP
    
  end type

  !</typeblock>
  !</types>

contains

  ! ***************************************************************************
  
!<subroutine>
  
  pure subroutine bgeom_initCoordSys2D(rcoordSys, Dorigin, drotation, &
                                       dscalingFactor)

!<description>
  ! Creates a 2D coordinate system from an optionally given origin, rotation
  ! and scaling factor.
!</description>

!<input>
  ! OPTIONAL: The origin of the coordinate system.
  ! Set to (/ 0.0_DP, 0.0_DP /) if not given.
  real(DP), dimension(:), optional,    intent(IN)  :: Dorigin
  
  ! OPTIONAL: The rotation angle of the coordinate system.
  ! Angle range: 0..2*PI
  ! Set to 0.0_DP if not given.
  real(DP), optional,                  intent(IN)  :: drotation
  
  ! The scaling factor. Should be != 0.0_DP
  ! Set to 1.0_DP if not given.
  real(DP), optional,                  intent(IN)  :: dscalingFactor

!</input>

!<output>
  ! A t_coordinateSystem2D structure to be written.
  type(t_coordinateSystem2D),          intent(OUT) :: rcoordSys

!</output>

!</subroutine>

    ! Set the origin, if given.
    if (present(Dorigin)) then
      rcoordSys%Dorigin = (/ Dorigin(1), Dorigin(2) /)
    else
      rcoordSys%Dorigin = (/ 0.0_DP, 0.0_DP /)
    end if
      
    ! Set the rotation, if given.
    if (present(drotation)) then
      
      rcoordSys%drotation = drotation
        
      ! Calculate SIN and COS of rotation angle
      rcoordSys%dsin_rotation = sin(drotation)
      rcoordSys%dcos_rotation = cos(drotation)
        
    else
      
      rcoordSys%drotation = 0.0_DP
        
      ! Set SIN and COS values
      rcoordSys%dsin_rotation = 0.0_DP
      rcoordSys%dcos_rotation = 1.0_DP
        
    end if
      
    ! Set the scaling factor, if given.
    if (present(dscalingFactor)) then
      rcoordSys%dscalingFactor = dscalingFactor
    else
      rcoordSys%dscalingFactor = 1.0_DP
    endif
  
    ! That's it
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  pure subroutine bgeom_transformPoint2D(rcoordSys, DpointIn, DpointOut)

!<description>
  ! Transforms a given 2D point with coordinates relative to the given
  ! coordinate system to world coordinates.
!</description>

!<input>
  ! The 2D coordinate system that is to be applied to the given point.
  type(t_coordinateSystem2D),  intent(IN)  :: rcoordSys

  ! The 2D point that is to be transformed.
  real(DP), dimension(:),      intent(IN)  :: DpointIn
!</input>

!<output>
  ! The transformed 2D point.
  real(DP), dimension(:),      intent(OUT) :: DpointOut
!</output>

!</subroutine>

    ! Check if the rotation of the coordinate system is non-zero.
    if (rcoordSys%drotation .ne. 0.0_DP) then

      ! Apply rotation to the input coordinates.
      DpointOut(1) = (rcoordSys%dcos_rotation * DpointIn(1)) - &
                     (rcoordSys%dsin_rotation * DpointIn(2))
      DpointOut(2) = (rcoordSys%dsin_rotation * DpointIn(1)) + &
                     (rcoordSys%dcos_rotation * DpointIn(2))

    else

      ! No rotation in the coordinate system, so simply copy the input coords.
      DpointOut(1) = DpointIn(1)
      DpointOut(2) = DpointIn(2)

    end if

    ! Now scale the given coordinates by the scaling factor and translate them
    ! by the origin of our coordinate system.
    DpointOut(1) = (rcoordSys%dscalingFactor * DpointOut(1)) + &
                    rcoordSys%Dorigin(1)
    DpointOut(2) = (rcoordSys%dscalingFactor * DpointOut(2)) + &
                    rcoordSys%Dorigin(2)

    ! That's it!

  end subroutine


  ! ***************************************************************************
  
!<subroutine>

  pure subroutine bgeom_transformBackPoint2D(rcoordSys, DpointIn, DpointOut)

!<description>
  ! Transform a 2D point given in world coordinates to coordinates relative
  ! to a given 2D coordinate system.
  ! This subroutine is the inverse routine of bgeom_transformPoint2D.
!</description>

!<input>
  ! The 2D coordinate system of our input point.
  type(t_coordinateSystem2D), intent(IN)  :: rcoordSys

  ! The 2D point that is to be transformed, relative to the given coordinate
  ! system.
  real(DP), dimension(:),     intent(IN)  :: DpointIn
!</input>

!<output>
  ! The transformed 2D point, in world coordinates.
  real(DP), dimension(:),     intent(OUT) :: DpointOut
!</output>

!</subroutine>

    ! local variables
    real(DP) :: X,Y

    ! Translate the given point by the negatives of our coordinate system
    ! origin.
    X = DpointIn(1) - rcoordSys%Dorigin(1)
    Y = DpointIn(2) - rcoordSys%Dorigin(2)

    ! Now scale the point by the inverse of our coordinate system sclaing
    ! factor, of course only if it is non-zero ...
    if (rcoordSys%dscalingFactor .ne. 0.0_DP) then

      X = X / rcoordSys%dscalingFactor
      Y = Y / rcoordSys%dscalingFactor

    else

      ! ... or set the result to zero and exit subroutine otherwise.
      DpointOut(1) = 0.0_DP
      DpointOut(2) = 0.0_DP
      return

    endif

    ! Check if the rotation of the coordinate system is non-zero.
    if (rcoordSys%drotation .ne. 0.0_DP) then

      ! Apply rotation to the input coordinates.
      DpointOut(1) = ( rcoordSys%dcos_rotation * X) + (rcoordSys%dsin_rotation * Y)
      DpointOut(2) = (-rcoordSys%dsin_rotation * X) + (rcoordSys%dcos_rotation * Y)

    else

      ! No rotation in the coordinate system, so simply copy the input coords.
      DpointOut(1) = X
      DpointOut(2) = Y

    end if

    ! That's it

  end subroutine

end module
