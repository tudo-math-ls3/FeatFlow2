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
!#  1.) bgeom_initCoordSys2D, bgeom_initCoordSys3D
!#      -> Creates a 2D/3D coordinate system from an optionally given origin,
!#         rotation and scaling factor.
!#
!#  2.) bgeom_transformPoint2D
!#      -> Transforms a 2D point with coordinates relative to a given
!#         coordinate system to world coordinates.
!#
!#  3.) bgeom_transformBackPoint2D
!#      -> Transforms a 2D point with world coordinates to coordinates relative
!#         to a given coordinate system.
!#
!#  4.) bgeom_transformPoint3D
!#      -> Transforms a 3D point with coordinates relative to a given
!#         coordinate system to world coordinates.
!#
!#  5.) bgeom_transformBackPoint3D
!#      -> Transforms a 3D point with world coordinates to coordinates relative
!#         to a given coordinate system.

!# </purpose>
!##############################################################################

module basicgeometry

!$use omp_lib
  use fsystem

  implicit none
  
  private

!<constants>

!<constantblock description="Dimension constants">
  
  ! Dimension constant for 1D triangulations.
  integer, parameter, public :: NDIM1D = 1

  ! Dimension constant for 2D triangulations.
  integer, parameter, public :: NDIM2D = 2

  ! Dimension constant for 3D triangulations.
  integer, parameter, public :: NDIM3D = 3

!</constantblock>

!<constantblock description="Element shape identifiers">
  ! Unknown shape
  integer(I32), parameter, public :: BGEOM_SHAPE_UNKNOWN = 0_I32
  
  ! Line shape (1D)
  integer(I32), parameter, public :: BGEOM_SHAPE_LINE  = 1000002_I32
  
  ! Triangle shape (2D)
  integer(I32), parameter, public :: BGEOM_SHAPE_TRIA  = 2000303_I32
  
  ! Quadrilateral shape (2D)
  integer(I32), parameter, public :: BGEOM_SHAPE_QUAD  = 2000404_I32
  
  ! Tetrahedron shape (3D)
  integer(I32), parameter, public :: BGEOM_SHAPE_TETRA = 3040604_I32
  
  ! Hexahedron shape (3D)
  integer(I32), parameter, public :: BGEOM_SHAPE_HEXA  = 3061208_I32
  
  ! Pyramid shape (3D)
  integer(I32), parameter, public :: BGEOM_SHAPE_PYRA  = 3050805_I32
  
  ! Prism shape (3D)
  integer(I32), parameter, public :: BGEOM_SHAPE_PRISM = 3050906_I32

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
  
  public :: t_point2D

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
  
  public :: t_point3D

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

  public :: t_coordinateSystem2D

  !</typeblock>
  
  !<typeblock>
  
  ! The point structure for regular 3D coordinate systems, consisting of
  ! a rotated X-, Y-axes and Z-axes
  type t_coordinateSystem3D
  
    ! Coordinates of the origin
    real(DP), dimension(3) :: Dorigin = (/0.0_dp,0.0_dp,0.0_dp/)
    
    ! Rotation angle; 0..2*PI
    real(DP) :: drotationX = 0.0_dp
    
    real(DP) :: drotationY = 0.0_dp
    
    real(DP) :: drotationZ = 0.0_dp
    
    ! scaling factor of the coordinate system; usually = 1.0
    real(DP) :: dscalingFactor = 1.0_dp
    
  end type
  
    public :: t_coordinateSystem3D

  !</typeblock>
  
  !</types>

  public :: bgeom_initCoordSys2D
  public :: bgeom_initCoordSys3D
  public :: bgeom_transformPoint2D
  public :: bgeom_transformBackPoint2D
  public :: bgeom_transformPoint3D
  public :: bgeom_transformBackPoint3D
  
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
  real(DP), dimension(:), optional,    intent(in)  :: Dorigin
  
  ! OPTIONAL: The rotation angle of the coordinate system.
  ! Angle range: 0..2*PI
  ! Set to 0.0_DP if not given.
  real(DP), optional,                  intent(in)  :: drotation
  
  ! The scaling factor. Should be != 0.0_DP
  ! Set to 1.0_DP if not given.
  real(DP), optional,                  intent(in)  :: dscalingFactor

!</input>

!<output>
  ! A t_coordinateSystem2D structure to be written.
  type(t_coordinateSystem2D),          intent(out) :: rcoordSys

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
  
    ! That is it
  
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
  type(t_coordinateSystem2D),  intent(in)  :: rcoordSys

  ! The 2D point that is to be transformed.
  real(DP), dimension(:),      intent(in)  :: DpointIn
!</input>

!<output>
  ! The transformed 2D point.
  real(DP), dimension(:),      intent(out) :: DpointOut
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

    ! That is it!

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
  type(t_coordinateSystem2D), intent(in)  :: rcoordSys

  ! The 2D point that is to be transformed, relative to the given coordinate
  ! system.
  real(DP), dimension(:),     intent(in)  :: DpointIn
!</input>

!<output>
  ! The transformed 2D point, in world coordinates.
  real(DP), dimension(:),     intent(out) :: DpointOut
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

    ! That is it

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  pure subroutine bgeom_initCoordSys3D(rcoordSys, Dorigin, drotationX, drotationY,&
                                       drotationZ, dscalingFactor)

!<description>
  ! Creates a 3D coordinate system from an optionally given origin, rotation angles
  ! and scaling factor.
!</description>

!<input>
  ! OPTIONAL: The origin of the coordinate system.
  ! Set to (/ 0.0_dp, 0.0_dp, 0.0_dp /) if not given.
  real(DP), dimension(:), optional,    intent(in)  :: Dorigin
  
  ! OPTIONAL: The rotation angles of the coordinate system.
  ! Angle range: 0..2*PI
  ! Set to 0.0_DP if not given.
  real(DP), optional,                  intent(in)  :: drotationX
  real(DP), optional,                  intent(in)  :: drotationY
  real(DP), optional,                  intent(in)  :: drotationZ
  
  ! The scaling factor. Should be != 0.0_DP
  ! Set to 1.0_DP if not given.
  real(DP), optional,                  intent(in)  :: dscalingFactor

!</input>

!<output>
  ! A t_coordinateSystem3D structure to be written.
  type(t_coordinateSystem3D),          intent(out) :: rcoordSys

!</output>

!</subroutine>

    ! Set the origin, if given.
    if (present(Dorigin)) then
      rcoordSys%Dorigin = (/ Dorigin(1), Dorigin(2), Dorigin(3) /)
    else
      rcoordSys%Dorigin = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    end if
      
    ! Set the rotation, if given.
    if (present(drotationX)) then
      rcoordSys%drotationX = drotationX
    else
      rcoordSys%drotationX = 0.0_DP
    end if

    ! Set the rotation, if given.
    if (present(drotationY)) then
      rcoordSys%drotationY = drotationY
    else
      rcoordSys%drotationY = 0.0_DP
    end if

    ! Set the rotation, if given.
    if (present(drotationZ)) then
      rcoordSys%drotationZ = drotationZ
    else
      rcoordSys%drotationZ = 0.0_DP
    end if
      
    ! Set the scaling factor, if given.
    if (present(dscalingFactor)) then
      rcoordSys%dscalingFactor = dscalingFactor
    else
      rcoordSys%dscalingFactor = 1.0_DP
    endif
  
    ! That is it
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  pure subroutine bgeom_transformPoint3D(rcoordSys, DpointIn, DpointOut)

!<description>
  ! Transforms a 3D point in model coordinates to world coordinates
  ! This routine may be inefficient if a large number of of points is to
  ! be transformed. Then a matrix based transformation would be more efficient.
  ! But this would require a rewriting of the mechanism how routines, that call
  ! this function(i.e. geom_isingeometry) work.
!</description>

!<input>
  ! The 3D coordinate system that is to be applied to the given point.
  type(t_coordinateSystem3D),  intent(in)  :: rcoordSys

  ! The 3D point that is to be transformed.
  real(DP), dimension(:),      intent(in)  :: DpointIn
!</input>

!<output>
  ! The transformed 3D point.
  real(DP), dimension(:),      intent(out) :: DpointOut
!</output>

!</subroutine>

  real(dp) :: sx,sy,sz,cx,cy,cz,x1,y1,z1,x2,y2,z2
  
  ! Check if the rotation of the coordinate system is non-zero.
  if (rcoordSys%drotationX .ne. 0.0_DP) then

    ! setup
    sx = sin(rcoordSys%drotationX)
    cx = cos(rcoordSys%drotationX)
    sy = sin(rcoordSys%drotationY)
    cy = cos(rcoordSys%drotationY)
    sz = sin(rcoordSys%drotationZ)
    cz = cos(rcoordSys%drotationZ)

    x1 =  DpointIn(1) * cz +  DpointIn(2) * sz
    y1 =  DpointIn(2) * cz -  DpointIn(1) * sz
    z1 =  DpointIn(3)
    x2 = x1 * cy + z1 * sy
    y2 = z1
    z2 = z1 * cy - x1 * sy
    DpointOut(1) = x2
    DpointOut(2) = y2 * cx + z1 * sx
    DpointOut(3) = z2 * cx - x1 * sx

  else

    ! No rotation in the coordinate system, so simply copy the input coords.
    DpointOut(1) = DpointIn(1)
    DpointOut(2) = DpointIn(2)
    DpointOut(3) = DpointIn(3)

  end if

  ! Now scale the given coordinates by the scaling factor and translate them
  ! by the origin of our coordinate system.
  DpointOut(1) = (rcoordSys%dscalingFactor * DpointOut(1)) + &
                  rcoordSys%Dorigin(1)
  DpointOut(2) = (rcoordSys%dscalingFactor * DpointOut(2)) + &
                  rcoordSys%Dorigin(2)
  DpointOut(3) = (rcoordSys%dscalingFactor * DpointOut(3)) + &
                  rcoordSys%Dorigin(3)

  ! That is it!

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  pure subroutine bgeom_transformBackPoint3D(rcoordSys, DpointIn, DpointOut)

!<description>
  ! Transform a 3d point from world coordinates to model coordinates
!</description>

!<input>
  ! The 3D coordinate system of our input point.
  type(t_coordinateSystem3D), intent(in)  :: rcoordSys

  ! The 3D point that is to be transformed, relative to the given coordinate
  ! system.
  real(DP), dimension(:),     intent(in)  :: DpointIn
!</input>

!<output>
  ! The transformed 3D point, in world coordinates.
  real(DP), dimension(:),     intent(out) :: DpointOut
!</output>

!</subroutine>

  ! local variables
  real(DP) :: X,Y,Z
    
  real(dp) :: sx,sy,sz,cx,cy,cz,x1,y1,z1,x2,y2,z2
  
  ! Translate the given point by the negatives of our coordinate system
  ! origin.
  X = DpointIn(1) - rcoordSys%Dorigin(1)
  Y = DpointIn(2) - rcoordSys%Dorigin(2)
  Z = DpointIn(3) - rcoordSys%Dorigin(3)

  ! Now scale the point by the inverse of our coordinate system sclaing
  ! factor, of course only if it is non-zero ...
  if (rcoordSys%dscalingFactor .ne. 0.0_DP) then

    X = X / rcoordSys%dscalingFactor
    Y = Y / rcoordSys%dscalingFactor
    Z = Z / rcoordSys%dscalingFactor

  else

    ! ... or set the result to zero and exit subroutine otherwise.
    DpointOut(1) = 0.0_DP
    DpointOut(2) = 0.0_DP
    DpointOut(3) = 0.0_DP
    return

  endif
  
  ! Check if the rotation of the coordinate system is non-zero.
  if (rcoordSys%drotationX .ne. 0.0_DP) then

    ! setup
    sx = sin(-rcoordSys%drotationX)
    cx = cos(-rcoordSys%drotationX)
    sy = sin(-rcoordSys%drotationY)
    cy = cos(-rcoordSys%drotationY)
    sz = sin(-rcoordSys%drotationZ)
    cz = cos(-rcoordSys%drotationZ)

    x1 =  X * cz + Y * sz
    y1 =  Y * cz - X * sz
    z1 =  Z
    x2 = x1 * cy + z1 * sy
    y2 = z1
    z2 = z1 * cy - x1 * sy
    DpointOut(1) = x2
    DpointOut(2) = y2 * cx + z1 * sx
    DpointOut(3) = z2 * cx - x1 * sx

  else

    ! No rotation in the coordinate system, so simply copy the input coords.
    DpointOut(1) = X
    DpointOut(2) = Y
    DpointOut(3) = Z

  end if

  end subroutine
  
  ! ***************************************************************************

end module
