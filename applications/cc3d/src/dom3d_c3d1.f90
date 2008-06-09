!##############################################################################
!# ****************************************************************************
!# <name> dom3d_c3d1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some auxiliary routines for the 3D domain which is
!# used in the "flow around a hexahedron" benchmark.
!#
!# The domain can be described as the [0, 2.5]x[0, 0.41]x[0, 0.41] hexahedron
!# without the [0.45, 0.55]x[0.15, 0.25]x[0, 0.41] hexahedron.
!#
!# To be (mathematically) more precise:
!#
!#   A := { (x,y,z) \in R^3 | 0 <= x <= 2.5 , 0 <= y, z <= 0.41 }
!#   B := { (x,y,z) \in R^3 | 0.45 < x < 0.55 , 0.15 < y < 0.25 
!#                          , 0 <= z <= 0.41 }
!#
!# Then the domain is A \ B.
!#
!# The domain is described by seven boundary regions:
!# 1. The bottom face, i.e. all cells whose Z-coordinate is 0.
!# 2. The front face,  i.e. all cells whose Y-coordinate is 0.
!# 3. The right face,  i.e. all cells whose X-coordinate is 2.5.
!# 4. The back face,   i.e. all cells whose Y-coordinate is 0.41.
!# 5. The left face,   i.e. all cells whose X-coordinate is 0.
!# 6. The top face,    i.e. all cells whose Z-coordinate is 0.41.
!# 7. The obstable,    i.e. all cells whose ...
!#
!# The following routines can be found here:
!#
!#  1.) dom3d_c3d1_calcMeshRegion
!#      -> Creates a mesh region describing a set of specific domain faces.
!#
!#  2.) dom3d_c3d1_HitTest
!#      -> A hit-test routine that is used by dom3d_c3d0_calcMeshRegion,
!#         to identify which cells of the triangulation belong to which
!#         domain face.
!#
!#  3.) dom3d_c3d1_calcParProfile
!#      -> Calculates the parabolic profile for a given domain region.
!#
!# </purpose>
!##############################################################################


MODULE dom3d_c3d1

  USE fsystem
  USE triangulation
  USE meshregion

  IMPLICIT NONE

!<constants>

!<constantblock description="Domain parametrisation">
  
  ! X-Coordinate ranges
  REAL(DP), PARAMETER :: DOM3D_C3D1_X_MIN  = 0.0_DP
  REAL(DP), PARAMETER :: DOM3D_C3D1_X_MAX  = 2.5_DP
  
  ! Y-Coordinate ranges
  REAL(DP), PARAMETER :: DOM3D_C3D1_Y_MIN  = 0.0_DP
  REAL(DP), PARAMETER :: DOM3D_C3D1_Y_MAX  = 0.41_DP
  
  ! Z-Coordinate ranges
  REAL(DP), PARAMETER :: DOM3D_C3D1_Z_MIN  = 0.0_DP
  REAL(DP), PARAMETER :: DOM3D_C3D1_Z_MAX  = 0.41_DP
  
  ! Obstacle midpoint coordinates
  REAL(DP), PARAMETER :: DOM3D_C3D1_X_MID  = 0.5_DP
  REAL(DP), PARAMETER :: DOM3D_C3D1_Y_MID  = 0.2_DP
  
  ! Edge length of the hexahedral obstacle
  REAL(DP), PARAMETER :: DOM3D_C3D1_EDGE   = 0.05_DP

!</constantblock>

!<constantblock description="Region identifiers">

  ! Identification flag for the bottom face:
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_BOTTOM   = 2**0

  ! Identification flag for the front face:
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_FRONT    = 2**1

  ! Identification flag for the right face:
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_RIGHT    = 2**2

  ! Identification flag for the back face:
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_BACK     = 2**3

  ! Identification flag for the left face:
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_LEFT     = 2**4

  ! Identification flag for the top face:
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_TOP      = 2**5
  
  ! Identification flag for the cylinder:
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_OBSTACLE = 2**6
  
  ! Frequently used: All regions except for the right face.
  ! This combination is often used in (Navier-)Stokes examples where
  ! every region (including the cylinder) is Dirichlet, except for the
  ! right face, which is left Neumann.
  INTEGER(I32), PARAMETER :: DOM3D_C3D1_REG_STOKES = DOM3D_C3D1_REG_BOTTOM &
                     + DOM3D_C3D1_REG_FRONT + DOM3D_C3D1_REG_BACK &
                     + DOM3D_C3D1_REG_LEFT + DOM3D_C3D1_REG_TOP &
                     + DOM3D_C3D1_REG_OBSTACLE

!</constantblock>

!</constants>

INTERFACE dom3d_c3d1_calcMeshRegion
  MODULE PROCEDURE dom3d_c3d1_calcMeshRegion_C
  MODULE PROCEDURE dom3d_c3d1_calcMeshRegion_I
END INTERFACE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE dom3d_c3d1_calcMeshRegion_C(rmeshRegion,rtriangulation,cregions)

!<description>
  ! This routine calculates a mesh region containing all vertices, edges and
  ! faces that belong to a specific set of domain regions.
!</description>

!<input>
  ! The triangulation on which the mesh region is to be defined.
  ! This trinagulation must represent the "flow around a cylinder" domain.
  TYPE(t_triangulation), TARGET, INTENT(IN)      :: rtriangulation
  
  ! A combination of DOM3D_C3D1_REG_XXXX contants defined above specifying
  ! which regions should be in the mesh region.
  INTEGER, INTENT(IN)                            :: cregions
!<input>

!<output>
  ! The mesh region that is to be created.
  TYPE(t_meshRegion), INTENT(OUT)                :: rmeshRegion
!<output>

!<subroutine>

  ! Some local variables
  INTEGER :: i
  INTEGER, DIMENSION(7) :: Iregions
  
    ! First of all, set up the Iregions array:
    DO i = 1, 7
      
      ! If the corresponding bit is set, we'll add the region, otherwise
      ! we will set the corresponding Iregions entry to -1.
      IF(IAND(cregions, 2**(i-1)) .NE. 0) THEN
        ! Add the region
        Iregions(i) = i
      ELSE
        ! Ignore the region
        Iregions(i) = -1
      END IF
      
    END DO

    ! Create a mesh-region holding all the regions, based on the hit-test
    ! routine below.
    CALL mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
         MSHREG_IDX_FACE, .TRUE., dom3d_c3d1_HitTest, Iregions)
    
    ! Now calculate the vertice and edge index arrays in the mesh region
    ! based on the face index array.
    CALL mshreg_recalcVerticesFromFaces(rmeshRegion)
    CALL mshreg_recalcEdgesFromFaces(rmeshRegion)
    
    ! That's it

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE dom3d_c3d1_calcMeshRegion_I(rmeshRegion,rtriangulation,Iregions)

!<description>
  ! This routine calculates a mesh region containing all vertices, edges and
  ! faces that belong to a specific set of domain regions.
!</description>

!<input>
  ! The triangulation on which the mesh region is to be defined.
  ! This trinagulation must represent the "flow around a cylinder" domain.
  TYPE(t_triangulation), TARGET, INTENT(IN)      :: rtriangulation
  
  ! An array holding the indices of the domain regions which should be in
  ! the mesh region.
  INTEGER, DIMENSION(:), INTENT(IN)              :: Iregions
!<input>

!<output>
  ! The mesh region that is to be created.
  TYPE(t_meshRegion), INTENT(OUT)                :: rmeshRegion
!<output>

!<subroutine>

    ! Create a mesh-region holding all the regions, based on the hit-test
    ! routine below.
    CALL mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
         MSHREG_IDX_FACE, .TRUE., dom3d_c3d1_HitTest, Iregions)
    
    ! Now calculate the vertice and edge index arrays in the mesh region
    ! based on the face index array.
    CALL mshreg_recalcVerticesFromFaces(rmeshRegion)
    CALL mshreg_recalcEdgesFromFaces(rmeshRegion)
    
    ! That's it

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE dom3d_c3d1_HitTest(inumCells,Dcoords,Ihit,rcollection)
  
  USE collection
  
!<description>
  ! This subroutine is called for hit-testing cells to decide whether a
  ! cell belongs to a specific mesh region or not.
  ! This routine corresponds to the interface defined in
  ! "intf_mshreghittest.inc".
!</description>
  
!<input>
  ! Number of cells (i.e. vertices,edges,faces,etc.) for which the
  ! hit-test is to be performed.
  INTEGER, INTENT(IN)                          :: inumCells
    
  ! Coordinates of the points for which the hit-test is to be performed.
  ! The dimension of the array is at least (1:idim,1:inumCells), where:
  ! -> idim is the dimension of the mesh, i.e. 1 for 1D, 2 for 2D, etc.
  ! -> inumCells is the parameter passed to this routine.
  REAL(DP), DIMENSION(:,:), INTENT(IN)         :: Dcoords
!</input>

!<output>
  ! An array that recieves the result of the hit-test.
  ! The dimension of the array is at least (1:inumCells).
  INTEGER, DIMENSION(:), INTENT(OUT)           :: Ihit
!</output>

!</inputoutput>
  ! OPTIONAL: A collection structure to provide additional information
  ! to the hit-test routine.
  TYPE(t_collection), INTENT(INOUT), OPTIONAL  :: rcollection
!</inputoutput>

!</subroutine>

  ! Some local variables
  INTEGER :: i
  REAL(DP) :: dl
  
  REAL(DP), PARAMETER :: tol = 0.0001_DP
  
    ! Loop through all cells
    DO i = 1, inumCells
      
      ! Bottom face?
      IF      ((Dcoords(3,i) - DOM3D_C3D1_Z_MIN) .LT. tol) THEN
        Ihit(i) = 1
      
      ! Front face?
      ELSE IF ((Dcoords(2,i) - DOM3D_C3D1_Y_MIN) .LT. tol) THEN
        Ihit(i) = 2
      
      ! Right face?
      ELSE IF ((Dcoords(1,i) - DOM3D_C3D1_X_MAX) .GT. -tol) THEN
        Ihit(i) = 3
      
      ! Back face?
      ELSE IF ((Dcoords(2,i) - DOM3D_C3D1_Y_MAX) .GT. -tol) THEN
        Ihit(i) = 4
      
      ! Left face?
      ELSE IF ((Dcoords(1,i) - DOM3D_C3D1_X_MIN) .LT. tol) THEN
        Ihit(i) = 5
      
      ! Top face?
      ELSE IF ((Dcoords(3,i) - DOM3D_C3D1_Z_MAX) .GT. -tol) THEN
        Ihit(i) = 6
      
      ! Hexahedral obstacle?
      ELSE IF ((ABS(Dcoords(1,i) - DOM3D_C3D1_X_MID) &
                  - DOM3D_C3D1_EDGE .LT. tol) .AND. &
               (ABS(Dcoords(2,i) - DOM3D_C3D1_Y_MID) &
                  - DOM3D_C3D1_EDGE .LT. tol)) THEN
        Ihit(i) = 7

      ELSE
        Ihit(i) = 0
      END IF
      
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE dom3d_c3d1_calcParProfile (dvalue, iregion, dx, dy, dz)
  
!<description>
  ! This routine calculates the parabolic profile for a specified domain
  ! boundary region.
!</description>

!<input>
  ! Specifies the index of the boundary region of which the parabolic profile
  ! is to be calculated.
  INTEGER, INTENT(IN)                            :: iregion
  
  ! The coordinates of the point in which the parabolic profile is to be
  ! evaluated.
  REAL(DP), INTENT(IN)                           :: dx, dy, dz
!</input>

!<output>
  ! The result of the evaluation of the parabolic profile.
  REAL(DP), INTENT(OUT)                          :: dvalue
!</output>

!</subroutine>

  ! Relative X-, Y- and Z-coordinates
  REAL(DP) :: x,y,z
    
    ! Convert coordinates to be in range [0,1]
    x = (dx - DOM3D_C3D1_X_MIN) / (DOM3D_C3D1_X_MAX - DOM3D_C3D1_X_MIN)
    y = (dy - DOM3D_C3D1_Y_MIN) / (DOM3D_C3D1_Y_MAX - DOM3D_C3D1_Y_MIN)
    z = (dz - DOM3D_C3D1_Z_MIN) / (DOM3D_C3D1_Z_MAX - DOM3D_C3D1_Z_MIN)

    SELECT CASE(iregion)
    CASE (1,6)
      ! X-Y-Profile
      dvalue = 16.0_DP*x*(1.0_DP-x)*y*(1.0_DP-y)
    
    CASE (2,4)
      ! X-Z-Profile
      dvalue = 16.0_DP*x*(1.0_DP-x)*z*(1.0_DP-z)
    
    CASE (3,5)
      ! Y-Z-Profile
      dvalue = 16.0_DP*y*(1.0_DP-y)*z*(1.0_DP-z)
    
    CASE (7)
      ! Z-Profile on the boundary of the hexahedral obstacle
      dvalue = 2.0_DP*z*(1.0_DP-z)

    CASE DEFAULT
      ! Invalid region
      PRINT *, 'ERROR: dom3d_c3d0_calcParProfile: Invalid region'
      CALL sys_halt()
    END SELECT
    
  END SUBROUTINE

END MODULE
