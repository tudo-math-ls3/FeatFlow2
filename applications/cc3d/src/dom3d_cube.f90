!##############################################################################
!# ****************************************************************************
!# <name> dom3d_cube </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some auxiliary routines for the 3D [0,1]x[0,1]x[0,1]
!# cube domain and its derived triangulations.
!#
!# The domain is described by six faces:
!# 1. The bottom face, i.e. all cells whose Z-coordinate is 0.
!# 2. The front face,  i.e. all cells whose Y-coordinate is 0.
!# 3. The right face,  i.e. all cells whose X-coordinate is 1.
!# 4. The back face,   i.e. all cells whose Y-coordinate is 1.
!# 5. The left face,   i.e. all cells whose X-coordinate is 0.
!# 6. The top face,    i.e. all cells whose Z-coordinate is 1.
!#
!# The following routines can be found here:
!#
!#  1.) dom3d_cube_calcMeshRegion
!#      -> Creates a mesh region describing a set of specific domain faces.
!#
!#  2.) dom3d_cube_HitTest
!#      -> A hit-test routine that is used by dom3d_cube_calcMeshRegion,
!#         to identify which cells of the triangulation belong to which
!#         domain face.
!#
!#  3.) dom3d_cube_calcParProfile
!#      -> Calculates the parabolic profile for a given domain region.
!#
!# </purpose>
!##############################################################################

MODULE dom3d_cube

  USE fsystem
  USE paramlist
  USE triangulation
  USE meshregion

  IMPLICIT NONE

!<constants>

!<constantblock description="Face identifiers">

  ! Identification flag for the bottom face:
  INTEGER(I32), PARAMETER :: DOM3D_CUBE_REG_BOTTOM = 2**0

  ! Identification flag for the front face:
  INTEGER(I32), PARAMETER :: DOM3D_CUBE_REG_FRONT  = 2**1

  ! Identification flag for the right face:
  INTEGER(I32), PARAMETER :: DOM3D_CUBE_REG_RIGHT  = 2**2

  ! Identification flag for the back face:
  INTEGER(I32), PARAMETER :: DOM3D_CUBE_REG_BACK   = 2**3

  ! Identification flag for the left face:
  INTEGER(I32), PARAMETER :: DOM3D_CUBE_REG_LEFT   = 2**4

  ! Identification flag for the top face:
  INTEGER(I32), PARAMETER :: DOM3D_CUBE_REG_TOP    = 2**5
  
  ! Frequently used: All faces except for the right one.
  ! This combination is often used in (Navier-)Stokes examples where
  ! every face is Dirichlet, except for the right one, which is left
  ! Neumann.
  INTEGER(I32), PARAMETER :: DOM3D_CUBE_REG_STOKES = DOM3D_CUBE_REG_BOTTOM &
                     + DOM3D_CUBE_REG_FRONT + DOM3D_CUBE_REG_BACK&
                     + DOM3D_CUBE_REG_LEFT + DOM3D_CUBE_REG_TOP

!</constantblock>

!</constants>

INTERFACE dom3d_cube_calcMeshRegion
  MODULE PROCEDURE dom3d_cube_calcMeshRegion_C
  MODULE PROCEDURE dom3d_cube_calcMeshRegion_I
END INTERFACE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE dom3d_cube_calcMeshRegion_C(rmeshRegion,rtriangulation,cfaces)

!<description>
  ! This routine calculates a mesh region containing all vertices, edges and
  ! faces that belong to a specific set of domain faces.
!</description>

!<input>
  ! The triangulation on which the mesh region is to be defined.
  ! This trinagulation must represent the 3D [0,1]^3 cube domain.
  TYPE(t_triangulation), TARGET, INTENT(IN)      :: rtriangulation
  
  ! A combination of DOM3D_CUBE_REG_XXXX contants defined above specifying
  ! which faces should be in the mesh region.
  INTEGER, INTENT(IN)                            :: cfaces
!<input>

!<output>
  ! The mesh region that is to be created.
  TYPE(t_meshRegion), INTENT(OUT)                :: rmeshRegion
!<output>

!<subroutine>

  ! Some local variables
  INTEGER :: i
  INTEGER, DIMENSION(6) :: Ifaces
  
    ! First of all, set up the Ifaces array:
    DO i = 1, 6
      
      ! If the corresponding bit is set, we'll add the face, otherwise
      ! we will set the corresponding Ifaces entry to -1.
      IF(IAND(cfaces, 2**(i-1)) .NE. 0) THEN
        ! Add the face
        Ifaces(i) = i
      ELSE
        ! Ignore the face
        Ifaces(i) = -1
      END IF
      
    END DO

    ! Create a mesh-region holding all the faces, based on the hit-test
    ! routine below.
    CALL mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
         MSHREG_IDX_FACE, .TRUE., dom3d_cube_HitTest, Ifaces)
    
    ! Now calculate the vertice and edge index arrays in the mesh region
    ! based on the face index array.
    CALL mshreg_recalcVerticesFromFaces(rmeshRegion)
    CALL mshreg_recalcEdgesFromFaces(rmeshRegion)
    
    ! That's it

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE dom3d_cube_calcMeshRegion_I(rmeshRegion,rtriangulation,Ifaces)

!<description>
  ! This routine calculates a mesh region containing all vertices, edges and
  ! faces that belong to a specific set of domain faces.
!</description>

!<input>
  ! The triangulation on which the mesh region is to be defined.
  ! This trinagulation must represent the 3D [0,1]^3 cube domain.
  TYPE(t_triangulation), TARGET, INTENT(IN)      :: rtriangulation
  
  ! An array holding the indices of the faces which should be in the mesh
  ! region.
  INTEGER, DIMENSION(:), INTENT(IN)              :: Ifaces
!<input>

!<output>
  ! The mesh region that is to be created.
  TYPE(t_meshRegion), INTENT(OUT)                :: rmeshRegion
!<output>

!<subroutine>

    ! Create a mesh-region holding all the faces, based on the hit-test
    ! routine below.
    CALL mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
         MSHREG_IDX_FACE, .TRUE., dom3d_cube_HitTest, Ifaces)
    
    ! Now calculate the vertice and edge index arrays in the mesh region
    ! based on the face index array.
    CALL mshreg_recalcVerticesFromFaces(rmeshRegion)
    CALL mshreg_recalcEdgesFromFaces(rmeshRegion)
    
    ! That's it

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE dom3d_cube_HitTest(inumCells,Dcoords,Ihit,rcollection)
  
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

  INTEGER :: i
  
  ! We will not check the coordinates to be <= 0 (or >= 1), but we will
  ! check them against a small tolerance to avoid that boundary cells
  ! are treated as "inside the domain" due to rounding errors.
  REAL(DP), PARAMETER :: dzero = 0.0001_DP
  REAL(DP), PARAMETER :: done = 0.9999_DP
  
    ! Loop through all cells
    DO i = 1, inumCells
      
      ! Bottom face?
      IF      (Dcoords(3,i) .LT. dzero) THEN
        Ihit(i) = 1
      
      ! Front face?
      ELSE IF (Dcoords(2,i) .LT. dzero) THEN
        Ihit(i) = 2
      
      ! Right face?
      ELSE IF (Dcoords(1,i) .GT. done) THEN
        Ihit(i) = 3
      
      ! Back face?
      ELSE IF (Dcoords(2,i) .GT. done) THEN
        Ihit(i) = 4
      
      ! Left face?
      ELSE IF (Dcoords(1,i) .LT. dzero) THEN
        Ihit(i) = 5
      
      ! Top face?
      ELSE IF (Dcoords(3,i) .GT. done) THEN
        Ihit(i) = 6
      
      ! Inner face? 
      ELSE
        Ihit(i) = 0
      END IF
      
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE dom3d_cube_calcParProfile (dvalue, iregion, dx, dy, dz)
  
!<description>
  ! This routine calculates the parabolic profile for a speicified domain
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

    SELECT CASE(iregion)
    CASE (1,6)
      ! X-Y-Profile
      dvalue = 4.0_DP*dx*(1.0_DP-dx)*dy*(1.0_DP-dy)
    
    CASE (2,4)
      ! X-Z-Profile
      dvalue = 4.0_DP*dx*(1.0_DP-dx)*dz*(1.0_DP-dz)
    
    CASE (3,5)
      ! Y-Z-Profile
      dvalue = 4.0_DP*dy*(1.0_DP-dy)*dz*(1.0_DP-dz)
      
    CASE DEFAULT
      ! Invalid region
      PRINT *, 'ERROR: dom3d_cube_calcParProfile: Invalid region'
      CALL sys_halt()
    END SELECT
    
  END SUBROUTINE

END MODULE
