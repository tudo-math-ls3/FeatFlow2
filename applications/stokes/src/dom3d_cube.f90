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
!# </purpose>
!##############################################################################

module dom3d_cube

  use fsystem
  use triangulation
  use meshregion

  implicit none

!<constants>

!<constantblock description="Face identifiers">

  ! Identification flag for the bottom face:
  integer(I32), parameter :: DOM3D_CUBE_REG_BOTTOM = 2**0

  ! Identification flag for the front face:
  integer(I32), parameter :: DOM3D_CUBE_REG_FRONT  = 2**1

  ! Identification flag for the right face:
  integer(I32), parameter :: DOM3D_CUBE_REG_RIGHT  = 2**2

  ! Identification flag for the back face:
  integer(I32), parameter :: DOM3D_CUBE_REG_BACK   = 2**3

  ! Identification flag for the left face:
  integer(I32), parameter :: DOM3D_CUBE_REG_LEFT   = 2**4

  ! Identification flag for the top face:
  integer(I32), parameter :: DOM3D_CUBE_REG_TOP    = 2**5
  
  ! Frequently used: All faces except for the right one.
  ! This combination is often used in (Navier-)Stokes examples where
  ! every face is Dirichlet, except for the right one, which is left
  ! Neumann.
  integer(I32), parameter :: DOM3D_CUBE_REG_STOKES = DOM3D_CUBE_REG_BOTTOM &
                     + DOM3D_CUBE_REG_FRONT + DOM3D_CUBE_REG_BACK&
                     + DOM3D_CUBE_REG_LEFT + DOM3D_CUBE_REG_TOP

!</constantblock>

!</constants>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine dom3d_cube_calcMeshRegion(rmeshRegion,rtriangulation,cfaces)

!<description>
  ! This routine calculates a mesh region containing all vertices, edges and
  ! faces that belong to a specific set of domain faces.
!</description>

!<input>
  ! The triangulation on which the mesh region is to be defined.
  ! This trinagulation must represent the 3D [0,1]^3 cube domain.
  type(t_triangulation), target, intent(in)      :: rtriangulation
  
  ! A combination of DOM3D_CUBE_REG_XXXX contants defined above specifying
  ! which faces should be in the mesh region.
  integer(I32), intent(in)                       :: cfaces
!<input>

!<output>
  ! The mesh region that is to be created.
  type(t_meshRegion), intent(out)                :: rmeshRegion
!<output>

!<subroutine>

  ! Some local variables
  integer :: i
  integer, dimension(6) :: Ifaces
  
    ! First of all, set up the Ifaces array:
    do i = 1, 6
      
      ! If the corresponding bit is set, we will add the face, otherwise
      ! we will set the corresponding Ifaces entry to -1.
      if(iand(cfaces, int(2**(i-1),I32)) .ne. 0) then
        ! Add the face
        Ifaces(i) = i
      else
        ! Ignore the face
        Ifaces(i) = -1
      end if
      
    end do

    ! Create a mesh-region holding all the faces, based on the hit-test
    ! routine below.
    call mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
         MSHREG_IDX_FACE, .true., dom3d_cube_HitTest, Ifaces)
    
    ! Now calculate the vertice and edge index arrays in the mesh region
    ! based on the face index array.
    call mshreg_recalcVerticesFromFaces(rmeshRegion)
    call mshreg_recalcEdgesFromFaces(rmeshRegion)
    
    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine dom3d_cube_HitTest(inumCells,Dcoords,Ihit,rcollection)
  
  use collection
  
!<description>
  ! This subroutine is called for hit-testing cells to decide whether a
  ! cell belongs to a specific mesh region or not.
  ! This routine corresponds to the interface defined in
  ! "intf_mshreghittest.inc".
!</description>
  
!<input>
  ! Number of cells (i.e. vertices,edges,faces,etc.) for which the
  ! hit-test is to be performed.
  integer, intent(in)                          :: inumCells
    
  ! Coordinates of the points for which the hit-test is to be performed.
  ! The dimension of the array is at least (1:idim,1:inumCells), where:
  ! -> idim is the dimension of the mesh, i.e. 1 for 1D, 2 for 2D, etc.
  ! -> inumCells is the parameter passed to this routine.
  real(DP), dimension(:,:), intent(in)         :: Dcoords
!</input>

!<output>
  ! An array that recieves the result of the hit-test.
  ! The dimension of the array is at least (1:inumCells).
  integer, dimension(:), intent(out)           :: Ihit
!</output>

!</inputoutput>
  ! OPTIONAL: A collection structure to provide additional information
  ! to the hit-test routine.
  type(t_collection), intent(inout), optional  :: rcollection
!</inputoutput>

!</subroutine>

  integer :: i
  
  ! We will not check the coordinates to be <= 0 (or >= 1), but we will
  ! check them against a small tolerance to avoid that boundary cells
  ! are treated as "inside the domain" due to rounding errors.
  real(DP), parameter :: dzero = 0.0001_DP
  real(DP), parameter :: done = 0.9999_DP
  
    ! Loop through all cells
    do i = 1, inumCells
      
      ! Bottom face?
      if      (Dcoords(3,i) .lt. dzero) then
        Ihit(i) = 1
      
      ! Front face?
      else if (Dcoords(2,i) .lt. dzero) then
        Ihit(i) = 2
      
      ! Right face?
      else if (Dcoords(1,i) .gt. done) then
        Ihit(i) = 3
      
      ! Back face?
      else if (Dcoords(2,i) .gt. done) then
        Ihit(i) = 4
      
      ! Left face?
      else if (Dcoords(1,i) .lt. dzero) then
        Ihit(i) = 5
      
      ! Top face?
      else if (Dcoords(3,i) .gt. done) then
        Ihit(i) = 6
      
      ! Inner face?
      else
        Ihit(i) = 0
      end if
      
    end do

  end subroutine

end module
