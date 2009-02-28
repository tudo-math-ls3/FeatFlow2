!##############################################################################
!# ****************************************************************************
!# <name> cc3dmini_aux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some auxiliary routines for the cc3d solvers.
!# </purpose>
!##############################################################################

module cc3dmini_aux

  use fsystem
  use storage
  use triangulation
  use meshregion

  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc3daux_calcCubeDirichletRegion(rmeshRegion)

!<description>
  ! This auxiliary routine extracts the faces from a mesh region describing
  ! the standard-cube's boundary which will have Dirichlet boundary
  ! conditions.
  !
  ! In the cc3d examples, we will prescribe Dirichlet boundary conditions
  ! on all boundary faces except the face, where the X-coordinate is 1.
  ! This one face will have Neumann boundary conditions, and therefore it must
  ! be excluded from the mesh region for the Dirichlet boundary conditions.
!</description>

!<inputoutput>
  ! The mesh region describing the cube's boundary.
  ! On entry: the complete cube boundary
  ! On exit: the cube boundary except the Neumann boundary face
  type(t_meshRegion), intent(INOUT)                 :: rmeshRegion

!</inputoutput>
!</subroutine>

   ! A hand full of local variables
   type(t_triangulation), pointer :: p_rtria
   integer, dimension(:), pointer :: p_IfaceIdx
   integer, dimension(:), allocatable :: IfacesInMR
   integer, dimension(:,:), pointer :: p_IvertsAtFace
   real(DP), dimension(:,:), pointer :: p_Dcoords
   integer :: i, iface, iNAT
   real(DP) :: dx

    ! Let's see if the mesh region has faces at all...
    if((rmeshRegion%h_IfaceIdx .eq. ST_NOHANDLE) .or. &
       (rmeshRegion%NAT .le. 0)) then
      print *, 'ERROR: st3daux_calcCubeRegion'
      print *, 'Mesh region does not have any faces!'
      stop
    end if
    
    ! Get the arrays from the triangulation
    p_rtria => rmeshRegion%p_rtriangulation
    call storage_getbase_int2D(p_rtria%h_IverticesAtFace, p_IvertsAtFace)
    call storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dcoords)
    
    ! Get the face indices from the mesh region
    call storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! Allocate a temporary array
    allocate(IfacesInMR(p_rtria%NAT))
    IfacesInMR = 0
    
    ! Now go through all faces in the input mesh region
    do i = 1, rmeshRegion%NAT
    
      ! Get the face index
      iface = p_IfaceIdx(i)
    
      ! Calculate midpoint x-coordinate of the face
      dx = 0.25_DP * (p_Dcoords(1,p_IvertsAtFace(1,iface)) &
                    + p_Dcoords(1,p_IvertsAtFace(2,iface)) &
                    + p_Dcoords(1,p_IvertsAtFace(3,iface)) &
                    + p_Dcoords(1,p_IvertsAtFace(4,iface)))
                    
      ! Do we have to add this face?
      if (dx .lt. 0.999_DP) IfacesInMR(iface) = 1
        
    end do
    
    ! Count how many faces there are in the new mesh region
    iNAT = 0
    do i = 1, ubound(IfacesInMR,1)
      iNAT = iNAT + IfacesInMR(i)
    end do
    
    ! At this point, we can release the old mesh region
    call mshreg_done(rmeshRegion)

    ! Make sure that there are any faces left
    if (iNAT .le. 0) then
      print *, 'ERROR: cc3daux_calcCubeRegion'
      print *, 'Could not find any suitable faces!'
      stop
    end if
    
    ! Call the storage to allocate the face index array
    call storage_new('cc3daux_calcCubeRegion', 'p_IfaceIdx', iNAT, &
                     ST_INT, rmeshRegion%h_IfaceIdx, ST_NEWBLOCK_NOINIT)

    ! Get the face indices from the mesh region
    call storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! And build up the face index array
    iface = 1
    do i = 1, ubound(IfacesInMR,1)
      if (IfacesInMR(i) .ne. 0) then
        p_IfaceIdx(iface) = i
        iface = iface + 1
      end if
    end do
    
    ! Copy the triangulation pointer to the new mesh region
    rmeshRegion%p_rtriangulation => p_rtria
    rmeshRegion%NAT = iNAT
    
    ! Now we can deallocate the work array
    deallocate(IfacesInMR)
    
    ! Finally, recalculate the edge and vertex index arrays from the faces
    call mshreg_recalcEdgesFromFaces(rmeshRegion)
    call mshreg_recalcVerticesFromFaces(rmeshRegion)
    
    ! That's it

  end subroutine

end module
