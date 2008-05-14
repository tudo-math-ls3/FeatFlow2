!##############################################################################
!# ****************************************************************************
!# <name> cc3dmini_aux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some auxiliary routines for the cc3d solvers.
!# </purpose>
!##############################################################################

MODULE cc3dmini_aux

  USE fsystem
  USE storage
  USE triangulation
  USE meshregion

  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc3daux_calcCubeDirichletRegion(rmeshRegion)

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
  TYPE(t_meshRegion), INTENT(INOUT)                 :: rmeshRegion

!</inputoutput>
!</subroutine>

   ! A hand full of local variables
   TYPE(t_triangulation), POINTER :: p_rtria
   INTEGER, DIMENSION(:), POINTER :: p_IfaceIdx
   INTEGER, DIMENSION(:), ALLOCATABLE :: IfacesInMR
   INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertsAtFace
   REAL(DP), DIMENSION(:,:), POINTER :: p_Dcoords
   INTEGER :: i, iface, iNAT
   REAL(DP) :: dx

    ! Let's see if the mesh region has faces at all...
    IF((rmeshRegion%h_IfaceIdx .EQ. ST_NOHANDLE) .OR. &
       (rmeshRegion%NAT .LE. 0)) THEN
      PRINT *, 'ERROR: st3daux_calcCubeRegion'
      PRINT *, 'Mesh region does not have any faces!'
      STOP
    END IF
    
    ! Get the arrays from the triangulation
    p_rtria => rmeshRegion%p_rtriangulation
    CALL storage_getbase_int2D(p_rtria%h_IverticesAtFace, p_IvertsAtFace)
    CALL storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dcoords)
    
    ! Get the face indices from the mesh region
    CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! Allocate a temporary array
    ALLOCATE(IfacesInMR(p_rtria%NAT))
    IfacesInMR = 0
    
    ! Now go through all faces in the input mesh region
    DO i = 1, rmeshRegion%NAT
    
      ! Get the face index
      iface = p_IfaceIdx(i)
    
      ! Calculate midpoint x-coordinate of the face
      dx = 0.25_DP * (p_Dcoords(1,p_IvertsAtFace(1,iface)) &
                    + p_Dcoords(1,p_IvertsAtFace(2,iface)) &
                    + p_Dcoords(1,p_IvertsAtFace(3,iface)) &
                    + p_Dcoords(1,p_IvertsAtFace(4,iface)))
                    
      ! Do we have to add this face?
      IF (dx .LT. 0.999_DP) IfacesInMR(iface) = 1
        
    END DO
    
    ! Count how many faces there are in the new mesh region
    iNAT = 0
    DO i = 1, UBOUND(IfacesInMR,1)
      iNAT = iNAT + IfacesInMR(i)
    END DO
    
    ! At this point, we can release the old mesh region
    CALL mshreg_done(rmeshRegion)

    ! Make sure that there are any faces left
    IF (iNAT .LE. 0) THEN
      PRINT *, 'ERROR: cc3daux_calcCubeRegion'
      PRINT *, 'Could not find any suitable faces!'
      STOP
    END IF
    
    ! Call the storage to allocate the face index array
    CALL storage_new('cc3daux_calcCubeRegion', 'p_IfaceIdx', iNAT, &
                     ST_INT, rmeshRegion%h_IfaceIdx, ST_NEWBLOCK_NOINIT)

    ! Get the face indices from the mesh region
    CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! And build up the face index array
    iface = 1
    DO i = 1, UBOUND(IfacesInMR,1)
      IF (IfacesInMR(i) .NE. 0) THEN
        p_IfaceIdx(iface) = i
        iface = iface + 1
      END IF
    END DO
    
    ! Copy the triangulation pointer to the new mesh region
    rmeshRegion%p_rtriangulation => p_rtria
    rmeshRegion%NAT = iNAT
    
    ! Now we can deallocate the work array
    DEALLOCATE(IfacesInMR)
    
    ! Finally, recalculate the edge and vertex index arrays from the faces
    CALL mshreg_recalcEdgesFromFaces(rmeshRegion)
    CALL mshreg_recalcVerticesFromFaces(rmeshRegion)
    
    ! That's it

  END SUBROUTINE

END MODULE