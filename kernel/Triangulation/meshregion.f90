!##############################################################################
!# ****************************************************************************
!# <name> meshregion </name>
!# ****************************************************************************
!#
!# <purpose>
!# TODO
!# </purpose>
!##############################################################################

MODULE meshregion

  USE fsystem
  USE storage
  USE triangulation

  IMPLICIT NONE

!<constants>

!<constantblock description="Share flags for mesh regions">

  ! A mesh region shares its IvertexIdx array with another structure
  INTEGER(I32), PARAMETER :: MESH_REGION_SHARE_VERTEX  = 2**0

  ! A mesh region shares its IedgeIdx array with another structure
  INTEGER(I32), PARAMETER :: MESH_REGION_SHARE_EDGE    = 2**1

  ! A mesh region shares its IfaceIdx array with another structure
  INTEGER(I32), PARAMETER :: MESH_REGION_SHARE_FACE    = 2**2

  ! A mesh region shares its IelementIdx array with another structure
  INTEGER(I32), PARAMETER :: MESH_REGION_SHARE_ELEMENT = 2**3

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Structure defining a mesh region.
  TYPE t_meshRegion
  
    ! A pointer to the triangulation on which this mesh region is defined.
    TYPE(t_triangulation), POINTER :: p_rtriangulation => NULL()
    
    ! Flags of the mesh region.
    INTEGER(I32) :: cflags = 0
    
    ! Number of vertices in the mesh region.
    INTEGER(PREC_VERTEXIDX) :: NVT = 0
    
    ! Number of edges in the mesh region (only 2D/3D).
    INTEGER(PREC_EDGEIDX) :: NMT = 0
    
    ! Number of faces in the mesh region (only 3D).
    INTEGER(PREC_FACEIDX) :: NAT = 0
    
    ! Number of elements in the mesh region.
    INTEGER(PREC_ELEMENTIDX) :: NEL = 0
    
    ! Handle to integer array holding the indices of all vertices which
    ! belong to this mesh region.
    INTEGER :: h_IvertexIdx = ST_NOHANDLE
    
    ! Handle to integer array holding the indices of all edges which
    ! belong to this mesh region (only 2D/3D).
    INTEGER :: h_IedgeIdx = ST_NOHANDLE
    
    ! Handle to integer array holding the indices of all faces which
    ! belong to this mesh region (only 3D).
    INTEGER :: h_IfaceIdx = ST_NOHANDLE
    
    ! Handle to integer array holding the indces of all cells which
    ! belong to this mesh region.
    INTEGER :: h_IelementIdx = ST_NOHANDLE
    
  END TYPE

!</typeblock>

!</types>

  CONTAINS

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_createFromNodalProp(rmeshRegion, p_rtriangulation, &
                                        inodalProperty)

!<description>
  ! Creates a new mesh region based on the nodal property array of
  ! the triangulation
!</description>

!<input>
  ! The triangulation that is to be used.
  TYPE(t_triangulation), TARGET, INTENT(IN) :: p_rtriangulation
  
  ! OPTIONAL: A nodal property for which the region is to be created.
  ! If present, all vertices, edges and faces which have the specified
  ! nodal property will be added to this mesh region.
  ! If not present, all vertices, edges and faces which have a positive
  ! nodal property will be added to this mesh region.
  INTEGER, OPTIONAL, INTENT(IN) :: inodalProperty
!</input>

!<output>
  ! The mesh region of the triangulation.
  TYPE(t_meshRegion), INTENT(OUT) :: rmeshRegion
!</output>

!</subroutine>

  ! Some local variables
  INTEGER :: iNVT, iNMT, iNAT, i, idim, ioff
  INTEGER, DIMENSION(:), POINTER :: p_InodalProperty, p_Idata
  
    ! Hang in the triangulation
    rmeshRegion%p_rtriangulation => p_rtriangulation
    
    ! Get the nodal property array of the triangulation
    CALL storage_getbase_int(p_rtriangulation%h_InodalProperty, p_InodalProperty)
    
    ! Get the dimension of the triangulation
    idim = p_rtriangulation%ndim
    
    ! First find out how many vertices, edges and faces belong to this
    ! mesh region.
    iNVT = 0
    iNMT = 0
    iNAT = 0
    IF(PRESENT(inodalProperty)) THEN
      
      ! Go through all vertices in the mesh.
      DO i = 1, p_rtriangulation%NVT
        
        IF (p_InodalProperty(i) .EQ. inodalProperty) iNVT = iNVT + 1
        
      END DO
      
      ! If this mesh is at least a 2D mesh, go through all edges
      IF (idim .GE. 2) THEN
      
        ioff = p_rtriangulation%NVT
        
        DO i = 1, p_rtriangulation%NMT
          
          IF (p_InodalProperty(ioff+i) .EQ. inodalProperty) iNMT = iNMT + 1
          
        END DO
        
      END IF
      
      ! If this mesh is a 3D mesh, go through all faces
      IF (idim .GE. 3) THEN
      
        ioff = p_rtriangulation%NVT + p_rtriangulation%NMT

        DO i = 1, p_rtriangulation%NAT
          
          IF (p_InodalProperty(ioff+i) .EQ. inodalProperty) iNAT = iNAT + 1
          
        END DO
        
      END IF
    
    ELSE
      ! No nodal property was given - so get all vertices, edges and faces
      ! which have positive nodal property...
    
      ! Go through all vertices in the mesh.
      DO i = 1, p_rtriangulation%NVT
        
        IF (p_InodalProperty(i) .GT. 0) iNVT = iNVT + 1
        
      END DO
      
      ! If this mesh is at least a 2D mesh, go through all edges
      IF (idim .GE. 2) THEN
      
        ioff = p_rtriangulation%NVT
        
        DO i = 1, p_rtriangulation%NMT
          
          IF (p_InodalProperty(ioff+i) .GT. 0) iNMT = iNMT + 1
          
        END DO
        
      END IF
      
      ! If this mesh is a 3D mesh, go through all faces
      IF (idim .GE. 3) THEN
      
        ioff = p_rtriangulation%NVT + p_rtriangulation%NMT

        DO i = 1, p_rtriangulation%NAT
          
          IF (p_InodalProperty(ioff+i) .GT. 0) iNAT = iNAT + 1
          
        END DO
        
      END IF

    END IF
    
    ! If the mesh region is empty, print a warning and exit here
    IF ((iNVT .EQ. 0) .AND. (iNMT .EQ. 0) .AND. (iNAT .EQ. 0)) THEN
    
      PRINT *, 'Warning: Created empty mesh region'
      RETURN
      
    ENDIF
    
    IF (iNVT .GT. 0) THEN
    
      ! Allocate the vertice index array
      CALL storage_new('mshreg_createFromNodalProp', 'p_IvertexIdx', iNVT, &
                       ST_INT, rmeshRegion%h_IvertexIdx, ST_NEWBLOCK_NOINIT)
      
      ! get the vertice index array
      CALL storage_getbase_int(rmeshRegion%h_IvertexIdx, p_Idata)
      
      ! Now go and map the vertices
      iNVT = 0
      IF(PRESENT(inodalProperty)) THEN

        DO i = 1, p_rtriangulation%NVT
          
          IF (p_InodalProperty(i) .EQ. inodalProperty) THEN
            iNVT = iNVT + 1
            p_Idata(iNVT) = i
          END IF
          
        END DO
      
      ELSE
      
        DO i = 1, p_rtriangulation%NVT
          
          IF (p_InodalProperty(i) .GT. 0) THEN
            iNVT = iNVT + 1
            p_Idata(iNVT) = i
          END IF
          
        END DO
        
      END IF
      
      ! store the vertice count
      rmeshRegion%NVT = iNVT

    END IF
    
    IF (iNMT .GT. 0) THEN
    
      ! Allocate the edge index array
      CALL storage_new('mshreg_createFromNodalProp', 'p_IedgeIdx', iNVT, &
                       ST_INT, rmeshRegion%h_IedgeIdx, ST_NEWBLOCK_NOINIT)
      
      ! get the edge index array
      CALL storage_getbase_int(rmeshRegion%h_IedgeIdx, p_Idata)
      
      ! Now go and map the edges
      iNMT = 0
      ioff = p_rtriangulation%NVT
      IF(PRESENT(inodalProperty)) THEN

        DO i = 1, p_rtriangulation%NMT
          
          IF (p_InodalProperty(ioff+i) .EQ. inodalProperty) THEN
            iNMT = iNMT + 1
            p_Idata(iNMT) = i
          END IF
          
        END DO
      
      ELSE
      
        DO i = 1, p_rtriangulation%NMT
          
          IF (p_InodalProperty(ioff+i) .GT. 0) THEN
            iNMT = iNMT + 1
            p_Idata(iNMT) = i
          END IF
          
        END DO
        
      END IF
      
      ! store the edge count
      rmeshRegion%NMT = iNMT

    END IF

    IF (iNAT .GT. 0) THEN
    
      ! Allocate the face index array
      CALL storage_new('mshreg_createFromNodalProp', 'p_IfaceIdx', iNVT, &
                       ST_INT, rmeshRegion%h_IfaceIdx, ST_NEWBLOCK_NOINIT)
      
      ! get the face index array
      CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_Idata)
      
      ! Now go and map the faces
      iNAT = 0
      ioff = p_rtriangulation%NVT + p_rtriangulation%NMT
      IF(PRESENT(inodalProperty)) THEN

        DO i = 1, p_rtriangulation%NAT
          
          IF (p_InodalProperty(ioff+i) .EQ. inodalProperty) THEN
            iNAT = iNAT + 1
            p_Idata(iNAT) = i
          END IF
          
        END DO
      
      ELSE
      
        DO i = 1, p_rtriangulation%NAT
          
          IF (p_InodalProperty(ioff+i) .GT. 0) THEN
            iNAT = iNAT + 1
            p_Idata(iNAT) = i
          END IF
          
        END DO
        
      END IF
      
      ! store the face count
      rmeshRegion%NAT = iNAT

    END IF
    
    ! That's it

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_done(rmeshRegion)

!<description>
  ! Releases a mesh region structure and frees all allocated memory from
  ! the heap.
!</description>

!<inputoutput>
  TYPE(t_meshRegion), INTENT(INOUT) :: rmeshRegion
!</inputoutput>

!</subroutine>

    ! If the region does not share its vertex index array, release it
    IF ((rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) .AND. &
        (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_VERTEX) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IvertexIdx)
        
    END IF
  
    ! If the region does not share its edge index array, release it
    IF ((rmeshRegion%h_IedgeIdx .NE. ST_NOHANDLE) .AND. &
        (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_EDGE) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IedgeIdx)
        
    END IF
  
    ! If the region does not share its face index array, release it
    IF ((rmeshRegion%h_IfaceIdx .NE. ST_NOHANDLE) .AND. &
        (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_FACE) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IfaceIdx)
        
    END IF
  
    ! If the region does not share its element index array, release it
    IF ((rmeshRegion%h_IelementIdx .NE. ST_NOHANDLE) .AND. &
        (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_ELEMENT) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IelementIdx)
        
    END IF
    
    ! Last but not least: reset the structure
    rmeshRegion%p_rtriangulation => NULL()
    rmeshRegion%cflags = 0
    rmeshRegion%NVT = 0
    rmeshRegion%NMT = 0
    rmeshRegion%NAT = 0
    rmeshRegion%NEL = 0
    ! Although usually the handles are reset by the storage_free() routine,
    ! we need to reset the handles by hand, too, because if a handle was given,
    ! and the array represented by the handle was marked as shared in the flags,
    ! the storage_free() routine is not called.
    rmeshRegion%h_IvertexIdx = ST_NOHANDLE
    rmeshRegion%h_IedgeIdx = ST_NOHANDLE
    rmeshRegion%h_IfaceIdx = ST_NOHANDLE
    rmeshRegion%h_IelementIdx = ST_NOHANDLE
    
    ! That's it

  END SUBROUTINE
  
END MODULE
