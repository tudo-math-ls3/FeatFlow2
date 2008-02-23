!##############################################################################
!# ****************************************************************************
!# <name> meshregion </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains structures and routines for managing a sepcific region
!# inside a given triangulation. These region may be used e.g. for boundary
!# condition assembly.
!#
!# The following routines can be found here:
!#
!#  1.) mshreg_createFromNodalProp
!#      -> Creates a mesh region based on the triangulation's
!#         nodal property array.
!#
!#  2.) mshred_done
!#      -> Releases a mesh region structure.
!#
!#  3.) mshreg_recalcEdgesFromFaces
!#      -> Recalculates the edge index array of the mesh region from the
!#         face index array of the mesh region.
!#
!#  4.) mshreg_recalcVerticesFromFaces
!#      -> Recalculates the vertex index array of the mesh region from the
!#         face index array of the mesh region.
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

!<constantblock description="Mask operator types">

  ! Kick the old index array
  INTEGER(I32), PARAMETER :: MESH_REGION_MASK_KICK     = 0
  
  ! Apply the OR-operator
  INTEGER(I32), PARAMETER :: MESH_REGION_MASK_OR       = 1
  
  ! Apply the AND-operator
  INTEGER(I32), PARAMETER :: MESH_REGION_MASK_AND      = 2

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

  PRIVATE :: mshreg_calcArrayFromMap

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
      CALL storage_new('mshreg_createFromNodalProp', 'p_IedgeIdx', iNMT, &
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
      CALL storage_new('mshreg_createFromNodalProp', 'p_IfaceIdx', iNAT, &
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

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_calcArrayFromMap(Imap,h_Iarray,inumEntries)

!<description>
  ! This auxiliary routine creates an index array from a map.
  !
  ! This subroutine is used by some other routines in this module and
  ! is not meant to be called from outside this module, therefore, it is
  ! declared as private.
!</description>

!<input>
  ! The map that is to be converted into an index array.
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: Imap
  
!</input>

!<output>
  ! A storage handle to the index array.
  INTEGER(I32), INTENT(OUT) :: h_Iarray
  
  ! The number of entries in the index array
  INTEGER, INTENT(OUT) :: inumEntries
!</output>

!</subroutine>

  ! Three nice local variables
  INTEGER :: i,ientry
  INTEGER, DIMENSION(:), POINTER :: p_Iarray

    h_Iarray = ST_NOHANDLE

    ! First of all, count the entries
    inumEntries = 0
    DO i = LBOUND(Imap,1), UBOUND(Imap,1)
      IF (Imap(i) .NE. 0) inumEntries = inumEntries + 1
    END DO
    
    ! If there are no entries, we can leave here
    IF (inumEntries .EQ. 0) RETURN
    
    ! Otherwise allocate an index array
    CALL storage_new('mshreg_calcArrayFromMap', 'p_Iarray', inumEntries, &
                     ST_INT, h_Iarray, ST_NEWBLOCK_NOINIT)
    
    ! Get the index array
    CALL storage_getbase_int(h_Iarray, p_Iarray)
    
    ! Go and calculate the indices
    ientry = 1
    DO i = LBOUND(Imap,1), UBOUND(Imap,1)
      IF (Imap(i) .NE. 0) THEN
        p_Iarray(ientry) = i
        ientry = ientry + 1
      END IF
    END DO
    
    ! That's it
    
  END SUBROUTINE
  
!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_recalcEdgesFromFaces(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the edge index array from the face index array.
  ! Every edge that belongs to a face in the mesh region will be added
  ! to the mesh region's edge index array.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old edge index array.
  ! One of the MESH_REGION_MASK_XXXX constants defined above.
  ! If not given, MESH_REGION_MASK_KICK is used.
  INTEGER(I32), INTENT(IN), OPTIONAL   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  TYPE(t_meshRegion), INTENT(INOUT)    :: rmeshRegion
!</inputoutput>

!</subroutine>
  
  INTEGER(I32), DIMENSION(:), ALLOCATABLE :: IedgeMap
  INTEGER, DIMENSION(:), POINTER :: p_IfaceIdx, p_IedgeIdx
  INTEGER, DIMENSION(:,:), POINTER :: p_IedgesAtFace
  INTEGER(I32) :: copMask
  INTEGER :: i, j, iface
  TYPE(t_triangulation), POINTER :: p_rtria

    ! Decide what operator to use
    copMask = MESH_REGION_MASK_KICK
    IF (PRESENT(coperatorMask)) copMask = coperatorMask
    
    ! Do we have any faces at all?
    IF ((rmeshRegion%NAT .LE. 0) .OR. (rmeshRegion%h_IfaceIdx .EQ. &
         ST_NOHANDLE)) THEN
    
      ! If the mesh region already had an edge index array and we do not
      ! apply the OR-operator, kick the old edge index array.
      IF (copMask .NE. MESH_REGION_MASK_OR) THEN
      
        IF ((rmeshRegion%h_IedgeIdx .NE. ST_NOHANDLE) .AND. &
            (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_EDGE) .EQ. 0)) THEN
            
            CALL storage_free(rmeshRegion%h_IedgeIdx)
            
        END IF
        
        ! Remove anything that has to do with edges from the mesh region
        rmeshRegion%NMT = 0
        rmeshRegion%h_IedgeIdx = ST_NOHANDLE
        rmeshRegion%cflags = IAND(rmeshRegion%cflags, &
                              NOT(MESH_REGION_SHARE_EDGE))
      
      END IF
      
      ! Return here
      RETURN
    
    END IF
    
    p_rtria => rmeshRegion%p_rtriangulation
    
    ! Get the face index array
    CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! Get the edges-at-face array from the triangulation
    CALL storage_getbase_int2D(p_rtria%h_IedgesAtFace,&
                               p_IedgesAtFace)
    
    ! Allocate an edge map
    ALLOCATE(IedgeMap(p_rtria%NMT))
    DO i=1, p_rtria%NMT
      IedgeMap(i) = 0
    END DO
    
    ! Go through all faces
    DO i=1, rmeshRegion%NAT
    
      ! Get the face index
      iface = p_IfaceIdx(i)
      
      ! Go through all edges adjacent to that face
      DO j=1, UBOUND(p_IedgesAtFace,1)
        
        IF (p_IedgesAtFace(j,iface) .GT. 0) THEN
          IedgeMap(p_IedgesAtFace(j,iface)-p_rtria%NVT) = 1
        END IF
        
      END DO
      
    END DO
    
    ! Now let's see what to do with the old edges
    IF (rmeshRegion%h_IedgeIdx .NE. ST_NOHANDLE) THEN
    
      ! Get the old edge index array
      CALL storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)
      
      SELECT CASE(copMask)
      CASE (MESH_REGION_MASK_OR)
        ! Go through all edges and apply an OR-operator.
        DO i=1, rmeshRegion%NMT
          IedgeMap(p_IedgeIdx(i)) = 1
        END DO
      
      CASE (MESH_REGION_MASK_AND)
        ! Go through all edges and apply an AND-operator.
        DO i=1, rmeshRegion%NMT
          IedgeMap(p_IedgeIdx(i)) = IAND(IedgeMap(p_IedgeIdx(i)),1)
        END DO
        
      END SELECT
      
      ! Now let's destroy the old edge index array
      IF (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_EDGE) .EQ. 0) THEN
          CALL storage_free(rmeshRegion%h_IedgeIdx)
      END IF
      
      ! Remove anything that has to do with edges from the mesh region
      rmeshRegion%NMT = 0
      rmeshRegion%h_IedgeIdx = ST_NOHANDLE
      rmeshRegion%cflags = IAND(rmeshRegion%cflags, &
                            NOT(MESH_REGION_SHARE_EDGE))
      
    END IF
    
    ! Now we have the edge index map, so create an index array from this.
    CALL mshreg_calcArrayFromMap(IedgeMap, rmeshRegion%h_IedgeIdx, &
                                 rmeshRegion%NMT)
    
    ! And release the edge map
    DEALLOCATE(IedgeMap)
    
    ! That's it
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_recalcVerticesFromFaces(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the vertice index array from the face index array.
  ! Every vertice that belongs to a face in the mesh region will be added
  ! to the mesh region's vertice index array.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old vertice index array.
  ! One of the MESH_REGION_MASK_XXXX constants defined above.
  ! If not given, MESH_REGION_MASK_KICK is used.
  INTEGER(I32), INTENT(IN), OPTIONAL   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  TYPE(t_meshRegion), INTENT(INOUT)    :: rmeshRegion
!</inputoutput>

!</subroutine>
  
  INTEGER(I32), DIMENSION(:), ALLOCATABLE :: IvertMap
  INTEGER, DIMENSION(:), POINTER :: p_IfaceIdx, p_IvertIdx
  INTEGER, DIMENSION(:,:), POINTER :: p_IvertsAtFace
  INTEGER(I32) :: copMask
  INTEGER :: i, j, iface
  TYPE(t_triangulation), POINTER :: p_rtria

    ! Decide what operator to use
    copMask = MESH_REGION_MASK_KICK
    IF (PRESENT(coperatorMask)) copMask = coperatorMask
    
    ! Do we have any faces at all?
    IF ((rmeshRegion%NAT .LE. 0) .OR. (rmeshRegion%h_IfaceIdx .EQ. &
         ST_NOHANDLE)) THEN
    
      ! If the mesh region already had an vertex index array and we do not
      ! apply the OR-operator, kick the old vertex index array.
      IF (copMask .NE. MESH_REGION_MASK_OR) THEN
      
        IF ((rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) .AND. &
            (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_VERTEX) .EQ. 0)) THEN
            
            CALL storage_free(rmeshRegion%h_IvertexIdx)
            
        END IF
        
        ! Remove anything that has to do with vertices from the mesh region
        rmeshRegion%NVT = 0
        rmeshRegion%h_IvertexIdx = ST_NOHANDLE
        rmeshRegion%cflags = IAND(rmeshRegion%cflags, &
                              NOT(MESH_REGION_SHARE_VERTEX))
      
      END IF
      
      ! Return here
      RETURN
    
    END IF
    
    p_rtria => rmeshRegion%p_rtriangulation
    
    ! Get the face index array
    CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! Get the vertices-at-face array from the triangulation
    CALL storage_getbase_int2D(p_rtria%h_IverticesAtFace,&
                               p_IvertsAtFace)
    
    ! Allocate an vertex map
    ALLOCATE(IvertMap(p_rtria%NVT))
    DO i=1, p_rtria%NVT
      IvertMap(i) = 0
    END DO
    
    ! Go through all faces
    DO i=1, rmeshRegion%NAT
    
      ! Get the face index
      iface = p_IfaceIdx(i)
      
      ! Go through all edges adjacent to that face
      DO j=1, UBOUND(p_IvertsAtFace,1)
        
        IF (p_IvertsAtFace(j,iface) .GT. 0) THEN
          IvertMap(p_IvertsAtFace(j,iface)) = 1
        END IF
        
      END DO
      
    END DO
    
    ! Now let's see what to do with the old vertices
    IF (rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) THEN
    
      ! Get the old vertices index array
      CALL storage_getbase_int(rmeshRegion%h_IvertexIdx, p_IvertIdx)
      
      SELECT CASE(copMask)
      CASE (MESH_REGION_MASK_OR)
        ! Go through all vertices and apply an OR-operator.
        DO i=1, rmeshRegion%NVT
          IvertMap(p_IvertIdx(i)) = 1
        END DO
      
      CASE (MESH_REGION_MASK_AND)
        ! Go through all vertices and apply an AND-operator.
        DO i=1, rmeshRegion%NVT
          IvertMap(p_IvertIdx(i)) = IAND(IvertMap(p_IvertIdx(i)),1)
        END DO
        
      END SELECT
      
      ! Now let's destroy the old vertex index array
      IF ((rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) .AND. &
          (IAND(rmeshRegion%cflags, MESH_REGION_SHARE_VERTEX) .EQ. 0)) THEN
          
          CALL storage_free(rmeshRegion%h_IvertexIdx)
          
      END IF
      
      ! Remove anything that has to do with vertices from the mesh region
      rmeshRegion%NVT = 0
      rmeshRegion%h_IvertexIdx = ST_NOHANDLE
      rmeshRegion%cflags = IAND(rmeshRegion%cflags, &
                            NOT(MESH_REGION_SHARE_VERTEX))
      
    END IF
    
    ! Now we have the vertex index map, so create an index array from this.
    CALL mshreg_calcArrayFromMap(IvertMap, rmeshRegion%h_IvertexIdx, &
                                 rmeshRegion%NVT)
    
    ! And release the vertex map
    DEALLOCATE(IvertMap)
    
    ! That's it
    
  END SUBROUTINE

END MODULE
