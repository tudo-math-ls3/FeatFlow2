!##############################################################################
!# ****************************************************************************
!# <name> meshadjacency </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the calculation of adjacency graphs
!# for the cells of a mesh. The adjacency graph is similar to the sparsity
!# pattern of a CSR matrix (aka matrix storage technique 7/9), describing
!# which cells are adjacent to each other in respect to a specific adjacency
!# relation. The adjacency graph can be seen as a more extended version of
!# the mesh's IneighboursAtElement array, which is used e.g. for the
!# generation of a coloured partitioning of the mesh's cells.
!#
!# Currently, this module supports up to 3 adjacency types (depending on the
!# mesh's dimension):
!# 1. vertice-based element adjacency (1D/2D/3D):
!#    -> Two elements of the mesh are adjacent if they share at least one
!#       common vertice.
!#
!# 2. edge-based element adjacency (2D/3D):
!#    -> Two elements of the mesh are adjacent if they share at least one
!#       common edge.
!#
!# 3. face-based element adjacency (3D):
!#    -> Two elements of the mesh are adjacent if they share at least one
!#       common face.
!#
!# Additionally, an dimension-dependend alias called "neighbour-based element
!# adjacency" is defined, which is equivalent to "vertice-based" in 1D,
!# "edge-based" in 2D and "face-based" in 3D.
!#
!# The following routines can be found here:
!#
!# 1.) mshadj_calcElementAdjacency
!#     -> Calculates the adjacency graph for a given mesh in repsect to a
!#        specific adjacency relation.
!# </purpose>
!##############################################################################

MODULE meshadjacency

  USE fsystem
  USE storage
  USE triangulation

  IMPLICIT NONE

!<constants>

!<constantblock description="mesh adjacency type identifiers">

  ! Element adjacency is calculated based on the IneighboursAtElement array
  ! from the triangulation.
  INTEGER(I32), PARAMETER :: MSHADJ_ADJ_NEIGHBOUR       = -1
  
  ! Element adjacency is calculated by vertex connectivity, i.e. two elements
  ! are adjacent to each other if they share at least one common vertice.
  INTEGER(I32), PARAMETER :: MSHADJ_ADJ_BY_VERTEX       = 1
  
  ! Element adjacency is calculated by edge connectivity, i.e. two elements
  ! are adjacent to each other if they share at least one common edge.
  ! In the case of a 1D grid, this is equal to MSHADJ_ADJ_BY_VERTEX.
  INTEGER(I32), PARAMETER :: MSHADJ_ADJ_BY_EDGE         = 2

  ! Element adjacency is calculated by face connectivity, i.e. two elements
  ! are adjacent to each other if they share at least one common face.
  ! In the case of a 1D grid, this is equal to MSHADJ_ADJ_BY_VERTEX.
  ! In the case of a 2D grid, this is equal to MSHADJ_ADJ_BY_EDGE.
  INTEGER(I32), PARAMETER :: MSHADJ_ADJ_BY_FACE         = 3

!</constantblock>

!</constants>

CONTAINS

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE mshadj_calcElementAdjacency(rtria, h_Iptr, h_Iidx, cadjacency)

!<description>
  ! This routine calculates the element adjacency arrays for a given
  ! triangulation.
  ! The triangulation is silently assumed to be a conformal standard mesh.
!</description>

!<input>
  ! The triangulation structure for which the adjencency is to be calculated.
  TYPE(t_triangulation), INTENT(IN) :: rtria
  
  ! OPTIONAL: One of the MSHADJ_ADJ_XXXX identifiers which specifies how the
  ! adjacency is to be calculated. If not given, MSHADJ_ADJ_NEIGHBOUR is used.
  INTEGER(I32), OPTIONAL, INTENT(IN) :: cadjacency
!</input>

!<output>
  ! A storage handle to the pointer-array of the adjacency graph.
  INTEGER(I32), INTENT(OUT) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  INTEGER(I32), INTENT(OUT) :: h_Iidx
!</output>

!</subroutine>

  INTEGER, DIMENSION(:), POINTER :: p_Iptr, p_Iidx, p_IadjIdx, p_Iadj
  INTEGER, DIMENSION(:,:), POINTER :: p_IneighboursAtElement,&
    p_IverticesAtElement, p_IedgesAtElement, p_Iprim
  INTEGER, DIMENSION(:), POINTER :: p_IelementsAtVertexIdx, p_IelementsAtVertex,&
    p_IelementsAtEdgeIdx3D, p_IelementsAtEdge3D
  INTEGER :: NEL,NNZE,NIAT,NIMT,i,j,k,iel,adj
  INTEGER(I32) :: cadj
  LOGICAL :: bFound
  
    ! Reset the handles
    h_Iptr = ST_NOHANDLE
    h_Iidx = ST_NOHANDLE
    
    ! Get the correct adjacency type. If the given adjacency type parameter
    ! is equivalent to MSHADJ_ADJ_NEIGHBOUR, we'll explicitly set it to
    ! MSHADJ_ADJ_NEIGHBOUR, as this adjacency is faster to compute than
    ! the more complex ones...
    cadj = MSHADJ_ADJ_NEIGHBOUR
    IF(PRESENT(cadjacency)) THEN
      SELECT CASE(cadjacency)
      CASE (MSHADJ_ADJ_NEIGHBOUR)
        cadj = MSHADJ_ADJ_NEIGHBOUR
      CASE (MSHADJ_ADJ_BY_VERTEX)
        IF(rtria%ndim .EQ. NDIM1D) THEN
          cadj = MSHADJ_ADJ_NEIGHBOUR
        ELSE
          cadj = MSHADJ_ADJ_BY_VERTEX
        END IF
      CASE (MSHADJ_ADJ_BY_EDGE)
        IF((rtria%ndim .EQ. NDIM1D) .OR. (rtria%ndim .EQ. NDIM2D)) THEN
          cadj = MSHADJ_ADJ_NEIGHBOUR
        ELSE
          cadj = MSHADJ_ADJ_BY_EDGE
        END IF
      CASE (MSHADJ_ADJ_BY_FACE)
        cadj = MSHADJ_ADJ_NEIGHBOUR
      END SELECT
    END IF
    
    ! Get all the arrays from the triangulation
    CALL storage_getbase_int2d(rtria%h_IneighboursAtElement,  p_IneighboursAtElement)
    CALL storage_getbase_int2d(rtria%h_IverticesAtElement, p_IverticesAtElement)
    IF (rtria%h_IedgesAtElement .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int2d(rtria%h_IedgesAtElement, p_IedgesAtElement)
    ELSE
      NULLIFY(p_IedgesAtElement)
    END IF
    IF (rtria%h_IelementsAtVertexIdx .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
      CALL storage_getbase_int(rtria%h_IelementsAtVertex, p_IelementsAtVertex)
    ELSE
      NULLIFY(p_IelementsAtVertexIdx)
      NULLIFY(p_IelementsAtVertex)
    END IF
    IF (rtria%h_IelementsAtEdgeIdx3D .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IelementsAtEdgeIdx3D, p_IelementsAtEdgeIdx3D)
      CALL storage_getbase_int(rtria%h_IelementsAtEdge3D, p_IelementsAtEdge3D)
    ELSE
      NULLIFY(p_IelementsAtEdgeIdx3D)
      NULLIFY(p_IelementsAtEdge3D)
    END IF

    ! Get the total number of elements from the mesh
    NEL = rtria%NEL

    ! Now comes the first interesting part - calculate the number of non-zero
    ! entries in our adjacency array...
    SELECT CASE(rtria%ndim)
    CASE (NDIM1D)
      ! In the 1D case every element is adjacent to 3 elements - except for the
      ! two at the boundary.
      NNZE = 3*NEL - 2
    
    CASE (NDIM2D)
    
      ! Calculate the total number of inner edges in the mesh
      NIMT = -rtria%NMT + 3*rtria%InelOfType(TRIA_NVETRI2D)&
                        + 4*rtria%InelOfType(TRIA_NVEQUAD2D)
      
      ! What type of adjacency do we have here?
      IF(cadj .EQ. MSHADJ_ADJ_BY_VERTEX) THEN
      
        ! Element adjacency is calculated based on the vertices.
        ! This is a bit tricky, but since we assume that the mesh is conformal,
        ! we are still able to calculate the number of non-zero entries.
        NNZE = 0
        
        ! Loop through the p_IelementsAtVertexIdx array
        DO j = 1, UBOUND(p_IelementsAtVertexIdx,1)-1
          
          ! Get the number of elements adjacent to vertice j
          k = p_IelementsAtVertexIdx(j+1) - p_IelementsAtVertexIdx(j)
          
          ! Every element adjacent to this vertice is coupled to all
          ! elements adjacent to this vertice...
          NNZE = NNZE + k*k
          
        END DO
        
        ! Now every triangle in the mesh was coupled three times to itself -
        ! one time for each vertice.
        NNZE = NNZE - 2*rtria%InelOfType(TRIA_NVETRI2D)
        
        ! And every quadrilateral was coupled four times to itself.
        NNZE = NNZE - 3*rtria%InelOfType(TRIA_NVEQUAD2D)
        
        ! Now if two elements share a common edge, then we have counted the
        ! adjacency twice - one time for each vertice of the common edge.
        NNZE = NNZE - NIMT
        
      ELSE
      
        ! Element adjacency is calculated based on the edges.
        NNZE = NEL + 2*NIMT
        
      END IF
    
    CASE(NDIM3D)
    
      ! Calculate the total number of inner faces in the mesh
      !NIAT = -rtria%NAT + 6*rtria%InelOfType(TRIA_NVEHEXA3D)
      NIAT = -rtria%NAT + 6*NEL
      
      ! What type of adjacency do we have here?
      IF(cadj .EQ. MSHADJ_ADJ_BY_VERTEX) THEN

        PRINT *, 'ERROR: mshadj_calcElementAdjacency'
        PRINT *, '3D vertice-based adjacency not yet implemented'
        CALL sys_halt()
 
      ELSE IF(cadj .EQ. MSHADJ_ADJ_BY_VERTEX) THEN

        PRINT *, 'ERROR: mshadj_calcElementAdjacency'
        PRINT *, '3D edge-based adjacency not yet implemented'
        CALL sys_halt()
      
      ELSE
        
        ! Element adjacency is calculated based on the faces.
        NNZE = NEL + 2*NIAT
        
      END IF
      
    END SELECT
    
    ! Now call the storage to allocate the memory for our adjacency arrays
    CALL storage_new ('mshadj_calcElementAdjacency', 'p_Iptr', NEL+1, ST_INT, &
                      h_Iptr, ST_NEWBLOCK_NOINIT)
    CALL storage_new ('mshadj_calcElementAdjacency', 'p_Iidx', NNZE, ST_INT, &
                      h_Iidx, ST_NEWBLOCK_NOINIT)
    
    ! And get the pointers
    CALL storage_getbase_int(h_Iptr, p_Iptr)
    CALL storage_getbase_int(h_Iidx, p_Iidx)
    
    ! Do we calculate the adjacency based on the neighbours?
    IF(cadj .EQ. MSHADJ_ADJ_NEIGHBOUR) THEN
    
      ! This is the easy case - basically we only need to copy the
      ! neighbours-at-element array...
      k = 1
      DO i = 1, NEL
      
        ! Set the pointer for this element
        p_Iptr(i) = k
        
        ! The element is always adjacent to itself
        p_Iidx(k) = i
        k = k+1
        
        ! Now go through the neighbourhood
        DO j = 1, UBOUND(p_IneighboursAtElement,1)
        
          IF(p_IneighboursAtElement(j,i) .GT. 0) THEN
            
            ! Add the neighbour to the adjacency
            p_Iidx(k) = p_IneighboursAtElement(j,i)
            k = k+1
            
          END IF
          
        END DO
      
      END DO
      p_Iptr(NEL+1) = k
      
      ! That's it
      RETURN
      
    END IF
    
    ! Now we need to get the 3 arrays to assemble the adjacency from.
    ! These 3 arrays are:
    ! Iprim   = I****atElement
    ! IadjIdx = IelementsAt****Idx
    ! Iadj    = IelementsAt****
    ! **** must be replaced with 'Vertices' or 'Edges', depending on whether
    ! we assemble the adjacency based on vertice or edge connectivity.
    SELECT CASE(cadj)
    CASE (MSHADJ_ADJ_BY_VERTEX)
      p_Iprim   => p_IverticesAtElement
      p_IadjIdx => p_IelementsAtVertexIdx
      p_Iadj    => p_IelementsAtVertex
    
    CASE (MSHADJ_ADJ_BY_EDGE)
      p_Iprim   => p_IedgesAtElement
      p_IadjIdx => p_IelementsAtEdgeIdx3D
      p_Iadj    => p_IelementsAtEdge3D
      
    END SELECT
    
    ! Now loop through all elements in the mesh
    NNZE = 0
    p_Iptr(1) = 1
    DO iel = 1, NEL
    
      ! First of all, the element is adjacent to itself
      NNZE = NNZE+1
      p_Iidx(NNZE) = iel
    
      ! Loop through all vertices/edges at that element
      DO i = 1, UBOUND(p_Iprim,1)
      
        ! Skip this vertice/edge if it is zero (might happen for mixed
        ! triangle/quadrilateral meshes)
        IF(p_Iprim(i,iel) .LE. 0) CYCLE
        
        ! And loop through all elements adjacent to that vertice/edge
        DO j = p_IadjIdx(p_Iprim(i,iel)), p_IadjIdx(p_Iprim(i,iel)+1)-1
        
          ! Get the element number
          adj = p_Iadj(j)
          
          ! We assume that we have to add the element adjacency
          bFound = .FALSE.
          
          ! Go through the list
          DO k = p_Iptr(iel), NNZE
          
            ! Is this the element we want to add the adjacency to?
            IF(p_Iidx(k) .EQ. adj) THEN
              ! Yes, so we don't need to add it anymore
              bFound = .TRUE.
              EXIT
            END IF
          
          END DO ! k
          
          ! Do we already have this adjacency?
          IF(bFound) CYCLE
          
          ! Okay, add the entry at the end of our list
          NNZE = NNZE+1
          p_Iidx(NNZE) = adj
          
        END DO ! j
      
      END DO ! i

      ! Set the offset for the next element
      p_Iptr(iel+1) = NNZE+1
    
    END DO ! iel
    
    ! That's it

  END SUBROUTINE

END MODULE