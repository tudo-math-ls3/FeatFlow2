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
!#  2.) mshreg_createFromHitTest
!#      -> Creates a mesh region based on a hit-test callback function.
!#
!#  3.) mshreg_done
!#      -> Releases a mesh region structure.
!#
!#  4.) mshreg_recalcVerticesFromEdges
!#      -> Recalculates the vertex index array of the mesh region from the
!#         edge index array of the mesh region.
!#
!#  5.) mshreg_recalcVerticesFromFaces
!#      -> Recalculates the vertex index array of the mesh region from the
!#         face index array of the mesh region.
!#
!#  6.) mshreg_recalcEdgesFromFaces
!#      -> Recalculates the edge index array of the mesh region from the
!#         face index array of the mesh region.
!#
!#  7.) mshreg_calcBoundaryNormals2D
!#      -> Calculates the inner normal vectors for the boundary edges of a 
!#         2D mesh region.
!#
!# </purpose>
!##############################################################################

MODULE meshregion

  USE fsystem
  USE storage
  USE triangulation
  USE collection
  USE geometryaux

  IMPLICIT NONE

!<constants>

!<constantblock description="Cell type identificator flags">

  ! Identification flag for the vertice index array
  INTEGER(I32), PARAMETER :: MSHREG_IDX_VERTEX  = 2**0

  ! Identification flag for the edge index array
  INTEGER(I32), PARAMETER :: MSHREG_IDX_EDGE    = 2**1

  ! Identification flag for the face index array
  INTEGER(I32), PARAMETER :: MSHREG_IDX_FACE    = 2**2

  ! Identification flag for the element index array
  INTEGER(I32), PARAMETER :: MSHREG_IDX_ELEMENT = 2**3
  
  ! Identification flag for no index arrays
  INTEGER(I32), PARAMETER :: MSHREG_IDX_NONE    = 0
  
  ! Identification flag for all index arrays
  INTEGER(I32), PARAMETER :: MSHREG_IDX_ALL     = 2**0 + 2**1 + 2**2 + 2**3

!</constantblock>

!<constantblock description="Mask operator types">

  ! Kick the old index array
  INTEGER(I32), PARAMETER :: MSHREG_MASK_KICK     = 0
  
  ! Apply the OR-operator
  INTEGER(I32), PARAMETER :: MSHREG_MASK_OR       = 1
  
  ! Apply the AND-operator
  INTEGER(I32), PARAMETER :: MSHREG_MASK_AND      = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Structure defining a mesh region.
  TYPE t_meshRegion
  
    ! A pointer to the triangulation on which this mesh region is defined.
    TYPE(t_triangulation), POINTER :: p_rtriangulation => NULL()
    
    ! Share flags of the mesh region. A combination of MSHREG_IDX_XXXX
    ! constants defined above. If a MSHREG_IDX_XXXX flag is defined,
    ! then the mesh region shares the corresponding index array with
    ! another structure.
    INTEGER(I32) :: cshareFlags = 0
    
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

  PRIVATE :: mshreg_aux_calcIdxArray1, mshreg_aux_calcIdxArray2

  CONTAINS

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_aux_calcIdxArray1(IpropArray,ioff,ilen,h_IidxArray,&
                                     inumEntries,iidx,IallowedProp)

!<description>
  ! This auxiliary routine creates an index array from a nodal property
  ! array.
  !
  ! This subroutine is used by some other routines in this module and
  ! is not meant to be called from outside this module, therefore, it is
  ! declared as private.
!</description>

!<input>
  ! The nodal property array from which the index array is to be
  ! calculated.
  INTEGER(I32), DIMENSION(:), INTENT(IN)           :: IpropArray
  
  ! The offset of the first entry in the nodal property array.
  INTEGER, INTENT(IN)                              :: ioff
  
  ! The number of entries in the nodal property array.
  INTEGER, INTENT(IN)                              :: ilen

  ! ...  
  INTEGER, INTENT(IN)                              :: iidx
  
  ! OPTIONAL: An array holding all allowed nodal properties that
  ! are to be added into the index array. If not given, all entries
  ! with nodal property > 0 are added.
  INTEGER(I32), DIMENSION(:), OPTIONAL, INTENT(IN) :: IallowedProp
!</input>

!<output>
  ! A storage handle to the created index array.
  INTEGER(I32), INTENT(OUT)                        :: h_IidxArray
  
  ! The number of entries in the created index array.
  INTEGER, INTENT(OUT)                             :: inumEntries
!</output>

!</subroutine>

  ! Three nice local variables
  INTEGER :: i,j,ientry
  INTEGER, DIMENSION(:), POINTER :: p_Iarray

    ! Reset the handle
    h_IidxArray = ST_NOHANDLE

    ! First of all, count the entries
    inumEntries = 0
    IF (PRESENT(IallowedProp)) THEN
      DO i = 1, ilen
        DO j = LBOUND(IallowedProp,1), UBOUND(IallowedProp,1)
          IF (IpropArray(ioff+i) .EQ. IallowedProp(j)) THEN
            inumEntries = inumEntries + 1
            EXIT
          END IF
        END DO
      END DO
    
    ELSE
      DO i = 1, ilen
        IF (IpropArray(ioff+i) .GT. 0) inumEntries = inumEntries + 1
      END DO
    END IF
    
    ! If there are no entries, we can leave here
    IF (inumEntries .EQ. 0) RETURN
    
    ! Otherwise allocate an index array
    CALL storage_new('mshreg_aux_calcIdxArray2', 'p_IidxArray', inumEntries, &
                     ST_INT, h_IidxArray, ST_NEWBLOCK_NOINIT)
    
    ! Get the index array
    CALL storage_getbase_int(h_IidxArray, p_Iarray)
    
    ! Go and calculate the indices
    ientry = 1
    IF (PRESENT(IallowedProp)) THEN

      DO i = 1, ilen
        DO j = LBOUND(IallowedProp,1), UBOUND(IallowedProp,1)
          IF (IpropArray(ioff+i) .EQ. IallowedProp(j)) THEN
            p_Iarray(ientry) = iidx + i
            ientry = ientry + 1
          END IF
        END DO
      END DO
    
    ELSE

      DO i = 1, ilen
        IF (IpropArray(ioff+i) .GT. 0) THEN
          p_Iarray(ientry) = iidx + i
          ientry = ientry + 1
        END IF
      END DO

    END IF
    
    ! That's it
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_aux_calcIdxArray2(IpropArray,ioff,ilen,h_IidxArray,&
                                     inumEntries,Iindex,IallowedProp)

!<description>
  ! This auxiliary routine creates an index array from a nodal property
  ! array.
  !
  ! This subroutine is used by some other routines in this module and
  ! is not meant to be called from outside this module, therefore, it is
  ! declared as private.
!</description>

!<input>
  ! The nodal property array from which the index array is to be
  ! calculated.
  INTEGER(I32), DIMENSION(:), INTENT(IN)           :: IpropArray
  
  ! The offset of the first entry in the nodal property array.
  INTEGER, INTENT(IN)                              :: ioff
  
  ! The number of entries in the nodal property array.
  INTEGER, INTENT(IN)                              :: ilen
  
  ! An array holding the indices of the entries in the
  ! nodal property array.
  INTEGER, DIMENSION(:), INTENT(IN)                :: Iindex
  
  ! OPTIONAL: An array holding all allowed nodal properties that
  ! are to be added into the index array. If not given, all entries
  ! with nodal property > 0 are added.
  INTEGER(I32), DIMENSION(:), OPTIONAL, INTENT(IN) :: IallowedProp
!</input>

!<output>
  ! A storage handle to the created index array.
  INTEGER(I32), INTENT(OUT)                        :: h_IidxArray
  
  ! The number of entries in the created index array.
  INTEGER, INTENT(OUT)                             :: inumEntries
!</output>

!</subroutine>

  ! Three nice local variables
  INTEGER :: i,j,ientry
  INTEGER, DIMENSION(:), POINTER :: p_Iarray

    ! Reset the handle
    h_IidxArray = ST_NOHANDLE

    ! First of all, count the entries
    inumEntries = 0
    IF (PRESENT(IallowedProp)) THEN
      DO i = 1, ilen
        DO j = LBOUND(IallowedProp,1), UBOUND(IallowedProp,1)
          IF (IpropArray(ioff+i) .EQ. IallowedProp(j)) THEN
            inumEntries = inumEntries + 1
            EXIT
          END IF
        END DO
      END DO
    
    ELSE
      DO i = 1, ilen
        IF (IpropArray(ioff+i) .GT. 0) inumEntries = inumEntries + 1
      END DO
    END IF
    
    ! If there are no entries, we can leave here
    IF (inumEntries .EQ. 0) RETURN
    
    ! Otherwise allocate an index array
    CALL storage_new('mshreg_aux_calcIdxArray2', 'p_IidxArray', inumEntries, &
                     ST_INT, h_IidxArray, ST_NEWBLOCK_NOINIT)
    
    ! Get the index array
    CALL storage_getbase_int(h_IidxArray, p_Iarray)
    
    ! Go and calculate the indices
    ientry = 1
    IF (PRESENT(IallowedProp)) THEN

      DO i = 1, ilen
        DO j = LBOUND(IallowedProp,1), UBOUND(IallowedProp,1)
          IF (IpropArray(ioff+i) .EQ. IallowedProp(j)) THEN
            p_Iarray(ientry) = Iindex(i)
            ientry = ientry + 1
          END IF
        END DO
      END DO
        
    ELSE

      DO i = 1, ilen
        IF (IpropArray(ioff+i) .GT. 0) THEN
          p_Iarray(ientry) = Iindex(i)
          ientry = ientry + 1
        END IF
      END DO

    END IF
    
    ! That's it
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_createFromNodalProp(rmeshRegion, rtriangulation, &
                                        cidxCalc, IallowedProp)

!<description>
  ! Creates a new mesh region based on the nodal property array of
  ! the triangulation.
!</description>

!<input>
  ! The triangulation that is to be used.
  TYPE(t_triangulation), TARGET, INTENT(IN)      :: rtriangulation
  
  ! A combination of MSHREG_IDX_XXXX constants defined above which
  ! specifies which index arrays should be calculated from the nodal
  ! property array. May also be MSHREG_IDX_NONE or MSHREG_IDX_ALL.
  INTEGER(I32), INTENT(IN)                       :: cidxCalc
  
  ! OPTIONAL: A nodal property array for which the region is to be created.
  ! If present, all vertices, edges and faces which have one of the nodal
  ! properties specifies in the array will be added to this mesh region.
  ! If not present, all vertices, edges and faces which have a positive
  ! nodal property will be added to this mesh region.
  INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)    :: IallowedProp
!</input>

!<output>
  ! The mesh region of the triangulation.
  TYPE(t_meshRegion), INTENT(OUT)                :: rmeshRegion
!</output>

!</subroutine>

  ! Some local variables
  INTEGER :: idim
  INTEGER, DIMENSION(:), POINTER :: p_InodalProperty
  
    ! Hang in the triangulation
    rmeshRegion%p_rtriangulation => rtriangulation
    
    ! Get the nodal property array of the triangulation
    CALL storage_getbase_int(rtriangulation%h_InodalProperty, p_InodalProperty)
    
    ! Get the dimension of the triangulation
    idim = rtriangulation%ndim
    
    ! Should we calculate the vertice index array?
    IF(IAND(cidxCalc, MSHREG_IDX_VERTEX) .NE. 0) THEN
    
      ! Calculate the vertice index array
      IF (PRESENT(IallowedProp)) THEN

        CALL mshreg_aux_calcIdxArray1(p_InodalProperty, 0, rtriangulation%NVT, &
                                      rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, &
                                      0, IallowedProp)
      ELSE
                                   
        CALL mshreg_aux_calcIdxArray1(p_InodalProperty, 0, rtriangulation%NVT, &
                                     rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0)
      END IF
    
    END IF
    
    ! Should we calculate the edge index array?
    IF((idim .GE. 2) .AND. (IAND(cidxCalc, MSHREG_IDX_EDGE) .NE. 0)) THEN
    
      ! Calculate the edge index array
      IF (PRESENT(IallowedProp)) THEN

        CALL mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT, &
                                      rtriangulation%NMT, rmeshRegion%h_IedgeIdx, &
                                      rmeshRegion%NMT, 0, IallowedProp)
      ELSE
                                   
        CALL mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT, &
                                      rtriangulation%NMT, rmeshRegion%h_IedgeIdx, &
                                      rmeshRegion%NMT, 0)
      END IF
    
    END IF
    
    ! Should we calculate the face index array?
    IF((idim .GE. 3) .AND. (IAND(cidxCalc, MSHREG_IDX_FACE) .NE. 0)) THEN
    
      ! Calculate the face index array
      IF (PRESENT(IallowedProp)) THEN

        CALL mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT + &
                                      rtriangulation%NMT, rtriangulation%NAT, &
                                      rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, &
                                      0, IallowedProp)
      ELSE
                                   
        CALL mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT + &
                                      rtriangulation%NMT, rtriangulation%NAT, &
                                      rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, 0)
      END IF
    
    END IF

    ! That's it

  END SUBROUTINE
  
!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
          cidxCalc, bonlyBndry, fmshregHitTest, IallowedHit, rcollection)

!<description>
  ! Creates a new mesh region by calling a hit-test callback function
  ! to decide whether a cell belongs to the mesh region or not.
!</description>

!<input>
  ! The triangulation that is to be used.
  TYPE(t_triangulation), TARGET, INTENT(IN)      :: rtriangulation
  
  ! A combination of MSHREG_IDX_XXXX constants defined above which
  ! specifies which index arrays should be calculated from the nodal
  ! property array. May also be MSHREG_IDX_NONE or MSHREG_IDX_ALL.
  INTEGER(I32), INTENT(IN)                       :: cidxCalc
  
  ! If .TRUE., then only the cells which belong to the triangulation's
  ! boundary are hit-tested, otherwise all cells are hit-tested.
  LOGICAL, INTENT(IN)                            :: bonlyBndry
  
  ! A callback to the hit-test function.
  INCLUDE 'intf_mshreghittest.inc'
  
  ! OPTIONAL: 
  INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)    :: IallowedHit
  
!</input>

!<output>
  ! The mesh region of the triangulation.
  TYPE(t_meshRegion), INTENT(OUT)                :: rmeshRegion
!</output>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional information
  ! to the hit-test routine.
  TYPE(t_collection), OPTIONAL, INTENT(INOUT)    :: rcollection
!</inputoutput>

!</subroutine>

  ! Some local variables
  INTEGER :: idim,i,idx
  INTEGER, DIMENSION(:), POINTER :: p_Ihit, p_Iindex
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertsAtCell
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dcoords, p_Dverts
  
    ! Hang in the triangulation
    rmeshRegion%p_rtriangulation => rtriangulation
    
    ! Get the dimension of the triangulation
    idim = rtriangulation%ndim
    
    ! Get the vertice coordinates
    CALL storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_Dverts)
    
    ! Do we calculate vertices?
    IF (IAND(cidxCalc, MSHREG_IDX_VERTEX) .NE. 0) THEN
    
      ! Do we process only boundary vertices?
      IF (bonlyBndry) THEN
      
        ! Allocate a map
        ALLOCATE(p_Ihit(rtriangulation%NVBD))
        
        ! And call the hit-test routine
        IF (PRESENT(rcollection)) THEN
          CALL fmshregHitTest(rtriangulation%NVBD, p_Dverts, p_Ihit, &
                              rcollection)
        ELSE
          CALL fmshregHitTest(rtriangulation%NVBD, p_Dverts, p_Ihit)
        END IF
        
        ! Get the vertices-at-boundary array
        CALL storage_getbase_int(rtriangulation%h_IverticesAtBoundary, p_Iindex)
        
        ! Calculate the index array
        IF (PRESENT(IallowedHit)) THEN
          CALL mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NVBD, &
             rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, p_Iindex, IallowedHit)
        ELSE
          CALL mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NVBD, &
             rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, p_Iindex)
        END IF
        
        ! Free the hit-map
        DEALLOCATE(p_Ihit)
      
      ELSE
      
        ! Allocate a map
        ALLOCATE(p_Ihit(rtriangulation%NVT))
        
        ! And call the hit-test routine
        IF (PRESENT(rcollection)) THEN
          CALL fmshregHitTest(rtriangulation%NVT, p_Dverts, p_Ihit, &
                              rcollection)
        ELSE
          CALL fmshregHitTest(rtriangulation%NVT, p_Dverts, p_Ihit)
        END IF
        
        ! Calculate the index array
        IF (PRESENT(IallowedHit)) THEN
          CALL mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NVT, &
             rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0, IallowedHit)
        ELSE
          CALL mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NVT, &
             rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0)
        END IF
        
        ! Free the hit-map
        DEALLOCATE(p_Ihit)

      END IF
    
    END IF

    ! Do we calculate edges?
    IF ((idim .GE. 2) .AND. (IAND(cidxCalc, MSHREG_IDX_EDGE) .NE. 0)) THEN
    
      ! Do we process only boundary edgess?
      IF (bonlyBndry) THEN
      
        ! Allocate a map
        ALLOCATE(p_Ihit(rtriangulation%NMBD))
        
        ! Allocate a coordinate array
        ALLOCATE(p_Dcoords(idim,rtriangulation%NMBD))

        ! Get the edges-at-boundary array
        CALL storage_getbase_int(rtriangulation%h_IedgesAtBoundary, p_Iindex)
        
        ! Get the vertices-at-edge array
        CALL storage_getbase_int2D(rtriangulation%h_IverticesAtEdge, &
                                   p_IvertsAtCell)
        
        ! Calculate the edge mid-points
        DO i=1, rtriangulation%NMBD
          
          ! Get the index of the edge
          idx = p_Iindex(i)
          
          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.5_DP * (p_Dverts(:, p_IvertsAtCell(1,idx)) &
                                  +  p_Dverts(:, p_IvertsAtCell(2,idx)))
        END DO
        
        ! And call the hit-test routine
        IF (PRESENT(rcollection)) THEN
          CALL fmshregHitTest(rtriangulation%NMBD, p_Dcoords, p_Ihit, &
                              rcollection)
        ELSE
          CALL fmshregHitTest(rtriangulation%NMBD, p_Dcoords, p_Ihit)
        END IF
        
        
        ! Calculate the index array
        IF (PRESENT(IallowedHit)) THEN
          CALL mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NMBD, &
             rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, p_Iindex, IallowedHit)
        ELSE
          CALL mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NMBD, &
             rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, p_Iindex)
        END IF
        
        ! Deallocate midpoint array
        DEALLOCATE(p_Dcoords)
        
        ! Free the hit-map
        DEALLOCATE(p_Ihit)
      
      ELSE
      
        ! Allocate a coordinate array
        ALLOCATE(p_Dcoords(idim,rtriangulation%NMT))

        ! Get the vertices-at-edge array
        CALL storage_getbase_int2D(rtriangulation%h_IverticesAtEdge, &
                                   p_IvertsAtCell)
        
        ! Calculate the edge mid-points
        DO i=1, rtriangulation%NMBD
          
          ! Get the index of the edge
          idx = rtriangulation%NVT + i
          
          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.5_DP * (p_Dverts(:, p_IvertsAtCell(1,idx)) &
                                  +  p_Dverts(:, p_IvertsAtCell(2,idx)))
        END DO
        
        ! And call the hit-test routine
        IF (PRESENT(rcollection)) THEN
          CALL fmshregHitTest(rtriangulation%NMT, p_Dcoords, p_Ihit, &
                              rcollection)
        ELSE
          CALL fmshregHitTest(rtriangulation%NMT, p_Dcoords, p_Ihit)
        END IF
        
        ! Calculate the index array
        IF (PRESENT(IallowedHit)) THEN
          CALL mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NMT, &
             rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, 0, IallowedHit)
        ELSE
          CALL mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NMT, &
             rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, 0)
        END IF
        
        ! Deallocate midpoint array
        DEALLOCATE(p_Dcoords)
        
        ! Free the hit-map
        DEALLOCATE(p_Ihit)
      
      END IF
    
    END IF

    ! Do we calculate face?
    IF ((idim .GE. 3) .AND. (IAND(cidxCalc, MSHREG_IDX_FACE) .NE. 0)) THEN
    
      ! Do we process only boundary edgess?
      IF (bonlyBndry) THEN
      
        ! Allocate a map
        ALLOCATE(p_Ihit(rtriangulation%NABD))
        
        ! Allocate a coordinate array
        ALLOCATE(p_Dcoords(idim,rtriangulation%NABD))

        ! Get the faces-at-boundary array
        CALL storage_getbase_int(rtriangulation%h_IfacesAtBoundary, p_Iindex)
        
        ! Get the vertices-at-face array
        CALL storage_getbase_int2D(rtriangulation%h_IverticesAtFace, &
                                   p_IvertsAtCell)
        
        ! Calculate the face mid-points
        DO i=1, rtriangulation%NABD
          
          ! Get the index of the face
          idx = p_Iindex(i)
          
          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.25_DP * (p_Dverts(:, p_IvertsAtCell(1,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(2,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(3,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(4,idx)))
        END DO
        
        ! And call the hit-test routine
        IF (PRESENT(rcollection)) THEN
          CALL fmshregHitTest(rtriangulation%NABD, p_Dcoords, p_Ihit, &
                              rcollection)
        ELSE
          CALL fmshregHitTest(rtriangulation%NABD, p_Dcoords, p_Ihit)
        END IF
        
        
        ! Calculate the index array
        IF (PRESENT(IallowedHit)) THEN
          CALL mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NABD, &
             rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, p_Iindex, IallowedHit)
        ELSE
          CALL mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NABD, &
             rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, p_Iindex)
        END IF
        
        ! Deallocate midpoint array
        DEALLOCATE(p_Dcoords)
        
        ! Free the hit-map
        DEALLOCATE(p_Ihit)
      
      ELSE
      
        ! Allocate a coordinate array
        ALLOCATE(p_Dcoords(idim,rtriangulation%NMT))

        ! Get the vertices-at-face array
        CALL storage_getbase_int2D(rtriangulation%h_IverticesAtEdge, &
                                   p_IvertsAtCell)
        
        ! Calculate the face mid-points
        DO i=1, rtriangulation%NMBD
          
          ! Get the index of the face
          idx = rtriangulation%NVT + rtriangulation%NMT + i
          
          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.25_DP * (p_Dverts(:, p_IvertsAtCell(1,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(2,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(3,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(4,idx)))
        END DO
        
        ! And call the hit-test routine
        IF (PRESENT(rcollection)) THEN
          CALL fmshregHitTest(rtriangulation%NAT, p_Dcoords, p_Ihit, &
                              rcollection)
        ELSE
          CALL fmshregHitTest(rtriangulation%NAT, p_Dcoords, p_Ihit)
        END IF
        
        ! Calculate the index array
        IF (PRESENT(IallowedHit)) THEN
          CALL mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NAT, &
             rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, 0, IallowedHit)
        ELSE
          CALL mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NAT, &
             rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, 0)
        END IF
        
        ! Deallocate midpoint array
        DEALLOCATE(p_Dcoords)
        
        ! Free the hit-map
        DEALLOCATE(p_Ihit)
      
      END IF
    
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
        (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IvertexIdx)
        
    END IF
  
    ! If the region does not share its edge index array, release it
    IF ((rmeshRegion%h_IedgeIdx .NE. ST_NOHANDLE) .AND. &
        (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_EDGE) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IedgeIdx)
        
    END IF
  
    ! If the region does not share its face index array, release it
    IF ((rmeshRegion%h_IfaceIdx .NE. ST_NOHANDLE) .AND. &
        (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_FACE) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IfaceIdx)
        
    END IF
  
    ! If the region does not share its element index array, release it
    IF ((rmeshRegion%h_IelementIdx .NE. ST_NOHANDLE) .AND. &
        (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_ELEMENT) .EQ. 0)) THEN
        
        CALL storage_free(rmeshRegion%h_IelementIdx)
        
    END IF
    
    ! Last but not least: reset the structure
    rmeshRegion%p_rtriangulation => NULL()
    rmeshRegion%cshareFlags = 0
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

  SUBROUTINE mshreg_recalcVerticesFromEdges(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the vertice index array from the edge index array.
  ! The operator mask parameter decides which vertices will be added into
  ! the mesh region's vertice index array:
  ! 1. coperatorMask = MSHREG_MASK_KICK:
  !    The old vertice index array is deleted (if it exists), and every
  !    vertice that belongs to an edge in the mesh region will be added to
  !    the vertice index array.
  !
  ! 2. coperatorMask = MSHREG_MASK_OR:
  ! -> Every vertice that already exists in the mesh region or that belongs
  !    to an edge in the mesh region will be added to vertice index array.
  !
  ! 3. coperatorMask = MSHREG_MASK_AND:
  ! -> Every vertice that already exists in the mesh region and that belongs
  !    to an edge in the mesh region will be added to vertice index array.
  !    This implies that if the mesh region's vertice index array is empty
  !    on entry, it will also be empty on exit.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old vertice index array.
  ! One of the MSHREG_MASK_XXXX constants defined above.
  ! If not given, MSHREG_MASK_KICK is used.
  INTEGER(I32), INTENT(IN), OPTIONAL   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  TYPE(t_meshRegion), INTENT(INOUT)    :: rmeshRegion
!</inputoutput>

!</subroutine>
  
  INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Imap
  INTEGER, DIMENSION(:), POINTER :: p_IedgeIdx, p_IvertIdx
  INTEGER, DIMENSION(:,:), POINTER :: p_IvertsAtEdge
  INTEGER, DIMENSION(1) :: IallowedProp
  INTEGER(I32) :: copMask
  INTEGER :: i, iedge
  TYPE(t_triangulation), POINTER :: p_rtria

    ! Decide what operator to use
    copMask = MSHREG_MASK_KICK
    IF (PRESENT(coperatorMask)) copMask = coperatorMask
    
    ! If we use the AND-operator, the allowed property is 2, otherwise 1
    IF (copMask .EQ. MSHREG_MASK_AND) THEN
      IallowedProp(1) = 2
    ELSE
      IallowedProp(1) = 1
    END IF
    
    ! Do we have any edges at all?
    IF ((rmeshRegion%NMT .LE. 0) .OR. (rmeshRegion%h_IedgeIdx .EQ. &
         ST_NOHANDLE)) THEN
    
      ! If the mesh region already had an vertex index array and we do not
      ! apply the OR-operator, kick the old vertex index array.
      IF (copMask .NE. MSHREG_MASK_OR) THEN
      
        IF ((rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) .AND. &
            (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .EQ. 0)) THEN
            
            CALL storage_free(rmeshRegion%h_IvertexIdx)
            
        END IF
        
        ! Remove anything that has to do with vertices from the mesh region
        rmeshRegion%NVT = 0
        rmeshRegion%h_IvertexIdx = ST_NOHANDLE
        rmeshRegion%cshareFlags = IAND(rmeshRegion%cshareFlags, &
                                       NOT(MSHREG_IDX_VERTEX))
      
      END IF
      
      ! Return here
      RETURN
    
    END IF
    
    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation
    
    ! Get the edge index array
    CALL storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)
    
    ! Get the vertices-at-edge array from the triangulation
    CALL storage_getbase_int2D(p_rtria%h_IverticesAtEdge, p_IvertsAtEdge)
    
    ! Allocate an vertex map
    ALLOCATE(Imap(p_rtria%NVT))
    DO i=1, p_rtria%NVT
      Imap(i) = 0
    END DO
    
    ! Go through all edges
    DO i=1, rmeshRegion%NMT
    
      ! Get the edge index
      iedge = p_IedgeIdx(i)
      
      ! Go through all vertices adjacent to that edge
      Imap(p_IvertsAtEdge(1,iedge)) = 1
      Imap(p_IvertsAtEdge(2,iedge)) = 1
      
    END DO
    
    ! Now let's see what to do with the old vertices
    IF (rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) THEN
    
      ! Get the old vertices index array
      CALL storage_getbase_int(rmeshRegion%h_IvertexIdx, p_IvertIdx)
      
      SELECT CASE(copMask)
      CASE (MSHREG_MASK_OR)
        ! Go through all vertices and set the map entry to 1
        DO i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = 1
        END DO
      
      CASE (MSHREG_MASK_AND)
        ! Go through all vertices and add 1 to the map entry
        DO i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = Imap(p_IvertIdx(i)) + 1
        END DO
        
        ! Now all vertices that have already been in the mesh region
        ! and which are adjacent to an edge in the mesh region have
        ! a value of 2 in the map.

      END SELECT
      
      ! Now let's destroy the old vertex index array
      IF ((rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) .AND. &
          (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .EQ. 0)) THEN
          
          CALL storage_free(rmeshRegion%h_IvertexIdx)
          
      END IF
      
      ! Remove anything that has to do with vertices from the mesh region
      rmeshRegion%NVT = 0
      rmeshRegion%h_IvertexIdx = ST_NOHANDLE
      rmeshRegion%cshareFlags = IAND(rmeshRegion%cshareFlags, &
                                     NOT(MSHREG_IDX_VERTEX))
      
    END IF
    
    ! Now we have the vertex index map, so create an index array from this.
    CALL mshreg_aux_calcIdxArray1(Imap, 0, p_rtria%NVT, &
        rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0, IallowedProp)
       
    ! And release the vertex map
    DEALLOCATE(Imap)
    
    ! That's it
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_recalcVerticesFromFaces(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the vertice index array from the face index array.
  ! The operator mask parameter decides which vertices will be added into
  ! the mesh region's vertice index array:
  ! 1. coperatorMask = MSHREG_MASK_KICK:
  !    The old vertice index array is deleted (if it exists), and every
  !    vertice that belongs to a face in the mesh region will be added to
  !    the vertice index array.
  !
  ! 2. coperatorMask = MSHREG_MASK_OR:
  ! -> Every vertice that already exists in the mesh region or that belongs
  !    to a face in the mesh region will be added to vertice index array.
  !
  ! 3. coperatorMask = MSHREG_MASK_AND:
  ! -> Every vertice that already exists in the mesh region and that belongs
  !    to a face in the mesh region will be added to vertice index array.
  !    This implies that if the mesh region's vertice index array is empty
  !    on entry, it will also be empty on exit.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old vertice index array.
  ! One of the MSHREG_MASK_XXXX constants defined above.
  ! If not given, MSHREG_MASK_KICK is used.
  INTEGER(I32), INTENT(IN), OPTIONAL   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  TYPE(t_meshRegion), INTENT(INOUT)    :: rmeshRegion
!</inputoutput>

!</subroutine>
  
  INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Imap
  INTEGER, DIMENSION(:), POINTER :: p_IfaceIdx, p_IvertIdx
  INTEGER, DIMENSION(:,:), POINTER :: p_IvertsAtFace
  INTEGER, DIMENSION(1) :: IallowedProp
  INTEGER(I32) :: copMask
  INTEGER :: i, j, iface
  TYPE(t_triangulation), POINTER :: p_rtria

    ! Decide what operator to use
    copMask = MSHREG_MASK_KICK
    IF (PRESENT(coperatorMask)) copMask = coperatorMask

    ! If we use the AND-operator, the allowed property is 2, otherwise 1
    IF (copMask .EQ. MSHREG_MASK_AND) THEN
      IallowedProp(1) = 2
    ELSE
      IallowedProp(1) = 1
    END IF
    
    ! Do we have any faces at all?
    IF ((rmeshRegion%NAT .LE. 0) .OR. (rmeshRegion%h_IfaceIdx .EQ. &
         ST_NOHANDLE)) THEN
    
      ! If the mesh region already had an vertex index array and we do not
      ! apply the OR-operator, kick the old vertex index array.
      IF (copMask .NE. MSHREG_MASK_OR) THEN
      
        IF ((rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) .AND. &
            (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .EQ. 0)) THEN
            
            CALL storage_free(rmeshRegion%h_IvertexIdx)
            
        END IF
        
        ! Remove anything that has to do with vertices from the mesh region
        rmeshRegion%NVT = 0
        rmeshRegion%h_IvertexIdx = ST_NOHANDLE
        rmeshRegion%cshareFlags = IAND(rmeshRegion%cshareFlags, &
                                  NOT(MSHREG_IDX_VERTEX))
      
      END IF
      
      ! Return here
      RETURN
    
    END IF
    
    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation
    
    ! Get the face index array
    CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! Get the vertices-at-face array from the triangulation
    CALL storage_getbase_int2D(p_rtria%h_IverticesAtFace,&
                               p_IvertsAtFace)
    
    ! Allocate an vertex map
    ALLOCATE(Imap(p_rtria%NVT))
    DO i=1, p_rtria%NVT
      Imap(i) = 0
    END DO
    
    ! Go through all faces
    DO i=1, rmeshRegion%NAT
    
      ! Get the face index
      iface = p_IfaceIdx(i)
      
      ! Go through all edges adjacent to that face
      DO j=1, UBOUND(p_IvertsAtFace,1)
        
        IF (p_IvertsAtFace(j,iface) .GT. 0) THEN
          Imap(p_IvertsAtFace(j,iface)) = 1
        END IF
        
      END DO
      
    END DO
    
    ! Now let's see what to do with the old vertices
    IF (rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) THEN
    
      ! Get the old vertices index array
      CALL storage_getbase_int(rmeshRegion%h_IvertexIdx, p_IvertIdx)
      
      SELECT CASE(copMask)
      CASE (MSHREG_MASK_OR)
        ! Go through all vertices and set the map entry to 1
        DO i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = 1
        END DO
      
      CASE (MSHREG_MASK_AND)
        ! Go through all vertices and add 1 to the map entry
        DO i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = Imap(p_IvertIdx(i)) + 1
        END DO
        
        ! Now all vertices that have already been in the mesh region
        ! and which are adjacent to a face in the mesh region have
        ! a value of 2 in the map.
        
      END SELECT
      
      ! Now let's destroy the old vertex index array
      IF ((rmeshRegion%h_IvertexIdx .NE. ST_NOHANDLE) .AND. &
          (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .EQ. 0)) THEN
          
          CALL storage_free(rmeshRegion%h_IvertexIdx)
          
      END IF
      
      ! Remove anything that has to do with vertices from the mesh region
      rmeshRegion%NVT = 0
      rmeshRegion%h_IvertexIdx = ST_NOHANDLE
      rmeshRegion%cshareFlags = IAND(rmeshRegion%cshareFlags, &
                                     NOT(MSHREG_IDX_VERTEX))
      
    END IF
    
    ! Now we have the vertex index map, so create an index array from this.
    CALL mshreg_aux_calcIdxArray1(Imap, 0, p_rtria%NVT, &
        rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0, IallowedProp)
    
    ! And release the vertex map
    DEALLOCATE(Imap)
    
    ! That's it
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_recalcEdgesFromFaces(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the edge index array from the face index array.
  ! The operator mask parameter decides which edges will be added into
  ! the mesh region's edge index array:
  ! 1. coperatorMask = MSHREG_MASK_KICK:
  !    The old edge index array is deleted (if it exists), and every
  !    edge that belongs to a face in the mesh region will be added to
  !    the edge index array.
  !
  ! 2. coperatorMask = MSHREG_MASK_OR:
  ! -> Every edge that already exists in the mesh region or that belongs
  !    to a face in the mesh region will be added to edge index array.
  !
  ! 3. coperatorMask = MSHREG_MASK_AND:
  ! -> Every edge that already exists in the mesh region and that belongs
  !    to a face in the mesh region will be added to edge index array.
  !    This implies that if the mesh region's edge index array is empty
  !    on entry, it will also be empty on exit.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old edge index array.
  ! One of the MSHREG_MASK_XXXX constants defined above.
  ! If not given, MSHREG_MASK_KICK is used.
  INTEGER(I32), INTENT(IN), OPTIONAL   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  TYPE(t_meshRegion), INTENT(INOUT)    :: rmeshRegion
!</inputoutput>

!</subroutine>
  
  INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Imap
  INTEGER, DIMENSION(:), POINTER :: p_IfaceIdx, p_IedgeIdx
  INTEGER, DIMENSION(:,:), POINTER :: p_IedgesAtFace
  INTEGER, DIMENSION(1) :: IallowedProp
  INTEGER(I32) :: copMask
  INTEGER :: i, j, iface
  TYPE(t_triangulation), POINTER :: p_rtria

    ! Decide what operator to use
    copMask = MSHREG_MASK_KICK
    IF (PRESENT(coperatorMask)) copMask = coperatorMask

    ! If we use the AND-operator, the allowed property is 2, otherwise 1
    IF (copMask .EQ. MSHREG_MASK_AND) THEN
      IallowedProp(1) = 2
    ELSE
      IallowedProp(1) = 1
    END IF
    
    ! Do we have any faces at all?
    IF ((rmeshRegion%NAT .LE. 0) .OR. (rmeshRegion%h_IfaceIdx .EQ. &
         ST_NOHANDLE)) THEN
    
      ! If the mesh region already had an edge index array and we do not
      ! apply the OR-operator, kick the old edge index array.
      IF (copMask .NE. MSHREG_MASK_OR) THEN
      
        IF ((rmeshRegion%h_IedgeIdx .NE. ST_NOHANDLE) .AND. &
            (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_EDGE) .EQ. 0)) THEN
            
            CALL storage_free(rmeshRegion%h_IedgeIdx)
            
        END IF
        
        ! Remove anything that has to do with edges from the mesh region
        rmeshRegion%NMT = 0
        rmeshRegion%h_IedgeIdx = ST_NOHANDLE
        rmeshRegion%cshareFlags = IAND(rmeshRegion%cshareFlags, &
                                       NOT(MSHREG_IDX_EDGE))
      
      END IF
      
      ! Return here
      RETURN
    
    END IF
    
    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation
    
    ! Get the face index array
    CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! Get the edges-at-face array from the triangulation
    CALL storage_getbase_int2D(p_rtria%h_IedgesAtFace, p_IedgesAtFace)
    
    ! Allocate an edge map
    ALLOCATE(Imap(p_rtria%NMT))
    DO i = 1, p_rtria%NMT
      Imap(i) = 0
    END DO
    
    ! Go through all faces
    DO i=1, rmeshRegion%NAT
    
      ! Get the face index
      iface = p_IfaceIdx(i)
      
      ! Go through all edges adjacent to that face
      DO j=1, UBOUND(p_IedgesAtFace,1)
        
        IF (p_IedgesAtFace(j,iface) .GT. 0) THEN
          Imap(p_IedgesAtFace(j,iface) - p_rtria%NVT) = 1
        END IF
        
      END DO
      
    END DO
    
    ! Now let's see what to do with the old edges
    IF (rmeshRegion%h_IedgeIdx .NE. ST_NOHANDLE) THEN
    
      ! Get the old edge index array
      CALL storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)
      
      SELECT CASE(copMask)
      CASE (MSHREG_MASK_OR)
        ! Go through all edges and set the map entry to 1
        DO i=1, rmeshRegion%NMT
          Imap(p_IedgeIdx(i)) = 1
        END DO
      
      CASE (MSHREG_MASK_AND)
        ! Go through all edges and add 1 to the map entry
        DO i=1, rmeshRegion%NMT
          Imap(p_IedgeIdx(i)) = Imap(p_IedgeIdx(i)) + 1
        END DO

        ! Now all edges that have already been in the mesh region
        ! and which are adjacent to a face in the mesh region have
        ! a value of 2 in the map.
        
      END SELECT
      
      ! Now let's destroy the old edge index array
      IF (IAND(rmeshRegion%cshareFlags, MSHREG_IDX_EDGE) .EQ. 0) THEN
          CALL storage_free(rmeshRegion%h_IedgeIdx)
      END IF
      
      ! Remove anything that has to do with edges from the mesh region
      rmeshRegion%NMT = 0
      rmeshRegion%h_IedgeIdx = ST_NOHANDLE
      rmeshRegion%cshareFlags = IAND(rmeshRegion%cshareFlags, &
                                     NOT(MSHREG_IDX_EDGE))
      
    END IF
    
    ! Now we have the edge index map, so create an index array from this.
    CALL mshreg_aux_calcIdxArray1(Imap, 0, p_rtria%NMT, &
        rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, 0, IallowedProp)
    
    ! And release the edge map
    DEALLOCATE(Imap)
    
    ! That's it
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_calcBoundaryNormals2D(rmeshRegion,Dnormals)

!<description>
  ! Calculates the inner normal vectors for the edges of a 2D mesh
  ! region, i.e. the normal vectors which point "into the domain".
  ! It is silently assumed that all edges of the mesh region lie on
  ! the boundary of the mesh.
!</description>

!<input>
  ! The mesh region for whose edges the normals are to be calcuated.
  TYPE(t_meshRegion), INTENT(IN)        :: rmeshRegion
!</input>

!<output>
  ! An array that recieves the normal vectors of the edges.
  ! Dimension: Dnormals(NDIM2D,rmeshRegion%NMT)
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dnormals
!</output>

!</subroutine>

  ! Some local variables
  TYPE(t_triangulation), POINTER :: p_rtria
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dcoords
  INTEGER, DIMENSION(:), POINTER :: p_IedgeIdx
  INTEGER, DIMENSION(:,:), POINTER :: p_IvertAtEdge
  REAL(DP) :: dx, dy, dt
  INTEGER :: imt, iedge
  
  
    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation
    
    ! Is this a 2D triangulation?
    IF (p_rtria%ndim .NE. 2) THEN
      PRINT *, 'ERROR: mshreg_calcBoundaryNormals2D'
      PRINT *, 'Triangulation must be 2D'
      CALL sys_halt()
    END IF
    
    ! Don't we have any edges here?
    IF (rmeshRegion%NMT .LE. 0) RETURN
    
    ! In 2D the task of calculating the normal vectors is quite easy,
    ! if we exploit some knowledge about the algorithm which calculates
    ! the edges in a 2D mesh.
    ! As we know that all edges on the boundary are oriented in mathematically
    ! positive direction, the only thing that we have to do is rotate the
    ! edge vector by 90° counter-clock-wise (and normalise, of course)
    ! to get the inner normal vector of a boundary edge...
        
    ! Get a pointer to the vertices-at-edge array
    CALL storage_getbase_int2D(p_rtria%h_IverticesAtEdge, p_IvertAtEdge)
    
    ! Get the vertice coordinate array
    CALL storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dcoords)
    
    ! And get the edge index array of the mesh region
    CALL storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)
    
    ! Now loop through all edges in the mesh region
    DO imt = 1, rmeshRegion%NMT
      
      ! Get the index of the edge
      iedge = p_IedgeIdx(imt)
      
      ! Calculate edge
      dx = p_Dcoords(1,p_IvertAtEdge(2,iedge)) &
         - p_Dcoords(1,p_IvertAtEdge(1,iedge))
      dy = p_Dcoords(2,p_IvertAtEdge(2,iedge)) &
         - p_Dcoords(2,p_IvertAtEdge(1,iedge))
      
      ! Calculate length of edge
      dt = 1.0_DP / SQRT(dx**2 + dy**2)

      ! Calculate normal
      Dnormals(1,imt) = -dy * dt
      Dnormals(2,imt) =  dx * dt
      
    END DO
    
    ! That's it
  
  END SUBROUTINE
  

!************************************************************************

!<subroutine>

  SUBROUTINE mshreg_calcBoundaryNormals3D(rmeshRegion,Dnormals)

!<description>
  ! Calculates the inner normal vectors for the faces of a 3D mesh
  ! region, i.e. the normal vectors which point "into the domain".
  ! It is silently assumed that all faces of the mesh region lie on
  ! the boundary of the mesh.
!</description>

!<input>
  ! The mesh region for whose edges the normals are to be calcuated.
  TYPE(t_meshRegion), INTENT(IN)        :: rmeshRegion
!</input>

!<output>
  ! An array that recieves the normal vectors of the faces.
  ! Dimension: Dnormals(NDIM3D,rmeshRegion%NAT)
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dnormals
!</output>

!</subroutine>

  ! Some local variables
  TYPE(t_triangulation), POINTER :: p_rtria
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dcoords
  INTEGER, DIMENSION(:), POINTER :: p_IfaceIdx
  INTEGER, DIMENSION(:,:), POINTER :: p_IvertAtElem,&
    p_IfaceAtElem, p_IelemAtFace
  REAL(DP), DIMENSION(3) :: Du,Dv,Dn
  REAL(DP), DIMENSION(3,8) :: Dcorners
  REAL(DP) :: dt
  INTEGER :: iat,ivt,iel,ifae,iface
  
  
    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation
    
    ! Is this a 3D triangulation?
    IF (p_rtria%ndim .NE. 3) THEN
      PRINT *, 'ERROR: mshreg_calcBoundaryNormals3D'
      PRINT *, 'Triangulation must be 3D'
      CALL sys_halt()
    END IF
    
    ! Don't we have any faces here?
    IF (rmeshRegion%NAT .LE. 0) RETURN
    
    ! Get a pointer to the faces-at-element array
    CALL storage_getbase_int2D(p_rtria%h_IfacesAtElement, p_IfaceAtElem)

    ! Get a pointer to the elements-at-face array
    CALL storage_getbase_int2D(p_rtria%h_IelementsAtFace, p_IelemAtface)
       
    ! Get a pointer to the vertices-at-element array
    CALL storage_getbase_int2D(p_rtria%h_IverticesAtElement, p_IvertAtElem)
    
    ! Get the vertice coordinate array
    CALL storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dcoords)
    
    ! And get the edge index array of the mesh region
    CALL storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
    
    ! Now loop through all face in the mesh region
    DO iat = 1, rmeshRegion%NAT
      
      ! Get the index of the face
      iface = p_IfaceIdx(iat)
      
      ! Get the index of the element that is adjacent to that face
      iel = p_IelemAtFace(1,iface)
      
      ! Now go through all faces of the element and search for this one
      DO ifae = 1, 6
        IF (p_IfaceAtElem(ifae,iel) .EQ. iface) EXIT
      END DO
      
      ! Get the eight corner vertices of the hexahedron
      DO ivt = 1, 8
        Dcorners(1:3,ivt) = p_Dcoords(1:3,p_IvertAtElem(ivt,iel))
      END DO
      
      ! Calculate the normal of the face
      SELECT CASE(ifae)
      CASE (1)
        ! (1->3) x (4->2)
        Du(:) = Dcorners(:,3) - Dcorners(:,1)
        Dv(:) = Dcorners(:,2) - Dcorners(:,4)
      CASE (2)
        ! (1->6) x (2->5)
        Du(:) = Dcorners(:,6) - Dcorners(:,1)
        Dv(:) = Dcorners(:,5) - Dcorners(:,2)
      CASE (3)
        ! (2->7) x (6->3)
        Du(:) = Dcorners(:,7) - Dcorners(:,2)
        Dv(:) = Dcorners(:,3) - Dcorners(:,6)
      CASE (4)
        ! (3->8) x (4->7)
        Du(:) = Dcorners(:,8) - Dcorners(:,3)
        Dv(:) = Dcorners(:,7) - Dcorners(:,4)
      CASE (5)
        ! (4->5) x (1->8)
        Du(:) = Dcorners(:,5) - Dcorners(:,4)
        Dv(:) = Dcorners(:,8) - Dcorners(:,1)
      CASE (6)
        ! (5->7) x (6->8)
        Du(:) = Dcorners(:,7) - Dcorners(:,5)
        Dv(:) = Dcorners(:,8) - Dcorners(:,6)
      END SELECT
      
      ! Calculate normal by 3D cross product
      Dn(1) = Du(2)*Dv(3) - Du(3)*Dv(2)
      Dn(2) = Du(3)*Dv(1) - Du(1)*Dv(3)
      Dn(3) = Du(1)*Dv(2) - Du(2)*Dv(1)
      
      ! Normalise it
      dt = 1.0_DP / SQRT(Dn(1)**2 + Dn(2)**2 + Dn(3)**2)
      Dn = dt * Dn
      
      ! Store the normal
      IF (gaux_isFlipped_hexa3D(Dcorners)) THEN
        Dnormals(:,iat) = Dn(:)
      ELSE
        Dnormals(:,iat) = -Dn(:)
      END IF
      
    END DO
    
    ! That's it
  
  END SUBROUTINE
  
END MODULE
