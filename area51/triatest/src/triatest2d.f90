!##############################################################################
!# ****************************************************************************
!# <name> triatest2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module reads in a 2D mesh, optionally refines it a few times and
!# creates a standard mesh from it. Afterwards (nearly) all arrays from the
!# triangulation are printed onto the screen and to the log file.
!# </purpose>
!##############################################################################
MODULE triatest2d

  USE fsystem
  USE genoutput
  USE boundary
  USE triangulation
  USE storage

  IMPLICIT NONE

CONTAINS

  SUBROUTINE tria_test2d(NLMAX)
  INTEGER, OPTIONAL, INTENT(IN) :: NLMAX
  
  TYPE(t_boundary) :: rbnd
  TYPE(t_triangulation) :: rtria
  INTEGER, DIMENSION(:), POINTER :: p_Idata,p_Iidx
  INTEGER, DIMENSION(:,:), POINTER :: p_Idata2
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata2
  INTEGER :: i,j
  
    ! Print a header
    CALL output_separator(OU_SEP_STAR)
    CALL output_line('T R I A _ T E S T 2 D')
    CALL output_separator(OU_SEP_STAR)
  
    ! Read in the boundary parametrisation
    CALL boundary_read_prm(rbnd, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rtria, './pre/QUAD.tri', rbnd)
    
    ! Optionally refine the mesh
    IF(PRESENT(NLMAX)) THEN
      CALL tria_quickRefine2LevelOrdering (NLMAX-1,rtria,rbnd)
    END IF

    ! Standardise the mesh
    CALL tria_initStandardMeshFromRaw (rtria, rbnd)
    
    ! First of all, print out the basic information
    CALL output_line('Triangulation Scalars')
    CALL output_line('----------------------------------')
    CALL output_line('NVT....................: ' // TRIM(sys_si(rtria%NVT,8)))
    CALL output_line('NMT....................: ' // TRIM(sys_si(rtria%NMT,8)))
    CALL output_line('NEL....................: ' // TRIM(sys_si(rtria%NEL,8)))
    CALL output_line('NBCT...................: ' // TRIM(sys_si(rtria%NBCT,8)))
    CALL output_line('NblindBCT..............: ' // TRIM(sys_si(rtria%NblindBCT,8)))
    CALL output_line('NVBD...................: ' // TRIM(sys_si(rtria%NVBD,8)))
    CALL output_line('NMBD...................: ' // TRIM(sys_si(rtria%NMBD,8)))
    CALL output_line('NNVE...................: ' // TRIM(sys_si(rtria%NNVE,8)))
    CALL output_line('NNEE...................: ' // TRIM(sys_si(rtria%NNEE,8)))
    CALL output_line('NNelAtVertex...........: ' // TRIM(sys_si(rtria%NNelAtVertex,8)))
    CALL output_line('NNelAtEdge.............: ' // TRIM(sys_si(rtria%NNelAtEdge,8)))
    CALL output_line('nverticesPerEdge.......: ' // TRIM(sys_si(rtria%nverticesPerEdge,8)))
    CALL output_line('nverticesOnAllEdges....: ' // TRIM(sys_si(rtria%nVerticesOnAllEdges,8)))
    CALL output_line('nverticesInEachElement.: ' // TRIM(sys_si(rtria%nverticesInEachElement,8)))
    CALL output_line('nverticesInAllElements.: ' // TRIM(sys_si(rtria%nverticesInAllElements,8)))
    CALL output_line('nadditionalVertices....: ' // TRIM(sys_si(rtria%NMT,8)))
    
    ! Print the vertice coordinates
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_double2d(rtria%h_DvertexCoords, p_Ddata2)
    CALL output_line('p_DvertexCoords')
    DO i = LBOUND(p_Ddata2,2), UBOUND(p_Ddata2,2)
      CALL output_line(TRIM(sys_si(i,8)) // ': (/ ' // TRIM(sys_sdEP(p_Ddata2(1,i),20,12)) //&
              ' , ' // TRIM(sys_sdEP(p_Ddata2(2,i),20,12)) // ' /)')
    END DO
    
    ! Print vertices-at-Element
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_int2d(rtria%h_IverticesAtElement,p_Idata2)
    CALL output_line('p_IverticesAtElement')
    DO i = LBOUND(p_Idata2,2), UBOUND(p_Idata2,2)
      CALL output_line(TRIM(sys_si(i,8)) // ': (/ ', bnolinebreak = .TRUE.)
      DO j = LBOUND(p_Idata2,1), UBOUND(p_Idata2,1)-1
        CALL output_line(TRIM(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .TRUE.)
      END DO
      CALL output_line(TRIM(sys_si(p_Idata2(UBOUND(p_Idata2,1),i),8)) // ' /)')
    END DO
    
    ! Print edges-at-Element
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_int2d(rtria%h_IedgesAtElement,p_Idata2)
    CALL output_line('p_IedgesAtElement')
    DO i = LBOUND(p_Idata2,2), UBOUND(p_Idata2,2)
      CALL output_line(TRIM(sys_si(i,8)) // ': (/ ', bnolinebreak = .TRUE.)
      DO j = LBOUND(p_Idata2,1), UBOUND(p_Idata2,1)-1
        CALL output_line(TRIM(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .TRUE.)
      END DO
      CALL output_line(TRIM(sys_si(p_Idata2(UBOUND(p_Idata2,1),i),8)) // ' /)')
    END DO
    
    ! Print neighbours-at-Element
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_int2d(rtria%h_IneighboursAtElement,p_Idata2)
    CALL output_line('h_IneighboursAtElement')
    DO i = LBOUND(p_Idata2,2), UBOUND(p_Idata2,2)
      CALL output_line(TRIM(sys_si(i,8)) // ': (/ ', bnolinebreak = .TRUE.)
      DO j = LBOUND(p_Idata2,1), UBOUND(p_Idata2,1)-1
        CALL output_line(TRIM(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .TRUE.)
      END DO
      CALL output_line(TRIM(sys_si(p_Idata2(UBOUND(p_Idata2,1),i),8)) // ' /)')
    END DO
    
    ! Print vertices-at-edge
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_int2d(rtria%h_IverticesAtEdge,p_Idata2)
    CALL output_line('h_IverticesAtEdge')
    DO i = LBOUND(p_Idata2,2), UBOUND(p_Idata2,2)
      CALL output_line(TRIM(sys_si(i,8)) // ': (/ ', bnolinebreak = .TRUE.)
      DO j = LBOUND(p_Idata2,1), UBOUND(p_Idata2,1)-1
        CALL output_line(TRIM(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .TRUE.)
      END DO
      CALL output_line(TRIM(sys_si(p_Idata2(UBOUND(p_Idata2,1),i),8)) // ' /)')
    END DO

    ! Print elements-at-edge
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_int2d(rtria%h_IelementsAtEdge,p_Idata2)
    CALL output_line('h_IelementsAtEdge')
    DO i = LBOUND(p_Idata2,2), UBOUND(p_Idata2,2)
      CALL output_line(TRIM(sys_si(i,8)) // ': (/ ', bnolinebreak = .TRUE.)
      DO j = LBOUND(p_Idata2,1), UBOUND(p_Idata2,1)-1
        CALL output_line(TRIM(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .TRUE.)
      END DO
      CALL output_line(TRIM(sys_si(p_Idata2(UBOUND(p_Idata2,1),i),8)) // ' /)')
    END DO

    ! Print edges at vertices
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_IedgesAtVertexIdx .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IedgesAtVertexIdx,p_Iidx)
      CALL storage_getbase_int(rtria%h_IedgesAtVertex,p_Idata)
      CALL output_line('h_IedgesAtVertexIdx / h_IedgesAtVertex')
      DO i = LBOUND(p_Iidx,1), UBOUND(p_Iidx,1)-1
        CALL output_line(TRIM(sys_si(i,8)) // ' , '// TRIM(sys_si(p_Iidx(i),8)) //&
                         ' : ' // TRIM(sys_si(p_Idata(p_Iidx(i)),8)))
        DO j = p_Iidx(i)+1, p_Iidx(i+1)-1
          CALL output_line('           ' // TRIM(sys_si(j,8)) // ' : ' // TRIM(sys_si(p_Idata(j),8)))
        END DO
      END DO
    ELSE
      CALL output_line('h_IedgesAtVertexIdx = ST_NOHANDLE')
    END IF
    
    ! Print elements at vertex
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_int(rtria%h_IelementsAtVertexIdx,p_Iidx)
    CALL storage_getbase_int(rtria%h_IelementsAtVertex,p_Idata)
    CALL output_line('h_IelementsAtVertexIdx / h_IelementsAtVertex')
    DO i = LBOUND(p_Iidx,1), UBOUND(p_Iidx,1)-1
      CALL output_line(TRIM(sys_si(i,8)) // ' , '// TRIM(sys_si(p_Iidx(i),8)) //&
                       ' : ' // TRIM(sys_si(p_Idata(p_Iidx(i)),8)))
      DO j = p_Iidx(i)+1, p_Iidx(i+1)-1
        CALL output_line('           ' // TRIM(sys_si(j,8)) // ' : ' // TRIM(sys_si(p_Idata(j),8)))
      END DO
    END DO
    
    ! Print nodal property
    CALL output_separator(OU_SEP_MINUS)
    CALL storage_getbase_int(rtria%h_InodalProperty,p_Idata)
    CALL output_line('h_InodalProperty')
    DO i = LBOUND(p_Idata,1), UBOUND(p_Idata,1)
      CALL output_line(TRIM(sys_si(i,8)) // ': ' // TRIM(sys_si(p_Idata(i),8)))
    END DO

    ! Print maco nodal property
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_ImacroNodalProperty .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_ImacroNodalProperty,p_Idata)
      CALL output_line('h_ImacroNodalProperty')
      DO i = LBOUND(p_Idata,1), UBOUND(p_Idata,1)
        CALL output_line(TRIM(sys_si(i,8)) // ': ' // TRIM(sys_si(p_Idata(i),8)))
      END DO
    ELSE
      CALL output_line('h_ImacroNodalProperty = ST_NOHANDLE')
    END IF

    ! Print element volume
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_DelementVolume .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_double(rtria%h_DelementVolume,p_Ddata)
      CALL output_line('h_DelementVolume')
      DO i = LBOUND(p_Ddata,1), UBOUND(p_Ddata,1)
        CALL output_line(TRIM(sys_si(i,8)) // ': ' // TRIM(sys_sdEP(p_Ddata(i),20,12)))
      END DO
    ELSE
      CALL output_line('h_DelementVolume = ST_NOHANDLE')
    END IF

    ! Print refinement patch
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_IrefinementPatchIdx .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IrefinementPatchIdx,p_Iidx)
      CALL storage_getbase_int(rtria%h_IrefinementPatch,p_Idata)
      CALL output_line('h_IrefinementPatchIdx / h_IrefinementPatch')
      DO i = LBOUND(p_Iidx,1), UBOUND(p_Iidx,1)-1
        CALL output_line(TRIM(sys_si(i,8)) // ' , '// TRIM(sys_si(p_Iidx(i),8)) //&
                         ' : ' // TRIM(sys_si(p_Idata(p_Iidx(i)),8)))
        DO j = p_Iidx(i)+1, p_Iidx(i+1)-1
          CALL output_line('           ' // TRIM(sys_si(j,8)) // ' : ' // TRIM(sys_si(p_Idata(j),8)))
        END DO
      END DO
    ELSE
      CALL output_line('h_IrefinementPatchIdx = ST_NOHANDLE')
    END IF
    
    ! Print coarse grid element
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_IcoarseGridElement .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IcoarseGridElement,p_Idata)
      CALL output_line('h_IcoarseGridElement')
      DO i = LBOUND(p_Idata,1), UBOUND(p_Idata,1)
        CALL output_line(TRIM(sys_si(i,8)) // ': ' // TRIM(sys_si(p_Idata(i),8)))
      END DO
    ELSE
      CALL output_line('h_IcoarseGridElement = ST_NOHANDLE')
    END IF
    
    ! Print boundary vertices
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_IboundaryCpIdx .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IboundaryCpIdx,p_Iidx)
      CALL storage_getbase_int(rtria%h_IverticesAtBoundary,p_Idata)
      CALL output_line('h_IboundaryCpIdx / h_IverticesAtBoundary')
      DO i = LBOUND(p_Iidx,1), UBOUND(p_Iidx,1)-1
        CALL output_line(TRIM(sys_si(i,8)) // ' , '// TRIM(sys_si(p_Iidx(i),8)) //&
                         ' : ' // TRIM(sys_si(p_Idata(p_Iidx(i)),8)))
        DO j = p_Iidx(i)+1, p_Iidx(i+1)-1
          CALL output_line('           ' // TRIM(sys_si(j,8)) // ' : ' // TRIM(sys_si(p_Idata(j),8)))
        END DO
      END DO
    ELSE
      CALL output_line('h_IboundaryCpIdx = ST_NOHANDLE')
    END IF
    
    ! Print boundary edges
    CALL output_separator(OU_SEP_MINUS)
    !IF (rtria%h_IboundaryCpEdgesIdx .NE. ST_NOHANDLE) THEN
    !  CALL storage_getbase_int(rtria%h_IboundaryCpEdgesIdx,p_Iidx)
    IF (rtria%h_IboundaryCpIdx .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IboundaryCpIdx,p_Iidx)
      CALL storage_getbase_int(rtria%h_IedgesAtBoundary,p_Idata)
      CALL output_line('h_IboundaryCpEdgesIdx / h_IedgesAtBoundary')
      DO i = LBOUND(p_Iidx,1), UBOUND(p_Iidx,1)-1
        CALL output_line(TRIM(sys_si(i,8)) // ' , '// TRIM(sys_si(p_Iidx(i),8)) //&
                         ' : ' // TRIM(sys_si(p_Idata(p_Iidx(i)),8)))
        DO j = p_Iidx(i)+1, p_Iidx(i+1)-1
          CALL output_line('           ' // TRIM(sys_si(j,8)) // ' : ' // TRIM(sys_si(p_Idata(j),8)))
        END DO
      END DO
    ELSE
      CALL output_line('h_IboundaryCpEdgesIdx = ST_NOHANDLE')
    END IF

    ! Print boundary elements
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_IelementsAtBoundary .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rtria%h_IelementsAtBoundary,p_Idata)
      CALL output_line('h_IelementsAtBoundary')
      DO i = LBOUND(p_Idata,1), UBOUND(p_Idata,1)
        CALL output_line(TRIM(sys_si(i,8)) // ': ' // TRIM(sys_si(p_Idata(i),8)))
      END DO
    ELSE
      CALL output_line('h_IelementsAtBoundary = ST_NOHANDLE')
    END IF

    ! Print vertex parameter values
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_DvertexParameterValue .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_double(rtria%h_DvertexParameterValue,p_Ddata)
      CALL output_line('h_DvertexParameterValue')
      DO i = LBOUND(p_Ddata,1), UBOUND(p_Ddata,1)
        CALL output_line(TRIM(sys_si(i,8)) // ': ' // TRIM(sys_sdEP(p_Ddata(i),20,12)))
      END DO
    ELSE
      CALL output_line('h_DvertexParameterValue = ST_NOHANDLE')
    END IF

    ! Print edge parameter values
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_DedgeParameterValue .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_double(rtria%h_DedgeParameterValue,p_Ddata)
      CALL output_line('h_DedgeParameterValue')
      DO i = LBOUND(p_Ddata,1), UBOUND(p_Ddata,1)
        CALL output_line(TRIM(sys_si(i,8)) // ': ' // TRIM(sys_sdEP(p_Ddata(i),20,12)))
      END DO
    ELSE
      CALL output_line('h_DedgeParameterValue = ST_NOHANDLE')
    END IF
    
    ! Print boundary vertex positions
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_IboundaryVertexPos .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int2d(rtria%h_IboundaryVertexPos, p_Idata2)
      CALL output_line('h_IboundaryVertexPos')
      DO i = LBOUND(p_Idata2,2), UBOUND(p_Idata2,2)
        CALL output_line(TRIM(sys_si(i,8)) // ': (/ ' // TRIM(sys_si(p_Idata2(1,i),8)) //&
                ' , ' // TRIM(sys_si(p_Idata2(2,i),8)) // ' /)')
      END DO
    ELSE
      CALL output_line('h_IboundaryVertexPos = ST_NOHANDLE')
    END IF
    
    ! Print boundary edge positions
    CALL output_separator(OU_SEP_MINUS)
    IF (rtria%h_IboundaryEdgePos .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int2d(rtria%h_IboundaryEdgePos, p_Idata2)
      CALL output_line('h_IboundaryEdgePos')
      DO i = LBOUND(p_Idata2,2), UBOUND(p_Idata2,2)
        CALL output_line(TRIM(sys_si(i,8)) // ': (/ ' // TRIM(sys_si(p_Idata2(1,i),8)) //&
                ' , ' // TRIM(sys_si(p_Idata2(2,i),8)) // ' /)')
      END DO
    ELSE
      CALL output_line('h_IboundaryEdgePos = ST_NOHANDLE')
    END IF

    ! Release the mesh
    CALL tria_done(rtria)
    
    ! Release the boundary
    CALL boundary_release(rbnd)
  
  END SUBROUTINE

END MODULE