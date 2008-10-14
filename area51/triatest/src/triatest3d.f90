!##############################################################################
!# ****************************************************************************
!# <name> triatest3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module reads in a 3D mesh, optionally refines it a few times and
!# creates a standard mesh from it. Afterwards (nearly) all arrays from the
!# triangulation are printed onto the screen and to the log file.
!# </purpose>
!##############################################################################
module triatest3d

  use fsystem
  use genoutput
  use boundary
  use triangulation
  use storage

  implicit none

contains

  subroutine tria_test3d(NLMAX)
  integer, optional, intent(IN) :: NLMAX
  
  type(t_triangulation) :: rtria
  integer, dimension(:), pointer :: p_Idata,p_Iidx
  integer, dimension(:,:), pointer :: p_Idata2
  real(DP), dimension(:), pointer :: p_Ddata
  real(DP), dimension(:,:), pointer :: p_Ddata2
  integer :: i,j
  
    ! Print a header
    call output_separator(OU_SEP_STAR)
    call output_line('T R I A _ T E S T 3 D')
    call output_separator(OU_SEP_STAR)
  
    ! Now read in the basic triangulation.
    call tria_readTriFile3D (rtria, './pre/CUBE.tri')
    
    ! Optionally refine the mesh
    if(present(NLMAX)) then
      call tria_quickRefine2LevelOrdering (NLMAX-1,rtria)
    end if

    ! Standardise the mesh
    call tria_initStandardMeshFromRaw (rtria)
    
    ! First of all, print out the basic information
    call output_line('Triangulation Scalars')
    call output_line('----------------------------------')
    call output_line('NVT....................: ' // trim(sys_si(rtria%NVT,8)))
    call output_line('NMT....................: ' // trim(sys_si(rtria%NMT,8)))
    call output_line('NAT....................: ' // trim(sys_si(rtria%NAT,8)))
    call output_line('NEL....................: ' // trim(sys_si(rtria%NEL,8)))
    call output_line('NBCT...................: ' // trim(sys_si(rtria%NBCT,8)))
    call output_line('NblindBCT..............: ' // trim(sys_si(rtria%NblindBCT,8)))
    call output_line('NVBD...................: ' // trim(sys_si(rtria%NVBD,8)))
    call output_line('NMBD...................: ' // trim(sys_si(rtria%NMBD,8)))
    call output_line('NABD...................: ' // trim(sys_si(rtria%NABD,8)))
    call output_line('NNVE...................: ' // trim(sys_si(rtria%NNVE,8)))
    call output_line('NNEE...................: ' // trim(sys_si(rtria%NNEE,8)))
    call output_line('NNAE...................: ' // trim(sys_si(rtria%NNAE,8)))
    call output_line('NNelAtVertex...........: ' // trim(sys_si(rtria%NNelAtVertex,8)))
    call output_line('NNelAtEdge.............: ' // trim(sys_si(rtria%NNelAtEdge,8)))
    call output_line('nverticesPerEdge.......: ' // trim(sys_si(rtria%nverticesPerEdge,8)))
    call output_line('nverticesOnAllEdges....: ' // trim(sys_si(rtria%nVerticesOnAllEdges,8)))
    call output_line('nverticesInEachElement.: ' // trim(sys_si(rtria%nverticesInEachElement,8)))
    call output_line('nverticesInAllElements.: ' // trim(sys_si(rtria%nverticesInAllElements,8)))
    call output_line('nadditionalVertices....: ' // trim(sys_si(rtria%NMT,8)))
    
    ! Print the vertice coordinates
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_double2d(rtria%h_DvertexCoords, p_Ddata2)
    call output_line('p_DvertexCoords')
    do i = lbound(p_Ddata2,2), ubound(p_Ddata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ' // trim(sys_sdEP(p_Ddata2(1,i),20,12)) //&
              ' , ' // trim(sys_sdEP(p_Ddata2(2,i),20,12)) // &
              ' , ' // trim(sys_sdEP(p_Ddata2(3,i),20,12)) // ' /)')
    end do
    
    ! Print vertices-at-Element
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IverticesAtElement,p_Idata2)
    call output_line('p_IverticesAtElement')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do
    
    ! Print edges-at-Element
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IedgesAtElement,p_Idata2)
    call output_line('p_IedgesAtElement')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do
    
    ! Print faces-at-Element
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IfacesAtElement,p_Idata2)
    call output_line('h_IfacesAtElement')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do
    
    ! Print neighbours-at-Element
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IneighboursAtElement,p_Idata2)
    call output_line('h_IneighboursAtElement')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do
    
    ! Print vertices-at-edge
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IverticesAtEdge,p_Idata2)
    call output_line('h_IverticesAtEdge')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do
    
    ! Print faces at edge
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IfacesAtEdgeIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IfacesAtEdgeIdx,p_Iidx)
      call storage_getbase_int(rtria%h_IfacesAtEdge,p_Idata)
      call output_line('h_IfacesAtEdgeIdx / h_IfacesAtEdge')
      do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
        call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                         ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
        do j = p_Iidx(i)+1, p_Iidx(i+1)-1
          call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
        end do
      end do
    else
      call output_line('h_IfacesAtEdgeIdx = ST_NOHANDLE')
    end if

    ! Print elements at edge
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IelementsAtEdgeIdx3d .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IelementsAtEdgeIdx3d,p_Iidx)
      call storage_getbase_int(rtria%h_IelementsAtEdge3d,p_Idata)
      call output_line('h_IelementsAtEdgeIdx3d / h_IelementsAtEdge3d')
      do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
        call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                         ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
        do j = p_Iidx(i)+1, p_Iidx(i+1)-1
          call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
        end do
      end do
    else
      call output_line('h_IelementsAtEdgeIdx3d = ST_NOHANDLE')
    end if
    
    ! Print vertices at face
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IverticesAtFace,p_Idata2)
    call output_line('h_IverticesAtFace')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do
    
    ! Print edges at face
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IedgesAtFace,p_Idata2)
    call output_line('h_IedgesAtFace')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do

    ! Print elements at face
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int2d(rtria%h_IelementsAtFace,p_Idata2)
    call output_line('h_IelementsAtFace')
    do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
      call output_line(trim(sys_si(i,8)) // ': (/ ', bnolinebreak = .true.)
      do j = lbound(p_Idata2,1), ubound(p_Idata2,1)-1
        call output_line(trim(sys_si(p_Idata2(j,i),8)) // ' , ', bnolinebreak = .true.)
      end do
      call output_line(trim(sys_si(p_Idata2(ubound(p_Idata2,1),i),8)) // ' /)')
    end do

    ! Print edges at vertices
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IedgesAtVertexIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IedgesAtVertexIdx,p_Iidx)
      call storage_getbase_int(rtria%h_IedgesAtVertex,p_Idata)
      call output_line('h_IedgesAtVertexIdx / h_IedgesAtVertex')
      do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
        call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                         ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
        do j = p_Iidx(i)+1, p_Iidx(i+1)-1
          call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
        end do
      end do
    else
      call output_line('h_IedgesAtVertexIdx = ST_NOHANDLE')
    end if
    
    ! Print elements at vertex
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int(rtria%h_IelementsAtVertexIdx,p_Iidx)
    call storage_getbase_int(rtria%h_IelementsAtVertex,p_Idata)
    call output_line('h_IelementsAtVertexIdx / h_IelementsAtVertex')
    do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
      call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                       ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
      do j = p_Iidx(i)+1, p_Iidx(i+1)-1
        call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
      end do
    end do

    ! Print nodal property
    call output_separator(OU_SEP_MINUS)
    call storage_getbase_int(rtria%h_InodalProperty,p_Idata)
    call output_line('h_InodalProperty')
    do i = lbound(p_Idata,1), ubound(p_Idata,1)
      call output_line(trim(sys_si(i,8)) // ': ' // trim(sys_si(p_Idata(i),8)))
    end do

    ! Print maco nodal property
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_ImacroNodalProperty .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_ImacroNodalProperty,p_Idata)
      call output_line('h_ImacroNodalProperty')
      do i = lbound(p_Idata,1), ubound(p_Idata,1)
        call output_line(trim(sys_si(i,8)) // ': ' // trim(sys_si(p_Idata(i),8)))
      end do
    else
      call output_line('h_ImacroNodalProperty = ST_NOHANDLE')
    end if

    ! Print element volume
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_DelementVolume .ne. ST_NOHANDLE) then
      call storage_getbase_double(rtria%h_DelementVolume,p_Ddata)
      call output_line('h_DelementVolume')
      do i = lbound(p_Ddata,1), ubound(p_Ddata,1)
        call output_line(trim(sys_si(i,8)) // ': ' // trim(sys_sdEP(p_Ddata(i),20,12)))
      end do
    else
      call output_line('h_DelementVolume = ST_NOHANDLE')
    end if


    ! Print refinement patch
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IrefinementPatchIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IrefinementPatchIdx,p_Iidx)
      call storage_getbase_int(rtria%h_IrefinementPatch,p_Idata)
      call output_line('h_IrefinementPatchIdx / h_IrefinementPatch')
      do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
        call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                         ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
        do j = p_Iidx(i)+1, p_Iidx(i+1)-1
          call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
        end do
      end do
    else
      call output_line('h_IrefinementPatchIdx = ST_NOHANDLE')
    end if
    
    ! Print coarse grid element
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IcoarseGridElement .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IcoarseGridElement,p_Idata)
      call output_line('h_IcoarseGridElement')
      do i = lbound(p_Idata,1), ubound(p_Idata,1)
        call output_line(trim(sys_si(i,8)) // ': ' // trim(sys_si(p_Idata(i),8)))
      end do
    else
      call output_line('h_IcoarseGridElement = ST_NOHANDLE')
    end if
    
    ! Print boundary vertices
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IboundaryCpIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IboundaryCpIdx,p_Iidx)
      call storage_getbase_int(rtria%h_IverticesAtBoundary,p_Idata)
      call output_line('h_IboundaryCpIdx / h_IverticesAtBoundary')
      do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
        call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                         ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
        do j = p_Iidx(i)+1, p_Iidx(i+1)-1
          call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
        end do
      end do
    else
      call output_line('h_IboundaryCpIdx = ST_NOHANDLE')
    end if
    
    ! Print boundary edges
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IboundaryCpEdgesIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IboundaryCpEdgesIdx,p_Iidx)
      call storage_getbase_int(rtria%h_IedgesAtBoundary,p_Idata)
      call output_line('h_IboundaryCpEdgesIdx / h_IedgesAtBoundary')
      do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
        call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                         ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
        do j = p_Iidx(i)+1, p_Iidx(i+1)-1
          call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
        end do
      end do
    else
      call output_line('h_IboundaryCpEdgesIdx = ST_NOHANDLE')
    end if
    
    ! Print boundary faces
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IboundaryCpFacesIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IboundaryCpFacesIdx,p_Iidx)
      call storage_getbase_int(rtria%h_IfacesAtBoundary,p_Idata)
      call output_line('h_IboundaryCpFacesIdx / h_IfacesAtBoundary')
      do i = lbound(p_Iidx,1), ubound(p_Iidx,1)-1
        call output_line(trim(sys_si(i,8)) // ' , '// trim(sys_si(p_Iidx(i),8)) //&
                         ' : ' // trim(sys_si(p_Idata(p_Iidx(i)),8)))
        do j = p_Iidx(i)+1, p_Iidx(i+1)-1
          call output_line('           ' // trim(sys_si(j,8)) // ' : ' // trim(sys_si(p_Idata(j),8)))
        end do
      end do
    else
      call output_line('h_IboundaryCpFacesIdx = ST_NOHANDLE')
    end if

    ! Print boundary elements
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IelementsAtBoundary .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IelementsAtBoundary,p_Idata)
      call output_line('h_IelementsAtBoundary')
      do i = lbound(p_Idata,1), ubound(p_Idata,1)
        call output_line(trim(sys_si(i,8)) // ': ' // trim(sys_si(p_Idata(i),8)))
      end do
    else
      call output_line('h_IelementsAtBoundary = ST_NOHANDLE')
    end if

    ! Print vertex parameter values
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_DvertexParameterValue .ne. ST_NOHANDLE) then
      call storage_getbase_double(rtria%h_DvertexParameterValue,p_Ddata)
      call output_line('h_DvertexParameterValue')
      do i = lbound(p_Ddata,1), ubound(p_Ddata,1)
        call output_line(trim(sys_si(i,8)) // ': ' // trim(sys_sdEP(p_Ddata(i),20,12)))
      end do
    else
      call output_line('h_DvertexParameterValue = ST_NOHANDLE')
    end if

    ! Print edge parameter values
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_DedgeParameterValue .ne. ST_NOHANDLE) then
      call storage_getbase_double(rtria%h_DedgeParameterValue,p_Ddata)
      call output_line('h_DedgeParameterValue')
      do i = lbound(p_Ddata,1), ubound(p_Ddata,1)
        call output_line(trim(sys_si(i,8)) // ': ' // trim(sys_sdEP(p_Ddata(i),20,12)))
      end do
    else
      call output_line('h_DedgeParameterValue = ST_NOHANDLE')
    end if
    
    ! Print boundary vertex positions
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IboundaryVertexPos .ne. ST_NOHANDLE) then
      call storage_getbase_int2d(rtria%h_IboundaryVertexPos, p_Idata2)
      call output_line('h_IboundaryVertexPos')
      do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
        call output_line(trim(sys_si(i,8)) // ': (/ ' // trim(sys_si(p_Idata2(1,i),8)) //&
                ' , ' // trim(sys_si(p_Idata2(2,i),8)) // ' /)')
      end do
    else
      call output_line('h_IboundaryVertexPos = ST_NOHANDLE')
    end if
    
    ! Print boundary edge positions
    call output_separator(OU_SEP_MINUS)
    if (rtria%h_IboundaryEdgePos .ne. ST_NOHANDLE) then
      call storage_getbase_int2d(rtria%h_IboundaryEdgePos, p_Idata2)
      call output_line('h_IboundaryEdgePos')
      do i = lbound(p_Idata2,2), ubound(p_Idata2,2)
        call output_line(trim(sys_si(i,8)) // ': (/ ' // trim(sys_si(p_Idata2(1,i),8)) //&
                ' , ' // trim(sys_si(p_Idata2(2,i),8)) // ' /)')
      end do
    else
      call output_line('h_IboundaryEdgePos = ST_NOHANDLE')
    end if

    ! Release the mesh
    call tria_done(rtria)
  
  end subroutine

end module
