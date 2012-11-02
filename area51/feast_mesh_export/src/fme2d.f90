module fme2d

  use fsystem
  use genoutput
  use paramlist
  use storage
  use basicgeometry
  use boundary
  use triangulation
  use io
  
  use fmecore

  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fme_2d(rparam)
  type(t_parlist), intent(inout) :: rparam

!</subroutine>

  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  integer :: nLevel
  character(len=SYS_STRLEN) :: sPrmFile, sTriFile, sMeshFile, sChartFile
  
    ! fetch parameters
    call parlst_getvalue_string(rparam, '', 'SPRMFILE', sPrmFile, '')
    call parlst_getvalue_string(rparam, '', 'STRIFILE', sTriFile, '')
    call parlst_getvalue_int(rparam, '', 'NLEVEL', nLevel, 0)

    ! Read in parametrisation
    call boundary_read_prm(rbnd, './mesh/' // trim(sPrmFile) // '.prm')

    ! Read in triangulation
    call tria_readTriFile2D (rtria, './mesh/' // trim(sTriFile) // '.tri', rbnd)

    ! Refine up to desired level
    call tria_quickRefine2LevelOrdering (nLevel, rtria, rbnd)
    call tria_initStandardMeshFromRaw (rtria, rbnd)
    
    ! write chart file
    if(parlst_queryvalue(rparam, '', 'SCHARTFILE') .gt. 0) then
      call parlst_getvalue_string(rparam, '', 'SCHARTFILE', sChartFile, '')
      call fma2d_writeChartFile('./out/' // trim(sChartFile), rbnd, rparam)
    end if

    ! write mesh file
    if(parlst_queryvalue(rparam, '', 'SMESHFILE') .gt. 0) then
      call parlst_getvalue_string(rparam, '', 'SMESHFILE', sMeshFile, '')
      call fma2d_writeMeshFile('./out/' // trim(sMeshFile), rtria, rbnd, rparam)
    end if

    ! Clean up the mess
    call tria_done (rtria)
    call boundary_release (rbnd)

  end subroutine

  ! ***********************************************************************************************

  subroutine fma2d_writeChartFile(sChartPath, rbnd, rparam)
  character(len=*), intent(in) :: sChartPath
  type(t_boundary), intent(inout) :: rbnd
  type(t_parlist), intent(inout) :: rparam

  character(len=256) :: sInfoLine
  integer :: iunit, nlines, i, ibct
  type(t_boundaryRegion) :: rrgn

    ! open file for writing
    call io_openFileForWriting(trim(sChartPath), iunit, SYS_REPLACE) !, bformatted=.true.)

    ! Fail?
    if(iunit .le. 0) then
      call output_line('Failed to open "' // trim(sChartPath) // '" for writing!', &
          OU_CLASS_ERROR, OU_MODE_STD)
      call sys_halt()
    end if
    
    ! write file starter
    write(iunit,'(A)') '<feast_chart_file>'
    
    ! write file header
    write(iunit,'(A)') '<header>'
    write(iunit,'(A)') ' version 1'
    write(iunit,'(A)') ' charts ' // trim(sys_sil(boundary_igetNBoundComp(rbnd),4))
    write(iunit,'(A)') '</header>'
    
    ! write file info (if present)
    nlines = parlst_querysubstrings(rparam, '', 'sChartFileInfo')
    if(nlines .gt. 0) then
      write(iunit,'(A)') '<info>'
      do i=1, nlines
        call parlst_getvalue_string (rparam, '', 'sChartFileInfo', sInfoLine, '', i)
        write(iunit,'(A)') ' ' // trim(sInfoLine)
      end do
      write(iunit,'(A)') '</info>'
    end if

    ! loop over all boundary components
    do ibct = 1, boundary_igetNBoundComp(rbnd)
    
      ! check for uniformity of boundary component
      call fma2d_checkBndComp(rbnd, ibct)
    
      ! write chart chunk
      call boundary_createRegion(rbnd, ibct, 1, rrgn)
      select case(rrgn%isegmentType)
      case (BOUNDARY_TYPE_LINE)
        call fma2d_writeChartPolyline(iunit, rbnd, ibct)
        
      case (BOUNDARY_TYPE_CIRCLE)
        call fma2d_writeChartCircle(iunit, rbnd, ibct)
      end select
      
    end do
    
    ! write file terminator
    write(iunit,'(A)') '</feast_chart_file>'
    
    ! close output file
    close(iunit)  

  end subroutine
  
  ! ***********************************************************************************************

  subroutine fma2d_checkBndComp(rbnd, ibct)
  type(t_boundary), intent(inout) :: rbnd
  integer, intent(in) :: ibct
  
  type(t_boundaryRegion) :: rrgn
  integer :: iseg, nseg, ctype
  
    nseg = boundary_igetNsegments(rbnd, ibct)
    call boundary_createRegion(rbnd, ibct, 1, rrgn)
    ctype = rrgn%isegmentType
    
    do iseg=2, nseg
      call boundary_createRegion(rbnd, ibct, iseg, rrgn)
      if(ctype .ne. rrgn%isegmentType) then
        call output_line('ERROR: boundary component ' // trim(sys_sil(ibct,4)) &
          // ' has mixed segment types!', OU_CLASS_ERROR, OU_MODE_STD)
        call sys_halt()
      end if
    end do

  end subroutine

  ! ***********************************************************************************************

  subroutine fma2d_writeChartPolyline(iunit, rbnd, ibct)
  integer, intent(in) :: iunit
  type(t_boundary), intent(inout) :: rbnd
  integer, intent(in) :: ibct
  
  integer :: nseg, i
  real(DP) :: dx,dy
  
    ! fetch the segment count
    nseg = boundary_igetNsegments(rbnd, ibct)
    
    ! write chart starter
    write(iunit,'(A)') '<chart>'

    ! write header
    write(iunit,'(A)') ' <header>'
    write(iunit,'(A)') '  name ' // trim(sys_sil(ibct,8))
    write(iunit,'(A)') '  type polyline'
    write(iunit,'(A)') '  img_dim 2'
    write(iunit,'(A)') ' </header>'
    
    ! write polyline chunk
    write(iunit,'(A)') ' <polyline>'
    write(iunit,'(A)') '  num_points ' // trim(sys_sil(nseg+1,8))
    write(iunit,'(A)') '  <points>'
    do i=0, nseg
      call boundary_getCoords(rbnd, ibct, real(i,DP), dx, dy)
      write(iunit,'(A)') '   ' // trim(sys_sil(i,8)) // ' ' // &
        trim(fmtCoord(dx)) // ' ' // trim(fmtCoord(dy))
    end do
    write(iunit,'(A)') '  </points>'
    write(iunit,'(A)') ' </polyline>'
    
    ! write chart terminator
    write(iunit,'(A)') '</chart>'

  end subroutine

  ! ***********************************************************************************************

  subroutine fma2d_writeChartCircle(iunit, rbnd, ibct)
  integer, intent(in) :: iunit
  type(t_boundary), intent(inout) :: rbnd
  integer, intent(in) :: ibct
  
  real(DP) :: dz, dx0,dy0, dx1, dy1, dr
  
    ! fetch maximum parameter value
    dz = boundary_dgetMaxParVal(rbnd, ibct)
    
    ! fetch point coords
    call boundary_getCoords(rbnd, ibct, 0.0_DP, dx0, dy0)
    call boundary_getCoords(rbnd, ibct, 0.5_DP*dz, dx1, dy1)
    
    ! compute radius
    dr = 0.5_DP * sqrt((dx1-dx0)**2 + (dy1-dy0)**2)
  
    ! write chart starter
    write(iunit,'(A)') '<chart>'

    ! write header
    write(iunit,'(A)') ' <header>'
    write(iunit,'(A)') '  name ' // trim(sys_sil(ibct,8))
    write(iunit,'(A)') '  type circle'
    write(iunit,'(A)') ' </header>'
    
    ! write circle chunk
    write(iunit,'(A)') ' <circle>'
    write(iunit,'(A)') '  radius ' // trim(fmtCoord(dr))
    write(iunit,'(A)') '  midpoint ' // &
      trim(fmtCoord(0.5_DP*(dx0+dx1))) // ' ' // &
      trim(fmtCoord(0.5_DP*(dy0+dy1)))
    write(iunit,'(A)') '  domain 0 ' // trim(fmtCoord(dz))
    write(iunit,'(A)') ' </circle>'
    
    ! write circle terminator
    write(iunit,'(A)') '</chart>'

  end subroutine
  
  ! ***********************************************************************************************

  subroutine fma2d_writeMeshFile(sMeshPath, rtria, rbnd, rparam)
  character(len=*), intent(in) :: sMeshPath
  type(t_triangulation), intent(inout) :: rtria
  type(t_boundary), intent(inout) :: rbnd
  type(t_parlist), intent(inout) :: rparam

  character(len=256) :: sInfoLine, sChartFile
  integer :: iunit, nlines, i, ibct
  logical :: bhaveChart
  
    bhaveChart = .false.

    ! open file for writing
    call io_openFileForWriting(sMeshPath, iunit, SYS_REPLACE) !, bformatted=.true.)

    ! Fail?
    if(iunit .le. 0) then
      call output_line('Failed to open "' // trim(sMeshPath) // '" for writing!', &
          OU_CLASS_ERROR, OU_MODE_STD)
      call sys_halt()
    end if
    
    ! write file starter
    write(iunit,'(A)') '<feast_mesh_file>'
    
    ! write file header
    write(iunit,'(A)') '<header>'
    write(iunit,'(A)') ' version 1'
    if(parlst_queryvalue(rparam, '', 'SCHARTFILE') .gt. 0) then
      ! write chart file name
      call parlst_getvalue_string(rparam, '', 'SCHARTFILE', sChartFile, '')
      write(iunit,'(A)') ' chart_file ' // trim(sChartFile)
      bhaveChart = .true.
    end if
    write(iunit,'(A)') ' submeshes ' // trim(sys_sil(rtria%NBCT,4))
    write(iunit,'(A)') ' cellsets 0'
    write(iunit,'(A)') '</header>'
    
    ! write file info (if present)
    nlines = parlst_querysubstrings(rparam, '', 'sMeshFileInfo')
    if(nlines .gt. 0) then
      write(iunit,'(A)') '<info>'
      do i=1, nlines
        call parlst_getvalue_string (rparam, '', 'sMeshFileInfo', sInfoLine, '', i)
        write(iunit,'(A)') ' ' // trim(sInfoLine)
      end do
      write(iunit,'(A)') '</info>'
    end if
    
    ! write mesh chunk
    call fma2d_writeMeshChunk(iunit, rtria, rbnd)
    
    ! write submeshes
    do ibct = 1, rtria%NBCT
      call fma2d_writeSubmeshChunk(iunit, rtria, rbnd, ibct, bhaveChart)
    end do
    
    ! write file terminator
    write(iunit,'(A)') '</feast_mesh_file>'
    
    ! close output file
    close(iunit)  
  
  end subroutine

  ! ***********************************************************************************************

  subroutine fma2d_writeMeshChunk(iunit, rtria, rbnd)
  integer, intent(in) :: iunit
  type(t_triangulation), intent(inout) :: rtria
  type(t_boundary), intent(inout) :: rbnd
  
    ! write mesh starter
    write(iunit,'(A)') '<mesh>'
    
    ! write mesh header
    write(iunit,'(A)') ' <header>'
    write(iunit,'(A)') '  type conformal'
    if(rtria%InelOfType(TRIA_NVETRI2D) .eq. 0) then
      write(iunit,'(A)') '  shape quad'
    else if(rtria%InelOfType(TRIA_NVEQUAD2D) .eq. 0) then
      write(iunit,'(A)') '  shape tria'
    else
      write(iunit,'(A)') '  shape tria-quad'
    end if
    write(iunit,'(A)') '  coords 2'
    write(iunit,'(A)') ' </header>'

    ! write counts
    write(iunit,'(A)') ' <counts>'
    write(iunit,'(A)') '  verts ' // trim(sys_sil(rtria%NVT, 16))
    write(iunit,'(A)') '  edges ' // trim(sys_sil(rtria%NMT, 16))
    if(rtria%InelOfType(TRIA_NVETRI2D) .ne. 0) &
      write(iunit,'(A)') '  trias ' // trim(sys_sil(rtria%InelOfType(TRIA_NVETRI2D), 16))
    if(rtria%InelOfType(TRIA_NVEQUAD2D) .ne. 0) &
      write(iunit,'(A)') '  quads ' // trim(sys_sil(rtria%InelOfType(TRIA_NVEQUAD2D), 16))
    write(iunit,'(A)') ' </counts>'
    
    ! write coords
    call fme_writeCoordsChunk(iunit, rtria)
    
    ! write adjacencies
    call fma_writeAdjacencyChunkEdge(iunit, rtria%h_IverticesAtEdge, rtria%NMT)
    if(rtria%InelOfType(TRIA_NVETRI2D) .eq. 0) then
      call fma_writeAdjacencyChunkQuad(iunit, rtria%h_IverticesAtElement, rtria%NEL)
    else if(rtria%InelOfType(TRIA_NVEQUAD2D) .eq. 0) then
      call fma_writeAdjacencyChunkTria(iunit, rtria%h_IverticesAtElement, rtria%NEL)
    else
      call fma_writeAdjacencyChunkTriaQuad(iunit, rtria%h_IverticesAtElement, rtria%NEL)
    end if
    
    ! write mesh terminator
    write(iunit,'(A)') '</mesh>'
    
  end subroutine
  
  ! ***********************************************************************************************

  subroutine fma2d_writeSubmeshChunk(iunit, rtria, rbnd, ibct, bhaveChart)
  integer, intent(in) :: iunit
  type(t_triangulation), intent(inout) :: rtria
  type(t_boundary), intent(inout) :: rbnd
  integer, intent(in) :: ibct
  logical, intent(in) :: bhaveChart
  
  integer :: i, n0, n1
  integer, dimension(:), pointer :: p_IbdCptIdx, p_Ivtx, p_Iedge
  real(DP), dimension(:), pointer :: p_Dvtx
  
    ! fetch vertex parameter array
    call storage_getbase_double(rtria%h_DvertexParameterValue, p_Dvtx)
    call storage_getbase_int(rtria%h_IboundaryCpIdx, p_IbdCptIdx)
    call storage_getbase_int(rtria%h_IverticesAtBoundary, p_Ivtx)
    call storage_getbase_int(rtria%h_IedgesAtBoundary, p_Iedge)
    
    ! fetch offsets
    n0 = p_IbdCptIdx(ibct)
    n1 = p_IbdCptIdx(ibct+1)
  
    ! write submesh starter
    write(iunit,'(A)') '<submesh>'
    
    ! write submesh header
    write(iunit,'(A)') ' <header>'
    write(iunit,'(A)') '  name ' // trim(sys_sil(ibct,8))
    write(iunit,'(A)') '  parent root'
    if(bhaveChart) &
      write(iunit,'(A)') '  chart ' // trim(sys_sil(ibct,8))
    write(iunit,'(A)') '  type conformal'
    write(iunit,'(A)') '  shape edge'
    write(iunit,'(A)') '  coords 1'
    write(iunit,'(A)') ' </header>'

    ! write counts
    write(iunit,'(A)') ' <counts>'
    write(iunit,'(A)') '  verts ' // trim(sys_sil(n1-n0+1, 16))
    write(iunit,'(A)') '  edges ' // trim(sys_sil(n1-n0, 16))
    write(iunit,'(A)') ' </counts>'
    
    ! write coords
    write(iunit, '(A)') ' <coords>'
    do i=n0, n1-1
      write(iunit, '(A)') '  ' // trim(fmtCoord(p_Dvtx(i)))
    end do
    ! add a last vertex containing the maximum parameter value
    write(iunit, '(A)') '  ' // trim(fmtCoord(boundary_dgetMaxParVal(rbnd, ibct)))
    write(iunit, '(A)') ' </coords>'
    
    ! write vert@edge
    write(iunit, '(A)') ' <vert@edge>'
    do i=1, n1-n0
      write(iunit, '(A)') '  ' // trim(sys_sil(i-1, 16)) // ' ' // trim(sys_sil(i, 16))
    end do
    write(iunit, '(A)') ' </vert@edge>'
    
    ! write vertex indices
    write(iunit, '(A)') ' <vert_idx>'
    do i=n0, n1-1
      write(iunit, '(A)') '  ' // trim(sys_sil(p_Ivtx(i)-1, 16))
    end do
    write(iunit, '(A)') '  ' // trim(sys_sil(p_Ivtx(n0)-1, 16))
    write(iunit, '(A)') ' </vert_idx>'
    
    ! write edge indices
    write(iunit, '(A)') ' <edge_idx>'
    do i=n0, n1-1
      write(iunit, '(A)') '  ' // trim(sys_sil(p_Iedge(i)-1, 16))
    end do
    write(iunit, '(A)') ' </edge_idx>'
    
    ! write submesh terminator
    write(iunit,'(A)') '</submesh>'

  end subroutine
  
end module
