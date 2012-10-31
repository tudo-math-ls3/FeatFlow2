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
  type(t_parlist), intent(INOUT) :: rparam

!</subroutine>

  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  integer :: nLevel, iunit, ibct, nlines, i
  character(len=SYS_STRLEN) :: sPrmFile, sTriFile, sMeshFile, sMeshPath
  character(len=256) :: sInfoLine
  
    ! fetch parameters
    call parlst_getvalue_string(rparam, '', 'SPRMFILE', sPrmFile, '')
    call parlst_getvalue_string(rparam, '', 'STRIFILE', sTriFile, '')
    call parlst_getvalue_string(rparam, '', 'SMESHFILE', sMeshFile, '')
    call parlst_getvalue_int(rparam, '', 'NLEVEL', nLevel, 0)

    ! Read in parametrisation
    call boundary_read_prm(rbnd, './mesh/' // trim(sPrmFile) // '.prm')

    ! Read in triangulation
    call tria_readTriFile2D (rtria, './mesh/' // trim(sTriFile) // '.tri', rbnd)

    ! Refine up to desired level
    call tria_quickRefine2LevelOrdering (nLevel, rtria, rbnd)
    call tria_initStandardMeshFromRaw (rtria, rbnd)

    ! open file for writing
    sMeshPath = './out/' // trim(sMeshFile) // '.txt'
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
    write(iunit,'(A)') ' submeshes ' // trim(sys_sil(rtria%NBCT,4))
    write(iunit,'(A)') ' cellsets 0'
    write(iunit,'(A)') '</header>'
    
    ! write file info (if present)
    nlines = parlst_querysubstrings(rparam, '', 'sFileInfo')
    if(nlines .gt. 0) then
      write(iunit,'(A)') '<info>'
      do i=1, nlines
        call parlst_getvalue_string (rparam, '', 'sFileInfo', sInfoLine, '', i)
        write(iunit,'(A)') ' ' // trim(sInfoLine)
      end do
      write(iunit,'(A)') '</info>'
    end if
    
    
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
    call fme_writeCoordsChunk(iunit, rtria, 8)
    
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
    
    ! write submeshes
    do ibct = 1, rtria%NBCT
      call fma2d_writeSubmeshChunk(iunit, rtria, rbnd, ibct, 8)
    end do
    
    ! write file terminator
    write(iunit,'(A)') '</feast_mesh_file>'
    
    ! close output file
    close(iunit)

    ! Clean up the mess
    call tria_done (rtria)
    call boundary_release (rbnd)

  end subroutine

  subroutine fma2d_writeSubmeshChunk(iunit, rtria, rbnd, ibct, ndigits)
  integer, intent(in) :: iunit
  type(t_triangulation), intent(inout) :: rtria
  type(t_boundary), intent(inout) :: rbnd
  integer, intent(in) :: ibct
  integer, intent(in) :: ndigits
  
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
      write(iunit, '(A)') '  ' // trim(sys_sdel(p_Dvtx(i), ndigits))
    end do
    ! add a last vertex containing the maximum parameter value
    write(iunit, '(A)') '  ' // trim(sys_sdel(boundary_dgetMaxParVal(rbnd, ibct), ndigits))
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
