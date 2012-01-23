module mmdump2d

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use linearsystemscalar
  use bilinearformevaluation
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use element
  use derivatives
  use stdoperators
  use paramlist
  use io
  use statistics
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine mmdump_2d(rparam)
  type(t_parlist), intent(inout) :: rparam
  
!<description>
!</description>

!</subroutine>

  ! Definitions of variables.
  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  type(t_spatialDiscretisation) :: rdisc
  type(t_scalarCubatureInfo) :: rcub
  type(t_matrixScalar) :: rmat
  type(t_timer) :: rtimer
  integer :: i, nlevel, casmMatrix, cdumpMesh, cdumpMatrix, ndigVertex, ndigMatrix
  integer(I32) :: celement, ccubature
  character(LEN=64) :: smesh, selement, scubature, sopName, selName, scuName
  character(LEN=256) :: smeshName, smatName
  real(DP) :: dasmTime

    ! read in the parameters
    call parlst_getvalue_string(rparam, '', 'SMESH', smesh, 'quad')
    call parlst_getvalue_int(rparam, '', 'NLEVEL', nlevel, 0)
    call parlst_getvalue_string(rparam, '', 'SELEMENT', selement, 'EL_Q1_2D')
    call parlst_getvalue_string(rparam, '', 'SCUBATURE', scubature, 'AUTO_G3')
    call parlst_getvalue_int(rparam, '', 'CASMMATRIX', casmMatrix, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPMESH', cdumpMesh, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPMATRIX', cdumpMatrix, 0)
    call parlst_getvalue_int(rparam, '', 'NDIGVERTEX', ndigVertex, 8)
    call parlst_getvalue_int(rparam, '', 'NDIGMATRIX', ndigMatrix, 8)

    ! parse element and cubature
    celement = elem_igetID(selement)
    ccubature = cub_igetID(scubature)

    ! build operator name
    select case(casmMatrix)
    case (0)
      sopName = 'none'
    case (1)
      sopName = 'mass'
    case (2)
      sopName = 'laplace'
    case default
      ! invalid matrix type
      call output_line('ERROR: invalid CASMMATRIX parameter')
      call sys_halt()
    end select

    ! build element name
    select case(celement)
    case (EL_Q1_2D)
      selName = 'q1'
    case (EL_Q2_2D)
      selName = 'q2'
    case default
      selName = elem_getName(celement)
    end select

    ! build cubature rule name
    select case(ccubature)
    case (CUB_GEN_AUTO_G1, CUB_G1_2D)
      scuName = 'g1'
    case (CUB_GEN_AUTO_G2, CUB_G2_2D)
      scuName = 'g2'
    case (CUB_GEN_AUTO_G3, CUB_G3_2D)
      scuName = 'g3'
    case (CUB_GEN_AUTO_G4, CUB_G4_2D)
      scuName = 'g4'
    case (CUB_GEN_AUTO_G5, CUB_G5_2D)
      scuName = 'g5'
    case default
      scuName = cub_getName(ccubature)
    end select

    ! build dumping filenames
    smeshName = trim(smesh) // '_lvl' // trim(sys_sil(nlevel,4)) // '_mesh.txt'
    smatName = trim(smesh) // '_lvl' //  trim(sys_sil(nlevel,4)) // '_' // trim(sopName) &
      // '_' // trim(selName) // '_' // trim(scuName) // '.txt'

    ! print what we have parsed onto the output stream
    call output_line('Parsed input from data file:')
    call output_line('----------------------------')
    call output_line('sMesh       = ' // trim(smesh))
    call output_line('nLevel      = ' // trim(sys_sil(nlevel,4)))
    call output_line('sElement    = ' // trim(elem_getName(celement)))
    call output_line('sCubature   = ' // trim(cub_getName(ccubature)))
    call output_line('cAmsMatrix  = ' // trim(sys_sil(casmMatrix,4)))
    call output_line('cDumpMesh   = ' // trim(sys_sil(cdumpMesh,4)))
    call output_line('cDumpMatrix = ' // trim(sys_sil(cdumpMatrix,4)))
    call output_line('nDigVertex  = ' // trim(sys_sil(ndigVertex,4)))
    call output_line('nDigMatrix  = ' // trim(sys_sil(ndigMatrix,4)))
    call output_separator(OU_SEP_STAR)

    ! Read the domain, read the mesh, refine
    call output_line('Reading coarse mesh...')
    call boundary_read_prm(rbnd, './mesh/' // trim(smesh) // '.prm')
    call tria_readTriFile2D (rtria, './mesh/' // trim(smesh) // '.tri', rbnd)

    ! Refine up to desired level
    call output_line('Refining mesh...')
    !                      123456789ABCDEF-123456789ABCDEF-123456789ABCDEF-
    call output_line('LEVEL          #VERTS          #EDGES          #QUADS')
    call output_line('-----------------------------------------------------')
    do i = 1, nlevel
      call output_line(trim(sys_si(i-1,5)) // trim(sys_si(rtria%NVT,16)) &
          // trim(sys_si(rtria%NMT,16)) // trim(sys_si(rtria%NEL,16)))
      call tria_quickRefine2LevelOrdering(1, rtria, rbnd)
    end do
    call tria_initStandardMeshFromRaw (rtria,rbnd)
    call output_line(trim(sys_si(nlevel,5)) // trim(sys_si(rtria%NVT,16)) &
        // trim(sys_si(rtria%NMT,16)) // trim(sys_si(rtria%NEL,16)))
    call output_lbrk()

    ! Dump mesh if desired
    if(cdumpMesh .ne. 0) then
      call output_line('Dumping mesh to "' // trim(smeshName) // '"...')
      call dumpMesh(rtria, './out/' // trim(smeshName), ndigVertex)
    end if
    
    ! Assemble matrix?
    if(casmMatrix .gt. 0) then

      ! Set up the discretisation
      call output_line('Setting up discretisation...')
      call spdiscr_initDiscr_simple (rdisc, EL_Q1, rtria, rbnd)

      ! Set up cubature
      call output_line('Setting up cubature rule...')
      call spdiscr_createDefCubStructure(rdisc, rcub, CUB_GEN_AUTO_G3)

      ! Assemble matrix structure
      call output_line('Assembling matrix structure...')
      call bilf_createMatrixStructure (rdisc, LSYSSC_MATRIX9, rmat)
      call output_line('NEQ : ' // trim(sys_sil(rmat%NEQ,16)))
      call output_line('NNZE: ' // trim(sys_sil(rmat%NA,16)))

      ! start timer
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer)

      ! Assemble matrix entries
      select case(casmMatrix)
      case (1)
        ! Assemble Mass matrix
        call output_line('Assembling Mass matrix...')
        call stdop_assembleSimpleMatrix (rmat, DER_FUNC, DER_FUNC, 1.0_DP, .true., rcub)

      case (2)
        ! Assemble Laplace matrix
        call output_line('Assembling Laplace matrix...')
        call stdop_assembleLaplaceMatrix (rmat, .true., 1.0_DP, rcub)

      end select

      ! stop timer and fetch time
      call stat_stopTimer(rtimer)
      dasmTime = rtimer%delapsedReal
      call output_line('Assembly Time: ' // trim(sys_sdl(dasmTime,12)) // ' seconds')

      ! Dump matrix if desired
      if(cdumpMatrix .ne. 0) then
        call output_line('Dumping matrix to "' // trim(smatName) // '"...')
        call dumpMatrix(rmat, './out/' // trim(smatName), ndigMatrix, dasmTime)
      end if

    end if
    
    ! Clean up the mess
    call output_line('Cleaning up the mess...')
    if(casmMatrix .gt. 0) then
      call lsyssc_releaseMatrix (rmat)
      call spdiscr_releaseCubStructure(rcub)
      call spdiscr_releaseDiscr(rdisc)
    end if
    call tria_done (rtria)
    call boundary_release (rbnd)
    
  end subroutine

 ! ****************************************************************************

  subroutine dumpMesh(rtria, sfilename, ndigits)
  type(t_triangulation), intent(in) :: rtria
  character(len=*), intent(in) :: sfilename
  integer, intent(in) :: ndigits

  integer :: iunit,i
  real(DP), dimension(:,:), pointer :: p_Dvtx
  integer, dimension(:,:), pointer :: p_Iidx

    ! open file for writing
    call io_openFileForWriting(sfilename, iunit, SYS_REPLACE, bformatted=.true.)

    ! Fail?
    if(iunit .le. 0) then
      call output_line('Failed to open "' // trim(sfilename) // '" for writing!', &
          OU_CLASS_ERROR, OU_MODE_STD)
      call sys_halt()
    end if

    ! write header
    write(iunit,'(A)') 'MESH'
    write(iunit,'(A)') &
      trim(sys_sil(rtria%NVT,16)) // ' ' // &
      trim(sys_sil(rtria%NMT,16)) // ' ' // &
      trim(sys_sil(rtria%NEL,16))

    ! fetch vertex coords
    call storage_getbase_double2d(rtria%h_DvertexCoords, p_Dvtx)

    ! dump vertex coords
    write(iunit,'(A)') 'VERTEX_COORDS'
    do i = 1, rtria%NVT
      write(iunit,'(A)') &
        trim(sys_sdel(p_Dvtx(1,i), ndigits)) // ' ' // &
        trim(sys_sdel(p_Dvtx(2,i), ndigits))
    end do

    ! fetch vertex indices
    call storage_getbase_int2d(rtria%h_IverticesAtElement, p_Iidx)

    ! dump vertex indices
    write(iunit,'(A)') 'VERTEX_INDICES'
    do i = 1, rtria%NEL
      write(iunit,'(A)') &
        trim(sys_sil(p_Iidx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(3,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(4,i)-1,16))
    end do

    ! fetch edge indices
    call storage_getbase_int2d(rtria%h_IedgesAtElement, p_Iidx)

    ! dump vertex indices
    write(iunit,'(A)') 'EDGE_INDICES'
    do i = 1, rtria%NEL
      write(iunit,'(A)') &
        trim(sys_sil(p_Iidx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(3,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(4,i)-1,16))
    end do

    ! close file
    write(iunit,'(A)') 'END'
    close(iunit)

  end subroutine

  ! ****************************************************************************

  subroutine dumpMatrix(rmatrix, sfilename, ndigits, dasmTime)
  type(t_matrixScalar), intent(in) :: rmatrix
  character(len=*), intent(in) :: sfilename
  integer, intent(in) :: ndigits
  real(DP), intent(in) :: dasmTime

  integer :: iunit,i
  real(DP), dimension(:), pointer :: p_DA
  integer, dimension(:), pointer :: p_Idx

    ! open file for writing
    call io_openFileForWriting(sfilename, iunit, SYS_REPLACE, bformatted=.true.)

    ! Fail?
    if(iunit .le. 0) then
      call output_line('Failed to open "' // trim(sfilename) // '" for writing!', &
          OU_CLASS_ERROR, OU_MODE_STD)
      call sys_halt()
    end if

    ! write header
    write(iunit, '(A)') 'MATRIX'
    write(iunit, '(A)') &
      trim(sys_sil(rmatrix%NEQ,16)) // ' ' // &
      trim(sys_sil(rmatrix%NA,16)) // ' ' // &
      trim(sys_sdL(dasmTime,12))

    ! fetch row pointer array
    call lsyssc_getbase_Kld(rmatrix, p_Idx)

    ! dump row pointer
    write(iunit,'(A)') 'ROW_PTR'
    do i = 1, rmatrix%NEQ+1
      write(iunit,'(A)') trim(sys_sil(p_Idx(i)-1,16))
    end do

    ! fetch column index array
    call lsyssc_getbase_Kcol(rmatrix, p_Idx)

    ! dump column index array
    write(iunit,'(A)') 'COL_IDX'
    do i = 1, rmatrix%NA
      write(iunit,'(A)') trim(sys_sil(p_Idx(i)-1,16))
    end do

    ! fetch data array
    call lsyssc_getbase_double(rmatrix, p_DA)

    ! dump matrix data array
    write(iunit,'(A)') 'DATA'
    do i = 1, rmatrix%NA
      write(iunit,'(A)') trim(sys_sdEL(p_DA(i),ndigits))
    end do

    ! close file
    write(iunit,'(A)') 'END'
    close(iunit)

  end subroutine

end module
