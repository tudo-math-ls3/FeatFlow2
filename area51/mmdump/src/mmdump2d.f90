module mmdump2d

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use linearsystemscalar
  use linearsystemblock
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
  use adjacency
  use blockmatassembly
  use blockmatassemblybase
  use blockmatassemblystdop
  
  implicit none

contains

  ! ***************************************************************************

  subroutine mmdump_2d(rparam)
  type(t_parlist), intent(inout) :: rparam

  ! Definitions of variables.
  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  type(t_spatialDiscretisation) :: rdisc
  type(t_blockDiscretisation) :: rdiscBlock
  type(t_scalarCubatureInfo) :: rcub
  type(t_matrixScalar) :: rmat
  type(t_matrixBlock) :: rmatBlock
  type(t_timer) :: rtimer
  integer :: i, nlevel, casmMatrix, cdumpMesh, cdumpMatrix, cdumpColor, ndigVertex, &
    ndigMatrix, casmType, ncubRefines
  integer :: h_Icpt, h_Iprm
  integer(I32) :: celement, ccubature
  character(LEN=64) :: smesh, selement, scubature, sopName, selName, scuName
  character(LEN=256) :: smeshName, smatName, scptName
  real(DP) :: dasmTime

    ! read in the parameters
    call parlst_getvalue_string(rparam, '', 'SMESH', smesh, 'quad')
    call parlst_getvalue_int(rparam, '', 'NLEVEL', nlevel, 0)
    call parlst_getvalue_string(rparam, '', 'SELEMENT', selement, 'EL_Q1_2D')
    call parlst_getvalue_string(rparam, '', 'SCUBATURE', scubature, 'AUTO_G3')
    call parlst_getvalue_int(rparam, '', 'NCUBREFINES', ncubRefines, 0)
    call parlst_getvalue_int(rparam, '', 'CASMMATRIX', casmMatrix, 0)
    call parlst_getvalue_int(rparam, '', 'CASMTYPE', casmType, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPMESH', cdumpMesh, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPMATRIX', cdumpMatrix, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPCOLOR', cdumpColor, 0)
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
    case (EL_EN30_2D)
      selName = 'q1t'
    case default
      selName = elem_getName(celement)
    end select

    ! build cubature rule name
    if(ncubRefines .gt. 0) then
      select case(ccubature)
      case (CUB_GEN_AUTO_G1, CUB_G1_2D)
        scuName = 'g1-r' // trim(sys_sil(ncubRefines,4))
      case (CUB_GEN_AUTO_G2, CUB_G2_2D)
        scuName = 'g2-r' // trim(sys_sil(ncubRefines,4))
      case (CUB_GEN_AUTO_G3, CUB_G3_2D)
        scuName = 'g3-r' // trim(sys_sil(ncubRefines,4))
      case (CUB_GEN_AUTO_G4, CUB_G4_2D)
        scuName = 'g4-r' // trim(sys_sil(ncubRefines,4))
      case (CUB_GEN_AUTO_G5, CUB_G5_2D)
        scuName = 'g5-r' // trim(sys_sil(ncubRefines,4))
      case default
        scuName = cub_getName(ccubature) // '-r' // trim(sys_sil(ncubRefines,4))
      end select
    else
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
    end if

    ! build dumping filenames
    smeshName = trim(smesh) // '_lvl' // trim(sys_sil(nlevel,4)) // '_mesh'
    smatName = trim(smesh) // '_lvl' //  trim(sys_sil(nlevel,4)) // '_' // &
      trim(sopName) // '_' // trim(selName) // '_' // trim(scuName)
    scptName = trim(smesh) // '_lvl' //  trim(sys_sil(nlevel,4)) // '_cpt'

    ! print what we have parsed onto the output stream
    call output_line('Parsed input from data file:')
    call output_line('----------------------------')
    call output_line('sMesh       = ' // trim(smesh))
    call output_line('nLevel      = ' // trim(sys_sil(nlevel,4)))
    call output_line('sElement    = ' // trim(elem_getName(celement)))
    call output_line('sCubature   = ' // trim(cub_getName(ccubature)))
    call output_line('nCubRefines = ' // trim(sys_sil(ncubRefines,4)))
    call output_line('cAsmMatrix  = ' // trim(sys_sil(casmMatrix,4)))
    call output_line('cAsmType    = ' // trim(sys_sil(casmType,4)))
    call output_line('cDumpMesh   = ' // trim(sys_sil(cdumpMesh,4)))
    call output_line('cDumpMatrix = ' // trim(sys_sil(cdumpMatrix,4)))
    call output_line('cDumpColor  = ' // trim(sys_sil(cdumpColor,4)))
    call output_line('nDigVertex  = ' // trim(sys_sil(ndigVertex,4)))
    call output_line('nDigMatrix  = ' // trim(sys_sil(ndigMatrix,4)))
    call output_separator(OU_SEP_STAR)

    ! Read the domain, read the mesh, refine
    call output_line('Reading coarse mesh...')
    call boundary_read_prm(rbnd, './mesh/' // trim(smesh) // '.prm')
    call tria_readTriFile2D (rtria, './mesh/' // trim(smesh) // '.tri', rbnd)

    ! Ensure that we use a quadrilateral mesh - triangular meshes are not supported.
    if(rtria%NNVE .ne. 4) then
      call output_line('ERROR: only quadrilateral meshes supported')
      call sys_halt()
    end if

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
    if(iand(cdumpMesh,1) .ne. 0) then
      call output_line('Dumping mesh to "' // trim(smeshName) // '.txt" in text format...')
      call dumpMeshText(rtria, './out/' // trim(smeshName) // '.txt', ndigVertex)
    end if
    if(iand(cdumpMesh,2) .ne. 0) then
      call output_line('Dumping mesh to "' // trim(smeshName) // '.bin" in binary format...')
      call dumpMeshBinary(rtria, './out/' // trim(smeshName) // '.bin')
    end if

    ! Assemble color partition table?
    if(cdumpColor .gt. 0) then

      ! calculate cpt
      call output_line('Calculating color partition table...')
      call calcColorPart(rtria, rbnd, h_Icpt, h_Iprm)

      ! Dump cpt if desired
      if(iand(cdumpColor,1) .ne. 0) then
        call output_line('Dumping coloring to "' // trim(scptName) // '.txt" in text format...')
        call dumpColorText(h_Icpt, h_Iprm, './out/' // trim(scptName) // '.txt')
      end if
      if(iand(cdumpColor,2) .ne. 0) then
        call output_line('Dumping coloring to "' // trim(scptName) // '.bin" in binary format...')
        call dumpColorBinary(h_Icpt, h_Iprm, './out/' // trim(scptName) // '.bin')
      end if

      ! free cpt
      call storage_free(h_Iprm)
      call storage_free(h_Icpt)

    end if

    ! Assemble matrix?
    if(casmMatrix .gt. 0) then

      ! Set up the discretisation
      call output_line('Setting up discretisation...')
      call spdiscr_initDiscr_simple (rdisc, celement, rtria, rbnd)

      ! Set up cubature
      call output_line('Setting up cubature rule...')
      call spdiscr_createDefCubStructure(rdisc, rcub, ccubature, ncubRefines)

      ! Assemble matrix structure
      call output_line('Assembling matrix structure...')
      call bilf_createMatrixStructure (rdisc, LSYSSC_MATRIX9, rmat)
      call output_line('NEQ : ' // trim(sys_sil(rmat%NEQ,16)))
      call output_line('NNZE: ' // trim(sys_sil(rmat%NA,16)))

      call lsyssc_allocEmptyMatrix(rmat, LSYSSC_SETM_ZERO)

      ! start timer
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer)

      ! Assemble matrix entries
      select case(casmMatrix)
      case (1)
        ! Assemble Mass matrix
        call output_line('Assembling Mass matrix...')
        select case(casmType)
        case (0)
          call stdop_assembleSimpleMatrix (rmat, DER_FUNC, DER_FUNC, 1.0_DP, .true., rcub)
        case (1)
          call spdiscr_createBlockDiscrInd(rdisc, rdiscBlock)
          call lsysbl_createMatFromScalar(rmat, rmatBlock, rdiscBlock)
          call bma_buildMatrix (rmatBlock, BMA_CALC_STANDARD, bma_fcalc_mass, rcubatureInfo=rcub)
        end select

      case (2)
        ! Assemble Laplace matrix
        call output_line('Assembling Laplace matrix...')
        select case(casmType)
        case (0)
          call stdop_assembleLaplaceMatrix (rmat, .true., 1.0_DP, rcub)
        case (1)
          call spdiscr_createBlockDiscrInd(rdisc, rdiscBlock)
          call lsysbl_createMatFromScalar(rmat, rmatBlock, rdiscBlock)
          call bma_buildMatrix (rmatBlock, BMA_CALC_STANDARD, bma_fcalc_laplace, rcubatureInfo=rcub)
        end select

      end select

      ! stop timer and fetch time
      call stat_stopTimer(rtimer)
      dasmTime = rtimer%delapsedReal
      call output_line('Assembly Time: ' // trim(sys_sdl(dasmTime,12)) // ' seconds')

      ! Dump matrix if desired
      if(iand(cdumpMatrix,1) .ne. 0) then
        call output_line('Dumping matrix to "' // trim(smatName) // '.txt" in text format...')
        call dumpMatrixText(rmat, './out/' // trim(smatName) // '.txt', ndigMatrix, dasmTime)
      end if
      if(iand(cdumpMatrix,2) .ne. 0) then
        call output_line('Dumping matrix to "' // trim(smatName) // '.bin" in binary format...')
        call dumpMatrixBinary(rmat, './out/' // trim(smatName) // '.bin', dasmTime)
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

  subroutine calcColorPart(rtria, rbnd, h_Icpt, h_Iperm)
  type(t_triangulation), intent(in) :: rtria
  type(t_boundary), intent(in) :: rbnd
  integer, intent(out) :: h_Icpt, h_Iperm

  type(t_spatialDiscretisation) :: rdiscQ0, rdiscQ1
  type(t_matrixScalar) :: rmatA, rmatB
  integer :: h_Iptr, h_Iidx

    ! create a Q0/Q1 discretisation
    call spdiscr_initDiscr_simple (rdiscQ0, EL_Q0_2D, rtria, rbnd)
    call spdiscr_initDiscr_simple (rdiscQ1, EL_Q1_2D, rtria, rbnd)

    ! assemble a Q0/Q1 matrix structure
    call bilf_createMatrixStructure (rdiscQ1, LSYSSC_MATRIX9, rmatA, rdiscQ0)

    ! transose it
    call lsyssc_transposeMatrix (rmatA, rmatB, LSYSSC_TR_STRUCTURE)

    ! compose adjacency graphs
    call adj_compose(rmatA%h_Kld, rmatA%h_Kcol, rmatB%h_Kld, rmatB%h_Kcol, h_Iptr, h_Iidx)

    ! calculate color partitioning
    call adj_calcColouring(h_Iptr, h_Iidx, h_Icpt, h_Iperm)

    ! clean up
    call storage_free(h_Iidx)
    call storage_free(h_Iptr)
    call lsyssc_releaseMatrix(rmatB)
    call lsyssc_releaseMatrix(rmatA)
    call spdiscr_releaseDiscr(rdiscQ1)
    call spdiscr_releaseDiscr(rdiscQ0)

  end subroutine

 ! ****************************************************************************

  subroutine dumpMeshText(rtria, sfilename, ndigits)
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
    write(iunit,'(A)') 'TMSH'
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

    ! fetch vertex-at-edge indices
    call storage_getbase_int2d(rtria%h_IverticesAtEdge, p_Iidx)

    ! dump vertex indices
    write(iunit,'(A)') 'VERT_AT_EDGE'
    do i = 1, rtria%NMT
      write(iunit,'(A)') &
        trim(sys_sil(p_Iidx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(2,i)-1,16))
    end do

    ! fetch vertex-at-element indices
    call storage_getbase_int2d(rtria%h_IverticesAtElement, p_Iidx)

    ! dump vertex indices
    write(iunit,'(A)') 'VERT_AT_ELEM'
    do i = 1, rtria%NEL
      write(iunit,'(A)') &
        trim(sys_sil(p_Iidx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(3,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(4,i)-1,16))
    end do

    ! fetch edge-at-element indices
    call storage_getbase_int2d(rtria%h_IedgesAtElement, p_Iidx)

    ! dump vertex indices
    write(iunit,'(A)') 'EDGE_AT_ELEM'
    do i = 1, rtria%NEL
      write(iunit,'(A)') &
        trim(sys_sil(p_Iidx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(3,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Iidx(4,i)-1,16))
    end do

    ! close file
    write(iunit,'(A)') 'END_OF_FILE'
    close(iunit)

  end subroutine

 ! ****************************************************************************

  subroutine dumpMeshBinary(rtria, sfilename)
  type(t_triangulation), intent(in) :: rtria
  character(len=*), intent(in) :: sfilename

  real(DP), dimension(:,:), pointer :: p_Dvtx
  integer, dimension(:,:), pointer :: p_IVertAtEdge, p_IVertAtElem, p_IEdgeAtElem

    ! fetch vertex coords and index arrays
    call storage_getbase_double2d(rtria%h_DvertexCoords, p_Dvtx)
    call storage_getbase_int2d(rtria%h_IverticesAtEdge, p_IVertAtEdge)
    call storage_getbase_int2d(rtria%h_IverticesAtElement, p_IVertAtElem)
    call storage_getbase_int2d(rtria%h_IedgesAtElement, p_IEdgeAtElem)

    call mmaux_dumpmesh(sfilename, rtria%NVT, rtria%NMT, rtria%NEL, &
        p_Dvtx, p_IVertAtEdge, p_IVertAtElem, p_IEdgeAtElem)

  end subroutine

  ! ****************************************************************************

  subroutine dumpMatrixText(rmatrix, sfilename, ndigits, dasmTime)
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
    write(iunit, '(A)') 'TMTX'
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
    write(iunit,'(A)') 'END_OF_FILE'
    close(iunit)

  end subroutine

  ! ****************************************************************************

  subroutine dumpMatrixBinary(rmatrix, sfilename, dasmTime)
  type(t_matrixScalar), intent(in) :: rmatrix
  character(len=*), intent(in) :: sfilename
  real(DP), intent(in) :: dasmTime

  real(DP), dimension(:), pointer :: p_DA
  integer, dimension(:), pointer :: p_IrowPtr, p_IcolIdx
  real(SP) :: fasmTime

    ! convert time to float
    fasmTime = real(dasmTime,SP)

    ! fetch row pointer array
    call lsyssc_getbase_Kld(rmatrix, p_IrowPtr)
    ! fetch column index array
    call lsyssc_getbase_Kcol(rmatrix, p_IcolIdx)
    ! fetch data array
    call lsyssc_getbase_double(rmatrix, p_DA)

    ! dump matrix
    call mmaux_dumpmatrix(sfilename, rmatrix%NEQ, rmatrix%NA, fasmTime, &
        p_IrowPtr, p_IcolIdx, p_DA)

  end subroutine


  ! ****************************************************************************

  subroutine dumpColorText(h_Icpt, h_Iperm, sfilename)
  integer, intent(in) :: h_Icpt, h_Iperm
  character(len=*), intent(in) :: sfilename

  integer :: iunit,i,ncol,nadj
  integer, dimension(:), pointer :: p_Icpt, p_Iperm

    ! open file for writing
    call io_openFileForWriting(sfilename, iunit, SYS_REPLACE, bformatted=.true.)

    ! Fail?
    if(iunit .le. 0) then
      call output_line('Failed to open "' // trim(sfilename) // '" for writing!', &
          OU_CLASS_ERROR, OU_MODE_STD)
      call sys_halt()
    end if

    ! fetch arrays
    call storage_getbase_int(h_Icpt, p_Icpt)
    call storage_getbase_int(h_Iperm, p_Iperm)

    ! fetch lengths
    ncol = ubound(p_Icpt, 1) - 1
    nadj = ubound(p_Iperm, 1)

    ! write header
    write(iunit, '(A)') 'TCPT'
    write(iunit, '(A)') &
      trim(sys_sil(ncol,16)) !// ' ' // &
      !trim(sys_sil(radj,16))

    ! dump colot pointer
    write(iunit,'(A)') 'COLOR'
    do i = 1, ncol+1
      write(iunit,'(A)') trim(sys_sil(p_Icpt(i)-1,16))
    end do

    ! dump element index array
    write(iunit,'(A)') 'ELEM_IDX'
    do i = 1, nadj
      write(iunit,'(A)') trim(sys_sil(p_Iperm(i)-1,16))
    end do

    ! close file
    write(iunit,'(A)') 'END_OF_FILE'
    close(iunit)

  end subroutine

  ! ****************************************************************************

  subroutine dumpColorBinary(h_Icpt, h_Iperm, sfilename)
  integer, intent(in) :: h_Icpt, h_Iperm
  character(len=*), intent(in) :: sfilename

  integer :: ncol, nadj
  integer, dimension(:), pointer :: p_Icpt, p_Iperm

    ! fetch arrays
    call storage_getbase_int(h_Icpt, p_Icpt)
    call storage_getbase_int(h_Iperm, p_Iperm)

    ! fetch lengths
    ncol = ubound(p_Icpt, 1) - 1
    nadj = ubound(p_Iperm, 1)

    ! dump matrix
    call mmaux_dumpcolor(sfilename, ncol, nadj, p_Icpt, p_Iperm)

  end subroutine

end module
