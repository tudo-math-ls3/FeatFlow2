!##############################################################################
!# ****************************************************************************
!# <name> gridtransdbg_test1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module gridtransdbg_test1

  use fsystem
  use genoutput
  use paramlist
  use basicgeometry
  use storage
  use cubature
  use basicgeometry
  use boundary
  use triangulation
  use element
  use spatialdiscretisation
  use multileveloperators

  use gridtransdbg_aux
  
  implicit none

contains

! ***************************************************************************

!<subroutine>

  subroutine gridtransdbg_1(rparam, sConfigSection, itest)
  type(t_parlist), intent(INOUT) :: rparam
  character(LEN=*), intent(IN) :: sConfigSection
  integer, intent(IN) :: itest
  
!<description>
!</description>

!</subroutine>
  
  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtriaC, rtriaF
  type(t_blockDiscretisation) :: rdiscrC, rdiscrF
  integer :: ilevel, itransferOperator, iprojTypeA, iprojTypeB, cavrgType
  integer(I32) :: celement, ccubature, cshape
  character(LEN=64) :: selement,scubature
  real(DP), dimension(:,:), pointer :: p_DA, p_DB, p_DE
  real(DP) :: dtol
  integer :: i,j,m,n
  
    dtol = SYS_EPSREAL_DP * 100.0_DP
      
    ! Fetch minimum and maximum levels
    call parlst_getvalue_int(rparam, sConfigSection, 'ILEVEL', ILEVEL, 1)
    
    ! Fetch element and cubature rule
    call parlst_getvalue_string(rparam, sConfigSection, 'SELEMENT', selement, '')
    call parlst_getvalue_string(rparam, sConfigSection, 'SCUBATURE', scubature, '')

    ! Fetch transfer operator
    call parlst_getvalue_int(rparam, sConfigSection, 'ITRANSFEROPERATOR', itransferOperator, -1)

    ! Fetch projection types
    call parlst_getvalue_int(rparam, sConfigSection, 'IPROJTYPE_A', iprojTypeA, -1)
    call parlst_getvalue_int(rparam, sConfigSection, 'IPROJTYPE_B', iprojTypeB, -1)

    ! Parse element and cubature
    celement = elem_igetID(selement)
    ccubature = cub_igetID(scubature)
    
    ! Get the shape of the element
    cshape = elem_igetShape(celement)
    if(cshape .eq. BGEOM_SHAPE_UNKNOWN) then
      call output_line('Element is not a valid element identifier!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()
    end if
    if(cshape .ne. cub_igetShape(ccubature)) then
      call output_line('Element and cubature formula incompatible!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()
    end if
    
    cavrgType = MLOP_AVRG_ARITHMETIC
    
    call output_separator(OU_SEP_STAR)
    call output_line('GRID-TRANSFER-DEBUGGER: TEST #1')
    call output_line('===============================')
    
    ! Print out that we are going to do:
    call output_line('ILEVEL.............: ' // trim(sys_siL(ILEVEL,4)))
    select case(itransferOperator)
    case(1)
      call output_line('Transfer operator..: Prolongation')
    case(2)
      call output_line('Transfer operator..: Restriction')
    case(3)
      call output_line('Transfer operator..: Interpolation')
    case default
      call output_line('Invalid ITRANSFEROPERATOR parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()
    end select
    select case(iprojTypeA)
    case(1)
      call output_line('Projection Type A..: hard-coded')
    case(2)
      call output_line('Projection Type A..: L2-projection')
    case(3)
      call output_line('Projection Type A..: sliced L2-projection')
    case default
      call output_line('Invalid IPROJTYPE_A parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()
    end select
    select case(iprojTypeB)
    case(1)
      call output_line('Projection Type B..: hard-coded')
    case(2)
      call output_line('Projection Type B..: L2-projection')
    case(3)
      call output_line('Projection Type B..: sliced L2-projection')
    case default
      call output_line('Invalid IPROJTYPE_B parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()
    end select
    call output_lbrk()
    
    ! Create the coarse-mesh triangulation
    select case(cshape)
    case (BGEOM_SHAPE_LINE)
      ! Create 1D line mesh
      call tria_createRawTria1D(rtriaC, 0.0_DP, 1.0_DP, 1)
      
    case (BGEOM_SHAPE_TRIA)
      ! Create 2D triangular mesh
      call boundary_read_prm(rbnd, './pre/TRIA.prm')
      call tria_readTriFile2D(rtriaC, './pre/TRIA.tri', rbnd)
    
    case (BGEOM_SHAPE_QUAD)
      ! Create 2D quadrilateral mesh
      call boundary_read_prm(rbnd, './pre/QUAD.prm')
      call tria_readTriFile2D(rtriaC, './pre/QUAD.tri', rbnd)
    
    case (BGEOM_SHAPE_TETRA)
      call output_line('Tetrahedral elements not supported!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()
    
    case (BGEOM_SHAPE_HEXA)
      ! Create 3D hexahedral mesh
      call tria_readTriFile3D(rtriaC, './pre/CUBE.tri', rbnd)

    case (BGEOM_SHAPE_PYRA)
      call output_line('Pyramid elements not supported!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()

    case (BGEOM_SHAPE_PRISM)
      call output_line('Prism elements not supported!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg_1')
      call sys_halt()
    end select
    
    ! Refine up to ilevel
    call tria_quickRefine2LevelOrdering (ilevel-1, rtriaC, rbnd)
    
    ! Create standard mesh
    call tria_initStandardMeshFromRaw (rtriaC)
    
    ! Refine once more
    call tria_refine2LevelOrdering(rtriaC, rtriaF, rbnd)
    
    ! Create standard mesh
    call tria_initStandardMeshFromRaw (rtriaF)
    
    ! Set up coarse and fine mesh discretisations
    call spdiscr_initBlockDiscr (rdiscrC, 1, rtriaC)
    call spdiscr_initDiscr_simple (rdiscrC%RspatialDiscr(1), &
        celement, ccubature, rtriaC)
    call spdiscr_initBlockDiscr (rdiscrF, 1, rtriaF)
    call spdiscr_initDiscr_simple (rdiscrF%RspatialDiscr(1), &
        celement, ccubature, rtriaF)
    
    ! Assemble matrices
    p_DA => null()
    p_DB => null()
    select case(itransferOperator)
    case(1)
      ! Prolongation
      select case(iprojTypeA)
      case(1)
        call gtaux_asmProlMat_hc(rdiscrC, rdiscrF, p_DA)
      case(2)
        call gtaux_asmProlMat_l2(rdiscrC, rdiscrF, p_DA)
      case(3)
        call gtaux_asmProlMat_sl2(rdiscrC, rdiscrF, p_DA, cavrgType)
      end select
      select case(iprojTypeB)
      case(1)
        call gtaux_asmProlMat_hc(rdiscrC, rdiscrF, p_DB)
      case(2)
        call gtaux_asmProlMat_l2(rdiscrC, rdiscrF, p_DB)
      case(3)
        call gtaux_asmProlMat_sl2(rdiscrC, rdiscrF, p_DB, cavrgType)
      end select
    
    case(2)
      ! Restriction
      select case(iprojTypeA)
      case(1)
        call gtaux_asmRestMat_hc(rdiscrC, rdiscrF, p_DA)
      case(2)
        call gtaux_asmRestMat_l2(rdiscrC, rdiscrF, p_DA)
      case(3)
        call gtaux_asmRestMat_sl2(rdiscrC, rdiscrF, p_DA, cavrgType)
      end select
      select case(iprojTypeB)
      case(1)
        call gtaux_asmRestMat_hc(rdiscrC, rdiscrF, p_DB)
      case(2)
        call gtaux_asmRestMat_l2(rdiscrC, rdiscrF, p_DB)
      case(3)
        call gtaux_asmRestMat_sl2(rdiscrC, rdiscrF, p_DB, cavrgType)
      end select
    
    case(3)
      ! Interpolation
      select case(iprojTypeA)
      case(1)
        call gtaux_asmInterpMat_hc(rdiscrC, rdiscrF, p_DA)
      case(2)
        call gtaux_asmInterpMat_l2(rdiscrC, rdiscrF, p_DA)
      case(3)
        call gtaux_asmInterpMat_sl2(rdiscrC, rdiscrF, p_DA, cavrgType)
      end select
      select case(iprojTypeB)
      case(1)
        call gtaux_asmInterpMat_hc(rdiscrC, rdiscrF, p_DB)
      case(2)
        call gtaux_asmInterpMat_l2(rdiscrC, rdiscrF, p_DB)
      case(3)
        call gtaux_asmInterpMat_sl2(rdiscrC, rdiscrF, p_DB, cavrgType)
      end select
    end select
    
    ! Allocate error matrix
    m = ubound(p_DA,1)
    n = ubound(p_DA,2)
    allocate(p_DE(m,n))
    
    ! Calculate error matrix
    do j = 1, n
      do i = 1, m
      
        p_DE(i,j) = p_DA(i,j) - p_DB(i,j)

        if(abs(p_DA(i,j)) .le. dtol) &
          p_DA(i,j) = 0.0_DP
        if(abs(p_DB(i,j)) .le. dtol) &
          p_DB(i,j) = 0.0_DP
        if(abs(p_DE(i,j)) .le. dtol) &
          p_DE(i,j) = 0.0_DP

      end do ! i
    end do ! j
    
    ! Print results to screen      -                      -                      -
    call output_separator(OU_SEP_MINUS)
    call output_line('Entry        Projection A          Projection B          Error')
    do i = 1, m
      do j = 1, n
        call output_line('(' // trim(sys_si(i,3)) // ',' // trim(sys_si(j,3)) &
          // ') = ' // trim(sys_sdEP(p_DA(i,j),20,13)) &
          // '  '   // trim(sys_sdEP(p_DB(i,j),20,13)) &
          // '  '   // trim(sys_sdEP(p_DE(i,j),20,13)))
      end do ! j
    end do ! i
    
    ! Deallocate matrices
    deallocate(p_DE)
    deallocate(p_DB)
    deallocate(p_DA)
        
    ! Release discretisations
    call spdiscr_releaseBlockDiscr(rdiscrF)
    call spdiscr_releaseBlockDiscr(rdiscrC)
    
    ! Release triangulations
    call tria_done (rtriaF)
    call tria_done (rtriaC)
    
    ! Release boundary
    call boundary_release (rbnd)

  end subroutine

end module