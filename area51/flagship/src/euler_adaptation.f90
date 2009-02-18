!##############################################################################
!# ****************************************************************************
!# <name> euler_adaptation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# grid adaptation for the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!# 1.) fs_performGridAdaptation
!#     -> perform adaptive grid refinement/coarsening
!#
!# 2.) fs_reinitConstOperators1D
!#     -> reinitialize the constant operators after grid adaptation in 1D
!#
!# 3.) fs_reinitConstOperators2D
!#     -> reinitialize the constant operators after grid adaptation in 2D
!#
!# 4.) fs_reinitConstOperators3D
!#     -> reinitialize the constant operators after grid adaptation in 3D
!#
!# The following auxiliary routines are available:
!#
!# 1.) fcb_hadaptCallback1D
!#      -> Callback function which is invoked by the grid adaptation
!#         procedure in 1D
!#
!# 2.) fcb_hadaptCallback2D
!#      -> Callback function which is invoked by the grid adaptation
!#         procedure in 2D
!#
!# 3.) fcb_hadaptCallback3D
!#      -> Callback function which is invoked by the grid adaptation
!#         procedure in 3D
!#
!# </purpose>
!##############################################################################

module euler_adaptation

  use afcstabilisation
  use collection
  use fsystem
  use genoutput
  use graph
  use hadaptaux1d
  use hadaptaux2d
  use hadaptaux3d
  use hadaptivity
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use triangulation

  use euler_basic
  use problem

  implicit none

  private

  public :: fs_performGridAdaptation
  public :: fs_reinitConstOperators1D
  public :: fs_reinitConstOperators2D
  public :: fs_reinitConstOperators3D
  public :: fcb_hadaptCallback1D
  public :: fcb_hadaptCallback2D
  public :: fcb_hadaptCallback3D

contains
  
  !*****************************************************************************

!<subroutine>

  subroutine fs_performGridAdaptation(rcollection, rhadapt, rproblemLevel,&
                                      rsolution, rgridIndicator, fcb_hadaptCallback)

!<description>
    ! This subroutine performs adaptive grid refinement/coarsening
!</description>

!<input>
    ! grid indicator
    type(t_vectorScalar), intent(IN) :: rgridIndicator

    ! callback routines
    include 'intf_hadaptcallback.inc'
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection

    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! Multigrid level for which grid adaptation should be performed
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! Solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution
!</inputoutput>
!</subroutine>

    ! local variables
    integer(I32), dimension(1) :: Ivalue


    ! Start time measurement for grid adaptivity
    call stat_startTimer(rtimer_adaptivity, STAT_TIMERSHORT)

    ! Check if adaptivity structure has been prepared.
    if (rhadapt%iSpec .eq. HADAPT_HAS_PARAMETERS) then
      
      ! Initialize adaptivity structure
      call hadapt_initFromTriangulation(rhadapt, rproblemLevel%rtriangulation)

      ! Add template matrix and solution vector to collection
      call collct_setvalue_matsca(rcollection, "MATRIX_TEMPLATE",&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE), .true.)
      call collct_setvalue_vec(rcollection, "VECTOR_SOLUTION",&
                               rsolution, .true.)

      ! Initialze callback routine
      Ivalue = 0
      call fcb_hadaptCallback(rcollection, HADAPT_OPR_INITCALLBACK,&
                              Ivalue, Ivalue)

    else

      ! Refresh adaptivity structure
      call hadapt_refreshAdaptation(rhadapt, rproblemLevel%rtriangulation)

    end if
        
    ! Perform grid adaptivity
    call hadapt_performAdaptation(rhadapt, rgridIndicator,&
                                  rcollection, fcb_hadaptCallback)
    
    ! Stop time measurement for grid adaptivity
    call stat_stopTimer(rtimer_adaptivity)
    
  end subroutine fs_performGridAdaptation
  
  !*****************************************************************************

!<subroutine>

  subroutine fs_reinitConstOperators1D(rcollection, rhadapt, rproblemLevel,&
                                       rsolution, ieltype, fcb_hadaptCallback)

!<description>
    ! This subroutine reinitializes the constant operators in 1D
!</description>

!<input>
    ! type of finite elements
    integer, intent(IN) :: ieltype
    
    ! callback routines
    include 'intf_hadaptcallback.inc'
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection

    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! Multigrid level for which grid adaptation should be performed
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! Solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution
!</inputoutput>
!</subroutine>

    ! local variables
    integer(I32), dimension(1) :: Ivalue
    integer :: neq,nedge,ivar,jvar

    ! Start time measurement for triangulation
    call stat_startTimer(rtimer_triangulation, STAT_TIMERSHORT)

    ! Generate raw mesh from adaptivity structure
    call hadapt_generateRawMesh(rhadapt, rproblemLevel%rtriangulation)

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw (rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
    
    ! Re-initialize the discretization structure   
    select case(ieltype)
    case (-1,1,11)
      ! P1 finite elements
      call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                    EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
                                    rproblemLevel%rtriangulation,&
                                    rproblemLevel%p_rproblem%rboundary)
          
    case DEFAULT
      call output_line('Unsupproted element type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators1D')
      call sys_halt()
    end select

    ! Stop time measurement for triangulation
    call stat_stopTimer(rtimer_triangulation)
    

    ! Start time measurement for coefficient matrices
    call stat_startTimer(rtimer_assembly_coeff, STAT_TIMERSHORT)

    ! Get sparsity pattern of the template matrix from collection
    Ivalue = 0
    call fcb_hadaptCallback(rcollection, -3, Ivalue, Ivalue)

    ! Resize scalar matrices MC, CX
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)

    ! Generate the global coefficient matrices for the FE method
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                    DER_FUNC1D, DER_FUNC1D)
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                    DER_DERIV1D_X, DER_FUNC1D)

    ! Perform mass lumping ML:=lumped(MC)
    call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_ML))
    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_ML), LSYSSC_LUMP_DIAG)

    ! Resize AFC stabilisation for inviscid part
    neq   = rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ
    nedge = int(0.5*(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NA-&
                     rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ),I32)
    
    call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID), neq, nedge)


    ! Adjust system matrix
    select case(isystemFormat)
      
    case (SYSTEM_INTERLEAVEFORMAT)

      ! Relase pseudo block matrix
      call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))
      
      ! Resize scalar matrix
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ,&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NCOLS,&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NA,&
                               .false., .false., bforce=.true.)

      ! Create pseudo block matrix from global operator
      call lsysbl_createMatFromScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                      rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                      rproblemLevel%rdiscretisation)

      ! Resize solution vector
      call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                       rsolution, .false., .true.)

    case (SYSTEM_BLOCKFORMAT)
      
      ! What kind of global operator should be adopted?
      select case(isystemCoupling)
      case (FLOW_SEGREGATED)
          
        ! Resize diagonal blocks
        do ivar = 1, NVAR1D
          call lsyssc_resizeMatrix(&
              rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,ivar),&
              rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE), .false., .false., .true.)
        end do

        ! Update global matrix structure
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))
        
        ! Resize solution vector
        call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                         rsolution, .false., .true.)

      case (FLOW_ALLCOUPLED)
        
        ! Resize all blocks
        do ivar = 1, NVAR1D
          do jvar = 1, NVAR1D
            call lsyssc_resizeMatrix(&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,jvar),&
                rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE), .false., .false., .true.)
          end do
        end do
        
        ! Update global matrix structure
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))

        ! Resize solution vector
        call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                         rsolution, .false., .true.)

      case DEFAULT
        call output_line('Unsupported block matrix format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators1D')
        call sys_halt()
      end select
      

    case DEFAULT
      call output_line('Unsupproted system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators1D')
      call sys_halt()
    end select

    ! Stop time measurement for coefficient matrices
    call stat_stopTimer(rtimer_assembly_coeff)

  end subroutine fs_reinitConstOperators1D

  !*****************************************************************************

!<subroutine>

  subroutine fs_reinitConstOperators2D(rcollection, rhadapt, rproblemLevel,&      
                                       rsolution, ieltype, fcb_hadaptCallback)

!<description>
    ! This subroutine reinitializes the constant operators in 2D
!</description>

!<input>
    ! type of finite elements
    integer, intent(IN) :: ieltype
    
    ! callback routines
    include 'intf_hadaptcallback.inc'
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection

    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! Multigrid level for which grid adaptation should be performed
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! Solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution
!</inputoutput>
!</subroutine>

    ! local variables
    integer(I32), dimension(1) :: Ivalue
    integer :: neq,nedge,ivar,jvar

    ! Start time measurement for triangulation
    call stat_startTimer(rtimer_triangulation, STAT_TIMERSHORT)

    ! Generate raw mesh from adaptivity structure
    call hadapt_generateRawMesh(rhadapt, rproblemLevel%rtriangulation)

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw (rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
    
    ! Re-initialize the discretization structure   
    select case(ieltype)
    case (1)
      ! P1 finite elements
      call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                    EL_E001, SPDISC_CUB_AUTOMATIC,&
                                    rproblemLevel%rtriangulation,&
                                    rproblemLevel%p_rproblem%rboundary)
      
    case (11)
      ! Q1 finite elements
      call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                    EL_E011, SPDISC_CUB_AUTOMATIC,&
                                    rproblemLevel%rtriangulation,&
                                    rproblemLevel%p_rproblem%rboundary)
      
    case (-1)
      ! mixed P1/Q1 finite elements
      call spdiscr_initDiscr_triquad(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                     EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC,&
                                     SPDISC_CUB_AUTOMATIC,&
                                     rproblemLevel%rtriangulation,&
                                     rproblemLevel%p_rproblem%rboundary)
      
    case DEFAULT
      call output_line('Unsupproted element type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators2D')
      call sys_halt()
    end select

    ! Stop time measurement for triangulation
    call stat_stopTimer(rtimer_triangulation)
    

    ! Start time measurement for coefficient matrices
    call stat_startTimer(rtimer_assembly_coeff, STAT_TIMERSHORT)

    ! Get sparsity pattern of the template matrix from collection
    Ivalue = 0
    call fcb_hadaptCallback(rcollection, -3, Ivalue, Ivalue)

    ! Resize scalar matrices MC, CX, and CY
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)

    ! Generate the global coefficient matrices for the FE method
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                    DER_FUNC2D, DER_FUNC2D)
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                    DER_DERIV2D_X, DER_FUNC2D)
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                                    DER_DERIV2D_Y, DER_FUNC2D)

    ! Perform mass lumping ML:=lumped(MC)
    call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_ML))
    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_ML), LSYSSC_LUMP_DIAG)

    ! Resize AFC stabilisation for inviscid part
    neq   = rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ
    nedge = int(0.5*(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NA-&
                     rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ),I32)
    
    call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID), neq, nedge)


    ! Adjust system matrix
    select case(isystemFormat)
      
    case (SYSTEM_INTERLEAVEFORMAT)

      ! Relase pseudo block matrix
      call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))
      
      ! Resize scalar matrix
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ,&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NCOLS,&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NA,&
                               .false., .false., bforce=.true.)

      ! Create pseudo block matrix from global operator
      call lsysbl_createMatFromScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                      rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                      rproblemLevel%rdiscretisation)

      ! Resize solution vector
      call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                       rsolution, .false., .true.)

    case (SYSTEM_BLOCKFORMAT)
      
      ! What kind of global operator should be adopted?
      select case(isystemCoupling)
      case (FLOW_SEGREGATED)
          
        ! Resize diagonal blocks
        do ivar = 1, NVAR2D
          call lsyssc_resizeMatrix(&
              rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,ivar),&
              rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE), .false., .false., .true.)
        end do

        ! Update global matrix structure
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))
        
        ! Resize solution vector
        call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                         rsolution, .false., .true.)

      case (FLOW_ALLCOUPLED)
        
        ! Resize all blocks
        do ivar = 1, NVAR2D
          do jvar = 1, NVAR2D
            call lsyssc_resizeMatrix(&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,jvar),&
                rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE), .false., .false., .true.)
          end do
        end do
        
        ! Update global matrix structure
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))

        ! Resize solution vector
        call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                         rsolution, .false., .true.)

      case DEFAULT
        call output_line('Unsupported block matrix format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators2D')
        call sys_halt()
      end select
      

    case DEFAULT
      call output_line('Unsupproted system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators2D')
      call sys_halt()
    end select

    ! Stop time measurement for coefficient matrices
    call stat_stopTimer(rtimer_assembly_coeff)

  end subroutine fs_reinitConstOperators2D

  !*****************************************************************************

!<subroutine>

  subroutine fs_reinitConstOperators3D(rcollection, rhadapt, rproblemLevel,&
                                       rsolution, ieltype, fcb_hadaptCallback)

!<description>
    ! This subroutine reinitializes the constant operators in 3D
!</description>

!<input>
    ! type of finite elements
    integer, intent(IN) :: ieltype
    
    ! callback routines
    include 'intf_hadaptcallback.inc'
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection

    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! Multigrid level for which grid adaptation should be performed
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! Solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution
!</inputoutput>
!</subroutine>

    ! local variables
    integer(I32), dimension(1) :: Ivalue
    integer :: neq,nedge,ivar,jvar

    ! Start time measurement for triangulation
    call stat_startTimer(rtimer_triangulation, STAT_TIMERSHORT)

    ! Generate raw mesh from adaptivity structure
    call hadapt_generateRawMesh(rhadapt, rproblemLevel%rtriangulation)

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw (rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
    
    ! Re-initialize the discretization structure   
    select case(ieltype)
    case (1)
      ! P1 finite elements
      call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                    EL_P1_3D, SPDISC_CUB_AUTOMATIC,&
                                    rproblemLevel%rtriangulation,&
                                    rproblemLevel%p_rproblem%rboundary)

    case (11)
      ! Q1 finite elements
      call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                    EL_Q1_3D, SPDISC_CUB_AUTOMATIC,&
                                    rproblemLevel%rtriangulation,&
                                    rproblemLevel%p_rproblem%rboundary)
      
    case (-1)
      ! mixed P1/Q1 finite elements
      call spdiscr_initDiscr_triquad(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                     EL_P1_3D, EL_Q1_3D, SPDISC_CUB_AUTOMATIC,&
                                     SPDISC_CUB_AUTOMATIC,&
                                     rproblemLevel%rtriangulation,&
                                     rproblemLevel%p_rproblem%rboundary)
      
    case DEFAULT
      call output_line('Unsupproted element type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators3D')
      call sys_halt()
    end select

    ! Stop time measurement for triangulation
    call stat_stopTimer(rtimer_triangulation)
    

    ! Start time measurement for coefficient matrices
    call stat_startTimer(rtimer_assembly_coeff, STAT_TIMERSHORT)

    ! Get sparsity pattern of the template matrix from collection
    Ivalue = 0
    call fcb_hadaptCallback(rcollection, -3, Ivalue, Ivalue)

    ! Resize scalar matrices MC, CX, CY, and CZ
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CZ),&
                             rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                             .false., .false., .true.)

    ! Generate the global coefficient matrices for the FE method
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                    DER_FUNC3D, DER_FUNC3D)
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                    DER_DERIV3D_X, DER_FUNC3D)
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                                    DER_DERIV3D_Y, DER_FUNC3D)
    call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_CZ),&
                                    DER_DERIV3D_Z, DER_FUNC3D)

    ! Perform mass lumping ML:=lumped(MC)
    call lsyssc_releaseMatrix  (rproblemLevel%Rmatrix(CNSE_MATRIX_ML))
    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                 LSYSSC_LUMP_DIAG)

    ! Resize AFC stabilisation for inviscid part
    neq   = rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ
    nedge = int(0.5*(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NA-&
                     rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ),I32)
    
    call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID), neq, nedge)


    ! Adjust system matrix
    select case(isystemFormat)
      
    case (SYSTEM_INTERLEAVEFORMAT)

      ! Relase pseudo block matrix
      call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))
      
      ! Resize scalar matrix
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NEQ,&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NCOLS,&
                               rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE)%NA,&
                               .false., .false., bforce=.true.)

      ! Create pseudo block matrix from global operator
      call lsysbl_createMatFromScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                      rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                      rproblemLevel%rdiscretisation)

      ! Resize solution vector
      call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                       rsolution, .false., .true.)

    case (SYSTEM_BLOCKFORMAT)
      
      ! What kind of global operator should be adopted?
      select case(isystemCoupling)
      case (FLOW_SEGREGATED)
          
        ! Resize diagonal blocks
        do ivar = 1, NVAR3D
          call lsyssc_resizeMatrix(&
              rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,ivar),&
              rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE), .false., .false., .true.)
        end do

        ! Update global matrix structure
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))
        
        ! Resize solution vector
        call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                         rsolution, .false., .true.)

      case (FLOW_ALLCOUPLED)
        
        ! Resize all blocks
        do ivar = 1, NVAR3D
          do jvar = 1, NVAR3D
            call lsyssc_resizeMatrix(&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,jvar),&
                rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE), .false., .false., .true.)
          end do
        end do
        
        ! Update global matrix structure
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))

        ! Resize solution vector
        call lsysbl_resizeVecBlockIndMat(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
            rsolution, .false., .true.)

      case DEFAULT
        call output_line('Unsupported block matrix format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators3D')
        call sys_halt()
      end select
      

    case DEFAULT
      call output_line('Unsupproted system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_reinitConstOperators3D')
      call sys_halt()
    end select

    ! Stop time measurement for coefficient matrices
    call stat_stopTimer(rtimer_assembly_coeff)

  end subroutine fs_reinitConstOperators3D
  
  !*****************************************************************************

!<subroutine>

  subroutine fcb_hadaptCallback1D(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), save :: rgraph
    type(t_matrixScalar), pointer, save :: rmatrix
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution

    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve matrix from collection and build sparsity-graph
      rmatrix => collct_getvalue_matsca(rcollection, "MATRIX_TEMPLATE")
      call grph_createGraphFromMatrix(rmatrix, rgraph)

      ! Retrieve solution vector from colletion and set pointer
      rsolution => collct_getvalue_vec(rcollection, "VECTOR_SOLUTION")
      call lsysbl_getbase_double(rsolution, p_Dsolution)


    case(HADAPT_OPR_DONECALLBACK)
      ! Release all data
      call grph_releaseGraph(rgraph)
      nullify(rmatrix, rsolution, p_Dsolution)

    case(-3)
      ! Retrieve matrix from collection and generate template matrix
      rmatrix => collct_getvalue_matsca(rcollection, "MATRIX_TEMPLATE")
      call grph_generateMatrix(rgraph, rmatrix)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      if (rsolution%NEQ .ne. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, Ivertices(1))

      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      p_Dsolution(Ivertices(1)) = 0.5_DP*&
          (p_Dsolution(Ivertices(2))+p_Dsolution(Ivertices(3)))


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph and solution
      if (Ivertices(2) .ne. 0) then
        call grph_removeVertex(rgraph,Ivertices(1), Ivertices(2))
        p_Dsolution(Ivertices(1)) = p_Dsolution(Ivertices(2))
      else
        call grph_removeVertex(rgraph, Ivertices(1))
        p_Dsolution(Ivertices(1)) = 0.0_DP
      end if

    case(HADAPT_OPR_REF_LINE2LINE)
      ! Delete broken edge (I1,I2) and add new edges (I1,I3),(I2,I3)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      
    case(HADAPT_OPR_CRS_2LINE1LINE)
      ! Delete broken edges (I1,I3) and (I2,I3) and add new edge (I1,I2)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      
    case DEFAULT
      call output_line('Invalid operation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fcb_hadaptCallback1D')
      call sys_halt()
    end select
  end subroutine fcb_hadaptCallback1D

  !*****************************************************************************

!<subroutine>

  subroutine fcb_hadaptCallback2D(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), save :: rgraph
    type(t_matrixScalar), pointer, save :: rmatrix
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: NEQ,ivar

    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve matrix from collection and build sparsity-graph
      rmatrix => collct_getvalue_matsca(rcollection, "MATRIX_TEMPLATE")
      call grph_createGraphFromMatrix(rmatrix, rgraph)

      ! Retrieve solution vector from colletion and set pointer
      rsolution => collct_getvalue_vec(rcollection, "VECTOR_SOLUTION")
      call lsysbl_getbase_double(rsolution, p_Dsolution)


    case(HADAPT_OPR_DONECALLBACK)
      ! Release all data
      call grph_releaseGraph(rgraph)
      nullify(rmatrix, rsolution, p_Dsolution)


    case(-3)
      ! Retrieve matrix from collection and generate template matrix
      rmatrix => collct_getvalue_matsca(rcollection, "MATRIX_TEMPLATE")
      call grph_generateMatrix(rgraph, rmatrix)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      if (rsolution%NEQ .ne. NVAR2D*Ivertices(1)) then
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1),&
                                        .false., .true.)
          call lsysbl_getbase_double(rsolution, p_Dsolution)

        case (SYSTEM_BLOCKFORMAT)
          call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1),&
                                        .false., .true.)
          call lsysbl_getbase_double(rsolution, p_Dsolution)

        case DEFAULT
          call output_line('Unsupported system format!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'fcb_hadaptCallback2D')
          call sys_halt()
        end select
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, Ivertices(1))

      ! Insert vertex into solution vector
      select case(isystemFormat)

      case (SYSTEM_INTERLEAVEFORMAT)
        if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
          call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1),&
                                        .false., .true.)
          call lsysbl_getbase_double(rsolution, p_Dsolution)
        end if

        do ivar = 1, NVAR2D
          p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) =&
              0.5_DP*( p_Dsolution((Ivertices(2)-1)*NVAR2D+ivar)+&
                       p_Dsolution((Ivertices(3)-1)*NVAR2D+ivar) )
        end do

      case (SYSTEM_BLOCKFORMAT)
        if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
          call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1),&
                                        .false., .true.)
          call lsysbl_getbase_double(rsolution, p_Dsolution)
        end if

        NEQ = rsolution%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_Dsolution((ivar-1)*NEQ+Ivertices(1)) =&
              0.5_DP*( p_Dsolution((ivar-1)*NEQ+Ivertices(2))+&
                       p_Dsolution((ivar-1)*NEQ+Ivertices(3)) )
        end do

      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'fcb_hadaptCallback2D')
        call sys_halt()
      end select
        

    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, Ivertices(1))

      ! Insert vertex into solution vector
      select case(isystemFormat)

      case (SYSTEM_INTERLEAVEFORMAT)
        if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
          call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1),&
                                        .false., .true.)
          call lsysbl_getbase_double(rsolution, p_Dsolution)
        end if

        do ivar = 1, NVAR2D
          p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) =&
              0.25_DP*( p_Dsolution((Ivertices(2)-1)*NVAR2D+ivar)+&
                        p_Dsolution((Ivertices(3)-1)*NVAR2D+ivar)+&
                        p_Dsolution((Ivertices(4)-1)*NVAR2D+ivar)+&
                        p_Dsolution((Ivertices(5)-1)*NVAR2D+ivar) )
        end do

      case (SYSTEM_BLOCKFORMAT)
        if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
          call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1),&
                                        .false., .true.)
          call lsysbl_getbase_double(rsolution, p_Dsolution)
        end if

        NEQ = rsolution%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_Dsolution((ivar-1)*NEQ+Ivertices(1)) =&
              0.25_DP*( p_Dsolution((ivar-1)*NEQ+Ivertices(2))+&
                        p_Dsolution((ivar-1)*NEQ+Ivertices(3))+&
                        p_Dsolution((ivar-1)*NEQ+Ivertices(4))+&
                        p_Dsolution((ivar-1)*NEQ+Ivertices(5)) )
        end do

      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'fcb_hadaptCallback2D')
        call sys_halt()
      end select


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph and solution
      if (Ivertices(2) .ne. 0) then
        call grph_removeVertex(rgraph,Ivertices(1), Ivertices(2))

        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          do ivar = 1, NVAR2D
            p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) =&
                p_Dsolution((Ivertices(2)-1)*NVAR2D+ivar)
          end do

        case (SYSTEM_BLOCKFORMAT)
          NEQ = rsolution%NEQ/NVAR2D
          do ivar = 1, NVAR2D
            p_Dsolution((ivar-1)*NEQ+Ivertices(1)) =&
                p_Dsolution((ivar-1)*NEQ+Ivertices(2))
          end do

        case DEFAULT
          call output_line('Unsupported system format!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'fcb_hadaptCallback2D')
          call sys_halt()
        end select

      else
        call grph_removeVertex(rgraph, Ivertices(1))

        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          do ivar = 1, NVAR2D
            p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) = 0.0_DP
          end do

        case (SYSTEM_BLOCKFORMAT)
          NEQ = rsolution%NEQ/NVAR2D
          do ivar = 1, NVAR2D
            p_Dsolution((ivar-1)*NEQ+Ivertices(1)) = 0.0_DP
          end do

        case DEFAULT
          call output_line('Unsupported system format!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'fcb_hadaptCallback2D')
          call sys_halt()
        end select
      end if


    case(HADAPT_OPR_REF_TRIA2TRIA)
      ! Delete broken edge (I1,I2) and add three new edges 
      ! (I1,I4), (I2,I4), and (I3,I4) if this is necessary
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))


    case(HADAPT_OPR_REF_TRIA3TRIA12)
      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Add new edges (I3,I4),(I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))


    case(HADAPT_OPR_REF_TRIA3TRIA23)
      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if

      ! Delete broken edges (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Add new edges (I1,I5),(I4,I5)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))


    case(HADAPT_OPR_REF_TRIA4TRIA)
      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Delete broken edge (I1,I3) and add new edges (I3,I6),(I1,I6)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(6))
      end if

      ! Add new edges I4,I5),(I4,I6), and (I5,I6)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))


    case(HADAPT_OPR_REF_QUAD2QUAD)
      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I6),(I4,I6)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      end if

      ! Add new edges (I1,I6),(I4,I5),(I2,I6),(I3,I5),(I5,I6)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))


    case(HADAPT_OPR_REF_QUAD3TRIA)
      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Add new edges (I3,I5),(I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))


    case(HADAPT_OPR_REF_QUAD4TRIA)
      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Add new edges (I4,I5),(I4,I6), and (I5,I6)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))


    case(HADAPT_OPR_REF_QUAD4QUAD)
      ! Delete broken edges (I1,I3) and (I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(7))
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
      end if

      ! Add new edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),
      ! (I2,I9),(I5,I6),(I3,I9),(I6,I7),(I4,I9), and (I7,I8)
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))


    case(HADAPT_OPR_CVT_TRIA2TRIA)
      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Delete broken edge (I1,I3) and add new edges (I1,I6),(I3,I6)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(6))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I5,I6),(I4,I6),(I4,I5)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))


    case(HADAPT_OPR_CVT_QUAD2QUAD)

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Delete broken edge (I1,I4) and add new edges  (I1,I8),(I4,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
      end if

      ! Delete broken edges (I5,I7),(I2,I7),(I3,I5),(I1,I7),(I4,I5) and
      ! add new edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I2,I9),
      ! (I3,I9),(I4,I9),(I5,I8),(I5,I6),(I6,I7),(I7,I8)
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))


    case(HADAPT_OPR_CVT_QUAD3TRIA)
      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(7))
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
      end if

      ! Delete broken edges (I5,I3) and (I5,I4) and add new edges
      ! (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),(I2,I9),(I5,I6)
      ! (I3,I9),I6,I7),(I4,I9), and (I7,I8)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))

    case(HADAPT_OPR_CVT_QUAD4TRIA)
      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(7))
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
      end if

      ! Delete broken edges (I4,I5),(I5,I6) and (I4,I6) and add new
      ! edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),(I2,I9),
      ! (I5,I6),(I3,I9),I6,I7),(I4,I9), and (I7,I8)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))   


    case(HADAPT_OPR_CRS_2TRIA1TRIA)
      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if
      
      ! Delete broken edge (I3,I4)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
      

    case(HADAPT_OPR_CRS_4TRIA1TRIA)
      ! Delete broken edges (I4,I5),(I4,I6), and (I5,I6)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA1)
      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I3,I4)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA2)
      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I1,I5)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA3)
      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I2,I6)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4QUAD1QUAD)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(5))
     
      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if

      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
 
      
    case(HADAPT_OPR_CRS_4QUAD2QUAD)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(5))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
      end if
      
      ! Add new edges (I5,I7),(I1,I7),(I4,I5),(I2,I7) and (I3,I5)
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      
      
    case(HADAPT_OPR_CRS_4QUAD3TRIA)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(5))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
      end if
      
      ! Add new edges (I3,I5) and (I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      
      
    case(HADAPT_OPR_CRS_4QUAD4TRIA)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7), and (I7,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I3,I5) and (I3,I8)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(8))
      
      
    case(HADAPT_OPR_CRS_2QUAD1QUAD)
      ! Delete broken edges (I5,I7), (I1,I7),(I4,I5),(I2,(I7) and (I3,I5)
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      
      
    case(HADAPT_OPR_CRS_2QUAD3TRIA)
      ! Delete broken edges (I5,I7),(I1,I7),(I4,I5),(I2,(I7) and (I3,I5)
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I3,I5) and (I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      
      
    case(HADAPT_OPR_CRS_3TRIA1QUAD)
      ! Delete broken edges (I3,I5) and (I4,I5)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      
      
    case(HADAPT_OPR_CRS_4TRIA1QUAD)
      ! Delete broken edges (I1,I6),(I1,I7) and (I6,I7)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      
      
    case(HADAPT_OPR_CRS_4TRIA3TRIA2)
      ! Delete broken edges (I1,I7) and (I6,I7)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edge (I4,I6)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      
      
    case(HADAPT_OPR_CRS_4TRIA3TRIA3)
      ! Delete broken edges (I1,I6) and (I6,I7)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Add new edge (I2,I7)
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(7))
      
      
    case DEFAULT
      call output_line('Invalid operation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fcb_hadaptCallback2D')
      call sys_halt()
    end select
  end subroutine fcb_hadaptCallback2D

  !*****************************************************************************

!<subroutine>

  subroutine fcb_hadaptCallback3D(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), save :: rgraph
    type(t_matrixScalar), pointer, save :: rmatrix
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution

    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve matrix from collection and build sparsity-graph
      rmatrix => collct_getvalue_matsca(rcollection, "MATRIX_TEMPLATE")
      call grph_createGraphFromMatrix(rmatrix, rgraph)

      ! Retrieve solution vector from colletion and set pointer
      rsolution => collct_getvalue_vec(rcollection, "VECTOR_SOLUTION")
      call lsysbl_getbase_double(rsolution, p_Dsolution)


    case(HADAPT_OPR_DONECALLBACK)
      ! Release all data
      call grph_releaseGraph(rgraph)
      nullify(rmatrix, rsolution, p_Dsolution)


    case(-3)
      ! Retrieve matrix from collection and generate template matrix
      rmatrix => collct_getvalue_matsca(rcollection, "MATRIX_TEMPLATE")
      call grph_generateMatrix(rgraph, rmatrix)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      if (rsolution%NEQ .ne. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
    

    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph and solution
      if (Ivertices(2) .ne. 0) then
        call grph_removeVertex(rgraph,Ivertices(1), Ivertices(2))
        p_Dsolution(Ivertices(1)) = p_Dsolution(Ivertices(2))
      else
        call grph_removeVertex(rgraph, Ivertices(1))
        p_Dsolution(Ivertices(1)) = 0.0_DP
      end if

    case DEFAULT
      call output_line('Invalid operation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fcb_hadaptCallback3D')
      call sys_halt()
    end select
  end subroutine fcb_hadaptCallback3D
end module euler_adaptation
