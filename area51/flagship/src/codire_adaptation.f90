!##############################################################################
!# ****************************************************************************
!# <name> codire_adaptation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# grid adaptation for a conservation law for a scalar variable
!#
!# The following routines are available:
!#
!# 1.) codire_performGridAdaptation
!#     -> perform adaptive grid refinement/coarsening
!#
!# 2.) codire_reinitConstOperators1D
!#     -> reinitialize the constant operators after grid adaptation in 1D
!#
!# 3.) codire_reinitConstOperators2D
!#     -> reinitialize the constant operators after grid adaptation in 2D
!#
!# 4.) codire_reinitConstOperators3D
!#     -> reinitialize the constant operators after grid adaptation in 3D
!#
!# </purpose>
!##############################################################################

module codire_adaptation

  use afcstabilisation
  use bilinearformevaluation
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

  use codire_basic
  use codire_callback
  use problem

  implicit none

  private

  public :: codire_performGridAdaptation
  public :: codire_reinitConstOperators1D
  public :: codire_reinitConstOperators2D
  public :: codire_reinitConstOperators3D
  
contains

  !*****************************************************************************
  
!<subroutine>

  subroutine codire_performGridAdaptation(rcollection, rhadapt, rproblemLevel, rsolution,&
                                       rgridIndicator, fcb_hadaptCallback)

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
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE), .true.)
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

  end subroutine codire_performGridAdaptation

  !*****************************************************************************
  
!<subroutine>

  subroutine codire_reinitConstOperators1D(rcollection, rhadapt, rproblemLevel, rsolution,&
                                        ieltype, fcb_hadaptCallback)

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
    integer :: neq,nedge,istabilisation
    logical :: bcompatible
    
    ! Start time measurement for triangulation
    call stat_startTimer(rtimer_triangulation, STAT_TIMERSHORT)

    ! Generate raw mesh from adaptivity structure
    call hadapt_generateRawMesh(rhadapt, rproblemLevel%rtriangulation)

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
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
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators1D')
      call sys_halt()
    end select
    
    ! Stop time measurement for triangulation
    call stat_stopTimer(rtimer_triangulation)
    

    ! Start time measurement for coefficient matrices
    call stat_startTimer(rtimer_assembly_coeff, STAT_TIMERSHORT)

    ! Get sparsity pattern of the template matrix from collection
    Ivalue = 0
    call fcb_hadaptCallback(rcollection, -3, Ivalue, Ivalue)

    
    ! Adjust Jacobian matrix accordingly (if required)
    if (rproblemLevel%Rmatrix(CDEQ_MATRIX_J)%cmatrixFormat .ne. LSYSSC_MATRIXUNDEFINED) then
            
      ! What kind of stabilization are we?
      istabilisation = max(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation,&
                           rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation)

      select case (istabilisation)
      case (AFCSTAB_GALERKIN,&
            AFCSTAB_UPWIND,&
            AFCSTAB_FEMFCT,&
            AFCSTAB_DMP)
        ! Adopt the standard sparsity pattern
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                 rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                 .false., .false., .true.)
        

      case (AFCSTAB_FEMGP, AFCSTAB_FEMTVD)
        if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iextendedJacobian .eq. 0) then
          ! Adopt the standard sparsity pattern
          call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                   .false., .false., .true.)
        else
          ! Extend the standard sparsity pattern
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
          call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                           rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
        end if
        
      case (AFCSTAB_SYMMETRIC)
        if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%iextendedJacobian .eq. 0) then
          ! Adopt the standard sparsity pattern
          call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                   .false., .false., .true.)
        else
          ! Extend the standard sparsity pattern
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
          call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                           rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
        end if
        

      case DEFAULT
        call output_line('Type of stabilization is not available',&
                         OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators1D')
        call sys_halt()
      end select
    end if


    ! Generate the global coefficient matrices for the FE method
    select case(abs(ivelocitytype))
    case (VELOCITY_ZERO)
      ! zero velocity, do nothing
      
    case (VELOCITY_CONSTANT,&
          VELOCITY_TIMEDEP,&
          VELOCITY_BURGERS1D,&
          VELOCITY_BUCKLEV1D)
      ! non-zero velocity, assemble coefficient matrices CX
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                      DER_DERIV1D_X, DER_FUNC1D)
      
      ! Resize AFC stabilisation for convective part
      neq   = rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ
      nedge = int(0.5*(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NA-&
                       rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ),I32)

      call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
                                       neq, nedge)

      ! Remove indicator for subdiagonal edges structure.
      ! If they are needed, then the have to be re-generated.
      rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iSpec = iand(&
          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iSpec,&
          not(AFCSTAB_SUBDIAGONALEDGES))

    case DEFAULT
      call output_line('Invalid type of velocity!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators1D')
    end select


    ! Generate the diffusion matrix
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! No diffusion, do nothing

    case (DIFFUSION_ISOTROPIC,&
          DIFFUSION_ANISOTROPIC)
      ! Isotropic diffusion, so generate the standard Laplace matrix and 
      ! scale it by the negative value of the diffusion coefficient.
      if (DdiffusionMatrix1D(1,1) > 0.0_DP) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                 rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                 .false., .false., .true.)

        call stdop_assembleLaplaceMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .true.,&
                                         -DdiffusionMatrix2D(1,1))
      end if
   
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators2d')
      call sys_halt()
    end select
    

    ! Adjust consistent mass matrix MC
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                             rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                             .false., .false., .true.)

    call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                     DER_FUNC2D, DER_FUNC2D)
    
    ! Perform mass lumping ML:=lumped(MC)
    call lsyssc_releaseMatrix  (rproblemLevel%Rmatrix(CDEQ_MATRIX_ML))
    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)

    call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                 LSYSSC_LUMP_DIAG)


    ! Adjust matrix A and low-order operator L
    call lsyssc_isMatrixCompatible(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_L), bcompatible)
    if (bcompatible) then
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)
    else
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                               .false., .false., .true.)

      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)
    end if

    ! Resize solution vector
    call lsysbl_resizeVectorBlock(rsolution,&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ, .false.)

    ! Mark velocity for update
    bvelocityUpdate = .true.
    
    ! Stop time measurement for coefficient matrices
    call stat_stopTimer(rtimer_assembly_coeff)

  end subroutine codire_reinitConstOperators1D

  !*****************************************************************************
  
!<subroutine>

  subroutine codire_reinitConstOperators2D(rcollection, rhadapt, rproblemLevel, rsolution,&
                                        ieltype, fcb_hadaptCallback)

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
    type(t_bilinearForm) :: rform
    integer(I32), dimension(1) :: Ivalue
    integer :: neq,nedge,istabilisation
    logical :: bcompatible
    
    ! Start time measurement for triangulation
    call stat_startTimer(rtimer_triangulation, STAT_TIMERSHORT)

    ! Generate raw mesh from adaptivity structure
    call hadapt_generateRawMesh(rhadapt, rproblemLevel%rtriangulation)

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
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
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators2D')
      call sys_halt()
    end select
    
    ! Stop time measurement for triangulation
    call stat_stopTimer(rtimer_triangulation)
    

    ! Start time measurement for coefficient matrices
    call stat_startTimer(rtimer_assembly_coeff, STAT_TIMERSHORT)

    ! Get sparsity pattern of the template matrix from collection
    Ivalue = 0
    call fcb_hadaptCallback(rcollection, -3, Ivalue, Ivalue)

    
    ! Adjust Jacobian matrix accordingly (if required)
    if (rproblemLevel%Rmatrix(CDEQ_MATRIX_J)%cmatrixFormat .ne. LSYSSC_MATRIXUNDEFINED) then
            
      ! What kind of stabilization are we?
      istabilisation = max(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation,&
                           rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation)

      select case (istabilisation)
      case (AFCSTAB_GALERKIN,&
            AFCSTAB_UPWIND,&
            AFCSTAB_FEMFCT,&
            AFCSTAB_DMP)
        ! Adopt the standard sparsity pattern
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                 rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                 .false., .false., .true.)
        

      case (AFCSTAB_FEMGP, AFCSTAB_FEMTVD)
        if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iextendedJacobian .eq. 0) then
          ! Adopt the standard sparsity pattern
          call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                   .false., .false., .true.)
        else
          ! Extend the standard sparsity pattern
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
          call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                           rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
        end if
        
      case (AFCSTAB_SYMMETRIC)
        if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%iextendedJacobian .eq. 0) then
          ! Adopt the standard sparsity pattern
          call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                   .false., .false., .true.)
        else
          ! Extend the standard sparsity pattern
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
          call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                           rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
        end if
        

      case DEFAULT
        call output_line('Type of stabilization is not available',&
                         OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators2D')
        call sys_halt()
      end select
    end if


    ! Generate the global coefficient matrices for the FE method
    select case(abs(ivelocitytype))
    case (VELOCITY_ZERO)
      ! zero velocity, do nothing
      
    case (VELOCITY_CONSTANT, VELOCITY_TIMEDEP,&
          VELOCITY_BURGERS_SPACETIME,&
          VELOCITY_BUCKLEV_SPACETIME,&
          VELOCITY_BURGERS2D)
      ! non-zero velocity, assemble coefficient matrices CX and CY
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                      DER_DERIV2D_X, DER_FUNC2D)
      
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                                      DER_DERIV2D_Y, DER_FUNC2D)

      ! Resize AFC stabilisation for convective part
      neq   = rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ
      nedge = int(0.5*(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NA-&
                       rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ),I32)

      call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
                                       neq, nedge)

      ! Remove indicator for subdiagonal edges structure.
      ! If they are needed, then the have to be re-generated.
      rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iSpec = iand(&
          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iSpec,&
          not(AFCSTAB_SUBDIAGONALEDGES))

    case DEFAULT
      call output_line('Invalid type of velocity!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators2D')
    end select


    ! Generate the diffusion matrix
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! No diffusion, do nothing

    case (DIFFUSION_ISOTROPIC)
      ! Isotropic diffusion, so generate the standard Laplace matrix and 
      ! scale it by the negative value of the diffusion coefficient.
      if (DdiffusionMatrix2D(1,1) > 0.0_DP) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                 rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                 .false., .false., .true.)

        call stdop_assembleLaplaceMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .true.,&
                                         -DdiffusionMatrix2D(1,1))
      end if
      
    case (DIFFUSION_ANISOTROPIC)
      
      ! For anisotropic diffusion, things are slightly more complicated.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      rform%itermCount = 4
      rform%Idescriptors(1,1) = DER_DERIV2D_X
      rform%Idescriptors(2,1) = DER_DERIV2D_X
      
      rform%Idescriptors(1,2) = DER_DERIV2D_X
      rform%Idescriptors(2,2) = DER_DERIV2D_Y
      
      rform%Idescriptors(1,3) = DER_DERIV2D_Y
      rform%Idescriptors(2,3) = DER_DERIV2D_X
      
      rform%Idescriptors(1,4) = DER_DERIV2D_Y
      rform%Idescriptors(2,4) = DER_DERIV2D_Y
      
      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -DdiffusionMatrix2D(1,1)
      rform%Dcoefficients(2)  = -DdiffusionMatrix2D(1,2)
      rform%Dcoefficients(3)  = -DdiffusionMatrix2D(2,1)
      rform%Dcoefficients(4)  = -DdiffusionMatrix2D(2,2)
      
      ! Now we can build the matrix entries
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call bilf_buildMatrixScalar(rform, .true.,&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_S))

      ! Resize AFC stabilisation for diffusive part
      neq   = rproblemLevel%Rmatrix(CDEQ_MATRIX_S)%NEQ
      nedge = int(0.5*(rproblemLevel%Rmatrix(CDEQ_MATRIX_S)%NA-&
                       rproblemLevel%Rmatrix(CDEQ_MATRIX_S)%NEQ),I32)

      call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION),&
                                       neq, nedge)

      ! Remove indicator for subdiagonal edges structure.
      ! If they are needed, then the have to be re-generated.
      rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%iSpec=iand(&
          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%iSpec, not(AFCSTAB_SUBDIAGONALEDGES))

    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators2D')
      call sys_halt()
    end select
    

    ! Adjust consistent mass matrix MC
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                             rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                             .false., .false., .true.)

    call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                     DER_FUNC2D, DER_FUNC2D)
    
    ! Perform mass lumping ML:=lumped(MC)
    call lsyssc_releaseMatrix  (rproblemLevel%Rmatrix(CDEQ_MATRIX_ML))
    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)

    call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                 LSYSSC_LUMP_DIAG)


    ! Adjust matrix A and low-order operator L
    call lsyssc_isMatrixCompatible(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_L),bcompatible)
    if (bcompatible) then
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)
    else
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                               .false., .false., .true.)

      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)
    end if

    ! Resize solution vector
    call lsysbl_resizeVectorBlock(rsolution,&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ, .false.)

    ! Mark velocity for update
    bvelocityUpdate = .true.
    
    ! Stop time measurement for coefficient matrices
    call stat_stopTimer(rtimer_assembly_coeff)

  end subroutine codire_reinitConstOperators2D

  !*****************************************************************************
  
!<subroutine>

  subroutine codire_reinitConstOperators3D(rcollection, rhadapt, rproblemLevel, rsolution,&
                                        ieltype, fcb_hadaptCallback)

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
    type(t_bilinearForm) :: rform
    integer(I32), dimension(1) :: Ivalue
    integer :: neq,nedge
    logical :: bcompatible
    
    ! Start time measurement for triangulation
    call stat_startTimer(rtimer_triangulation, STAT_TIMERSHORT)

    ! Generate raw mesh from adaptivity structure
    call hadapt_generateRawMesh(rhadapt, rproblemLevel%rtriangulation)

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
                                      rproblemLevel%p_rproblem%rboundary)
    
    ! Re-initialize the discretization structure   
    select case(ieltype)
    case (1)
      ! P1 finite elements
      call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                    EL_E001_3D, SPDISC_CUB_AUTOMATIC,&
                                    rproblemLevel%rtriangulation,&
                                    rproblemLevel%p_rproblem%rboundary)
      
    case (11)
      ! Q1 finite elements
      call spdiscr_initDiscr_simple(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                    EL_E011_3D, SPDISC_CUB_AUTOMATIC,&
                                    rproblemLevel%rtriangulation,&
                                    rproblemLevel%p_rproblem%rboundary)
      
    case (-1)
      ! mixed P1/Q1 finite elements
      call spdiscr_initDiscr_triquad(rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                     EL_E001_3D, EL_E011_3D, SPDISC_CUB_AUTOMATIC,&
                                     SPDISC_CUB_AUTOMATIC,&
                                     rproblemLevel%rtriangulation,&
                                     rproblemLevel%p_rproblem%rboundary)
      
    case DEFAULT
      call output_line('Unsupproted element type!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators3D')
      call sys_halt()
    end select
    
    ! Stop time measurement for triangulation
    call stat_stopTimer(rtimer_triangulation)
    

    ! Start time measurement for coefficient matrices
    call stat_startTimer(rtimer_assembly_coeff, STAT_TIMERSHORT)

    ! Get sparsity pattern of the template matrix from collection
    Ivalue = 0
    call fcb_hadaptCallback(rcollection, -3, Ivalue, Ivalue)

    
    ! Adjust Jacobian matrix accordingly (if required)
    if (rproblemLevel%Rmatrix(CDEQ_MATRIX_J)%cmatrixFormat .ne. LSYSSC_MATRIXUNDEFINED) then
      
      ! What kind of stabilization are we?
      select case (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
      case (AFCSTAB_GALERKIN,&
            AFCSTAB_UPWIND,&
            AFCSTAB_FEMFCT)
        ! Adopt the standard sparsity pattern
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                 rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                 .false., .false., .true.)
        
      case (AFCSTAB_FEMGP,&
            AFCSTAB_FEMTVD)
        if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iextendedJacobian .eq. 0) then
          ! Adopt the standard sparsity pattern
          call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                   .false., .false., .true.)
          
        else
          ! Extend the standard sparsity pattern
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
          call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                           rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
        end if
        
      case DEFAULT
        call output_line('Type of stabilization is not available',&
                         OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators3D')
        call sys_halt()
      end select
    end if


    ! Generate the global coefficient matrices for the FE method
    select case(abs(ivelocitytype))
    case (VELOCITY_ZERO)
      ! zero velocity, do nothing
      
    case (VELOCITY_CONSTANT,&
          VELOCITY_TIMEDEP,&
          VELOCITY_BURGERS_SPACETIME,&
          VELOCITY_BUCKLEV_SPACETIME,&
          VELOCITY_BURGERS2D)
      ! non-zero velocity, assemble coefficient matrices CX and CY
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                      DER_DERIV3D_X, DER_FUNC3D)
      
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                                      DER_DERIV3D_Y, DER_FUNC3D)

      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CZ),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CZ),&
                                      DER_DERIV3D_Z, DER_FUNC3D)

      ! Resize AFC stabilisation for convective part
      neq   = rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ
      nedge = int(0.5*(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NA-&
                       rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ),I32)

      call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
                                       neq, nedge)

    case DEFAULT
      call output_line('Invalid type of velocity!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators3D')
    end select


    ! Generate the diffusion matrix
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! No diffusion, do nothing

    case (DIFFUSION_ISOTROPIC)
      ! Isotropic diffusion, so generate the standard Laplace matrix and 
      ! scale it by the negative value of the diffusion coefficient.
      if (DdiffusionMatrix3D(1,1) > 0.0_DP) then
        call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                 rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                 .false., .false., .true.)

        call stdop_assembleLaplaceMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .true.,&
                                         -DdiffusionMatrix3D(1,1))
      end if
      
    case (DIFFUSION_ANISOTROPIC)
      
      ! For anisotropic diffusion, things are slightly more complicated.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 3D.
      rform%itermCount = 9
      rform%Idescriptors(1,1) = DER_DERIV3D_X
      rform%Idescriptors(2,1) = DER_DERIV3D_X
      
      rform%Idescriptors(1,2) = DER_DERIV3D_X
      rform%Idescriptors(2,2) = DER_DERIV3D_Y
     
      rform%Idescriptors(1,3) = DER_DERIV3D_X
      rform%Idescriptors(2,3) = DER_DERIV3D_Z
      
      rform%Idescriptors(1,4) = DER_DERIV3D_Y
      rform%Idescriptors(2,4) = DER_DERIV3D_X

      rform%Idescriptors(1,5) = DER_DERIV3D_Y
      rform%Idescriptors(2,5) = DER_DERIV3D_Y

      rform%Idescriptors(1,6) = DER_DERIV3D_Y
      rform%Idescriptors(2,6) = DER_DERIV3D_Z

      rform%Idescriptors(1,7) = DER_DERIV3D_Z
      rform%Idescriptors(2,7) = DER_DERIV3D_X

      rform%Idescriptors(1,8) = DER_DERIV3D_Z
      rform%Idescriptors(2,8) = DER_DERIV3D_Y

      rform%Idescriptors(1,9) = DER_DERIV3D_Z
      rform%Idescriptors(2,9) = DER_DERIV3D_Z

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -DdiffusionMatrix3D(1,1)
      rform%Dcoefficients(2)  = -DdiffusionMatrix3D(1,2)
      rform%Dcoefficients(3)  = -DdiffusionMatrix3D(1,3)
      rform%Dcoefficients(4)  = -DdiffusionMatrix3D(2,1)
      rform%Dcoefficients(5)  = -DdiffusionMatrix3D(2,2)
      rform%Dcoefficients(6)  = -DdiffusionMatrix3D(2,3)
      rform%Dcoefficients(7)  = -DdiffusionMatrix3D(3,1)
      rform%Dcoefficients(8)  = -DdiffusionMatrix3D(3,2)
      rform%Dcoefficients(9)  = -DdiffusionMatrix3D(3,3)
      
      ! Now we can build the matrix entries
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call bilf_buildMatrixScalar(rform, .true.,&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_S))

      ! Resize AFC stabilisation for diffusive part
      neq   = rproblemLevel%Rmatrix(CDEQ_MATRIX_S)%NEQ
      nedge = int(0.5*(rproblemLevel%Rmatrix(CDEQ_MATRIX_S)%NA-&
                       rproblemLevel%Rmatrix(CDEQ_MATRIX_S)%NEQ),I32)

      call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION),&
                                       neq, nedge)

    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'codire_reinitConstOperators3D')
      call sys_halt()
    end select
    

    ! Adjust consistent mass matrix MC
    call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                             rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                             .false., .false., .true.)

    call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                     DER_FUNC3D, DER_FUNC3D)
    
    ! Perform mass lumping ML:=lumped(MC)
    call lsyssc_releaseMatrix  (rproblemLevel%Rmatrix(CDEQ_MATRIX_ML))
    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)

    call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                 LSYSSC_LUMP_DIAG)


    ! Adjust matrix A and low-order operator L
    call lsyssc_isMatrixCompatible(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_L),bcompatible)
    if (bcompatible) then
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)

      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)
    else
      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                               .false., .false., .true.)

      call lsyssc_resizeMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                               .false., .false., .true.)
    end if

    ! Resize solution vector
    call lsysbl_resizeVectorBlock(rsolution,&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE)%NEQ, .false.)

    ! Mark velocity for update
    bvelocityUpdate = .true.
    
    ! Stop time measurement for coefficient matrices
    call stat_stopTimer(rtimer_assembly_coeff)

  end subroutine codire_reinitConstOperators3D
end module codire_adaptation
