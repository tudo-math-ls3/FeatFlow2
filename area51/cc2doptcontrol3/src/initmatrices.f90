!##############################################################################
!# ****************************************************************************
!# <name> initmatrices </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to allocate and assemble template matrices
!# for the optimal control problem.
!#
!# Routines in this module:
!#
!# 1.) inmat_initSpaceLevel
!#     -> Basic initialisaion of a template assembly structure
!#        with information for the optimal control problem.
!#
!# 2.) inmat_allocStaticMatrices
!#     -> Allocate memory for matrices a template assembly structure
!#
!# 3.) inmat_releaseStaticMatrices
!#     -> Release matrices in a template assembly structure
!#
!# 4.) inmat_generateStaticMatrices
!#     -> Calculate matrices in a template assembly structure
!#
!# 5.) inmat_initStaticAsmTemplHier
!#     -> Create a level hierarchy of template assembly structures
!#
!# 6.) inmat_doneStaticAsmTemplHier
!#     -> Release a level hierarchy of template assembly structures
!#
!# 7.) inmat_calcStaticLevelAsmHier
!#     -> Calculates the matrices in a template assembly hierarchy
!# </purpose>
!##############################################################################

module initmatrices

  use fsystem
  use storage
  use genoutput
  use element
  use linearalgebra
  use boundary
  use triangulation
  use spatialdiscretisation
  use collection
  use derivatives
  use bilinearformevaluation
  use stdoperators
  use linearsystemscalar
  use linearsystemblock
  use convection
  use matrixmodification

  use fespacehierarchybase
  use fespacehierarchy
  
  use constantsdiscretisation
  use structuresdiscretisation
  use assemblytemplates
  
  implicit none
  
  private
  
  public :: inmat_initSpaceLevel
  public :: inmat_doneSpaceLevel
  public :: inmat_allocStaticMatrices
  public :: inmat_releaseStaticMatrices
  public :: inmat_generateStaticMatrices
  public :: inmat_initStaticAsmTemplHier
  public :: inmat_doneStaticAsmTemplHier
  public :: inmat_calcStaticLevelAsmHier
  public :: imnat_getL2PrjMatrix
  public :: inmat_printTmplHierStatistic
  
contains

  ! ***************************************************************************

!<subrotine>
  
  subroutine inmat_initSpaceLevel(rstaticAsmTemplates,rtriangulation,&
      rsettingsSpaceDiscr,rphysics,&
      rdiscretisationPrimal,rdiscretisationDual,rdiscretisationControl)

!<description>
  ! Basic initialisation of a t_staticSpaceAsmTemplates structure.
!</description>

!<input>
  ! Underlying triangulation
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! Settings controlling the spatial discretisation (cubature)
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Physics of the problem
  type(t_settings_physics), intent(in) :: rphysics

  ! Discretisation of the primal space
  type(t_blockDiscretisation), intent(in), target :: rdiscretisationPrimal

  ! Discretisation of the dual space
  type(t_blockDiscretisation), intent(in), target :: rdiscretisationDual

  ! Discretisation of the control space
  type(t_blockDiscretisation), intent(in), target :: rdiscretisationControl
!</input>

!<output>
  ! The Assembly template structure to initialise.
  type(t_staticSpaceAsmTemplates), intent(out) :: rstaticAsmTemplates
!</output>

!</subroutine>

    ! Just set pointers.
    rstaticAsmTemplates%p_rtriangulation => rtriangulation
    rstaticAsmTemplates%p_rdiscr => rdiscretisationPrimal%RspatialDiscr(1)
    
    ! Initialise the cubature information structures for the assembly
    ! of the matrices and vectors.
    
    call spdiscr_createDefCubStructure(&
        rstaticAsmTemplates%p_rdiscr,&
        rstaticAsmTemplates%rcubatureInfo,&
        rsettingsSpaceDiscr%icubStokes)

    call spdiscr_createDefCubStructure(&
        rstaticAsmTemplates%p_rdiscr,&
        rstaticAsmTemplates%rcubatureInfoMass,&
        rsettingsSpaceDiscr%icubMass)

    call spdiscr_createDefCubStructure(&
        rstaticAsmTemplates%p_rdiscr,&
        rstaticAsmTemplates%rcubatureInfoRHS,&
        rsettingsSpaceDiscr%icubF)

    select case (rphysics%cequation)
    
    ! ---------------------------------------------------------------
    ! Stokes/Navier Stokes.
    ! ---------------------------------------------------------------
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      
      ! Additional initialisation
      rstaticAsmTemplates%p_rdiscrPressure => rdiscretisationPrimal%RspatialDiscr(3)
      
      call spdiscr_createDefCubStructure(&
          rstaticAsmTemplates%p_rdiscrPressure,&
          rstaticAsmTemplates%rcubatureInfoDiv,&
          rsettingsSpaceDiscr%icubB)
      
      call spdiscr_createDefCubStructure(&
          rstaticAsmTemplates%p_rdiscrPressure,&
          rstaticAsmTemplates%rcubatureInfoMassPressure,&
          rsettingsSpaceDiscr%icubMass)

      call spdiscr_createDefCubStructure(&
          rstaticAsmTemplates%p_rdiscrPressure,&
          rstaticAsmTemplates%rcubatureInfoRHScontinuity,&
          rsettingsSpaceDiscr%icubF)

    end select
    
  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine inmat_doneSpaceLevel(rstaticAsmTemplates)

!<description>
  ! Final cleanup of a t_staticSpaceAsmTemplates structure.
!</description>

!<inputoutput>
  ! The Assembly template structure to clean up
  type(t_staticSpaceAsmTemplates), intent(inout) :: rstaticAsmTemplates
!</inputoutput>

!</subroutine>

    ! Remove the pointers
    nullify(rstaticAsmTemplates%p_rtriangulation)
    nullify(rstaticAsmTemplates%p_rdiscr)
    nullify(rstaticAsmTemplates%p_rdiscrPressure)
    
    ! Release cubature information structures
    call spdiscr_releaseCubStructure(rstaticAsmTemplates%rcubatureInfo)
    call spdiscr_releaseCubStructure(rstaticAsmTemplates%rcubatureInfoMass)
    call spdiscr_releaseCubStructure(rstaticAsmTemplates%rcubatureInfoRHS)
    
    if (rstaticAsmTemplates%rcubatureInfoDiv%ninfoBlockCount .ne. 0) &
      call spdiscr_releaseCubStructure(rstaticAsmTemplates%rcubatureInfoDiv)
    if (rstaticAsmTemplates%rcubatureInfoMassPressure%ninfoBlockCount .ne. 0) &
      call spdiscr_releaseCubStructure(rstaticAsmTemplates%rcubatureInfoMassPressure)
    if (rstaticAsmTemplates%rcubatureInfoRHScontinuity%ninfoBlockCount .ne. 0) &
      call spdiscr_releaseCubStructure(rstaticAsmTemplates%rcubatureInfoRHScontinuity)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine inmat_allocStaticMatrices (rstaticAsmTemplates,rphysics)
  
!<description>
  ! Allocates memory and generates the structure of all static matrices
  ! in rstaticAsmTemplates.
!</description>

!<input>
  ! Physics of the problem
  type(t_settings_physics), intent(in) :: rphysics
!</input>

!<inputoutput>
  ! A t_staticLevelInfo structure. The static matrices in this structure are generated.
  type(t_staticSpaceAsmTemplates), intent(inout) :: rstaticAsmTemplates
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: cmatBuildType
  
    ! When the jump stabilisation is used, we have to create an extended
    ! matrix stencil for the velocity template matrices!
    cmatBuildType = BILF_MATC_ELEMENTBASED
    
    select case (rphysics%cequation)
    
    ! ---------------------------------------------------------------
    ! Stokes/Navier Stokes.
    ! ---------------------------------------------------------------
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      
      ! -----------------------------------------------------------------------
      ! General matrix templates
      ! -----------------------------------------------------------------------
      ! Get a pointer to the template FEM matrix. This is used for the
      ! Laplace/Stokes matrix and probably for the mass matrix.
      ! It defines KCOL/KLD.
      
      ! Create the matrix structure for the velocity FEM space.
      call bilf_createMatrixStructure (&
          rstaticAsmTemplates%p_rdiscr,LSYSSC_MATRIX9,&
          rstaticAsmTemplates%rmatrixTemplateFEM,cconstrType=cmatBuildType)
          
      ! In case the element-based routine is used to create the matrices,
      ! the "offdiagonal" matrices have the same structure. If we used
      ! the edge-based construction routine, the "offdiagonal" matrices
      ! can still be constructed with a smaller stencil.
      if (cmatBuildType .ne. BILF_MATC_ELEMENTBASED) then
        call bilf_createMatrixStructure (&
            rstaticAsmTemplates%p_rdiscr,LSYSSC_MATRIX9,&
            rstaticAsmTemplates%rmatrixTemplateFEMOffdiag,cconstrType=BILF_MATC_ELEMENTBASED)
      else
        call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEM,&
            rstaticAsmTemplates%rmatrixTemplateFEMOffdiag,LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
      end if

      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
          rstaticAsmTemplates%p_rdiscrPressure,LSYSSC_MATRIX9,&
          rstaticAsmTemplates%rmatrixTemplateFEMPressure)

      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
          rstaticAsmTemplates%p_rdiscrPressure,LSYSSC_MATRIX9,&
          rstaticAsmTemplates%rmatrixTemplateGradient,&
          rstaticAsmTemplates%p_rdiscr)
                
      ! Transpose the B-structure to get the matrix template for the
      ! divergence matrices.
      call lsyssc_transposeMatrix (rstaticAsmTemplates%rmatrixTemplateGradient,&
          rstaticAsmTemplates%rmatrixTemplateDivergence,LSYSSC_TR_STRUCTURE)
      
      ! Ok, now we use the matrices from above to create the actual submatrices.
      
      ! -----------------------------------------------------------------------
      ! Laplace/Stokes matrix
      ! -----------------------------------------------------------------------
      ! Connect the Stokes matrix to the template FEM matrix such that they
      ! use the same structure.
      !
      ! Don"t create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEM,&
          rstaticAsmTemplates%rmatrixLaplace,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      
      ! Allocate memory for the entries; don"t initialise the memory.
      call lsyssc_allocEmptyMatrix (&
          rstaticAsmTemplates%rmatrixLaplace,LSYSSC_SETM_UNDEFINED)
      
      ! -----------------------------------------------------------------------
      ! B-matrices
      ! -----------------------------------------------------------------------
      ! Create matrices for the gradient term.
      ! Don"t create a content array yet, it will be created by
      ! the assembly routines later.
      ! Allocate memory for the entries; don"t initialise the memory.
      
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateGradient,&
          rstaticAsmTemplates%rmatrixB1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateGradient,&
          rstaticAsmTemplates%rmatrixB2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixB1,LSYSSC_SETM_UNDEFINED)
      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixB2,LSYSSC_SETM_UNDEFINED)

      ! -----------------------------------------------------------------------
      ! D-matrices
      ! -----------------------------------------------------------------------
      ! Set up memory for the divergence matrices D1 and D2.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateDivergence,&
          rstaticAsmTemplates%rmatrixD1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateDivergence,&
          rstaticAsmTemplates%rmatrixD2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixD1,LSYSSC_SETM_UNDEFINED)
      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixD2,LSYSSC_SETM_UNDEFINED)
      
      ! -----------------------------------------------------------------------
      ! D^T-matrices matrix
      ! -----------------------------------------------------------------------
      ! The D1^T and D2^T matrices are by default the same as B1 and B2.
      ! These matrices may be different for special VANCA variants if
      ! B1 and B2 is different from D1 and D2 (which is actually a rare case).
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixB1,&
          rstaticAsmTemplates%rmatrixD1T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                  
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixB2,&
          rstaticAsmTemplates%rmatrixD2T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      ! -----------------------------------------------------------------------
      ! Mass matrices
      ! -----------------------------------------------------------------------
      ! Generate mass matrix. The matrix has basically the same structure as
      ! our template FEM matrix, so we can take that.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEM,&
          rstaticAsmTemplates%rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixMass,LSYSSC_SETM_UNDEFINED)

      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEMPressure,&
          rstaticAsmTemplates%rmatrixMassPressure,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixMassPressure,&
          LSYSSC_SETM_UNDEFINED)

      ! Mass matrix with extended structure for pure Dirichlet problems
      ! and UMFPACK.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEMPressure,&
          rstaticAsmTemplates%rmatrixMassPressureExtStruc,LSYSSC_DUP_COPY,LSYSSC_DUP_REMOVE)
      call mmod_expandToFullRow (rstaticAsmTemplates%rmatrixMassPressureExtStruc,1)
      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixMassPressureExtStruc,&
          LSYSSC_SETM_UNDEFINED)
          
    ! ---------------------------------------------------------------
    ! Heat equation
    ! ---------------------------------------------------------------
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)

      ! -----------------------------------------------------------------------
      ! General matrix templates
      ! -----------------------------------------------------------------------
      ! Get a pointer to the template FEM matrix. This is used for the
      ! Laplace/Stokes matrix and probably for the mass matrix.
      ! It defines KCOL/KLD.
      
      ! Create the matrix structure for the main FEM space
      call bilf_createMatrixStructure (&
          rstaticAsmTemplates%p_rdiscr,LSYSSC_MATRIX9,&
          rstaticAsmTemplates%rmatrixTemplateFEM,cconstrType=cmatBuildType)

      ! Expand the first row to full to support pure Neumann problems
      ! call mmod_expandToFullRow (rstaticAsmTemplates%rmatrixTemplateFEM,1)

      ! Ok, now we use the matrices from above to create the actual submatrices.
      
      ! -----------------------------------------------------------------------
      ! Laplace/Stokes matrix
      ! -----------------------------------------------------------------------
      ! Connect the Stokes matrix to the template FEM matrix such that they
      ! use the same structure.
      !
      ! Don"t create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEM,&
          rstaticAsmTemplates%rmatrixLaplace,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      
      ! Allocate memory for the entries; don"t initialise the memory.
      call lsyssc_allocEmptyMatrix (&
          rstaticAsmTemplates%rmatrixLaplace,LSYSSC_SETM_UNDEFINED)
      
      ! -----------------------------------------------------------------------
      ! Mass matrix
      ! -----------------------------------------------------------------------
      ! Generate mass matrix. The matrix has basically the same structure as
      ! our template FEM matrix, so we can take that.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEM,&
          rstaticAsmTemplates%rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_allocEmptyMatrix (rstaticAsmTemplates%rmatrixMass,LSYSSC_SETM_UNDEFINED)

    end select
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine inmat_releaseStaticMatrices (rstaticAsmTemplates)
  
!<description>
  ! Releases all static matrices from rstaticAsmTemplates structure.
!</description>

!<inputoutput>
  ! A t_staticLevelInfo structure to be cleaned up.
  type(t_staticSpaceAsmTemplates), intent(inout) :: rstaticAsmTemplates
!</inputoutput>

!</subroutine>

    ! If there is an existing mass matrix, release it.
    if (rstaticAsmTemplates%rmatrixMass%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixMass)
    if (rstaticAsmTemplates%rmatrixMassPressure%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixMassPressure)
    if (rstaticAsmTemplates%rmatrixMassPressureExtStruc%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixMassPressureExtStruc)
    if (rstaticAsmTemplates%rmatrixMassPressureLumpInt%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixMassPressureLumpInt)
    if (rstaticAsmTemplates%rmatrixMassLumpInt%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixMassLumpInt)

    ! Release Stokes, B1, B2,... matrices
    if (rstaticAsmTemplates%rmatrixD2T%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixD2T)
    if (rstaticAsmTemplates%rmatrixD1T%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixD1T)
    if (rstaticAsmTemplates%rmatrixD2%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixD2)
    if (rstaticAsmTemplates%rmatrixD1%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixD1)
    if (rstaticAsmTemplates%rmatrixB2%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixB2)
    if (rstaticAsmTemplates%rmatrixB1%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixB1)
    if (rstaticAsmTemplates%rmatrixLaplace%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixLaplace)
    
    ! Release the template matrices. This is the point, where the
    ! memory of the matrix structure is released.
    if (rstaticAsmTemplates%rmatrixTemplateDivergence%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixTemplateDivergence)
    if (rstaticAsmTemplates%rmatrixTemplateGradient%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixTemplateGradient)
    if (rstaticAsmTemplates%rmatrixTemplateFEM%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixTemplateFEM)
    if (rstaticAsmTemplates%rmatrixTemplateFEMPressure%NEQ .ne. 0) &
      call lsyssc_releaseMatrix (rstaticAsmTemplates%rmatrixTemplateFEMPressure)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine inmat_generateStaticMatrices (rstaticAsmTemplates,rsettingsSpaceDiscr,rphysics)
  
!<description>
  ! Calculates entries of all static matrices (Stokes, B-matrices,...)
  ! in the specified problem structure, i.e. the entries of all matrices
  ! that do not change during the computation or which serve as a template for
  ! generating other matrices.
  !
  ! Memory for those matrices must have been allocated before with
  ! inmat_allocStaticMatrices!
!</description>

!<input>
  ! Settings controlling the spatial discretisation (stabilisation parameters).
  ! This must coincide with the structure passed to inmat_initSpaceLevel.
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Physics of the problem
  type(t_settings_physics), intent(in) :: rphysics
!</input>

!<inputoutput>
  ! A t_staticSpaceAsmTemplates structure. The static matrices in this structure are generated.
  type(t_staticSpaceAsmTemplates), intent(inout), target :: rstaticAsmTemplates
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: nel

    select case (rphysics%cequation)
    
    ! ---------------------------------------------------------------
    ! Stokes/Navier Stokes.
    ! ---------------------------------------------------------------
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
    
      ! ---------------------------------
      ! Laplace matrix
      ! ---------------------------------
      
      ! Assemble the Laplace operator:
      call stdop_assembleLaplaceMatrix (rstaticAsmTemplates%rmatrixLaplace,&
          .true.,1.0_DP,rstaticAsmTemplates%rcubatureInfo)

      ! ---------------------------------
      ! B/D-matrices
      ! ---------------------------------

      ! Build the first pressure matrix B1.
      call stdop_assembleSimpleMatrix (rstaticAsmTemplates%rmatrixB1,&
          DER_FUNC,DER_DERIV_X,-1.0_DP,.true.,rstaticAsmTemplates%rcubatureInfoDiv)

      ! Build the second pressure matrix B2.
      call stdop_assembleSimpleMatrix (rstaticAsmTemplates%rmatrixB2,&
          DER_FUNC,DER_DERIV_Y,-1.0_DP,.true.,rstaticAsmTemplates%rcubatureInfoDiv)
      
      ! Set up the matrices D1 and D2 by transposing B1 and B2.
      call lsyssc_transposeMatrix (rstaticAsmTemplates%rmatrixB1,&
          rstaticAsmTemplates%rmatrixD1,LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rstaticAsmTemplates%rmatrixB2,&
          rstaticAsmTemplates%rmatrixD2,LSYSSC_TR_CONTENT)
          
      ! -----------------------------------------------------------------------
      ! Mass matrices. They are used in so many cases, is it better we always
      ! have them available.
      ! -----------------------------------------------------------------------

      ! Call the standard matrix setup routine to build the mass matrices.
      call stdop_assembleSimpleMatrix (rstaticAsmTemplates%rmatrixMass,&
          DER_FUNC,DER_FUNC,1.0_DP,.true.,rstaticAsmTemplates%rcubatureInfoMass)

      call stdop_assembleSimpleMatrix (rstaticAsmTemplates%rmatrixMassPressure,&
          DER_FUNC,DER_FUNC,1.0_DP,.true.,rstaticAsmTemplates%rcubatureInfoMassPressure)

      call stdop_assembleSimpleMatrix (rstaticAsmTemplates%rmatrixMassPressureExtStruc,&
          DER_FUNC,DER_FUNC,1.0_DP,.true.,rstaticAsmTemplates%rcubatureInfoMassPressure)
          
      ! Create a lumped mass matrix that represents the integral over the domain
      ! if being multiplied to a FEM function and summed up.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixMassPressure,&
          rstaticAsmTemplates%rmatrixMassPressureLumpInt,LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
          
      select case (rsettingsSpaceDiscr%ielementType)
      case (4)
        ! QP1 pressure. The matrix is created in a slightly different way
        ! to represent the integral over the domain
        call lsyssc_lumpMatrixScalar (rstaticAsmTemplates%rmatrixMassPressureLumpInt,&
            LSYSSC_LUMP_STD,.true.)
        nel = rstaticAsmTemplates%rmatrixMassPressureLumpInt%p_rspatialDiscrTest%p_rtriangulation%NEL
        call lsyssc_getbase_double (rstaticAsmTemplates%rmatrixMassPressureLumpInt,p_Ddata)
        call lalg_clearVector (p_Ddata(NEL+1:))
      case default
        call lsyssc_lumpMatrixScalar (rstaticAsmTemplates%rmatrixMassPressureLumpInt,&
            LSYSSC_LUMP_DIAG,.true.)
      end select
          
    ! ---------------------------------------------------------------
    ! Heat equation
    ! ---------------------------------------------------------------
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)

      ! ---------------------------------
      ! Laplace matrix
      ! ---------------------------------
      
      ! Assemble the Laplace operator:
      call stdop_assembleLaplaceMatrix (rstaticAsmTemplates%rmatrixLaplace,&
          .true.,1.0_DP,rstaticAsmTemplates%rcubatureInfo)

      ! -----------------------------------------------------------------------
      ! Mass matrices. They are used in so many cases, is it better we always
      ! have them available.
      ! -----------------------------------------------------------------------

      call stdop_assembleSimpleMatrix (rstaticAsmTemplates%rmatrixMass,&
          DER_FUNC,DER_FUNC,1.0_DP,.true.,rstaticAsmTemplates%rcubatureInfoMass)

      ! Create a lumped mass matrix that represents the integral over the domain
      ! if being multiplied to a FEM function and summed up.
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixMass,&
          rstaticAsmTemplates%rmatrixMassLumpInt,LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_lumpMatrixScalar (rstaticAsmTemplates%rmatrixMassLumpInt,&
          LSYSSC_LUMP_DIAG,.true.)

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine inmat_generateEOJmatrix (rphysics,rsettingsSpaceDiscr,rstaticAsmTemplates)
  
!<description>
  ! Calculates a matrix fort EOJ stabilisation.
!</description>

!<input>
  ! Definition of the underlying equation
  type(t_settings_physics), intent(in) :: rphysics

  ! Settings controlling the spatial discretisation (stabilisation parameters).
  ! This must coincide with the structure passed to inmat_initSpaceLevel.
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr
!</input>

!<inputoutput>
  ! A t_staticSpaceAsmTemplates structure. The static matrices in this structure are generated.
  type(t_staticSpaceAsmTemplates), intent(inout), target :: rstaticAsmTemplates
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_jumpStabilisation) :: rjumpStabil

    ! -------------------------------------------------------------------------
    ! EOJ matrix
    ! -------------------------------------------------------------------------
    
    select case (rphysics%cequation)
    
    ! ---------------------------------------------------------------
    ! Heat/Stokes/Navier Stokes.
    ! ---------------------------------------------------------------
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D,CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
    
      ! Create an empty matrix
      call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEM,&
          rstaticAsmTemplates%rmatrixEOJPrimal,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_clearMatrix (rstaticAsmTemplates%rmatrixEOJPrimal)

      ! Set up the jump stabilisation structure.
      ! There"s not much to do. Viscosity is set to 1.0 here.
      ! The operator is linear, so scaling the matrix by nu, one
      ! obtains the corresponding EOJ matrix.
      rjumpStabil%dnu = 1.0_DP
      
      ! Set stabilisation parameter
      rjumpStabil%dgamma = rsettingsSpaceDiscr%rstabilConvecPrimal%dupsam
      
      ! Matrix weight
      rjumpStabil%dtheta = 1.0_DP

      ! Call the jump stabilisation technique to stabilise that stuff.
      ! We can assemble the jump part any time as is it independent of any
      ! convective parts...
      call conv_jumpStabilisation2d (&
          rjumpStabil, CONV_MODMATRIX, rstaticAsmTemplates%rmatrixEOJPrimal)
          
      ! Subtract the boundary operator.
      !if (rsettingsSpaceDiscr%rstabilConvecPrimal%ceojStabilOnBoundary .eq. 0) then
      !  call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rstaticAsmTemplates%rmatrixEOJPrimal)
      !end if
      
      ! Primal and dual matrices identical?
      if ((rsettingsSpaceDiscr%rstabilConvecPrimal%dupsam .eq. &
          rsettingsSpaceDiscr%rstabilConvecDual%dupsam) &
          ! .and. &
          !(rsettingsSpaceDiscr%rstabilConvecPrimal%ceojStabilOnBoundary .eq. &
          !rsettingsSpaceDiscr%rstabilConvecDual%ceojStabilOnBoundary)
          ) then
      else
        ! Create another one for the dual space.
        rjumpStabil%dgamma = rsettingsSpaceDiscr%rstabilConvecDual%dupsam

        call lsyssc_duplicateMatrix (rstaticAsmTemplates%rmatrixTemplateFEM,&
            rstaticAsmTemplates%rmatrixEOJDual,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix (rstaticAsmTemplates%rmatrixEOJDual)

        call conv_jumpStabilisation2d (&
            rjumpStabil, CONV_MODMATRIX, rstaticAsmTemplates%rmatrixEOJDual)

        ! Subtract the boundary operator.
        !if (rsettingsSpaceDiscr%rstabilConvecDual%ceojStabilOnBoundary .eq. 0) then
        !  call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rstaticAsmTemplates%rmatrixEOJDual)
        !end if
      end if
      
    end select
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine inmat_initStaticAsmTemplHier(rhierarchy,rsettingsSpaceDiscr,&
      rfeHierarchyPrimal,rfeHierarchyDual,rfeHierarchyControl,rphysics)
  
!<description>
  ! Initialises the static matrices on all levels.
!</description>

!<input>
  ! A hierarchy of space levels for the primal space
  type(t_feHierarchy), intent(in) :: rfeHierarchyPrimal

  ! A hierarchy of space levels for the dual space
  type(t_feHierarchy), intent(in) :: rfeHierarchyDual

  ! A hierarchy of space levels for the control space
  type(t_feHierarchy), intent(in) :: rfeHierarchyControl

  ! Settings controlling the spatial discretisation (cubature)
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Physics of the problem
  type(t_settings_physics), intent(in) :: rphysics
!</input>

!<inputoutput>
  ! Level info hierarchy to initialise
  type(t_staticSpaceAsmHierarchy), intent(out) :: rhierarchy
!</inputoutput>

!</subroutine>

    integer :: ilevel

    ! Create the hierarchy.
    call astmpl_createSpaceAsmHier (rhierarchy,rfeHierarchyPrimal%nlevels)

    ! Calculate the structures
    do ilevel = 1,rfeHierarchyPrimal%nlevels
      call inmat_initSpaceLevel(&
          rhierarchy%p_RasmTemplList(ilevel),&
          rfeHierarchyPrimal%rmeshHierarchy%p_Rtriangulations(ilevel),rsettingsSpaceDiscr,&
          rphysics,&
          rfeHierarchyPrimal%p_rfeSpaces(ilevel)%p_rdiscretisation,&
          rfeHierarchyDual%p_rfeSpaces(ilevel)%p_rdiscretisation,&
          rfeHierarchyControl%p_rfeSpaces(ilevel)%p_rdiscretisation)
          
      call inmat_allocStaticMatrices (rhierarchy%p_RasmTemplList(ilevel),rphysics)
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine inmat_printTmplHierStatistic(rhierarchy)
  
!<description>
  ! Prints a statistic about the matrices in the hierarchy.
!</description>

!<input>
  ! Level info hierarchy
  type(t_staticSpaceAsmHierarchy), intent(in) :: rhierarchy
!</input>

!</subroutine>

    integer :: ilevel

    ! Print statistic data.
    call output_line ("Lv.       NA V-Mat       NA B-mat")
    call output_line ("---------------------------------")
    do ilevel = 1,rhierarchy%nlevels
      call output_line (&
          trim(sys_si(ilevel,3)) //" "// &
          trim(sys_si(rhierarchy%p_RasmTemplList(ilevel)%rmatrixTemplateFEM%NA,14)) //" "// &
          trim(sys_si(rhierarchy%p_RasmTemplList(ilevel)%rmatrixTemplateFEMPressure%NA,14)))
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine inmat_calcStaticLevelAsmHier(rhierarchy,rsettingsSpaceDiscr,rphysics,bprint)
  
!<description>
  ! Calculates the static matrices on all levels.
!</description>

!<input>
  ! OPTIONAL: Whether to print the current state of the assembly to
  ! the terminal. If set to TRUE, numbers "2 3 4..:" will be printed
  ! to the terminal for every level currently being processed.
  logical, intent(in), optional :: bprint
!</input>

!<input>
  ! Settings controlling the spatial discretisation (stabilisation parameters).
  ! This must coincide with the structure passed to inmat_initSpaceLevel.
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Physics of the problem
  type(t_settings_physics), intent(in) :: rphysics
!</input>

!<inputoutput>
  ! Level info hierarchy to initialise
  type(t_staticSpaceAsmHierarchy), intent(inout) :: rhierarchy
!</inputoutput>

!</subroutine>

    integer :: ilevel

    ! Calculate structures
    do ilevel = 1,rhierarchy%nlevels
      if (bprint) then
        ! Print current state.
        if (ilevel .eq. 1) then
          call output_line (trim(sys_siL(ilevel,10)),bnolinebreak=.true., &
              cdateTimeLogPolicy = OU_DTP_NONE)
        else
          call output_line (" "//trim(sys_siL(ilevel,10)),bnolinebreak=.true.,&
              cdateTimeLogPolicy = OU_DTP_NONE)
        end if
      end if

      call inmat_generateStaticMatrices (&
          rhierarchy%p_RasmTemplList(ilevel),rsettingsSpaceDiscr,rphysics)

    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine inmat_doneStaticAsmTemplHier(rhierarchy)
  
!<description>
  ! Cleans up the static matrices on all levels.
!</description>

!<inputoutput>
  ! Level info hierarchy to initialise
  type(t_staticSpaceAsmHierarchy), intent(inout) :: rhierarchy
!</inputoutput>

!</subroutine>

    integer :: ilevel

    ! Clean up the levels
    do ilevel = 1,rhierarchy%nlevels
      call inmat_releaseStaticMatrices(rhierarchy%p_RasmTemplList(ilevel))
      call inmat_doneSpaceLevel(rhierarchy%p_RasmTemplList(ilevel))
    end do
    
    ! Release the hierarchy
    call astmpl_releaseSpaceAsmHier (rhierarchy)
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine imnat_getL2PrjMatrix(rphysics,rasmTempl,rblockDiscr,rmassMatrix)
  
!<description>
  ! Creates a block mass matrix for the desired space which can be used
  ! for L2 projection
!</description>

!<input>
  ! Definition of the underlying equation
  type(t_settings_physics), intent(in) :: rphysics

  ! Assembly template structure
  type(t_staticSpaceAsmTemplates), intent(in) :: rasmTempl
  
  ! Discretisation structure of the space.
  type(t_blockDiscretisation), intent(in), target :: rblockDiscr
!</input>

!<inputoutput>
  ! Mass matrix.
  type(t_matrixBlock), intent(out) :: rmassMatrix
!</inputoutput>

!</subroutine>

    select case (rphysics%cequation)

    ! *************************************************************
    ! Stokes/Navier Stokes.
    ! *************************************************************
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      ! Create a block matrix
      call lsysbl_createMatBlockByDiscr (rblockDiscr,rmassMatrix)
      
      ! Plug in the mass matrices
      call lsyssc_duplicateMatrix (rasmTempl%rmatrixMass,rmassMatrix%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsyssc_duplicateMatrix (rasmTempl%rmatrixMass,rmassMatrix%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsyssc_duplicateMatrix (rasmTempl%rmatrixMassPressure,rmassMatrix%RmatrixBlock(3,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    ! *************************************************************
    ! Heat equation
    ! *************************************************************
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
      ! Create a block matrix
      call lsysbl_createMatBlockByDiscr (rblockDiscr,rmassMatrix)
      
      ! Plug in the mass matrices
      call lsyssc_duplicateMatrix (rasmTempl%rmatrixMass,rmassMatrix%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
    end select
    
  end subroutine

end module
