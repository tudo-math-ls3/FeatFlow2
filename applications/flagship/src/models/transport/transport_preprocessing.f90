!##############################################################################
!# ****************************************************************************
!# <Name> transport_preprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all preprocessing routines which are required to
!# solve scalar conservation laws in arbitrary spatial dimensions.
!#
!# The following routines are available:
!#
!# 1.) transp_initSolvers
!#     -> Initialises the solve structures from the parameter list.
!#
!# 2.) transp_initProblemDescriptor
!#     -> Initialises the abstract problem descriptor based on the
!#        parameter settings given by the parameter list.
!#
!# 3.) transp_initProblemLevel
!#     -> Initialises the individual problem level based on the
!#        parameter settings given by the parameter list.
!#        This routine is called repeatedly by the global
!#        initialisation routine transp_initAllProblemLevels.
!#
!# 4.) transp_initAllProblemLevels
!#     -> Initialises ALL problem levels attached to the global
!#        problem structure based on the parameter settings
!#        given by the parameter list.
!#
!# 5.) transp_initSolution
!#     -> Initialises the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# 6.) transp_initRHS
!#     -> Initialises the right-hand side vector based on the
!#        parameter settings given by the application descriptor
!#
!# 7.) transp_initTargetFunc
!#     -> Initialises the target functional for the dual problem
!# </purpose>
!##############################################################################

module transport_preprocessing

#include "../../flagship.h"

!$use omp_lib
  use afcstabbase
  use afcstabscalar
  use afcstabscalarfct
  use basicgeometry
  use bilinearformevaluation
  use boundary
  use boundarycondaux
  use collection
  use cubature
  use derivatives
  use dofmapping
  use element
  use fparser
  use fsystem
  use genoutput
  use groupfembase
  use linearalgebra
  use linearformevaluation
  use lineariser
  use linearsystemblock
  use linearsystemscalar
  use meshmodification
  use paramlist
  use pprocsolution
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use stdoperators
  use storage
  use timestep
  use timestepaux
  use triangulation

  ! Modules from transport model
  use transport_basic
  use transport_callback

  implicit none

  private

  public :: transp_initAllProblemLevels
  public :: transp_initProblemDescriptor
  public :: transp_initProblemLevel
  public :: transp_initRHS
  public :: transp_initSolution
  public :: transp_initSolvers
  public :: transp_initTargetFunc

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_initSolvers(rparlist, ssectionName,&
      rtimestep, rsolver)

!<description>
    ! This subroutine initialises the time-stepping structure
    ! and the top-level solver structure using the
    ! parameter settings defined in the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<output>
    ! time-stepping structure
    type(t_timestep), intent(out) :: rtimestep

    ! solver struchture
    type(t_solver), intent(out) :: rsolver
!</output>
!</subroutine>

    ! section name for the top-level solver
    character(LEN=SYS_STRLEN) :: ssolverName

    ! section name for time-stepping scheme
    character(LEN=SYS_STRLEN) :: stimestepName

    ! local variables
    integer :: nlmin, nlmax


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'timestep', stimestepName)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'solver',   ssolverName)

    ! Initialise time-stepping
    call tstep_createTimestep(rparlist, stimestepName, rtimestep)

    ! Initialise solver structure
    call solver_createSolver(rparlist, ssolverName, rsolver)
    if (rsolver%csolverType .eq. SV_FMG) then
      nlmin = rsolver%p_rsolverMultigrid%nlmin
      nlmax = rsolver%p_rsolverMultigrid%nlmax
      call solver_adjustHierarchy(rsolver, nlmin, nlmax)
    else
      call solver_adjustHierarchy(rsolver)
    end if
    call solver_updateStructure(rsolver)

  end subroutine transp_initSolvers

  !*****************************************************************************

!<subroutine>
  
  subroutine transp_initProblemDescriptor(rparlist, ssectionName,&
      nlmin, nlmax, rproblemDescriptor)

!<description>
    ! This subroutine initialises the abstract problem descriptor
    ! using the parameters settings defined in the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! minimum/maximum problem level
    integer, intent(in) :: nlmin, nlmax
!</input>

!<output>
    ! problem descriptor
    type(t_problemDescriptor), intent(out) :: rproblemDescriptor
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ddisturbmesh
    integer :: discretisation,velocityfield,dofCoords
    integer :: massAFC,templateGFEM
    integer :: convectionAFC,convectionGFEM
    integer :: diffusionAFC,diffusionGFEM
    integer :: primalBdrGFEM,dualBdrGFEM
    integer :: templateMatrix
    integer :: systemMatrix
    integer :: jacobianMatrix
    integer :: transportMatrix
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: coeffMatrix_S
    integer :: iconvToTria
    integer :: ncubatureInfo
        
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iconvtotria', iconvToTria, TRI_CONVERT_NONE)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dofCoords', dofCoords, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'velocityfield', velocityfield, 0)
    call parlst_getvalue_double(rparlist,&
        ssectionName, 'ddisturbmesh', ddisturbmesh, 0.0_DP)

    ncubatureInfo = parlst_querysubstrings(rparlist,&
        ssectionName, 'CubatureInfo')

    ! Get global positions of matrices
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateMatrix', templateMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemMatrix', systemMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianMatrix', jacobianMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'transportMatrix', transportMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'consistentMassMatrix', consistentMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedMassMatrix', lumpedMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CX', coeffMatrix_CX, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CY', coeffMatrix_CY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CZ', coeffMatrix_CZ, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_S', coeffMatrix_S, 0)

    ! Default is no stabilization
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC, 0)

    ! By default the same identifier is used for the group finite
    ! element formulation and the stabilization structure
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateGFEM', templateGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionGFEM', convectionGFEM, convectionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionGFEM', diffusionGFEM, diffusionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'primalbdrGFEM', primalBdrGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dualbdrGFEM', dualBdrGFEM, 0)

    ! Consistency check
    if (templateGFEM .eq. 0) then
      convectionGFEM = 0
      diffusionGFEM  = 0
      primalBdrGFEM  = 0
      dualBdrGFEM    = 0
    end if

    ! Set additional problem descriptor
    rproblemDescriptor%nblockdiscretisation = max(0, discretisation)
    rproblemDescriptor%ncubatureInfo   = max(0, ncubatureInfo)
    rproblemDescriptor%nafcstab        = max(0, massAFC,&
                                                convectionAFC,&
                                                diffusionAFC)
    rproblemDescriptor%ngroupfemBlock  = max(0, templateGFEM,&
                                                convectionGFEM,&
                                                diffusionGFEM,&
                                                primalBdrGFEM,&
                                                dualBdrGFEM)
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = max(0, templateMatrix,&
                                                systemMatrix,&
                                                jacobianMatrix,&
                                                transportMatrix,&
                                                consistentMassMatrix,&
                                                lumpedMassMatrix,&
                                                coeffMatrix_CX,&
                                                coeffMatrix_CY,&
                                                coeffMatrix_CZ,&
                                                coeffMatrix_S)
    rproblemDescriptor%nmatrixBlock    = 0
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = max(0, velocityfield, dofCoords)
    rproblemDescriptor%iconvStrategy   = iconvToTria
    rproblemDescriptor%ddisturbmesh    = ddisturbmesh

  end subroutine transp_initProblemDescriptor

  !*****************************************************************************

!<subroutine>

  subroutine transp_initProblemLevel(rparlist, ssectionName,&
      rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

!<description>
    ! This subroutine initielises the individual problem level. It
    ! generates the discretisation, the template matrix and the
    ! coefficient matrices as duplicates of the template matrix.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: boundary condition for primal problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondPrimal

    ! OPTIONAL: boundary condition for dual problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondDual
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</output>
!</subroutine>

    ! local variables
    integer :: templateMatrix
    integer :: systemMatrix
    integer :: jacobianMatrix
    integer :: transportMatrix
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: coeffMatrix_S
    integer :: massAFC,templateGFEM
    integer :: convectionAFC,convectionGFEM
    integer :: diffusionAFC,diffusionGFEM
    integer :: discretisation
    integer :: dofCoords
    integer :: ijacobianFormat
    integer :: imatrixFormat
    integer :: primalbdrGFEM,dualbdrGFEM

    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_triangulation) , pointer :: p_rtriangulation
    type(t_groupFEMSet), pointer :: p_rgroupFEMSet
    type(t_fparser), pointer :: p_rfparser
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_matrixScalar) :: rmatrixBdrSx,rmatrixBdrSy,rmatrixBdrSz
    integer(I32), dimension(:), allocatable :: Celement
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    real(DP) :: dmeshdisturb
    integer :: nlevel,ncubatureInfo
    integer :: i,ibdc,isegment,nmatrices,nsubstrings,neq
    character(len=SYS_STRLEN) :: selemName,smass,sconvection,sdiffusion
    character(len=SYS_STRLEN) :: scubatureInfo,scubType

    ! Retrieve application specific parameters from the parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templatematrix', templateMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemmatrix', systemMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianmatrix', jacobianMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'transportmatrix', transportMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_s', coeffMatrix_S, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cx', coeffMatrix_CX, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cy', coeffMatrix_CY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cz', coeffMatrix_CZ, 0)

    ! Default is no stabilization
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC, 0)

    ! By default the same identifier is used for the group finite
    ! element formulation and the stabilization structure
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateGFEM', templateGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionGFEM', convectionGFEM, convectionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionGFEM', diffusionGFEM, diffusionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'primalBdrGFEM', primalBdrGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dualBdrGFEM', dualBdrGFEM, 0)

    ! Default no summed cubature
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dofCoords', dofCoords, 0)

    ! Default is empty section, i.e. no configuration
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'mass', smass, '')
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'diffusion', sdiffusion, '')
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'convection', sconvection, '')

    ! Set pointers to triangulation and boundary structure
    p_rtriangulation => rproblemLevel%rtriangulation

    ! Disturb the mesh?
    call parlst_getvalue_double(rparlist,&
        ssectionName, 'dmeshdisturb', dmeshdisturb, 0.0_DP)
    if (dmeshdisturb .gt. 0.0_DP)&
        call meshmod_disturbMesh(p_rtriangulation, dmeshdisturb)
   
    !---------------------------------------------------------------------------
    ! Create discretisation structure
    !---------------------------------------------------------------------------
    if (discretisation > 0) then

      ! Initialise the block discretisation structure
      p_rdiscretisation => rproblemLevel%RblockDiscretisation(discretisation)
      if (p_rdiscretisation%ndimension .eq. 0) then
        call spdiscr_initBlockDiscr(p_rdiscretisation, 1, p_rtriangulation,&
            rproblemLevel%p_rproblem%rboundary)
      end if

      ! Allocate temporal memory
      nsubstrings = max(1, parlst_querysubstrings(rparlist,&
                           ssectionName, 'celement'))
      allocate(Celement(nsubstrings))

      ! Get IDs of element types
      do i = 1, nsubstrings
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'celement', selemName, isubstring=i)
        Celement(i) = elem_igetID(selemName)
      end do
      
      ! Get spatial dimension
      select case(p_rdiscretisation%ndimension)
      case (NDIM1D)
        call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1),&
            Celement(1), p_rtriangulation, rproblemLevel%p_rproblem%rboundary)

      case (NDIM2D)
        if (size(Celement) .eq. 1) then
          call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1),&
              Celement(1), p_rtriangulation, rproblemLevel%p_rproblem%rboundary)
        else
          call spdiscr_initDiscr_triquad(&
              p_rdiscretisation%RspatialDiscr(1), Celement(1), Celement(2),&
              p_rtriangulation, rproblemLevel%p_rproblem%rboundary)
        end if

      case (NDIM3D)
        call spdiscr_initDiscr_simple(p_rdiscretisation%RspatialDiscr(1),&
            Celement(1), p_rtriangulation, rproblemLevel%p_rproblem%rboundary)
      
      case default
        call output_line('Invalid number of spatial dimensions',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_initProblemLevel')
        call sys_halt()
      end select
      
      ! Deallocate temporal memory
      deallocate(Celement)

      !-------------------------------------------------------------------------
      ! Create cubature info structures
      !-------------------------------------------------------------------------
      ncubatureInfo = parlst_querysubstrings(rparlist,&
          ssectionName, 'CubatureInfo')
      
      do i = 1, ncubatureInfo
        ! Get section name of cubatureinfo structire
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'CubatureInfo', scubatureinfo, isubstring=i)
        
        ! Read cubatureinfo from parameter file
        call parlst_getvalue_string(rparlist,&
            scubatureInfo, 'scubType', scubType)
        call parlst_getvalue_int(rparlist,&
            scubatureInfo, 'nlevels', nlevel)
        
        ! Create cubatureinfo structure
        call spdiscr_createDefCubStructure(&  
            p_rdiscretisation%RspatialDiscr(1),&
            rproblemLevel%RcubatureInfo(i),&
            cub_igetID(scubType), nlevel)
      end do
      
      !-------------------------------------------------------------------------
      ! Calculate coordinates of the global DOF`s
      if (dofCoords > 0) then
        ! Check if block vector has been initialised
        if (rproblemLevel%RvectorBlock(dofCoords)%nblocks .eq. 0) then
          call lsysbl_createVectorBlock(p_rdiscretisation,&
              p_rdiscretisation%ndimension,&
              rproblemLevel%RvectorBlock(dofCoords), .false.)
        else
          neq = dof_igetNDofGlobBlock(p_rdiscretisation)
          if (rproblemLevel%RvectorBlock(dofCoords)%NEQ .ne. neq) then
            call lsysbl_resizeVectorBlock(p_rdiscretisation,&
                rproblemLevel%RvectorBlock(dofCoords), .false.)
          end if
        end if
        ! Calculate coordinates
        call lin_calcDofCoordsBlock(p_rdiscretisation,&
            rproblemLevel%RvectorBlock(dofCoords))
      end if

      !-------------------------------------------------------------------------
      ! Create finite element matrices
      !-------------------------------------------------------------------------
      
      ! If the template matrix has no structure data then generate the
      ! finite element matrix sparsity structure based on the spatial
      ! descretisation and store it as the template matrix. Otherwise we
      ! assume that the template matrix has been generated externally.
      if (.not.lsyssc_hasMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix))) then
        call parlst_getvalue_int(rparlist, ssectionName, 'imatrixFormat', imatrixFormat)
        call bilf_createMatrixStructure(p_rdiscretisation%RspatialDiscr(1),&
            imatrixFormat, rproblemLevel%RmatrixScalar(templateMatrix))
      end if
      
      !-------------------------------------------------------------------------
      ! Create system matrix as duplicate of the template matrix
      if (systemMatrix > 0)&
          call initMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix),&
                                   rproblemLevel%RmatrixScalar(systemMatrix))

      !-------------------------------------------------------------------------
      ! Create transport matrix as duplicate of the template matrix
      if (transportMatrix > 0)&
          call initMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix),&
                                   rproblemLevel%RmatrixScalar(transportMatrix))

      !-------------------------------------------------------------------------
      ! Create Jacobian matrix. This is a little bit tricky. If the
      ! Jacobian matrix has the same sparsity pattern as the template
      ! matrix, we can just create the Jacobian matrix as a duplicate of
      ! the template matrix. If the Jacobian matrix has an extended
      ! sparsity pattern we must create it by using the template matrix
      if (jacobianMatrix > 0) then
        
        ! What format do we have for the Jacobian matrix?
        call parlst_getvalue_int(rparlist, ssectionName,&
            'ijacobianFormat', ijacobianFormat)
        
        if (lsyssc_hasMatrixStructure(rproblemLevel%RmatrixScalar(jacobianMatrix))) then
          if (ijacobianFormat .eq. 0) then
            call lsyssc_resizeMatrix(&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(templateMatrix),&
                .false., .false., bforce=.true.)
          else
            call afcstab_genExtSparsity(&
                rproblemLevel%RmatrixScalar(templateMatrix),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if
        else
          if (ijacobianFormat .eq. 0) then
            call lsyssc_duplicateMatrix(&
                rproblemLevel%RmatrixScalar(templateMatrix),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
          else
            call afcstab_genExtSparsity(&
                rproblemLevel%RmatrixScalar(templateMatrix),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if
        end if
      end if

      !-------------------------------------------------------------------------
      ! Create consistent mass matrix as duplicate of the template matrix
      if (consistentMassMatrix > 0) then
        call initMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix),&
                                 rproblemLevel%RmatrixScalar(consistentMassMatrix))
        ! Do we have a precomputed cubatureinfo structure?
        i = -1; nsubstrings =&
            parlst_querysubstrings(rparlist,&
            ssectionName, 'consistentmassmatrix')
        if (nsubstrings .eq. 1) then
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'consistentmassmatrix', scubatureInfo, isubstring=1)
          i = parlst_findvalue(rparlist,&
              ssectionName, 'CubatureInfo', trim(adjustl(scubatureInfo)))
        end if

        if (i .eq. -1) then
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              DER_FUNC, DER_FUNC)
        else
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              DER_FUNC, DER_FUNC, rcubatureInfo=rproblemLevel%RcubatureInfo(i))
        end if
                
        ! Create lumped mass matrix
        if (lumpedMassMatrix > 0) then
          call lsyssc_duplicateMatrix(&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
          call lsyssc_lumpMatrix(&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
        end if
      elseif (lumpedMassMatrix > 0) then
        ! Create lumped mass matrix
        call lsyssc_duplicateMatrix(&
            rproblemLevel%RmatrixScalar(templateMatrix),&
            rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

        ! Do we have a precomputed cubatureinfo structure?
        i = -1; nsubstrings =&
            parlst_querysubstrings(rparlist,&
            ssectionName, 'lumpedmassmatrix')
        if (nsubstrings .eq. 1) then
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'lumpedmassmatrix', scubatureInfo, isubstring=1)
          i = parlst_findvalue(rparlist,&
              ssectionName, 'CubatureInfo', trim(adjustl(scubatureInfo)))
        end if

        if (i .eq. -1) then
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              DER_FUNC, DER_FUNC)
        else
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              DER_FUNC, DER_FUNC, rcubatureInfo=rproblemLevel%RcubatureInfo(i))
        end if

        call lsyssc_lumpMatrix(&
            rproblemLevel%RmatrixScalar(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
      end if
      
      !-------------------------------------------------------------------------
      ! Create diffusion matrix as duplicate of the template matrix
      if (coeffMatrix_S > 0) then
        call initMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix),&
                                 rproblemLevel%RmatrixScalar(coeffMatrix_S))
        
        ! Get function parser from collection
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)

        ! Do we have a precomputed cubatureinfo structure?
        i = -1; nsubstrings =&
            parlst_querysubstrings(rparlist,&
            ssectionName, 'consistentmassmatrix')
        if (nsubstrings .eq. 1) then
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'coeffMatrix_S', scubatureInfo, isubstring=1)
          i = parlst_findvalue(rparlist,&
              ssectionName, 'CubatureInfo', trim(adjustl(scubatureInfo)))
        end if
        
        if (i .eq. -1) then
          ! Assemble diffusion matrix without precomputed cubatureinfo
          select case(p_rtriangulation%ndim)
          case (NDIM1D)
            call initDiffusionMatrix1D(p_rfparser,&
                rproblemLevel%RmatrixScalar(coeffMatrix_S))
          case (NDIM2D)
            call initDiffusionMatrix2D(p_rfparser,&
                rproblemLevel%RmatrixScalar(coeffMatrix_S))
          case (NDIM3D)
            call initDiffusionMatrix3D(p_rfparser,&
                rproblemLevel%RmatrixScalar(coeffMatrix_S))
          case default
            call lsyssc_releaseMatrix(rproblemLevel%RmatrixScalar(coeffMatrix_S))
          end select
        else
          ! Assemble diffusion matrix with precomputed cubatureinfo
          select case(p_rtriangulation%ndim)
          case (NDIM1D)
            call initDiffusionMatrix1D(p_rfparser,&
                rproblemLevel%RmatrixScalar(coeffMatrix_S),&
                rproblemLevel%RcubatureInfo(i))
          case (NDIM2D)
            call initDiffusionMatrix2D(p_rfparser,&
                rproblemLevel%RmatrixScalar(coeffMatrix_S),&
                rproblemLevel%RcubatureInfo(i))
          case (NDIM3D)
            call initDiffusionMatrix3D(p_rfparser,&
                rproblemLevel%RmatrixScalar(coeffMatrix_S),&
                rproblemLevel%RcubatureInfo(i))
          case default
            call lsyssc_releaseMatrix(rproblemLevel%RmatrixScalar(coeffMatrix_S))
          end select
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dx) as duplicate of the
      ! template matrix
      if (coeffMatrix_CX > 0) then
        call initMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix),&
                                 rproblemLevel%RmatrixScalar(coeffMatrix_CX))
        ! Do we have a precomputed cubatureinfo structure?
        i = -1; nsubstrings =&
            parlst_querysubstrings(rparlist,&
            ssectionName, 'coeffMatrix_CX')
        if (nsubstrings .eq. 1) then
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'coeffMatrix_CX', scubatureInfo, isubstring=1)
          i = parlst_findvalue(rparlist,&
              ssectionName, 'CubatureInfo', trim(adjustl(scubatureInfo)))
        end if

        if (i .eq. -1) then
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CX),&
              DER_DERIV3D_X, DER_FUNC)
        else
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CX),&
              DER_DERIV3D_X, DER_FUNC, rcubatureInfo=rproblemLevel%RcubatureInfo(i))
        end if
      end if

      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dy) as duplicate of the
      ! template matrix
      if (coeffMatrix_CY > 0) then
        call initMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix),&
                                 rproblemLevel%RmatrixScalar(coeffMatrix_CY))
        ! Do we have a precomputed cubatureinfo structure?
        i = -1; nsubstrings =&
            parlst_querysubstrings(rparlist,&
            ssectionName, 'coeffMatrix_CY')
        if (nsubstrings .eq. 1) then
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'coeffMatrix_CY', scubatureInfo, isubstring=1)
          i = parlst_findvalue(rparlist,&
              ssectionName, 'CubatureInfo', trim(adjustl(scubatureInfo)))
        end if

        if (i .eq. -1) then
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CY),&
              DER_DERIV3D_Y, DER_FUNC)
        else
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CY),&
              DER_DERIV3D_Y, DER_FUNC, rcubatureInfo=rproblemLevel%RcubatureInfo(i))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dz) as duplicate of the
      ! template matrix
      if (coeffMatrix_CZ > 0) then
        call initMatrixStructure(rproblemLevel%RmatrixScalar(templateMatrix),&
                                 rproblemLevel%RmatrixScalar(coeffMatrix_CZ))
        ! Do we have a precomputed cubatureinfo structure?
        i = -1; nsubstrings =&
            parlst_querysubstrings(rparlist,&
            ssectionName, 'coeffMatrix_CZ')
        if (nsubstrings .eq. 1) then
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'coeffMatrix_CZ', scubatureInfo, isubstring=1)
          i = parlst_findvalue(rparlist,&
              ssectionName, 'CubatureInfo', trim(adjustl(scubatureInfo)))
        end if

        if (i .eq. -1) then
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CZ),&
              DER_DERIV3D_Z, DER_FUNC)
        else
          call stdop_assembleSimpleMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CZ),&
              DER_DERIV3D_Z, DER_FUNC, rcubatureInfo=rproblemLevel%RcubatureInfo(i))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Create group finite element structures and AFC-stabilisations
      !-------------------------------------------------------------------------
      
      ! Initialise/resize template group finite element structure and
      ! generate the edge structure derived from the template matrix
      if (templateGFEM > 0) then
        ! Check if structure has been initialised
        if (rproblemLevel%RgroupFEMBlock(templateGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(templateGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialise first group finite element set for edge-based assembly
          call gfem_initGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(templateMatrix), 0, 0, 0, GFEM_EDGEBASED)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(templateMatrix))
        end if
        
        ! Generate diagonal and edge structure derived from template matrix
        call gfem_genDiagList(rproblemLevel%RmatrixScalar(templateMatrix),&
            p_rgroupFEMSet)
        call gfem_genEdgeList(rproblemLevel%RmatrixScalar(templateMatrix),&
            p_rgroupFEMSet)
      else
        convectionGFEM = 0
        diffusionGFEM  = 0
        primalBdrGFEM  = 0
        dualBdrGFEM    = 0
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/resize group finite element structure as duplicate of
      ! the template group finite element structure and fill it with the
      ! precomputed matrix coefficients for the convective term
      if (convectionGFEM > 0) then
        ! Check if structure has been initialised
        if (rproblemLevel%RgroupFEMBlock(convectionGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(convectionGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialise first group finite element set for edge-based
          ! assembly as aduplicate of the template structure
          call gfem_duplicateGroupFEMSet(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              p_rgroupFEMSet, GFEM_DUP_STRUCTURE, .false.)
          
          ! Compute number of matrices to be copied
          nmatrices = 0
          if (coeffMatrix_CX > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CY > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CZ > 0) nmatrices = nmatrices+1
          
          ! Allocate memory for matrix entries
          call gfem_allocCoeffs(p_rgroupFEMSet, nmatrices, 0, nmatrices)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(templateMatrix))
        end if
        
        ! Duplicate edge-based structure from template
        call gfem_duplicateGroupFEMSet(&
            rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
            p_rgroupFEMSet, GFEM_DUP_DIAGLIST+GFEM_DUP_EDGELIST, .true.)
        
        ! Copy constant coefficient matrices to group finite element set
        nmatrices = 0
        if (coeffMatrix_CX > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(coeffMatrix_CX), nmatrices)
        end if
        if (coeffMatrix_CY > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(coeffMatrix_CY), nmatrices)
        end if
        if (coeffMatrix_CZ > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(coeffMatrix_CZ), nmatrices)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the convective
      ! term by duplicating parts of the template group finite element set
      if (convectionAFC > 0) then
        if (rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          call afcstab_initFromParameterlist(rparlist, sconvection,&
              rproblemLevel%Rafcstab(convectionAFC))
          call afcsc_initStabilisation(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rafcstab(convectionAFC), p_rdiscretisation)
        else
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/resize group finite element structure as duplicate of
      ! the template group finite element structure and fill it with the
      ! precomputed matrix coefficients for the diffusife term
      if (diffusionGFEM > 0) then
        ! Check if structure has been initialised
        if (rproblemLevel%RgroupFEMBlock(diffusionGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(diffusionGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(diffusionGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialise first group finite element set for edge-based
          ! assembly as aduplicate of the template structure
          call gfem_duplicateGroupFEMSet(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              p_rgroupFEMSet, GFEM_DUP_STRUCTURE, .false.)
          
          ! Compute number of matrices to be copied
          nmatrices = 0
          if (coeffMatrix_S > 0) nmatrices = nmatrices+1
          
          ! Allocate memory for matrix entries
          call gfem_allocCoeffs(p_rgroupFEMSet, nmatrices, 0, nmatrices)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(templateMatrix))
        end if
        
        ! Duplicate edge-based structure from template
        call gfem_duplicateGroupFEMSet(&
            rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
            p_rgroupFEMSet, GFEM_DUP_DIAGLIST+GFEM_DUP_EDGELIST, .true.)
        
        ! Copy constant coefficient matrices to group finite element set
        nmatrices = 0
        if (coeffMatrix_S > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%RmatrixScalar(coeffMatrix_S), nmatrices)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the diffusive
      ! term by duplicating parts of the template group finite element set
      if (diffusionAFC > 0) then
        if (rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          call afcstab_initFromParameterlist(rparlist, sdiffusion,&
              rproblemLevel%Rafcstab(diffusionAFC))
          rproblemLevel%Rafcstab(diffusionAFC)%bisSymmetricOperator = .true.
          call afcsc_initStabilisation(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rafcstab(diffusionAFC), p_rdiscretisation)
        else
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(diffusionAFC),&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the mass matrix by
      ! duplicating parts of the template group finite element set
      if (massAFC > 0) then
        if (rproblemLevel%Rafcstab(massAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          call afcstab_initFromParameterlist(rparlist, smass,&
              rproblemLevel%Rafcstab(massAFC))
          call afcsc_initStabilisation(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rafcstab(massAFC), p_rdiscretisation)
        else
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(massAFC),&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1))
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Create group finite element structures on the boundary
      !-------------------------------------------------------------------------
      
      if (((primalBdrGFEM > 0) .and. present(rbdrCondPrimal)) .or.&
          ((dualBdrGFEM   > 0) .and. present(rbdrCondDual))) then
        
        if (coeffMatrix_CX > 0) then
          call lsyssc_duplicateMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CX), rmatrixBdrSx,&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_COPYOVERWRITE)
          call lsyssc_createMatrixSymmPart(rmatrixBdrSx, 1.0_DP)
        end if
        
        if (coeffMatrix_CY > 0) then
          call lsyssc_duplicateMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CY), rmatrixBdrSy,&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_COPYOVERWRITE)
          call lsyssc_createMatrixSymmPart(rmatrixBdrSy, 1.0_DP)
        end if

        if (coeffMatrix_CZ > 0) then
          call lsyssc_duplicateMatrix(&
              rproblemLevel%RmatrixScalar(coeffMatrix_CZ), rmatrixBdrSz,&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_COPYOVERWRITE)
          call lsyssc_createMatrixSymmPart(rmatrixBdrSz, 1.0_DP)
        end if
      end if

      !-------------------------------------------------------------------------
      ! Primal boundary condition
      !-------------------------------------------------------------------------
      
      if ((primalBdrGFEM > 0) .and. present(rbdrCondPrimal)) then
        
        ! Check if structure has been initialised
        if (rproblemLevel%RgroupFEMBlock(primalBdrGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
            bdrc_getNumberOfRegions(rbdrCondPrimal))
        
        ! Compute number of matrices to be copied
        nmatrices = 0
        if (coeffMatrix_CX > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CY > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CZ > 0) nmatrices = nmatrices+1
        
        ! Set pointer
        call storage_getbase_int(rbdrCondPrimal%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
        
        ! Loop over all boundary components
        do ibdc = 1, rbdrCondPrimal%iboundarycount
      
          ! Loop over all boundary segments
          do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
            
            ! Create boundary region
            call bdrc_createRegion(rbdrCondPrimal, ibdc,&
                isegment-p_IbdrCondCpIdx(ibdc)+1, rboundaryRegion)
            
            ! Set pointer to group finite element set
            p_rgroupFEMSet =>&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM)%RgroupFEMBlock(isegment)
            
            ! Initialise group finite element
            call initGroupFEMSetBoundary(rboundaryRegion,&
                rproblemLevel%RmatrixScalar(templateMatrix), nmatrices, p_rgroupFEMSet)
            
            ! Copy constant coefficient matrices to group finite element set
            if (coeffMatrix_CX > 0) then
              call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixBdrSx, 1)
            end if
            if (coeffMatrix_CY > 0) then
              call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixBdrSy, 2)
            end if
            if (coeffMatrix_CZ > 0) then
              call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixBdrSz, 3)
            end if
            
          end do
        end do
      end if
      
      !-------------------------------------------------------------------------
      ! Dual boundary condition
      !-------------------------------------------------------------------------
      
      if ((dualBdrGFEM > 0) .and. present(rbdrCondDual)) then
        
        ! Check if structure has been initialised
        if (rproblemLevel%RgroupFEMBlock(dualBdrGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
            bdrc_getNumberOfRegions(rbdrCondDual))
        
        ! Compute number of matrices to be copied
        nmatrices = 0
        if (coeffMatrix_CX > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CY > 0) nmatrices = nmatrices+1
        if (coeffMatrix_CZ > 0) nmatrices = nmatrices+1
        
        ! Set pointer
        call storage_getbase_int(rbdrCondDual%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
        
        ! Loop over all boundary components
        do ibdc = 1, rbdrCondDual%iboundarycount
          
          ! Loop over all boundary segments
          do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
            
            ! Create boundary region
            call bdrc_createRegion(rbdrCondDual, ibdc,&
                isegment-p_IbdrCondCpIdx(ibdc)+1, rboundaryRegion)
            
            ! Set pointer to group finite element set
            p_rgroupFEMSet =>&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM)%RgroupFEMBlock(isegment)
            
            ! Initialise group finite element
            call initGroupFEMSetBoundary(rboundaryRegion,&
                rproblemLevel%RmatrixScalar(templateMatrix), nmatrices, p_rgroupFEMSet)
            
            ! Copy constant coefficient matrices to group finite element set
            if (coeffMatrix_CX > 0) then
              call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixBdrSx, 1)
            end if
            if (coeffMatrix_CY > 0) then
              call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixBdrSy, 2)
            end if
            if (coeffMatrix_CZ > 0) then
              call gfem_initCoeffsFromMatrix(p_rgroupFEMSet, rmatrixBdrSz, 3)
            end if
          end do
        end do
      end if

      ! Clean auxiliary matrices
      call lsyssc_releaseMatrix(rmatrixBdrSx)
      call lsyssc_releaseMatrix(rmatrixBdrSy)
      call lsyssc_releaseMatrix(rmatrixBdrSz)
      
    end if   ! discretisation > 0

    !---------------------------------------------------------------------------
    ! Set update notifiers for the discrete transport operator and the
    ! preconditioned in the problem level structure
    !---------------------------------------------------------------------------
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_TROPER_UPDATE)
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_PRECOND_UPDATE)

  contains

    !**************************************************************
    ! Initialise the matrix structure by duplicating the template matrix

    subroutine initMatrixStructure(rmatrixTemplate, rmatrix)

      type(t_matrixScalar), intent(in) :: rmatrixTemplate
      type(t_matrixScalar), intent(inout) :: rmatrix

      if (lsyssc_isMatrixStructureShared(rmatrix, rmatrixTemplate)) then
        call lsyssc_resizeMatrix(rmatrix, rmatrixTemplate,&
            .false., .false., bforce=.true.)
      else
        call lsyssc_duplicateMatrix(rmatrixTemplate, rmatrix,&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      end if

    end subroutine initMatrixStructure
    
    !**************************************************************
    ! Initialise the diffusion matrix in 1D

    subroutine initDiffusionMatrix1D(rfparser, rmatrix, rcubatureinfo)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix
      type(t_scalarCubatureInfo), intent(in), optional :: rcubatureInfo

      ! local variables
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC,&
            DIFFUSION_ANISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist, ssectionName, 'sdiffusionname',&
                                    sdiffusionName, isubString=1)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, sdiffusionName, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha, rcubatureInfo)

      case default
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix1D

    !**************************************************************
    ! Initialise the diffusion matrix in 2D

    subroutine initDiffusionMatrix2D(rfparser, rmatrix, rcubatureInfo)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix
      type(t_scalarCubatureInfo), intent(in), optional :: rcubatureInfo

      ! local variables
      type(t_bilinearform) :: rform
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: i,idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, sdiffusionName, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha, rcubatureInfo)

      case (DIFFUSION_ANISOTROPIC)

        do i = 1, 4
          ! Retrieve name/number of expression describing the diffusion coefficient
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'sdiffusionname', sdiffusionName, isubString=i)

          ! Evaluate the constant coefficient from the function parser
          call fparser_evalFunction(rfparser, sdiffusionName, Dunity,&
              rform%Dcoefficients(i))
        end do

        ! We have constant coefficients
        rform%ballCoeffConstant   = .true.
        rform%Dcoefficients(1:4)  = -rform%Dcoefficients(1:4)

        ! Initialise the bilinear form
        rform%itermCount = 4
        rform%Idescriptors(1,1) = DER_DERIV2D_X
        rform%Idescriptors(2,1) = DER_DERIV2D_X

        rform%Idescriptors(1,2) = DER_DERIV2D_X
        rform%Idescriptors(2,2) = DER_DERIV2D_Y

        rform%Idescriptors(1,3) = DER_DERIV2D_Y
        rform%Idescriptors(2,3) = DER_DERIV2D_X

        rform%Idescriptors(1,4) = DER_DERIV2D_Y
        rform%Idescriptors(2,4) = DER_DERIV2D_Y

        ! Assemble the anisotropic diffusion matrix
        call bilf_buildMatrixScalar(rform, .true., rmatrix, rcubatureInfo)

      case default
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix2D

    !**************************************************************
    ! Initialise the diffusion matrix in 3D

    subroutine initDiffusionMatrix3D(rfparser, rmatrix, rcubatureInfo)

      type(t_fparser), intent(in) :: rfparser
      type(t_matrixScalar), intent(inout) :: rmatrix
      type(t_scalarCubatureInfo), intent(in), optional :: rcubatureInfo

      ! local variables
      type(t_bilinearform) :: rform
      character(LEN=SYS_STRLEN) :: sdiffusionName
      real(DP), dimension(1) :: Dunity = (/1.0_DP/)
      real(DP) :: dalpha
      integer :: i,idiffusiontype

      ! Retrieve data from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'idiffusiontype', idiffusiontype)

      select case(idiffusiontype)
      case (DIFFUSION_ISOTROPIC)

        ! Retrieve name/number of expression describing the diffusion coefficient
        call parlst_getvalue_string(rparlist,&
            ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)

        ! Evaluate the constant coefficient from the function parser
        call fparser_evalFunction(rfparser, sdiffusionName, Dunity, dalpha)

        ! Assemble the Laplace matrix multiplied by the negative value
        ! of the physical diffusion parameter alpha
        call stdop_assembleLaplaceMatrix(rmatrix, .true., -dalpha, rcubatureInfo)

      case (DIFFUSION_ANISOTROPIC)

        do i = 1, 9
          ! Retrieve name/number of expression describing the diffusion coefficient
          call parlst_getvalue_string(rparlist,&
              ssectionName, 'sdiffusionname', sdiffusionName, isubString=i)

          ! Evaluate the constant coefficient from the function parser
          call fparser_evalFunction(rfparser, sdiffusionName, Dunity,&
              rform%Dcoefficients(i))
        end do

        ! We have constant coefficients
        rform%ballCoeffConstant   = .true.
        rform%Dcoefficients(1:9)  = -rform%Dcoefficients(1:9)

        ! Initialise the bilinear form
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

        ! Assemble the anisotropic diffusion matrix
        call bilf_buildMatrixScalar(rform, .true., rmatrix, rcubatureInfo)

      case default
        call lsyssc_clearMatrix(rmatrix)
      end select

    end subroutine initDiffusionMatrix3D

    !**************************************************************
    ! Initialise the group finite element set for evaluating the
    ! bilinear and linear forms on the boundary node-by-node.

    subroutine initGroupFEMSetBoundary(rregion, rmatrix, nmatrices, rgroupFEMSet)

      ! input parameters
      type(t_boundaryRegion), intent(in) :: rregion
      type(t_matrixScalar), intent(in) :: rmatrix
      integer, intent(in) :: nmatrices

      ! input/output parameters
      type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
      
      if (rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
        ! Initialise finite element set for node-based assembly
        call gfem_initGroupFEMSetBoundary(rgroupFEMSet, rmatrix,&
            0, 0, 0, GFEM_NODEBASED, rregionTrial=rregion)
        
        ! Allocate memory for matrix entries
        call gfem_allocCoeffs(rgroupFEMSet, 0, nmatrices, 0)
      else
        ! Resize first group finite element set
        call gfem_resizeGroupFEMSetBoundary(rgroupFEMSet, rmatrix,&
            rregionTrial=rregion)
      end if

!!$      ! Generate diagonal and edge structure derived from template matrix
!!$      call gfem_genDiagList(rmatrix, rgroupFEMSet)
!!$      call gfem_genEdgeList(rmatrix, rgroupFEMSet)

      ! Generate node structure derived from template matrix
      call gfem_genNodeList(rmatrix, rgroupFEMSet)

    end subroutine initGroupFEMSetBoundary
    
  end subroutine transp_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine transp_initAllProblemLevels(rparlist, ssectionName,&
      rproblem, rcollection, rbdrCondPrimal, rbdrCondDual)

!<description>
    ! This subroutine initialises the all problem levels attached to
    ! the global problem structure. It generates the discretisation,
    ! the template matrix and the coefficient matrices as duplicates
    ! of the template matrix.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: boundary condition for primal problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondPrimal

    ! OPTIONAL: boundary condition for dual problem
    type(t_boundaryCondition), intent(in), optional :: rbdrCondDual
!</input>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</intputoutput>
!</subroutine>

    ! pointer to the problem level
    type(t_problemLevel), pointer :: p_rproblemLevel


    ! loop over all problem levels
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(p_rproblemLevel))

      ! Initialise individual problem level
      call transp_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine transp_initAllProblemLevels

  !*****************************************************************************

!<subroutine>

  subroutine transp_initSolution(rparlist, ssectionName,&
    rproblemLevel, dtime, rvector, rcollection)

!<description>
    ! This subroutine initialises the solution vector
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! problem level
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! time for solution evaluation
    real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_afcstab) :: rafcstab
    type(t_collection) :: rcollectionTmp
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    type(t_matrixScalar), pointer :: p_rlumpedMassMatrix,p_rConsistentMassMatrix
    type(t_matrixScalar), target :: rlumpedMassMatrix,rconsistentMassMatrix
    type(t_pgm) :: rpgm
    type(t_vectorBlock) :: rvectorHigh,rvectorAux
    real(DP), dimension(:), pointer :: p_Ddata,p_DdofCoords
    real(DP) :: depsAbsSolution,depsRelSolution,dnorm0,dnorm
    character(LEN=SYS_STRLEN) :: ssolutionname
    integer :: iter,isolutiontype,nmaxIterationsSolution
    integer :: lumpedMassMatrix,consistentMassMatrix,systemMatrix,dofCoords

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isolutiontype', isolutiontype)

    ! How should the solution be initialised?
    select case(isolutionType)
    case (SOLUTION_ZERO)

      !-------------------------------------------------------------------------
      ! Initialise solution by zeros
      !-------------------------------------------------------------------------

      call lsysbl_clearVector(rvector)

      
    case (SOLUTION_GRAYMAP)
      
      !-------------------------------------------------------------------------
      ! Initialise the nodal values by the data of a graymap image
      !-------------------------------------------------------------------------

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
      
      if (dofCoords > 0) then
        ! Set pointers
        call lsyssc_getbase_double(rvector%RvectorBlock(1), p_Ddata)
        call lsyssc_getbase_double(&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1), p_DdofCoords)
              
        ! Initialise solution from portable graymap image
        call ppsol_readPGM(0, ssolutionName, rpgm)
        
        ! Initialise the solution by the image data
        call ppsol_initArrayPGM(rpgm, p_DdofCoords, p_Ddata)
        
        ! Release portable graymap image
        call ppsol_releasePGM(rpgm)
      else
        call output_line('Coordinates of DOFs not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_initSolution')
        call sys_halt()
      end if


    case (SOLUTION_ANALYTIC_POINTVALUE)

      !-------------------------------------------------------------------------
      ! Initialise the nodal values by the data of an analytical expression
      !-------------------------------------------------------------------------
      
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
      
      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)
      
      if (dofCoords > 0) then
        ! Set pointers
        call lsyssc_getbase_double(rvector%RvectorBlock(1), p_Ddata)
        call lsyssc_getbase_double(&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1), p_DdofCoords)

        ! Evaluate solution values in the positions of the degrees of freedom
        call fparser_evalFuncBlockByName2(p_rfparser, ssolutionname,&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1)%NVAR,&
            rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1)%NEQ,&
            p_DdofCoords, rvector%RvectorBlock(1)%NEQ, p_Ddata, (/dtime/))
      else
        call output_line('Coordinates of DOFs not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_initSolution')
        call sys_halt()
      end if

      
    case (SOLUTION_ANALYTIC_L2_CONSISTENT,&
          SOLUTION_ANALYTIC_L2_LUMPED)

      !-------------------------------------------------------------------------
      ! Initialise the FE-function by the L2-projection of the analytical data
      !-------------------------------------------------------------------------

      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist,&
          ssectionName, 'ssolutionname', ssolutionName)

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)

      ! Retrieve the lumped and consistent mass matrices from the
      ! problem level structure or recompute them on-the-fly.
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'consistentMassMatrix', consistentMassMatrix)
      
      if (consistentMassMatrix .gt. 0) then
        p_rconsistentMassMatrix => rproblemLevel%RmatrixScalar(consistentMassMatrix)
      else
        call bilf_createMatrixStructure(&
            rvector%p_rblockDiscr%RspatialDiscr(1),&
            LSYSSC_MATRIX9, rconsistentMassMatrix)
        call stdop_assembleSimpleMatrix(&
            rconsistentMassMatrix, DER_FUNC, DER_FUNC, 1.0_DP, .true.)
        p_rconsistentMassMatrix => rconsistentMassMatrix
      end if
      
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
      
      if (lumpedMassMatrix .gt. 0) then
        p_rlumpedMassMatrix => rproblemLevel%RmatrixScalar(lumpedMassMatrix)
      else
        call lsyssc_duplicateMatrix(p_rconsistentMassMatrix,&
            rlumpedMassMatrix, LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_TEMPLATE)
        call lsyssc_lumpMatrix(rlumpedMassMatrix, LSYSSC_LUMP_DIAG)
        p_rlumpedMassMatrix => rlumpedMassMatrix
      end if
      
      ! Initialise temporal collection structure
      call collct_init(rcollectionTmp)
      
      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, ssolutionname)
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)

      ! Set up the linear form
      rform%itermCount = 1
      rform%Idescriptors(1) = DER_FUNC
      
      ! Assemble the linear form for the scalar subvector
      call linf_buildVectorScalar(rform, .true., rvector%RvectorBlock(1),&
          fcoeff_buildVectorSc_sim=transp_coeffVectorAnalytic,&
          rcollection=rcollectionTmp)

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)

      ! Store norm of load vector (if required)
      dnorm0 = lsyssc_vectorNorm(&
          rvector%RvectorBlock(1), LINALG_NORML2)

      ! Compute the lumped L2-projection
      call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
          rvector%RvectorBlock(1), 1.0_DP, rvector%RvectorBlock(1))

      !-------------------------------------------------------------------------
      ! Restore contribution of the consistent mass matrix of the L2-projection
      !-------------------------------------------------------------------------
      if (isolutionType .eq. SOLUTION_ANALYTIC_L2_CONSISTENT) then

        ! Get configuration from parameter list
        call parlst_getvalue_double(rparlist,&
            ssectionName, 'depsAbsSolution', depsAbsSolution, 1e-6_DP)
        call parlst_getvalue_double(rparlist,&
            ssectionName, 'depsRelSolution', depsRelSolution, 1e-4_DP)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'nmaxIterationsSolution', nmaxIterationsSolution, 100)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'systemMatrix', systemMatrix)

        ! Compute auxiliary vectors for high-order solution and increment
        call lsysbl_duplicateVector(rvector, rvectorHigh,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_COPY)
        call lsysbl_duplicateVector(rvector, rvectorAux,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
        
        ! Compute the consistent L2-projection by Richardson iteration
        richardson: do iter = 1, nmaxIterationsSolution
          ! Compute the increment for each scalar subvector
          call lsyssc_matVec(p_rconsistentMassMatrix,&
              rvectorHigh%RvectorBlock(1),&
              rvectorAux%RvectorBlock(1), 1.0_DP, 0.0_DP)
          call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
              rvectorAux%RvectorBlock(1), 1.0_DP,&
              rvectorAux%RvectorBlock(1))
          call lsyssc_vectorLinearComb(rvector%RvectorBlock(1),&
              rvectorAux%RvectorBlock(1), 1.0_DP, -1.0_DP)
          
          ! Update the scalar subvector of thesolution
          call lsyssc_vectorLinearComb(rvectorAux%RvectorBlock(1),&
              rvectorHigh%RvectorBlock(1), 1.0_DP, 1.0_DP)
          
          ! Check for convergence
          dnorm = lsyssc_vectorNorm(&
              rvectorAux%RvectorBlock(1), LINALG_NORML2)
          if ((dnorm .le. depsAbsSolution) .or.&
              (dnorm .le. depsRelSolution*dnorm0)) exit richardson
        end do richardson
        
        ! Initialise stabilisation structure by hand
        rafcstab%istabilisationSpec = AFCSTAB_UNDEFINED
        rafcstab%cprelimitingType   = AFCSTAB_PRELIMITING_NONE
        rafcstab%cafcstabType = AFCSTAB_LINFCT_MASS
        call afcsc_initStabilisation(rproblemLevel%RmatrixScalar(systemMatrix), rafcstab)

        ! Compute the raw antidiffusive mass fluxes
        call afcsc_buildFluxFCT(rafcstab, rvectorHigh,&
            0.0_DP, 0.0_DP, 1.0_DP, .true., .true.,&
            AFCSTAB_FCTFLUX_EXPLICIT,&
            rmatrix=p_rconsistentMassMatrix, rxTimeDeriv=rvectorHigh)

        ! Apply flux correction to solution profile
        call afcsc_buildVectorFCT(rafcstab,&
            p_rlumpedMassMatrix, rvector, 1.0_DP, .false.,&
            AFCSTAB_FCTALGO_STANDARD+AFCSTAB_FCTALGO_SCALEBYMASS, rvector)

        ! Release stabilisation structure
        call afcstab_releaseStabilisation(rafcstab)

        ! Release auxiliary vectors
        call lsysbl_releaseVector(rvectorHigh)
        call lsysbl_releaseVector(rvectorAux)
      end if

      ! Release temporal matrices (if any)
      call lsyssc_releaseMatrix(rconsistentMassMatrix)
      call lsyssc_releaseMatrix(rlumpedMassMatrix)


    case default
      call output_line('Invalid type of solution profile!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'transp_initSolution')
      call sys_halt()
    end select

  end subroutine transp_initSolution

  !*****************************************************************************

!<subroutine>

  subroutine transp_initRHS(rparlist, ssectionName, rproblemLevel,&
      dtime, rvector, rcollection)

!<description>
    ! This subroutine initialises the right-hand side vector for the
    ! primal problem based on the data from the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! time for right-hand side evaluation
    real(DP), intent(in) :: dtime
!</input>

!<intputoutput>
    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    type(t_collection) :: rcollectionTmp
    character(LEN=SYS_STRLEN) :: srhsname
    integer :: irhstype


    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        'irhstype', irhstype)

    ! How should the right-hand side be initialised?
    select case(irhstype)
    case (RHS_ZERO)
      ! Initialise right-hand side by zeros
      call lsysbl_clearVector(rvector)


    case (RHS_ANALYTIC)
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, ssectionName,&
          'srhsname', srhsname)

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)

      ! Initialise temporal collection structure
      call collct_init(rcollectionTmp)

      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, srhsname)
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)
      
      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretised right-hand side vector
      call linf_buildVectorScalar(rform, .true., rvector%RvectorBlock(1),&
          fcoeff_buildVectorSc_sim=transp_coeffVectorAnalytic,&
          rcollection=rcollectionTmp)

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)


    case default
      call output_line('Invalid type of target functional!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_initRHS')
      call sys_halt()
    end select

  end subroutine transp_initRHS

  !*****************************************************************************

!<subroutine>

  subroutine transp_initTargetFunc(rparlist, ssectionName,&
      rproblemLevel, dtime, rvector, rcollection)

!<description>
    ! This subroutine calculates the target functional which serves as
    ! right-hand side vector for the dual problem in the framework of
    ! goal-oriented error estimation.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! time for target function evaluation
    real(DP), intent(in) :: dtime
!</input>

!<intputoutput>
    ! target function vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_linearForm) :: rform
    type(t_collection) :: rcollectionTmp
    character(LEN=SYS_STRLEN) :: stargetfuncname
    integer :: itargetfunctype

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        'itargetfunctype', itargetfunctype)

    ! How should the target functional be initialised?
    select case(itargetfunctype)
    case (TFUNC_ZERO,&
          TFUNC_SURFINTG)

      ! Initialise target functional by zeros. The contribution of
      ! the surfae integral comes in be weakly impose boundary conditions
      call lsysbl_clearVector(rvector)


    case (TFUNC_VOLINTG,&
          TFUNC_MIXINTG)
      ! Get global configuration from parameter list
      call parlst_getvalue_string(rparlist, ssectionName,&
          'stargetfuncname', stargetfuncname)
      
      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)
      
      ! Initialise temporal collection structure
      call collct_init(rcollectionTmp)
      
      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, stargetfuncname)
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)


      ! Set up the corresponding linear form
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the discretised target functional. The contribution of
      ! the surfae integral comes in be weakly impose boundary conditions
      call linf_buildVectorScalar(rform, .true., rvector%RvectorBlock(1),&
          fcoeff_buildVectorSc_sim=transp_coeffVectorAnalytic,&
          rcollection=rcollectionTmp)

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)


    case default
      call output_line('Invalid type of target functional!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_initTargetFunc')
      call sys_halt()
    end select

  end subroutine transp_initTargetFunc

end module transport_preprocessing
