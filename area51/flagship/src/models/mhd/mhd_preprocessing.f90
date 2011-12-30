!##############################################################################
!# ****************************************************************************
!# <Name> mhd_preprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all preprocessing routines which are required to
!# solve the compressible MHD equations in arbitrary spatial dimensions.
!#
!# The following routines are available:
!#
!# 1.) mhd_initSolvers
!#     -> Initialises the solve structures from the parameter list.
!#
!# 2.) mhd_initProblemDescriptor
!#     -> Initialises the abstract problem descriptor based on the
!#        parameter settings given by the parameter list.
!#
!# 3.) mhd_initProblemLevel
!#     -> Initialises the individual problem level based on the
!#        parameter settings given by the parameter list.
!#        This routine is called repeatedly by the global
!#        initialisation routine mhd_initAllProblemLevels.
!#
!# 4.) mhd_initAllProblemLevels
!#     -> Initialises ALL problem levels attached to the global
!#        problem structure based on the parameter settings
!#        given by the parameter list.
!#
!# 5.) mhd_initSolution
!#     -> Initialises the solution vector based on the parameter
!#        settings given by the parameter list
!#
!# </purpose>
!##############################################################################

module mhd_preprocessing

#include "mhd.h"

  use afcstabbase
  use afcstabsystem
  use basicgeometry
  use bilinearformevaluation
  use boundary
  use boundarycondaux
  use collection
  use cubature
  use derivatives
  use dofmapping
  use element
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use groupfembase
  use groupfemsystem
  use linearalgebra
  use linearformevaluation
  use lineariser
  use linearsystemblock
  use linearsystemscalar
  use meshmodification
  use paramlist
  use pprocerror
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use stdoperators
  use timestep
  use timestepaux
  use triangulation

  ! Modules from MHD model
  use mhd_basic
  use mhd_callback

  implicit none

  private

  public :: mhd_initAllProblemLevels
  public :: mhd_initProblemDescriptor
  public :: mhd_initProblemLevel
  public :: mhd_initSolution
  public :: mhd_initSolvers

contains

  !*****************************************************************************

!<subroutine>

  subroutine mhd_initSolvers(rparlist, ssectionName,&
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
      nlmin = rsolver%p_solverMultigrid%nlmin
      nlmax = rsolver%p_solverMultigrid%nlmax
      call solver_adjustHierarchy(rsolver, nlmin, nlmax)
    else
      call solver_adjustHierarchy(rsolver)
    end if
    call solver_updateStructure(rsolver)

  end subroutine mhd_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine mhd_initProblemDescriptor(rparlist, ssectionName,&
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
    integer :: discretisation,dofCoords
    integer :: massAFC,templateGFEM
    integer :: inviscidAFC,inviscidGFEM
    integer :: viscousAFC,viscousGFEM
    integer :: primalBdrGFEM,dualBdrGFEM
    integer :: templateMatrix
    integer :: systemMatrix
    integer :: jacobianMatrix
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: coeffMatrix_CXX
    integer :: coeffMatrix_CYY
    integer :: coeffMatrix_CZZ
    integer :: coeffMatrix_CXY
    integer :: coeffMatrix_CXZ
    integer :: coeffMatrix_CYZ
    integer :: iconvToTria

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'trifile', rproblemDescriptor%trifile)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'prmfile', rproblemDescriptor%prmfile, '')
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ndimension', rproblemDescriptor%ndimension)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iconvtotria', iconvToTria, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dofCoords', dofCoords, 0)

    ! Get global positions of matrices
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateMatrix', templateMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemMatrix', systemMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianMatrix', jacobianMatrix, 0)
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
        ssectionName, 'coeffMatrix_CXX', coeffMatrix_CXX, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CYY', coeffMatrix_CYY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CZZ', coeffMatrix_CZZ, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CXY', coeffMatrix_CXY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CXZ', coeffMatrix_CXZ, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffMatrix_CYZ', coeffMatrix_CYZ, 0)

    ! Default is no stabilization
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'viscousAFC', viscousAFC, 0)

    ! By default the same identifier is used for the group finite
    ! element formulation and the stabilization structure
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateGFEM', templateGFEM)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'inviscidGFEM', inviscidGFEM, inviscidAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'viscousGFEM', viscousGFEM, viscousAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'primalbdrGFEM', primalbdrGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dualbdrGFEM', dualbdrGFEM, 0)

    ! Consistency check
    if (templateGFEM .eq. 0) then
      inviscidGFEM = 0
      viscousGFEM  = 0
      primalBdrGFEM  = 0
      dualBdrGFEM    = 0
    end if

    ! Set additional problem descriptor
    rproblemDescriptor%ndiscretisation = max(0, discretisation)
    rproblemDescriptor%nafcstab        = max(0, massAFC,&
                                                inviscidAFC,&
                                                viscousAFC)
    rproblemDescriptor%ngroupfemBlock  = max(0, templateGFEM,&
                                                inviscidGFEM,&
                                                viscousGFEM,&
                                                primalBdrGFEM,&
                                                dualBdrGFEM)
    rproblemDescriptor%nlmin           = nlmin
    rproblemDescriptor%nlmax           = nlmax
    rproblemDescriptor%nmatrixScalar   = max(0, templateMatrix,&
                                                consistentMassMatrix,&
                                                lumpedMassMatrix,&
                                                coeffMatrix_CX,&
                                                coeffMatrix_CY,&
                                                coeffMatrix_CZ,&
                                                coeffMatrix_CXX,&
                                                coeffMatrix_CYY,&
                                                coeffMatrix_CZZ,&
                                                coeffMatrix_CXY,&
                                                coeffMatrix_CXZ,&
                                                coeffMatrix_CYZ)
    rproblemDescriptor%nmatrixBlock    = max(0, systemMatrix,&
                                                jacobianMatrix)
    rproblemDescriptor%nvectorScalar   = 0
    rproblemDescriptor%nvectorBlock    = max(0, dofCoords)

    ! Check if quadrilaterals should be converted to triangles
    if (iconvToTria .ne. 0) then
      rproblemDescriptor%iproblemSpec = rproblemDescriptor%iproblemSpec &
                                      + PROBDESC_MSPEC_CONVTRIANGLES
    end if

  end subroutine mhd_initProblemDescriptor

  !*****************************************************************************

!<subroutine>

  subroutine mhd_initProblemLevel(rparlist, ssectionName,&
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
    integer :: consistentMassMatrix
    integer :: lumpedMassMatrix
    integer :: coeffMatrix_CX
    integer :: coeffMatrix_CY
    integer :: coeffMatrix_CZ
    integer :: coeffMatrix_CXX
    integer :: coeffMatrix_CYY
    integer :: coeffMatrix_CZZ
    integer :: coeffMatrix_CXY
    integer :: coeffMatrix_CXZ
    integer :: coeffMatrix_CYZ
    integer :: massAFC,templateGFEM
    integer :: inviscidAFC,inviscidGFEM
    integer :: viscousAFC,viscousGFEM
    integer :: discretisation
    integer :: dofCoords
    integer :: isystemFormat
    integer :: isystemCoupling
    integer :: imatrixFormat
    integer :: primalbdrGFEM,dualbdrGFEM

    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_triangulation) , pointer :: p_rtriangulation
    type(t_groupFEMSet), pointer :: p_rgroupFEMSet
    character(len=SYS_STRLEN) :: slimitingvariable
    type(t_matrixScalar) :: rmatrixSX,rmatrixSY,rmatrixSZ
    integer, dimension(:), allocatable :: Celement
    real(DP) :: dmeshdisturb
    integer :: i,j,ivar,jvar,ivariable,nvariable,nvartransformed,neq
    integer :: nsumcubRefBilForm,nsumcubRefLinForm,nsumcubRefEval
    integer :: nmatrices,nsubstrings,ccubType
    character(len=SYS_STRLEN) :: selemName,smass,sinviscid,sviscous

    ! Retrieve application specific parameters from the parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templatematrix', templateMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'systemmatrix', systemMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'jacobianmatrix', jacobianMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cx', coeffMatrix_CX, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cy', coeffMatrix_CY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cz', coeffMatrix_CZ, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cxx', coeffMatrix_CXX, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cyy', coeffMatrix_CYY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_czz', coeffMatrix_CZZ, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cxy', coeffMatrix_CXY, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cxz', coeffMatrix_CXZ, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cyz', coeffMatrix_CYZ, 0)

    ! Default is no stabilization
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'viscousAFC', viscousAFC, 0)

    ! By default the same identifier is used for the group finite
    ! element formulation and the stabilization structure
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'templateGFEM', templateGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'inviscidGFEM', inviscidGFEM, inviscidAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'viscousGFEM', viscousGFEM, viscousAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'primalBdrGFEM', primalBdrGFEM, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dualBdrGFEM', dualBdrGFEM, 0)

    ! Default no summed cubature
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'discretisation', discretisation, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dofCoords', dofCoords, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'imatrixFormat', imatrixFormat)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemFormat', isystemFormat)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemCoupling', isystemCoupling)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefBilForm', nsumcubRefBilForm, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefLinForm', nsumcubRefLinForm, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'nsumcubRefEval', nsumcubRefEval, 0)
    
    ! Default is empty section, i.e. no configuration
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'mass', smass, '')
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'inviscid', sinviscid, '')
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'viscous', sviscous, '')

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

      ! Initialise the discretisation structure
      p_rdiscretisation => rproblemLevel%Rdiscretisation(discretisation)
      if (p_rdiscretisation%ndimension .eq. 0) then
        select case(isystemFormat)
        case (SYSTEM_INTERLEAVEFORMAT)
          call spdiscr_initBlockDiscr(p_rdiscretisation, 1,&
              rproblemLevel%rtriangulation)

        case (SYSTEM_BLOCKFORMAT)
          call spdiscr_initBlockDiscr(p_rdiscretisation,&
              mhd_getNVAR(rproblemLevel), p_rtriangulation)

        case default
          call output_line('Unsupported system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
          call sys_halt()
        end select
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
        call spdiscr_initDiscr_simple(&
            p_rdiscretisation%RspatialDiscr(1),&
            Celement(1), SPDISC_CUB_AUTOMATIC,&
            p_rtriangulation)

      case (NDIM2D)
        if (size(Celement) .eq. 1) then
          call spdiscr_initDiscr_simple(&
              p_rdiscretisation%RspatialDiscr(1),&
              Celement(1), SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, rproblemLevel%p_rproblem%rboundary)
        else
          call spdiscr_initDiscr_triquad(&
              p_rdiscretisation%RspatialDiscr(1), Celement(1), Celement(2),&
              SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC,&
              p_rtriangulation, rproblemLevel%p_rproblem%rboundary)
        end if

      case (NDIM3D)
        call spdiscr_initDiscr_simple(&
            p_rdiscretisation%RspatialDiscr(1),&
            Celement(1), SPDISC_CUB_AUTOMATIC,&
            p_rtriangulation)
      
      case default
        call output_line('Invalid number of spatial dimensions',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
        call sys_halt()
      end select
      
      ! Deallocate temporal memory
      deallocate(Celement)
      
      !-------------------------------------------------------------------------
      ! Configure evaluation of (bi-)linear formes
      !-------------------------------------------------------------------------

      if (parlst_queryvalue(rparlist, ssectionName, 'ccubTypeBilForm') .ne. 0) then
        ! Check if special cubature formula for evaluating integral
        ! terms of the bilinear form are requested by the user
        nsubstrings = max(1, parlst_querysubstrings(rparlist,&
            ssectionName, 'ccubTypeBilForm'))
        if (nsubstrings .ne. 0) then
          if (nsubstrings .eq. p_rdiscretisation%RspatialDiscr(1)%inumFESpaces) then
            do i = 1, nsubstrings
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'ccubTypeBilForm', ccubType, 0, i)
              if (ccubType .ne. 0)&
                  p_rdiscretisation%RspatialDiscr(1)%RelementDistr(i)%ccubTypeBilForm = ccubType
            end do
          else
            call output_line('Number of substrings does not match number of FE spaces!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
            call sys_halt()
          end if
        end if
      end if
      
      if (parlst_queryvalue(rparlist, ssectionName, 'ccubTypeLinForm') .ne. 0) then
        ! Check if special cubature formula for evaluating integral
        ! terms of the linear form are requested by the user
        nsubstrings = max(1, parlst_querysubstrings(rparlist,&
                             ssectionName, 'ccubTypeLinForm'))
        if (nsubstrings .ne. 0) then
          if (nsubstrings .eq. p_rdiscretisation%RspatialDiscr(1)%inumFESpaces) then
            do i = 1, nsubstrings
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'ccubTypeLinForm', ccubType, 0, i)
              if (ccubType .ne. 0)&
                  p_rdiscretisation%RspatialDiscr(1)%RelementDistr(i)%ccubTypeLinForm = ccubType
            end do
          else
            call output_line('Number of substrings does not match number of FE spaces!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
            call sys_halt()
          end if
        end if
      end if

      if (parlst_queryvalue(rparlist, ssectionName, 'ccubTypeEval') .ne. 0) then
        ! Check if special cubature formula for evaluating integral
        ! terms of allother terms are requested by the user
        nsubstrings = max(1, parlst_querysubstrings(rparlist,&
                             ssectionName, 'ccubTypeEval'))
        if (nsubstrings .ne. 0) then
          if (nsubstrings .eq. p_rdiscretisation%RspatialDiscr(1)%inumFESpaces) then
            do i = 1, nsubstrings
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'ccubTypeEval', ccubType, 0, i)
              if (ccubType .ne. 0)&
                  p_rdiscretisation%RspatialDiscr(1)%RelementDistr(i)%ccubTypeEval = ccubType
            end do
          else
            call output_line('Number of substrings does not match number of FE spaces!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
            call sys_halt()
          end if
        end if
      end if

      ! Duplicate scalar discretisation structure for block matrix format
      if (isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, mhd_getNVAR(rproblemLevel)
          call spdiscr_duplicateDiscrSc(&
              p_rdiscretisation%RspatialDiscr(1),&
              p_rdiscretisation%RspatialDiscr(ivar), .true.)
        end do
      end if

      ! Enforce using summed cubature formula (if any)
      do i = 1, p_rdiscretisation%ncomponents
        do j = 1, p_rdiscretisation%RspatialDiscr(i)%inumFESpaces
          p_rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%ccubTypeBilForm =&
              cub_getSummedCubType(p_rdiscretisation%RspatialDiscr(i)&
              %RelementDistr(j)%ccubTypeBilForm, nsumcubRefBilForm)
          p_rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%ccubTypeLinForm =&
              cub_getSummedCubType(p_rdiscretisation%RspatialDiscr(i)&
              %RelementDistr(j)%ccubTypeLinForm, nsumcubRefLinForm)
          p_rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%ccubTypeEval =&
              cub_getSummedCubType(p_rdiscretisation%RspatialDiscr(i)&
              %RelementDistr(j)%ccubTypeEval, nsumcubRefEval)
        end do
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
            call lsysbl_resizeVectorBlock(&
                rproblemLevel%RvectorBlock(dofCoords), neq, .false.)
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
      if (.not.lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(templateMatrix))) then
        call bilf_createMatrixStructure(&
            p_rdiscretisation%RspatialDiscr(1), imatrixFormat,&
            rproblemLevel%Rmatrix(templateMatrix))
      end if
      
      !-------------------------------------------------------------------------
      ! Create system matrix
      if (systemMatrix > 0) then
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          ! The global operator is stored as an interleave matrix with
          ! NVAR components. However, the row and column structure of
          ! the template matrix can be adopted without modification
          if (lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(systemMatrix))) then
            
            ! Release pseudo block matrix
            call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
            
            ! Resize scalar matrix
            call lsyssc_resizeMatrix(&
                rproblemLevel%Rmatrix(systemMatrix),&
                rproblemLevel%Rmatrix(templateMatrix)%NEQ,&
                rproblemLevel%Rmatrix(templateMatrix)%NCOLS,&
                rproblemLevel%rmatrix(templateMatrix)%NA,&
                .false., .false., bforce=.true.)
            
          else   ! System matrix has no structure
            
            call lsyssc_duplicateMatrix(&
                rproblemLevel%Rmatrix(templateMatrix),&
                rproblemLevel%Rmatrix(systemMatrix),&
                LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
            
            ! Set number of variables per node
            rproblemLevel%Rmatrix(systemMatrix)%NVAR = mhd_getNVAR(rproblemLevel)
            
            ! What matrix format should be used?
            select case(imatrixFormat)
            case (LSYSSC_MATRIX7)
              rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX7INTL
              
            case (LSYSSC_MATRIX9)
              rproblemLevel%Rmatrix(systemMatrix)%cmatrixFormat = LSYSSC_MATRIX9INTL
              
            case default
              call output_line('Unsupported matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
              call sys_halt()
            end select
            
            ! What kind of global operator should be adopted?
            select case(isystemCoupling)
            case (SYSTEM_SEGREGATED)
              rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIXD
              
            case (SYSTEM_ALLCOUPLED)
              rproblemLevel%Rmatrix(systemMatrix)%cinterleavematrixFormat = LSYSSC_MATRIX1
              
            case default
              call output_line('Unsupported interleave matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
              call sys_halt()
            end select
            
            ! Create global operator physically
            call lsyssc_allocEmptyMatrix(&
                rproblemLevel%Rmatrix(systemMatrix), LSYSSC_SETM_UNDEFINED)
            
          end if
          
          ! Create pseudo block matrix from global operator
          call lsysbl_createMatFromScalar(&
              rproblemLevel%Rmatrix(systemMatrix),&
              rproblemLevel%RmatrixBlock(systemMatrix), p_rdiscretisation)
          
          
        case (SYSTEM_BLOCKFORMAT)
          ! The global operator is stored as a block matrix with
          ! NVARxNVAR blocks made up from scalar matrices
          
          if ((rproblemLevel%RmatrixBlock(systemMatrix)%nblocksPerRow .ne. 0) .and.&
              (rproblemLevel%RmatrixBlock(systemMatrix)%nblocksPerCol .ne. 0)) then
            
            ! What kind of global operator should be adopted?
            select case(isystemCoupling)
              
            case (SYSTEM_SEGREGATED)
              ! Create only NVAR diagonal blocks
              do ivar = 1, mhd_getNVAR(rproblemLevel)
                call lsyssc_resizeMatrix(&
                    rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                    rproblemLevel%Rmatrix(templateMatrix), .false., .false., bforce=.true.)
              end do
              
            case (SYSTEM_ALLCOUPLED)
              ! Create all NVAR x NVAR blocks
              do ivar = 1, mhd_getNVAR(rproblemLevel)
                do jvar = 1, mhd_getNVAR(rproblemLevel)
                  call lsyssc_resizeMatrix(&
                      rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar ,jvar),&
                      rproblemLevel%Rmatrix(templateMatrix), .false., .false., bforce=.true.)
                end do
              end do
              
            case default
              call output_line('Unsupported block matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
              call sys_halt()
            end select
            
          else   ! System matrix has no structure
            
            ! Create empty NVARxNVAR block matrix directly
            call lsysbl_createEmptyMatrix(&
                rproblemLevel%RmatrixBlock(systemMatrix),&
                mhd_getNVAR(rproblemLevel),&
                mhd_getNVAR(rproblemLevel))
            
            ! Specify matrix as 'group matrix'
            rproblemLevel%RmatrixBlock(systemMatrix)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX
            
            ! What kind of global operator should be adopted?
            select case(isystemCoupling)
              
            case (SYSTEM_SEGREGATED)
              ! Create only NVAR diagonal blocks
              do ivar = 1, mhd_getNVAR(rproblemLevel)
                call lsyssc_duplicateMatrix(&
                    rproblemLevel%Rmatrix(templateMatrix),&
                    rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
              end do
              
            case (SYSTEM_ALLCOUPLED)
              ! Create all NVAR x NVAR blocks
              do ivar = 1, mhd_getNVAR(rproblemLevel)
                do jvar = 1, mhd_getNVAR(rproblemLevel)
                  call lsyssc_duplicateMatrix(&
                      rproblemLevel%Rmatrix(templateMatrix),&
                      rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                      LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
                end do
              end do
              
            case default
              call output_line('Unsupported block matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
              call sys_halt()
            end select
            
          end if
          
          ! Update internal structure of block matrix
          call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(systemMatrix))
          
        case default
          call output_line('Unsupported system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_initProblemLevel')
          call sys_halt()
        end select
      end if
      
      !-------------------------------------------------------------------------
      ! Create consistent (and lumped) mass matrix as duplicate of the template matrix
      if (consistentMassMatrix > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(consistentMassMatrix))
        call stdop_assembleSimpleMatrix(&
            rproblemLevel%Rmatrix(consistentMassMatrix), DER_FUNC, DER_FUNC)
        
        ! Create lumped mass matrix
        if (lumpedMassMatrix > 0) then
          call lsyssc_duplicateMatrix(&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
          call lsyssc_lumpMatrixScalar(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
        end if
      elseif (lumpedMassMatrix > 0) then
        ! Create lumped mass matrix
        call lsyssc_duplicateMatrix(&
            rproblemLevel%Rmatrix(templateMatrix),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), DER_FUNC, DER_FUNC)
        call lsyssc_lumpMatrixScalar(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dx) as duplicate of the
      ! template matrix
      if (coeffMatrix_CX > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CX))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CX),&
                                        DER_DERIV3D_X, DER_FUNC)
      end if

      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dy) as duplicate of the
      ! template matrix
      if (coeffMatrix_CY > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CY))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CY),&
                                        DER_DERIV3D_Y, DER_FUNC)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (phi, dphi/dz) as duplicate of the
      ! template matrix
      if (coeffMatrix_CZ > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CZ))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CZ),&
                                        DER_DERIV3D_Z, DER_FUNC)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (dphi/dx, dphi/dx) as duplicate of the
      ! template matrix
      if (coeffMatrix_CXX > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CXX))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CXX),&
                                        DER_DERIV3D_X, DER_DERIV3D_X)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (dphi/dy, dphi/dy) as duplicate of the
      ! template matrix
      if (coeffMatrix_CYY > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CYY))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CYY),&
                                        DER_DERIV3D_Y, DER_DERIV3D_Y)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (dphi/dz, dphi/dz) as duplicate of the
      ! template matrix
      if (coeffMatrix_CZZ > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CZZ))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CZZ),&
                                        DER_DERIV3D_Z, DER_DERIV3D_Z)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (dphi/dx, dphi/dy) as duplicate of the
      ! template matrix
      if (coeffMatrix_CXY > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CXY))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CXY),&
                                        DER_DERIV3D_X, DER_DERIV3D_Y)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (dphi/dx, dphi/dz) as duplicate of the
      ! template matrix
      if (coeffMatrix_CXZ > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CXZ))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CXZ),&
                                        DER_DERIV3D_X, DER_DERIV3D_Z)
      end if
      
      !-------------------------------------------------------------------------
      ! Create coefficient matrix (dphi/dy, dphi/dz) as duplicate of the
      ! template matrix
      if (coeffMatrix_CYZ > 0) then
        call initMatrixStructure(rproblemLevel%Rmatrix(templateMatrix),&
                                 rproblemLevel%Rmatrix(coeffMatrix_CYZ))
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(coeffMatrix_CYZ),&
                                        DER_DERIV3D_Y, DER_DERIV3D_Z)
      end if
      
      !-------------------------------------------------------------------------
      ! Create group finite element structures and AFC-stabilisations
      !-------------------------------------------------------------------------
      
      ! Initialise/resize template group finite elment structure and
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
              rproblemLevel%Rmatrix(templateMatrix), 0, 0, 0, GFEM_EDGEBASED)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(templateMatrix))
        end if
        
        ! Generate diagonal and edge structure derived from template matrix
        call gfem_genDiagList(rproblemLevel%Rmatrix(templateMatrix),&
            p_rgroupFEMSet)
        call gfem_genEdgeList(rproblemLevel%Rmatrix(templateMatrix),&
            p_rgroupFEMSet)
      else
        inviscidGFEM = 0
        viscousGFEM  = 0
        primalBdrGFEM  = 0
        dualBdrGFEM    = 0
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/resize group finite element structure as duplicate of
      ! the template group finite element structure and fill it with the
      ! precomputed matrix coefficients for the inviscid term
      if (inviscidGFEM > 0) then
        ! Check if structure has been initialised
        if (rproblemLevel%RgroupFEMBlock(inviscidGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(inviscidGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialise first group finite element set for edge-based
          ! assembly as aduplicate of the template structure
          call gfem_duplicateGroupFEMSet(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              p_rgroupFEMSet, GFEM_DUP_STRUCTURE, .false.)
          
          ! Adjust number of variables
          p_rgroupFEMSet%NVAR = mhd_getNVAR(rproblemLevel)
          
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
              rproblemLevel%Rmatrix(templateMatrix))
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
              rproblemLevel%Rmatrix(coeffMatrix_CX), nmatrices)
        end if
        if (coeffMatrix_CY > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CY), nmatrices)
        end if
        if (coeffMatrix_CZ > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CZ), nmatrices)
        end if
      end if

      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the inviscid term
      ! by duplicating parts of the corresponding group finite element set
      if (inviscidAFC > 0) then
        if (rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          
          ! Initialise stabilisation structure from parameter list
          call afcstab_initFromParameterlist(rparlist, sinviscid,&
              rproblemLevel%Rafcstab(inviscidAFC))

          ! Determine the number of transformed variables
          select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
          case (AFCSTAB_GALERKIN, AFCSTAB_UPWIND)
            ! No variable transformation needs to be performed
            NVARtransformed = 0
            
          case (AFCSTAB_TVD, AFCSTAB_GP)
            ! Transformation to characteristic variables
            NVARtransformed = mhd_getNVAR(rproblemLevel)
            
          case default
            ! Get number of expressions for limiting variables
            nvariable = max(1,&
                parlst_querysubstrings(rparlist,&
                ssectionName, 'slimitingvariable'))
            
            ! Initialise number of limiting variables
            NVARtransformed = 1
            
            ! Determine maximum number of limiting variables in a single set
            do ivariable = 1, nvariable
              call parlst_getvalue_string(rparlist,&
                  ssectionName, 'slimitingvariable',&
                  slimitingvariable, isubstring=ivariable)
              NVARtransformed = max(NVARtransformed,&
                  mhd_getNVARtransformed(rproblemLevel, slimitingvariable))
            end do
          end select
          
          ! Initialise stabilisation structure
          call afcsys_initStabilisation(p_rgroupFEMSet,&
              rproblemLevel%Rafcstab(inviscidAFC), p_rdiscretisation, NVARtransformed)
        else
          ! Resize stabilisation structure
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(inviscidAFC),&
              p_rgroupFEMSet)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/resize group finite element structure as duplicate of
      ! the template group finite element structure and fill it with the
      ! precomputed matrix coefficients for the viscous term
      if (viscousGFEM > 0) then
        ! Check if structure has been initialised
        if (rproblemLevel%RgroupFEMBlock(viscousGFEM)%nblocks .eq. 0)&
            call gfem_initGroupFEMBlock(rproblemLevel%RgroupFEMBlock(viscousGFEM), 1)
        
        ! Set pointer to first group finite element set of this block
        p_rgroupFEMSet =>&
            rproblemLevel%RgroupFEMBlock(viscousGFEM)%RgroupFEMBlock(1)
        
        if (p_rgroupFEMSet%isetSpec .eq. GFEM_UNDEFINED) then
          ! Initialise first group finite element set for edge-based
          ! assembly as aduplicate of the template structure
          call gfem_duplicateGroupFEMSet(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              p_rgroupFEMSet, GFEM_DUP_STRUCTURE, .false.)
          
          ! Adjust number of variables
          p_rgroupFEMSet%NVAR = mhd_getNVAR(rproblemLevel)
          
          ! Compute number of matrices to be copied
          nmatrices = 0
          if (coeffMatrix_CXX > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CYY > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CZZ > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CXY > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CXZ > 0) nmatrices = nmatrices+1
          if (coeffMatrix_CYZ > 0) nmatrices = nmatrices+1
          
          ! Allocate memory for matrix entries
          call gfem_allocCoeffs(p_rgroupFEMSet, nmatrices, 0, nmatrices)
        else
          ! Resize first group finite element set
          call gfem_resizeGroupFEMSet(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(templateMatrix))
        end if
        
        ! Duplicate edge-based structure from template
        call gfem_duplicateGroupFEMSet(&
            rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
            p_rgroupFEMSet, GFEM_DUP_DIAGLIST+GFEM_DUP_EDGELIST, .true.)
        
        ! Copy constant coefficient matrices to group finite element set
        nmatrices = 0
        if (coeffMatrix_CXX > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CXX), nmatrices)
        end if
        if (coeffMatrix_CYY > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CYY), nmatrices)
        end if
        if (coeffMatrix_CZZ > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CZZ), nmatrices)
        end if
        if (coeffMatrix_CXY > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CXY), nmatrices)
        end if
        if (coeffMatrix_CXZ > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CXZ), nmatrices)
        end if
        if (coeffMatrix_CYZ > 0) then
          nmatrices = nmatrices+1
          call gfem_initCoeffsFromMatrix(p_rgroupFEMSet,&
              rproblemLevel%Rmatrix(coeffMatrix_CYZ), nmatrices)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the inviscid term
      ! by duplicating parts of the corresponding group finite element set
      if (viscousAFC > 0) then
        if (rproblemLevel%Rafcstab(viscousAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          
          ! Get number of expressions for limiting variables
          nvariable = max(1,&
              parlst_querysubstrings(rparlist,&
              ssectionName, 'slimitingvariable'))
          
          ! Initialise number of limiting variables
          NVARtransformed = 1
          
          ! Determine maximum number of limiting variables in a single set
          do ivariable = 1, nvariable
            call parlst_getvalue_string(rparlist,&
                ssectionName, 'slimitingvariable',&
                slimitingvariable, isubstring=ivariable)
            NVARtransformed = max(NVARtransformed,&
                mhd_getNVARtransformed(rproblemLevel, slimitingvariable))
          end do
          
          ! Initialise stabilisation structure
          call afcstab_initFromParameterlist(rparlist, sviscous,&
              rproblemLevel%Rafcstab(viscousAFC))
          call afcsys_initStabilisation(p_rgroupFEMSet,&
              rproblemLevel%Rafcstab(viscousAFC), p_rdiscretisation, NVARtransformed)
        else
          ! Resize stabilisation structure
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(viscousAFC),&
              p_rgroupFEMSet)
        end if
      end if
      
      !-------------------------------------------------------------------------
      ! Initialise/Resize stabilisation structure for the mass matrix by
      ! duplicating parts of the template group finite element set
      if (massAFC > 0) then
        if (rproblemLevel%Rafcstab(massAFC)%istabilisationSpec&
            .eq. AFCSTAB_UNDEFINED) then
          
          ! Get number of expressions for limiting variables
          nvariable = max(1,&
              parlst_querysubstrings(rparlist,&
              ssectionName, 'slimitingvariable'))
          
          ! Initialise number of limiting variables
          NVARtransformed = 1
          
          ! Determine maximum number of limiting variables in a single set
          do ivariable = 1, nvariable
            call parlst_getvalue_string(rparlist,&
                ssectionName, 'slimitingvariable',&
                slimitingvariable, isubstring=ivariable)
            NVARtransformed = max(NVARtransformed,&
                mhd_getNVARtransformed(rproblemLevel, slimitingvariable))
          end do
          
          ! Initialise stabilisation structure
          call afcstab_initFromParameterlist(rparlist, smass,&
              rproblemLevel%Rafcstab(massAFC))
          call afcsys_initStabilisation(&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rafcstab(massAFC), p_rdiscretisation, NVARtransformed)
          
          ! Adjust number of variables
          rproblemLevel%Rafcstab(massAFC)%NVAR = mhd_getNVAR(rproblemLevel)
        else
          ! Resize stabilisation structure
          call afcstab_resizeStabilisation(rproblemLevel%Rafcstab(massAFC),&
              rproblemLevel%RgroupFEMBlock(templateGFEM)%RgroupFEMBlock(1))
        end if
      end if
      
      ! The auxiliary matrices are not needed any more
      if (coeffMatrix_CX > 0)&
          call lsyssc_releaseMatrix(rmatrixSX)
      if (coeffMatrix_CY > 0)&
          call lsyssc_releaseMatrix(rmatrixSY)
      if (coeffMatrix_CZ > 0)&
          call lsyssc_releaseMatrix(rmatrixSZ)
      
    end if   ! discretisation > 0

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
            0, 0, 0, GFEM_NODEBASED, rregionTest=rregion,&
            brestrictToBoundary=.true.)
        
        ! Allocate memory for matrix entries
        call gfem_allocCoeffs(rgroupFEMSet, 0, nmatrices,0)
      else
        ! Resize first group finite element set
        call gfem_resizeGroupFEMSetBoundary(rgroupFEMSet, rmatrix,&
            rregionTest=rregion, rregionTrial=rregion)
      end if

      ! Generate node structure derived from template matrix
      call gfem_genNodeList(rmatrix, rgroupFEMSet)

    end subroutine initGroupFEMSetBoundary

  end subroutine mhd_initProblemLevel

  !*****************************************************************************

!<subroutine>

  subroutine mhd_initAllProblemLevels(rparlist, ssectionName,&
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
      call mhd_initProblemLevel(rparlist, ssectionName,&
          p_rproblemLevel, rcollection, rbdrCondPrimal, rbdrCondDual)

      ! Switch to next coarser level
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine mhd_initAllProblemLevels

  !*****************************************************************************

!<subroutine>

  subroutine mhd_initSolution(rparlist, ssectionName,&
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
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_afcstab) :: rafcstab
    type(t_linearForm) :: rform
    type(t_vectorBlock) :: rvectorBlock,rvectorHigh,rvectorAux
    type(t_vectorScalar) :: rvectorScalar
    type(t_matrixScalar), target :: rlumpedMassMatrix,rconsistentMassMatrix
    type(t_matrixScalar), pointer :: p_rlumpedMassMatrix,p_rConsistentMassMatrix
    type(t_fparser), pointer :: p_rfparser
    type(t_collection) :: rcollectionTmp
    real(DP), dimension(:), pointer :: p_Ddata,p_DdofCoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP), dimension(:), pointer :: Dnorm0
    real(DP) :: depsAbsSolution,depsRelSolution,dnorm,dmin,dmax
    character(len=SYS_STRLEN), dimension(:), pointer :: SsolutionFailsafeVariables
    character(LEN=SYS_STRLEN) :: ssolutionName
    integer :: isolutiontype,nexpression,nsolutionfailsafe
    integer :: iblock,ivar,ieq,idim,iter
    integer :: lumpedMassMatrix,consistentMassMatrix,systemMatrix,dofCoords
    integer :: nmaxIterationsSolution,ivariable,nvariable
    logical :: bisAccepted


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

      call output_line('Initialisation if solution by graymap image is not yet supported!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'mhd_initSolution')
      call sys_halt()


    case (SOLUTION_ANALYTIC_POINTVALUE)

      !-------------------------------------------------------------------------
      ! Initialise the nodal values by the data of an analytical expression
      !-------------------------------------------------------------------------

      ! Initialise total number of expressions
      nexpression = 0

      ! Compute total number of expressions
      do iblock = 1, rvector%nblocks
        nexpression = nexpression + rvector%RvectorBlock(iblock)%NVAR
      end do

      ! Check if array of solution names is available
      if (parlst_querysubstrings(rparlist, ssectionName,&
          'ssolutionname') .lt. nexpression) then
        call output_line('Invalid number of expressions!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'mhd_initSolution')
        call sys_halt()
      end if

      ! Get global configuration from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'dofCoords', dofCoords)
      
      if (dofCoords > 0) then
        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)
        
        ! Initialise variable values
        Dvalue           = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        nexpression      = 0
        
        ! Loop over all blocks of the global solution vector
        do iblock = 1, rvector%nblocks
          
          ! Set pointers
          call lsyssc_getbase_double(rvector%RvectorBlock(iblock), p_Ddata)
          call lsyssc_getbase_double(&
              rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(iblock),&
              p_DdofCoords)
          
          ! Loop over all variables of the solution vector
          do ivar = 1, rvector%RvectorBlock(iblock)%NVAR
            
            ! Get the function name of the component used for evaluating the initial solution.
            call parlst_getvalue_string(rparlist, ssectionName,&
                'ssolutionName', ssolutionName, isubstring=nexpression+ivar)
            
            if (rvector%RvectorBlock(iblock)%NVAR .eq. 1) then
              ! Evaluate the function parser for all coordinates
              call fparser_evalFuncBlockByName2(p_rfparser, ssolutionName,&
                  rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(iblock)%NVAR,&
                  rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(iblock)%NEQ,&
                  p_DdofCoords, rvector%RvectorBlock(iblock)%NEQ, p_Ddata, (/dtime/))
            else
              
              ! Loop over all equations of the scalar subvector
              do ieq = 1, rvector%RvectorBlock(iblock)%NEQ
                ! Set coordinates and evalution time
                do idim = 1, rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(iblock)%NVAR
                  Dvalue(idim) = p_DdofCoords((ieq-1)*&
                      rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(iblock)%NVAR+idim)
                end do
                
                ! Evaluate the function parser
                call fparser_evalFunction(p_rfparser, ssolutionName,&
                    Dvalue, p_Ddata((ieq-1)*rvector%RvectorBlock(iblock)%NVAR+ivar))
              end do   ! ieq
            end if
            
          end do   ! ivar
          
          ! Increase number of processed expressions
          nexpression = nexpression + rvector%RvectorBlock(iblock)%NVAR
          
        end do   ! iblock

      else
        call output_line('Coordinates of DOFs not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_initSolution')
        call sys_halt()
      end if
   

    case (SOLUTION_ANALYTIC_L2_CONSISTENT,&
          SOLUTION_ANALYTIC_L2_LUMPED)

      !-------------------------------------------------------------------------
      ! Initialise the FE-function by the L2-projection of the analytical data
      !-------------------------------------------------------------------------

      ! Initialise total number of expressions
      nexpression = 0

      ! Compute total number of expressions
      do iblock = 1, rvector%nblocks
        nexpression = nexpression + rvector%RvectorBlock(iblock)%NVAR
      end do

      ! Check if array of solution names is available
      if (parlst_querysubstrings(rparlist, ssectionName,&
          'ssolutionname') .lt. nexpression) then
        call output_line('Invalid number of expressions!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'mhd_initSolution')
        call sys_halt()
      end if

      ! Get function parser from collection structure
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)

      ! Retrieve the lumped and consistent mass matrices from the
      ! problem level structure or recompute them on-the-fly.
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'consistentMassMatrix', consistentMassMatrix)
      
      if (consistentMassMatrix .gt. 0) then
        p_rconsistentMassMatrix => rproblemLevel%Rmatrix(consistentMassMatrix)
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
        p_rlumpedMassMatrix => rproblemLevel%Rmatrix(lumpedMassMatrix)
      else
        call lsyssc_duplicateMatrix(p_rconsistentMassMatrix,&
            rlumpedMassMatrix, LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_TEMPLATE)
        call lsyssc_lumpMatrixScalar(rlumpedMassMatrix, LSYSSC_LUMP_DIAG)
        p_rlumpedMassMatrix => rlumpedMassMatrix
      end if
      
      ! Initialise temporal collection structure
      call collct_init(rcollectionTmp)
      
      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)

      ! Set up the linear form
      rform%itermCount = 1
      rform%Idescriptors(1) = DER_FUNC
      
      ! Initialise number of expressions
      nexpression = 0

      ! Loop over all blocks of the global solution vector
      do iblock = 1, rvector%nblocks
                
        ! Scalar vectors in interleaved format have to be treated differently
        if (rvector%RvectorBlock(iblock)%NVAR .eq. 1) then

          ! Get the function name of the component used for evaluating the initial solution.
          call parlst_getvalue_string(rparlist, ssectionName,&
              'ssolutionName', ssolutionName, isubstring=nexpression+1)

          ! Set the number of the component used for evaluating the initial solution
          rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, ssolutionname)
          
          ! Assemble the linear form for the scalar subvector
          call linf_buildVectorScalar2(rform, .true.,&
              rvector%RvectorBlock(iblock), mhd_coeffVectorAnalytic, rcollectionTmp)

          ! Increase number of processed expressions
          nexpression = nexpression + 1
        else

          ! Convert scalar vector in interleaved format to true block vector
          call lsysbl_convertScalarBlockVector(&
              rvector%RvectorBlock(iblock), rvectorBlock)

          ! Loop over all blocks
          do ivar = 1, rvectorBlock%nblocks
            
            ! Get the function name of the component used for evaluating the initial solution.
            call parlst_getvalue_string(rparlist, ssectionName,&
                'ssolutionName', ssolutionName, isubstring=nexpression+ivar)
            
            ! Set the number of the component used for evaluating the initial solution
            rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, ssolutionname)

            ! Assemble the linear form for the scalar subvector
            call linf_buildVectorScalar2(rform, .true.,&
                rvectorBlock%RvectorBlock(ivar), mhd_coeffVectorAnalytic, rcollectionTmp)
          end do

          ! Convert block vector back to scalar vector in interleaved format
          call lsysbl_convertBlockScalarVector(&
              rvectorBlock,rvector%RvectorBlock(iblock))
          
          ! Increase number of processed expressions
          nexpression = nexpression + rvectorBlock%nblocks
          
          ! Release temporal block vector
          call lsysbl_releaseVector(rvectorBlock)
        end if

      end do

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)

      ! Store norm of load vector (if required)
      if (isolutionType .eq. SOLUTION_ANALYTIC_L2_CONSISTENT) then
        allocate(Dnorm0(rvector%nblocks))
        do iblock = 1, rvector%nblocks
          Dnorm0(iblock) = lsyssc_vectorNorm(&
              rvector%RvectorBlock(iblock), LINALG_NORML2)
        end do
      end if

      ! Compute the lumped L2-projection
      do iblock = 1, rvector%nblocks
        call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
            rvector%RvectorBlock(iblock), 1.0_DP, rvector%RvectorBlock(iblock))
      end do

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
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'nsolutionfailsafe', nsolutionfailsafe, 0)

        ! Compute auxiliary vectors for high-order solution and increment
        call lsysbl_duplicateVector(rvector, rvectorHigh,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_COPY)
        call lsysbl_duplicateVector(rvector, rvectorAux,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
        
        ! Compute the consistent L2-projection by Richardson iteration
        do iblock = 1, rvector%nblocks
          richardson: do iter = 1, nmaxIterationsSolution
            ! Compute the increment for each scalar subvector
            call lsyssc_scalarMatVec(p_rconsistentMassMatrix,&
                rvectorHigh%RvectorBlock(iblock),&
                rvectorAux%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
            call lsyssc_invertedDiagMatVec(p_rlumpedMassMatrix,&
                rvectorAux%RvectorBlock(iblock), 1.0_DP,&
                rvectorAux%RvectorBlock(iblock))
            call lsyssc_vectorLinearComb(rvector%RvectorBlock(iblock),&
                rvectorAux%RvectorBlock(iblock), 1.0_DP, -1.0_DP)

            ! Update the scalar subvector of thesolution
            call lsyssc_vectorLinearComb(rvectorAux%RvectorBlock(iblock),&
                rvectorHigh%RvectorBlock(iblock), 1.0_DP, 1.0_DP)
            
            ! Check for convergence
            dnorm = lsyssc_vectorNorm(&
                rvectorAux%RvectorBlock(iblock), LINALG_NORML2)
            if ((dnorm .le. depsAbsSolution) .or.&
                (dnorm .le. depsRelSolution*Dnorm0(iblock))) exit richardson
          end do richardson
        end do
        
        ! Initialise stabilisation structure by hand
        rafcstab%istabilisationSpec = AFCSTAB_UNDEFINED
        rafcstab%cprelimitingType   = AFCSTAB_PRELIMITING_NONE
        rafcstab%cafcstabType = AFCSTAB_LINFCT_MASS
        call afcsys_initStabilisation(&
            rproblemLevel%RmatrixBlock(systemMatrix), rafcstab)

        ! Compute the raw antidiffusive mass fluxes. Note that we may supply any
        ! callback function for assembling the antidiffusive fluxes since it
        ! will not be used for assembling antidiffusive mass fluxes !!!
        call afcsys_buildFluxFCT(rafcstab, rvectorHigh,&
            0.0_DP, 0.0_DP, 1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
            rmatrix=p_rconsistentMassMatrix, rcollection=rcollection)
        
        if (nsolutionfailsafe .gt. 0) then

          ! Get number of failsafe variables
          nvariable = max(1,&
              parlst_querysubstrings(rparlist,&
              ssectionName, 'ssolutionfailsafevariable'))
          
          ! Allocate character array that stores all failsafe variable names
          allocate(SsolutionFailsafeVariables(nvariable))
          
          ! Initialise character array with failsafe variable names
          do ivariable = 1, nvariable
            call parlst_getvalue_string(rparlist,&
                ssectionName, 'ssolutionfailsafevariable',&
                SsolutionFailsafevariables(ivariable), isubstring=ivariable)
          end do

          ! Compute and apply FEM-FCT correction
          call mhd_calcCorrectionFCT(rproblemLevel, rvector, 1.0_DP,&
              .false., AFCSTAB_FCTALGO_STANDARD-AFCSTAB_FCTALGO_CORRECT,&
              rvector, ssectionName, rcollection,&
              rafcstab, 'ssolutionconstrainvariable')
          
          ! Apply failsafe flux correction
          call afcsys_failsafeFCT(rafcstab, p_rlumpedMassMatrix,&
              rvector, 1.0_DP, 1e-8_DP, AFCSTAB_FAILSAFEALGO_STANDARD,&
              bisAccepted, nsteps=nsolutionfailsafe,&
              CvariableNames=SsolutionFailsafeVariables,&
              fcb_extractVariable=mhd_getVariable,&
              rcollection=rcollection)
          
          ! Deallocate temporal memory
          deallocate(SsolutionFailsafeVariables)

        else
          
          ! Compute and apply FEM-FCT correction
          call mhd_calcCorrectionFCT(rproblemLevel, rvector, 1.0_DP,&
              .false., AFCSTAB_FCTALGO_STANDARD+AFCSTAB_FCTALGO_SCALEBYMASS,&
              rvector, ssectionName, rcollection,&
              rafcstab, 'ssolutionconstrainvariable')
        end if
        
        ! Release stabilisation structure
        call afcstab_releaseStabilisation(rafcstab)

        ! Release auxiliary vectors
        call lsysbl_releaseVector(rvectorHigh)
        call lsysbl_releaseVector(rvectorAux)

        ! Release temporal memory
        deallocate(Dnorm0)
      end if

      ! Release temporal matrices (if any)
      call lsyssc_releaseMatrix(rconsistentMassMatrix)
      call lsyssc_releaseMatrix(rlumpedMassMatrix)
        

    case default
      call output_line('Invalid type of solution profile!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'mhd_initSolution')
      call sys_halt()
    end select

    ! Output statistics
    call output_lbrk()
    call output_separator(OU_SEP_AT)
    call output_line('Initial solution')
    call output_separator(OU_SEP_AT)

    ! Density
    call mhd_getVariable(rvector, 'density', rvectorScalar)
    call pperr_scalar(PPERR_L1ERROR, dnorm, rvectorScalar)
    call lsyssc_getbase_double(rvectorScalar, p_Ddata)
    dmin = SYS_MAXREAL_DP; dmax = -SYS_MAXREAL_DP
    
    do ieq = 1, rvectorScalar%NEQ
      dmin = min(dmin, p_Ddata(ieq))
      dmax = max(dmax, p_Ddata(ieq))
    end do

    call output_line('Total mass:       '//trim(sys_sdEL(dnorm,5)))
    call output_line('min/max density:  '//trim(sys_sdEL(dmin,5))//' / '&
                                         //trim(sys_sdEL(dmax,5)))
    
    ! Total energy
    call mhd_getVariable(rvector, 'total_energy', rvectorScalar)
    call pperr_scalar(PPERR_L1ERROR, dnorm, rvectorScalar)
    call lsyssc_getbase_double(rvectorScalar, p_Ddata)
    dmin = SYS_MAXREAL_DP; dmax = -SYS_MAXREAL_DP
    
    do ieq = 1, rvectorScalar%NEQ
      dmin = min(dmin, p_Ddata(ieq))
      dmax = max(dmax, p_Ddata(ieq))
    end do

    call output_line('Total energy:     '//trim(sys_sdEL(dnorm,5)))
    call output_line('min/max density:  '//trim(sys_sdEL(dmin,5))//' / '&
                                         //trim(sys_sdEL(dmax,5)))
    
    ! Pressure
    call mhd_getVariable(rvector, 'pressure', rvectorScalar)
    call pperr_scalar(PPERR_L1ERROR, dnorm, rvectorScalar)
    call lsyssc_getbase_double(rvectorScalar, p_Ddata)
    dmin = SYS_MAXREAL_DP; dmax = -SYS_MAXREAL_DP
    
    do ieq = 1, rvectorScalar%NEQ
      dmin = min(dmin, p_Ddata(ieq))
      dmax = max(dmax, p_Ddata(ieq))
    end do

    call output_line('min/max pressure: '//trim(sys_sdEL(dmin,5))//' / '&
                                         //trim(sys_sdEL(dmax,5)))

    call output_separator(OU_SEP_AT)
    call output_lbrk()

    ! Release temporal memory
    call lsyssc_releaseVector(rvectorScalar)
  end subroutine mhd_initSolution

end module mhd_preprocessing
