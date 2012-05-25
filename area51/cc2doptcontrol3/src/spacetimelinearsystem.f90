!##############################################################################
!# ****************************************************************************
!# <name> spacetimelinearsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module spacetimelinearsystem

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use basicgeometry
  use boundary
  use linearalgebra
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use linearsystemscalar
  use linearsystemblock
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use paramlist
  use timestepping
  use discretebc
  use discretefbc
  use domainintegration
  use analyticprojection
  use element
  
  use collection
  use convection
    
  use constantsdiscretisation
  use assemblytemplates
  use assemblytemplatesoptc
  use structuresoptc

  use analyticsolution
  use spatialbcdef
  use spacetimevectors
  use timediscretisation
  use timeboundaryconditions
  use spacediscretisation
  !use spacetimediscretisation
  use spacematvecassembly
  use fespacehierarchybase
  use fespacehierarchy
  
  use spacetimehierarchy
  
  use spacetimeneumannbc
  use spacetimedirichletbcc

  !use spacepreconditioner
  !use spacepreconditionerinit

  !use timeanalysis
  !use spatialbc
  
  !use spacetimediscretisation
  !use timeboundaryconditions
  use spacematvecassembly
  use structuresoptflow
  use user_callback
  
  use matrixio
  use vectorio
    
  implicit none

  private

!<constants>

!<constantblock description="Filter-ID's of the filters to be applied in a matrix vector multiplication">

  ! No filter
  integer(I32), parameter, public :: SPTID_FILTER_NONE  = 0

  ! Filter the output vector for boundary conditions in the defect.
  integer(I32), parameter, public :: SPTID_FILTER_BCDEF = 2**0
  
  ! Filter the output vector for initial conditions in the defect
  integer(I32), parameter, public :: SPTID_FILTER_ICDEF = 2**1

  ! Filter the output vector for terminal conditions in the defect
  integer(I32), parameter, public :: SPTID_FILTER_TCDEF = 2**2
  
  ! Apply all filters that are typically applied to a defect vector
  integer(I32), parameter, public :: SPTID_FILTER_DEFECT = SPTID_FILTER_BCDEF + &
                                                   SPTID_FILTER_ICDEF + &
                                                   SPTID_FILTER_TCDEF
                                                   
!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! All discretisation and level related data structures that are necessary
  ! to discretise equations in the KKT system.
  type t_spaceTimeDiscrData
  
    ! Pointer to the underlying physics.
    type(t_settings_physics), pointer :: p_rphysicsPrimal => null()
    
    ! Pointer to the stabilisation structures of the primal and dual
    ! equation.
    type(t_settings_stabil), pointer :: p_rstabilPrimal => null()
    type(t_settings_stabil), pointer :: p_rstabilDual => null()
    
    ! Possible constraints in space/time. May point to NULL if there are
    ! no constraints.
    type(t_optcconstraintsSpaceTime), pointer :: p_rconstraints => null()
    
    ! Discretisation of the level, primal space.
    type(t_blockDiscretisation), pointer :: p_rdiscrPrimal => null()

    ! Discretisation of the level, dual space.
    type(t_blockDiscretisation), pointer :: p_rdiscrDual => null()

    ! Discretisation of the level, control space.
    type(t_blockDiscretisation), pointer :: p_rdiscrDual => null()
    
    ! Discretisation of the time mesh.
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()
    
    ! Pointer to optimal control parameters.
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl => null()
    
    ! Pointer to the static assembly data on the corresponding space level.
    type(t_staticSpaceAsmTemplates), pointer :: p_rstaticSpaceAsmTempl => null()

    ! Pointer to the static optimal control assembly data on the corresponding space level.
    type(t_staticSpaceAsmTemplatesOptC), pointer :: p_rstaticSpaceAsmTemplOptC => null()
    
  end type

!</typeblock>

  public :: t_spaceTimeMatrixDiscrData


!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  type t_ccoptSpaceTimeMatrix
  
    ! Type of the matrix. One of the MATT_xxxx constants
    integer :: cmatrixType = MATT_OPTCONTROL
    
    ! Number of unknowns in time.
    integer :: NEQtime = 0
    
    ! <!--
    !  This matrix realises a nonlinear operator in time! For the assembly it is
    !  therefore necessary to carry with all information which is necessary
    !  to handle the nonlinearity -- physics, stabilisation, discretisation, etc.
    ! -->
    
    ! Discretisation related data.
    type(t_spaceTimeMatrixDiscrData) :: rdiscrData

    ! Pointer to global data.
    type(t_globalData), pointer :: p_rglobalData => null()
  
    ! Pointer to program wide debug flags.
    type(t_optcDebugFlags), pointer :: p_rdebugFlags => null()
    
    ! Pointer to a space-time solution vector that defines the point
    ! where the nonlinearity is evaluated when applying or calculating
    ! matrices.
    ! If this points to NULL(), there is no global solution specified,
    ! so the matrix assembly/application routines use the vector that
    ! is specified in their input parameters.
    type(t_spacetimeVector), pointer :: p_rsolution => null()
    
    ! A t_sptiNeumannBoundary structure defining the Neumann boundary
    type(t_sptiNeumannBoundary), pointer :: p_rneumannBoundary => null()

    ! A t_sptiDirichletBCCBoundary structure defining the Dirichlet 
    ! boundary control
    type(t_sptiDirichletBCCBoundary), pointer :: p_rdirichletBCCBoundary => null()
    
  end type

!</typeblock>

!</types>

  public :: t_ccoptSpaceTimeMatrix
  public :: stlin_initSpaceTimeMatrix
  public :: stlin_releaseSpaceTimeMatrix
  public :: stlin_setupMatrixWeights
  public :: stlin_getFullMatrixDummy
  public :: stlin_initSpaceConstraints
  public :: stlin_doneSpaceConstraints
  public :: stlin_initSpaceAssembly
  public :: stlin_initSpaceAssemblyFromGl
  public :: stlin_doneSpaceAssembly
  public :: stlin_spaceTimeMatVec
  public :: stlin_dampDualSolution

contains

  ! ***************************************************************************

    subroutine ffunctionConstraint (cderivative, rdiscretisation, &
                                   nelements, npointsPerElement, Dpoints, &
                                   IdofsTest, rdomainIntSubset, &
                                   Dvalues, rcollection)
    
  !<description>
    ! Reference function used for the L2 projetion of analytically given
    ! constraints.
  !</description>
    
    integer, intent(in) :: cderivative
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    integer, intent(in) :: nelements
    integer, intent(in) :: npointsPerElement
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    integer, dimension(:,:), intent(in) :: IdofsTest
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
    type(t_collection), intent(inout), optional :: rcollection
    
    real(DP), dimension(:,:), intent(out) :: Dvalues
  
      integer :: ierror, idim
  
      idim = rcollection%IquickAccess(1)
  
      ! Evaluate the function "SOL" from the collection
      call ansol_evaluate (rcollection,"SOL",idim,Dvalues,&
          npointsPerElement,nelements,Dpoints,ierror=ierror)
          
      if (ierror .ne. 0) then
        call output_line ("Error evaluating constraint for projection!");
      end if
  
    end subroutine
    
  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_initSpaceConstraints (&
      rspaceTimeConstr,dcontrolConstrTime,dstateConstrTime,rspaceDiscrPrimal,rspaceConstr)

!<description>
  ! Creates a discrete version of the constraints at a specific point in time.
  ! The routine can be called multiple times; every new call updates the
  ! data in rspaceConstr.
!</description>

!<input>
  ! Space-time assembly structure.
  type(t_optcconstraintsSpaceTime), intent(in) :: rspaceTimeConstr
  
  ! Current time where the control constraints should be applied
  real(DP), intent(in) :: dcontrolConstrTime

  ! Current time where the state constraints should be applied
  real(DP), intent(in) :: dstateConstrTime
  
  ! Spatial discretisation structure of the primal space.
  type(t_blockDiscretisation), intent(in) :: rspaceDiscrPrimal
!</input>

!<output>
  ! Space constraints structure.
  type(t_optcconstraintsSpace), intent(out) :: rspaceConstr
!</output>

!</subroutine>

    ! local variables
    type(t_collection), target :: rcollection
    integer :: idim

    ! Transfer all possible parameters.
    !
    ! Control constraints.
    rspaceConstr%dumin1 = rspaceTimeConstr%dumin1
    rspaceConstr%dumax1 = rspaceTimeConstr%dumax1
    rspaceConstr%dumin2 = rspaceTimeConstr%dumin2
    rspaceConstr%dumax2 = rspaceTimeConstr%dumax2
    rspaceConstr%ccontrolConstraints = &
        rspaceTimeConstr%ccontrolConstraints

    rspaceConstr%p_rumin1 => rspaceTimeConstr%p_rumin1
    rspaceConstr%p_rumax1 => rspaceTimeConstr%p_rumax1
    rspaceConstr%p_rumin2 => rspaceTimeConstr%p_rumin2
    rspaceConstr%p_rumax2 => rspaceTimeConstr%p_rumax2
    rspaceConstr%ccontrolConstraintsType = &
        rspaceTimeConstr%ccontrolConstraintsType

    rspaceConstr%dcontrolConstrTime = dcontrolConstrTime

    select case (rspaceConstr%ccontrolConstraints)
    case (1:)
      select case (rspaceConstr%ccontrolConstraintsType)
      case (1)
        if (.not. associated(rspaceConstr%p_rvectorumin)) then
          allocate(rspaceConstr%p_rvectorumin)
          call lsysbl_createVectorBlock(rspaceDiscrPrimal,rspaceConstr%p_rvectorumin,.true.)
        end if
        
        if (.not. associated(rspaceConstr%p_rvectorumax)) then
          allocate(rspaceConstr%p_rvectorumax)
          call lsysbl_createVectorBlock(rspaceDiscrPrimal,rspaceConstr%p_rvectorumax,.true.)
        end if
      
        call collct_init (rcollection)

        ! -----
        ! U1-min
        call ansol_prepareEval (rspaceConstr%p_rumin1,rcollection,"SOL",dcontrolConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 1
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorumin%RvectorBlock(1),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        ! -----
        ! U1-max
        call ansol_prepareEval (rspaceConstr%p_rumax1,rcollection,"SOL",dcontrolConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 1
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorumax%RvectorBlock(1),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        ! -----
        ! U2-min
        call ansol_prepareEval (rspaceConstr%p_rumin2,rcollection,"SOL",dcontrolConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 2
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorumin%RvectorBlock(2),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        ! -----
        ! U2-max
        call ansol_prepareEval (rspaceConstr%p_rumax2,rcollection,"SOL",dcontrolConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 2
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorumin%RvectorBlock(2),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        call collct_done (rcollection)
        
      end select
    end select
    
    ! State constraints.
    rspaceConstr%dymin1 = rspaceTimeConstr%dymin1
    rspaceConstr%dymax1 = rspaceTimeConstr%dymax1
    rspaceConstr%dymin2 = rspaceTimeConstr%dymin2
    rspaceConstr%dymax2 = rspaceTimeConstr%dymax2
    rspaceConstr%cstateConstraints = &
        rspaceTimeConstr%cstateConstraints

    rspaceConstr%p_rymin1 => rspaceTimeConstr%p_rymin1
    rspaceConstr%p_rymax1 => rspaceTimeConstr%p_rymax1
    rspaceConstr%p_rymin2 => rspaceTimeConstr%p_rymin2
    rspaceConstr%p_rymax2 => rspaceTimeConstr%p_rymax2
    rspaceConstr%cstateConstraintsType = &
        rspaceTimeConstr%cstateConstraintsType

    rspaceConstr%dstateConstrReg = rspaceTimeConstr%dstateConstrReg
    rspaceConstr%dstateConstrTime = dstateConstrTime
    
    select case (rspaceConstr%cstateConstraints)
    case (1:)
      select case (rspaceConstr%cstateConstraintsType)
      case (1)
        if (.not. associated(rspaceConstr%p_rvectorymin)) then
          allocate(rspaceConstr%p_rvectorymin)
          call lsysbl_createVectorBlock(rspaceDiscrPrimal,rspaceConstr%p_rvectorymin,.true.)
        end if

        if (.not. associated(rspaceConstr%p_rvectorymax)) then
          allocate(rspaceConstr%p_rvectorymax)
          call lsysbl_createVectorBlock(rspaceDiscrPrimal,rspaceConstr%p_rvectorymax,.true.)
        end if
      
        call collct_init (rcollection)

        ! -----
        ! U1-min
        call ansol_prepareEval (rspaceConstr%p_rymin1,rcollection,"SOL",dstateConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 1
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorymin%RvectorBlock(1),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        ! -----
        ! U1-max
        call ansol_prepareEval (rspaceConstr%p_rymax1,rcollection,"SOL",dstateConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 1
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorymax%RvectorBlock(1),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        ! -----
        ! U2-min
        call ansol_prepareEval (rspaceConstr%p_rymin2,rcollection,"SOL",dstateConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 2
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorymin%RvectorBlock(2),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        ! -----
        ! U2-max
        call ansol_prepareEval (rspaceConstr%p_rymax2,rcollection,"SOL",dstateConstrTime)

        ! Set current dimension for the callback routine
        rcollection%IquickAccess(1) = 2
      
        ! Prepare a projection to a vector and evaluate.
        call anprj_discrDirect (rspaceConstr%p_rvectorymin%RvectorBlock(2),&
            ffunctionConstraint, rcollection)

        call ansol_doneEval (rcollection,"SOL")

        call collct_done (rcollection)
        
      end select
    end select

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_doneSpaceConstraints (rspaceConstr)

!<description>
  ! Releases memory allocated by stlin_initSpaceAssembly.
!</description>

!<inputoutput>
  ! Space constraints structure.
  type(t_optcconstraintsSpace), intent(inout) :: rspaceConstr
!</inputoutput>

!</subroutine>

    if (associated(rspaceConstr%p_rvectorumin)) then
      
      ! Analytical constraints active. Release the discrete constraints.
      call lsysbl_releaseVector(rspaceConstr%p_rvectorumin)
      call lsysbl_releaseVector(rspaceConstr%p_rvectorumax)
      
      deallocate(rspaceConstr%p_rvectorumin)
      deallocate(rspaceConstr%p_rvectorumax)
    end if

    if (associated(rspaceConstr%p_rvectorymin)) then
      
      ! Analytical constraints active. Release the discrete constraints.
      call lsysbl_releaseVector(rspaceConstr%p_rvectorymin)
      call lsysbl_releaseVector(rspaceConstr%p_rvectorymax)
      
      deallocate(rspaceConstr%p_rvectorymin)
      deallocate(rspaceConstr%p_rvectorymax)
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_initSpaceAssembly (rspaceTimeDiscr,rdebugFlags,&
      dcontrolConstrTime,dstateConstrTime,rspaceDiscr)

!<description>
  ! Creates a space-assembly data structure from the space-time assembly
  ! data structure. If analytical constraints are active, the routine
  ! invokes an L2 projection to discretise the comnstraints!
  ! The routine can be called multiple times; every new call updates the
  ! data in rspaceDiscr.
!</description>

!<input>
  ! Space-time assembly structure.
  type(t_spaceTimeMatrixDiscrData), intent(in) :: rspaceTimeDiscr
  
  ! Current time where the control constraints should be applied
  real(DP), intent(in) :: dcontrolConstrTime

  ! Current time where the state constraints should be applied
  real(DP), intent(in) :: dstateConstrTime
  
  ! Debug flags
  type(t_optcDebugFlags), target :: rdebugFlags
!</input>

!<output>
  ! Space-assembly structure.
  type(t_spatialMatrixDiscrData), intent(out) :: rspaceDiscr
!</output>

!</subroutine>

    ! Transfer all possible parameters.
    rspaceDiscr%rphysicsPrimal = rspaceTimeDiscr%p_rphysicsPrimal
    rspaceDiscr%rstabilPrimal = rspaceTimeDiscr%p_rstabilPrimal
    rspaceDiscr%rstabilDual = rspaceTimeDiscr%p_rstabilDual
    
    rspaceDiscr%p_rdiscrPrimal => rspaceTimeDiscr%p_rdiscrPrimal
    rspaceDiscr%p_rdiscrPrimalDual => rspaceTimeDiscr%p_rdiscrPrimalDual
    
    rspaceDiscr%p_rstaticAsmTemplates => rspaceTimeDiscr%p_rstaticSpaceAsmTempl
    rspaceDiscr%p_rstaticAsmTemplatesOptC => rspaceTimeDiscr%p_rstaticSpaceAsmTemplOptC

    rspaceDiscr%p_rdebugFlags => rdebugFlags
    rspaceDiscr%p_DobservationArea => rspaceTimeDiscr%p_DobservationArea
    
    ! Initialise the constraints.
    call stlin_initSpaceConstraints (rspaceTimeDiscr%p_rconstraints,&
        dcontrolConstrTime,dstateConstrTime,&
        rspaceTimeDiscr%p_rdiscrPrimal,rspaceDiscr%rconstraints)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stlin_initSpaceAssemblyFromGl (&
      rsettings, ilevel, dcontrolConstrTime, dstateConstrTime, rdiscrData, rphysics)
  
!<description>
  ! Fetches all discretisation data from the main program structure that is
  ! necessary to set up nonlinear matrices at a certain point in time.
!</description>

!<input>
  ! The structure of the main solver
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Level where the discretisation takes place.
  integer, intent(in) :: ilevel

  ! Current time where the control constraints should be applied
  real(DP), intent(in) :: dcontrolConstrTime

  ! Current time where the state constraints should be applied
  real(DP), intent(in) :: dstateConstrTime

  ! OPTIONAL: Alternative physics definition to use.
  ! If not present, the standard global physics settings are used.
  type(t_settings_physics), intent(in), optional :: rphysics
!</input>

!<output>
  ! Discretisation related data, can be used to generate nonlinear matrices.
  type(t_spatialMatrixDiscrData), intent(out) :: rdiscrData
!</output>

!</subroutine>

    ! Get level-independent data.
    rdiscrData%rphysicsPrimal = rsettings%rphysicsPrimal
    if (present(rphysics)) then
      rdiscrData%rphysicsPrimal = rphysics
    end if
    rdiscrData%rsettingsSpaceDiscr = rsettings%rsettingsSpaceDiscr
    rdiscrData%rstabilPrimal = rsettings%rstabilPrimal
    rdiscrData%rstabilDual = rsettings%rstabilDual

    rdiscrData%p_rdiscrPrimal => &
        rsettings%rfeHierPrimal%p_rfeSpaces(ilevel)%p_rdiscretisation

    rdiscrData%p_rdiscrPrimalDual => &
        rsettings%rfeHierPrimalDual%p_rfeSpaces(ilevel)%p_rdiscretisation
    
    rdiscrData%p_rstaticAsmTemplates => &
        rsettings%rspaceAsmHierarchy%p_RasmTemplList(ilevel)

    rdiscrData%p_rstaticAsmTemplatesOptC => &
        rsettings%rspaceAsmHierarchyOptC%p_RasmTemplList(ilevel)

    rdiscrData%p_rdebugFlags => rsettings%rdebugFlags
    
    rdiscrData%p_DobservationArea => rsettings%rsettingsOptControl%p_DobservationArea

    ! Initialise the constraints.
    call stlin_initSpaceConstraints (rsettings%rsettingsOptControl%rconstraints,&
        dcontrolConstrTime, dstateConstrTime,rdiscrData%p_rdiscrPrimal,rdiscrData%rconstraints)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_doneSpaceAssembly (rspaceDiscr)

!<description>
  ! Releases memory allocated by stlin_initSpaceAssembly.
!</description>

!<inputoutput>
  ! Space-assembly structure.
  type(t_spatialMatrixDiscrData), intent(out) :: rspaceDiscr
!</inputoutput>

!</subroutine>

    call stlin_doneSpaceConstraints (rspaceDiscr%rconstraints)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_initSpaceTimeMatrix (&
      rspaceTimeMatrix,cmatrixType,rdiscrData,rsolution,&
      rneumannBoundary,rdirichletBCCBoundary,rglobalData,rdebugFlags)

!<description>
  ! Initialises a space-time matrix.
!</description>

!<input>
  ! Type of the matrix. One of the MATT_xxxx constants.
  integer, intent(in) :: cmatrixType
  
  ! Structure containing all discretisation related data.
  type(t_spaceTimeMatrixDiscrData), intent(in) :: rdiscrData

  ! Space-time solution vector that specifies the evaluation
  ! point of the nonlinearity.
  type(t_spacetimeVector), intent(in), target :: rsolution
  
  ! Definition of the Neumann boundary conditions for all timesteps.
  ! Used for nonlinear boundary conditions.
  type(t_sptiNeumannBoundary), intent(in), target :: rneumannBoundary

  ! Definition of the Dirichlet boundary control boundary conditions for all timesteps.
  ! Used for nonlinear boundary conditions.
  type(t_sptiDirichletBCCBoundary), intent(in), target :: rdirichletBCCBoundary

  ! Global data for callback routines.
  type(t_globalData), intent(in), target :: rglobalData

  ! All debug flags used by the application
  type(t_optcDebugFlags), intent(in), target :: rdebugFlags
!</input>

!<output>
  ! The space-time matrix to initialise.
  type(t_ccoptSpaceTimeMatrix), intent(out) :: rspaceTimeMatrix
!</output>

!</subroutine>

    ! Fetch all necessary information.
    rspaceTimeMatrix%cmatrixType = cmatrixType
    
    rspaceTimeMatrix%NEQtime = rdiscrData%p_rtimeDiscr%nintervals+1
    
!    if (present(rsolution)) then
!      rspaceTimeMatrix%p_rsolution => rsolution
!    else
!      nullify(rspaceTimeMatrix%p_rsolution)
!    end if
    
    ! Remember data about how to discretise.
    rspaceTimeMatrix%rdiscrData = rdiscrData
    rspaceTimeMatrix%p_rsolution => rsolution
    rspaceTimeMatrix%p_rneumannBoundary => rneumannBoundary
    rspaceTimeMatrix%p_rdirichletBCCBoundary => rdirichletBCCBoundary
    rspaceTimeMatrix%p_rdebugFlags => rdebugFlags
    rspaceTimeMatrix%p_rglobalData => rglobalData

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_releaseSpaceTimeMatrix (rspaceTimeMatrix)

!<description>
  ! Cleans up a space-time matrix.
!</description>

!<inputoutput>
  ! The space-time matrix to clean up.
  type(t_ccoptSpaceTimeMatrix), intent(inout) :: rspaceTimeMatrix
!</inputoutput>

!</subroutine>

    ! Nothing to do at the moment. Just basic clean up.
    rspaceTimeMatrix%cmatrixType = MATT_OPTCONTROL
    rspaceTimeMatrix%NEQtime = 0
    
    nullify(rspaceTimeMatrix%p_rneumannBoundary)
    nullify(rspaceTimeMatrix%p_rsolution)
    nullify(rspaceTimeMatrix%p_rdebugFlags)
    nullify(rspaceTimeMatrix%p_rglobalData)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_setupMatrixWeights (rspaceTimeMatrix,&
      ieqTime,irelpos,rnonlinearSpatialMatrix)

!<description>
  ! This routine sets up the matrix weights in the rnonlinearSpatialMatrix structure
  ! according to the position of the corresponding submatrix.
  ! ieqTime defines the timestep which is to be tackled -- which corresponds
  ! to the row in the supermatrix.
  ! irelpos specifies the column in the supermatrix relative to the diagonal --
  ! thus =0 means the diagonal, =1 the submatrix above the diagonal and =-1
  ! the submatrix below the diagonal.
  !
  ! The matrix weights in rnonlinearSpatialMatrix are initialised based on this
  ! information such that the 'assembleMatrix' and 'assembleDefect' routines of
  ! the core equation will build the system matrix / defect at that position
  ! in the supermatrix.
  !
  ! The routine does not initialise the pointers to the basic matrices/discretisation
  ! structures in the rnonlinearSpatialMatrix structure. This has to be done by
  ! the caller!
!</description>

!<input>
  ! A t_ccoptSpaceTimeMatrix structure defining the discretisation of the
  ! coupled space-time matrix.
  type(t_ccoptSpaceTimeMatrix), intent(IN), target :: rspaceTimeMatrix

  ! Substep in the time-dependent simulation = row in the supermatrix.
  ! Range 1..neqTime
  integer, intent(IN) :: ieqTime
  
  ! Specifies the column in the supermatrix relative to the diagonal.
  ! =0: set matrix weights for the diagonal.
  integer, intent(IN) :: irelpos
!</input>

!<inputoutput>
  ! A t_nonlinearSpatialMatrix structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6
  ! block matrix system are initialised depending on the position
  ! specified by ieqTime and neqTime.
  !
  ! The structure must have been initialised with smva_initNonlinMatrix!
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>

!</subroutine>

!    ! If the following constant is set from 1.0 to 0.0, the primal system is
!    ! decoupled from the dual system!
!    real(DP), parameter :: dprimalDualCoupling = 1.0_DP
!
!    ! If the following constant is set from 1.0 to 0.0, the dual system is
!    ! decoupled from the primal system!
!    real(DP), parameter :: ddualPrimalCoupling = 1.0_DP
!
!    ! If the following parameter is set from 1.0 to 0.0, the terminal
!    ! condition between the primal and dual equation is decoupled, i.e.
!    ! the dual equation gets independent from the primal one.
!    real(DP), parameter :: dterminalCondDecoupled = 1.0_DP
!
!    ! If the following parameter is set from 1.0 to 0.0, the time coupling
!    ! is disabled, resulting in a stationary simulation in every timestep.
!    real(DP), parameter :: dtimeCoupling = 1.0_DP

    real(DP) :: dprimalDualCoupling, ddualPrimalCoupling
    real(DP) :: dterminalCondDecoupled, dtimeCoupling

    ! This constant defines the type of equation. There are two equivalent
    ! formulations of the dual equation which only differs in the sign
    ! of the dual velocity.
    ! A constant of "1" here means to use the formulation with "z-y" in
    ! the RHS of the dual equation, while a constant of "-1" means to use the
    ! formulation with "-(z-y)" there.
    ! Note: If this is changed, a "-" sign must be implemented / removed from
    ! the RHS, too!
    real(DP) :: dequationType
    
    ! other local variables
    real(DP) :: dnewton
    real(DP) :: dtstep,dtheta,dtmp,dtimePrimal,dtimeDual
    integer :: ipressureFullyImplicit,cthetaschemetype
    logical :: bconvectionExplicit
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    integer :: neqtime
    
    dprimalDualCoupling = rspaceTimeMatrix%p_rdebugFlags%dprimalDualCoupling
    ddualPrimalCoupling = rspaceTimeMatrix%p_rdebugFlags%ddualPrimalCoupling
    dterminalCondDecoupled = rspaceTimeMatrix%p_rdebugFlags%dterminalCondDecoupled
    dtimeCoupling = rspaceTimeMatrix%p_rdebugFlags%dtimeCoupling
    ipressureFullyImplicit = rspaceTimeMatrix%p_rdebugFlags%ipressureFullyImplicit
    
    dequationType = 1.0_DP
    if (rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%ispaceTimeFormulation .ne. 0) &
      dequationType = -1.0_DP
      
    ! Treat the convection explicitely?
    bconvectionExplicit = &
        rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%iconvectionExplicit .ne. 0

    ! What's the matrix type we should set up? If we have to set up a
    ! Newton matrix, we put dnewton to 1.0 so the Newton part in the
    ! primal velocity is assembled.
    dnewton = 0.0_DP
    if ((rspaceTimeMatrix%cmatrixType .eq. MATT_LINOPTCONTROL) .or. &
        (rspaceTimeMatrix%cmatrixType .eq. MATT_LINPRIMAL) .or. &
        (rspaceTimeMatrix%cmatrixType .eq. MATT_DUAL) .or. &
        (rspaceTimeMatrix%cmatrixType .eq. MATT_LINDUAL)) then
      if (rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation .eq. 0) then
        ! Newton is only to be assembled in Navier-Stokes!
        dnewton = 1.0_DP
      end if
    end if
 
    ! Get the corresponding time level
    !p_rtimeDiscr => rspaceTimeMatrix%rdiscrData%p_rspaceTimeDiscr%p_rtimeDiscr
    p_rtimeDiscr => rspaceTimeMatrix%rdiscrData%p_rtimeDiscr
   
    dtheta = p_rtimeDiscr%dtheta
    neqTime = tdiscr_igetNDofGlob(p_rtimeDiscr)
    
    ! The timestep length is directly taken from the structure.
    dtstep = p_rtimeDiscr%dtstep

    ! Get timestep and current time
    call tdiscr_getTimestep(rspaceTimeMatrix%rdiscrData%p_rtimeDiscr,ieqTime+irelpos,&
        dtimePrimal)
    dtimeDual = dtimePrimal - (1.0_DP-dtheta)*dtstep

    ! The 'old' one step scheme (indicated by ITAG=1) sets up a time discretisation
    ! with all solutions located in the points in time, not inbetween.
    ! This does not correspond to any known analytical minimisation problem.
    cthetaschemetype = p_rtimeDiscr%itag
   
    ! Clear the coefficients
    call smva_clearMatrix (rnonlinearSpatialMatrix)

    select case (rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D

      ! Specify which vectors should be used for evaluating the nonlinearity.
      ! Aii-1 = Aii+1(x_i-1,lambda_i-1).
      ! Aii = Aii(x_i,lambda_i).
      ! Aii+1 = Aii+1(x_i,lambda_i+1)

      if (irelpos .eq. -1) then
      
        ! Matrix on the left of the diagonal.

        rnonlinearSpatialMatrix%iprimalSol = 1
        rnonlinearSpatialMatrix%idualSol = 1

        ! Switch off any stabilisation
        rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal%dupsam = 0.0_DP
        rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam = 0.0_DP

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rnonlinearSpatialMatrix%Dmass(1,1) = dtimeCoupling * (-1.0_DP)/dtstep
        
        rnonlinearSpatialMatrix%Dstokes(1,1) = dtimeCoupling * (1.0_DP-dtheta)
        
        if (ipressureFullyImplicit .eq. 0) then
          rnonlinearSpatialMatrix%DBmat(1,1) = dtimeCoupling * (1.0_DP-dtheta)
        end if

        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dygrad(1,1) = dtimeCoupling * &
              (1.0_DP-dtheta) * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
          
          rnonlinearSpatialMatrix%Dgrady(1,1) = dtimeCoupling * &
              (1.0_DP-dtheta) * dnewton

        else
        
          rnonlinearSpatialMatrix%Dygrad(1,1) = dtimeCoupling * &
              dtheta * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
          
          rnonlinearSpatialMatrix%Dgrady(1,1) = dtimeCoupling * &
              dtheta * dnewton
        
        end if

        if (cthetaschemetype .eq. 0) then
          ! Activate the old timestep scheme. The coupling to the dual must be
          ! weighted by 1-dtheta.
          !
          ! Switch off if alpha <= 0.
          if (rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dalphaC .gt. 0.0_DP) then
            rnonlinearSpatialMatrix%Dmass(1,2) = dtimeCoupling * dprimalDualCoupling * &
                (-dequationType) * (1.0_DP-dtheta) &
                            / rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dalphaC
          end if
        end if

      else if (irelpos .eq. 0) then

        ! The diagonal matrix.
        
        if (ieqTime .lt. neqTime) then
          ! Not the last row
          if (.not. bconvectionExplicit) then
            rnonlinearSpatialMatrix%iprimalSol = 2
            rnonlinearSpatialMatrix%idualSol = 2
            rnonlinearSpatialMatrix%idualSol2 = 3
          else
            rnonlinearSpatialMatrix%iprimalSol = 2
            rnonlinearSpatialMatrix%idualSol = 3
            rnonlinearSpatialMatrix%idualSol2 = 2
          end if
        else
          ! Last row.
          rnonlinearSpatialMatrix%iprimalSol = 2
          rnonlinearSpatialMatrix%idualSol = 2
          rnonlinearSpatialMatrix%idualSol2 = 2
        end if

        rnonlinearSpatialMatrix%Dmass(1,1) = dtimeCoupling * 1.0_DP/dtstep

        if (ieqTime .lt. neqTime) then
          ! Not the last row. Probably deactivate the time coupling.
          rnonlinearSpatialMatrix%Dmass(2,2) = dtimeCoupling * 1.0_DP/dtstep
        else
          ! Last row. The mass matrix must be there.
          rnonlinearSpatialMatrix%Dmass(2,2) = 1.0_DP/dtstep
        end if
        
        rnonlinearSpatialMatrix%Dstokes(1,1) = dtheta
        rnonlinearSpatialMatrix%Dstokes(2,2) = dtheta

        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dygrad(1,1) = &
              dtheta * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
          rnonlinearSpatialMatrix%DygradAdj(2,2) = &
              dtheta * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
          
            rnonlinearSpatialMatrix%Dgrady(1,1) = dtheta * dnewton
            rnonlinearSpatialMatrix%DgradyT(2,2) = &
                  dtheta * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
                
        end if

        rnonlinearSpatialMatrix%DBmat(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%DBmat(2,2) = 1.0_DP

        if (ipressureFullyImplicit .eq. 0) then
          rnonlinearSpatialMatrix%DBmat(1,1) = dtheta
          rnonlinearSpatialMatrix%DBmat(2,2) = dtheta
        end if
        
        rnonlinearSpatialMatrix%DBTmat(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%DBTmat(2,2) = 1.0_DP
        
        ! No coupling to the dual in the first timestep.
        if (ieqTime .gt. 1) then
          ! Switch off if alpha <= 0.
          if (rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dalphaC .gt. 0.0_DP) then
            rnonlinearSpatialMatrix%Dmass(1,2) = dprimalDualCoupling * &
                (-dequationType) * 1.0_DP &
                              / rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dalphaC
          end if
        end if
        
        ! Dirichlet boundary control

        ! Multiply all weigths by the penalty parameter.
        dtmp = 1.0_DP !rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%ddirichletBCPenalty
        rnonlinearSpatialMatrix%DdirichletBCCY(1,1) = dtmp
        
        if (rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dbetaC .gt. 0.0_DP) then
          ! In the first timestep, the dual equation is not connected to the primal one.
          ! Initial condition...
          if (ieqTime .gt. 1) then
            rnonlinearSpatialMatrix%DdirichletBCCLambda(1,2) = &
                (-dequationType) * (-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%dnuConst * &
                dtmp/rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dbetaC)
            rnonlinearSpatialMatrix%DdirichletBCCXi(1,2) = &
                (-dequationType) * (dtmp/rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dbetaC)
          end if
        end if

        select case (cthetaschemetype)
        case (0)
          ! Old 1-step scheme, solutions located at the time interval endpoints
          if (ieqTime .gt. 1) then
            ! Activate the old timestep scheme. The coupling to the dual must be
            ! weighted by dtheta.
            !
            ! Switch off if alpha <= 0.
            if (rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dalphaC .gt. 0.0_DP) then
              rnonlinearSpatialMatrix%Dmass(1,2) = dprimalDualCoupling * &
                  (-dequationType) * dtheta &
                              / rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dalphaC
            end if
          end if

          if (ieqTime .lt. neqTime) then
            ! Not the last row. Standard coupling between primal and dual
            rnonlinearSpatialMatrix%Dmass(2,1) = ddualPrimalCoupling * &
                ( dequationType) * dtheta
          else
            ! Last row. Realise the terminal condition.
            rnonlinearSpatialMatrix%Dmass(2,1) = dterminalCondDecoupled * &
                ( dequationType) * &
                (1.0_DP + rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dgammaC / dtstep)
          end if
          
        case (1)
          ! New 1-step scheme, solutions located inbetween the timesteps
          if (ieqTime .eq. 1) then
            rnonlinearSpatialMatrix%Dmass(2,1) = ddualPrimalCoupling * &
                ( dequationType) * (1.0_DP - dtheta)
                
          else if (ieqTime .lt. neqTime-1) then
            ! Row inbetween. Standard coupling between primal and dual
            rnonlinearSpatialMatrix%Dmass(2,1) = ddualPrimalCoupling * &
                ( dequationType) * 1.0_DP
                
          else if (ieqTime .eq. neqTime-1) then
            ! Last but one row. The terminal condition has influence here.
            rnonlinearSpatialMatrix%Dmass(2,1) = ddualPrimalCoupling * &
                ( dequationType) * (1.0_DP+dterminalCondDecoupled*(1.0_DP-dtheta)*&
                rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dgammaC/dtstep)
                
          else
            ! Last row. Realise the terminal condition.
            rnonlinearSpatialMatrix%Dmass(2,1) = dterminalCondDecoupled * &
                ( dequationType) * &
                (dtheta + dtheta*rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dgammaC / dtstep)
          end if
        end select
            
        if (ieqTime .lt. neqTime) then
          
          if (.not. bconvectionExplicit) then

            if (dnewton .ne. 0.0_DP) then
              rnonlinearSpatialMatrix%DygradT(2,1) = ddualPrimalCoupling * dtheta
              rnonlinearSpatialMatrix%DgradyAdj(2,1) = ddualPrimalCoupling * dtheta

              ! For Crank-Nicolson there appears a 2nd reactive term
              ! stemming from the next timestep.

              rnonlinearSpatialMatrix%DygradT2(2,1) = ddualPrimalCoupling * (1.0_DP-dtheta)
              rnonlinearSpatialMatrix%DgradyAdj2(2,1) = ddualPrimalCoupling * (1.0_DP-dtheta)
            end if
            
          else
          
            if (dnewton .ne. 0.0_DP) then
              rnonlinearSpatialMatrix%DygradT(2,1) = ddualPrimalCoupling * dtheta
              rnonlinearSpatialMatrix%DgradyAdj(2,1) = ddualPrimalCoupling * dtheta
            end if
          
          end if
          
        else
          
          if (.not. bconvectionExplicit) then
          
            ! Weight the mass matrix by GAMMA instead of delta(T).
            ! That's the only difference to the implementation above!
            if (dnewton .ne. 0.0_DP) then
              rnonlinearSpatialMatrix%DygradT(2,1) = dterminalCondDecoupled * dtheta
              rnonlinearSpatialMatrix%DgradyAdj(2,1) = dterminalCondDecoupled * dtheta
            end if

          end if
        end if

      else if (irelpos .eq. 1) then
            
        ! Matrix on the right of the diagonal.
        !
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        rnonlinearSpatialMatrix%iprimalSol = 2
        rnonlinearSpatialMatrix%idualSol = 2
        rnonlinearSpatialMatrix%idualSol2 = 3
        
        ! WARNING: For a very strange situation, taking xi here (which is said
        ! to be the correct evaluation point from the theory) does not lead
        ! to quadratic convergence in the Newton. Taking xi+1 does!?!
        ! So we take xi+1, although the theory tells us to take xi!
        ! ...
        ! No, that does not to be right. Commented out since the above works
        ! as well and should be correct due to the theory.
        ! rnonlinearSpatialMatrix%idualSol = 3

        ! Switch off any stabilisation
        rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal%dupsam = 0.0_DP
        rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam = 0.0_DP

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        rnonlinearSpatialMatrix%Dmass(2,2) = dtimeCoupling * (-1.0_DP)/dtstep
        
        rnonlinearSpatialMatrix%Dstokes(2,2) = dtimeCoupling * (1.0_DP-dtheta)
        
        if (ipressureFullyImplicit .eq. 0) then
          rnonlinearSpatialMatrix%DBmat(2,2) = dtimeCoupling * (1.0_DP-dtheta)
        end if
        
        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%DygradAdj(2,2) = dtimeCoupling * &
              (1.0_DP-dtheta) * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
          
          rnonlinearSpatialMatrix%DgradyT(2,2) = dtimeCoupling * &
              (1.0_DP-dtheta) * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
              
        else
        
          rnonlinearSpatialMatrix%DygradAdj(2,2) = dtimeCoupling * &
                dtheta * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
          
          rnonlinearSpatialMatrix%DgradyT(2,2) = dtimeCoupling * &
                dtheta * real(1-rspaceTimeMatrix%rdiscrData%p_rphysicsPrimal%cequation,DP)
        
        end if

        if (cthetaschemetype .eq. 0) then
          ! Activate the old timestep scheme. The coupling to the dual must be
          ! weighted by 1-dtheta.
          rnonlinearSpatialMatrix%Dmass(2,1) = dtimeCoupling * ddualPrimalCoupling * &
              ( dequationType) * (1.0_DP-dtheta)
        end if
            
        if (bconvectionExplicit) then

          ! DON'T KNOW IF THIS IS CORRECT!!!
        
          if (dnewton .ne. 0.0_DP) then
            rnonlinearSpatialMatrix%DygradTAdj(2,1) = dtimeCoupling * ddualPrimalCoupling * &
                (1.0_DP-dtheta)
            rnonlinearSpatialMatrix%Dgrady(2,1) = -dtimeCoupling * ddualPrimalCoupling * &
                (1.0_DP-dtheta)
          end if
          
        end if
            
      end if
      
      ! Initialise/Update the constraints.
      call stlin_initSpaceConstraints (rspaceTimeMatrix%rdiscrData%p_rconstraints,&
          dtimeDual,dtimePrimal,rspaceTimeMatrix%rdiscrData%p_rdiscrPrimal,&
          rnonlinearSpatialMatrix%rdiscrData%rconstraints)

    end select
    
    ! General parameters.
    rnonlinearSpatialMatrix%dalphaC = rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dalphaC
    rnonlinearSpatialMatrix%dbetaC = rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%dbetaC
    rnonlinearSpatialMatrix%ddirichletBCPenalty = &
        rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%ddirichletBCPenalty
    rnonlinearSpatialMatrix%cmatrixType = rspaceTimeMatrix%cmatrixType
    
    ! Probably scale the system entries.
    if (rspaceTimeMatrix%rdiscrData%p_rsettingsOptControl%csystemScaling .ne. 0) then
    
      ! Do not scale identity matrices.
      ! rnonlinearSpatialMatrix%DidentityY(:,:)    = dtstep * rnonlinearSpatialMatrix%DidentityY(:,:)
      
      rnonlinearSpatialMatrix%Dmass(:,:)   = dtstep * rnonlinearSpatialMatrix%Dmass(:,:)
      rnonlinearSpatialMatrix%Dstokes(:,:)   = dtstep * rnonlinearSpatialMatrix%Dstokes(:,:)
      rnonlinearSpatialMatrix%Dygrad(:,:)   = dtstep * rnonlinearSpatialMatrix%Dygrad(:,:)
      rnonlinearSpatialMatrix%DygradAdj(:,:)   = dtstep * rnonlinearSpatialMatrix%DygradAdj(:,:)
      rnonlinearSpatialMatrix%Dgrady(:,:)  = dtstep * rnonlinearSpatialMatrix%Dgrady(:,:)
      rnonlinearSpatialMatrix%DgradyAdj(:,:)  = dtstep * rnonlinearSpatialMatrix%DgradyAdj(:,:)
      rnonlinearSpatialMatrix%DygradT(:,:)  = dtstep * rnonlinearSpatialMatrix%DygradT(:,:)
      rnonlinearSpatialMatrix%DygradTAdj(:,:)  = dtstep * rnonlinearSpatialMatrix%DygradTAdj(:,:)
      rnonlinearSpatialMatrix%Dgrady2(:,:) = dtstep * rnonlinearSpatialMatrix%Dgrady2(:,:)
      rnonlinearSpatialMatrix%DgradyAdj2(:,:) = dtstep * rnonlinearSpatialMatrix%DgradyAdj2(:,:)
      rnonlinearSpatialMatrix%DygradT2(:,:) = dtstep * rnonlinearSpatialMatrix%DygradT2(:,:)
      rnonlinearSpatialMatrix%DygradTAdj2(:,:) = dtstep * rnonlinearSpatialMatrix%DygradTAdj2(:,:)
      rnonlinearSpatialMatrix%DgradyT(:,:) = dtstep * rnonlinearSpatialMatrix%DgradyT(:,:)
      rnonlinearSpatialMatrix%DBmat(:,:)     = dtstep * rnonlinearSpatialMatrix%DBmat(:,:)
      rnonlinearSpatialMatrix%DdirichletBCCY(:,:) = &
          dtstep * rnonlinearSpatialMatrix%DdirichletBCCY(:,:)
      rnonlinearSpatialMatrix%DdirichletBCCLambda(:,:) = &
          dtstep * rnonlinearSpatialMatrix%DdirichletBCCLambda(:,:)
      rnonlinearSpatialMatrix%DdirichletBCCXi(:,:) = &
          dtstep * rnonlinearSpatialMatrix%DdirichletBCCXi(:,:)
      
      ! Do not scale the pressure block.
      ! rnonlinearSpatialMatrix%DidentityP(:,:)   = dtstep * rnonlinearSpatialMatrix%DidentityP(:,:)

      ! Do not scale the divergence equation.
      ! rnonlinearSpatialMatrix%DBTmat(:,:)     = dtstep * rnonlinearSpatialMatrix%DBTmat(:,:)
      
    end if
    
    ! DEBUG!!!
    !rnonlinearSpatialMatrix%DgradyT(:,:)  = 0.0_DP
    

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine stlin_getFullMatrixDummy(rphysics,rnonlinearSpatialMatrix)
  
!<description>
  ! Configures the parameters in rnonlinearSpatialMatrix with dummy
  ! coefficients so that the configuration resembles the matrix structure
  ! in the current situation. This is used during the allocation of memory
  ! to figure out the matrix structure and allocate necessary submatrices.
!</description>

!<input>
  ! Current physics structure of the primal equation
  type(t_settings_physics), intent(in) :: rphysics
!</input>

!<inputoutput>
  ! A t_nonlinearSpatialMatrix structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps. The coefficients in this
  ! structure are either set to 1.0 or 0.0 depending on whether a
  ! submatrix is active or not.
  !
  ! The structure must have been initialised with smva_initNonlinMatrix!
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>

!</subroutine>

    select case (rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D

      rnonlinearSpatialMatrix%Dstokes(1,1) = 1.0_DP   ! A velocity block
      rnonlinearSpatialMatrix%DBmat(1,1) = 1.0_DP     ! A gradient block
      rnonlinearSpatialMatrix%DBTmat(1,1) = 1.0_DP     ! A divergence block
      rnonlinearSpatialMatrix%DidentityP(1,1) = 1.0_DP   ! Pressure block
      ! A Newton block, if we have Navier-Stokes. For the case
      ! that we use Newton
      rnonlinearSpatialMatrix%Dgrady(1,1) = real(1-rphysics%cequation,DP)
      
      rnonlinearSpatialMatrix%Dstokes(2,2) = 1.0_DP   ! A velocity block
      rnonlinearSpatialMatrix%DBmat(2,2) = 1.0_DP     ! A gradient block
      rnonlinearSpatialMatrix%DBTmat(2,2) = 1.0_DP     ! A divergence block
      ! A Newton block, if we have Navier-Stokes
      rnonlinearSpatialMatrix%Dgrady(2,2) = real(1-rphysics%cequation,DP)
      rnonlinearSpatialMatrix%DidentityP(2,2) = 1.0_DP   ! Pressure block
      
      rnonlinearSpatialMatrix%Dmass(1,1) = 1.0_DP
      rnonlinearSpatialMatrix%Dmass(2,1) = 1.0_DP
      rnonlinearSpatialMatrix%Dmass(1,2) = 1.0_DP
      rnonlinearSpatialMatrix%Dmass(2,2) = 1.0_DP
      rnonlinearSpatialMatrix%Dgrady(2,1) = real(1-rphysics%cequation,DP)
      
      ! For Dirichlet boundary control, coupling matrix to xi
      ! in the primal equation.
      rnonlinearSpatialMatrix%DdirichletBCCXi(1,2) = 1.0_DP
      
    end select
      
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_spaceTimeMatVec (rspaceTimeMatrix, rx, rd, cx, cy, &
      cfilter, dnorm, rboundaryConditions, bprintRes,&
      cx2, rx2)

!<description>
  ! This routine performs a matrix-vector multiplication with the
  ! system matrix A defined by rspaceTimeMatrix.
  !    rd  :=  cx A(p_rsolution + cx2 * rx2) rx  +  cy rd
  ! If rspaceTimeMatrix does not specify an evaluation point for th nonlinearity
  ! in A(.), the routine calculates
  !    rd  :=  cx A(rx + cx2 * rx2) rx  +  cy rd
  ! rd is overwritten by the result.
!</description>

!<input>
  ! A t_ccoptSpaceTimeMatrix structure defining the space-time matrix.
  type(t_ccoptSpaceTimeMatrix), intent(inout) :: rspaceTimeMatrix

  ! A space-time vector defining the current solution.
  type(t_spacetimeVector), intent(in) :: rx
  
  ! Type of filter to apply to the vectors. A combination of SPTID_FILTER_xxxx
  ! flags that specifies which type of filter is to be applied to the
  ! output vector during the matriux vector multiplication.
  integer(I32), intent(in) :: cfilter

  ! OPTIONAL: Space-time boundary conditions. Must be specified if
  ! the boundary-condition filter should be applied.
  type(t_optcBDC), intent(in), optional :: rboundaryConditions
  
  ! If set to TRUE and dnorm is present, too, the residuals are printed to the
  ! terminal.
  logical, intent(in), optional :: bprintRes
!</input>

!<inputoutput>
  ! A second space-time vector.
  type(t_spacetimeVector), intent(inout) :: rd
  
  ! Multiplication factor for rx
  real(DP), intent(in) :: cx
  
  ! Multiplication factor for rd
  real(DP), intent(in) :: cy
  
  ! OPTIONAL: Multiplication factor for rx2
  real(DP), intent(in), optional :: cx2
  
  ! OPTIONAL: Correction vector
  type(t_spacetimeVector), intent(in), optional :: rx2
!</inputoutput>

!<output>
  ! OPTIONAL: If specified, returns the $l_2$-norm of rd.
  real(DP), intent(OUT), optional :: dnorm
!<output>

!</subroutine>

    ! local variables
    integer :: ieqTime,icp,neqTime
    logical :: blocalPrintRes
    type(t_vectorBlock) :: rtempVectorD
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    type(t_vectorBlock), dimension(3) :: rtempVector
    type(t_vectorBlock), dimension(3), target :: rtempVectorEval
    type(t_vectorBlock) :: rtempVector2
    type(t_blockDiscretisation), pointer :: p_rdiscr
    real(DP) :: dnormpart
    type(t_matrixBlock) :: rblockTemp
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    real(DP) :: dtstep,dtimePrimal,dtimeDual
    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    type(t_spatialMatrixDiscrData) :: rspaceDiscr
    type(t_spatialMatrixNonlinearData), target :: rnonlinearity
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dx3,p_Db
    real(DP), dimension(:), pointer :: p_DxE1,p_DxE2,p_DxE3
!#ifndef TESTCODE
!    type(t_matrixAssemblyFlags) :: rflags
!    type(t_matrixBlock) :: rmatrix
!#endif
    
    blocalPrintRes = .false.
    if (present(bprintRes)) blocalPrintRes = bprintRes
    
    ! Get the current space- and time discretisation
    p_rdiscr => rspaceTimeMatrix%rdiscrData%p_rspaceDiscr
    p_rtimeDiscr => rspaceTimeMatrix%rdiscrData%p_rtimeDiscr
    
    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)

    ! Create a temp vector that contains the part of rd which is to be modified.
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorD,.false.)
    
    ! The vector will be a defect vector. Assign the boundary conditions so
    ! that we can implement them.
    call lsysbl_assignDiscreteBC(rtempVectorD,rdiscreteBC)
    call lsysbl_assignDiscreteFBC(rtempVectorD,rdiscreteFBC)
    
    ! Create a temp vector for the X-vectors at timestep i-1, i and i+1.
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(1),.false.)
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(2),.false.)
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(3),.false.)
    
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector2,.false.)
    
    ! Now what's the evaluation point where to evaluate the nonlinearity/ies?
    ! If the structure defines an evaluation point, we take that one, so
    ! we need another temp vector that holds the evaluation point.
    ! If not, we take the vector rx. This is archived by creating new
    ! vectors that share their information with rx.
    if (associated(rspaceTimeMatrix%p_rsolution)) then
      call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(1),.false.)
      call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(2),.false.)
      call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(3),.false.)
    else
      call lsysbl_duplicateVector (rtempVector(1),rtempVectorEval(1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsysbl_duplicateVector (rtempVector(2),rtempVectorEval(2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsysbl_duplicateVector (rtempVector(3),rtempVectorEval(3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end if
    
    if (.not. associated(rspaceTimeMatrix%p_rsolution)) then
      
      ! Get the parts of the X-vector which are to be modified at first --
      ! subvector 1, 2 and 3. rtempVector(1) contains the 'previous' solution,
      ! rtempVector(2) the 'current' and rtempVector(3) the 'next' one.
      call sptivec_getTimestepData(rx, 1+0, rtempVector(2))
      
      ! If necessary, include the correction
      if (present(rx2)) then
        call sptivec_getTimestepData(rx2,1+0, rtempVector2)
        if (.not. present(cx2)) then
          call output_line ("Multiplication factor cx2 not specified!", &
              OU_CLASS_ERROR,OU_MODE_STD,"stlin_spaceTimeMatVec")
          call sys_halt()
        end if
        
        call lsysbl_vectorLinearComb (rtempVector2,rtempVector(2),cx2,1.0_DP)
      end if
      
      ! If necesary, multiply the rtempVectorX. We have to take a -1 into
      ! account as the actual matrix multiplication routine smva_assembleDefect
      ! introduces another -1!
      if (cx .ne. -1.0_DP) then
        call lsysbl_scaleVector (rtempVector(2),-cx)
      end if

    else
    
      ! Get the parts of the X-vector which are to be modified at first --
      ! subvector 1, 2 and 3. rtempVector(1) contains the 'previous' solution,
      ! rtempVector(2) the 'current' and rtempVector(3) the 'next' one.
      call sptivec_getTimestepData(rx, 1+0, rtempVector(2))
      
      ! If necesary, multiply the rtempVectorX. We have to take a -1 into
      ! account as the actual matrix multiplication routine smva_assembleDefect
      ! introduces another -1!
      if (cx .ne. -1.0_DP) then
        call lsysbl_scaleVector (rtempVector(2),-cx)
      end if

      ! If a nonlinearity is involved, rtempVectorEval(1) contains the 'previous',
      ! rtempVectorEval(2) the 'current' and rtempVectorEval(3) the 'next'
      ! solution (relative to the current time step) where to evaluate
      ! the nonlinearity.

      call sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, &
          1+0, rtempVectorEval(2))

      ! If necessary, include the correction.
      ! This is only incorporated into the nonlinearity, not into rx.
      if (present(rx2)) then
        call sptivec_getTimestepData(rx2,1+0, rtempVector2)
        if (.not. present(cx2)) then
          call output_line ("Multiplication factor cx2 not specified!", &
              OU_CLASS_ERROR,OU_MODE_STD,"stlin_spaceTimeMatVec")
          call sys_halt()
        end if
        
        call lsysbl_vectorLinearComb (rtempVector2,rtempVectorEval(2),cx2,1.0_DP)
      end if
      
    end if
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVector(1),p_Dx1)
    call lsysbl_getbase_double (rtempVector(2),p_Dx2)
    call lsysbl_getbase_double (rtempVector(3),p_Dx3)
    call lsysbl_getbase_double (rtempVectorD,p_Db)
    call lsysbl_getbase_double (rtempVectorEval(1),p_DxE1)
    call lsysbl_getbase_double (rtempVectorEval(2),p_DxE2)
    call lsysbl_getbase_double (rtempVectorEval(3),p_DxE3)
    
    ! If dnorm is specified, clear it.
    if (present(dnorm)) then
      dnorm = 0.0_DP
    end if
    
    ! Loop through the substeps
    neqTime = tdiscr_igetNDofGlob(p_rtimeDiscr)
    
    do ieqTime = 1,neqTime
            
      ! Get the current time.
      call tdiscr_getTimestep(p_rtimeDiscr,ieqTime-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal - (1.0_DP-p_rtimeDiscr%dtheta)*dtstep
      
      ! Get a space-assembly structure from our space-time assembly structure.
      ! Necessary for assembling matrices.
      call stlin_initSpaceAssembly (rspaceTimeMatrix%rdiscrData,&
          rspaceTimeMatrix%p_rdebugFlags,dtimeDual,dtimePrimal,rspaceDiscr)
      
      ! Get the part of rd which is to be modified.
      if (cy .ne. 0.0_DP) then
        call sptivec_getTimestepData(rd, ieqTime, rtempVectorD)
        
        ! If cy <> 1, multiply rtempVectorD by that.
        if (cy .ne. 1.0_DP) then
          call lsysbl_scaleVector (rtempVectorD,cy)
        end if
      else
        call lsysbl_clearVector (rtempVectorD)
      end if
      
      ! Form a t_spatialMatrixNonlinearData structure that encapsules the nonlinearity
      ! of the spatial matrix.
      call smva_initNonlinearData (rnonlinearity,&
          rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3),&
          dtimePrimal,dtimeDual,&
          rspaceTimeMatrix%p_rneumannBoundary%p_RbdRegion(ieqTime),&
          rspaceTimeMatrix%p_rdirichletBCCBoundary%p_RbdRegion(ieqTime))
      
      if (ieqTime .ne. neqTime) then
      
        ! Read the solution of the 'next' timestep and the 'next' evaluation
        ! point.

        if (.not. associated(rspaceTimeMatrix%p_rsolution)) then
          
          call sptivec_getTimestepData(rx, ieqTime+1, rtempVector(3))

          ! If necessary, include the correction
          if (present(rx2)) then
            call sptivec_getTimestepData(rx2,ieqTime+1, rtempVector2)
            call lsysbl_vectorLinearComb (rtempVector2,rtempVector(3),cx2,1.0_DP)
          end if

          ! If necesary, multiply the rtempVectorX. We have to take a -1 into
          ! account as the actual matrix multiplication routine smva_assembleDefect
          ! introduces another -1!
          if (cx .ne. -1.0_DP) then
            call lsysbl_scaleVector (rtempVector(3),-cx)
          end if
          
        else
        
          ! Nonlinearity specified. The correction is only included into the
          ! nonlinearity.

          call sptivec_getTimestepData(rx, ieqTime+1, rtempVector(3))

          ! If necesary, multiply the rtempVectorX. We have to take a -1 into
          ! account as the actual matrix multiplication routine smva_assembleDefect
          ! introduces another -1!
          if (cx .ne. -1.0_DP) then
            call lsysbl_scaleVector (rtempVector(3),-cx)
          end if

          ! Get the nonlinearity          
          call sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, &
              ieqTime+1, rtempVectorEval(3))

          ! If necessary, include the correction
          if (present(rx2)) then
            call sptivec_getTimestepData(rx2,ieqTime+1, rtempVector2)
            call lsysbl_vectorLinearComb (rtempVector2,rtempVectorEval(3),cx2,1.0_DP)
          end if
          
        end if

      end if

      ! The first and last substep is a little bit special concerning
      ! the matrix!
      if (ieqTime .eq. 1) then
        
        ! We are in the first substep. Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( A11 A12   0   0 ... )  ( x1 )  =  ( f1 )
        !  ( ... ... ... ... ... )  ( .. )     ( .. )
        !
        ! So we have to compute:
        !
        !  d1  :=  f1  -  A11 x1  -  A12 x2
        !
        ! -----
        
        ! The diagonal matrix.
      
        ! Set up the matrix weights of that submatrix.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
        call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,0,rnonlinearSpatialMatrix)
          
        ! Subtract: rd = rd - A11 x1
        call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(2),rtempVectorD,1.0_DP)

        ! -----
      
        ! Create the matrix
        !   A12  :=  -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and subtract A12 x2 from rd.

        ! Set up the matrix weights of that submatrix.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
        call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,1,rnonlinearSpatialMatrix)

        ! Subtract: rd = rd - A12 x2
        call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(3),rtempVectorD,1.0_DP)

        ! Release the block mass matrix.
        call lsysbl_releaseMatrix (rblockTemp)

      else if (ieqTime .lt. neqTime) then

        ! We are sonewhere in the middle of the matrix. There is a substep
        ! ieqTime+1 and a substep ieqTime-1!  Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
        !  ( ... Aii-1 Aii Aii+1 ... )  ( xi )  =  ( fi )
        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
        !
        ! So we have to compute:
        !
        !  dn  :=  fn  -  Aii-1 xi-1  -  Aii xi  -  Aii+1 xi+1
        
        ! -----
        
        ! Create the matrix
        !   Aii-1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the primal velocity.

        ! Set up the matrix weights of that submatrix.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
        call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,-1,rnonlinearSpatialMatrix)
            
        ! Subtract: rd = rd - Aii-1 xi-1.
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(1),rtempVectorD,1.0_DP)

        ! Release the block mass matrix.
        call lsysbl_releaseMatrix (rblockTemp)

        ! -----

        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
        call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,0,rnonlinearSpatialMatrix)

        ! Subtract: rd = rd - Aii xi
        call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(2),rtempVectorD,1.0_DP)
            
        ! -----
        
        ! Create the matrix
        !   Aii+1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
        call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,1,rnonlinearSpatialMatrix)
          
        ! Subtract: rd = rd - Aii+1 xi+1
        call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(3),rtempVectorD,1.0_DP)
        
        ! Release the block mass matrix.
        call lsysbl_releaseMatrix (rblockTemp)

      else
        
        ! We are in the last substep. Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( ... ... ... ...   ... )  ( .. )     ( .. )
        !  ( ... ...   0 Ann-1 Ann )  ( xn )  =  ( fn )
        !
        ! So we have to compute:
        !
        !  dn  :=  fn  -  Ann-1 xn-1  -  Ann xn
        
        ! -----
        
        ! Create the matrix
        !   Ann-1 = -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
        call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,-1,rnonlinearSpatialMatrix)
          
        ! Subtract: rd = rd - Ann-1 xn-1
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xn due to the discretisation scheme!!!
        call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(1),rtempVectorD,1.0_DP)
     
        ! -----
        
        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
            rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
        call stlin_setupMatrixWeights (rspaceTimeMatrix,&
          ieqTime,0,rnonlinearSpatialMatrix)

        ! Subtract: rd = rd - Ann xn
        call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(2),rtempVectorD,1.0_DP)
      
      end if

      ! Implement the BC`s?
      if (iand(cfilter,SPTID_FILTER_BCDEF) .ne. 0) then
      
        if (.not. present(rboundaryConditions)) then
          call output_line ('Boundary conditions not specified!', &
              OU_CLASS_ERROR,OU_MODE_STD,'stlin_spaceTimeMatVec')
          call sys_halt()
        end if

        ! Discretise the boundary conditions at the new point in time.
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call sbc_assembleBDconditions (rboundaryConditions,dtimePrimal,dtimeDual,&
            CCSPACE_PRIMALDUAL,rspaceTimeMatrix%p_rglobalData,SBC_DIRICHLETBC,&
            p_rtimediscr,p_rdiscr,rdiscreteBC)

        ! Implement the boundary conditions into the defect.
        call vecfil_discreteBCdef (rtempVectorD)
        call vecfil_discreteFBCdef (rtempVectorD)
      end if

      if ((ieqTime .eq. 1) .and. (iand(cfilter,SPTID_FILTER_ICDEF) .ne. 0)) then
        ! Implement the initial conditions into the defect.
        call tbc_implementInitCondDefSingle (rtempVectorD)
      end if

      if ((ieqTime .eq. neqTime) .and. &
          (iand(cfilter,SPTID_FILTER_TCDEF) .ne. 0)) then
        ! Implement the initial conditions into the defect.
        call tbc_implementTermCondDefSingle (rtempVectorD)
      end if
      
      ! If dnorm is specified, calculate the norm of the sub-defect vector and
      ! add it to dnorm.
      if (present(dnorm) .or. blocalPrintRes) then
        dnormpart = lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)**2
        if (present(dnorm)) &
          dnorm = dnorm + dnormpart
        
        if (blocalPrintRes) then
          call output_line ('||D_'//trim(sys_siL(ieqTime,10))//'|| = '//&
              trim(sys_sdEL(sqrt(dnormpart),10)) )
          do icp=1,rtempVectorD%nblocks
            call output_line ('  ||D_'//&
                trim(sys_siL(ieqTime,10))//'^'//trim(sys_siL(icp,2))&
                //'|| = '//&
                trim(sys_sdEL(&
                    lsyssc_vectorNorm(rtempVectorD%RvectorBlock(icp),LINALG_NORML2),10)) )
          end do

!#ifndef TESTCODE
!          if (ieqTime .eq. 2) then
!            ! Test
!            call smva_initNonlinMatrix (rnonlinearSpatialMatrix,&
!                rspaceTimeMatrix%p_rglobalData,rspaceDiscr,rnonlinearity)
!            call stlin_setupMatrixWeights (rspaceTimeMatrix,&
!              ieqTime,0,rnonlinearSpatialMatrix)
!
!            ! Subtract: rd = rd - Ann xn
!            call sptivec_getTimestepData(rd, ieqTime, rtempVectorD)
!            call smva_assembleDefect (rnonlinearSpatialMatrix,rtempVector(2),rtempVectorD,1.0_DP)
!            if (iand(cfilter,SPTID_FILTER_BCDEF) .ne. 0) then
!              call vecfil_discreteBCdef (rtempVectorD)
!            end if
!            print *,lsyssc_vectorNorm(rtempVectorD%RvectorBlock(1),LINALG_NORML2)
!
!            call sptivec_getTimestepData(rd, ieqTime, rtempVectorD)
!            call smva_assembleMatrix (CCMASM_ALLOCMEM+CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,rflags,&
!                rnonlinearSpatialMatrix,rmatrix)
!
!            call vecio_writeBlockVectorHR (rtempVectorD, "rb",&
!                .true., 0, "rb2.txt", "(E15.5)")
!            call vecio_writeBlockVectorHR (rtempVector(2), "rx",&
!                .true., 0, "rx2.txt", "(E15.5)")
!
!            call lsysbl_blockmatvec (rmatrix,rtempvector(2),rtempvectord,-1.0_dp,1.0_dp)
!
!            if (iand(cfilter,sptid_filter_bcdef) .ne. 0) then
!              call vecfil_discretebcdef (rtempvectord)
!              call matfil_discretebc (rmatrix,rtempvectord%p_rdiscretebc)
!            end if
!
!            call vecio_writeblockvectorhr (rtempvectord, "rd",&
!                .true., 0, "rd2.txt", "(e15.5)")
!
!            call matio_writeBlockMatrixHR (rmatrix, "matrix",&
!              .true., 0, "matrix2.txt", "(E15.5)", dthreshold=1E-12_DP)
!
!            print *,lsyssc_vectorNorm(rtempVectorD%RvectorBlock(1),LINALG_NORML2)
!
!            call lsysbl_releaseMatrix (rmatrix)
!          end if
!#endif
                  
        end if
      end if
      
      ! Save the defect vector back to rd.
      call sptivec_setTimestepData(rd, ieqTime, rtempVectorD)
            
      ! Cycle the solution vectors and the evaluation points: 1 <- 2 <- 3.
      call lsysbl_copyVector (rtempVector(2),rtempVector(1))
      call lsysbl_copyVector (rtempVector(3),rtempVector(2))

      if (associated(rspaceTimeMatrix%p_rsolution)) then
        call lsysbl_copyVector (rtempVectorEval(2),rtempVectorEval(1))
        call lsysbl_copyVector (rtempVectorEval(3),rtempVectorEval(2))
      end if
            
      ! Now, the 3rd subvector is again free to take the next vector.
      
      ! Release the space assembly structure again.
      call stlin_doneSpaceAssembly (rspaceDiscr)
      
    end do
    
    ! If dnorm is specified, normalise it.
    ! It was calculated from rspaceTimeDiscr%niterations+1 subvectors.
    if (present(dnorm)) then
      dnorm = sqrt(dnorm / real(neqTime,DP))
    end if
    
    ! Release the temp vectors.
    call lsysbl_releaseVector (rtempVectorEval(3))
    call lsysbl_releaseVector (rtempVectorEval(2))
    call lsysbl_releaseVector (rtempVectorEval(1))
    call lsysbl_releaseVector (rtempVector(3))
    call lsysbl_releaseVector (rtempVector(2))
    call lsysbl_releaseVector (rtempVector(1))
    call lsysbl_releaseVector (rtempVectorD)
    call lsysbl_releaseVector (rtempVector2)
    
    ! Release the BC's again.
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
    call bcasm_releaseDiscreteBC(rdiscreteBC)
    
  end subroutine
   
  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_dampDualSolution (rphysics,rx,cx)

!<description>
  ! Applies damping to the dual solution. Multiplies all components of the dual
  ! solution by cx.
!</description>

!<input>
  ! Current physics structure of the primal equation
  type(t_settings_physics), intent(in) :: rphysics

  ! Damping factor
  real(DP), intent(in) :: cx
!</input>

!<inputoutput>
  ! Space-time vector to apply damping to.
  type(t_spacetimeVector), intent(inout) :: rx
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_vectorBlock) :: rvector
    
    if (cx .eq. 1.0_DP) return
    call lsysbl_createVectorBlock (rx%p_rspaceDiscr,rvector)

    select case (rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes
      do i=1,rx%NEQtime
        ! Fetch the vector, scale, write back.
        call sptivec_getTimestepData(rx,i,rvector)
        
        call lsyssc_scaleVector (rvector%RvectorBlock(4),cx)
        call lsyssc_scaleVector (rvector%RvectorBlock(5),cx)
        call lsyssc_scaleVector (rvector%RvectorBlock(6),cx)
        
        call sptivec_setTimestepData(rx,i,rvector)
      end do
    end select
    
    call lsysbl_releaseVector (rvector)

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
!      rtempvector1, rtempvector2, rtempvector3, bimplementBC)
!
!!<description>
!  ! Assembles the space-time RHS vector rb. Bondary conditions are NOT
!  ! implemented!
!  !
!  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
!  ! this routine!
!!</description>
!
!!<input>
!  ! A problem structure that provides information on all
!  ! levels as well as temporary vectors.
!  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!
!  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
!  ! coupled space-time system.
!  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!!</input>
!
!!<inputoutput>
!  ! A temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector1
!
!  ! A second temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector2
!
!  ! A third temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector3
!
!  ! A space-time vector that receives the RHS.
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
!
!  ! Whether to implement boundary conditions into the RHS or not.
!  LOGICAL, INTENT(IN) :: bimplementBC
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    INTEGER :: isubstep
!    REAL(DP) :: dtheta
!    TYPE(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
!
!    ! A temporary vector for the creation of the RHS.
!    TYPE(t_vectorBlock) :: rtempVectorRHS
!
!    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd, p_Drhs
!
!    ! Theta-scheme identifier.
!    ! =1: impliciz Euler.
!    ! =0.5: Crank Nicolson
!    dtheta = rproblem%rtimedependence%dtimeStepTheta
!
!    ! ----------------------------------------------------------------------
!    ! Generate the global RHS vector
!
!    CALL lsysbl_getbase_double (rtempVector1,p_Dx)
!    CALL lsysbl_getbase_double (rtempVector2,p_Db)
!    CALL lsysbl_getbase_double (rtempVector3,p_Dd)
!
!    ! Assemble 1st RHS vector in X temp vector.
!    CALL generateRHS (rproblem,0,rb%ntimesteps,rspaceTimeDiscr%rtimeDiscr%dtstep,&
!        rtempVector1, .TRUE., .FALSE.)
!
!    ! Assemble the 2nd RHS vector in the RHS temp vector
!    CALL generateRHS (rproblem,1,rb%ntimesteps,rspaceTimeDiscr%rtimeDiscr%dtstep,&
!        rtempVector2, .TRUE., .FALSE.)
!
!    ! Assemble the 3rd RHS vector in the defect temp vector
!    IF (rspaceTimeDiscr%rtimeDiscr%nintervals .GE. 2) THEN
!      CALL generateRHS (rproblem,2,rb%ntimesteps,rspaceTimeDiscr%rtimeDiscr%dtstep,&
!          rtempVector3, .TRUE., .FALSE.)
!    ELSE
!      CALL lsysbl_copyVector (rtempVector2,rtempVector3)
!    END IF
!
!    ! Create a copy of the X temp vector (RHS0). That vector will be
!    ! our destination vector for assembling the RHS in all timesteps.
!    CALL lsysbl_copyVector (rtempVector1,rtempVectorRHS)
!
!    ! DEBUG!!!
!    CALL lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
!
!    ! RHS 0,1,2 -> 1-2-3
!
!    DO isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals
!
!      IF (isubstep .EQ. 0) THEN
!
!        ! Primal RHS comes from rtempVector1. The dual from the
!        ! isubstep+1'th RHS in rtempVector2.
!        !
!        ! primal RHS(0) = PRIMALRHS(0)
!        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
!
!        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(1),rtempVectorRHS%RvectorBlock(1))
!        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(2),rtempVectorRHS%RvectorBlock(2))
!        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(4),rtempVector2%RvectorBlock(4),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(4))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(5),rtempVector2%RvectorBlock(5),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(5))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(6),rtempVector2%RvectorBlock(6),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(6))
!
!        ! In the 0'th timestep, there is no RHS in the dual equation!
!        ! That is because of the initial condition, which fixes the primal solution
!        ! => dual solution has no influence on the primal one
!        ! => setting up a dual RHS in not meaningful as the dual RHS cannot
!        !    influence the primal solution
!        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(4))
!        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(5))
!        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(6))
!
!      ELSE IF (isubstep .LT. rspaceTimeDiscr%rtimeDiscr%nintervals) THEN
!
!        ! We are somewhere 'in the middle'.
!        !
!        ! Dual RHS comes from rtempVector3. The primal from the
!        ! isubstep-1'th RHS.
!        !
!        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(1),rtempVector2%RvectorBlock(1),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(1))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(2),rtempVector2%RvectorBlock(2),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(2))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(3),rtempVector2%RvectorBlock(3),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(3))
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(4),rtempVector3%RvectorBlock(4),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(4))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(5),rtempVector3%RvectorBlock(5),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(5))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(6),rtempVector3%RvectorBlock(6),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(6))
!
!        IF (isubstep .LT. rspaceTimeDiscr%rtimeDiscr%nintervals-1) THEN
!          ! Shift the RHS vectors and generate the RHS for the next time step.
!          ! (Yes, I know, this could probably be solved more elegant without copying anything
!          ! using a ring buffer ^^)
!          CALL lsysbl_copyVector(rtempVector2,rtempVector1)
!          CALL lsysbl_copyVector(rtempVector3,rtempVector2)
!          CALL generateRHS (rproblem,isubstep+2,&
!              rspaceTimeDiscr%rtimeDiscr%nintervals,&
!              rspaceTimeDiscr%rtimeDiscr%dtstep,rtempVector3, .TRUE., .FALSE.)
!        END IF
!
!      ELSE
!
!        ! We are 'at the end'.
!        !
!        ! Dual RHS comes from rtempVector3. The primal from the
!        ! isubstep-1'th RHS and rtempVector3.
!        !
!        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!        ! dual RHS(0)   = DUALRHS(0)
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(1),rtempVector3%RvectorBlock(1),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(1))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(2),rtempVector3%RvectorBlock(2),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(2))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(3),rtempVector3%RvectorBlock(3),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(3))
!
!        !CALL generateRHS (rproblem,isubstep+1,rspaceTimeDiscr%niterations,&
!        !    rtempVector3, .TRUE., .FALSE.)
!
!        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
!        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
!        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(6),rtempVectorRHS%RvectorBlock(6))
!
!        ! Multiply the last RHS of the dual equation -z by gamma, that's it.
!        CALL lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),rspaceTimeDiscr%dgammaC)
!        CALL lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),rspaceTimeDiscr%dgammaC)
!
!      END IF
!
!      ! Implement the boundary conditions into the RHS vector
!      IF (bimplementBC) THEN
!        CALL generateRHS (rproblem,isubstep,&
!            rspaceTimeDiscr%rtimeDiscr%nintervals,&
!            rspaceTimeDiscr%rtimeDiscr%dtstep, rtempVectorRHS, .FALSE., .TRUE.)
!      END IF
!
!      ! Save the RHS.
!      CALL sptivec_setTimestepData(rb, isubstep, rtempVectorRHS)
!
!    END DO
!
!    ! Release the temp vector for generating the RHS.
!    CALL lsysbl_releaseVector (rtempVectorRHS)
!
!  CONTAINS
!
!    SUBROUTINE generateRHS (rproblem,isubstep,nsubsteps,dtstep,&
!        rvector, bgenerate, bincludeBC)
!
!    ! Generate the RHS vector of timestep isubstep and/or include boundary
!    ! conditions.
!
!    ! Problem structure.
!    TYPE(t_problem), INTENT(INOUT) :: rproblem
!
!    ! Number of the substep where to generate the RHS vector
!    INTEGER, INTENT(IN) :: isubstep
!
!    ! Total number of substeps
!    INTEGER, INTENT(IN) :: nsubsteps
!
!    ! Length od one timestep
!    REAL(DP), INTENT(IN) :: dtstep
!
!    ! Destination vector
!    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!
!    ! Whether to generate the RHS vector or not
!    LOGICAL, INTENT(IN) :: bgenerate
!
!    ! Whether to include boundary conditions
!    LOGICAL, INTENT(IN) :: bincludeBC
!
!    ! DEBUG!!!
!    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!
!      ! DEBUG!!!
!      CALL lsysbl_getbase_double (rvector,p_Ddata)
!
!      ! Set the time where we are at the moment
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*dtstep
!      rproblem%rtimedependence%itimestep = isubstep
!
!      ! Assemble the RHS?
!      IF (bgenerate) THEN
!
!        ! Generate the RHS of that point in time into the vector.
!        CALL cc_generateBasicRHS (rproblem,rvector)
!
!      END IF
!
!      ! Include BC's?
!      IF (bincludeBC) THEN
!
!        ! Initialise the collection for the assembly process with callback routines.
!        ! Basically, this stores the simulation time in the collection if the
!        ! simulation is nonstationary.
!        CALL user_initCollectForAssembly (rproblem,rproblem%rcollection)
!
!        ! Discretise the boundary conditions at the new point in time --
!        ! if the boundary conditions are nonconstant in time!
!        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!          CALL cc_updateDiscreteBC (rproblem)
!        END IF
!
!        ! Implement the boundary conditions into the RHS.
!        ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
!        ! as Dirichlet values mustn't be multiplied with GAMMA!
!        CALL vecfil_discreteBCsol (rvector)
!        CALL vecfil_discreteFBCsol (rvector)
!
!        ! Clean up the collection (as we are done with the assembly, that's it.
!        CALL user_doneCollectForAssembly (rproblem,rproblem%rcollection)
!
!      END IF
!
!    END SUBROUTINE
!
!  END SUBROUTINE

!betrachte die GMV's. Entkoppelt gibt's 2. Ordnung in der Geschwindigkeit
!(Druck=0), Gekoppelt nicht :(

  
end module
