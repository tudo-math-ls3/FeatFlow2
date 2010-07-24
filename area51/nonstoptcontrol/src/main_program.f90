!##############################################################################
!# ****************************************************************************
!# <name> main_program </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module solves an optimal control problem for the stationary and
!# nonstationary Navier-Stokes optimal control problem 
!#
!#  $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alphga/2||u||^2 $$
!#
!#  $$- \nu Delta(y) + y*\Nabla(y) + \Nabla p = f $$
!#  $$ \Nabla \cdot y = 0$$
!#  $$- \nu Delta(\lambda) - y*\Nabla(\lambda) + \lambda\Nabla y + \Nabla \xi = y-z $$
!#  $$ \Nabla \cdot \lambda = 0$$
!#              
!#
!# on a 2D domain for a 2D function $y=(y_1,y_2)$, a pressure $p$,
!# a dual velocity $\lambda$ and a dual pressure $\xi$. $u$ is the control
!# and $z$ a desired flow field.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# For the nonlinearity, the nonlinear solver is invoked. The
!# defect that is setted up there is preconditioned by a linear Multigrid
!# solver with a simple-VANKA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type. As coarse grid solver,
!# UMFPACK is used.
!# </purpose>
!##############################################################################

module main_program

  use fsystem
  use genoutput
  use storage
  use basicgeometry
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use element
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use fparser
  use statistics
  
  use collection
  use convection
    
  use multilevelprojection
  use externalstorage
  use paramlist
  
  use meshhierarchy
  use fespacehierarchy
  use timescalehierarchy
  
  use physics
  use spacetimevectors
  use spacetimerhs
  use spacetimehierarchy
  use spacetimelinsol
  use spatialoperators
  use postprocessing
  use spacetimeinterlevelprojection

  implicit none

!<globals>

  ! Directory containing the data files.
  character(LEN=SYS_STRLEN), save :: DIR_DATA = "./data";

!</globals>

  ! Main data of the solver
  type t_maindata
    
    ! The domain
    type(t_boundary) :: rboundary
    
    ! Mesh hierarchy.
    type(t_meshHierarchy) :: rmeshHierarchy
    
    ! FE space hierarchy.
    type(t_feHierarchy) :: rfeHierarchy
    
    ! Time coarse mesh.
    type(t_timeDiscretisation) :: rtimecoarse
    
    ! Time hierarchy
    type(t_timescaleHierarchy) :: rtimeHierarchy
    
    ! Space-time hierarchy, the solver works on.
    type(t_spacetimeHierarchy) :: rspaceTimeHierarchy
    
    ! Physics of the equation
    type(t_physics) :: rphysics
    
    ! Boundary conditions on all levels
    type(t_spacetimeBC), dimension(:), pointer :: p_RspaceTimeBC
    
    ! Template matrices on all levels.
    type(t_matvecTemplates), dimension(:), pointer :: p_RmatvecTempl
    
  end type
  
  ! A solver encapsuling the total linear space-time solver.
  type t_linearSpaceTimeSolver
    
    ! Type of the solver.
    ! =0: Defect correction, Block-Jacobi preconditioning
    integer :: csolverType
    
    ! The main solver node.
    type(t_spacetimelinsol) :: rsolver

    ! Preconditioner
    type(t_spacetimelinsol) :: rpreconditioner
    
    ! Coarse grid solver for a MG solver
    type(t_spacetimelinsol) :: rcoarseGridSolver

    ! coarse grid preconditioner
    type(t_spacetimelinsol) :: rcoarsePreconditioner
    
    ! Smoothers for a MG solver
    type(t_spacetimelinsol), dimension(:), pointer :: p_Rsmoothers

    ! Preconditioner of the smoothers 
    type(t_spacetimelinsol), dimension(:), pointer :: p_RsmootherPrecond
    
    ! DEBUG!!! Preconditioners on each level
    type(t_spacetimelinsol), dimension(:), pointer :: p_Rpreconditioners
    
    ! Projection hierarchy for the MG solver
    type(t_sptiProjHierarchy) :: rprojection

    ! Projection hierarchy in space for the MG solver
    type(t_interlevelProjectionHier) :: rprojHierarchySpace
    
  end type
  
contains

  ! ***************************************************************************

  subroutine main_getDiscr(rtriangulation,rdiscr,rboundary,rcollection)
  
  ! Returns a discretisation structure
  
  type(t_triangulation), intent(in) :: rtriangulation
  type(t_blockDiscretisation), intent(out) :: rdiscr
  type(t_collection), intent(inout), optional :: rcollection
  type(t_boundary), intent(in), optional :: rboundary
  
  ! local varianbles
  integer :: cequation, cfespace, creferenceProblem
  
    cequation = rcollection%IquickAccess (1)
    creferenceProblem = rcollection%IquickAccess (2)
    cfespace = rcollection%IquickAccess (3)
    
    ! cfespace: -1=automatic.
    ! =1: 1st order (Q1, Q1~/Q0)
    ! =2: 2nd order (Q2, Q2/QP1)
    
    select case (cequation)
    case (0)
      ! Heat equation
      call spdiscr_initBlockDiscr (rdiscr,2,rtriangulation, rboundary)
      
      select case (cfespace)
      case (-1,1)
        call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q1_2D, CUB_G4X4,&
            rtriangulation, rboundary)
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_Q1_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(2))

      case (2)
        call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q2_2D, CUB_G4X4,&
            rtriangulation, rboundary)
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_Q2_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(2))

      case default
      
        call output_line ("FE-space not supported.")
        call sys_halt()

      end select
    
    case (1)
      ! Stokes equation
      call spdiscr_initBlockDiscr (rdiscr,6,rtriangulation, rboundary)
      
      select case (cfespace)
      case (1)
        call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_EM30_2D, CUB_G4X4,&
            rtriangulation, rboundary)
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_EM30_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(2))
            
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_EM30_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(4))
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_EM30_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(5))

        call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(3),EL_Q0, CUB_G4X4,&
            rtriangulation, rboundary)
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(3), EL_Q0, CUB_G4X4, &
            rdiscr%RspatialDiscr(6))

      case (-1,2)
        call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q2_2D, CUB_G4X4,&
            rtriangulation, rboundary)
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_Q2_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(2))
            
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_Q2_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(4))
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_Q2_2D, CUB_G4X4, &
            rdiscr%RspatialDiscr(5))

        call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(3),EL_QP1, CUB_G4X4,&
            rtriangulation, rboundary)
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(3), EL_QP1, CUB_G4X4, &
            rdiscr%RspatialDiscr(6))

      case default
      
        call output_line ("FE-space not supported.")
        call sys_halt()

      end select

    case (2)
      ! Heat equation in 1D
      call spdiscr_initBlockDiscr (rdiscr,2,rtriangulation, rboundary)
      
      select case (cfespace)
      case (-1,1)
        call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_P1_1D, CUB_G3_1D,&
            rtriangulation, rboundary)
        call spdiscr_deriveSimpleDiscrSc (rdiscr%RspatialDiscr(1), EL_P1_1D, CUB_G3_1D, &
            rdiscr%RspatialDiscr(2))

      case default
      
        call output_line ("FE-space not supported.")
        call sys_halt()

      end select

    case default
    
      call output_line ("Equation not supported.")
      call sys_halt()

    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine main_initLinearSolver (rparlist,rparams,nlevels,rsolver)
  
!<description>
  ! Calculates the optimisation problem based on the parameter
  ! in the parameter list.
!</description>

!<input>
  ! Parameter list containing the DAT file parameters
  type(t_parlist), intent(in) :: rparlist
  
  ! Parameter block of the problem
  type(t_maindata), intent(in) :: rparams
  
  ! Number of levels
  integer, intent(in) :: nlevels
!</input>

!<output>
  ! The solver.
  type(t_linearSpaceTimeSolver), intent(out) :: rsolver
!</output>

!</subroutine>

    ! local variables
    integer :: csolverType
    integer :: iadcgcorr,niteReinit
    real(DP) :: dadcgcorrMin,dadcgcorrMax
    integer :: ilev,nsmoothingSteps,itypeProjection,nmaxiterations
    integer :: ismoother,icoarsegridsolver,ifullcouplingFBGS2,ioutputLevel
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    real(DP) ::  ddampingCoarseGridCorrection,ddamping,drelax,depsrel,depsabs
    integer :: cspaceSolver, cspacePreconditioner
    
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "csolverType", csolverType)
    
    if (csolverType .eq. -1) then
      if (nlevels .eq. 1) then
        csolverType = 0
      else
        csolverType = 1
      end if
    end if
    
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "nsmoothingSteps", nsmoothingSteps)
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "nmaxiterations", nmaxiterations)
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "ismoother", ismoother)
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "icoarsegridsolver", icoarsegridsolver)
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "itypeProjection", itypeProjection)
    call parlst_getvalue_double (rparlist, "SPACETIME-LINEARSOLVER", &
        "ddampingCoarseGridCorrection", ddampingCoarseGridCorrection)
    call parlst_getvalue_double (rparlist, "SPACETIME-LINEARSOLVER", &
        "depsrel", depsrel)
    call parlst_getvalue_double (rparlist, "SPACETIME-LINEARSOLVER", &
        "ddamping", ddamping)
    call parlst_getvalue_double (rparlist, "SPACETIME-LINEARSOLVER", &
        "drelax", drelax)
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "ifullcouplingFBGS2", ifullcouplingFBGS2)
    call parlst_getvalue_int (rparlist, "SPACETIME-LINEARSOLVER", &
        "iadcgcorr", iadcgcorr)
    call parlst_getvalue_double (rparlist, "SPACETIME-LINEARSOLVER", &
        "dadcgcorrMin", dadcgcorrMin,0.5_DP)
    call parlst_getvalue_double (rparlist, "SPACETIME-LINEARSOLVER", &
        "dadcgcorrMax", dadcgcorrMax,2.0_DP)
    
    call parlst_getvalue_int (rparlist, "SPACE-LINEARSOLVER", &
        "csolverType", cspaceSolver)
    call parlst_getvalue_int (rparlist, "SPACE-LINEARSOLVER", &
        "cpreconditioner", cspacePreconditioner)
    
    rsolver%csolverType = csolverType
    
    ! Rematch some variables to solver constanrs
    select case (cspaceSolver)
    case (0)
      cspaceSolver = LINSOL_ALG_UMFPACK4
    case (1)
      cspaceSolver = LINSOL_ALG_MULTIGRID2
    end select

    select case (cspacePreconditioner)
    case (0)
      select case (rparams%rphysics%cequation)
      case (0,2)
        ! Jacobi
        cspacePreconditioner = 0
      case (1)
        ! VANKA
        cspacePreconditioner = 1
      end select
    end select
    
    select case (rsolver%csolverType)
    case (0)
      ! Defect correction with Block Jacobi preconditioning

      ! Create the solver.
      call stls_initBlockJacobi (rsolver%rpreconditioner,rparams%rspacetimeHierarchy,&
          nlevels,1.0_DP,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)
      call stls_initDefCorr (rsolver%rsolver,rparams%rspacetimeHierarchy,nlevels,&
          rsolver%rpreconditioner)
      
      rsolver%rpreconditioner%domega = ddamping
      rsolver%rsolver%nmaxIterations = nmaxiterations
          
    case (2)
      ! Defect correction with Block FBGS preconditioning

      ! Create the solver.
      call stls_initBlockFBGS2 (rsolver%rpreconditioner,rparams%rspacetimeHierarchy,&
          nlevels,drelax,cspaceSolver,cspacePreconditioner,0,rparams%p_RmatvecTempl)
      call stls_initDefCorr (rsolver%rsolver,rparams%rspacetimeHierarchy,nlevels,&
          rsolver%rpreconditioner)
          
      rsolver%rpreconditioner%domega = ddamping
      rsolver%rsolver%nmaxIterations = nmaxiterations
          
    case (1)
      ! MG-Solver with block Jacobi preconditioning.
      
      ! Create the coarse grid solver.

      call parlst_getvalue_double (rparlist, "SPACETIME-COARSEGRIDSOLVER", &
          "ddamping", ddamping)
      call parlst_getvalue_double (rparlist, "SPACETIME-COARSEGRIDSOLVER", &
          "depsrel", depsrel)
      call parlst_getvalue_double (rparlist, "SPACETIME-COARSEGRIDSOLVER", &
          "depsabs", depsabs)
      call parlst_getvalue_double (rparlist, "SPACETIME-COARSEGRIDSOLVER", &
          "drelax", drelax)
      call parlst_getvalue_int (rparlist, "SPACETIME-COARSEGRIDSOLVER", &
          "ifullcouplingFBGS2", ifullcouplingFBGS2)
      call parlst_getvalue_int (rparlist, "SPACETIME-COARSEGRIDSOLVER", &
          "ioutputlevel", ioutputlevel)
      
      if (icoarsegridsolver .eq. 0) then

        call stls_initBlockJacobi (rsolver%rcoarsePreconditioner,rparams%rspacetimeHierarchy,&
            1,1.0_DP,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)

        rsolver%rcoarsePreconditioner%domega = ddamping
        
        call stls_initDefCorr (rsolver%rcoarseGridSolver,rparams%rspacetimeHierarchy,1,&
            rsolver%rcoarsePreconditioner)

      else if (icoarsegridsolver .eq. 1) then

        call stls_initBlockFBGS (rsolver%rcoarsePreconditioner,rparams%rspacetimeHierarchy,&
            1,drelax,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)

        rsolver%rcoarsePreconditioner%domega = ddamping
        
        call stls_initDefCorr (rsolver%rcoarseGridSolver,rparams%rspacetimeHierarchy,1,&
            rsolver%rcoarsePreconditioner)

      else if (icoarsegridsolver .eq. 2) then

        call stls_initBlockFBGS2 (rsolver%rcoarsePreconditioner,rparams%rspacetimeHierarchy,&
            1,drelax,cspaceSolver,cspacePreconditioner,ifullcouplingFBGS2,rparams%p_RmatvecTempl)

        rsolver%rcoarsePreconditioner%domega = ddamping
        
        call stls_initDefCorr (rsolver%rcoarseGridSolver,rparams%rspacetimeHierarchy,1,&
            rsolver%rcoarsePreconditioner)

      else if (icoarsegridsolver .eq. 3) then

        call stls_initBlockJacobi (rsolver%rcoarsePreconditioner,rparams%rspacetimeHierarchy,&
            1,1.0_DP,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)

        rsolver%rcoarsePreconditioner%domega = ddamping
        
        call stls_initBiCGStab (rsolver%rcoarseGridSolver,rparams%rspacetimeHierarchy,1,&
            rsolver%rcoarsePreconditioner)

      else if (icoarsegridsolver .eq. 4) then

        call stls_initBlockFBGS (rsolver%rcoarsePreconditioner,rparams%rspacetimeHierarchy,&
            1,drelax,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)

        rsolver%rcoarsePreconditioner%domega = ddamping
        
        call stls_initBiCGStab (rsolver%rcoarseGridSolver,rparams%rspacetimeHierarchy,1,&
            rsolver%rcoarsePreconditioner)

      else if (icoarsegridsolver .eq. 5) then

        call stls_initBlockFBGS2 (rsolver%rcoarsePreconditioner,rparams%rspacetimeHierarchy,&
            1,drelax,cspaceSolver,cspacePreconditioner,ifullcouplingFBGS2,rparams%p_RmatvecTempl)

        rsolver%rcoarsePreconditioner%domega = ddamping
        
        call stls_initBiCGStab (rsolver%rcoarseGridSolver,rparams%rspacetimeHierarchy,1,&
            rsolver%rcoarsePreconditioner)

      end if
      
      rsolver%rcoarseGridSolver%ioutputLevel = ioutputlevel
      ! rsolver%rcoarseGridSolver%domega = 0.7_DP
      rsolver%rcoarseGridSolver%domega = ddampingCoarseGridCorrection
      rsolver%rcoarseGridSolver%depsrel = depsrel
      rsolver%rcoarseGridSolver%depsabs = depsabs
      
      if (rsolver%rcoarseGridSolver%domega .eq. 0.0_DP) then
        rsolver%rcoarseGridSolver%nmaxIterations = 0
      end if
          
      ! Create the smoothers.
      
      call parlst_getvalue_double (rparlist, "SPACETIME-SMOOTHER", &
          "ddamping", ddamping)
      call parlst_getvalue_double (rparlist, "SPACETIME-SMOOTHER", &
          "drelax", drelax)
      call parlst_getvalue_int (rparlist, "SPACETIME-SMOOTHER", &
          "ifullcouplingFBGS2", ifullcouplingFBGS2)
      call parlst_getvalue_int (rparlist, "SPACETIME-SMOOTHER", &
          "ioutputlevel", ioutputlevel)
      call parlst_getvalue_int (rparlist, "SPACETIME-SMOOTHER", &
          "niteReinit", niteReinit)
      
      allocate (rsolver%p_RsmootherPrecond(nlevels))
      allocate (rsolver%p_Rsmoothers(nlevels))
      allocate (rsolver%p_Rpreconditioners(nlevels))
      do ilev = 2,nlevels
        if (ismoother .eq. 0) then
        
          call stls_initBlockJacobi (rsolver%p_RsmootherPrecond(ilev),rparams%rspacetimeHierarchy,&
              ilev,1.0_DP,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)
        
          rsolver%p_RsmootherPrecond(ilev)%domega = ddamping
          
          call stls_initDefCorr (rsolver%p_Rsmoothers(ilev),rparams%rspacetimeHierarchy,ilev,&
              rsolver%p_RsmootherPrecond(ilev))
          call stls_convertToSmoother (rsolver%p_Rsmoothers(ilev),nsmoothingSteps)

        else if (ismoother .eq. 1) then
        
          call stls_initBlockFBGS (rsolver%p_RsmootherPrecond(ilev),rparams%rspacetimeHierarchy,&
              ilev,drelax,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)
          !rsolver%p_RsmootherPrecond(ilev)%domega = 0.5_DP
          !rsolver%p_RsmootherPrecond(ilev)%domega = 1.0_DP
        
          rsolver%p_RsmootherPrecond(ilev)%domega = ddamping
          
          call stls_initDefCorr (rsolver%p_Rsmoothers(ilev),rparams%rspacetimeHierarchy,ilev,&
              rsolver%p_RsmootherPrecond(ilev))
          call stls_convertToSmoother (rsolver%p_Rsmoothers(ilev),nsmoothingSteps)

        else if (ismoother .eq. 2) then
        
          call stls_initBlockFBGS2 (rsolver%p_RsmootherPrecond(ilev),rparams%rspacetimeHierarchy,&
              ilev,drelax,cspaceSolver,cspacePreconditioner,ifullcouplingFBGS2,rparams%p_RmatvecTempl)
          !rsolver%p_RsmootherPrecond(ilev)%domega = 0.5_DP
        
          rsolver%p_RsmootherPrecond(ilev)%domega = ddamping
          
          call stls_initDefCorr (rsolver%p_Rsmoothers(ilev),rparams%rspacetimeHierarchy,ilev,&
              rsolver%p_RsmootherPrecond(ilev))
          call stls_convertToSmoother (rsolver%p_Rsmoothers(ilev),nsmoothingSteps)

        else if (ismoother .eq. 3) then
        
          call stls_initBlockJacobi (rsolver%p_RsmootherPrecond(ilev),rparams%rspacetimeHierarchy,&
              ilev,1.0_DP,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)
        
          rsolver%p_RsmootherPrecond(ilev)%domega = drelax
          
          call stls_initBiCGStab (rsolver%p_Rsmoothers(ilev),rparams%rspacetimeHierarchy,ilev,&
              rsolver%p_RsmootherPrecond(ilev))
          call stls_convertToSmoother (rsolver%p_Rsmoothers(ilev),nsmoothingSteps)
          rsolver%p_Rsmoothers(ilev)%niteReinit = niteReinit
          rsolver%p_Rsmoothers(ilev)%domega = ddamping

        else if (ismoother .eq. 4) then
        
          call stls_initBlockFBGS (rsolver%p_RsmootherPrecond(ilev),rparams%rspacetimeHierarchy,&
              ilev,drelax,cspaceSolver,cspacePreconditioner,rparams%p_RmatvecTempl)
        
          rsolver%p_RsmootherPrecond(ilev)%domega = ddamping
          
          call stls_initBiCGStab (rsolver%p_Rsmoothers(ilev),rparams%rspacetimeHierarchy,ilev,&
              rsolver%p_RsmootherPrecond(ilev))
          call stls_convertToSmoother (rsolver%p_Rsmoothers(ilev),nsmoothingSteps)

        else if (ismoother .eq. 5) then
        
          call stls_initBlockFBGS2 (rsolver%p_RsmootherPrecond(ilev),rparams%rspacetimeHierarchy,&
              ilev,drelax,cspaceSolver,cspacePreconditioner,ifullcouplingFBGS2,rparams%p_RmatvecTempl)
        
          rsolver%p_RsmootherPrecond(ilev)%domega = ddamping
          
          call stls_initBiCGStab (rsolver%p_Rsmoothers(ilev),rparams%rspacetimeHierarchy,ilev,&
              rsolver%p_RsmootherPrecond(ilev))
          call stls_convertToSmoother (rsolver%p_Rsmoothers(ilev),nsmoothingSteps)

        end if

        rsolver%p_Rsmoothers(ilev)%ioutputlevel = ioutputlevel
        
        !call stls_initDefCorr (rsolver%p_Rpreconditioners(ilev),rparams%rspacetimeHierarchy,ilev,&
        !    rsolver%p_RsmootherPrecond(ilev))
        !rsolver%p_Rpreconditioners(ilev)%nmaxIterations = 50
      end do
      
      ! Initialise the interlevel projection in space
      call mlprj_initPrjHierarchy(rsolver%rprojHierarchySpace,1,rparams%rfeHierarchy%nlevels)
      do ilev = 1,rparams%rfeHierarchy%nlevels

        p_rfeSpaceLevel => rparams%rfeHierarchy%p_rfeSpaces(ilev)
        call mlprj_initPrjHierarchyLevel(rsolver%rprojHierarchySpace,ilev,&
            p_rfeSpaceLevel%p_rdiscretisation)
            
      end do

      call mlprj_commitPrjHierarchy (rsolver%rprojHierarchySpace)
      
      ! Initialise the interlevel projection in space/time
      call sptipr_initProjection (rsolver%rprojection,rparams%rspacetimeHierarchy,&
          rsolver%rprojHierarchySpace,rparams%rphysics,itypeProjection)
      
      call parlst_getvalue_double (rparlist, "SPACETIME-LINEARSOLVER", &
          "depsrel", depsrel)

      ! Create the solver.
      call stls_initMultigrid (rsolver%rsolver,rparams%rspacetimeHierarchy,nlevels,&
          rsolver%rprojection, rsolver%rcoarseGridSolver, RpostSmoothers=rsolver%p_Rsmoothers)
          !,&
          !Rpreconditioners=rsolver%p_Rpreconditioners)
      rsolver%rsolver%nmaxIterations = nmaxiterations
      rsolver%rsolver%depsrel = depsrel
      rsolver%rsolver%iadcgcorr = iadcgcorr
      rsolver%rsolver%dadcgcorrMin = dadcgcorrMin
      rsolver%rsolver%dadcgcorrMax = dadcgcorrMax
          
    end select
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine main_doneLinearSolver (rsolver)
  
!<description>
  ! Releases the linear solver.
!</description>

!<inputoutput>
  ! The solver to release.
  type(t_linearSpaceTimeSolver),intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev

    select case (rsolver%csolverType)
    case (0)
      ! Defect correction with Block Jacobi preconditioning

      ! Release the solver.
      call stls_done(rsolver%rpreconditioner)
      call stls_done(rsolver%rsolver)
      
    case (2)
      ! Defect correction with Block FBGS preconditioning

      ! Release the solver.
      call stls_done(rsolver%rpreconditioner)
      call stls_done(rsolver%rsolver)

    case (1)
      ! MG-Solver with block Jacobi preconditioning.
      call mlprj_releasePrjHierarchy (rsolver%rprojHierarchySpace)
      call sptipr_doneProjection(rsolver%rprojection)

      call stls_done(rsolver%rcoarsePreconditioner)
      call stls_done(rsolver%rcoarseGridSolver)
          
      ! Release the smoothers.
      do ilev = 2,rsolver%rsolver%ilevel
        call stls_done(rsolver%p_RsmootherPrecond(ilev))
        call stls_done(rsolver%p_Rsmoothers(ilev))
      end do
      deallocate (rsolver%p_RsmootherPrecond)
      deallocate (rsolver%p_Rsmoothers)
      
    case default
    
      call output_line ("main_doneLinearSolver: Unknown solver!")
      call sys_halt()
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine main_calculate (rparlist)
  
!<description>
  ! Calculates the optimisation problem based on the parameter
  ! in the parameter list.
!</description>

!<input>
  ! Parameter list containing the DAT file parameters
  type(t_parlist), intent(in) :: rparlist
!</input>

!</subroutine>

    ! Local variables
    type(t_maindata) :: rparams
    type(t_collection) :: rcollection
    integer :: nspacelevels,ntimelevels,nminlevelspace,nminleveltime,ilev,ntstepscoarse
    integer :: nmaxminleveltime,csolverStatus,nminminleveltime
    real(DP) :: dtheta
    type(t_spaceTimeVector) :: rrhs, rsolution, rtemp
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    type(t_spaceTimeMatrix), dimension(:), pointer :: p_Rmatrices
    type(t_triangulation) :: rtria1D
    integer :: iwriteUCD, cfespace, icalcError
    character(LEN=SYS_STRLEN) :: smesh, sboundary
    integer :: ithetaschemetype
    
    type(t_linearSpaceTimeSolver) :: rlinearSolver
    
    call parlst_getvalue_int (rparlist, "SPACETIME-DISCRETISATION", &
        "nspacelevels", nspacelevels)
    call parlst_getvalue_int (rparlist, "SPACETIME-DISCRETISATION", &
        "ntimelevels", ntimelevels)
    call parlst_getvalue_int (rparlist, "SPACETIME-DISCRETISATION", &
        "nminlevelspace", nminlevelspace)
    call parlst_getvalue_int (rparlist, "SPACETIME-DISCRETISATION", &
        "nmaxminleveltime", nmaxminleveltime)
    call parlst_getvalue_int (rparlist, "SPACETIME-DISCRETISATION", &
        "nminminleveltime", nminminleveltime)
    call parlst_getvalue_double (rparlist, "SPACETIME-DISCRETISATION", &
        "dtheta", dtheta)
    call parlst_getvalue_int (rparlist, "SPACETIME-DISCRETISATION", &
        "ithetaschemetype", ithetaschemetype)
        
    call parlst_getvalue_string (rparlist, "PARAMTRIANG", &
        "sboundary", sboundary,bdequote=.true.)    
    call parlst_getvalue_string (rparlist, "PARAMTRIANG", &
        "smesh", smesh,bdequote=.true.)    

    call parlst_getvalue_int (rparlist, "POSTPROC", &
        "iwriteUCD", iwriteUCD)

    call parlst_getvalue_int (rparlist, "POSTPROC", &
        "icalcError", icalcError)

    call parlst_getvalue_int (rparlist, "SPACETIME-DISCRETISATION", &
        "cfespace", cfespace, -1)

    do nminleveltime = nminminleveltime,nmaxminleveltime

      call output_separator (OU_SEP_MINUS)
      call output_line ("Level "//trim(sys_siL(nminleveltime+ntimelevels-1,10)))
      call output_lbrk

      ntstepscoarse = 5*2**(nminleveltime-1) !1 !5*2**(nminleveltime-1)
    
      ! Get parameters
      call cb_getPhysicsHeatEqn(rparlist,rparams%rphysics)
      
      ! Read boundary and generate mesh hierarchy
      select case (rparams%rphysics%cequation)
      case (0,1)
        
        ! 2D mesh
        call boundary_read_prm(rparams%rboundary, sboundary)
        call mshh_initHierarchy (rparams%rmeshHierarchy,nspacelevels,smesh,&
            NDIM2D,nminlevelspace-1,rparams%rboundary)

      case (2)

        ! 1D mesh
        call tria_createRawTria1D(rtria1D, 0.0_DP, 1.0_DP, 2**(nminlevelspace-1))
        call tria_initStandardMeshFromRaw (rtria1D)
        call mshh_initHierarchy (rparams%rmeshHierarchy,rtria1D,0,nspacelevels,&
            cdupFlag=TR_SHARE_NONE)
        call tria_done (rtria1D)

      case default
          
        call output_line ("Equation not supported.")
        call sys_halt()
    
      end select

      call mshh_refineHierarchy2lv(rparams%rmeshHierarchy,nspacelevels,&
          rboundary=rparams%rboundary)
      
      ! Space hierarchy
      rcollection%IquickAccess (1) = rparams%rphysics%cequation
      rcollection%IquickAccess (2) = rparams%rphysics%creferenceProblem
      rcollection%IquickAccess (3) = cfespace
      call fesph_createHierarchy (rparams%rfeHierarchy,nspacelevels,&
          rparams%rmeshHierarchy,main_getDiscr,rcollection,rparams%rboundary)

      ! Time coarse mesh and time hierarchy
      call tdiscr_initOneStepTheta (rparams%rtimecoarse, 0.0_DP, 1.0_DP, ntstepscoarse, dtheta)
      
      ! The ithetaschemetype flag is saved as tag.
      rparams%rtimecoarse%itag = ithetaschemetype

      call tmsh_createHierarchy (rparams%rtimecoarse,rparams%rtimeHierarchy,&
          0,ntimelevels)
      
      ! Space-time hierarchy.
      call sth_initHierarchy (rparams%rspacetimeHierarchy,&
          rparams%rfeHierarchy,rparams%rtimeHierarchy)
      call sth_defineHierarchyByCoarsening (rparams%rspacetimeHierarchy,&
          nspacelevels,nspacelevels,1,ntimelevels)
      
      ! Boundary conditions on all levels.
      allocate (rparams%p_RspaceTimeBC(rparams%rspacetimeHierarchy%nlevels))
      allocate (rparams%p_RmatvecTempl(nspacelevels))
      allocate (p_Rmatrices(rparams%rspacetimeHierarchy%nlevels))
      
      do ilev=1,rparams%rspacetimeHierarchy%nlevels
      
        ! Boundary conditions
        call spop_createBC(ilev,rparams%rspacetimeHierarchy,rparams%rphysics,&
            rparams%p_RspaceTimeBC(ilev))

        ! Matrix
        call stmv_createMatrix (0,ilev,rparams%rspacetimeHierarchy,rparams%rphysics,&
            rparams%p_RspaceTimeBC(ilev),rparams%p_RmatvecTempl,p_Rmatrices(ilev))

      end do
      
      do ilev=1,nspacelevels
        p_rspaceDiscr => rparams%rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation
        ! Template matrices
        call spop_calcMatrices (p_rspaceDiscr,rparams%rphysics,rparams%p_RmatvecTempl(ilev))
      end do

      ! Get the discretisation info of the max. level.
      call sth_getLevel (rparams%rspacetimeHierarchy,rparams%rspacetimeHierarchy%nlevels,&
        p_rfeSpaceLevel,p_rtimeDiscr)
      p_rspaceDiscr => p_rfeSpaceLevel%p_rdiscretisation
      
      ! Create the right hand side.
      call sptivec_initVector (rrhs,p_rtimeDiscr,p_rspaceDiscr)
      call strhs_assembleRHS (rparams%rphysics,rrhs)
      
      ! Create the solver.
      call main_initLinearSolver (rparlist,rparams,ntimelevels,rlinearSolver)
      
      ! Attach matrices
      do ilev = 1,rparams%rspacetimeHierarchy%nlevels
        call stls_setMatrix (rlinearSolver%rsolver,ilev,p_Rmatrices(ilev))
      end do

      ! Apply boundary conditions, prepare to solve.    
      call sptivec_initVector (rtemp,p_rtimeDiscr,p_rspaceDiscr)
      call sptivec_initVector (rsolution,p_rtimeDiscr,p_rspaceDiscr)
      call sptivec_clearVector (rsolution)
      
      call spop_applyBC (rparams%p_RspaceTimeBC(rparams%rspacetimeHierarchy%nlevels), &
          SPOP_RHS, rrhs)
      call spop_applyBC (rparams%p_RspaceTimeBC(rparams%rspacetimeHierarchy%nlevels), &
          SPOP_SOLUTION, rsolution)
      
      !call sptivec_saveToFileSequence (rrhs,"('ns/rhs_"//&
      !    trim(sys_siL(nminleveltime,2))//"-"//&
      !    trim(sys_siL(ntimelevels,2))//"lv.txt.',I5.5)",.true.)

      ! Initial defect    
      !call stpp_printDefectSubnorms (p_Rmatrices(rparams%rspacetimeHierarchy%nlevels),&
      !    rsolution,rrhs,rtemp)
      
      ! Solve
      call stls_initData (rlinearSolver%rsolver)
      call stls_solveAdaptively (rlinearSolver%rsolver,rsolution,rrhs,rtemp)
      csolverStatus = rlinearSolver%rsolver%csolverStatus
      call stls_doneData (rlinearSolver%rsolver)
      
      call main_doneLinearSolver (rlinearSolver)

      call spop_applyBC (rparams%p_RspaceTimeBC(rparams%rspacetimeHierarchy%nlevels), &
          SPOP_SOLUTION, rsolution)
      
      ! Postprocessing
      !call stpp_printDefectSubnorms (p_Rmatrices(rparams%rspacetimeHierarchy%nlevels),&
      !    rsolution,rrhs,rtemp)
      !call sptivec_saveToFileSequence (rsolution,"('ns/solution_"//&
      !    trim(sys_siL(nminleveltime,2))//"-"//&
      !    trim(sys_siL(ntimelevels,2))//"lv.txt.',I5.5)",.true.)
      call stpp_postproc (rparams%rphysics,rsolution,iwriteUCD .ne. 0,icalcError .ne. 0)
      
      do ilev=1,rparams%rspacetimeHierarchy%nlevels
        call stmv_releaseMatrix(p_Rmatrices(ilev))
        call spop_releaseBC(rparams%p_RspaceTimeBC(ilev))
      end do

      do ilev=1,nspacelevels
        call spop_releaseMatrices(rparams%p_RmatvecTempl(ilev))
      end do

      deallocate(p_Rmatrices)
      deallocate(rparams%p_RmatvecTempl)
      deallocate(rparams%p_RspaceTimeBC)
      
      call sptivec_releaseVector(rtemp)
      call sptivec_releaseVector(rsolution)
      call sptivec_releaseVector(rrhs)
      
      call sth_doneHierarchy(rparams%rspacetimeHierarchy)
      
      call tmsh_releaseHierarchy(rparams%rtimeHierarchy)
      call tdiscr_done(rparams%rtimecoarse)
      
      call fesph_releaseHierarchy (rparams%rfeHierarchy)
          
      call mshh_releaseHierarchy(rparams%rmeshHierarchy)
      call boundary_release(rparams%rboundary)

      if (csolverStatus .eq. 2) then
        call output_line ("Stopping iteration due to divergence.")
        exit
      end if

    end do    

  end subroutine

  ! ***************************************************************************

  subroutine main_optc
    
    ! Program parameters
    type(t_parlist) :: rparlist
    character(len=SYS_STRLEN) :: soutput
    character(LEN=SYS_STRLEN) :: smaster
    logical :: bexists
    
    ! The very first thing in every application: 
    ! Initialise system-wide settings:
    call system_init()
    
    ! Read the program parameters.
    call parlst_init (rparlist)
    
    ! Check if a command line parameter specifies the master.dat file.
    call sys_getcommandLineArg(1,smaster,sdefault='./data/nonstoptcontrol.dat')

    ! Read the file 'master.dat'.
    ! If that does not exist, try to manually read files with parameters from a
    ! couple of files.
    inquire(file=smaster, exist=bexists)
    
    if (.not. bexists) then
      call output_line ("Cannot find INI file: "//trim(smaster))
      call sys_halt()
    end if
    
    call parlst_readfromfile(rparlist,smaster)
    
    ! Initialise log file for output.
    call parlst_getvalue_string (rparlist, "OUTPUT", &
        "soutputFile", soutput,bdequote=.true.)    
    call output_init (soutput)
    OU_LINE_LENGTH = 132
    cdefaultDateTimeLogPolicy = OU_DTP_ADDDATETIME
    cdatetimeLogFormat = 1
    
    call parlst_info(rparlist)
    
    ! Now we can really start!
    !
    ! Initialise the storage management: 
    call storage_init(999, 100)
    call exstor_init (999,100)
    
    ! Initialise the parser
    call fparser_init ()
    
    ! Call the problem to solve. 
    call output_lbrk ()
    call output_line ('Calculating problem')
    
    call main_calculate (rparlist)

    ! Release the parser
    call fparser_done ()

    ! Information about external storage usage
    call output_lbrk ()
    call exstor_info (bprintHandles=.true.)
    call exstor_done ()

    ! Print out heap statistics - just to check if everything
    ! is cleaned up.
    ! This should display 'Handles in use=0' and 'Memory in use=0'!
    call output_lbrk ()
    call storage_info(.true.)
    
    ! Clean up the storage management, parameter list, finish
    call storage_done()
    call parlst_done (rparlist)
    
  end subroutine

end module
