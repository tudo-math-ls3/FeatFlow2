!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2nonlinearcore </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the routines to solve the stationary core equation
!# of the problem with a nonlinear solver.
!# (This is comparable to the NSDEF routine in older FEAT versions)
!#
!# The discretised core equation reads at the moment:
!#
!#   alpha*M*u + theta*Stokes*u + gamma*N(u)u + B*p = f1
!#                                            B^T*u = f2
!#
!# with:
!#   alpha = 0/1    - for time dependence; 
!#                    =0 for stationary problem
!#   theta in [0,1] - for time dependent THETA scheme; 
!#                    =1 for stationary problem
!#   gamma = 0/1    - for Stokes/Navier-Stokes
!#
!# "Stokes" is the Stokes matrix (=nu*Laplace) in deformation or gradient
!# formulation.
!# N(u)u is the nonlinearity with stabilisation. This may also include things
!# like Newton matrices, etc.
!#
!# The core equation is abstractly written a nonlinear system of the form
!#
!#  $$ A(x)x = b $$
!#
!# and solved with the nonlinear solver from the kernel, using the defect
!# correction approach
!#
!#  $$  x_{n+1}  =  x_n  +  \omega_n C^{-1} ( b - A(x_n) x_n )  $$
!#
!# where $C^{-1}$ means to apply a suitable preconditioner (inverse mass
!# matrix, apply the linearised $A(x_n)^-1$ with multigrid, apply Newton or 
!# do something similar). 
!#
!# The following routines can be found here:
!#
!#  1.) c2d2_getNonlinearSolver
!#      -> Initialise nonlinear solver configuration with information 
!#         from INI/DAT files
!#
!#  2.) c2d2_createNonlinearLoop
!#      -> Creates a nonlinear iteration structure, allocates memory.
!#
!#  3.) c2d2_setupCoreEquation
!#      -> Initialises the parameters of the core equation and creates links to
!#         matrices on all levels.
!#
!#  4.) c2d2_solveCoreEquation
!#      -> Starts the nonlinear iteration to solve the core equation.
!#
!#  5.) c2d2_doneNonlinearLoop
!#      -> Releases allocated memory in the nonlinear iteration structure.
!#
!# Callback routines for the nonlinear solver:
!#
!#  6.) c2d2_getDefect
!#      -> Callback routine. Calculate nonlinear defect
!#
!#  7.) c2d2_getOptimalDamping
!#     -> Auxiliary routine. Calculate optimal damping parameter
!#
!#  8.) c2d2_precondDefect
!#      -> Callback routine. Preconditioning of nonlinear defect
!#
!#  9.) c2d2_resNormCheck
!#      -> Callback routine: Check the residuals for convergence
!#
!# Auxiliary routines:
!#
!# 10.) c2d2_saveNonlinearLoop
!#      -> Saves the parameters of a nonlinear iteration structure to a
!#         collection
!#
!# 11.) c2d2_restoreNonlinearLoop
!#      -> Restores a nonlinear iteration structure from a collection
!#
!# 12.) c2d2_removeNonlinearLoop
!#      -> Releases the parameters of a nonlinear iteration structure
!#         from a collection.
!#
!# 13.) c2d2_assembleLinearisedMatrices
!#      -> Assembles the linearised nonlinear matrix/matrices.
!#
!# 14.) c2d2_assembleConvDiffDefect
!#      -> Assemble the convection-diffusion part of the nonlinear defect.
!#
!# 15.) c2d2_assembleNonlinearDefect
!#      -> Assembles the nonlinear defect to a given RHS and a given solution
!#         vector.
!#
!# 16.) c2d2_getDefectNorm
!#      -> Calculate a couple of norms of a residual vector
!#
!# To solve a system with the core equation, one has to deal with two
!# main structures. On one hand, one has a nonlinear iteration structure t_nlsolNode
!# from the kernel; this is initialised by c2d2_getNonlinearSolver.
!# On the other hand, one has to maintain a 'nonlinear iteration structure'
!# of type t_ccNonlinearIteration, which configures the core equation and
!# specifies parameters for the solver how to work.
!#
!# The basic procedure for dealing with the core equation is as follows:
!#
!#  a) c2d2_getNonlinearSolver  -> Basic initialisation on the nonlinear solver
!#
!#  b) c2d2_createNonlinearLoop -> Basic initialisation of the core equation
!#                                 structures
!#
!#  e) Initialise further parameters in the core equation structure manually
!#     (e.g. preconditioner, pointer to matrices, ...).
!#     It's important, that the 'outer' application initialises pointers to
!#     matrices, otherwise nothing will work!
!#     This all has to be done with the nonlinear-iteration-structure directly.
!#
!#  d) c2d2_setupCoreEquation -> Initialises constants in the core equation
!#
!#  e) c2d2_solveCoreEquation -> Solve the core equation with the nonlinear
!#                               solver structure from c2d2_getNonlinearSolver.
!#
!#  f) c2d2_doneNonlinearLoop -> Release core equation structure
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2nonlinearcore

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  USE linearsolverautoinitialise
  USE matrixrestriction
  USE trilinearformevaluation
  
  USE collection
  USE convection
    
  USE cc2dmedium_callback
  
  IMPLICIT NONE
  
!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  INTEGER, PARAMETER :: CCPREC_NONE         = -1

  ! Preconditioning with inverse mass matrix (not yet implemented)
  INTEGER, PARAMETER :: CCPREC_INVERSEMASS   = 0

  ! Preconditioning by linear solver, solving the linearised system
  INTEGER, PARAMETER :: CCPREC_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  INTEGER, PARAMETER :: CCPREC_NEWTON        = 2

  ! Preconditioning by dynamic Newton-Iteration (uses defect correction
  ! and switches automatically to Newton if the error is small enough)
  INTEGER, PARAMETER :: CCPREC_NEWTONDYNAMIC = 3

!</constantblock>

!</constants>

  
!<types>

!<typeblock>

  ! This structure controls the Newton iteration -- i.e. the preconditioning
  ! with the Frechet derivative of the Navier--Stokes equation, which
  ! can lead to quadratic covergence of the nonlinear solver.
  ! As Newton works only in the basin of attraction of the solution,
  ! the parameters in this structure allow to define a switching criterion
  ! when to use Newton. In the first couple of iterations, defect correction
  ! is used, while the iteration switches to Newton if the residuum is small
  ! enough.
  ! This block is used if CCPREC_NEWTONDYNAMIC is used as preconditioner.
  TYPE t_ccDynamicNewtonControl
  
    ! Minimum number of usul fix point iteration before to switch to
    ! preconfitioning with the Newton matrix. (IFIXMIN)

    INTEGER :: nminFixPointIterations = 0

    ! Maximum number of usul fix point iteration before to switch to
    ! preconfitioning with the Newton matrix. (IFIXMAX)

    INTEGER :: nmaxFixPointIterations = 999

    ! Norm of absolute residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Stamndard value = 1E-5.

    REAL(DP) :: depsAbsNewton = 1.0E-5_DP

    ! Norm of relative residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Standard value = 1E99 -> The absolute residuum counts.

    REAL(DP) :: depsRelNewton = 1.0E99_DP
  
  END TYPE
  
!</typeblock>


!<typeblock>

  ! Preconditioner structure for CCxD. This structure saves the configuration of the
  ! preconditioner that is used during the nonlinear iteration. 
  
  TYPE t_ccPreconditioner
  
    ! Type of preconditioner.
    ! This is one of the CCPREC_xxxx flags as defined above (CCPREC_INVERSEMASS for
    ! preconditioning with inverse mass matrix, CCPREC_LINEARSOLVER for solving a linear
    ! system, CCPREC_NEWTON for a Newton iteration,...)
    INTEGER :: ctypePreconditioning = CCPREC_NONE
    
    ! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock), POINTER :: p_rprojection
    
    ! A filter chain that is used for implementing boundary conditions or other
    ! things when invoking the linear solver.
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! Configuration block for the adaptive Newton preconditioner.
    ! Is only valid if ctypePreconditioning=CCPREC_NEWTONDYNAMIC!
    TYPE(t_ccDynamicNewtonControl) :: radaptiveNewton

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    TYPE(t_vectorScalar), POINTER :: p_rtempVectorSc

    ! Temporary scalar vector; used for calculating the optimal damping
    ! parameter.
    TYPE(t_vectorScalar), POINTER :: p_rtempVectorSc2

  END TYPE

!</typeblock>

!<typeblock>

  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it's noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It's more or less
  ! a collection if different flags.
  TYPE t_ccFinalAssemblyInfo
    ! This flag is set to YES if the B matrices must be assembled
    ! transposedly. This may be necessary for special VANCA type smoothers
    ! if a linear solver is used as preconditioner.
    INTEGER :: iBmatricesTransposed = 0
    
    ! Whether to use 'adaptive matrices', i.e. set up coarse grid matrices
    ! with the help of fine grid matrices. This is used for very special
    ! discretisations only (e.g. Q1~/Q0). =0: deactivate
    INTEGER :: iadaptiveMatrices    = 0
    
    ! A configuration parameter for adaptive matrices.
    REAL(DP) :: dadMatThreshold     = 0.0_DP
  END TYPE

!</typeblock>

!<typeblock>

  ! Represents the core equation on one level of the discretisation.
  ! Collects all information that are necessary to assemble the 
  ! (linearised) system matrix and RHS vector.
  TYPE t_cccoreEquationOneLevel
  
    ! The (linearised) system matrix for that specific level. 
    TYPE(t_matrixBlock), POINTER :: p_rmatrix => NULL()

    ! Stokes matrix for that specific level (=nu*Laplace)
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes => NULL()

    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB1 => NULL()

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB2 => NULL()

    ! Mass matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixMass => NULL()

    ! Block matrix, which is used in the defect correction / Newton
    ! algorithm as preconditioner matrix of the correspnding underlying
    ! linear sytem. Is usually the (linearise) system matrix or
    ! a Newton matrix. This matrix is changed during the
    ! nonlinear iteration and used e.g. if a linear solver (Multigrid) is
    ! used for preconditioning.
    TYPE(t_matrixBlock), POINTER :: p_rmatrixPreconditioner
  
    ! Velocity coupling matrix $A_{12}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixVelocityCoupling12 => NULL()

    ! Velocity coupling matrix $A_{12}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixVelocityCoupling21 => NULL()
    
  END TYPE

!</typeblock>

!<typeblock>

  ! This configuration block configures all parameters that are needed
  ! by the callback routines to perform the nonlinear iteration.
  ! On start of the program, a structure of this type is initialised.
  ! The entries of this structure aresaved to the collection structure.
  ! When a callback routine is called, the structure is rebuild from
  ! the collection. When he nonlinear iteration is finished, the
  ! parameters are removed from the colletion again.
  TYPE t_ccNonlinearIteration
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    REAL(DP) :: dalpha = 0.0_DP
    
    ! THETA-parameter that controls the weight of the Stokes matrix
    ! in the core equation. =1.0 for stationary simulations.
    REAL(DP) :: dtheta = 0.0_DP
    
    ! GAMMA-parameter that controls the weight in front of the
    ! nonlinearity. =1.0 for Navier-Stokes, =0.0 for Stokes equation.
    REAL(DP) :: dgamma = 0.0_DP
    
    ! Minimum allowed damping parameter; OMGMIN
    REAL(DP) :: domegaMin = 0.0_DP
    
    ! Maximum allowed damping parameter; OMGMAX
    REAL(DP) :: domegaMax = 2.0_DP
    
    ! Minimum discretisation level
    INTEGER :: NLMIN = 1
    
    ! Maximum discretisation level
    INTEGER :: NLMAX = 1
    
    ! Pointer to the solution vector that is changed during the iteration
    TYPE(t_vectorBlock), POINTER :: p_rsolution => NULL()
    
    ! Pointer to the RHS vector
    TYPE(t_vectorBlock), POINTER :: p_rrhs => NULL()
    
    ! A t_ccPreconditioner saving information about the preconditioner.
    TYPE(t_ccPreconditioner) :: rpreconditioner
    
    ! A t_ccFinalAssemblyInfo structure that saves information about
    ! special 'tweaks' in matrices such that everything works.
    TYPE(t_ccFinalAssemblyInfo) :: rfinalAssembly
    
    ! An array of t_cccoreEquationOneLevel structures for all levels
    ! of the discretisation.
    TYPE(t_cccoreEquationOneLevel), DIMENSION(:), POINTER :: RcoreEquation => NULL()
    
    ! Auxiliary variable: Saves the initial defect in the nonlinear iteration
    REAL(DP), DIMENSION(2) :: DresidualInit = 0.0_DP
    
    ! Auxiliary variable: Saves the last defect in the nonlinear iteration
    REAL(DP), DIMENSION(2) :: DresidualOld = 0.0_DP

    ! Auxiliary variable: Norm of the relative change = norm of the 
    ! preconditioned residual in the nonlinear iteration
    REAL(DP), DIMENSION(3) :: DresidualCorr = 0.0_DP
    
    ! Auxiliary variable: Convergence criteria of the nonlinear solver
    REAL(DP), DIMENSION(5) :: DepsNL = 0.0_DP
    
    ! Auxiliary variable: Last calculated damping parameter
    REAL(DP) :: domegaNL = 0.0_DP
    
    ! Auxiliary variable: Convergence rate of linear solver (if this is
    ! applied as preconditioner).
    REAL(DP) :: drhoLinearSolver = 0.0_DP
  END TYPE

!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************
  ! Routines to create a nonlinear iteration structure, to save it
  ! to a collection, to rebuild it from there and to clean it up.
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_createNonlinearLoop (rnonlinearIteration,NLMIN,NLMAX)
  
!<description>
  ! This routine creates a nonlinear iteration structure. The structure is
  ! initialised to handle NLMAX-NLMIN+1 discretisation levels.
!</description>

!<input>
  ! Minimum discretisation level to be maintained
  INTEGER, INTENT(IN) :: NLMIN
  
  ! Maximum discretisation level to be maintained. The maximum level coincides
  ! with the level where to solve the system.
  INTEGER, INTENT(IN) :: NLMAX
!</input>

!<output>
  ! A nonlinear iteration structure to be initialised.
  TYPE(t_ccNonlinearIteration), INTENT(OUT) :: rnonlinearIteration
!</output>

!</subroutine>

    rnonlinearIteration%NLMIN = NLMIN
    rnonlinearIteration%NLMAX = NLMAX

    ! Initialise the matrix pointers on all levels that we have to maintain.
    ALLOCATE(rnonlinearIteration%RcoreEquation(NLMIN:NLMAX))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneNonlinearLoop (rnonlinearIteration)
  
!<description>
  ! Releases allocated memory in the nonlinear iteration structure.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that should be cleaned up.
  TYPE(t_ccNonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>
    
    IF (ASSOCIATED(rnonlinearIteration%RcoreEquation)) &
      DEALLOCATE(rnonlinearIteration%RcoreEquation)

    rnonlinearIteration%NLMIN = 0
    rnonlinearIteration%NLMAX = 0

  END SUBROUTINE

  ! ***************************************************************************
  ! Routines to save the nonlinear iteration structure to a collection and
  ! rebuild it from there.
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_saveNonlinearLoop (rnonlinearIteration,rcollection)
  
!<description>
  ! Saves the parameters in the rnonlinearIteration structure to a given
  ! collection.
!</description>

!<input>
  ! A nonlinear iteration structure with data.
  TYPE(t_ccNonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</input>

!<inputoutput>
  ! The collection structure where to save the data of rnonlinearIteration to.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ilevel

    ! Save the structure members
    CALL collct_setvalue_real(rcollection,'CCNL_OMEGAMIN',&
        rnonlinearIteration%domegaMin,.TRUE.)
    
    CALL collct_setvalue_real(rcollection,'CCNL_OMEGAMAX',&
        rnonlinearIteration%domegaMax,.TRUE.)
    
    CALL collct_setvalue_vec(rcollection,'CCNL_RHS',&
        rnonlinearIteration%p_rrhs,.TRUE.)
        
    CALL collct_setvalue_vec(rcollection,'CCNL_SOLUTION',&
        rnonlinearIteration%p_rsolution,.TRUE.)

    CALL collct_setvalue_real(rcollection,'CCNL_ALPHA',&
        rnonlinearIteration%dalpha,.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_THETA',&
        rnonlinearIteration%dtheta,.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_GAMMA',&
        rnonlinearIteration%dgamma,.TRUE.)
        
    ! Save the preconditioner
    CALL collct_setvalue_int(rcollection,'CCNL_PRECONDITIONER',&
        rnonlinearIteration%rpreconditioner%ctypePreconditioning,.TRUE.)

    ! Which preconditioner is chosen?
    SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Save temporary vectors
      CALL collct_setvalue_vecsca(rcollection,'CCNL_RTEMPSCALAR',&
          rnonlinearIteration%rpreconditioner%p_rtempVectorSc,.TRUE.)
     
      CALL collct_setvalue_vecsca(rcollection,'CCNL_RTEMP2SCALAR',&
          rnonlinearIteration%rpreconditioner%p_rtempVectorSc2,.TRUE.)

      ! Save the filter chain
      CALL collct_setvalue_fchn(rcollection,'CCNL_FILTERCHAIN',&
          rnonlinearIteration%rpreconditioner%p_RfilterChain,.TRUE.)

      ! Add the interlevel projection structure to the collection; we can
      ! use it later for setting up nonlinear matrices.
      CALL collct_setvalue_ilvp(rcollection,'CCNL_ILVPROJECTION',&
          rnonlinearIteration%rpreconditioner%p_rprojection,.TRUE.)

      ! Put the prepared solver node to the collection for later use.
      ! Remember, it's our preconditioner we need during the nonlinear
      ! iteration!
      CALL collct_setvalue_linsol(rcollection,'CCNL_LINSOLVER',&
          rnonlinearIteration%rpreconditioner%p_rsolverNode,.TRUE.)
      
      ! Remember the Newton parameters
      CALL collct_setvalue_int(rcollection,'CCNL_NEWTONMINFP',&
          rnonlinearIteration%rpreconditioner%radaptiveNewton%nminFixPointIterations,&
          .TRUE.)

      CALL collct_setvalue_int(rcollection,'CCNL_NEWTONMAXFP',&
          rnonlinearIteration%rpreconditioner%radaptiveNewton%nmaxFixPointIterations,&
          .TRUE.)

      CALL collct_setvalue_real(rcollection,'CCNL_NEWTONEPSABS',&
          rnonlinearIteration%rpreconditioner%radaptiveNewton%depsAbsNewton,&
          .TRUE.)

      CALL collct_setvalue_real(rcollection,'CCNL_NEWTONEPSREL',&
          rnonlinearIteration%rpreconditioner%radaptiveNewton%depsRelNewton,&
          .TRUE.)
    END SELECT

    CALL collct_setvalue_int(rcollection,'CCNL_IBMATTRANSPOSED',&
        rnonlinearIteration%rfinalAssembly%iBmatricesTransposed,.TRUE.)

    CALL collct_setvalue_int(rcollection,'CCNL_IADAPTIVEMATRIX',&
        rnonlinearIteration%rfinalAssembly%iadaptiveMatrices,.TRUE.)
                              
    CALL collct_setvalue_real(rcollection,'CCNL_DADMATTHRESHOLD',&
        rnonlinearIteration%rfinalAssembly%dadMatThreshold,.TRUE.)

    ! Min/Max level
    CALL collct_setvalue_int(rcollection,'CCNL_NLMIN',&
        rnonlinearIteration%NLMIN,.TRUE.)
        
    CALL collct_setvalue_int(rcollection,'CCNL_NLMAX',&
        rnonlinearIteration%NLMAX,.TRUE.)
    
    DO ilevel = rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
    
      ! Save the structure members on every level
      CALL collct_setvalue_mat(rcollection,'CCNL_SYSTEMMAT',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrix,.TRUE.,ilevel)

      CALL collct_setvalue_matsca(rcollection,'CCNL_STOKES',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixStokes,.TRUE.,ilevel)

      CALL collct_setvalue_matsca(rcollection,'CCNL_MATB1',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB1,.TRUE.,ilevel)

      CALL collct_setvalue_matsca(rcollection,'CCNL_MATB2',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB2,.TRUE.,ilevel)

      CALL collct_setvalue_matsca(rcollection,'CCNL_MATMASS',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixMass,.TRUE.,ilevel)

      CALL collct_setvalue_mat(rcollection,'CCNL_MATPREC',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixPreconditioner,&
          .TRUE.,ilevel)
      
      CALL collct_setvalue_matsca(rcollection,'CCNL_MATA12',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixVelocityCoupling12,&
          .TRUE.,ilevel)

      CALL collct_setvalue_matsca(rcollection,'CCNL_MATA21',&
          rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixVelocityCoupling21,&
          .TRUE.,ilevel)
    END DO
    
    ! Save auxiliary variables for the nonlinear iteration
    CALL collct_setvalue_real(rcollection,'CCNL_RESINIT1',&
        rnonlinearIteration%DresidualInit(1),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_RESINIT2',&
        rnonlinearIteration%DresidualInit(2),.TRUE.)

    CALL collct_setvalue_real(rcollection,'CCNL_RESOLD1',&
        rnonlinearIteration%DresidualOld(1),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_RESOLD2',&
        rnonlinearIteration%DresidualOld(2),.TRUE.)

    CALL collct_setvalue_real(rcollection,'CCNL_RESDEL1',&
        rnonlinearIteration%DresidualCorr(1),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_RESDEL2',&
        rnonlinearIteration%DresidualCorr(2),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_RESDEL3',&
        rnonlinearIteration%DresidualCorr(3),.TRUE.)

    ! Save convergence criteria of the nonlinear iteration
    CALL collct_setvalue_real(rcollection,'CCNL_EPSNL1',&
        rnonlinearIteration%DepsNL(1),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_EPSNL2',&
        rnonlinearIteration%DepsNL(2),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_EPSNL3',&
        rnonlinearIteration%DepsNL(3),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_EPSNL4',&
        rnonlinearIteration%DepsNL(4),.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_EPSNL5',&
        rnonlinearIteration%DepsNL(5),.TRUE.)

    ! Other statistical data
    CALL collct_setvalue_real(rcollection,'CCNL_OMEGANL',&
        rnonlinearIteration%domegaNL,.TRUE.)
    CALL collct_setvalue_real(rcollection,'CCNL_RHOMG',&
        rnonlinearIteration%drhoLinearSolver,.TRUE.)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_restoreNonlinearLoop (rnonlinearIteration,rcollection)
  
!<description>
  ! Restores the values of rnonlinearIteration from the given collection
  ! rcollection.
!</description>

!<input>
  ! The nonlinear iteration structure that should be reconstructed.
  TYPE(t_ccNonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</input>

!<output>
  ! The collection structure that contains the data in rnonlinearIteration.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</output>

!</subroutine>
 
    ! local variables
    INTEGER :: ilevel

    ! Restore the structure members

    rnonlinearIteration%domegaMin = collct_getvalue_real(rcollection,&
        'CCNL_OMEGAMIN')
    
    rnonlinearIteration%domegaMax = collct_getvalue_real(rcollection,&
        'CCNL_OMEGAMAX')
    
    rnonlinearIteration%p_rrhs => collct_getvalue_vec(rcollection,&
        'CCNL_RHS')
        
    rnonlinearIteration%p_rsolution => collct_getvalue_vec(rcollection,&
        'CCNL_SOLUTION')

    rnonlinearIteration%dalpha = collct_getvalue_real(rcollection,'CCNL_ALPHA')
    rnonlinearIteration%dtheta = collct_getvalue_real(rcollection,'CCNL_THETA')
    rnonlinearIteration%dgamma = collct_getvalue_real(rcollection,'CCNL_GAMMA')
        
    ! Get information about the preconditioner
    rnonlinearIteration%rpreconditioner%ctypePreconditioning = &
        collct_getvalue_int(rcollection,'CCNL_PRECONDITIONER')

    ! Which preconditioner is chosen?
    SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Get temporary vectors.
      rnonlinearIteration%rpreconditioner%p_rtempVectorSc => &
          collct_getvalue_vecsca(rcollection,'CCNL_RTEMPSCALAR')
      
      rnonlinearIteration%rpreconditioner%p_rtempVectorSc2 => &
          collct_getvalue_vecsca(rcollection,'CCNL_RTEMP2SCALAR')
          
      ! Interlevel projection structure
      rnonlinearIteration%rpreconditioner%p_rprojection => & 
          collct_getvalue_ilvp(rcollection,'CCNL_ILVPROJECTION')

      ! Restore the filter chain
      rnonlinearIteration%rpreconditioner%p_RfilterChain => &
          collct_getvalue_fchn(rcollection,'CCNL_FILTERCHAIN')
      
      ! Put the prepared solver node to the collection for later use.
      ! Remember, it's our preconditioner we need during the nonlinear
      ! iteration!
      rnonlinearIteration%rpreconditioner%p_rsolverNode => &
          collct_getvalue_linsol(rcollection,'CCNL_LINSOLVER')
      
      ! Get the Newton parameters
      rnonlinearIteration%rpreconditioner%radaptiveNewton%nminFixPointIterations = &
          collct_getvalue_int(rcollection,'CCNL_NEWTONMINFP')
      
      rnonlinearIteration%rpreconditioner%radaptiveNewton%nmaxFixPointIterations = &
          collct_getvalue_int(rcollection,'CCNL_NEWTONMAXFP')
      
      rnonlinearIteration%rpreconditioner%radaptiveNewton%depsAbsNewton = &
          collct_getvalue_real(rcollection,'CCNL_NEWTONEPSABS')
      
      rnonlinearIteration%rpreconditioner%radaptiveNewton%depsRelNewton = &
          collct_getvalue_real(rcollection,'CCNL_NEWTONEPSREL')
      
    END SELECT

    ! Reconstruct the final-assembly structure
    rnonlinearIteration%rfinalAssembly%iBmatricesTransposed = &
        collct_getvalue_int(rcollection,'CCNL_IBMATTRANSPOSED')

    rnonlinearIteration%rfinalAssembly%iadaptiveMatrices = &
        collct_getvalue_int(rcollection,'CCNL_IADAPTIVEMATRIX')

    rnonlinearIteration%rfinalAssembly%dadMatThreshold = &
        collct_getvalue_real(rcollection,'CCNL_DADMATTHRESHOLD')

    ! Get information about the levels
    rnonlinearIteration%NLMIN = &
        collct_getvalue_int(rcollection,'CCNL_NLMIN')

    rnonlinearIteration%NLMAX = &
        collct_getvalue_int(rcollection,'CCNL_NLMAX')

    ALLOCATE(rnonlinearIteration%RcoreEquation( &
        rnonlinearIteration%NLMIN:rnonlinearIteration%NLMAX))
    DO ilevel = rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
    
      ! Get the structure members on every level
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrix => &
          collct_getvalue_mat(rcollection,'CCNL_SYSTEMMAT',ilevel)
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixStokes => &
          collct_getvalue_matsca(rcollection,'CCNL_STOKES',ilevel)
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB1 => &
          collct_getvalue_matsca(rcollection,'CCNL_MATB1',ilevel)
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB2 => &
          collct_getvalue_matsca(rcollection,'CCNL_MATB2',ilevel)
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixMass => &
          collct_getvalue_matsca(rcollection,'CCNL_MATMASS',ilevel)
          
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixPreconditioner => &
          collct_getvalue_mat(rcollection,'CCNL_MATPrec',ilevel)

      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixVelocityCoupling12 => &
          collct_getvalue_matsca(rcollection,'CCNL_MATA12',ilevel)
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixVelocityCoupling21 => &
          collct_getvalue_matsca(rcollection,'CCNL_MATA21',ilevel)
    END DO
    
    ! Reconstruct residual information
    rnonlinearIteration%DresidualInit(1) = &
        collct_getvalue_real(rcollection,'CCNL_RESINIT1')
    rnonlinearIteration%DresidualInit(2) = &
        collct_getvalue_real(rcollection,'CCNL_RESINIT2')
    rnonlinearIteration%DresidualOld(1) = &
        collct_getvalue_real(rcollection,'CCNL_RESOLD1')
    rnonlinearIteration%DresidualOld(2) = &
        collct_getvalue_real(rcollection,'CCNL_RESOLD2')
    rnonlinearIteration%DresidualCorr(1) = &
        collct_getvalue_real(rcollection,'CCNL_RESDEL1')
    rnonlinearIteration%DresidualCorr(2) = &
        collct_getvalue_real(rcollection,'CCNL_RESDEL2')
    rnonlinearIteration%DresidualCorr(3) = &
        collct_getvalue_real(rcollection,'CCNL_RESDEL3')

    ! Convergence criteria of the nonlinear iteration
    rnonlinearIteration%DepsNL(1) = &
        collct_getvalue_real(rcollection,'CCNL_EPSNL1')
    rnonlinearIteration%DepsNL(2) = &
        collct_getvalue_real(rcollection,'CCNL_EPSNL2')
    rnonlinearIteration%DepsNL(3) = &
        collct_getvalue_real(rcollection,'CCNL_EPSNL3')
    rnonlinearIteration%DepsNL(4) = &
        collct_getvalue_real(rcollection,'CCNL_EPSNL4')
    rnonlinearIteration%DepsNL(5) = &
        collct_getvalue_real(rcollection,'CCNL_EPSNL5')
    
    ! Other statistical data
    rnonlinearIteration%domegaNL = &
        collct_getvalue_real(rcollection,'CCNL_OMEGANL')
    rnonlinearIteration%drhoLinearSolver = &
        collct_getvalue_real(rcollection,'CCNL_RHOMG')
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_removeNonlinearLoop (rnonlinearIteration,rcollection)
  
!<description>
  ! Removes the nonlinear iteration structure from the collection rcollection
  ! if it was saved there by c2d2_saveNonlinearLoop.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that should be cleaned up.
  TYPE(t_ccNonlinearIteration), INTENT(INOUT) :: rnonlinearIteration

  ! The collection structure which is to be cleaned up if
  ! the nonlinear iteratino structure is saved to that by c2d2_saveNonlinearLoop.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</inputoutput>

!</subroutine>

    INTEGER :: itypePreconditioner,ilevel

    ! Delete statistical data
    CALL collct_deletevalue(rcollection,'CCNL_OMEGANL')
    CALL collct_deletevalue(rcollection,'CCNL_RHOMG')

    ! Delete min./max. damping parameters from the collection
    CALL collct_deletevalue (rcollection,'CCNL_OMEGAMAX')
    CALL collct_deletevalue (rcollection,'CCNL_OMEGAMIN')

    ! Delete parameters that specify the core equation
    CALL collct_deletevalue (rcollection,'CCNL_ALPHA')
    CALL collct_deletevalue (rcollection,'CCNL_THETA')
    CALL collct_deletevalue (rcollection,'CCNL_GAMMA')

    ! Delete solution/RHS from the collection
    CALL collct_deletevalue (rcollection,'CCNL_RHS')
    CALL collct_deletevalue (rcollection,'CCNL_SOLUTION')
    
    ! Get information about the preconditioner
    itypePreconditioner = &
        collct_getvalue_int(rcollection,'CCNL_PRECONDITIONER')

    ! Which preconditioner is chosen?
    SELECT CASE (itypePreconditioner)
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)

      ! Remove the solver node from the collection - not needed anymore there
      CALL collct_deletevalue(rcollection,'CCNL_LINSOLVER')
      
      ! Remove the temporary vector from the collection
      CALL collct_deletevalue(rcollection,'CCNL_RTEMPSCALAR')
      CALL collct_deletevalue(rcollection,'CCNL_RTEMP2SCALAR')
      
      ! Remove the interlevel projection structure
      CALL collct_deletevalue(rcollection,'CCNL_ILVPROJECTION')
      
      ! Remove the filter chain
      CALL collct_deletevalue(rcollection,'CCNL_FILTERCHAIN')
      
      ! Remove parameters for adaptive matrix generation
      CALL collct_deletevalue(rcollection,'CCNL_IBMATTRANSPOSED')
      CALL collct_deletevalue(rcollection,'CCNL_IADAPTIVEMATRIX')
      
      ! Remove Newton parameters
      CALL collct_deletevalue(rcollection,'CCNL_NEWTONMINFP')
      CALL collct_deletevalue(rcollection,'CCNL_NEWTONMAXFP')
      CALL collct_deletevalue(rcollection,'CCNL_NEWTONEPSABS')
      CALL collct_deletevalue(rcollection,'CCNL_NEWTONEPSREL')
      
    END SELECT
  
    CALL collct_deletevalue(rcollection,'CCNL_PRECONDITIONER')
    
    ! Release the final-assembly structure
    CALL collct_deletevalue(rcollection,'CCNL_IBMATTRANSPOSED')
    CALL collct_deletevalue(rcollection,'CCNL_IADAPTIVEMATRIX')
    CALL collct_deletevalue(rcollection,'CCNL_DADMATTHRESHOLD')
    CALL collct_deletevalue(rcollection,'CCNL_INEUMANN')
    
    ! Delete residual information
    CALL collct_deletevalue (rcollection,'CCNL_RESINIT1')
    CALL collct_deletevalue (rcollection,'CCNL_RESINIT2')
    CALL collct_deletevalue (rcollection,'CCNL_RESOLD1')
    CALL collct_deletevalue (rcollection,'CCNL_RESOLD2')
    CALL collct_deletevalue (rcollection,'CCNL_RESPREC1')
    CALL collct_deletevalue (rcollection,'CCNL_RESPREC2')
    CALL collct_deletevalue (rcollection,'CCNL_RESPREC3')
    CALL collct_deletevalue (rcollection,'CCNL_RESDEL1')
    CALL collct_deletevalue (rcollection,'CCNL_RESDEL2')
    CALL collct_deletevalue (rcollection,'CCNL_RESDEL3')
    
    ! Delete convergence criteria of nonlinear iteration
    CALL collct_deletevalue (rcollection,'CCNL_EPSNL1')
    CALL collct_deletevalue (rcollection,'CCNL_EPSNL2')
    CALL collct_deletevalue (rcollection,'CCNL_EPSNL3')
    CALL collct_deletevalue (rcollection,'CCNL_EPSNL4')
    CALL collct_deletevalue (rcollection,'CCNL_EPSNL5')

    ! Release the level-array in the nonlinear iteration structure
    CALL collct_deletevalue(rcollection,'CCNL_NLMAX')
    CALL collct_deletevalue(rcollection,'CCNL_NLMIN')
    
    DO ilevel = rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
      CALL collct_deletevalue (rcollection,'CCNL_SYSTEMMAT',ilevel)
      CALL collct_deletevalue (rcollection,'CCNL_STOKES',ilevel)
      CALL collct_deletevalue (rcollection,'CCNL_MATB1',ilevel)
      CALL collct_deletevalue (rcollection,'CCNL_MATB2',ilevel)
      CALL collct_deletevalue (rcollection,'CCNL_MATMASS',ilevel)
      CALL collct_deletevalue (rcollection,'CCNL_MATA12',ilevel)
      CALL collct_deletevalue (rcollection,'CCNL_MATA21',ilevel)
      CALL collct_deletevalue (rcollection,'CCNL_MATPREC',ilevel)
    END DO
      
  END SUBROUTINE

  ! ***************************************************************************
  ! Routines to set up matrices and (nonlinear) defect vectors.
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_assembleLinearisedMatrices (rnonlinearIteration,rcollection,&
      blevelHirarchy,bboundaryConditions,bboundaryConditionsNonlin,&
      bassemblePreconditioner,bassembleNewton,rx)

  USE linearsystemblock
  USE collection
  
!<description>
  ! This routine assembles the linearised nonlinear system matrices.
  ! Pointers to these matrices must have been attached the rnonlinearIteration
  ! by a previous c2d2_setupCoreEquation.
  !
  ! The system matrices are assembled into the p_rmatrix variable in the
  ! nonlinear iteration structure.
!</description>

!<inputoutput>
  ! Reference to a collection structure that contains all parameters of the
  ! discretisation (for nonlinearity, etc.).
  TYPE(t_collection), INTENT(INOUT)                :: rcollection

  ! Nonlinear iteration structure where to write the linearised system matrix to.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!<input>
  ! Whether to initialise the whole level hirarchy.
  ! TRUE  = calculate the system matrices on all levels.
  ! FALSE = calculate the system matrix only on the maximum level where to solve
  !         the problem
  LOGICAL, INTENT(IN) :: blevelHirarchy
  
  ! TRUE  = include (linear) boundary conditions into the system 
  !   matrix/matrices after assembly
  ! FASLE = don't incorporate any boundary conditions
  LOGICAL, INTENT(IN) :: bboundaryConditions

  ! TRUE  = include nonlinear boundary conditions into the system 
  !   matrix/matrices after assembly
  ! FASLE = don't incorporate any boundary conditions
  LOGICAL, INTENT(IN) :: bboundaryConditionsNonlin

  ! TRUE  = Assemble the preconditioner matrices p_rmatrixPreconditioner
  !         on every level in rnonlinearIteration.
  ! FALSE = Assemble the linearised system matrices p_rmatrix
  !         on every level in rnonlinearIteration.
  LOGICAL, INTENT(IN) :: bassemblePreconditioner
  
  ! Applies only if bassemblePreconditioner=TRUE:
  ! TRUE  = Assemble the Newton preconditioner.
  ! FALSE = Assemble the standard defect correction preconditioner
  !         (i.e. the linearised system matrix).
  LOGICAL, INTENT(IN) :: bassembleNewton
  
  ! OPTIONAL: Current iteration vector. Must be specified if the nonlinearity
  ! should be discretised.
  TYPE(t_vectorBlock), INTENT(IN), TARGET, OPTIONAL       :: rx
 
!</input>

!</subroutine>

  ! local variables
  INTEGER :: iupwind

  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  TYPE(t_matrixBlock), POINTER :: p_rmatrix,p_rmatrixFine
  TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes,p_rmatrixMass
  TYPE(t_vectorScalar), POINTER :: p_rvectorTemp
  INTEGER :: NLMAX,NLMIN, ilev
  TYPE(t_convUpwind) :: rupwind
  TYPE(t_convStreamlineDiffusion) :: rstreamlineDiffusion
  TYPE(t_jumpStabilisation) :: rjumpStabil
  TYPE(t_vectorBlock), POINTER :: p_rvectorFine,p_rvectorCoarse
  
  TYPE(t_interlevelProjectionBlock), POINTER :: p_rprojection

  ! A filter chain for the linear solver
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
  ! DEBUG!!!
!    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    NLMAX = collct_getvalue_int (p_rcollection,'NLMAX')

    ! Get minimum and maximum level from the collection
    NLMIN = rnonlinearIteration%NLMIN
    NLMAX = rnonlinearIteration%NLMAX
    
    IF (.NOT. blevelHirarchy) THEN
      ! Discretise only maximum level
      NLMIN = NLMAX
    END IF
    
    ! Get the interlevel projection structure and the temporary vector
    ! from the collection.
    ! Our 'parent' prepared there how to interpolate the solution on the
    ! fine grid to coarser grids.
    p_rprojection => rnonlinearIteration%rpreconditioner%p_rprojection
    p_rvectorTemp => rnonlinearIteration%rpreconditioner%p_rtempVectorSc

    p_RfilterChain => rnonlinearIteration%rpreconditioner%p_RfilterChain

    ! If we discretise Navier-Stokes, we have to set up a stabilisation
    ! structure dor upwind / streamline diffusion / ...
    
    IF (rnonlinearIteration%dgamma .NE. 0.0_DP) THEN
    
      iupwind = collct_getvalue_int (rcollection,'IUPWIND')
    
      SELECT CASE (iupwind)
      CASE (0)
        ! Set up the SD structure for the creation of the defect.
        ! There's not much to do, only initialise the viscosity...
        rstreamlineDiffusion%dnu = collct_getvalue_real (rcollection,'NU')
        
        ! Set stabilisation parameter
        rstreamlineDiffusion%dupsam = collct_getvalue_real (rcollection,'UPSAM')
        
        ! Matrix weight
        rstreamlineDiffusion%dtheta = rnonlinearIteration%dgamma

      CASE (1)
        ! Set up the upwind structure for the creation of the defect.
        ! There's not much to do, only initialise the viscosity...
        rupwind%dnu = collct_getvalue_real (rcollection,'NU')
        
        ! Set stabilisation parameter
        rupwind%dupsam = collct_getvalue_real (rcollection,'UPSAM')

        ! Matrix weight
        rupwind%dtheta = rnonlinearIteration%dgamma

      CASE (2)
        ! Set up the jump stabilisation structure for the creation of the defect.
        ! There's not much to do, only initialise the viscosity...
        rjumpStabil%dnu = collct_getvalue_real (rcollection,'NU')
        
        ! Set stabilisation parameter
        rjumpStabil%dgammastar = collct_getvalue_real (rcollection,'UPSAM')
        rjumpStabil%dgamma = rjumpStabil%dgammastar
        
        ! Set up the SD structure for the creation of the defect.
        ! Initialise the viscosity...
        rstreamlineDiffusion%dnu = rjumpStabil%dnu
        
        ! Set upsam=0 to obtain a central-difference-type discretisation
        rstreamlineDiffusion%dupsam = 0.0_DP

        ! Matrix weight
        rstreamlineDiffusion%dtheta = rnonlinearIteration%dgamma
        rjumpStabil%dtheta = rnonlinearIteration%dgamma

      CASE DEFAULT
        PRINT *,'Don''t know how to set up nonlinearity!?!'
        STOP
      
      END SELECT
      
    ELSE
      iupwind = collct_getvalue_int (rcollection,'IUPWIND')
    
      SELECT CASE (iupwind)
      CASE (2)
        ! Jump stabilisation; can also be used with Stokes.
        !
        ! Set up the jump stabilisation structure for the creation of the defect.
        ! There's not much to do, only initialise the viscosity...
        rjumpStabil%dnu = collct_getvalue_real (rcollection,'NU')
        
        ! Set stabilisation parameter
        rjumpStabil%dgammastar = collct_getvalue_real (rcollection,'UPSAM')
        rjumpStabil%dgamma = rjumpStabil%dgammastar
      END SELECT      
      
    END IF
    
    ! On all levels, we have to set up the nonlinear system matrix,
    ! so that the linear solver can be applied to it.

    NULLIFY(p_rmatrix)
    DO ilev=NLMAX,NLMIN,-1
    
      ! Get the system matrix and the Stokes matrix
      p_rmatrixFine => p_rmatrix
      ! What should we assemble? Linearised matrix or preconditioner?
      IF (bassemblePreconditioner) THEN
        p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner
      ELSE
        p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrix
      END IF
      p_rmatrixStokes => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixStokes
      p_rmatrixMass => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixMass
      
      ! On the highest level, we use rx as solution to build the nonlinear
      ! matrix. On lower levels, we have to create a solution
      ! on that level from a fine-grid solution before we can use
      ! it to build the matrix!
      IF (ilev .EQ. NLMAX) THEN
        p_rvectorCoarse => rx
      ELSE
        ! Get the temporary vector on level i. Will receive the solution
        ! vector on that level. 
        p_rvectorCoarse => collct_getvalue_vec (rcollection,PAR_TEMPVEC,ilev)
        
        ! Get the solution vector on level i+1. This is either the temporary
        ! vector on that level, or the solution vector on the maximum level.
        IF (ilev .LT. NLMAX-1) THEN
          p_rvectorFine => collct_getvalue_vec (rcollection,PAR_TEMPVEC,ilev+1)
        ELSE
          p_rvectorFine => rx
        END IF

        ! Interpolate the solution from the finer grid to the coarser grid.
        ! The interpolation is configured in the interlevel projection
        ! structure we got from the collection.
        CALL mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                         p_rvectorFine,p_rvectorTemp)

        ! Apply the filter chain to the temp vector.
        ! This implements the boundary conditions that are attached to it.
        ! NOTE: Deactivated for standard CC2D compatibility -- and because
        ! it has to be checked whether the correct boundary conditions
        ! are attached to that vector!
        ! CALL filter_applyFilterChainVec (p_rvectorCoarse, p_RfilterChain)

      END IF
      
      ! Should we discretise the Navier-Stokes nonlinearity?
      IF (rnonlinearIteration%dgamma .NE. 0.0_DP) THEN
      
        IF (.NOT. PRESENT(rx)) THEN
          PRINT *,'c2d2_assembleLinearisedMatrices: rx not specified!'
          STOP
        END IF
      
        ! Which type of stabilisation/strategy for setting up the nonlinearity
        ! do we use?
        SELECT CASE (iupwind)
        CASE (0)
      
          ! The system matrix looks like:
          !   (  A    0   B1 )
          !   (  0    A   B2 )
          !   ( B1^T B2^T 0  )
          !
          ! The A-matrix consists of Mass+Stokes+Convection.
          ! We build them separately and add together.
          !
          ! So at first, initialise the A-matrix with the Mass+Stokes contribution.
          ! We ignore the structure and simply overwrite the content of the
          ! system submatrices with the Stokes matrix.
          
          IF ((rnonlinearIteration%dalpha .NE. 0.0_DP) .OR. &
              (rnonlinearIteration%dtheta .NE. 0.0_DP)) THEN
          
            ! Copy the Stokes matrix, overwrite p_rmatrix in any case.
            CALL lsyssc_duplicateMatrix (p_rmatrixStokes,p_rmatrix%RmatrixBlock(1,1),&
                                         LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
                                        
            ! Mass matrix?
            IF (rnonlinearIteration%dalpha .NE. 0.0_DP) THEN
              CALL lsyssc_matrixLinearComb (&
                  p_rmatrixMass,rnonlinearIteration%dalpha,&
                  p_rmatrix%RmatrixBlock(1,1),rnonlinearIteration%dtheta,&
                  p_rmatrix%RmatrixBlock(1,1),&
                  .FALSE.,.FALSE.,.TRUE.,.TRUE.)
            ELSE
              ! In this case, we may have to scale the Stokes matrix according to 
              ! theta; note that if the mass matrix is involved, this is done
              ! implicitely.
              IF (rnonlinearIteration%dtheta .NE. 1.0_DP) THEN
                CALL lsyssc_scaleMatrix (p_rmatrix%RmatrixBlock(1,1),&
                    rnonlinearIteration%dtheta)
              END IF
            END IF
          END IF
          
          ! Assemble either the standard linearised matrix or the Newton matrix.

          IF ((.NOT. bassemblePreconditioner) .OR. (.NOT. bassembleNewton)) THEN

            ! If A22=A11, that's all for the diagonals of the matrix.
            ! If A22 has a separate content in memory, we copy the content
            ! of A11 to A22 so that they are the same.
            IF (.NOT. lsyssc_isMatrixContentShared( &
              p_rmatrix%RmatrixBlock(1,1),p_rmatrix%RmatrixBlock(2,2))) THEN
              CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                  p_rmatrix%RmatrixBlock(2,2),&
                  LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
            END IF

            ! Make sure not to calculate the Newton part.
            rstreamlineDiffusion%dnewton = 0.0_DP

            ! Call the SD method to calculate the nonlinear part in A11 and A22.
            CALL conv_streamlineDiffusionBlk2d (&
                                p_rvectorCoarse, p_rvectorCoarse, &
                                1.0_DP, 0.0_DP,&
                                rstreamlineDiffusion, CONV_MODMATRIX, &
                                p_rmatrix)
                                
            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there's content at all),
            ! these matrices are ignored then by the kernel.
            
            p_rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            p_rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
                                
            ! Alternative iplementation with the scalar version of SD:
            !$
            ! Call the SD method to calculate the nonlinear matrix.
            !CALL conv_streamlineDiffusion2d (&
            !                    p_rvectorCoarse, p_rvectorCoarse, &
            !                    1.0_DP, 0.0_DP,&
            !                    rstreamlineDiffusion, CONV_MODMATRIX, &
            !                    p_rmatrix%RmatrixBlock(1,1))
            !
            ! If A22=A11, that's all for the diagonals of the matrix.
            ! If A22 has a separate content in memory, we copy the content
            ! of A11 to A22 so that they are the same.
            !IF (.NOT. lsyssc_isMatrixContentShared( &
            !  p_rmatrix%RmatrixBlock(1,1),p_rmatrix%RmatrixBlock(2,2))) THEN
            !  CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
            !      p_rmatrix%RmatrixBlock(2,2),&
            !      LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
            !END IF
            
          ELSE
            
            ! Copy A11 to A22 and clear A12 and A21. Then, assemble the
            ! convectve operator + Newton into A11, A12, A21 and A22.
            CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                p_rmatrix%RmatrixBlock(2,2),&
                LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
            
            CALL lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(1,2))
            CALL lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they aren't.
            p_rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            p_rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP

            ! Call the extended SD method to set up the velocity part
            ! including the Newton matrix.
            ! Configure SD to calculate Newton.
            rstreamlineDiffusion%dnewton = 1.0_DP
            CALL conv_streamlineDiffusionBlk2d (&
                                p_rvectorCoarse, p_rvectorCoarse, &
                                1.0_DP, 0.0_DP,&
                                rstreamlineDiffusion, CONV_MODMATRIX, &
                                p_rmatrix)

!            ! Alternative implementation of Newton with the trilinear form:
!            ! 
!            ! First assemble the nonlinearity into A11,
!            CALL conv_streamlineDiffusion2d (&
!                                p_rvectorCoarse, p_rvectorCoarse, &
!                                1.0_DP, 0.0_DP,&
!                                rstreamlineDiffusion, CONV_MODMATRIX, &
!                                p_rmatrix%RmatrixBlock(1,1))
!            ! Copy A11 to A22
!            CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
!                p_rmatrix%RmatrixBlock(2,2),&
!                LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
!
!            ! Add the Newton matrix to A11, A12, A21 and A22 using the
!            ! trilinear form.
!
!            TYPE(t_trilinearForm) :: rform
!
!            rform%itermCount = 1
!            rform%BconstantCoeff = .TRUE.          
!            rform%BallCoeffConstant = .TRUE.      
!            rform%Dcoefficients(1)  = rnonlinearIteration%dgamma
!            rform%Idescriptors(2,1) = DER_FUNC
!            rform%Idescriptors(3,1) = DER_FUNC
!
!            CALL lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(1,2))
!            CALL lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(2,1))
!
!            rform%Idescriptors(1,1) = DER_DERIV_X
!            CALL trilf_buildMatrix9d_conf2 (rform,.FALSE., &
!                p_rmatrix%RmatrixBlock(1,1),p_rvectorCoarse%RvectorBlock(1))
!            
!            rform%Idescriptors(1,1) = DER_DERIV_Y
!            CALL trilf_buildMatrix9d_conf2 (rform,.FALSE., &
!                p_rmatrix%RmatrixBlock(1,2),p_rvectorCoarse%RvectorBlock(1))
!
!            rform%Idescriptors(1,1) = DER_DERIV_X
!            CALL trilf_buildMatrix9d_conf2 (rform,.FALSE., &
!                p_rmatrix%RmatrixBlock(2,1),p_rvectorCoarse%RvectorBlock(2))
!
!            rform%Idescriptors(1,1) = DER_DERIV_Y
!            CALL trilf_buildMatrix9d_conf2 (rform,.FALSE., &
!                p_rmatrix%RmatrixBlock(2,2),p_rvectorCoarse%RvectorBlock(2))
                
          END IF
          
        CASE (1)
      
          ! The system matrix looks like:
          !   (  A    0   B1 )
          !   (  0    A   B2 )
          !   ( B1^T B2^T 0  )
          !
          ! The A-matrix consists of MassStokes+Convection.
          ! We build them separately and add together.
          !
          ! So at first, initialise the A-matrix with the Stokes contribution.
          ! We ignore the structure and simply overwrite the content of the
          ! system submatrices with the Stokes matrix.

          IF ((rnonlinearIteration%dalpha .NE. 0.0_DP) .OR. &
              (rnonlinearIteration%dtheta .NE. 0.0_DP)) THEN

            ! Copy the Stokes matrix, overwrite p_rmatrix in any case.
            CALL lsyssc_duplicateMatrix (p_rmatrixStokes,p_rmatrix%RmatrixBlock(1,1),&
                                        LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
            ! Mass matrix?
            IF (rnonlinearIteration%dalpha .NE. 0.0_DP) THEN
              CALL lsyssc_matrixLinearComb (&
                  p_rmatrixMass,rnonlinearIteration%dalpha,&
                  p_rmatrix%RmatrixBlock(1,1),rnonlinearIteration%dtheta,&
                  p_rmatrix%RmatrixBlock(1,1),&
                  .FALSE.,.FALSE.,.TRUE.,.TRUE.)
            ELSE
              ! In this case, we may have to scale the Stokes matrix according to 
              ! theta; note that if the mass matrix is involved, this is done
              ! implicitely.
              IF (rnonlinearIteration%dtheta .NE. 1.0_DP) THEN
                CALL lsyssc_scaleMatrix (p_rmatrix%RmatrixBlock(1,1),&
                    rnonlinearIteration%dtheta)
              END IF
            END IF
            
          END IF

          ! Call the upwind method to calculate the nonlinear matrix.
          CALL conv_upwind2d (p_rvectorCoarse, p_rvectorCoarse, &
                              1.0_DP, 0.0_DP,&
                              rupwind, CONV_MODMATRIX, &
                              p_rmatrix%RmatrixBlock(1,1))      
                              
          ! If A22=A11, that's all for the diagonals of the matrix.
          ! If A22 has a separate content in memory, we copy the content
          ! of A11 to A22 so that they are the same.
          IF (.NOT. lsyssc_isMatrixContentShared(p_rmatrix%RmatrixBlock(2,2))) THEN
            CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                                        p_rmatrix%RmatrixBlock(2,2),&
                                        LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
          END IF

        CASE (2)

          ! The system matrix looks like:
          !   (  A    0   B1 )
          !   (  0    A   B2 )
          !   ( B1^T B2^T 0  )
          !
          ! The A-matrix consists of Stokes+Convection.
          ! We build them separately and add together.
          !
          ! At first, initialise the A-matrix with the Stokes contribution.
          ! We ignore the structure and simply overwrite the content of the
          ! system submatrices with the Stokes matrix.
          IF ((rnonlinearIteration%dalpha .NE. 0.0_DP) .OR. &
              (rnonlinearIteration%dtheta .NE. 0.0_DP)) THEN

            ! Copy the Stokes matrix, overwrite p_rmatrix in any case.
            CALL lsyssc_duplicateMatrix (p_rmatrixStokes,p_rmatrix%RmatrixBlock(1,1),&
                                        LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
            ! Mass matrix?
            IF (rnonlinearIteration%dalpha .NE. 0.0_DP) THEN
              CALL lsyssc_matrixLinearComb (&
                  p_rmatrixMass,rnonlinearIteration%dalpha,&
                  p_rmatrix%RmatrixBlock(1,1),rnonlinearIteration%dtheta,&
                  p_rmatrix%RmatrixBlock(1,1),&
                  .FALSE.,.FALSE.,.TRUE.,.TRUE.)
            ELSE
              ! In this case, we may have to scale the Stokes matrix according to 
              ! theta; note that if the mass matrix is involved, this is done
              ! implicitely.
              IF (rnonlinearIteration%dtheta .NE. 1.0_DP) THEN
                CALL lsyssc_scaleMatrix (p_rmatrix%RmatrixBlock(1,1),&
                    rnonlinearIteration%dtheta)
              END IF
            END IF
            
          END IF

          ! Call the SD method to calculate the nonlinear matrix.
          ! We use the streamline-diffusion discretisation routine
          ! which uses a central-difference-like discretisation.
          ! As the parameter rstreamlineDiffusion%dbeta is =0 by default, the
          ! Stokes part is not calculated in this routine.
          CALL conv_streamlineDiffusion2d (&
                              p_rvectorCoarse, p_rvectorCoarse, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              p_rmatrix%RmatrixBlock(1,1))   
                              
          ! Call the jump stabilisation technique to stabilise that stuff.   
          CALL conv_jumpStabilisation2d (&
                              p_rvectorCoarse, p_rvectorCoarse, 1.0_DP, 0.0_DP,&
                              rjumpStabil, CONV_MODMATRIX, &
                              p_rmatrix%RmatrixBlock(1,1))   

          ! If A22=A11, that's all for the diagonals of the matrix.
          ! If A22 has a separate content in memory, we copy the content
          ! of A11 to A22 so that they are the same.
          IF (.NOT. lsyssc_isMatrixContentShared(p_rmatrix%RmatrixBlock(2,2))) THEN
            CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                                        p_rmatrix%RmatrixBlock(2,2),&
                                        LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
          END IF

        CASE DEFAULT
          PRINT *,'Don''t know how to set up nonlinearity!?!'
          STOP
          
        END SELECT
      
      ELSE
        ! The system matrix looks like:
        !   (  A    0   B1 )
        !   (  0    A   B2 )
        !   ( B1^T B2^T 0  )
        !
        ! The A-matrix is a simple Mass+Stokes.
        !
        ! Copy the Stokes matrix to A and filter it according to the boundary
        ! conditions - this gives then the system matrix.

        IF ((rnonlinearIteration%dalpha .NE. 0.0_DP) .OR. &
            (rnonlinearIteration%dtheta .NE. 0.0_DP)) THEN

          ! Copy the Stokes matrix, overwrite p_rmatrix in any case.
          CALL lsyssc_duplicateMatrix (p_rmatrixStokes,p_rmatrix%RmatrixBlock(1,1),&
                                      LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
          ! Mass matrix?
          IF (rnonlinearIteration%dalpha .NE. 0.0_DP) THEN
            CALL lsyssc_matrixLinearComb (&
                p_rmatrixMass,rnonlinearIteration%dalpha,&
                p_rmatrix%RmatrixBlock(1,1),rnonlinearIteration%dtheta,&
                p_rmatrix%RmatrixBlock(1,1),&
                .FALSE.,.FALSE.,.TRUE.,.TRUE.)
          ELSE
            ! In this case, we may have to scale the Stokes matrix according to 
            ! theta; note that if the mass matrix is involved, this is done
            ! implicitely.
            IF (rnonlinearIteration%dtheta .NE. 1.0_DP) THEN
              CALL lsyssc_scaleMatrix (p_rmatrix%RmatrixBlock(1,1),&
                  rnonlinearIteration%dtheta)
            END IF
          END IF

          SELECT CASE (iupwind)
          CASE (2)
            ! Jump stabilisation; can also be used with Stokes.
            !
            ! Call the jump stabilisation technique to stabilise that stuff.   
            CALL conv_jumpStabilisation2d (&
                                p_rvectorCoarse, p_rvectorCoarse, 1.0_DP, 0.0_DP,&
                                rjumpStabil, CONV_MODMATRIX, &
                                p_rmatrix%RmatrixBlock(1,1))   
          END SELECT
          
        END IF

        ! If A22=A11, that's all.
        ! If A22 has a separate content in memory, we copy the content
        ! of A11 to A22 so that they are the same.
        IF (.NOT. lsyssc_isMatrixContentShared(p_rmatrix%RmatrixBlock(2,2))) THEN
          CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                                      p_rmatrix%RmatrixBlock(2,2),&
                                      LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)
        END IF

      END IF

      ! For the construction of matrices on lower levels, call the matrix
      ! restriction. In case we have a uniform discretisation with Q1~,
      ! iadaptivematrix is <> 0 and so this will rebuild some matrix entries
      ! by a Galerkin approach using constant prolongation/restriction.
      ! This helps to stabilise the solver if there are elements in the
      ! mesh with high aspect ratio.
      IF (ilev .LT. NLMAX) THEN
        CALL mrest_matrixRestrictionEX3Y (p_rmatrixFine%RmatrixBlock(1,1), &
            p_rmatrix%RmatrixBlock(1,1), &
            rnonlinearIteration%rfinalAssembly%iadaptiveMatrices, &
            rnonlinearIteration%rfinalAssembly%dadmatthreshold)
            
        IF (.NOT. lsyssc_isMatrixContentShared(p_rmatrix%RmatrixBlock(2,2))) THEN
          CALL mrest_matrixRestrictionEX3Y (p_rmatrixFine%RmatrixBlock(2,2), &
              p_rmatrix%RmatrixBlock(1,1), &
              rnonlinearIteration%rfinalAssembly%iadaptiveMatrices, &
              rnonlinearIteration%rfinalAssembly%dadmatthreshold)
        END IF
      END IF
      
      IF (bboundaryConditions) THEN
    
        ! Apply the filter chain to the matrix.
        ! As the filter consists only of an implementation filter for
        ! boundary conditions, this implements the boundary conditions
        ! into the system matrix.
        CALL filter_applyFilterChainMat (p_rmatrix, p_RfilterChain)
        
      END IF
        
      IF (bboundaryConditionsNonlin) THEN

        ! 'Nonlinear' boundary conditions like slip boundary conditions
        ! are not implemented with a filter chain into a matrix.
        ! Call the appropriate matrix filter of 'nonlinear' boundary
        ! conditions manually:
        CALL matfil_discreteNLSlipBC (p_rmatrix,.TRUE.)
        
      END IF
      
    END DO
      
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_assembleConvDiffDefect (rnonlinearIteration,rx,rd,rcollection)

  USE linearsystemblock
  USE collection
  
!<description>
  ! Assembles the convection-diffusion-defect of the core equation.
  ! That means, from the core equation
  !   $alpha*M*u + theta*nu*Laplace*u + gamma*N(u)u + B*p = f1$
  !                                                $B^T*u = f2$
  ! we pick the linear+nonlinear part which is weighted by alpha,
  ! theta and gamma and assemble only that part of the defect:
  !   $d = f1 - alpha*mass - theta*nu*Laplace*u - gamma*N(u)u$
  !
  ! The routine does not implement any boundary conditions!
!</description>

!<input>
  ! Current iteration vector
  TYPE(t_vectorBlock), INTENT(IN),TARGET        :: rx
!</input>
              
!<inputoutput>
  ! Nonlinear iteration structure that specifies the core equation.
  TYPE(t_ccnonlinearIteration), INTENT(IN)      :: rnonlinearIteration

  ! Collection structure of the application. 
  TYPE(t_collection), INTENT(INOUT)             :: rcollection

  ! IN: Right hand side vector of the equation.
  ! OUT: Defect vector b-A(x)x. 
  TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: ilvmax,i
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes,p_rmatrixMass
  TYPE(t_matrixBlock) :: rmatrixTmpBlock
  TYPE(t_convUpwind) :: rupwind
  TYPE(t_convStreamlineDiffusion) :: rstreamlineDiffusion
  TYPE(T_jumpStabilisation) :: rjumpStabil
  
  ! DEBUG!!!
  !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

  ! A filter chain for the linear solver
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
    ! DEBUG!!!
    !CALL lsysbl_getbase_double (rx,p_Ddata)
    !CALL lsysbl_getbase_double (rd,p_Ddata2)
  
    ! Rebuild the nonlinear-iteration structure from the collection
    p_RfilterChain => rnonlinearIteration%rpreconditioner%p_RfilterChain
  
    ! Get minimum/maximum level from the collection
    ilvmax = rnonlinearIteration%NLMAX
    
    ! Get the system and the Stokes matrix on the maximum level
    p_rmatrix => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrix
    p_rmatrixStokes => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixStokes
    p_rmatrixMass => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass

    ! Now, in the first step, we build the linear part of the nonlinear defect:
    !     d_lin = rhs - theta*(-nu * Laplace(.))*solution
    !
    ! rd already contains the rhs.
    IF (rnonlinearIteration%dtheta .NE. 0.0_DP) THEN
      ! Build a temporary 3x3 block matrix rmatrixStokes with Stokes matrices 
      ! on the main diagonal.
      !
      ! (  L    0   B1 )
      ! (  0    L   B2 )
      ! ( B1^T B2^T 0  )
      
      CALL lsysbl_duplicateMatrix (p_rmatrix,rmatrixTmpBlock, &
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsyssc_duplicateMatrix (p_rmatrixStokes,rmatrixTmpBlock%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsyssc_duplicateMatrix (p_rmatrixStokes,rmatrixTmpBlock%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      ! Remove all rows and columns from the 3rd line on. For setting up the defect, we
      ! will treat these separately. For stationary simulations, this is the same as
      ! setting up the defect with the complete matrix as above. For nonstationary
      ! simulations, there is a difference, as the divergence part of the equation
      ! is not scaled with the time step! (see below).
      !
      ! (  L    0   0  )
      ! (  0    L   0  )
      ! (  0    0   0  )
      !
      DO i=NDIM2D+1,rmatrixTmpBlock%ndiagBlocks
        CALL lsysbl_releaseMatrixRow (rmatrixTmpBlock,i)
        CALL lsysbl_releaseMatrixColumn (rmatrixTmpBlock,i)
      END DO
    
      ! Build the defect in the velocity equations
      CALL lsysbl_blockMatVec (rmatrixTmpBlock, rx, rd, &
          -1.0_DP*rnonlinearIteration%dtheta, 1.0_DP)
          
      CALL lsysbl_releaseMatrix (rmatrixTmpBlock)
      
    END IF
      
    ! Let's see, if the weight of the mass matrix is <> 0, subtract that
    ! contribution, too:
    !     d_lin = dlin - dalpha*mass*solution
    IF (rnonlinearIteration%dalpha .NE. 0.0_DP) THEN
      ! Set up a new block mass matrix that is completely zero.
      ! Copy references to our scalar mass matrices to the diagonal blocks (1,1)
      ! and (2,2).
      ! Multiply with the mass matrix and release the block matrix again, as we
      ! don't need it afterwards.
      ! Note that no matrix entries will be copied, so this is a quick operation!
    
      CALL lsysbl_createMatBlockByDiscr (rx%p_rblockDiscretisation,rmatrixTmpBlock)
      CALL lsyssc_duplicateMatrix (p_rmatrixMass,rmatrixTmpBlock%RmatrixBlock(1,1),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
      CALL lsyssc_duplicateMatrix (p_rmatrixMass,rmatrixTmpBlock%RmatrixBlock(2,2),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
                                    
      CALL lsysbl_blockMatVec (rmatrixTmpBlock, rx, rd, &
          -1.0_DP*rnonlinearIteration%dalpha, 1.0_DP)
          
      CALL lsysbl_releaseMatrix (rmatrixTmpBlock)
    END IF
    
    ! Should we discretise the Navier-Stokes nonlinearity?
    IF (rnonlinearIteration%dgamma .NE. 0.0_DP) THEN
    
      ! Which type of stabilisation/strategy for setting up the nonlinearity
      ! do we use?
      SELECT CASE (collct_getvalue_int (rcollection,'IUPWIND'))
      CASE (0)
      
        ! Set up the SD structure for the creation of the defect.
        ! There's not much to do, only initialise the viscosity...
        rstreamlineDiffusion%dnu = collct_getvalue_real (rcollection,'NU')
        
        ! Set stabilisation parameter
        rstreamlineDiffusion%dupsam = collct_getvalue_real (rcollection,'UPSAM')
        
        ! Matrix weight
        rstreamlineDiffusion%dtheta = rnonlinearIteration%dgamma

        ! Call the SD method to calculate the nonlinear defect.
        ! As we calculate only the defect, the matrix is ignored!
        CALL conv_streamlinediffusion2d ( &
                          rx, rx, 1.0_DP, 0.0_DP,&
                          rstreamlineDiffusion, CONV_MODDEFECT, &
                          p_rmatrix%RmatrixBlock(1,1), rx, rd)
                
      CASE (1)
    
        ! Set up the upwind structure for the creation of the defect.
        ! There's not much to do, only initialise the viscosity...
        rupwind%dnu = collct_getvalue_real (rcollection,'NU')
        
        ! Set stabilisation parameter
        rupwind%dupsam = collct_getvalue_real (rcollection,'UPSAM')
        
        ! Matrix weight
        rupwind%dtheta = rnonlinearIteration%dgamma
        
        ! Call the upwind method to calculate the nonlinear defect.
        ! As we calculate only the defect, the matrix is ignored!
        CALL conv_upwind2d (rx, rx, 1.0_DP, 0.0_DP,&
                            rupwind, CONV_MODDEFECT, &
                            p_rmatrix%RmatrixBlock(1,1), rx, rd)      
                
      CASE (2)
      
        ! Set up the SD structure for the creation of the defect.
        ! Initialise the viscosity...
        rstreamlineDiffusion%dnu = collct_getvalue_real (rcollection,'NU')
        
        ! Set stabilisation parameter to 0.0 to get a central difference
        ! like discretisation.
        rstreamlineDiffusion%dupsam = 0.0_DP
        
        ! Matrix weight
        rstreamlineDiffusion%dtheta = rnonlinearIteration%dgamma

        ! Call the SD method to calculate the nonlinear defect.
        ! As we calculate only the defect, the matrix is ignored!
        CALL conv_streamlinediffusion2d ( &
                          rx, rx, 1.0_DP, 0.0_DP,&
                          rstreamlineDiffusion, CONV_MODDEFECT, &
                          p_rmatrix%RmatrixBlock(1,1), rx, rd)
                          
        ! Initialise the Jump stabilisation structure.
        rjumpStabil%dnu = rstreamlineDiffusion%dnu
        
        ! Initialise the jump stabilisation parameter
        rjumpStabil%dgamma = collct_getvalue_real (rcollection,'UPSAM')
        rjumpStabil%dgammastar = rjumpStabil%dgamma 
        
        ! Matrix weight
        rjumpStabil%dtheta = rnonlinearIteration%dgamma

        ! Call the Jump stabilisation to stabilise
        CALL conv_jumpStabilisation2d ( &
                          rx, rx, rnonlinearIteration%dgamma, 0.0_DP,&
                          rjumpStabil, CONV_MODDEFECT, &
                          p_rmatrix%RmatrixBlock(1,1), rx, rd)
                
      CASE DEFAULT
        PRINT *,'Don''t know how to set up nonlinearity!?!'
        STOP
        
      END SELECT
      
    ELSE
      ! Which type of stabilisation/strategy for setting up the nonlinearity
      ! do we use?
      SELECT CASE (collct_getvalue_int (rcollection,'IUPWIND'))
      CASE (2)
        ! Jump stabilisation. Can also be used with Stokes.
        !
        ! Initialise the Jump stabilisation structure.
        rjumpStabil%dnu = rstreamlineDiffusion%dnu
        
        ! Initialise the jump stabilisation parameter
        rjumpStabil%dgamma = collct_getvalue_real (rcollection,'UPSAM')
        rjumpStabil%dgammastar = rjumpStabil%dgamma 
        
        ! Matrix weight
        rjumpStabil%dtheta = rnonlinearIteration%dgamma

        ! Call the Jump stabilisation to stabilise
        CALL conv_jumpStabilisation2d ( &
                          rx, rx, rnonlinearIteration%dgamma, 0.0_DP,&
                          rjumpStabil, CONV_MODDEFECT, &
                          p_rmatrix%RmatrixBlock(1,1), rx, rd)

      END SELECT
      
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_assembleNonlinearDefect (rnonlinearIteration,rx,rd,&
      bboundaryConditions,bboundaryConditionsNonlin,rcollection)

  USE linearsystemblock
  USE collection
  
!<description>
  ! Assembles the nonlinear defect of the core equation.
!</description>

!<input>
  ! Current iteration vector
  TYPE(t_vectorBlock), INTENT(IN),TARGET        :: rx

  ! TRUE  = include (linear) boundary conditions into the system 
  !   matrix/matrices after assembly
  ! FASLE = don't incorporate any boundary conditions
  LOGICAL, INTENT(IN) :: bboundaryConditions

  ! TRUE  = include nonlinear boundary conditions into the system 
  !   matrix/matrices after assembly
  ! FASLE = don't incorporate any boundary conditions
  LOGICAL, INTENT(IN) :: bboundaryConditionsNonlin

!</input>
              
!<inputoutput>
  ! Nonlinear iteration structure that specifies the core equation.
  TYPE(t_ccnonlinearIteration), INTENT(IN)      :: rnonlinearIteration

  ! Collection structure of the application. 
  TYPE(t_collection), INTENT(INOUT)             :: rcollection

  ! IN: Right hand side vector of the equation.
  ! OUT: Defect vector b-A(x)x. 
  TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: ilvmax,i,j
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes
  TYPE(t_matrixBlock) :: rmatrixTmpBlock
  
  ! DEBUG!!!
  !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

  ! A filter chain for the linear solver
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
    ! DEBUG!!!
    !CALL lsysbl_getbase_double (rx,p_Ddata)
    !CALL lsysbl_getbase_double (rd,p_Ddata2)
  
    ! Rebuild the nonlinear-iteration structure from the collection
    p_RfilterChain => rnonlinearIteration%rpreconditioner%p_RfilterChain
  
    ! Get minimum/maximum level from the collection
    ilvmax = rnonlinearIteration%NLMAX
    
    ! Get the system and the Stokes matrix on the maximum level
    p_rmatrix => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrix
    p_rmatrixStokes => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixStokes

    ! Call the "assemble Convection-Diffusion-defect" routine to assemble
    ! the defect that arises from the Convection/Diffusion part and of those
    ! parts that depend on theta, gamma and/or alpha. These parts have a special role
    ! in time-dependent schemes, so the assembly routine for this part
    ! of the defect is separate.
    !
    ! rd already contains the rhs and is modified now:
    CALL c2d2_assembleConvDiffDefect (rnonlinearIteration,rx,rd,rcollection)

    ! To complete the setup of the defect, we have to treat the 'other' parts
    ! of the equation, which usually means that we have to tackle the divergence-
    ! and pressure-parts.
    IF (rnonlinearIteration%dtheta .NE. 0.0_DP) THEN
      ! Set up a block system for the divergence and pressure by duplicating 
      ! the system matrix and removing the Convection/Diffusion parts:
      !
      ! (  0    0   B1 )
      ! (  0    0   B2 )
      ! ( B1^T B2^T 0  )
      !
      CALL lsysbl_duplicateMatrix (p_rmatrix,rmatrixTmpBlock, &
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      DO j=1,NDIM2D
        DO i=1,NDIM2D
          CALL lsyssc_releaseMatrix (rmatrixTmpBlock%RmatrixBlock(i,j))
        END DO
      END DO

      ! Build the rest of the defect. Here, no time step is involved!
      CALL lsysbl_blockMatVec (rmatrixTmpBlock, rx, rd, &
          -1.0_DP, 1.0_DP)
          
      CALL lsysbl_releaseMatrix (rmatrixTmpBlock)
    END IF
      
    ! Apply the filter chain to the defect vector -- if this is desired.
    IF (bboundaryConditions) &
      CALL filter_applyFilterChainVec (rd, p_RfilterChain)
    
    ! Filter the resulting defect vector through the slip-boundary-
    ! condition vector filter for implementing nonlinear slip boundary
    ! conditions into a defect vector. This changes the vector only IF we
    ! have slip boundary conditions!
    IF (bboundaryConditionsNonlin) &
      CALL vecfil_discreteNLSlipBCdef (rd)

  END SUBROUTINE

  ! ***************************************************************************
  ! Callback routines for the nonlinear solver
  ! ***************************************************************************

  !<subroutine>
  
    SUBROUTINE c2d2_getDefect (ite,rx,rb,rd,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd. The routine accepts a pointer to a collection structure 
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
    INTEGER, INTENT(IN)                           :: ite

    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN),TARGET        :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN)               :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled with data by the callback routine.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
  !</inputoutput>
  
  !</subroutine>
  
      ! The nonlinear iteration structure
      TYPE(t_ccNonlinearIteration) :: rnonlinearIteration
      
      ! DEBUG!!!
      !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

      ! Rebuild the nonlinear-iteration structure from the collection
      CALL c2d2_restoreNonlinearLoop (rnonlinearIteration,p_rcollection)
    
      ! Build the nonlinear defect
      CALL lsysbl_copyVector (rb,rd)
      
      ! DEBUG!!!
      !CALL lsysbl_getbase_double (rx,p_Ddata)
      !CALL lsysbl_getbase_double (rd,p_Ddata2)
      
      CALL c2d2_assembleNonlinearDefect (rnonlinearIteration,rx,rd,&
          .TRUE.,.TRUE.,p_rcollection)      
          
      ! Release the nonlinear-iteration structure, that's it.
      CALL c2d2_doneNonlinearLoop (rnonlinearIteration)
      
    END SUBROUTINE
    
  ! ***************************************************************************

  !<subroutine>

    SUBROUTINE c2d2_getOptimalDamping (rd,rx,rb,rtemp1,rtemp2,domega,&
        rnonlinearIteration,p_rcollection)
  
  !<description>
    ! This subroutine is called inside of the nonlinear loop, to be precise,
    ! inside of c2d2_precondDefect. It calculates an optiman damping parameter
    ! for the nonlinear defect correction.
    !
    ! The nonlinear loop reads:
    !
    !     $$ u_(n+1) = u_n + OMEGA * C^{-1}d_n $$
    !
    ! with $d_n$ the nonlinear defect and $C^{-1}$ a preconditioner (usually
    ! the linearised system).
    ! Based on the current solution $u_n$, the defect vector $d_n$, the RHS 
    ! vector $f_n$ and the previous parameter OMEGA, a new 
    ! OMEGA=domega value is calculated.
    !
    ! The nonlinear system matrix on the finest level in the collection is 
    ! overwritten by $A(u_n+domega_{old}*C^{-1}d_n)$.
  !</description>

  !<input>
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN)               :: rx

    ! Current RHS vector of the nonlinear equation
    TYPE(t_vectorBlock), INTENT(IN)               :: rb

    ! Defect vector b-A(x)x. 
    TYPE(t_vectorBlock), INTENT(IN)               :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
  !</input>

  !<inputoutput>
    ! Nonlinear iteration structure that describes the equation and 
    ! discretisation on all levels.
    TYPE(t_ccnonlinearIteration), INTENT(INOUT)   :: rnonlinearIteration

    ! A temporary vector in the structure of rx
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rtemp1

    ! A 2nd temporary vector in the structure of rx
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rtemp2

    ! Damping parameter. On entry: an initial value given e.g. by the
    ! previous step.
    ! On return: The new damping parameter.
    REAL(DP), INTENT(INOUT)                       :: domega
  !</inputoutput>
  
  !</subroutine>

    ! local variables
    INTEGER :: ilvmax
    REAL(DP) :: dskv1,dskv2
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes,p_rmatrixMass

    ! A filter chain for the linear solver
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

!    ! DEBUG!!!:
!    real(dp), dimension(:), pointer :: p_vec,p_def,p_temp1,p_temp2,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    call lsysbl_getbase_double (rtemp1,p_temp1)
!    call lsysbl_getbase_double (rtemp2,p_temp2)
!    ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Is there anything to do?
      IF (rnonlinearIteration%domegaMin .GE. rnonlinearIteration%domegaMax) THEN
        ! No - cancel.
        domega = rnonlinearIteration%domegaMin
        RETURN
      END IF
      
      ! Get minimum/maximum level from the collection
      ilvmax = rnonlinearIteration%NLMAX
      
      p_RfilterChain => rnonlinearIteration%rpreconditioner%p_RfilterChain
      
      ! Get the system and the Stokes matrix on the maximum level
      p_rmatrix => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrix
      p_rmatrixStokes => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixStokes
      p_rmatrixMass => rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass

      ! We now want to calculate a new OMEGA parameter
      ! with OMGMIN < OMEGA < OMGMAX.
      !
      ! The defect correction for a problem like T(u)u=f has the form
      !
      !       u_(n+1)  =  u_n  +  OMEGA * C * ( f - T(u_n)u_n )
      !                =  u_n  +  OMEGA * d_n
      !
      ! with an appropriate preconditioner C, which we don't care here.
      ! In our case, this iteration system can be written as:
      !
      ! (u1)     (u1)                     ( (f1)   [ A         B1] (u1) )
      ! (u2)  := (u2)  + OMEGA * C^{-1} * ( (f2) - [      A    B2] (u2) )
      ! (p )     (p )                     ( (fp)   [ B1^T B2^T 0 ] (p ) )
      !
      !                                   |------------------------------|
      !                                              = d_n
      !                            |-------------------------------------|
      !                                        = Y = (y1,y2,yp) = rd
      !
      ! with KST1=KST1(u1,u2,p) and Y=rd being the solution from
      ! the Oseen equation with 
      !
      !                  [ A         B1 ]
      !    C = T(u_n) =  [      A    B2 ]
      !                  [ B1^T B2^T 0  ]
      !
      ! The parameter OMEGA is calculated as the result of the 1D
      ! minimization problem:
      !
      !   OMEGA = min_omega || T(u^l+omega*Y)*(u^l+omega*Y) - f ||_E
      !
      !           < T(u^l+omegaold*Y)Y , f - T(u^l+omegaold*Y)u^l >
      !        ~= -------------------------------------------------
      !              < T(u^l+omegaold*Y)Y , T(u^l+omegaold*Y)Y >
      !
      ! when choosing omegaold=previous omega, which is a good choice
      ! as one can see by linearization (see p. 170, Turek's book).
      !
      ! Here, ||.||_E denotes the the Euclidian norm to the Euclidian 
      ! scalar product <.,.>.
      
      ! ==================================================================
      ! First term of scalar product in the nominator
      !
      ! Calculate the new nonlinear block A at the
      ! point rtemp1 = u_n + omegaold*Y
      ! ==================================================================
      !
      ! At first, calculate the point rtemp1 = u_n+omegaold*Y where
      ! to evaluate the matrix.

      CALL lsysbl_copyVector(rd,rtemp1)
      CALL lsysbl_vectorLinearComb (rx,rtemp1,1.0_DP,domega)

      ! Should we discretise the Navier-Stokes nonlinearity?
      ! Would mean to rebuild the system matrix. Otherwise we can reuse
      ! our previous diffusion matrix without rebuilding it.
      
      IF (rnonlinearIteration%dgamma .NE. 0.0_DP) THEN
      
        ! Re-assemble the nonlinear system matrix on the maximum level
        ! at the point rtemp1.
        ! Don't incorporate nonlinear (slip) boundary conditions.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIteration,p_rcollection,&
          .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,rtemp1)      
      
      END IF
        
      ! ==================================================================
      ! Second term of the scalar product in the nominator
      ! Calculate the defect rtemp2 = F-T*u_n.
      ! ==================================================================

      CALL lsysbl_copyVector (rb,rtemp2)
      CALL lsysbl_blockMatVec (p_rmatrix, rx, rtemp2, -1.0_DP, 1.0_DP)
      
      ! This is a defect vector - filter it! This e.g. implements boundary
      ! conditions.
      CALL filter_applyFilterChainVec (rtemp2, p_RfilterChain)
      
      ! Filter the resulting defect vector through the slip-boundary-
      ! condition vector filter for implementing nonlinear slip boundary
      ! conditions into a defect vector.
      !
      ! Note: implementation commented out!
      ! Seems to work better without!
      ! CALL vecfil_discreteNLSlipBCdef (rtemp2)
      
      ! ==================================================================
      ! For all terms in the fraction:
      ! Calculate the value  rtemp1 = T*Y
      ! ==================================================================

      CALL lsysbl_blockMatVec (p_rmatrix, rd, rtemp1, 1.0_DP, 0.0_DP)
      
      ! This is a defect vector against 0 - filter it! This e.g. 
      ! implements boundary conditions.
      CALL filter_applyFilterChainVec (rtemp1, p_RfilterChain)
      
      ! Filter the resulting defect vector through the slip-boundary-
      ! condition vector filter for implementing nonlinear slip boundary
      ! conditions into a defect vector.
      !
      ! Note: implementation commented out!
      ! Seems to work better without!
      ! CALL vecfil_discreteNLSlipBCdef (rtemp1)

      ! ==================================================================
      ! Calculation of the fraction terms.
      ! Calculate nominator:    dskv1:= (T*Y,D)   = (rtemp1,rtemp2)
      ! Calculate denominator:  dskv2:= (T*Y,T*Y) = (rtemp1,rtemp1)
      ! ==================================================================
      
      dskv1 = lsysbl_scalarProduct (rtemp1, rtemp2)
      dskv2 = lsysbl_scalarProduct (rtemp1, rtemp1)
      
      IF (dskv2 .LT. 1.0E-40_DP) THEN
        PRINT *,'Error in c2d2_getOptimalDamping. dskv2 nearly zero.'
        PRINT *,'Optimal damping parameter singular.'
        PRINT *,'Is the triangulation ok??? .tri-file destroyed?'
        STOP
      END IF
      
      ! Ok, we have the nominator and the denominator. Divide them
      ! by each other to calculate the new OMEGA.
      
      domega = dskv1 / dskv2
      
      ! And make sure it's in the allowed range:
      
      domega = max(rnonlinearIteration%domegamin, &
                   min(rnonlinearIteration%domegamax,domega))
      
      ! That's it, we have our new Omega.
  
    END SUBROUTINE

  ! ***************************************************************************

  !<subroutine>

    SUBROUTINE c2d2_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd. The routine accepts a pointer to a collection structure 
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<inputoutput>
    ! Number of current iteration. 
    INTEGER, INTENT(IN)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    REAL(DP), INTENT(INOUT)                       :: domega

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    LOGICAL, INTENT(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Current right hand side of the nonlinear system
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_vectorScalar), POINTER :: p_rvectorTemp,p_rvectorTemp2
    TYPE(t_vectorBlock) :: rtemp1,rtemp2
    INTEGER :: ierror
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    INTEGER, DIMENSION(3) :: Cnorms
    
    INTEGER :: i
    REAL(DP) :: dresInit,dres
    LOGICAL :: bassembleNewton
    TYPE(t_matrixBlock), DIMENSION(:), ALLOCATABLE :: Rmatrices
    TYPE(t_ccDynamicNewtonControl), POINTER :: p_rnewton

    ! The nonlinear iteration structure
    TYPE(t_ccNonlinearIteration), TARGET :: rnonlinearIteration

    ! DEBUG!!!
    !real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !call lsysbl_getbase_double (rd,p_def)
    !call lsysbl_getbase_double (rx,p_vec)
!    NLMAX = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Reconstruct the nonlinear iteration structure from the collection.
      ! We need it for some parameters.
      CALL c2d2_restoreNonlinearLoop (rnonlinearIteration,p_rcollection)

      SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
      CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
        ! Preconditioning with a linear solver.
        !
        ! At first, assemble the preconditioner matrices on all levels
        ! and incorporate all boundary conditions.
        !
        ! Should we assemble the Newton part?
        
        bassembleNewton = .FALSE.
        
        IF (rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
            CCPREC_NEWTON) THEN
            
          ! Use Newton in any case.
          bassembleNewton = .TRUE.
          
        ELSE IF (rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
            CCPREC_NEWTONDYNAMIC) THEN
            
          ! Adaptive Newton. Check the iteration and the residual whether to use
          ! Newton or not. But at first, get a shortcut to the parameter structure
          ! of the adaptive Newton...
            
          p_rnewton => rnonlinearIteration%rpreconditioner%radaptiveNewton
          
          IF (ite .GT. p_rnewton%nmaxFixPointIterations) THEN
            ! Force Newton to be used.
            bassembleNewton = .TRUE.
          ELSE
            IF (ite .GT. p_rnewton%nminFixPointIterations) THEN
              ! In this case, the residuum of the last iterate decides on 
              ! whether to use Newton or not.
              dresInit = SQRT(rnonlinearIteration%DresidualInit(1)**2 + &
                            rnonlinearIteration%DresidualInit(2)**2)
              dres    = SQRT(rnonlinearIteration%DresidualOld(1)**2 + &
                            rnonlinearIteration%DresidualOld(2)**2)
              IF ((dres .LT. p_rnewton%depsAbsNewton) .AND. &
                  (dres .LT. p_rnewton%depsRelNewton * dresInit)) THEN
                bassembleNewton = .TRUE.
              END IF
            END IF
            
            ! Otherwise: Use fixpoint iteration...
            
          END IF
          
        END IF
        
        ! Ok, now assemble the preconditioner matrices in  rnonlinearIteration.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIteration,p_rcollection,&
            .TRUE.,.TRUE.,.TRUE.,.TRUE.,bassembleNewton,rx)
          
        ! Our 'parent' (the caller of the nonlinear solver) has prepared
        ! a preconditioner node for us (a linear solver with symbolically
        ! factorised matrices). Get this from the collection.
      
        p_rsolverNode => rnonlinearIteration%rpreconditioner%p_rsolverNode

        ! Re-attach the system matrices to the solver.
        ! Note that no pointers and no handles are changed, so we can savely do
        ! that without calling linsol_doneStructure/linsol_doneStructure.
        ! This simply informs the solver about possible new scaling factors
        ! in the matrices in case they have changed...
        ALLOCATE(Rmatrices(rnonlinearIteration%NLMIN:rnonlinearIteration%NLMAX))
        DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          CALL lsysbl_duplicateMatrix ( &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        END DO
        
        CALL linsol_setMatrices(&
            rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(:))
            
        ! The solver got the matrices; clean up Rmatrices, it was only of temporary
        ! nature...
        DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          CALL lsysbl_releaseMatrix (Rmatrices(i))
        END DO
        DEALLOCATE(Rmatrices)

        ! Initialise data of the solver. This in fact performs a numeric
        ! factorisation of the matrices in UMFPACK-like solvers.
        CALL linsol_initData (p_rsolverNode, ierror)
        IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
        
        ! Finally solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        CALL linsol_precondDefect (p_rsolverNode,rd)
        
        ! Remember convergence rate for output
        rnonlinearIteration%drhoLinearSolver = p_rsolverNode%dconvergenceRate

        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
        ! for the next iteration.
        CALL linsol_doneData (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .EQ. 0
        
      END SELECT
      
      ! Finally calculate a new damping parameter domega.
      !
      ! For this purpose, we need two temporary vectors.
      ! On one hand, we have p_rvectorTemp.
      ! Get the second temporary vector from the collection as it was
      ! prepared by our 'parent' that invoked the nonlinear solver.
      p_rvectorTemp => rnonlinearIteration%rpreconditioner%p_rtempVectorSc
      p_rvectorTemp2 => rnonlinearIteration%rpreconditioner%p_rtempVectorSc2
      
      ! Both temp vectors are scalar, but we need block-vectors in the
      ! structure of rx/rb. Derive block vectors in that structure that
      ! share their memory with the scalar temp vectors. Note that the
      ! temp vectors are created large enough by our parent!
      CALL lsysbl_createVecFromScalar (p_rvectorTemp,rtemp1)
      CALL lsysbl_enforceStructure (rb,rtemp1)

      CALL lsysbl_createVecFromScalar (p_rvectorTemp2,rtemp2)
      CALL lsysbl_enforceStructure (rb,rtemp2)

      ! Calculate the omega
      CALL c2d2_getOptimalDamping (rd,rx,rb,rtemp1,rtemp2,&
                                  domega,rnonlinearIteration,p_rcollection)

      ! Remember damping parameter for output
      rnonlinearIteration%domegaNL = domega

      ! Release the temp block vectors. This only cleans up the structure.
      ! The data is not released from heap as it belongs to the
      ! scalar temp vectors.
      CALL lsysbl_releaseVector (rtemp2)
      CALL lsysbl_releaseVector (rtemp1)
      
      ! Calculate the max-norm of the correction vector.
      ! This is used for the stopping criterium in c2d2_resNormCheck!
      Cnorms(:) = LINALG_NORMMAX
      rnonlinearIteration%DresidualCorr(:) = lsysbl_vectorNormBlock(rd,Cnorms)
      
      ! Save the nonlinear-iteration structure.
      CALL c2d2_saveNonlinearLoop (rnonlinearIteration,p_rcollection)
      
      IF ((.NOT. bsuccess) .AND. (domega .GE. 0.001_DP)) THEN
        ! The preconditioner did actually not work, but the solution is not
        ! 'too bad'. So we accept it all the same.
        bsuccess = .TRUE.
      END IF

      ! Release the nonlinear-iteration structure, that's it.
      CALL c2d2_doneNonlinearLoop (rnonlinearIteration)
      
    END SUBROUTINE

  ! ***************************************************************************

    SUBROUTINE c2d2_resNormCheck (ite,rx,rb,rd,bconvergence,bdivergence,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! Residual norm calculation & printing routine.
    ! This routine is called each time the norm of the residuum was calculated.
    ! It has to check the current residuum for convergence and/or divergence
    ! and can print the residuum to screen.
  !</description>

  !<inputoutput>
    ! Number of current iteration. Is set to 0 when the callback routine
    ! is called the first time. In this situation, rd describes the initial
    ! defect vector.
    INTEGER, INTENT(IN)                           :: ite
  
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb

    ! Defect vector b-A(x)x.
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
  !</inputoutput>
  
  !<output>
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is within a desired tolerance, so that the solver should treat
    ! the iteration as 'converged'.
    LOGICAL, INTENT(OUT)                        :: bconvergence
  
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is out of a desired tolerance, so that the solver should treat
    ! the iteration as 'diverged'.
    LOGICAL, INTENT(OUT)                        :: bdivergence
  !</output>

      ! local variables
      REAL(DP), DIMENSION(3) :: Dresiduals
      REAL(DP) :: dresOld,drhoNL,ddelP,ddelU,dtmp,dresU,dresDIV,dres,dresINIT
      REAL(DP) :: depsD,depsDiv,depsUR,depsPR,depsRES
      INTEGER, DIMENSION(3) :: Cnorms

      ! The nonlinear iteration structure
      TYPE(t_ccNonlinearIteration), TARGET :: rnonlinearIteration

      ! Reconstruct the nonlinear iteration structure from the collection.
      ! We need it for some parameters.
      CALL c2d2_restoreNonlinearLoop (rnonlinearIteration,p_rcollection)

      ! Calculate norms of the solution/defect vector
      CALL c2d2_getDefectNorm (rx,rb,rd,Dresiduals)
      Dresiduals(3) = SQRT(Dresiduals(1)**2 + Dresiduals(2)**2)

      ! Replace the 'old' residual by the current one
      rnonlinearIteration%DresidualOld(1:2) = Dresiduals(1:2)

      ! In the first iteration (initial defect), print the norm of the defect
      ! and save the norm of the initial residuum to the structure
      IF (ite .EQ. 0) THEN
      
        CALL output_separator (OU_SEP_MINUS)     
        CALL output_line (' IT  RELU     RELP     DEF-U    DEF-DIV'// &
                          '  DEF-TOT  RHONL    OMEGNL   RHOMG')
        CALL output_separator (OU_SEP_MINUS)     
        CALL output_line ('  0                   '// &
            TRIM(sys_sdEP(Dresiduals(1),9,2))//&
            TRIM(sys_sdEP(Dresiduals(2),9,2))//&
            TRIM(sys_sdEP(Dresiduals(3),9,2)))
        CALL output_separator (OU_SEP_MINUS)     

        rnonlinearIteration%DresidualInit (1:2) = Dresiduals(1:2)
        rnonlinearIteration%DresidualOld (1:2) = Dresiduals(1:2)

        bconvergence = .FALSE.
        bdivergence = .FALSE.
      
      ELSE
        ! In the other iterations, calculate the relative change and
        ! test the convergence criteria.
      
        ! Old defect:
        dresOld = SQRT(rnonlinearIteration%DresidualOld(1)**2 + &
                       rnonlinearIteration%DresidualOld(2)**2)
        
        ! Nonlinear convergence rate
        drhoNL = (Dresiduals(3)/dresOld) ** (1.0_DP/REAL(ite,DP))
        
        ! Calculate norms of the solution/defect vector, calculated above
        dresU   = Dresiduals(1)
        dresDIV = Dresiduals(2)
        dres    = SQRT(dresU**2 + dresDIV**2)
        
        ! Calculate relative maximum changes 
        ! This simply calculates some postprocessing values of the relative
        ! change in the solution.
        !
        ! Maximum norm of solution vector:
        !
        !              || YP ||_max       || Pnew - Pold ||_max
        !   DELP := ------------------- = ---------------------
        !              || P ||_max           || Pnew ||_max
        !
        !
        ! Relative change of solution vector:
        !
        !            || (Y1,Y2) ||_max    || Unew - Uold ||_max
        !   DELU := ------------------- = --------------------- 
        !           || (KU1,KU2) ||_max       || Unew ||_max
        !
        ! The norms || YP ||_max, || Yi ||_max are saved in the nonlinear
        ! iteration structure from the last preconditioning!
        ! The other MAX-norms have to be calculated from U...
        
        Cnorms(:) = LINALG_NORMMAX
        Dresiduals(:) = lsysbl_vectorNormBlock (rx,Cnorms)

        dtmp = MAX(Dresiduals(1),Dresiduals(2))
        IF (dtmp .LT. 1.0E-8_DP) dtmp = 1.0_DP
        ddelU = MAX(rnonlinearIteration%DresidualCorr(1),&
                    rnonlinearIteration%DresidualCorr(2))/dtmp
        
        dtmp = Dresiduals(3)
        IF (dtmp .LT. 1.0E-8_DP) dtmp = 1.0_DP
        ddelP = rnonlinearIteration%DresidualCorr(3)/dtmp
        
        ! Check if the nonlinear iteration can prematurely terminate.
        !        
        ! Get the stopping criteria from the parameters.
        ! Use the DepsNL data according to the initialisation above.
        dresINIT = SQRT(rnonlinearIteration%DresidualInit(1)**2 + &
                        rnonlinearIteration%DresidualInit(2)**2)
                        
        depsD   = rnonlinearIteration%DepsNL(1)
        depsDiv = rnonlinearIteration%DepsNL(2)
        depsUR  = rnonlinearIteration%DepsNL(3)
        depsPR  = rnonlinearIteration%DepsNL(4)
        depsRES = rnonlinearIteration%DepsNL(5)*dresINIT
        
        ! All residual information calculated.
        ! Check for divergence; use a 'NOT' for better NaN handling.
        bdivergence = .NOT. (dres/dresOld .LT. 1E2)
        
        ! Check for convergence
        IF((ddelU .LE. depsUR).AND.(ddelP .LE. depsPR)   .AND. &
           (dresU .LE. depsD) .AND.(depsDiv .LE. depsDiv).AND. &
           (dres .LE. depsRES)) THEN
          bconvergence = .TRUE.
        ELSE
          bconvergence = .FALSE.
        END IF
        
        ! Print residual information
        CALL output_line ( &
            TRIM(sys_si(ite,3))//' '//&
            TRIM(sys_sdEP(ddelU,9,2))// &
            TRIM(sys_sdEP(ddelP,9,2))// &
            TRIM(sys_sdEP(dresU,9,2))// &
            TRIM(sys_sdEP(dresDIV,9,2))// &
            TRIM(sys_sdEP(dres,9,2))// &
            TRIM(sys_sdEP(drhoNL,9,2))// &
            TRIM(sys_sdEP(rnonlinearIteration%domegaNL,9,2))// &
            TRIM(sys_sdEP(rnonlinearIteration%drhoLinearSolver,9,2)) &
            )
        
      END IF
      
      ! Save and release the nonlinear-iteration structure, that's it.
      CALL c2d2_saveNonlinearLoop (rnonlinearIteration,p_rcollection)
      CALL c2d2_doneNonlinearLoop (rnonlinearIteration)

    END SUBROUTINE

! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_getDefectNorm (rvector,rrhs,rdefect,Dresiduals)
  
!<description>
  ! Calculates a couple of norms from a given solution, defect and RHS vector
  ! of the nonlinear system. This can be used to check convergence criteria
  ! etc.
!</description>

!<input>
  ! The solution vector which is modified later during the nonlinear iteration.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
  
  ! A defect vector calculated with rvector and rrhs
  TYPE(t_vectorBlock), INTENT(IN) :: rdefect
!</input>

!<output>
  ! An array receiving different defect norms calculated by the above vectors.
  ! Dresiduals(1) = RESU   = ||defect_u|| / ||rhs|| = velocity residual
  ! Dresiduals(2) = RESDIV = ||p|| / ||u||          = divergence residual
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dresiduals
!</output>

!</subroutine>

    ! local variables
    REAL(DP) :: dresF,DresTmp(2),dnormU
    INTEGER, DIMENSION(2) :: Cnorms

    !-----------------------------------------------------------------------
    !     Compute the relative l2-norms  RESU,RESDIV
    !-----------------------------------------------------------------------

    Cnorms(:) = LINALG_NORMEUCLID

    ! RESF := max ( ||F1||_E , ||F2||_E )

    DresTmp = lsysbl_vectorNormBlock (rrhs,Cnorms)
    dresF = MAX(DresTmp(1),DresTmp(2))
    IF (dresF .LT. 1.0E-8_DP) dresF = 1.0_DP

    !               || (D1,D2) ||_E
    ! RESU = -----------------------------
    !        max ( ||F1||_E , ||F2||_E )

    DresTmp = lsysbl_vectorNormBlock (rdefect,Cnorms)
    Dresiduals(1) = SQRT(DresTmp(1)**2+DresTmp(2)**2)/dresF

    ! DNORMU = || (U1,U2) ||_l2 

    DresTmp = lsysbl_vectorNormBlock (rvector,Cnorms)
    dnormU = SQRT(DresTmp(1)**2+DresTmp(2)**2)
    IF (dnormU .LT. 1.0E-8_DP) dnormU = 1.0_DP

    !             || DP ||_E
    ! RESDIV = ----------------
    !          || (U1,U2) ||_E

    Dresiduals(2) = &
        lsyssc_vectorNorm (rdefect%RvectorBlock(3),LINALG_NORMEUCLID) / dnormU
        
  END SUBROUTINE

! ***************************************************************************
! Routines to solve the core equation
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_getNonlinearSolver (rnlSolver, rparamList, sname)
  
!<description>
  ! Creates a nonlinear solver node rnlSolver and initialises it with parameters
  ! from the INI/DAT files given in the parameter list rparamList.
  ! sname is the name of a section in rparamList that configures the
  ! parameter of the nonlinear solver.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  TYPE(t_parlist), INTENT(IN) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  CHARACTER(LEN=*), INTENT(IN) :: sname
!</input>

!<output>
  ! A t_nlsolNode structure that contains the configuration of the nonlinear
  ! solver. The parameters are initialised according to the information
  ! in the section sname of the parameter list rparamList
  TYPE(t_nlsolNode) :: rnlSolver
!</output>

    ! local variables
    TYPE(t_parlstSection), POINTER :: p_rsection

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    CALL parlst_querysection(rparamList, sname, p_rsection) 

    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      PRINT *,'Cannot create nonlinear solver; no section '''&
              //TRIM(sname)//'''!'
      STOP
    END IF
    
    ! Parse the given parameters now to initialise the solver node.
    ! This is now CCxD-specific!
    
    CALL parlst_getvalue_int (p_rsection, 'nminIterations', &
                              rnlSolver%nminIterations, rnlSolver%nminIterations)

    CALL parlst_getvalue_int (p_rsection, 'nmaxIterations', &
                              rnlSolver%nmaxIterations, rnlSolver%nmaxIterations)

    CALL parlst_getvalue_int (p_rsection, 'ioutputLevel', &
                              rnlSolver%ioutputLevel, rnlSolver%ioutputLevel)

    CALL parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnlSolver%DepsRel(1), rnlSolver%DepsRel(1))
    rnlSolver%DepsRel(2) = rnlSolver%DepsRel(1)

    CALL parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnlSolver%DepsRel(3), rnlSolver%DepsRel(3))

    CALL parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnlSolver%DepsAbs(1), rnlSolver%DepsAbs(1))
    rnlSolver%DepsAbs(2) = rnlSolver%DepsAbs(1)

    CALL parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnlSolver%DepsAbs(3), rnlSolver%DepsAbs(3))

    ! Initial damping parameter.
    CALL parlst_getvalue_double (p_rsection, 'domegaIni', &
                                 rnlSolver%domega, rnlSolver%domega)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_setupCoreEquation (rnonlinearIteration,&
      dalpha,dtheta,dgamma)
  
!<description>
  ! Initialises the coefficients in front of the terms of the core equation.
  ! Saves the pointers of all matrices to rnonlinearIteration
  ! that are needed during the assembly process of the linearised 
  ! nonlinear matrix.
!</description>

!<input>
  ! Weight in front of the mass matrix in the core equation
  REAL(DP), INTENT(IN) :: dalpha
  
  ! Weight in front of the Stokes matrix in the core equation
  REAL(DP), INTENT(IN) :: dtheta
  
  ! Weight in front of the nonlinearity in the core equation
  REAL(DP), INTENT(IN) :: dgamma

!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! Initialise the core equation parameters on all levels
    rnonlinearIteration%dalpha = dalpha
    rnonlinearIteration%dtheta = dtheta
    rnonlinearIteration%dgamma = dgamma
    
  END SUBROUTINE

! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_solveCoreEquation (rnlSolver,rnonlinearIteration,&
      rvector,rrhs,rcollection,rtempBlock)
  
!<description>
  ! This routine invokes the nonlinear solver to solve the core equation
  ! as configured in the core equation structure.
!</description>

!<input>
  ! The right-hand-side vector to use in the equation.
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
!</input>
  
!<inputoutput>
  ! A nonlinear solver configuration.
  ! Can be initialised e.g. by using c2d2_getNonlinearSolver.
  TYPE(t_nlsolNode), INTENT(INOUT) :: rnlSolver
  
  ! A nonlinear-iteration structure that configures the core equation to solve.
  ! Can be initialised e.g. by c2d2_createNonlinearLoop + c2d2_setupCoreEquation.
  TYPE(t_ccNonlinearIteration) :: rnonlinearIteration

  ! Initial solution vector. Is replaced by the new solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! A collection structure that is used to save parameters during the nonlinear
  ! iteration.
  TYPE(t_collection), INTENT(INOUT) :: rcollection

  ! A temporary block vector for the nonlinear iteration.
  ! OPTIONAL: If not specified, a temporary vector is automatically allocated.
  TYPE(t_vectorBlock), INTENT(INOUT), TARGET, OPTIONAL :: rtempBlock
  
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_vectorBlock), POINTER :: p_rtempBlock

    IF (.NOT. PRESENT(rtempBlock)) THEN
      ! Create a temporary vector we need for the nonlinear iteration.
      ALLOCATE (p_rtempBlock)
      CALL lsysbl_createVecBlockIndirect (rrhs, p_rtempBlock, .FALSE.)
    ELSE 
      p_rtempBlock => rtempBlock
    END IF

    ! Save the nonlinear-iteration structure to the collection.
    ! It's rebuild from the collection in the callback routines.
    CALL c2d2_saveNonlinearLoop (rnonlinearIteration,rcollection)

    ! Call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    CALL nlsol_performSolve(rnlSolver,rvector,rrhs,p_rtempBlock,&
                            c2d2_getDefect,c2d2_precondDefect,c2d2_resNormCheck,&
                            rcollection=rcollection)
        
    ! Remove parameters of the nonlinear loop from the collection.
    CALL c2d2_removeNonlinearLoop (rnonlinearIteration,rcollection)
             
    IF (.NOT. PRESENT(rtempBlock)) THEN
      ! Release the temporary vector
      CALL lsysbl_releaseVector (p_rtempBlock)
      DEALLOCATE (p_rtempBlock)
    END IF
        
  END SUBROUTINE

END MODULE
