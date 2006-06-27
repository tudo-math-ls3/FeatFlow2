!##############################################################################
!# ****************************************************************************
!# <name> linearsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a collection of linear solvers of various types.
!# The solvers are highly integrated and allow to be nested into each
!# other in a very general order. On the other hand, the solvers are
!# kept "black box like", so they can be configured from outside for
!# a large varienty of problems.
!#
!# The solvers basically deal with block systems of the following form:
!#
!# 1.) General linear systems of the form
!#         $$Ax = b$$
!#     with A being an arbitrary block matrix
!# 2.) Saddle-point block matrices of the form
!#           [A  B1] (x) = (f1)
!#           [B2  0] (y)   (f2)
!#     with a regular block matrix A and general block matrices B1 and B2
!# 3.) Saddle-point-like block matrices of the form
!#           [A  B1] (x) = (f1)
!#           [B2  C] (y)   (f2)
!#     with a regular block matrix A, general block matrices B1 and B2
!#     and a 'stabilisation' matrix $C \sim 0$.
!#
!# Cases 2.) + 3.) are special cases of 1.) in that sense, that specialised
!# solvers are available for these kind of systems which allow faster solving.
!#
!# To solve a problem, one has basically to call the following routines
!# 1.) linsol_initXXXX                - Initialises a solver, returns a solver
!#                                      structure identifying the solver
!# 2.) linsol_setMatrices             - Attach the system matrix/matrices to the
!#                                      solver
!# 3.) linsol_initStructure           - Allow the solvers to perform problem-
!#                                      structure specific initialisation
!#                                      (e.g. symbolical factorisation).
!#                                      During this phase, each solver allocates
!#                                      temporary memory it needs for the solution
!#                                      process
!# 3.) linsol_initData                - Allow the solvers to perform problem-
!#                                      data specific initialisation
!#                                      (e.g. numerical factorisation)
!# 4.) linsol_precondDefect           - Precondition a defect vector with
!#                                      a linear solver
!# or  linsol_solveAdaptively         - Solve the problem with an initial 
!#                                      solution vector
!# 5.) linsol_doneData                - Release problem-data specific information
!# 6.) linsol_doneStructure           - Release problem-structure specific 
!#                                      information. Release temporary memory.
!# 7.) linsol_releaseSolver           - Clean up solver structures, remove the
!#                                      solver structure from the heap.
!#
!# For initialising a multigrid solver, this sequenc changes a bit,
!# since the data of multiple levels have to be associated to
!# the solver before invoking it:
!#
!# 1.) linsol_initMultigrid           - Initialises multigrid  solver, returns a 
!#                                      solver structure identifying the solver
!# 2.) a) On the coarse grid:
!#          linsol_initXXXX           - Initialise coarse grid solver structure
!#
!#        On finer grids up to the maximum level:
!#          linsol_initXXXX
!#        + linsol_convertToSmoother  - Initialise a solver structure that can be 
!#                                      attached as smoother to one level of multigrid.
!#                                      *Not on the coars grid!'
!#     b) linsol_addMultigridLevel    - Add a multigrid level. First the coarse grid
!#                                      must be added, then level 2, then level 3 etc.
!#
!#     Steps a), b) must be repeated for all levels.
!#
!#     Btw., levels can be removed from the solver with linsol_removeMultigridLevel.
!#     The level stucture can be obtained with linsol_getMultigridLevel.
!#
!# 3.) linsol_setMatrices             - Attach the system matrices of all levels to the
!#                                      solver
!# 4.) linsol_initStructure           - Allow the solvers to perform problem-
!#                                      structure specific initialisation
!#                                      (e.g. symbolical factorisation).
!#                                      During this phase, each solver allocates
!#                                      temporary memory it needs for the solution
!#                                      process
!# 5.) linsol_initData                - Allow the solvers to perform problem-
!#                                      data specific initialisation
!#                                      (e.g. numerical factorisation)
!# 6.) linsol_precondDefect           - Precondition a defect vector 
!# or  linsol_solveAdaptively         - Solve the problem with an initial 
!#                                      solution vector
!# 7.) linsol_doneData                - Release problem-data specific information
!# 8.) linsol_doneStructure           - Release problem-structure specific 
!#                                      information. Release temporary memory.
!# 9.) linsol_releaseSolver           - Clean up solver structures, remove the
!#                                      solver structure from the heap. Remove all
!#                                      attached subsolvers (smoothers, coarse grid 
!#                                      solver) from the heap.
!#
!# Remark: THE SYSTEM MATRIX PRESCRIBES THE SPATIAL DISCRETISATION AND 
!#         THE BOUNDARY CONDITIONS!
!#
!# This allows the solver to create temporary vectors if necessary. In 
!# consequence, when an application calls a solver, the RHS/solution
!# vectors passed to the solver must match in their discretisation/
!# boundary conditions to the matrices attached previously!
!#
!# The following routines serve as auxiliary routines for the application to
!# maintain a Multigrid solver node:
!# 
!# 1.) linsol_removeMultigridLevel
!#     -> Deletes a level and attached solver structures of that level from 
!#        Multigrid
!#
!# 2.) linsol_cleanMultigridLevels
!#     -> Deletes all levels and attached solver structures
!#
!# 3.) linsol_getMultigridLevel
!#     -> Returns the level structure of a specific level
!#
!# 4.) linsol_getMultigridLevelCount
!#     -> Returns number of currently attached levels
!#
!#
!# Implementational details / structure of the solver library
!# ----------------------------------------------------------
!# When going through this library, a newcomer would think: 
!#
!#       'Where are the solvers? I can only find preconditioners!?!
!#        And what should that thing with the defect vectors mean?'
!# 
!# The reason is simple: 'Everything is a preconditioner for defect vectors!'
!#
!# A short explaination: Let's assume we want to solve a linear system:
!#
!#                  $$ Ax=b $$
!#
!# For this we have a linear solver $P^{-1} \approx A^{-1}$ (Gauss elimination
!# or whatever), so we get our solution vector x by saying:
!#
!#                 $$ x = P^{-1} b $$
!#
!# This is the usual way how to deal with a linear system in Numrics 1,
!# but it's a little bit hard to deal with in computer science. But we can
!# deal with the linear system more easily if we perform a simple
!# modification to this equation:
!# 
!#                          $$ Ax = b $$
!# $$ \Rightarrow        P^{-1}Ax = P^{-1}b             $$ 
!# $$ \Rightarrow               0 = P^{-1} (b-Ax)       $$
!#
!# $$ \Rightarrow         x_{n+1} = x_n  +  P^{-1} (b-Ax) $$
!#
!# So the linear solver $P^{-1}$ can equivalently be used in a defect
!# correction approach as a preconditioner for the defect $(b-Ax)$!
!# This can obviously be done for every linear solver.
!# So, this allows us to make a normalisation: 
!#
!# ! All linear solvers can be formulated as a preconditioner to be applied
!#   to defect vectors, we don't have to take care of solution vectors !
!#
!# (and with the problems related to that like implementing boundary
!# conditions...)
!#
!# So everything that can be found in this library is a preconditioner.
!# The xxx_solve routines are all formulated that way, that they accept
!# the defect vector $d := (b-Ax)$ and overwrite it by the preconditioned
!# defect: 
!#             $$ d := P^{-1}d $$
!#
!# which is actually the same as solving the linear system 
!'
!#             $$ Pd_{new} = d $$ 
!#
!# with the right hand side $d$ being a defect vector given from outside. 
!# Solving $Pd_{new} = d$ can then be done by an arbitrary linear solver, 
!# e.g. BiCGStab, UMFPACK (which takes $P=A$ and solves directly using Gauss),
!# or can even be the application pf a single Jacobi or ILU(s) 
!# preconditioner - everything is the same!
!#
!# Remark: As sometimes ( :-) ) the user wants to solve a system $Ax=b$ with
!#   a given solution vector and a given RHS which is not a defect vector,
!#   the routine 'linsol_solveAdaptively' can be called. This builds the
!#   defect, performs a preconditioning and corrects the solution vector -
!#   which is indeed solving the problem.
!#   When calling this routine with x:=0, the defect is exactly the RHS, so
!#   this routine falls back to the solution method known from Numerics 1:
!#
!#      $$  x_{new}  =  x + P^{-1}(b-Ax)  =  P^{-1}b  $$
!#
!#  FAQ - Frequently asked Questions
!# ----------------------------------
!# Ok, the solver library seems to be a little bit complicated. Therefore
!# here a small chapter about how to set up specific tasks for newbies:
!#
!# 1.) I'd like to set up a simple BiCGStab solver with ILU(0) preconditioner.
!#     How to do that?
!#
!#     Declare two solver structures, create the ILU(0) solver in the first,
!#     BiCGStab solver in the second and hang in the ILU(0) solver as
!#     preconditioner. Use NULL() for the filter chain as long as
!#     you don't need it:
!#
!#      TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
!#      TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
!#
!#      NULLIFY(p_RfilterChain)
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
!#
!#      -> p_rsolverNode is the solver identifying your BiCGStab(ILU(0)) 
!#         solver.
!#
!# 2.) In numerics 1 I learned the pure Richardson iteration. Where can I
!#     find it?
!#
!#     Pure Richardson iteration ($x_{n+1} = x_n + \omega (b-Ax)$) is called
!#     'defect correction' here. Set up a defect correction solver and
!#     modify the domega parameter in the solver node:
!#
!#      TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
!#      TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
!#
!#      NULLIFY(p_RfilterChain)
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
!#      p_rsolverNode%domega = 0.001   ! or whatever you like
!#
!# 3.) And where's the good old Jacobi iteration hidden? There's only that
!#     preconditioner!?!
!#
!#     Yes, i have only the Jacobi preconditioner. To get a Jacobi iteration,
!#     you must use it as a preconditioner in a defect correction loop:
!#
!#       $$ x_{n+1}  =  x_n  +  \omega D^{-1}  (b-Ax) $$
!#          ------------------  ^^^^^^^^^^^^^  ------  Defect correction loop
!#                              Jacobi preconditioner
!#     So,
!#      TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
!#      TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
!#
!#      NULLIFY(p_RfilterChain)
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initJacobi (p_rpreconditioner,0.7_DP)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
!#
!#    Alternatively, you can use the domega parameter from the defect
!#    correction to specify the $\omega$, while using no damping parameter
!#    in the actual Jacobi preconditioner:
!#
!#      NULLIFY(p_RfilterChain)
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initJacobi (p_rpreconditioner)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
!#      p_rsolverNode%domega = 0.7_DP
!#
!#  4.) I want to plug in my own solver XY. What do I have to do?
!#
!#      Well, the all solvers are build up in the same way with the same
!#      interfaces - that's the reason why every solver can be used as
!#      a preconditioner in another one. If you want to create your
!#      own solver and plug it into here, use the following guidline
!#      of what to do:
!#      a) Add new solver identifier LINSOL_ALG_xxxx
!#      b) Create solver-specific substructure (t_linsolSubnodeXXXX) if 
!#         necessary and add a pointer to the subnode in the main solver 
!#         structure.
!#      c) Copy/Paste solver-specific routines and rename them:
!#          - linsol_initXXXX
!#          - linsol_setMatrixXXXX
!#          - linsol_initStructureXXXX
!#          - linsol_initDataXXXX
!#          - linsol_doneDataXXXX
!#          - linsol_doneStructureXXXX
!#          - linsol_doneXXXX
!#         Don't forget to deallocate in linsol_doneXXXX everything that
!#         was allocated in linsol_initXXXX!!!
!#      d) Modify the routines so that they do their work.
!#      e) Alter the central solver routines to make your solver
!#         interface public:
!#          - linsol_initStructure
!#          - linsol_initData
!#          - linsol_doneData
!#          - linsol_doneStructure
!#          - linsol_releaseSolver
!#          - linsol_precondDefect
!#          - linsol_setMatrices
!#         Add calls to the new solver routines there.
!#      That's it, your solver should now be able to work ;-)
!#
!#  5.) And what if I need temporary vectors in my solver?
!#
!#      Allocate temporary vectors in linsol_initStructureXXXX,
!#      release them in linsol_doneStructureXXXX.
!# 
!#  6.) But if I make that, I often get messages like
!#      "Vector/Matrix not compatible"!?!
!#      How can I make sure, that the matrix/vector are compatible?
!#
!#      One remark to that: THE SYSTEM MATRIX PRESCRIBES THE SPATIAL
!#      DISCRETISATION AND THE BOUNDARY CONDITIONS!
!#
!#      Therefore, temporary vectors can be created with 
!#      lsysbl_createVecBlockIndMat using the system matrix as template.
!#      This also sets the boundary conditions etc. correctly, i.e.
!#      produce a compatible vector.
!#      Inside of the solver, you can also transfer properties of the 
!#      vector to your local temporary vectors using 
!#      lsysbl_assignDiscretIndirect.
!#
!#  7.) I want to use sorted vectors/matrices. Should I plug them into
!#      the solver sorted or not?
!#
!#      Plug them in sorted. Most of the solvers work only by matrix/vector
!#      multiplication and therefore directly need them sorted.
!#      Only some algorithms (like Multigrid) work internally on some places
!#      with unsorted vectors and therefore unsort them before using them.
!#      
!# </purpose>
!##############################################################################

MODULE linearsolver

  USE fsystem
  USE storage
  USE linearsystemblock
  USE multilevelprojection
  USE filtersupport
  USE coarsegridcorrection
  USE vanca
  USE globalsystem
  
  IMPLICIT NONE

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<constants>

!<constantblock description="Algorithm identifiers">

  ! Undefined algorithm
  INTEGER, PARAMETER :: LINSOL_ALG_UNDEFINED     = 0
  
  ! Preconditioned defect correction (Richardson iteration);
  ! $x_{n+1} = x_n + \omega P^{-1} (b-Ax)$
  INTEGER, PARAMETER :: LINSOL_ALG_DEFCORR       = 1
  
  ! Jacobi iteration $x_1 = x_0 + \omega D^{-1} (b-Ax_0)$
  INTEGER, PARAMETER :: LINSOL_ALG_JACOBI        = 2
  
  ! SOR/GS iteration $x_1 = x_0 + (L+\omega D)^{-1}(b-Ax_0)$
  INTEGER, PARAMETER :: LINSOL_ALG_SOR           = 4
  
  ! SSOR iteration
  INTEGER, PARAMETER :: LINSOL_ALG_SSOR          = 5
  
  ! CG iteration (preconditioned) 
  INTEGER, PARAMETER :: LINSOL_ALG_CG            = 6
  
  ! BiCGStab iteration (preconditioned) 
  INTEGER, PARAMETER :: LINSOL_ALG_BICGSTAB      = 7

  ! GMRES iteration (preconditioned) 
  INTEGER, PARAMETER :: LINSOL_ALG_GMRES         = 8
  
  ! Multigrid iteration
  INTEGER, PARAMETER :: LINSOL_ALG_MULTIGRID     = 9

  ! UMFPACK2
  INTEGER, PARAMETER :: LINSOL_ALG_UMFPACK2      = 10

  ! UMFPACK4
  INTEGER, PARAMETER :: LINSOL_ALG_UMFPACK4      = 11
  
  ! ILU(0) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_ILU01x1       = 50
  
  ! (M)ILU(s) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_MILUS1x1      = 51
  
  ! SPAI(0) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_SPAI01x1      = 52

  ! SPAI(k) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_SPAIK1x1      = 53
  
  ! VANCA iteration
  INTEGER, PARAMETER :: LINSOL_ALG_VANCA         = 54
  
!</constantblock>

!<constantblock description="Identifiers for stopping criteria">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. stop if both criteria hold
  INTEGER, PARAMETER :: LINSOL_STOP_STANDARD     = 0
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Bitfield identifiers for the ability of a linear solver">

  ! Solver can handle scalar systems
  INTEGER, PARAMETER :: LINSOL_ABIL_SCALAR       = 2**0

  ! Solver can handle block systems
  INTEGER, PARAMETER :: LINSOL_ABIL_BLOCK        = 2**1

  ! Solver can handle multiple levels
  INTEGER, PARAMETER :: LINSOL_ABIL_MULTILEVEL   = 2**2
  
  ! Solver allows checking the defect during the iteration.
  ! Solvers not capable of this perform only a fixed number of solution
  ! steps (e.g. UMFPACK performs always one step).
  INTEGER, PARAMETER :: LINSOL_ABIL_CHECKDEF     = 2**3
  
  ! Solver is a direct solver (e.g. UMFPACK, ILU).
  ! Otherwise the solver is of iterative nature and might perform
  ! multiple steps to solve the problem.
  INTEGER, PARAMETER :: LINSOL_ABIL_DIRECT       = 2**4
  
  ! Solver might use subsolvers (preconditioners, smoothers,...)
  INTEGER, PARAMETER :: LINSOL_ABIL_USESUBSOLVER = 2**5
  
  ! Solver supports filtering
  INTEGER, PARAMETER :: LINSOL_ABIL_USEFILTER    = 2**6

!</constantblock>

! *****************************************************************************

!<constantblock description="Error constants returned by initialisation routines">

  ! Initialisation routine went fine
  INTEGER, PARAMETER :: LINSOL_ERR_NOERROR       = 0

  ! Warning: Singular matrix
  INTEGER, PARAMETER :: LINSOL_ERR_SINGULAR      = -1

  ! Error during the initialisation: Not enough memory
  INTEGER, PARAMETER :: LINSOL_ERR_NOMEMORY      = 1

  ! General error during the initialisation
  INTEGER, PARAMETER :: LINSOL_ERR_INITERROR     = 2
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Variants of the VANCA solver">

  ! General VANCA solver
  INTEGER, PARAMETER :: LINSOL_VANCA_GENERAL     = 0   
  
  ! Simple Jacobi-like VANCA, 2D saddle-point problem, $\tilde Q_1/P_0$
  ! discretisation
  INTEGER, PARAMETER :: LINSOL_VANCA_2DSPQ1TQ0   = 1

!</constantblock>

!</constants>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<types>

!<typeblock>
  
  ! This is the central structure which repesents a solver.
  ! It collects all information for the solution process.
  !
  ! The caller who wants to solve a problem has to initialise
  ! such a structure with all information that is necessary
  ! to solve the linear system. For more complicated solvers
  ! which contain sub-solvers like preconditioners or smoothers
  ! (like in MG), such a structure must be created also for the
  ! sub-solver and attached to the main solver. This way, the
  ! structures form a ''solver tree'' with one solver-structure
  ! forming the root and the other structures, which are attached
  ! to this, forming the branches/nodes.
  !
  ! The structure contains general solver parameters which appear is
  ! many solvers. Whether or not a special parameter is respected
  ! by a solver depends on the solver (e.g. UMFPACK won't respect
  ! the ''maximum iterations'' parameter). Further solver-specific
  ! parameters which fit not into the general scheme are kept in
  ! solver specific structures. 
  !
  ! The structure contains INPUT, OUTPUT and STATUS parameters. INPUT
  ! parameters can be changed by the caller prior to solving.
  ! OUTPUT parameters indicate the final status of the solution process.
  ! STATUS parameters are only valid during the solution
  ! process. They should not be accessed from outside, as they are
  ! maintained internally by the solver.
  ! Parameters marked as READ ONLY are set during the initialisation phase
  ! of the solver and must not be changed by the caller.
  
  TYPE t_linsolNode
    
    ! OUTPUT: Result
    ! The result of the solution process.
    ! =0: success. =1: iteration broke down, diverging, =2: error in the parameters,
    ! <0: algorithm-specific error
    INTEGER                    :: iresult
    
    ! OUTPUT: Number of performed iterations, if the solver
    ! is of iterative nature.
    ! Is to 1 by the solver if not used (indicating at least 1 performed 
    ! iteration, which is always the case).
    INTEGER                    :: iiterations
    
    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Norm of initial residuum
    REAL(DP)                        :: dinitialDefect

    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Norm of final residuum
    REAL(DP)                        :: dfinalDefect

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK: 
    ! Convergence rate
    REAL(DP)                        :: dconvergenceRate

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK: 
    ! Asymptotic convergence rate
    REAL(DP)                        :: dasymptoticConvergenceRate

    ! OUTPUT PARAMETER:
    ! Total time for solver
    REAL(DP)                        :: dtimeTotal

    ! OUTPUT PARAMETER FOR SOLVERS THAT SUPPORT FILTERING:
    ! Total time for filtering
    REAL(DP)                        :: dtimeFiltering
    
    ! INPUT PARAMETER:
    ! General solver parameter; solver specific use.
    ! Standard value = 1.0 (corresponds to 'no damping' e.g. with the defect
    ! correction iteration)
    REAL(DP)                        :: domega  = 1.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    REAL(DP)                        :: depsRel = 1E-5_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    REAL(DP)                        :: depsAbs = 1E-5_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Relative divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY disables the relative divergence check.
    ! standard = 1E3
    REAL(DP)                        :: ddivRel = 1E3_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY disables the absolute divergence check.
    ! standard = SYS_INFINITY
    REAL(DP)                        :: ddivAbs = SYS_INFINITY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! RHS-vector is treated as zero if max(defect) < drhsZero
    REAL(DP)                        :: drhsZero = 1E-12_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Type of stopping criterion to use. One of the
    ! LINSOL_STOP_xxxx constants.
    INTEGER                    :: istoppingCriterion = LINSOL_STOP_STANDARD

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Minimum number of iterations top perform
    INTEGER                    :: nminIterations = 1

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Maximum number of iterations top perform
    INTEGER                    :: nmaxIterations = 50
    
    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Perform residual checks
    ! YES: check residuals (default), 
    ! NO: don't check residuals, simply perform as many iterations as 
    ! configured by nminIterations. Is typically set to NO for smoothing
    ! with a solver.
    INTEGER                    :: iresCheck = YES

    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Type of norm to use in the residual checking (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    INTEGER                    :: iresNorm = 2
    
    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Number of iterations that should be used to calculate
    ! the asymptotic convergtence rate, if the solver supports
    ! to generate an asymptotic convergence rate of the last
    ! couple of steps. =0: don't determine asymptotic convergence rate 
    INTEGER                    :: niteAsymptoticCVR = 0
    
    ! INPUT PARAMETER: Output level
    ! This determines the output level of the solver.
    ! =0: no output, =1: basic output, =2, extended output
    INTEGER                    :: ioutputLevel = 0

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK:
    ! Number of iterations to perform before printing out the
    ! norm of the residual to screen.
    ! =1: Print residual in every iteration
    INTEGER                    :: niteResOutput = 1
    
    ! INPUT PARAMETER: Default data type.
    ! This parameter prescribes the default data type that is used
    ! for temporary verctors inside the solvers. Solvers that need
    ! temporary vectors allocate them using this data type.
    ! The type identifier (either ST_SINGLE or ST_DOUBLE (default))
    ! should match the data type of the RHS/solution vectors,
    ! otherwise there might be performance loss!
    INTEGER                    :: cdefaultDataType = ST_DOUBLE

    ! INPUT PARAMETER: Solver subgroup
    ! By default, every solver in a soler tree belongs to solver
    ! subgroup 0. This means, the solver is initialised  the default
    ! call to linsol_initProblemStructure and linsol_initProblemData.
    ! By assigning a solver to a different subgroup , the solver will
    ! not initialise itself by a call to these routines. Instead,
    ! it will initialise itself if these routines are called for
    ! the specified subgroup number.
    ! Example: Multigrid solver, UMFPACK coarse grid solver,
    !  ILU(3) smoother.
    !  MG-Solver, UMFPACK solver: subgroup number := 0 (default)
    !  ILU(3) smoother:           subgroup number := 1
    ! linsol_initProblemStructure and linsol_initProblemData will
    ! initialise only MG and UMFPACK, not ILU(3). Instead,
    ! ILU(3) is initialised if linsol_initProblemStructure and 
    ! linsol_initProblemData are called with 'isolverSubgroup=1'.
    ! So, this mechanism provides a possibility to prevent some
    ! solvers from being initialised or allows to initialise
    ! some solvers manually.
    INTEGER                     :: isolverSubgroup = 0
    
    ! INPUT PARAMETER: t_matrixBlock structure that holds
    ! information about the linear system where to solve 
    ! (i.e. usually the finest level). 
    ! For multilevel-capable algorithms, information about the other 
    ! levels are saved in local solver-specific structures.
    TYPE(t_matrixBlock)                :: rsystemMatrix

    ! READ ONLY: Algorithm identifier. 
    ! One of the LINSOL_ALG_xxxx constants.
    ! Depending on the value, a solver-specific structure may be 
    ! assigned to this structure below.
    INTEGER                    :: calgorithm = LINSOL_ALG_UNDEFINED
    
    ! READ ONLY: Solver ability tag. 
    ! Bitfield. A combination of LINSOL_ABIL_xxxx
    ! flags that specify the ability of the actual solver (e.g.
    ! whether it can handle block-matrices or only scalar matrices,...).
    ! Calling a solver that is not able to handle a specific problem
    ! will cause an error by the framework.
    INTEGER(I32)                    :: ccapability = 0
    
    ! STATUS FOR ITERATIVE SOLVERS: Current iteration
    INTEGER                    :: icurrentIteration

    ! Pointer to a structure for the VANCA solver; NULL() if not set
    TYPE (t_linsolSubnodeVANCA), POINTER          :: p_rsubnodeVANCA       => NULL()
    
    ! Pointer to a structure for the Defect correction solver; NULL() if not set
    TYPE (t_linsolSubnodeDefCorr), POINTER        :: p_rsubnodeDefCorr     => NULL()

    ! Pointer to a structure for the BiCGStab solver; NULL() if not set
    TYPE (t_linsolSubnodeBiCGStab), POINTER       :: p_rsubnodeBiCGStab    => NULL()

    ! Pointer to a structure for the UMFPACK4 solver; NULL() if not set
    TYPE (t_linsolSubnodeUMFPACK4), POINTER       :: p_rsubnodeUMFPACK4    => NULL()

    ! Pointer to a structure for the Multigrid solver; NULL() if not set
    TYPE (t_linsolSubnodeMultigrid), POINTER      :: p_rsubnodeMultigrid   => NULL()

    ! Pointer to a structure for the ILU0 1x1 solver; NULL() if not set
    TYPE (t_linsolSubnodeILU01x1), POINTER        :: p_rsubnodeILU01x1     => NULL()
    
    ! Pointer to a structure for (M)ILUs 1x1 solver; NULL() if not set
    TYPE (t_linsolSubnodeMILUs1x1), POINTER       :: p_rsubnodeMILUs1x1    => NULL()

  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the scalar ILU(0) solver in
  ! FEAT style.
  
  TYPE t_linsolSubnodeILU01x1
  
    ! A pointer to the ILU0-decomposition of the main matrix.
    ! Shares the structure with the system matrix.
    TYPE(t_matrixScalar)              :: p_DiluDecomposition      
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the BiCGStab solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_linsolSubnodeDefCorr
  
    ! Temporary vector to use during the solution process
    TYPE(t_vectorBlock) :: rtempVector

    ! 2nd Temporary vector to use during the solution process
    TYPE(t_vectorBlock) :: rtempVector2

    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    TYPE(t_linsolNode), POINTER       :: p_rpreconditioner            => NULL()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:), POINTER      :: p_RfilterChain => NULL()
  
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the VANCA solver.
  
  TYPE t_linsolSubnodeVANCA
  
    ! Type of VANCA subsolver. One of the LINSOL_VANCA_xxxx constants.
    ! Specifies either thge general VANCA solver or a specialised version
    ! for improved speed.
    INTEGER             :: csubtypeVANCA
  
    ! For general VANCA (csubtypeVANCA=LINSOL_VANCA_GENERAL):
    ! Algorithm-specific structure.
    TYPE(t_vancaGeneral) :: rvancaGeneral
    
    ! For simple 2D-Saddle-Point Q1~/Q0 VANCA
    TYPE(t_vancaPointer2DSPQ1TQ0) :: rvanca2DSPQ1TQ0
  
    ! Temporary vector to use during the solution process
    TYPE(t_vectorBlock) :: rtempVector

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the BiCGStab solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_linsolSubnodeBiCGStab
  
    ! Temporary vectors to use during the solution process
    TYPE(t_vectorBlock), DIMENSION(6) :: RtempVectors

    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    TYPE(t_linsolNode), POINTER       :: p_rpreconditioner            => NULL()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:), POINTER      :: p_RfilterChain => NULL()
  
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the UMFPACK4 solver.
  
  TYPE t_linsolSubnodeUMFPACK4
  
    ! Control structure for UMFPACK4; contains parameter for the solver
    REAL(DP), DIMENSION(20) :: Dcontrol

    ! Handle for symbolic factorisation.
    ! This is not a FEAT-Handle!
    INTEGER(I32) :: isymbolic = 0

    ! Handle for numeric factorisation
    ! This is not a FEAT-Handle!
    INTEGER(I32) :: inumeric = 0
    
    ! Handle to a temporary vector for storing the solution
    TYPE(t_vectorBlock) :: rtempVector

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This saves information about ILU(k) matrices as defined in SPLIB.
  
  TYPE t_linsolSubnodeMILUs1x1

    ! Fill-in level for decomposition
    INTEGER :: ifill
    
    ! Relaxation factor for (M)ILU(s)
    REAL(DP) :: drelax

    ! Handle to array [1..*] of integer.
    ! Workspace containing the decomposed matrix.
    INTEGER :: h_Idata = ST_NOHANDLE

    ! Parameter 1 returned by SPLIB for factorised matrix
    INTEGER :: lu
    
    ! Parameter 2 returned by SPLIB for factorised matrix
    INTEGER :: jlu
    
    ! Parameter 3 returned by SPLIB for factorised matrix
    INTEGER :: ilup

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure collects timing information for the multigrid solver.

  TYPE t_linsolMGTiming
    
    ! Total time for prolongation
    REAL(DP)                           :: d_tmProlongation

    ! Total time for restriction
    REAL(DP)                           :: d_tmRestriction
  
    ! Total time for defect calculation
    REAL(DP)                           :: d_tmDefCalc
    
    ! Total time for smoothing
    REAL(DP)                           :: d_tmSmooth
    
    ! Total time for coarse grid solving
    REAL(DP)                           :: d_tmCoarseGridSolve

    ! Total time for coarse grid corrections
    REAL(DP)                           :: d_tmCorrection
  
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_linsolSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  TYPE t_linsolMGLevelInfo
  
    ! A level identifier. Can be set to identify the global number of the
    ! level, this structure refers to. Only for output purposes.
    INTEGER                        :: ilevel                 = 0
    
    ! Pointer to previous (lower) level or NULL() if the node is the 
    ! head of the list
    TYPE(t_linsolMGLevelInfo), POINTER  :: p_rprevLevel           => NULL()

    ! Pointer to next (higher) level or NULL() if the node is the 
    ! tail of the list
    TYPE(t_linsolMGLevelInfo), POINTER  :: p_rnextLevel           => NULL()
  
    ! t_matrixBlock structure that holds the system matrix.
    TYPE(t_matrixBlock)                :: rsystemMatrix
    
    ! A temporary vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_vectorBlock)                 :: rtempVector

    ! A RHS vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_vectorBlock)                 :: rrhsVector

    ! A solution vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_vectorBlock)                 :: rsolutionVector

    ! A pointer to the solver node for the presmoothing or NULL(),
    ! if no presmoother is used.
    TYPE(t_linsolNode), POINTER         :: p_rpresmoother         => NULL()

    ! A pointer to the solver node for the postsmoothing or NULL(),
    ! if no presmoother is used.
    TYPE(t_linsolNode), POINTER         :: p_rpostsmoother        => NULL()

    ! A pointer to the solver node for the coarse grid solver if
    ! the level corresponding to this structure represents the coarse
    ! grid. NULL() otherwise.
    TYPE(t_linsolNode), POINTER         :: p_rcoarseGridSolver    => NULL()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:),POINTER      :: p_RfilterChain => NULL()
    
    ! An interlevel projection structure that configures the transport
    ! of defect/solution vectors from one level to another.
    ! For the coarse grid, this structure is ignored.
    ! For a finer grid (e.g. level 4), this defines the grid transfer
    ! between the current and the lower level (so here between level 3 and
    ! level 4).
    TYPE(t_interlevelProjectionBlock)   :: rinterlevelProjection
    
    ! STATUS/INTERNAL: MG cycle information.
    ! Number of cycles to perform on this level.
    INTEGER                        :: ncycles

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of remaining cycles to perform on this level.
    INTEGER                        :: ncyclesRemaining
    
  END TYPE
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  TYPE t_linsolSubnodeMultigrid
  
    ! INPUT PARAMETER: Cycle identifier. 0=F-cycle, 1=V-cycle, 2=W-cycle
    INTEGER                       :: icycle                   = 0
    
    ! INPUT PARAMETER: Number of levels in the linked list of multigrid levels
    INTEGER                       :: nlevels                  = 0
    
    ! INPUT PARAMETER: Coarse grid correction structure for step length control. 
    ! Defines the algorithm for computing the optimal correction as well as the
    ! minimum and maximum step length ALPHAMIN/ALPHAMAX.
    ! The standard setting/initialisation is suitable for conforming elements.
    TYPE(t_coarseGridCorrection)  :: rcoarseGridCorrection
    
    ! Pointer to the head of the linked list of levels; corresponds
    ! to the lowest level.
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfoHead         => NULL()

    ! Pointer to the tail of the linked list of levels; corresponds
    ! to the highest level, i.e. the level where the system should be
    ! solved.
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfoTail         => NULL()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:),POINTER :: p_RfilterChain     => NULL()
  
    ! A temp vector for the prolongation/restriction.
    ! Memory for this vector is allocated in initStructure and released
    ! in doneStructure.
    TYPE(t_vectorScalar) :: rprjTempVector
    
  END TYPE
  
!</typeblock>

!</types>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

CONTAINS

  ! ***************************************************************************

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrices (rsolverNode,Rmatrices)
  
!<description>
  
  ! Initialises the system matrices in the solver node rsolverNode and
  ! in all nodes attached to it (preconditioners,...). For this purpose,
  ! the corresponding initialisation routine of each solver is called.
  !
  ! The last element in Rmatrices defines the system matrix of the
  ! linear system that is to solve. If there's a solver involved
  ! that needs matrices on multiple grids (e.g. multigrid), Rmatrices
  ! must contain the matrices of all levels. In this case, Rmatrices(1)
  ! defines the matrix on the coarsest level to solve, Rmatrices(2)
  ! the matrix of the first refinement etc. up to the level where
  ! to solve the problem.
  !
  ! The matrix structures passed in Rmatrices are copied to the system
  ! matrix structures in the solver node. If the system matrix changes
  ! (or more precisely, if handles or structural data of the system matrix
  ! changes), linsol_setMatrices has to be called again to update all
  ! the system matrix/matrices in the solver!
  
!</description>
  
!<input>
  
  ! Array of system matrices on all levels of the discretisation.
  ! This is passed through all initialisation routines, but actually used 
  ! only by the multigrid initialisation routine.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN) :: Rmatrices
  
!</input>
  
!<inputoutput>
  
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)             :: rsolverNode
  
!</inputoutput>
  
!</subroutine>

  INTEGER :: nblocks

  ! Copy the matrix structure on the finest level to the rsolverNode
  ! structure. This corresponds to the system we want to solve.
  rsolverNode%rsystemMatrix = Rmatrices(UBOUND(Rmatrices,1))
  
  ! Note in the system matrix structure that this matrix actually belongs to
  ! the application and not to the solver.
  ! (Don't care about empty sub-matrices here...)
  nblocks = rsolverNode%rsystemMatrix%ndiagBlocks
  rsolverNode%rsystemMatrix%RmatrixBlock(1:nblocks,1:nblocks)%imatrixSpec = &
    IOR(rsolverNode%rsystemMatrix%RmatrixBlock(1:nblocks,1:nblocks)%imatrixSpec,&
        LSYSSC_MSPEC_ISCOPY)
  
  ! Depending on the solver type, call the corresponding initialisation
  ! routine. For all single-grid solvers, pass the matrix on the highest
  ! level in RlinsysInfo. For multigrid we pass all matrices.
  ! That way the solvers can decide to do anything else with the matrix
  ! or to initialise their subsolvers.
  !
  ! For single grid solvers, the whole Rmatrices array is usually also
  ! passed to the initialisation routine. This allows to initialise
  ! the sub-nodes of each solver.
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (LINSOL_ALG_DEFCORR)
    CALL linsol_setMatrixDefCorr (rsolverNode,Rmatrices)
  CASE (LINSOL_ALG_UMFPACK4)
    ! UMFPACK needs no matrix initialisation routine, as it does not
    ! contain subsolvers. An attached matrix is processed in the 
    ! symbolical/numerical factorisation.
  CASE (LINSOL_ALG_VANCA)
    ! VANCA needs no matrix initialisation routine, as it does not
    ! contain subsolvers. An attached matrix is processed in the 
    ! symbolical/numerical factorisation.
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_setMatrixBiCGStab (rsolverNode, Rmatrices)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_setMatrixMultigrid (rsolverNode, Rmatrices)
  CASE DEFAULT
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_setOnelevelMatrixDirect (rsolverNode,rmatrix)
  
!<description>
  ! This is a special-case subroutine for the case that there is only
  ! one level to solve (e.g. when solving directly with UMFPACK).
  ! rmatrix is exactly one matrix of the level where to solve.
!</description>
  
!<input>
  ! System matrices of the discretisation.
  TYPE(t_matrixBlock), INTENT(IN),TARGET :: rmatrix
!</input>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>

!</subroutine>

    ! local variables
    
    ! t_linearSolverInfo structure array that receives the matrix pointer;
    ! it has length 1 here in this special situation.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    
    ! Copy the matrix pointer.
    Rmatrices(1) = rmatrix
    
    ! The rest of RlinsysInfo(.) will stay uninitialised - ok, there's no rest
    ! up to now :-)
    
    ! Call the standard initialisation routine
    CALL linsol_setMatrices (rsolverNode,Rmatrices)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructure (rsolverNode, ierror, isolverSubgroup)
  
!<description>
  
  ! Initialises the problem structure in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  ! The initialisation of the problem structure allowes the solver component
  ! to perform some 'precalculation', e.g. the UMFPACK4 or ILU solver can 
  ! perform a symbolical factorisation. The problem structure usually does
  ! not change during a simulation, except when the grid moves e.g..
  
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are initialised.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the structure-init routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_VANCA)
      CALL linsol_initStructureVANCA (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_initStructureDefCorr (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_initStructureUMFPACK4 (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MILUS1X1)
      ! No structure routine for (M)ILU(s)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_initStructureBiCGStab (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_initStructureMultigrid (rsolverNode,ierror,isubgroup)
    END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initData (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Initialises the problem data in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  ! The initialisation of the problem structure allowes the solver component
  ! to perform some 'precalculation', e.g. the UMFPACK4 or ILU solver can 
  ! perform a numerical factorisation. The problem structure usually does
  ! not change during a simulation, except when the grid moves e.g.
!</description>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are initialised.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the data-init routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_VANCA)
      CALL linsol_initDataVANCA (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_initDataDefCorr (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_initDataUMFPACK4 (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MILUS1X1)
      CALL linsol_initDataMILUs1x1 (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_initDataBiCGStab (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_initDataMultigrid (rsolverNode,ierror,isubgroup)
    END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_updateStructure (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  
  ! Reinitialises the problem structure in the solver node rsolverNode
  
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! For now, call the structure-done- and structure-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    
    CALL linsol_doneStructure (rsolverNode,isubgroup)
    CALL linsol_initStructure (rsolverNode,ierror,isubgroup)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_updateData (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Reinitialises the problem data in the solver node rsolverNode.
!</description>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<inputoutput>
  ! The solver node containing the solver confuguration
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>
  
!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! For now, call the data-done- and data-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    
    CALL linsol_doneData (rsolverNode, isubgroup)
    CALL linsol_initData (rsolverNode, isubgroup,ierror)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneData (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Releases the problem data in the solver node rsolverNode.
  
!</description>
  
!<inputoutput>
  
  ! The solver node containing the solver confuguration
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup

!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! by default, handle solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the data-done routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_VANCA)
      ! VANCA: no done-data-routine.
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_doneDataDefCorr (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_doneDataUMFPACK4 (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MILUS1X1)
      CALL linsol_doneDataMILUs1x1 (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_doneDataBiCGStab (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_doneDataMultigrid (rsolverNode,isubgroup)
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructure (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases the problem structure in the solver node rsolverNode.
  ! This is done by calling the appropriate routine of the
  ! actual solver.
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>
  
!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! by default, handle solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the data-done routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_VANCA)
      CALL linsol_doneStructureVANCA (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_doneStructureDefCorr (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_doneStructureUMFPACK4 (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MILUS1X1)
      ! No structure routine for (M)ILU(s)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_doneStructureBiCGStab (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_doneStructureMultigrid (rsolverNode,isubgroup)
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_releaseSolver (p_rsolverNode,bkeepSolverNode)
  
!<description>
  ! This routine releases the solver node rsolverNode and all attached
  ! sub-solver nodes (for preconditioning, smoothing etc. if specified)
  ! from the heap.
!</description>
  
!<input>
  ! OPTIONAL: If set to TRUE, the structure p_rsolverNode is not released
  ! from memory. If set to FALSE or not existent (the usual setting), the 
  ! structure p_rsolverNode will also be removed from the heap after 
  ! cleaning up.
  ! Remark: The subnodes of the solver (if there are any) are always
  ! removed from the heap!
  LOGICAL, INTENT(IN), OPTIONAL :: bkeepSolverNode
!</input>

!<inputoutput>
  ! The solver node which is to be released. If the node contains attached
  ! subsolvers (preconditioners, smoothers,...) they are also released.
  ! On return, rsolverNode is NULL().
  TYPE(t_linsolNode), POINTER     :: p_rsolverNode
!</inputoutput>
  
!</subroutine>

    IF (.NOT. ASSOCIATED(p_rsolverNode)) THEN
      ! Warning message, return
      PRINT *,'linsol_releaseSolver warning: Solver note not assigned!'
      RETURN
    END IF

    ! Depending on the solver type, call the corresponding done-routine
    SELECT CASE(p_rsolverNode%calgorithm)
    CASE (LINSOL_ALG_VANCA)
      CALL linsol_doneVANCA (p_rsolverNode)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_doneDefCorr (p_rsolverNode)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_doneUMFPACK4 (p_rsolverNode)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_doneBiCGStab (p_rsolverNode)
    CASE (LINSOL_ALG_MILUS1X1)
      CALL linsol_doneMILUs1x1 (p_rsolverNode)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_doneMultigrid (p_rsolverNode)
    CASE DEFAULT
    END SELECT
    
    ! Finally release the node itself (if we are allowed to).
    ! Deallocate the structure (if we are allowed to), finish.
    IF (.NOT. PRESENT(bkeepSolverNode)) THEN
      DEALLOCATE(p_rsolverNode)
    ELSE
      IF (.NOT. bkeepSolverNode) THEN
        DEALLOCATE(p_rsolverNode)
      END IF
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<function>
  
  LOGICAL FUNCTION linsol_testConvergence (rsolverNode, dvecNorm, rdef) RESULT(loutput)
  
!<description>
  ! Tests a defect vector rdef whether it is in a defined tolerance configured 
  ! in the solver node, so the iteration of an iterative solver can be
  ! can be treated as 'converged'.
  ! The iteration is treated as 'converged' if both, the relative and the
  ! absolute convergence criterion are fulfilled (or only one of them,
  ! respectively, if the other is switched off).
  !
  ! The solver must have initialised the 'dinitialDefect' variable
  ! of the solver structure for this routine to work properly!
!</description>
  
!<result>
  ! Boolean value. =TRUE if the convergence criterion is reached; 
  ! =FALSE otherwise.
!</result>
  
!<input>
  ! The solver node that contains the convergence criterion
  TYPE(t_linsolNode), INTENT(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks convergence depending on dvecNorm.
  TYPE(t_vectorBlock), INTENT(IN), OPTIONAL :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using dvecNorm.
  REAL(DP), INTENT(INOUT) :: dvecNorm
!</inputoutput>
  
!</function>

    ! Calculate the norm of the vector or take the one given
    ! as parameter
    IF (PRESENT(rdef)) THEN
      dvecNorm = lsysbl_vectorNorm (rdef,rsolverNode%iresNorm)
    END IF
    
    loutput = .TRUE.
    
    ! Absolute convergence criterion? Check the norm directly.
    IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
      IF (dvecNorm .GT. rsolverNode%depsAbs) THEN
        loutput = .FALSE.
      END IF
    END IF
    
    ! Relative convergence criterion? Multiply with initial residuum
    ! and check the norm. 
    IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
      IF (dvecNorm .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect) THEN
        loutput = .FALSE.
      END IF
    END IF
    
  END FUNCTION
  
  ! ***************************************************************************

!<function>
  
  LOGICAL FUNCTION linsol_testDivergence (rsolverNode, dvecNorm, rdef) RESULT(loutput)
  
!<description>
  ! Tests a defect vector rx whether it is out of a defined tolerance configured 
  ! in the solver node, so the iteration of an iterative solver can be
  ! can be treated as 'diverged'.
  ! The iteration is treated as 'diverged' if one criterion, the relative or the
  ! absolute divergence criterion is fulfilled.
  !
  ! The solver must have initialised the 'dinitialDefect' variable
  ! of the solver structure for this routine to work properly!
!</description>
  
!<result>
  ! Boolean value. =TRUE if the divergence criterion is reached; 
  ! =FALSE otherwise.
!</result>
  
!<input>
  ! The solver node that contains the divergence criterion
  TYPE(t_linsolNode), INTENT(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks divergence depending on dvecNorm.
  TYPE(t_vectorBlock), INTENT(IN), OPTIONAL :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and divergence is tested using dvecNorm.
  REAL(DP), INTENT(INOUT) :: dvecNorm
!</inputoutput>

!</function>

    ! Calculate the norm of the vector if not given
    ! as parameter
    IF (PRESENT(rdef)) THEN
      dvecNorm = lsysbl_vectorNorm (rdef,rsolverNode%iresNorm)
    END IF
    
    loutput = .FALSE.
    
    ! Absolute divergence criterion? Check the norm directly.
    IF (rsolverNode%ddivAbs .NE. SYS_INFINITY) THEN
     
      ! use NOT here - gives a better handling of special cases like NaN!
      IF ( .NOT. (dvecNorm .LE. rsolverNode%ddivAbs)) THEN
        loutput = .TRUE.
      END IF
      
    END IF
    
    ! Relative divergence criterion? Multiply with initial residuum
    ! and check the norm. 
    IF (rsolverNode%depsRel .NE. SYS_INFINITY) THEN
      IF ( .NOT. (dvecNorm .LE. rsolverNode%dinitialDefect*rsolverNode%ddivRel) ) THEN
        loutput = .TRUE.
      END IF
    END IF
  
  END FUNCTION
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precondDefect (rsolverNode,rd)
  
!<description>
  ! This routine applies the linear solver configured in rsolverNode
  ! as preconditioner to the defect vector rd.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
  !
  ! The vector rd must be compatible to the system matrix (same number
  ! of equations, same sorting strategy,...).
!</description>
  
!<inputoutput>
  ! The solver node containing the solver configuration
  TYPE(t_linsolNode), INTENT(INOUT)                :: rsolverNode
  
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)               :: rd
!</inputoutput>
  
!</subroutine>

    ! The only condition to this routine is that matrix and vector are compatible!
    CALL lsysbl_isMatrixCompatible(rd,rsolverNode%rsystemMatrix)

    ! Select the solver as configured in rsolverNode and let it perform
    ! the actual preconditioning task.
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_VANCA)
      CALL linsol_precVANCA (rsolverNode,rd)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_precDefCorr (rsolverNode,rd)
    CASE (LINSOL_ALG_JACOBI)
      CALL linsol_precJacobi (rsolverNode,rd)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_precUMFPACK4 (rsolverNode,rd)
    CASE (LINSOL_ALG_MILUS1x1)
      CALL linsol_precMILUS1x1 (rsolverNode,rd)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_precBiCGStab (rsolverNode,rd)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_precMultigrid (rsolverNode,rd)
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_solveAdaptively (rsolverNode,rx,rb,rtemp)
  
!<description>
  ! This routine starts the solution process to solve a linear system of the
  ! form $Ax=b$. The solver configuration must be given by rsolverNode.
  ! rb is the right hand side vector of the system. rx is an
  ! initial solution vector, which will be overwritten by the solution.
  !
  ! The solver does not implement any boundary conditions of the real
  ! system into the solution vector, this has to be done by the caller
  ! if necessary!
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<input>
  ! The RHS vector of the system
  TYPE(t_vectorBlock), INTENT(IN), TARGET            :: rb
!</input>
  
!<inputoutput>
  ! The solver node containing the solver configuration
  TYPE(t_linsolNode), INTENT(INOUT)                  :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  TYPE(t_vectorBlock), INTENT(INOUT)                 :: rx
  
  ! A temporary vector of the same size and structure as rx.
  TYPE(t_vectorBlock), INTENT(INOUT)                 :: rtemp
!</inputoutput>
  
!</subroutine>
    
    ! Method-specific remarks:
    ! The linear system $Ax=b$ is reformulated into a one-step defect-correction 
    ! approach
    !     $$ x  =  x_0  +  A^{-1}  ( b - A x_0 ) $$
    ! The standard solver P configured in rsovlverNode above is then used to 
    ! solve
    !     $$ Ay = b-Ax_0 $$
    ! and the solution $x$ is then calculated by $x=x+y$. 
    
    ! local variables
    
    ! Calculate the defect:
    ! To build (b-Ax), copy the RHS to the temporary vector
    ! and make a matrix/vector multiplication.
    CALL lsysbl_copyVector (rb,rtemp)
    CALL lsysbl_blockMatVec (rsolverNode%rsystemMatrix, rx, rtemp, -1.0_DP, 1.0_DP)
    
    ! Call linsol_precondDefect to solve the subproblem $Ay = b-Ax$.
    ! This overwrites rtemp with the correction vector.
    CALL linsol_precondDefect (rsolverNode,rtemp)

    ! Add the correction vector to the solution vector and release the memory.
    ! In case we have one... If the initial vector was assumed as zero, we don't
    ! have a correction vector, the result is directly in rx.
    
    ! Correct the solution vector: x=x+y
    CALL lsysbl_vectorLinearComb (rtemp,rx,1.0_DP,1.0_DP)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_convertToSmoother (rsolverNode,nsmoothingSteps,domega)
  
!<description>
  ! Converts a t_linsolNode to a smoother structure. A smoother is a solver
  ! that performs a fixed number of iterations without respecting any
  ! residuum. nsmoothingSteps is the number of steps the smoother should
  ! perform.
!</description>
  
!<input>
  ! Number of steps the smoother should perform
  INTEGER, INTENT(IN)          :: nsmoothingSteps
  
  ! OPTIONAL: Damping parameter.
  ! The parameter rsolverNode%domega is set to this value in order
  ! to set up the damping.
  REAL(DP), INTENT(IN), OPTIONAL :: domega
!</input>
  
!<inputoutput>
  ! Solver node which should be configured as smoother.
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>
  
!</subroutine>

    rsolverNode%depsRel = 0.0_DP
    rsolverNode%depsAbs = 0.0_DP
    rsolverNode%nminIterations = nsmoothingSteps
    rsolverNode%nmaxIterations = nsmoothingSteps
    rsolverNode%iresCheck      = NO
    
    IF (PRESENT(domega)) rsolverNode%domega = domega

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_initSolverGeneral (p_rsolverNode)
  
!<description>
  ! Creates a new solver node p_rsolverNode on the heap with default
  ! settings and returns a pointer to it.
!</description>
  
!<output>
  ! A pointer to a new solver node on the heap.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Allocate the node - the default initialisation will set most of the
    ! parameters in the structure.
    ALLOCATE(p_rsolverNode)

  END SUBROUTINE
  
! *****************************************************************************
! Routines for the Defect Correction iteration
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_Rfilter)
  
!<description>
  ! Creates a t_linsolNode solver structure for the defect correction iteration.
  ! The node can be used to directly solve a problem or to be attached 
  ! as solver or preconditioner to another solver structure. The node can be 
  ! deleted by linsol_releaseSolver.
  !
  ! The defect correction performs nmaxIterations iterations of the type
  !    $$ x_{n+1}  =  x_n  +  (b-Ax_n} $$
  ! with $x_0:=0$. 
  ! It's possible to include a damping parameter to this operation by 
  ! changing rsolverNode%domega to a value $\not =1$. In this case, the
  ! defect correction iteration changes to the Richardson iteration
  !
  !    $$ x_{n+1}  =  x_n  +  \omega(b-Ax_n} $$
  !
  ! By specifying an additional preconditioner, it's possible to perform
  ! the preconditioned defect correction iteration
  !
  !    $$ x_{n+1}  =  x_n  +  \omega P^{-1} (b-Ax_n} $$
  !
  ! By specifying an additional filter, the defect vector is filtered before
  ! each preconditioner step.
!</description>
  
!<input>
  ! OPTIONAL: A pointer to the solver structure of a solver that should be 
  ! used for preconditioning. If not given or set to NULL(), no preconditioning 
  ! will be used.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rpreconditioner
  
  ! OPTIONAL: A pointer to a filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the 
  ! iteration. If not given or set to NULL(), no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  TYPE(t_filterChain), DIMENSION(:), POINTER, OPTIONAL   :: p_Rfilter
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    CALL linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_DEFCORR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK    + &
                                LINSOL_ABIL_CHECKDEF + &
                                LINSOL_ABIL_USESUBSOLVER + &
                                LINSOL_ABIL_USEFILTER
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeDefCorr)
    
    ! Attach the preconditioner if given. 
    
    IF (PRESENT(p_rpreconditioner)) THEN 
      p_rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner => p_rpreconditioner
    END IF

    ! Attach the filter if given. 
    
    IF (PRESENT(p_Rfilter)) THEN
      p_rsolverNode%p_rsubnodeDefCorr%p_RfilterChain => p_Rfilter
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDefCorr (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Check if there's a preconditioner attached. If yes, release it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_releaseSolver(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    END IF
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeDefCorr)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixDefCorr (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_setMatrices (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                              Rmatrices)
    END IF

  END SUBROUTINE
  
! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureDefCorr (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_initStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                                isubgroup,ierror)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be initialised
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
    
    ! Intialisation. In our case: allocate temporary vectors for our data
    ! by using the associated matrix as template.
    ! That vectors are used in the defect correction so save the intermediate
    ! 'solution' vector.
    
    CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeDefCorr%rtempVector,.FALSE.,&
          rsolverNode%cdefaultDataType)

    CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeDefCorr%rtempVector2,.FALSE.,&
          rsolverNode%cdefaultDataType)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataDefCorr (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_initData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                            isubgroup,ierror)
    END IF
    
    ! Nothing to do here.
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataDefCorr (rsolverNode, isolverSUbgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_doneData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                            isubgroup)
    END IF
    
    ! Nothing to do here.
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureDefCorr (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_doneStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                                isubgroup)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be initialised
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

    ! Ok, we are the one to be released. We have temp-vectors to be released!
    IF (rsolverNode%p_rsubnodeDefCorr%rtempVector%NEQ .NE. 0) THEN
      CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector)
      CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector2)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precDefCorr (rsolverNode,rd)
  
!<description>
  ! Applies the Defect Correction preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: ireslength,ite,i
  REAL(DP) :: dres,dfr
  TYPE(t_vectorBlock), POINTER :: p_rx,p_rdef
  TYPE(t_linsolNode), POINTER :: p_rprecSubnode
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
  ! Damping parameter
  REAL(DP) :: domega

  ! The queue saves the current residual and the two previous residuals.
  REAL(DP), DIMENSION(32) :: Dresqueue
  
  ! The system matrix
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  INTEGER :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  LOGICAL bprec,bfilter
  
  ! The local subnode
  TYPE(t_linsolSubnodeDefCorr), POINTER :: p_rsubnode

    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeDefCorr
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    IF ((rd%NEQ .EQ. 0) .OR. (p_rmatrix%NEQ .EQ. 0) .OR. &
        (p_rmatrix%NEQ .NE. rd%NEQ) ) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate

    ireslength = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
    
    ! Damping parameter
    domega = rsolverNode%domega
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    bfilter = ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_RfilterChain)
    
    IF (bprec) THEN
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    END IF
    IF (bfilter) THEN
      p_RfilterChain => p_rsubnode%p_RfilterChain
    END IF

    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! We want to perform <= nmaxIterations of the form
    !
    !    $$ x_{n+1}  =  x_n  +  \omega P^{-1} (b-Ax) $$
    !
    ! At first, set up our temporary vectors, which holds the current
    ! 'solution' and the current 'defect'. We already allocated them 
    ! during the initialisation phase - now we want to use them!
    
    p_rx   => rsolverNode%p_rsubnodeDefCorr%rtempVector
    p_rdef => rsolverNode%p_rsubnodeDefCorr%rtempVector2
    
    ! The vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    CALL lsysbl_assignDiscretIndirect (rd,p_rx)
    CALL lsysbl_assignDiscretIndirect (rd,p_rdef)
  
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    CALL lsysbl_clearVector (p_rx)
  
    ! Copy our RHS rd to p_rdef. As the iteration vector is 0, this
    ! is also our initial defect.

    CALL lsysbl_copyVector(rd,p_rdef)
    IF (bfilter) THEN
      ! Apply the filter chain to the vector
      CALL filter_applyFilterChainVec (p_rdef, p_RfilterChain)
    END IF
    IF (bprec) THEN
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      CALL linsol_precondDefect (p_rprecSubnode,p_rdef)
    END IF
    
    ! Get the norm of the residuum
    dres = lsysbl_vectorNorm (p_rdef,rsolverNode%iresNorm)
    IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
              (dres .LE. 1E99_DP))) dres = 0.0_DP

    ! Initialize starting residuum
      
    rsolverNode%dinitialDefect = dres

    ! initialize the queue of the last residuals with RES

    Dresqueue = dres

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    IF ( rsolverNode%dinitialDefect .LT. rsolverNode%drhsZero ) THEN
     
      ! final defect is 0, as initialised in the output variable above

      CALL lsysbl_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    ELSE

      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        PRINT *,&
          'DefCorr: Iteration ',0,',  !!RES!! = ',rsolverNode%dinitialDefect
      END IF

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite
        
        ! In p_rdef, we now have the current residuum $P^{-1} (b-Ax)$.
        ! Add it (damped) to the current iterate p_x to get
        !   $$ x  :=  x  +  \omega P^{-1} (b-Ax) $$

        CALL lsysbl_vectorLinearComb (p_rdef ,p_rx,domega,1.0_DP)

        ! Calculate the residuum for the next step : (b-Ax)
        CALL lsysbl_copyVector (rd,p_rdef)
        CALL lsysbl_blockMatVec (p_rmatrix, p_rx,p_rdef, -1.0_DP,1.0_DP)
        IF (bfilter) THEN
          ! Apply the filter chain to the vector
          CALL filter_applyFilterChainVec (p_rdef, p_RfilterChain)
        END IF
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL linsol_precondDefect (p_rprecSubnode,p_rdef)
        END IF
        
        ! Get the norm of the new (final?) residuum
        dfr = lsysbl_vectorNorm (p_rdef,rsolverNode%iresNorm)
     
        ! Shift the queue with the last residuals and add the new
        ! residual to it. Check length of ireslength to be larger than
        ! 0 as some compilers might produce Floating exceptions
        ! otherwise! (stupid pgf95)
        IF (ireslength .GT. 0) &
          dresqueue(1:ireslength) = EOSHIFT(dresqueue(1:ireslength),1,dfr)

        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        IF (linsol_testDivergence(rsolverNode,dfr)) THEN
          PRINT *,'DefCorr: Solution diverging!'
          rsolverNode%iresult = 1
          EXIT
        END IF
     
        ! At least perform nminIterations iterations
        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          IF (linsol_testConvergence(rsolverNode,dfr)) EXIT
          
        END IF

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          !WRITE (MTERM,'(A,I7,A,D25.16)') 
          PRINT *,'DefCorr: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
        END IF

      END DO

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) &
        ite = rsolverNode%nmaxIterations

      ! Finish - either with an error or if converged.
      ! Print the last residuum.

      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        !WRITE (MTERM,'(A,I7,A,D25.16)') 
        PRINT *,'DefCorr: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
      END IF

    END IF

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    CALL lsysbl_copyVector (p_rx,rd)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
    
      ! Calculate asymptotic convergence rate
    
      IF (dresqueue(1) .GE. 1E-70_DP) THEN
        I = MAX(1,MIN(rsolverNode%iiterations,ireslength-1))
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(I,DP))
      END IF

      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :

      IF (rsolverNode%dfinalDefect .GT. rsolverNode%drhsZero) THEN
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/REAL(rsolverNode%iiterations,DP))
      END IF
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A)') 
        PRINT *,'DefCorr statistics:'
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A,I5)')     'Iterations              : ',
    !*            IPARAM(OITE)
        PRINT *,'Iterations              : ',rsolverNode%iiterations
        !WRITE (MTERM,'(A,D24.12)') '!!INITIAL RES!!         : ',
    !*            DPARAM(ODEFINI)
        PRINT *,'!!INITIAL RES!!         : ',rsolverNode%dinitialDefect
        !WRITE (MTERM,'(A,D24.12)') '!!RES!!                 : ',
    !*            DPARAM(ODEFFIN)
        PRINT *,'!!RES!!                 : ',rsolverNode%dfinalDefect
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
    !*              rparam%dfinalDefect / rparam%dinitialDefect
          PRINT *,'!!RES!!/!!INITIAL RES!! : ',&
                  rsolverNode%dfinalDefect / rsolverNode%dinitialDefect
        ELSE
          !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
    !*              0D0
          PRINT*,'!!RES!!/!!INITIAL RES!! : ',0.0_DP
        END IF
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A,D24.12)') 'Rate of convergence     : ',
    !*            DPARAM(ORHO)
        PRINT *,'Rate of convergence     : ',rsolverNode%dconvergenceRate

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'DefCorr: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'DefCorr: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
        PRINT *,&
              'DefCorr: Iterations/Rate of convergence: ',&
              rsolverNode%iiterations,' /',rsolverNode%dconvergenceRate
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the Jacobi solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initJacobi (p_rsolverNode, domega)
  
!<description>
  ! Creates a t_linsolNode solver structure for the Jacobi solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  !
  ! This Jacobi solver has no done-routine as there is no dynamic information
  ! allocated. The given damping parameter domega is saved to the solver
  ! structure rsolverNode%domega.
!</description>
  
!<input>
  ! OPTIONAL: Damping parameter. Is saved to rsolverNode%domega if specified.
  REAL(DP), OPTIONAL :: domega
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    
    CALL linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_JACOBI 
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
    
      ! No subnode for Jacobi. Only save domega to the structure.
    IF (PRESENT(domega)) THEN
      p_rsolverNode%domega = domega
    END IF
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precJacobi (rsolverNode,rd)
  
!<description>
  ! Applies the Jacobi preconditioner $D \in A$ to the defect 
  ! vector rd and solves $Dd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Jacobi solver
  TYPE(t_linsolNode), INTENT(INOUT),TARGET  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER iblock
    INTEGER(PREC_VECIDX) :: ieq
    
    TYPE (t_matrixScalar), POINTER :: p_rmatrix
    INTEGER (PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiag
    REAL(DP) :: dlocOmega
    REAL(SP) :: flocOmega
    REAL(DP), DIMENSION(:), POINTER :: p_Dvector, p_Dmatrix
    REAL(SP), DIMENSION(:), POINTER :: p_Fvector, p_Fmatrix
    
    ! Loop through all blocks. Each block corresponds to one
    ! diagonal block in the matrix.
    DO iblock = 1,rd%nblocks
      ! Get the matrix
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
      
      ! Now we have to make some decisions. At first, which matrix
      ! structure do we have?
      SELECT CASE (p_rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
        ! In case of structure 9, Kdiagonal points to the diagonal elements.
        ! In case of structure 7, Kld points to the diagonal elements
        IF (p_rmatrix%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
          CALL storage_getbase_int(p_rmatrix%h_Kdiagonal,p_Kdiag)
        ELSE
          CALL storage_getbase_int(p_rmatrix%h_Kld,p_Kdiag)
        END IF
      
        ! Get the omega parameter for the matrix. Don't forget the scaling 
        ! factor in the matrix!
        dlocomega = rsolverNode%domega * p_rmatrix%dScaleFactor

        ! Now, which data format do we have? Single or double?
        SELECT CASE (p_rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Get the matrix data arrays
          CALL storage_getbase_double (p_rmatrix%h_DA,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            DO ieq = 1,rd%NEQ
              p_Dvector(ieq) = dlocOmega * p_Dvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
            END DO
            
          CASE (ST_SINGLE)
            ! Get the data array
            CALL lsysbl_getbase_single (rd,p_Fvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            DO ieq = 1,rd%NEQ
              p_Dvector(ieq) = dlocOmega * p_Fvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
            END DO

          CASE DEFAULT
            PRINT *,'Jacobi: Unsupported vector format.'
            STOP
          END SELECT
          
        CASE (ST_SINGLE)

          ! Get the matrix data arrays
          CALL storage_getbase_double (p_rmatrix%h_Da,p_Fmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            DO ieq = 1,rd%NEQ
              p_Dvector(ieq) = dlocOmega * p_Dvector(ieq) / p_Fmatrix(p_Kdiag(ieq)) 
            END DO
            
          CASE (ST_SINGLE)
            ! Get the data array
            CALL lsysbl_getbase_single (rd,p_Fvector)
            
            ! Multiplication with Omega can be speeded up as we use
            ! sigle-precision only.
            flocOmega = REAL(dlocOmega,SP)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            DO ieq = 1,rd%NEQ
              p_Dvector(ieq) = flocOmega * p_Fvector(ieq) / p_Fmatrix(p_Kdiag(ieq))
            END DO

          CASE DEFAULT
            PRINT *,'Jacobi: Unsupported vector format.'
            STOP
          END SELECT

        CASE DEFAULT
          PRINT *,'Jacobi: Unsupported matrix format.'
          STOP
        END SELECT
      
      CASE DEFAULT
        PRINT *,'Jacobi: Unsupported matrix format.'
        STOP
      END SELECT
      
    END DO
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the VANCA CC2D solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initVANCA (p_rsolverNode,domega,csubtypeVANCA)
  
!<description>
  ! Creates a t_linsolNode solver structure for the VANCA solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  !
  ! VANCA is somehopw a special type of solver. There is one general VANCA
  ! solver, which applies to many situations, but which is quite slow.
  ! For higher performance, system-specific VANCA subsolvers must be created
  ! with hardcoded treatment of matrices/vectors!
  !
  ! For this purpose, there's the csubtypeVANCA flag. If this flag is present,
  ! it indicates a special VANCA variant for a special type of situation.
  ! the application must ensure then that the correct matrix/vector structure
  ! is used when solving the system with this special subtype, as each subtype
  ! implements a hardwired processing of matrix/vector data!
!</description>
  
!<input>
  ! OPTIONAL: Damping parameter. Is saved to rsolverNode%domega if specified.
  REAL(DP), OPTIONAL :: domega
  
  ! OPTIONAL: VANCA subtype.
  ! If not present, the standard VANCA solver is used.
  ! If present, this is one of the LINSOL_VANCA_xxxx flags that indicate a 
  ! special VANCA variant for higher performance.
  INTEGER, INTENT(IN), OPTIONAL       :: csubtypeVANCA
!</input>  
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>
  
  ! Create a default solver structure
  
  CALL linsol_initSolverGeneral(p_rsolverNode)
  
  ! Initialise the type of the solver
  p_rsolverNode%calgorithm = LINSOL_ALG_VANCA
  
  ! Initialise the ability bitfield with the ability of this solver:
  p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                              LINSOL_ABIL_DIRECT
  
  ! Allocate the subnode for VANCA.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(p_rsolverNode%p_rsubnodeVANCA)
  
  ! Initialise the VANCA subtype.
  p_rsolverNode%p_rsubnodeVANCA%csubtypeVANCA = LINSOL_VANCA_GENERAL
  IF (PRESENT(csubtypeVANCA)) THEN
    p_rsolverNode%p_rsubnodeVANCA%csubtypeVANCA = csubtypeVANCA
  END IF
  
  IF (PRESENT(domega)) p_rsolverNode%domega = domega
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureVANCA (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Memory allocation for the VANCA solver.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANCA solver
  TYPE(t_linsolNode), INTENT(INOUT),TARGET  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup

  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Stop if there's no matrix assigned
  
  IF (rsolverNode%rsystemMatrix%NEQ .EQ. 0) THEN
    PRINT *,'Error: No matrix associated!'
    STOP
  END IF
  
  ! If isubgroup does not coincide with isolverSubgroup from the solver
  ! structure, skip the rest here.
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

  ! Allocate a temporary vector
  CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
        rsolverNode%p_rsubnodeVANCA%rtempVector,.FALSE.,&
        rsolverNode%cdefaultDataType)
        
  ! Which VANCA solver do we actually have?
  SELECT CASE (rsolverNode%p_rsubnodeVANCA%csubtypeVANCA)
  CASE (LINSOL_VANCA_GENERAL)
    ! Special initialisation of the VANCA-general structure.
    CALL vanca_initGeneralVanca(rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANCA%rvancaGeneral)
  END SELECT
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataVANCA (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs final preparation of the VANCA solver subtype.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANCA solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup

  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Stop if there's no matrix assigned
  
  IF (rsolverNode%rsystemMatrix%NEQ .EQ. 0) THEN
    PRINT *,'Error: No matrix associated!'
    STOP
  END IF
  
  ! Which VANCA subtype do we have? Only in some variants there's
  ! something to do...
  SELECT CASE (rsolverNode%p_rsubnodeVANCA%csubtypeVANCA)
  CASE (LINSOL_VANCA_2DSPQ1TQ0)
    ! Special initialisation of a VANCA structure.
    CALL vanca_init2DSPQ1TQ0simple(rsolverNode%rsystemMatrix,&
                                   rsolverNode%p_rsubnodeVANCA%rvanca2DSPQ1TQ0)
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureVANCA (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases temporary memory of the VANCA solver.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANCA solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! If isubgroup does not coincide with isolverSubgroup from the solver
  ! structure, skip the rest here.
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
  
  ! Release the temp vector and associated data if associated
  IF (rsolverNode%p_rsubnodeVANCA%rtempVector%NEQ .NE. 0) THEN
  
    CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeVANCA%rtempVector)

    ! Which VANCA solver do we actually have?
    SELECT CASE (rsolverNode%p_rsubnodeVANCA%csubtypeVANCA)
    CASE (LINSOL_VANCA_GENERAL)
      ! Special cleanup of the VANCA-general structure.
      CALL vanca_doneGeneralVanca(rsolverNode%p_rsubnodeVANCA%rvancaGeneral)
    END SELECT

  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneVANCA (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the VANCA solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of VANCA which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  INTEGER :: isubgroup
  
  ! The done routine always releases everything.
  isubgroup = rsolverNode%isolverSubgroup

  ! Release temporary data.
  CALL linsol_doneStructureVANCA (rsolverNode, isubgroup)
  
  DEALLOCATE(rsolverNode%p_rsubnodeVANCA)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precVANCA (rsolverNode,rd)
  
!<description>
  ! Applies the VANCA preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANCA solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    TYPE (t_vectorBlock), POINTER :: p_rvector
    TYPE (t_matrixBlock), POINTER :: p_rmatrix
    REAL(DP) :: domega

    ! Status reset
    rsolverNode%iresult = 0

    ! Getch some information
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    IF ((rd%NEQ .EQ. 0) .OR. (p_rmatrix%NEQ .EQ. 0) .OR. &
        (p_rmatrix%NEQ .NE. rd%NEQ) ) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Damping parameter
    domega = rsolverNode%domega
      
    ! Get our temporary vector.
    p_rvector => rsolverNode%p_rsubnodeVANCA%rtempVector

    ! The vector share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vector.
    CALL lsysbl_assignDiscretIndirect (rd,p_rvector)
  
    ! Clear our solution vector
    CALL lsysbl_clearVector (p_rvector)

    ! Choose the correct (sub-)type of VANCA to call.
    SELECT CASE (rsolverNode%p_rsubnodeVANCA%csubtypeVANCA)
    
    CASE (LINSOL_VANCA_GENERAL)
      CALL vanca_general (rsolverNode%p_rsubnodeVANCA%rvancaGeneral, &
                          p_rvector, rd, domega)

    CASE (LINSOL_VANCA_2DSPQ1TQ0)
      CALL vanca_2DSPQ1TQ0simple (rsolverNode%p_rsubnodeVANCA%rvanca2DSPQ1TQ0, &
                                  p_rvector, rd, domega)
      
    CASE DEFAULT
      PRINT *,'Unknown VANCA variant!'
      STOP
    END SELECT
    
    ! Copy the solution vector to rd - it's our preconditioned defect now.
    CALL lsysbl_copyVector (p_rvector,rd)
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the UMFPACK4 solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initUMFPACK4 (rsolverNode)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the UMFPACK4 solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  
!</description>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: rsolverNode
   
!</output>
  
!</subroutine>

  ! Create a default solver structure
  
  CALL linsol_initSolverGeneral(rsolverNode)
  
  ! Initialise the type of the solver
  rsolverNode%calgorithm = LINSOL_ALG_UMFPACK4
  
  ! Initialise the ability bitfield with the ability of this solver:
  rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_DIRECT
  
  ! Allocate the subnode for UMFPACK4.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(rsolverNode%p_rsubnodeUMFPACK4)
  
  ! Initialize the UMFPACK4 control structure:
        
  CALL UMF4DEF(rsolverNode%p_rsubnodeUMFPACK4%Dcontrol)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureUMFPACK4 (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Performs a symbolic factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT),TARGET  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  TYPE(t_matrixScalar), POINTER :: p_rmatrix
  TYPE(t_matrixScalar) :: rtempMatrix
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  TYPE(t_matrixBlock), TARGET :: rmatrixLocal

  ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
  ! of a call to the solver routines.
  REAL(DP), DIMENSION(90) :: Dinfo
  
  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Stop if there's no matrix assigned
  
  IF (rsolverNode%rsystemMatrix%NEQ .EQ. 0) THEN
    PRINT *,'Error: No matrix associated!'
    STOP
  END IF
  
  ! If isubgroup does not coincide with isolverSubgroup from the solver
  ! structure, skip the rest here.
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

  ! Check out that we can handle the matrix.
  ! Check out that we can handle the matrix.
  IF (rsolverNode%rsystemMatrix%ndiagBlocks .NE. 1) THEN
    ! We have to create a global matrix first!
    CALL glsys_assembleGlobal (rsolverNode%rsystemMatrix,rmatrixLocal, &
                               .TRUE.,.TRUE.)
    p_rmatrix => rmatrixLocal%RmatrixBlock(1,1)
  ELSE
    p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock (1,1)
  END IF

  IF (p_rmatrix%cdataType .NE. ST_DOUBLE) THEN
    PRINT *,'UMFPACK can only handle double precision matrices!'
    STOP
  END IF

  IF (p_rmatrix%dScaleFactor .NE. 1.0_DP) THEN
    PRINT *,'UMFPACK cannot handle scaled matrices!'
    STOP
  END IF
  
  SELECT CASE (p_rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
    ! Format 9 is exactly the UMFPACK matrix.
    ! Make a copy of the matrix structure, but use the same matrix entries.
    CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                  LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
  CASE (LSYSSC_MATRIX7)
    ! For format 7, we have to modify the matrix slightly.
    ! Make a copy of the whole matrix:
    CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                  LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    ! Resort the entries to put the diagonal entry to the correct position.
    ! This means: Convert the structure-7 matrix to a structure-9 matrix:
    PRINT *,'UMFPACK: Convert 7->9 matrix not implemented!'
    STOP
  END SELECT
  
  ! Modify Kcol/Kld of the matrix. Subtract 1 to get the 0-based.
  CALL lsyssc_addIndex (rtempMatrix%h_Kcol,-1)
  CALL lsyssc_addIndex (rtempMatrix%h_Kld,-1)
  
  ! Get the data arrays.
  CALL storage_getbase_int (rtempMatrix%h_Kcol,p_Kcol)
  CALL storage_getbase_int (rtempMatrix%h_Kld,p_Kld)
  CALL storage_getbase_double (rtempMatrix%h_DA,p_DA)
  
  ! Perform a symbolic factorization...
  CALL UMF4SYM(rtempMatrix%NEQ,rtempMatrix%NEQ,p_Kld,p_Kcol,p_Da, &
               rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
               rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
               Dinfo)
               
  ! Check Dinfo(1) if there is an error
  SELECT CASE (INT(Dinfo(1)))
  CASE (0)
    ! ok.
  CASE (1)
    ! Singular matrix
    ierror = LINSOL_ERR_SINGULAR
  CASE (-1)
    ! no memory
    ierror = LINSOL_ERR_NOMEMORY
  CASE DEFAULT
    ! don't know what went wrong
    ierror = LINSOL_ERR_INITERROR
  END SELECT

  ! Throw away the temporary matrix/matrices
  CALL lsyssc_releaseMatrix (rtempMatrix)
  IF (rsolverNode%rsystemMatrix%ndiagBlocks .NE. 1) THEN
    CALL lsysbl_releaseMatrix (rmatrixLocal)
  END IF
  
  ! Allocate a temporary vector
  CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
        rsolverNode%p_rsubnodeUMFPACK4%rtempVector,.FALSE.,&
        rsolverNode%cdefaultDataType)
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataUMFPACK4 (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs a numeric factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  TYPE(t_matrixScalar), POINTER :: p_rmatrix
  TYPE(t_matrixBlock), TARGET :: rmatrixLocal
  TYPE(t_matrixScalar) :: rtempMatrix
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  REAL(DP), DIMENSION(:), POINTER :: p_DA

  ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
  ! of a call to the solver routines.
  REAL(DP), DIMENSION(90) :: Dinfo
  
  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Stop if there's no matrix assigned
  
  IF (rsolverNode%rsystemMatrix%NEQ .EQ. 0) THEN
    PRINT *,'Error: No matrix associated!'
    STOP
  END IF
  
  ! If isubgroup does not coincide with isolverSubgroup from the solver
  ! structure, skip the rest here.
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
  
  ! Check out that we can handle the matrix.
  IF (rsolverNode%rsystemMatrix%ndiagBlocks .NE. 1) THEN
    ! We have to create a global matrix first!
    CALL glsys_assembleGlobal (rsolverNode%rsystemMatrix,rmatrixLocal, &
                               .TRUE.,.TRUE.)
    p_rmatrix => rmatrixLocal%RmatrixBlock(1,1)
  ELSE
    p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock (1,1)
  END IF

  IF (p_rmatrix%cdataType .NE. ST_DOUBLE) THEN
    PRINT *,'UMFPACK can only handle double precision matrices!'
    STOP
  END IF

  IF (p_rmatrix%dScaleFactor .NE. 1.0_DP) THEN
    PRINT *,'UMFPACK cannot handle scaled matrices!'
    STOP
  END IF

  SELECT CASE (p_rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
    ! Format 9 is exactly the UMFPACK matrix.
    ! Make a copy of the matrix structure, but use the same matrix entries.
    CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                  LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
  CASE (LSYSSC_MATRIX7)
    ! For format 7, we have to modify the matrix slightly.
    ! Make a copy of the whole matrix:
    CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                 LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    ! Resort the entries to put the diagonal entry to the correct position.
    ! This means: Convert the structure-7 matrix to a structure-9 matrix:
    CALL lsyssc_convertMatrix (p_rmatrix,LSYSSC_MATRIX9)
  END SELECT

  ! Modify Kcol/Kld of the matrix. Subtract 1 to get the 0-based.
  CALL lsyssc_addIndex (rtempMatrix%h_Kcol,-1)
  CALL lsyssc_addIndex (rtempMatrix%h_Kld,-1)
  
  ! Get the data arrays.
  CALL storage_getbase_int (rtempMatrix%h_Kcol,p_Kcol)
  CALL storage_getbase_int (rtempMatrix%h_Kld,p_Kld)
  CALL storage_getbase_double (rtempMatrix%h_DA,p_DA)
  
  ! Perform a numeric factorization...
  CALL UMF4NUM(p_Kld,p_Kcol,p_Da, &
          rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
          rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
          rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
          Dinfo)
          
  ! Check Dinfo(1) if there is an error
  SELECT CASE (INT(Dinfo(1)))
  CASE (0)
    ! ok.
  CASE (1)
    ! Singular matrix
    ierror = LINSOL_ERR_SINGULAR
  CASE (-1)
    ! no memory
    ierror = LINSOL_ERR_NOMEMORY
  CASE DEFAULT
    ! don't know what went wrong
    ierror = LINSOL_ERR_INITERROR
  END SELECT

  ! Throw away the temporary matrix/matrices
  CALL lsyssc_releaseMatrix (rtempMatrix)
  IF (rsolverNode%rsystemMatrix%ndiagBlocks .NE. 1) THEN
    CALL lsysbl_releaseMatrix (rmatrixLocal)
  END IF
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataUMFPACK4 (rsolverNode,isolverSubgroup)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! If isubgroup does not coincide with isolverSubgroup from the solver
  ! structure, skip the rest here.
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
  
  ! Release the numerical factorisation if associated
  IF (rsolverNode%p_rsubnodeUMFPACK4%inumeric .NE. 0) THEN
    CALL UMF4FNUM(rsolverNode%p_rsubnodeUMFPACK4%inumeric)
    rsolverNode%p_rsubnodeUMFPACK4%inumeric = 0
  END IF  
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureUMFPACK4 (rsolverNode, isolverSUbgroup)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! If isubgroup does not coincide with isolverSubgroup from the solver
  ! structure, skip the rest here.
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
  
  ! Release the symbolical factorisation if associated
  IF (rsolverNode%p_rsubnodeUMFPACK4%isymbolic .NE. 0) THEN
    CALL UMF4FSYM(rsolverNode%p_rsubnodeUMFPACK4%isymbolic)
    rsolverNode%p_rsubnodeUMFPACK4%isymbolic = 0
  END IF
  
  ! Release the temp vector if associated
  IF (rsolverNode%p_rsubnodeUMFPACK4%rtempVector%NEQ .NE. 0) THEN
    CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeUMFPACK4%rtempVector)
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneUMFPACK4 (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  INTEGER :: isubgroup
  
  ! The done routine always releases everything.
  isubgroup = rsolverNode%isolverSubgroup

  ! Release symbolical and numerical factorisation if still associated...
  CALL linsol_doneDataUMFPACK4 (rsolverNode, isubgroup)
  CALL linsol_doneStructureUMFPACK4 (rsolverNode, isubgroup)
  
  DEALLOCATE(rsolverNode%p_rsubnodeUMFPACK4)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precUMFPACK4 (rsolverNode,rd)
  
!<description>
  ! Applies UMFPACK preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>
 
  ! local variables
  INTEGER KSYS
  REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Db
  TYPE(t_vectorBlock), POINTER :: p_rb
  TYPE(t_linsolSubnodeUMFPACK4), POINTER :: p_rsubnode

  ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
  ! of a call to the solver routines.
  REAL(DP), DIMENSION(90) :: Dinfo
  
    ! Check that our RHS db is double precision - UMFPACK supports only
    ! this.
    IF (rd%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'UMFPACK only supports double precision vectors!'
      STOP
    END IF

    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeUMFPACK4
    p_rb   => p_rsubnode%rtempVector
    
    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    CALL lsysbl_assignDiscretIndirect (rd,p_rb)
    
    ! Copy the RHS rd to the temp vector; it will be overwritten
    ! by the solution vector
    CALL lsysbl_copyVector (rd,p_rb)

    ! Get the RHS and solution vector data
    CALL lsysbl_getbase_double(rd,p_Dx)
    CALL lsysbl_getbase_double(p_rb,p_Db)

    ! Solve the system
    ! Solve the system. Note that UMFPACK expects the matrix in
    ! CSR format, which is transposed to our matrix format 9 --
    ! So we solve the transposed system:
    KSYS = 1
    CALL UMF4SOL(KSYS,p_Dx,p_Db,rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
                 rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,Dinfo)
                
    ! Check the solver status
    SELECT CASE (INT(Dinfo(1)))
    CASE (0) 
      ! All ok.
      rsolverNode%iiterations = 1
      rsolverNode%dconvergenceRate = 0.0_DP
    CASE DEFAULT
      ! We had an error. Don't know which one.
      rsolverNode%iresult = -1
    END SELECT
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the (M)ILUs 1x1 solver
! *****************************************************************************
! This is a 1-step solver usually used for preconditioning: x_1 = C^-1 x_0.
! As we use SPLIB, we have only a numerical factorisation, no symbolical
! factorisation.

!<subroutine>
  
  SUBROUTINE linsol_initMILUs1x1 (rsolverNode,ifill,drelax)
  
!<description>
  ! Creates a t_linsolNode solver structure for the (M)ILUs 1x1 solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! Fill-in level for factorisation.
  ! 0=ILU(0), 1=ILU(1),...
  INTEGER, INTENT(IN) :: ifill
  
  ! Relaxation parameter.
  ! 0.0=ILU(s), 1.0=MILU(s), between 0.0 and 1.0=multiplier to use before
  ! adding discarded fill to the diagonal.
  REAL(DP), INTENT(IN) :: drelax
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: rsolverNode
!</output>

!</subroutine>

    ! Create a default solver structure
    
    CALL linsol_initSolverGeneral(rsolverNode)
    
    ! Normally this solver iterates only once (used e.g. in a preconditioner).
    ! This is the default value. The application can set it higher to
    ! apply it multiple times.
    rsolverNode%nminIterations = 1
    rsolverNode%nmaxIterations = 1
    
    ! Initialise the type of the solver
    rsolverNode%calgorithm = LINSOL_ALG_MILUS1x1
    
    ! Initialise the ability bitfield with the ability of this solver.
    ! The solver only solves scalar (1x1) systems.
    ! The solver is a 1-step solver, i.e. does not support to solve
    ! by iteration (it basically perfoms only a forward-backward-substitution
    ! with the (M)ILUs matrix).
    rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_DIRECT
    
    ! Allocate the subnode for MILUs.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(rsolverNode%p_rsubnodeMILUs1x1)
    
    ! Save the (M)ILU(s) parameters to the structure
    rsolverNode%p_rsubnodeMILUs1x1%ifill = ifill
    rsolverNode%p_rsubnodeMILUs1x1%drelax = drelax
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataMILUs1x1 (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs a numeric factorisation on the assigned matrix, i.e.
  ! computes the (M)ILU(s) decomposition.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT),TARGET  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup,mneed,ierr,maxstr
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_matrixScalar), POINTER :: p_rmatrixSc
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld, p_Kcol
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  INTEGER(PREC_MATIDX) :: lu,jlu,ilup
  INTEGER, DIMENSION(:), POINTER :: p_Iwork,p_Iwork2
  INTEGER :: h_Iwork,h_Iwork2
  INTEGER, DIMENSION(1) :: IworkTemp
  INTEGER :: ifill
  REAL(DP) :: drelax
  
  ! Declare our ILUS-routine from SPLIB as interface to be sure, parameter
  ! interfaces are checked by the compiler.
  INTERFACE
    SUBROUTINE ilus(n,a,colind,rwptr,&
                  s,relax,&
                  lu,jlu,ilup,&
                  iwork,maxstr,&
                  ierr,mneed)
      INTEGER   n, iwork(*), s,  ierr, rwptr(*), colind(*)
      DOUBLE PRECISION a(*), relax
      INTEGER  mneed, maxstr, nzlu, remain
      LOGICAL milu
      INTEGER lu, jlu, ilup
    END SUBROUTINE
  END INTERFACE
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there's no matrix assigned
    
    IF (rsolverNode%rsystemMatrix%NEQ .EQ. 0) THEN
      PRINT *,'Error: No matrix associated!'
      STOP
    END IF
    
    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
    
    p_rmatrix => rsolverNode%rsystemMatrix

    ! We only support scalar 1x1 matrices in structure 7 and 9.
    IF (p_rmatrix%ndiagBlocks .NE. 1) THEN
      PRINT *,'(M)ILU(s) supports only 1x1 matrices!'
      STOP
    END IF

    p_rmatrixSc => p_rmatrix%RmatrixBlock(1,1)

    ! We only support scalar 1x1 matrices in structure 7 and 9.
    IF (p_rmatrixSc%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'(M)ILU(s) supports only double precision matrices!'
      STOP
    END IF
    
    IF ((p_rmatrixSc%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND.&
        (p_rmatrixSc%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
      PRINT *,'(M)ILU(s) supports only structure 7 and 9 matrices!'
      STOP
    END IF
    
    ! Get the matrix description
    CALL storage_getbase_int (p_rmatrixSc%h_Kld,p_Kld)
    CALL storage_getbase_int (p_rmatrixSc%h_Kcol,p_Kcol)
    CALL storage_getbase_double (p_rmatrixSc%h_DA,p_DA)

    ! Calculate a memory guess for how much memory the matrix needs.
    NULLIFY(p_Iwork)
    maxstr = 0

    ! Now calculate the decomposition. Probably allocate more memory if it's
    ! not enough.
    ! Calculate the (M)ILUs matrix using SPLIB.
    ! Get the parameters from the solver node; they were set during the
    ! initialisatin of the solver.
    ifill = rsolverNode%p_rsubnodeMILUs1x1%ifill
    drelax = rsolverNode%p_rsubnodeMILUs1x1%drelax
    
    ! Pass IworkTemp in the first call, not p_Iwork - as p_Iwork points
    ! to NULL. This would be a pill that makes some compilers die ^^.
    CALL ilus(p_rmatrixSc%NEQ,p_DA,p_Kcol,p_Kld,&
              ifill,drelax,&
              lu,jlu,ilup,&
              IworkTemp,maxstr,&
              ierr,mneed)
              
    maxstr = MAX(mneed,3*p_rmatrixSc%NA+(3+4*ifill)*p_rmatrixSc%NEQ)+10000
    !maxstr = MAX(mneed,3*p_rmatrixSc%NA+3*p_rmatrixSc%NEQ)+10000
    DO
      ! Allocate the memory
      CALL storage_new1D ('linsol_initDataMILUs1x1', 'Iwork', maxstr, &
                          ST_INT, h_Iwork, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int(h_Iwork,p_Iwork)
    
      ! Calculate the (M)ILUs matrix using SPLIB.
      CALL ilus(p_rmatrixSc%NEQ,p_DA,p_Kcol,p_Kld,&
                ifill,drelax,&
                lu,jlu,ilup,&
                p_Iwork,maxstr,&
                ierr,mneed)

      ! Error?
      SELECT CASE (ierr)
      CASE (:-1)
        PRINT *,'Warning: (M)ILU(s) decomposition singular!'
        EXIT
      CASE (0)
        EXIT   ! all ok
      CASE (1)
        ! Reallocate memory, we need more
        maxstr = maxstr + MAX(mneed,p_rmatrixSc%NEQ)
      CASE DEFAULT
        ierror = LINSOL_ERR_INITERROR
        RETURN
      END SELECT
      
      ! Try again...
      CALL storage_free(h_Iwork)
    
    END DO
    
    IF (h_Iwork .NE. ST_NOHANDLE) THEN
      ! If less than the half of the memory ILU wanted to have is used,
      ! we reallocate the memory. It does not make sense to have that much
      ! waste!
      IF (mneed .LT. MAXSTR/2) THEN
        CALL storage_new1D ('linsol_initDataMILUs1x1', 'Iwork', mneed, &
                            ST_INT, h_Iwork2, ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int(h_Iwork2,p_Iwork2)
        CALL lalg_copyVectorInt (p_Iwork(1:SIZE(p_Iwork2)),p_Iwork2)
        CALL storage_free (h_Iwork)
        h_Iwork = h_Iwork2
      END IF
      ! Save the handle and the matrix parameters to the MILU structure
      rsolverNode%p_rsubnodeMILUs1x1%h_Idata = h_Iwork
      rsolverNode%p_rsubnodeMILUs1x1%lu = lu
      rsolverNode%p_rsubnodeMILUs1x1%jlu = jlu
      rsolverNode%p_rsubnodeMILUs1x1%ilup = ilup
    END IF
      
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataMILUs1x1 (rsolverNode,isolverSubgroup)
  
!<description>
  ! Releases the memory of teh numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
    
    ! Release the ILU matrix.
    IF (rsolverNode%p_rsubnodeMILUs1x1%h_Idata .NE. ST_NOHANDLE) THEN
      CALL storage_free (rsolverNode%p_rsubnodeMILUs1x1%h_Idata)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneMILUs1x1 (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  INTEGER :: isubgroup
  
    ! The done-routine always releases everything.
    isubgroup = rsolverNode%isolverSubgroup
    
    ! Release symbolical and numerical factorisation if still associated...
    CALL linsol_doneDataMILUs1x1 (rsolverNode, isubgroup)
    
    DEALLOCATE(rsolverNode%p_rsubnodeMILUs1x1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precMILUS1x1 (rsolverNode,rd)
  
!<description>
  ! Applies (M)ILU(s) preconditioner $LU \approx A$ to the defect 
  ! vector rd and solves $LUd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the (M)ILU(s) solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Dd
    INTEGER(PREC_MATIDX) :: lu,jlu,ilup
    INTEGER, DIMENSION(:), POINTER :: p_Iwork
    INTEGER :: h_Iwork

    ! Declare SPLIB-routine as interface to make sure, procedure interfaces
    ! are checked by the compiler
    INTERFACE  
      SUBROUTINE lusolt (n, x, lu, jlu, uptr)
        INTEGER jlu(*),uptr(*),n
        DOUBLE PRECISION x(n)
        ! Note that we changed the interface here in contrast to the original
        ! LUSOLT routine - to make it possible to pass an integer array
        ! as double precision array. Bad practise, but SPLIB is set up 
        ! this way :(
        INTEGER lu(*)
      END SUBROUTINE
    END INTERFACE
    
    IF (rd%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'(M)ILU(s) only supports double precision vectors!'
      STOP
    END IF
    
    ! Get the data array of rd
    CALL lsysbl_getbase_double (rd,p_Dd)
    
    ! Get MILUs information from the parameter block
    h_Iwork = rsolverNode%p_rsubnodeMILUs1x1%h_Idata 
    lu = rsolverNode%p_rsubnodeMILUs1x1%lu 
    jlu = rsolverNode%p_rsubnodeMILUs1x1%jlu 
    ilup = rsolverNode%p_rsubnodeMILUs1x1%ilup 
    CALL storage_getbase_int (h_Iwork,p_Iwork)

    ! Solve the system. Call SPLIB, this overwrites the defect vector
    ! with the preconditioned one.
    CALL lusolt (SIZE(p_Dd),p_Dd, p_Iwork(lu:), p_Iwork(jlu:), p_Iwork(ilup:))
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the BiCGStab solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_Rfilter)
  
!<description>
  ! Creates a t_linsolNode solver structure for the BiCGStab solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! OPTIONAL: A pointer to the solver structure of a solver that should be 
  ! used for preconditioning. If not given or set to NULL(), no preconditioning 
  ! will be used.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rpreconditioner
  
  ! OPTIONAL: A pointer to a filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the 
  ! iteration. If not given or set to NULL(), no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  TYPE(t_filterChain), DIMENSION(:), POINTER, OPTIONAL   :: p_Rfilter
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    CALL linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_BICGSTAB
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK    + &
                                LINSOL_ABIL_CHECKDEF + &
                                LINSOL_ABIL_USESUBSOLVER + &
                                LINSOL_ABIL_USEFILTER
    
    ! Allocate the subnode for BiCGStab.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeBiCGStab)
    
    ! Attach the preconditioner if given. 
    
    IF (PRESENT(p_rpreconditioner)) THEN 
      p_rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner => p_rpreconditioner
    END IF

    ! Attach the filter if given. 
    
    IF (PRESENT(p_Rfilter)) THEN
      p_rsolverNode%p_rsubnodeBiCGStab%p_RfilterChain => p_Rfilter
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixBiCGStab (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of BiCGStab
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_setMatrices (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                              Rmatrices)
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureBiCGStab (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup,i
    TYPE(t_linsolSubnodeBiCGStab), POINTER :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_initStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                                isubgroup)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be initialised
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

    ! BiCGStab needs 5 temporary vectors + 1 for preconditioning. 
    ! Allocate that here! Use the default data type prescribed in the solver 
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    DO i=1,6
      CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%RtempVectors(i),.FALSE.,rsolverNode%cdefaultDataType)
    END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataBiCGStab (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_initData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                            isubgroup,ierror)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_doneData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                            isubgroup)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup,i
    TYPE(t_linsolSubnodeBiCGStab), POINTER :: p_rsubnode
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_doneStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                                isubgroup)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be released
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    IF (p_rsubnode%RtempVectors(1)%NEQ .NE. 0) THEN
      DO i=6,1,-1
        CALL lsysbl_releaseVector (p_rsubnode%RtempVectors(i))
      END DO
    END IF
      
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneBiCGStab (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the BiCGStab solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it's also released from the heap by calling 
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the BiCGStab solver.
  TYPE(t_linsolNode), POINTER         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there's a preconditioner attached. If yes, release it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_releaseSolver(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    END IF
    
    ! Release memory if still associated
    CALL linsol_doneDataBiCGStab (rsolverNode, rsolverNode%isolverSubgroup)
    CALL linsol_doneStructureBiCGStab (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the BiCGStab subnode
    DEALLOCATE(rsolverNode%p_rsubnodeBiCGStab)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precBiCGStab (rsolverNode,rd)
  
!<description>
  ! Applies BiCGStab preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
  !
  ! The implementation follows the original paper introducing BiCGStab:
  !   van der Vorst, H.A.; BiCGStab: A Fast and Smoothly Converging
  !   Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems;
  !   SIAM J. Sci. Stat. Comput. 1992, Vol. 13, No. 2, pp. 631-644
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the BiCGStab solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: dalpha,dbeta,domega0,domega1,domega2,dres
  REAL(DP) :: drho1,drho0,dfr
  INTEGER :: ireslength,ite,i

  ! The queue saves the current residual and the two previous residuals.
  REAL(DP), DIMENSION(32) :: Dresqueue
  
  ! The system matrix
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  INTEGER :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  LOGICAL bprec,bfilter
  
  ! Our structure
  TYPE(t_linsolSubnodeBiCGStab), POINTER :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  TYPE(t_vectorBlock), POINTER :: p_DR,p_DR0,p_DP,p_DPA,p_DSA,p_rx
  TYPE(t_linsolNode), POINTER :: p_rprecSubnode
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    IF ((rd%NEQ .EQ. 0) .OR. (p_rmatrix%NEQ .EQ. 0) .OR. &
        (p_rmatrix%NEQ .NE. rd%NEQ) ) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate

    ireslength = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    bfilter = ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_RfilterChain)
    
    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! Set pointers to the temporary vectors
    p_DR   => p_rsubnode%RtempVectors(1)
    p_DR0  => p_rsubnode%RtempVectors(2)
    p_DP   => p_rsubnode%RtempVectors(3)
    p_DPA  => p_rsubnode%RtempVectors(4)
    p_DSA  => p_rsubnode%RtempVectors(5)
    p_rx   => p_rsubnode%RtempVectors(6)
    
    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    CALL lsysbl_assignDiscretIndirect (rd,p_DR )
    CALL lsysbl_assignDiscretIndirect (rd,p_DR0)
    CALL lsysbl_assignDiscretIndirect (rd,p_DP )
    CALL lsysbl_assignDiscretIndirect (rd,p_DPA)
    CALL lsysbl_assignDiscretIndirect (rd,p_DSA)
    CALL lsysbl_assignDiscretIndirect (rd,p_rx)
    
    IF (bprec) THEN
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    END IF
    IF (bfilter) THEN
      p_RfilterChain => p_rsubnode%p_RfilterChain
    END IF
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    CALL lsysbl_clearVector (p_rx)
      
    ! Initialize used vectors with zero
      
    CALL lsysbl_clearVector(p_DP)
    CALL lsysbl_clearVector(p_DPA)
    
    ! Initialise the iteration vector with zero.

    ! Initialization

    drho0  = 1.0_DP
    dalpha = 1.0_DP
    domega0 = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.

    CALL lsysbl_copyVector(rd,p_DR)
    IF (bfilter) THEN
      ! Apply the filter chain to the vector
      CALL filter_applyFilterChainVec (p_DR, p_RfilterChain)
    END IF
    IF (bprec) THEN
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      CALL linsol_precondDefect (p_rprecSubnode,p_DR)
    END IF
    
    ! Get the norm of the residuum
    dres = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
    IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
              (dres .LE. 1E99_DP))) dres = 0.0_DP

    ! Initialize starting residuum
      
    rsolverNode%dinitialDefect = dres

    ! initialize the queue of the last residuals with RES

    Dresqueue = dres

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    IF ( rsolverNode%dinitialDefect .LT. rsolverNode%drhsZero ) THEN
     
      ! final defect is 0, as initialised in the output variable above

      CALL lsysbl_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    ELSE

      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        PRINT *,&
          'BiCGStab: Iteration ',0,',  !!RES!! = ',rsolverNode%dinitialDefect
      END IF

      CALL lsysbl_copyVector(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        drho1 = lsysbl_scalarProduct (p_DR0,p_DR) 

        IF (drho0*domega0 .EQ. 0.0_DP) THEN
          ! Should not happen
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            PRINT *,&
      'BiCGStab: Iteration prematurely stopped! Correction vector is zero!'
          END IF

          ! Some tuning for the output, then cancel.

          rsolverNode%iresult = -1
          rsolverNode%iiterations = ITE-1
          EXIT
          
        END IF

        dbeta=(drho1*dalpha)/(drho0*domega0)
        drho0 = drho1

        CALL lsysbl_vectorLinearComb (p_DR ,p_DP,1.0_DP,dbeta)
        CALL lsysbl_vectorLinearComb (p_DPA ,p_DP,-dbeta*domega0,1.0_DP)

        CALL lsysbl_blockMatVec (p_rmatrix, p_DP,p_DPA, 1.0_DP,0.0_DP)
        
        IF (bfilter) THEN
          ! Apply the filter chain to the vector
          CALL filter_applyFilterChainVec (p_DPA, p_RfilterChain)
        END IF
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL linsol_precondDefect (p_rprecSubnode,p_DPA)
        END IF

        dalpha = lsysbl_scalarProduct (p_DR0,p_DPA)
        
        IF (dalpha .EQ. 0.0_DP) THEN
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            PRINT *,'BiCGStab: Convergence failed, ALPHA=0!'
            rsolverNode%iresult = -2
            EXIT
          END IF
        END IF
        
        dalpha = drho1/dalpha

        CALL lsysbl_vectorLinearComb (p_DPA,p_DR,-dalpha,1.0_DP)

        CALL lsysbl_blockMatVec (p_rmatrix, p_DR,p_DSA, 1.0_DP,0.0_DP)
        
        IF (bfilter) THEN
          ! Apply the filter chain to the vector
          CALL filter_applyFilterChainVec (p_DSA, p_RfilterChain)
        END IF
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL linsol_precondDefect (p_rprecSubnode,p_DSA)
        END IF
        
        domega1 = lsysbl_scalarProduct (p_DSA,p_DR)
        domega2 = lsysbl_scalarProduct (p_DSA,p_DSA)
        
        IF (domega1 .EQ. 0.0_DP) THEN
          domega0 = 0.0_DP
        ELSE
          IF (domega2 .EQ. 0.0_DP) THEN
            IF (rsolverNode%ioutputLevel .GE. 2) THEN
              PRINT *,'BiCGStab: Convergence failed: omega=0!'
              rsolverNode%iresult = -2
              EXIT
            END IF
          END IF
          domega0 = domega1/domega2
        END IF

        CALL lsysbl_vectorLinearComb (p_DP ,p_rx,dalpha,1.0_DP)
        CALL lsysbl_vectorLinearComb (p_DR ,p_rx,domega0,1.0_DP)

        CALL lsysbl_vectorLinearComb (p_DSA,p_DR,-domega0,1.0_DP)

        ! Get the norm of the new (final?) residuum
        dfr = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        ! Shift the queue with the last residuals and add the new
        ! residual to it. Check length of ireslength to be larger than
        ! 0 as some compilers might produce Floating exceptions
        ! otherwise! (stupid pgf95)
        IF (ireslength .GT. 0) &
          dresqueue(1:ireslength) = EOSHIFT(dresqueue(1:ireslength),1,dfr)

        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        IF (linsol_testDivergence(rsolverNode,dfr)) THEN
          PRINT *,'BiCGStab: Solution diverging!'
          rsolverNode%iresult = 1
          EXIT
        END IF
     
        ! At least perform nminIterations iterations
        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          IF (linsol_testConvergence(rsolverNode,dfr)) EXIT
          
        END IF

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          !WRITE (MTERM,'(A,I7,A,D25.16)') 
          PRINT *,'BiCGStab: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
        END IF

      END DO

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) &
        ite = rsolverNode%nmaxIterations

      ! Finish - either with an error or if converged.
      ! Print the last residuum.


      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        !WRITE (MTERM,'(A,I7,A,D25.16)') 
        PRINT *,'BiCGStab: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
      END IF

    END IF

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    CALL lsysbl_copyVector (p_rx,rd)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
    
      ! Calculate asymptotic convergence rate
    
      IF (dresqueue(1) .GE. 1E-70_DP) THEN
        I = MAX(1,MIN(rsolverNode%iiterations,ireslength-1))
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(I,DP))
      END IF

      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :

      IF (rsolverNode%dfinalDefect .GT. rsolverNode%drhsZero) THEN
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/REAL(rsolverNode%iiterations,DP))
      END IF
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A)') 
        PRINT *,'BiCGStab statistics:'
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A,I5)')     'Iterations              : ',
    !*            IPARAM(OITE)
        PRINT *,'Iterations              : ',rsolverNode%iiterations
        !WRITE (MTERM,'(A,D24.12)') '!!INITIAL RES!!         : ',
    !*            DPARAM(ODEFINI)
        PRINT *,'!!INITIAL RES!!         : ',rsolverNode%dinitialDefect
        !WRITE (MTERM,'(A,D24.12)') '!!RES!!                 : ',
    !*            DPARAM(ODEFFIN)
        PRINT *,'!!RES!!                 : ',rsolverNode%dfinalDefect
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
    !*              rparam%dfinalDefect / rparam%dinitialDefect
          PRINT *,'!!RES!!/!!INITIAL RES!! : ',&
                  rsolverNode%dfinalDefect / rsolverNode%dinitialDefect
        ELSE
          !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
    !*              0D0
          PRINT*,'!!RES!!/!!INITIAL RES!! : ',0.0_DP
        END IF
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A,D24.12)') 'Rate of convergence     : ',
    !*            DPARAM(ORHO)
        PRINT *,'Rate of convergence     : ',rsolverNode%dconvergenceRate

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'BiCGStab: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'BiCGStab: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
        PRINT *,&
              'BiCGStab: Iterations/Rate of convergence: ',&
              rsolverNode%iiterations,' /',rsolverNode%dconvergenceRate
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the Multigrid solver
! *****************************************************************************

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_addMultigridLevel (p_rlevelInfo,rsolverNode, &
                    rinterlevelProjection, &
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver, &
                    iappend)
                    
!<description>
  ! This routine adds a new level to the linked list of levels in the multigrid
  ! solver identified by rsolverNode. 
  ! The given coarse-grid solver and smoother-structures are attached
  ! to the level. The routine returns a pointer to the new level-info structure
  ! in p_rlevelInfo to allow the caller to modify the standard settings of
  ! that level if necessary.
  !
  ! It's allowed to use p_rpresmoother=p_rpostsmoother. A value of
  ! NULL() is also permissable, this deactivates the corresponding smoother.
  !
  ! It's allowed to use the same rinterlevelProjection on all levels,
  ! since the structure is level independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway).
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver, where the level
  ! should be added to. Must already be initialised for the multigrid
  ! solver.
  TYPE(t_linsolNode)                             :: rsolverNode
!</inputoutput>
  
!<input>
  ! An interlevel projection structure that configures the projection
  ! between the solutions on a finer and a coarser grid. The structure
  ! must have been initialised with mlprj_initProjection.
  !
  ! Note that this structure is level-independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway), so the same structure can be used
  ! to initialise all levels!
  TYPE(t_interlevelProjectionBlock) :: rinterlevelProjection
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for presmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rpresmoother
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for postsmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rpostsmoother

  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for coarse grid solving. 
  ! Should only be given for the very first level.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rcoarseGridSolver
  
  ! Optional: Position where to put the new structure.
  ! YES or not specified: append the structure to the end of the list as new
  ! higher level.
  ! NO: Insert the structure as new lowest level. The caller has to make sure
  ! that the coarse grid solver on the previous lowest level is removed!
  ! (i.e. removed from the p_rlevelInfo node corresponding to the lowest
  ! level!)
  INTEGER, INTENT(IN), OPTIONAL                  :: iappend
!</input>
  
!<output>  
  ! The t_levelInfo structure for the new level that was added to the
  ! multigrid solver.
  TYPE(t_linsolMGLevelInfo), POINTER     :: p_rlevelInfo
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF
  
  ! Create a new level-info structure
  NULLIFY(p_rlevelInfo)
  ALLOCATE(p_rlevelInfo)
  
  ! Attach the sub-solvers
  IF (ASSOCIATED(p_rcoarseGridSolver)) THEN
    p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
  END IF

  IF (ASSOCIATED(p_rpresmoother)) THEN
    p_rlevelInfo%p_rpresmoother => p_rpresmoother
  END IF

  IF (ASSOCIATED(p_rpostsmoother)) THEN
    p_rlevelInfo%p_rpostsmoother => p_rpostsmoother
  END IF
  
  ! Initialise the interlevel projection structure,
  ! copy the data of our given template.
  p_rlevelInfo%rinterlevelProjection = rinterlevelProjection

  ! Attach the level-info structure to the linked list.
  ! Does the list exist?
  IF (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)) THEN
  
    ! List does not exist. Create a new one.
    rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
    rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
    rsolverNode%p_rsubnodeMultigrid%nlevels = 1
  
  ELSE
  
    ! Do we have to insert it as new head?
    i = YES
    IF (PRESENT(iappend)) i=iappend
    
    IF (i .EQ. NO) THEN
    
      ! New lowest level
      p_rlevelInfo%p_rnextLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead%p_rprevLevel => p_rlevelInfo
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
    
    ELSE
    
      ! New highest level
      p_rlevelInfo%p_rprevLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail%p_rnextLevel => p_rlevelInfo
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
      
    END IF
    
    ! Increase the number of existing levels
    rsolverNode%p_rsubnodeMultigrid%nlevels = &
      rsolverNode%p_rsubnodeMultigrid%nlevels + 1
  
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_removeMultigridLevel (rsolverNode,bcoarseLevel)
                    
!<description>
  ! This routine removes the highest or lowest multigrid level from the 
  ! multigrid structure and from the heap. The associated solver
  ! structures of subsolvers (coarse grid solver, smoother) are
  ! released from the heap.
  !
  ! If the coarse grid is removed, the caller has to make sure that
  ! the new coarse grid receives a coarse grid solver!
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver.
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!<input>
  ! OPTIONAL: Defines which level to delete.
  ! .TRUE. = Delete the coarse level.
  ! .FALSE. or not defined = Delete the finest level.
  LOGICAL, INTENT(IN), OPTIONAL :: bcoarseLevel
!</input>  
  
!</subroutine>

  ! local variables
  LOGICAL :: bcoarse
  
  ! The t_levelInfo structure to delete
  TYPE(t_linsolMGLevelInfo), POINTER     :: p_rlevelInfo
  
  bcoarse = .FALSE.
  IF (PRESENT(bcoarseLevel)) bcoarse = bcoarseLevel
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF
  
  ! Where's the levelInfo-structure to delete?
  IF (.NOT. bcoarse) THEN
    ! Remove the fine level node - if it exists.
    p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail

    IF (ASSOCIATED(p_rlevelInfo)) THEN
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo%p_rprevLevel
      
      IF (rsolverNode%p_rsubnodeMultigrid%nlevels .EQ. 1) THEN
        ! We remove the last node.
        NULLIFY(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)
      END IF
      
      IF (ASSOCIATED(p_rlevelInfo%p_rprevLevel)) THEN
        NULLIFY(p_rlevelInfo%p_rprevLevel%p_rnextLevel)
      END IF
    END IF
  ELSE
    ! Remove the coarse level node - if it exists.
    p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    IF (ASSOCIATED(p_rlevelInfo)) THEN
      
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo%p_rnextLevel
      
      IF (rsolverNode%p_rsubnodeMultigrid%nlevels .EQ. 1) THEN
        ! We remove the last node.
        NULLIFY(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail)
      END IF

      IF (ASSOCIATED(p_rlevelInfo%p_rnextLevel)) THEN
        NULLIFY(p_rlevelInfo%p_rnextLevel%p_rprevLevel)
      END IF
    END IF
  END IF
  
  IF (ASSOCIATED(p_rlevelInfo)) THEN
    ! Decrease the number of existing levels
    rsolverNode%p_rsubnodeMultigrid%nlevels = &
      rsolverNode%p_rsubnodeMultigrid%nlevels - 1

    ! Release the temporary vectors of this level. But be careful!
    !
    ! If the node is a coarse level node...
    IF (.NOT. ASSOCIATED(p_rlevelInfo%p_rprevLevel)) THEN
      ! ...don't release RHS/temp vector, as there is none. Instead,
      ! release RHS/solution vector of the next higher level as that
      ! this one gets the new coarse grid.
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)) THEN
        ! The vector may not exist - if the application has not called
        ! initStructure!
        IF (rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead%rrhsVector%NEQ .NE. 0) &
          CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeMultigrid% &
                                    p_rlevelInfoHead%rrhsVector)
        IF (rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead%rtempVector%NEQ .NE. 0) &
          CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeMultigrid% &
                                    p_rlevelInfoHead%rtempVector)
      END IF
      
      ! But there may be a solution vector in p_rlevelInfo that
      ! is to be deleted - if the level is not also the coarse grid
      IF (rsolverNode%p_rsubnodeMultigrid%nlevels .NE. 0) THEN
        IF (p_rlevelInfo%rsolutionVector%NEQ .NE. 0) &
          CALL lsysbl_releaseVector(p_rlevelInfo%rsolutionVector)
      END IF

    END IF

    ! If the node is a fine level node...
    IF (.NOT. ASSOCIATED(p_rlevelInfo%p_rnextLevel)) THEN
      ! ...don't release solution vector, as there is none (it shares
      ! its memory with the vector rd to be preconditioned). Instead,
      ! release solution vector of the next lower level as that
      ! this one gets the new fine grid.
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail)) THEN
        IF (rsolverNode%p_rsubnodeMultigrid% &
            p_rlevelInfoTail%rsolutionVector%NEQ .NE. 0) &
          CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeMultigrid% &
                                    p_rlevelInfoTail%rsolutionVector)
      END IF
      
      ! But there may be a RHS/temporary vector in p_rlevelInfo that
      ! is to be deleted - if the level is not also the coarse grid
      IF (rsolverNode%p_rsubnodeMultigrid%nlevels .NE. 0) THEN
        IF (p_rlevelInfo%rrhsVector%NEQ .NE. 0) &
          CALL lsysbl_releaseVector(p_rlevelInfo%rrhsVector)
        IF (p_rlevelInfo%rtempVector%NEQ .NE. 0) &
          CALL lsysbl_releaseVector(p_rlevelInfo%rtempVector)
      END IF
    END IF
  
    ! Release sub-solvers if associated.
    !
    ! Caution: Pre- and postsmoother may be identical!
    IF (ASSOCIATED(p_rlevelInfo%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rlevelInfo%p_rpreSmoother, &
                        p_rlevelInfo%p_rpostSmoother))) THEN
      CALL linsol_releaseSolver(p_rlevelInfo%p_rpostsmoother)
    END IF

    IF (ASSOCIATED(p_rlevelInfo%p_rpresmoother)) THEN
      CALL linsol_releaseSolver(p_rlevelInfo%p_rpresmoother)
    END IF

    IF (ASSOCIATED(p_rlevelInfo%p_rcoarseGridSolver)) THEN
      CALL linsol_releaseSolver(p_rlevelInfo%p_rcoarseGridSolver)
    END IF

    ! Remove the structure
    DEALLOCATE(p_rlevelInfo)
    
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_cleanMultigridLevels (rsolverNode)
                    
!<description>
  ! This routine removes all level information from the MG solver and releases
  ! all attached solver nodes on all levels (smoothers, coarse grid solvers,
  ! ...). 
  ! It can be used to clean up all levels before rebuilding the level
  ! structure by addMultigridLevel.
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver.
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!</subroutine>

  ! Remove all the levels
  DO WHILE (rsolverNode%p_rsubnodeMultigrid%nlevels .GT. 0) 
    CALL linsol_removeMultigridLevel (rsolverNode)
  END DO

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_getMultigridLevel (rsolverNode,ilevel,p_rlevelInfo)
                    
!<description>
  ! Searches inside of the multigrid structure for level ilevel and returns
  ! a pointer to the corresponding p_rlevelInfo structure
  ! (or NULL() if the level does not exist).
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!<input>
  ! Number of the level to fetch.
  INTEGER, INTENT(IN) :: ilevel
!</input>  
  
!<output>
  ! A pointer to the corresponding t_levelInfo structure or NULL()
  ! if the level does not exist.
  TYPE(t_linsolMGLevelInfo), POINTER     :: p_rlevelInfo
!</output>  
  
!</subroutine>

  ! local variables
  INTEGER :: i

  ! Do we have the level?
  IF ((ilevel .LT. 0) .OR. &
      (ilevel .GT. rsolverNode%p_rsubnodeMultigrid%nlevels)) THEN
    NULLIFY(p_rlevelInfo)
    RETURN
  END IF
  
  ! Go through the linear list to search for the level
  p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO i=1,ilevel-1
    p_rlevelInfo => p_rlevelInfo%p_rnextLevel
  END DO

  END SUBROUTINE
  
  ! ***************************************************************************

!<function>
  
  INTEGER FUNCTION linsol_getMultigridLevelCount (rsolverNode)
                    
!<description>
  ! Returns the number of levels currently attached to the multigrid solver.
!</description>
  
!<input>
  ! The solver structure of the multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
!</input>

!<return>
  ! The number of levels in the MG solver.
!</return>
  
!</function>

    linsol_getMultigridLevelCount = rsolverNode%p_rsubnodeMultigrid%nlevels

  END FUNCTION
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initMultigrid (p_rsolverNode,p_Rfilter)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the Multigrid solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  ! Before the solver can used for solving the problem, the caller must add
  ! information about all levels (matrices,...) to the solver. 
  ! This can be done by linsol_addMultigridLevel.
  
!</description>
  
!<input>
  
  ! Optional: A pointer to a filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the 
  ! iteration. If not given or set to NULL(), no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  TYPE(t_filterChain), DIMENSION(:), POINTER, OPTIONAL   :: p_Rfilter
  
!</input>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER             :: p_rsolverNode
   
!</output>
  
!</subroutine>

  ! Create a default solver structure
  CALL linsol_initSolverGeneral(p_rsolverNode)
  
  ! Initialise the type of the solver
  p_rsolverNode%calgorithm = LINSOL_ALG_MULTIGRID
  
  ! Initialise the ability bitfield with the ability of this solver:
  p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR     + LINSOL_ABIL_BLOCK        + &
                              LINSOL_ABIL_MULTILEVEL + LINSOL_ABIL_CHECKDEF     + &
                              LINSOL_ABIL_USESUBSOLVER + &
                              LINSOL_ABIL_USEFILTER
  
  ! Allocate the subnode for Multigrid.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(p_rsolverNode%p_rsubnodeMultigrid)
  
  ! Attach the filter if given. 
  
  IF (PRESENT(p_Rfilter)) THEN
    p_rsolverNode%p_rsubnodeMultigrid%p_RfilterChain => p_Rfilter
  END IF
  
  END SUBROUTINE
  
! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixMultigrid (rsolverNode,Rmatrices)
  
!<description>
  
  ! Assigns the system matrix rsystemMatrix to the Multigrid solver on 
  ! all levels.
  
!</description>
  
!<input>
  
  ! An array of system matrices on all levels.
  ! Each level in multigrid is initialised separately.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices

!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
!</inputoutput>
  
!</subroutine>

  ! local variables
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER                       :: ilevel,nlmin
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF
  
  ! Make sure we have the right amount of matrices
  IF (SIZE(Rmatrices) .NE. rsolverNode%p_rsubnodeMultigrid%nlevels) THEN
    PRINT *,'Error: Wrong number of matrices'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. If yes, attach the matrix
  ! to it.
  ! For this purpose, call the general initialisation routine, but
  ! pass only that part of the Rmatrices array that belongs
  ! to the range of levels between the coarse grid and the current grid.
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  nlmin  = LBOUND(Rmatrices,1)
  ilevel = LBOUND(Rmatrices,1)
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rpreSmoother, &
                                Rmatrices(nlmin:ilevel))
    END IF
    ! Pre- and postsmoother may be identical!
    ! Take care not to initialise the same smoother twice!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rpostSmoother, &
                                Rmatrices(nlmin:ilevel))
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rcoarseGridSolver, &
                                Rmatrices(nlmin:ilevel))
    END IF
    p_rcurrentLevel%rsystemMatrix = Rmatrices(ilevel)
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
    ilevel = ilevel + 1
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureMultigrid (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  INTEGER(PREC_VECIDX) :: imemmax
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_vectorBlock), POINTER :: p_rtemplVect
  
  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.

  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    ! Pre- and postsmoother may be identical!
    ! Take care not to initialise the same smoother twice!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  ! Cancel here, if we don't belong to the subgroup to be initialised
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

  imemmax = 0

  ! Multigrid needs temporary vectors for the RHS, solution and general
  ! data. Memory for that is allocated on all levels for temporary, RHS and 
  ! solution vectors.
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    p_rmatrix => p_rcurrentLevel%rsystemMatrix
    NULLIFY(p_rtemplVect)
    IF (ASSOCIATED(p_rcurrentLevel%p_rnextLevel)) THEN
      ! On the maximum level we don't need additional memory for the 
      ! solution vector - as solution vector on the fine grid, 
      ! the vector which comes as parameter
      ! to the multigrid  preconditioner is used.
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rsolutionVector,.FALSE.,&
            rsolverNode%cdefaultDataType)
      p_rtemplVect => p_rcurrentLevel%rsolutionVector
    END IF
    ! On the coarse grid, we don't need memory for temporary/RHS
    ! vectors. The RHS on the coarse grid is replaced in-situ
    ! by the solution vector.
    IF (ASSOCIATED(p_rcurrentLevel%p_rprevLevel)) THEN
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rrhsVector,.FALSE.,&
            rsolverNode%cdefaultDataType)
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rtempVector,.FALSE.,&
            rsolverNode%cdefaultDataType)
      p_rtemplVect => p_rcurrentLevel%rtempVector
    END IF

    ! Calculate the memory that is probably necessary for resorting
    ! vectors during the iteration.
    IF (ASSOCIATED(p_rtemplVect)) THEN
      imemmax = MAX(imemmax,p_rtemplVect%NEQ)
    END IF
          
    ! Ask the prolongation/restriction routines how much memory is necessary
    ! for prolongation/restriction on that level.
    IF (ASSOCIATED(p_rcurrentLevel%p_rprevLevel)) THEN ! otherwise coarse grid
      ! Calculate the memory that is necessary for prolongation/restriction -
      ! in case there is temporary memory needed.
      ! The system matrix on the fine/coarse grid specifies the discretisation.
      imemmax = MAX(imemmax,mlprj_getTempMemoryMat ( &
                    p_rcurrentLevel%rinterlevelProjection, &
                    p_rcurrentLevel%p_rprevLevel%rsystemMatrix,&
                    p_rcurrentLevel%rsystemMatrix))
    END IF

    ! All temporary vectors are marked as 'unsorted'. We can set the
    ! sorting flag directly here without using the resorting routines
    ! as the vectors have just been created.
    p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy = &
      -ABS(p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy)
    p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy = &
      -ABS(p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy)
    p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy = &
      -ABS(p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy)
    
    ! And the next level...
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  ! Do we have to allocate memory for prolongation/restriction/resorting?
  IF (imemmax .GT. 0) THEN
    CALL lsyssc_createVector (rsolverNode%p_rsubnodeMultigrid%rprjTempVector,&
                              imemmax,.FALSE.,rsolverNode%cdefaultDataType)
  ELSE
    rsolverNode%p_rsubnodeMultigrid%rprjTempVector%NEQ = 0
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataMultigrid (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  
  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup
  
  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rpreSmoother,isubgroup,ierror)
    END IF
    ! Pre- and postsmoother may be identical!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rpostSmoother,isubgroup,ierror)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup,ierror)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataMultigrid (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_doneData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    ! Pre- and postsmoother may be identical!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_doneData(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_doneData(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureMultigrid (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    ! Pre- and postsmoother may be identical!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
  END DO

  ! Cancel here, if we don't belong to the subgroup to be initialised
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

  ! Release temporary memory.
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rprevLevel)) THEN
      IF (p_rcurrentLevel%rtempVector%NEQ .NE. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rtempVector)
      END IF
      
      IF (p_rcurrentLevel%rrhsVector%NEQ .NE. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rrhsVector)
      END IF

    END IF

    IF (ASSOCIATED(p_rcurrentLevel%p_rnextLevel)) THEN
      IF (p_rcurrentLevel%rsolutionVector%NEQ .NE. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rsolutionVector)
      END IF
    END IF

    p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
  END DO

  ! Check if we have to release our temporary vector for prolongation/
  ! restriction.
  ! Do we have to allocate memory for prolongation/restriction?
  IF (rsolverNode%p_rsubnodeMultigrid%rprjTempVector%NEQ .GT. 0) THEN
    CALL lsyssc_releaseVector (rsolverNode%p_rsubnodeMultigrid%rprjTempVector)
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneMultigrid (rsolverNode)
  
!<description>
  
  ! This routine releases all temporary memory for the multigrid solver from
  ! the heap. In particular, this releases the solver structures of all
  ! subsolvers (smoother, coarse grid solver).
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
  
!</description>
  
!<inputoutput>
  ! A pointer to a t_linsolNode structure of the Multigrid solver.
  TYPE(t_linsolNode), INTENT(INOUT)                :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! local variables
  !TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

!  ! Check for every level, if there's a presmoother, postsmoother or
!  ! coarse grid solver structure attached. If yes, release it.
!  
!  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
!  DO WHILE(ASSOCIATED(p_rcurrentLevel))
!    ! Pre- and postsmoother may be identical!
!    ! Take care not to initialise the same smoother twice!
!    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
!        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
!                          p_rcurrentLevel%p_rpostSmoother))) THEN
!      CALL linsol_releaseSolver(p_rcurrentLevel%p_rpostSmoother)
!    END IF
!    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
!      CALL linsol_releaseSolver(p_rcurrentLevel%p_rpreSmoother)
!    END IF
!    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
!      CALL linsol_releaseSolver(p_rcurrentLevel%p_rcoarseGridSolver)
!    END IF
!    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
!  END DO

  ! Release data and structure if this has not happened yet.
  CALL linsol_doneDataMultigrid(rsolverNode)
  CALL linsol_doneStructureMultigrid(rsolverNode)
  
  ! Remove all the levels
  CALL linsol_cleanMultigridLevels (rsolverNode)
  
  ! Release MG substructure, that's it.
  DEALLOCATE(rsolverNode%p_rsubnodeMultigrid)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_smoothCorrection (rsolverNode,&
                       rmatrix,rx,rb,rtemp,p_RfilterChain)
  
!<description>
  ! This routine performs a smoothing process on the vector rx
  ! belonging to a linear system $Ax=b$.
  ! rsolverNode identifies a solver structure that is converted to a
  ! smoother using linsol_convertToSmoother: The number of smoothing
  ! steps is written to rsolverNode%nmaxIterations and 
  ! rsolverNode%nmaxIterations so that the solver performs a definite
  ! number of iterations regardless of the residual.
  !
  ! rx is assumed to be of 'defect' type. There are two types of smoothing
  ! processes, depending on which type of preconditioner rsolverNode is:
  ! 1.) If rsolverNode is an iterative solver, linsol_smoothCorrection
  !     calls the associated solver P to compute $x=P^{-1}b$ using
  !     rsolverNode%nmaxIterations iterations.
  ! 2.) If rsolverNode is a 1-step solver, linsol_smoothCorrection
  !     performs rsolverNode%nmaxIterations defect correction steps
  !      $$ x := x + P^{-1} (b-Ax) $$
  !     with $x_0=0$.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<input>
  ! The RHS vector of the system
  TYPE(t_vectorBlock), INTENT(IN), TARGET            :: rb

  ! The system matrix on the current level.
  TYPE(t_matrixBlock), INTENT(IN)                    :: rmatrix
  
  ! Pointer to the filter chain to use if rsolverNode is a direct
  ! solver. NULL() if no filtering is active.
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
!</input>
  
!<inputoutput>
  ! The solver node containing the solver configuration
  TYPE(t_linsolNode), INTENT(INOUT)                  :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  TYPE(t_vectorBlock), INTENT(INOUT)                 :: rx
  
  ! A temporary vector of the same size and structure as rx.
  TYPE(t_vectorBlock), INTENT(INOUT)                 :: rtemp
!</inputoutput>
  
!</subroutine>

  INTEGER :: i
  LOGICAL :: bfilter
  !DEBUG: REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
  
  ! Cancel if nmaxIterations = number of smoothing steps is =0.
  IF (rsolverNode%nmaxIterations .LE. 0) RETURN

  ! Do we have an iterative or one-step solver given?
  IF (IAND(rsolverNode%ccapability,LINSOL_ABIL_DIRECT) .EQ. 0) THEN
  
    ! Copy rd to rx and call the preconditioner to precondition rx.
    ! This is equivalent to:
    !
    !   $$   x_1   =   P^{-1} b   =   x_0 + P^{-1} (b-Ax_0)   $$
    !
    ! with $x_0 = 0$ and $P^{-1}$ performing nmaxIterations steps
    ! to approximate $A^{-1} b$.
    CALL lsysbl_copyVector(rb,rx)
    CALL linsol_precondDefect(rsolverNode,rx)
    
  ELSE
    
    ! This is a 1-step solver, we have to emulate the smoothing
    ! iterations. Perform rsolverNode%nmaxIterations steps of the
    ! form
    !     $$ x_{n+1} = x_n + P^{-1}(b-Ax_n) $$
    ! with $x_0 = 0$.
    
    bfilter = ASSOCIATED(p_RfilterChain)
    
    !DEBUG: CALL lsysbl_getbase_double (rx,p_Ddata)
    !DEBUG: CALL lsysbl_getbase_double (rtemp,p_Ddata2)
    
    DO i=1,rsolverNode%nmaxIterations
    
      CALL lsysbl_copyVector(rb,rtemp)
      CALL lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
      
      ! Apply the filter to this defect before preconditioning
      IF (bfilter) THEN
        CALL filter_applyFilterChainVec (rtemp, p_RfilterChain)
      END IF
      
      CALL linsol_precondDefect(rsolverNode,rtemp)
      CALL lsysbl_vectorLinearComb (rtemp,rx,1.0_DP,1.0_DP)
      
    END DO
    
  END IF
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precMultigrid (rsolverNode,rd)
  
!<description>
  ! Applies Multigrid preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The structure of the levels (pre/postsmoothers, interlevel projection
  ! structures) must have been prepared with linsol_addMultigridLevel
  ! before calling this routine.
  ! The matrices $A$ on all levels must be attached to the solver previously 
  ! by linsol_setMatrices.
  
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: nminIterations,nmaxIterations,ireslength,niteResOutput
  INTEGER :: ite,ilev,nlmax,i,nblocks
  REAL(DP) :: dres,dstep
  LOGICAL :: bfilter,bsort
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
  ! The queue saves the current residual and the two previous residuals.
  REAL(DP), DIMENSION(32) :: Dresqueue
  
  ! The system matrix on the current level
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  
  ! Our MG structure
  TYPE(t_linsolSubnodeMultigrid), POINTER :: p_rsubnode
  
  ! The current level and the next lower one.
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel,p_rlowerLevel
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeMultigrid
    
    bfilter = ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_RfilterChain)
    p_RfilterChain => rsolverNode%p_rsubnodeMultigrid%p_RfilterChain
    
    ! Get the system matrix on the finest level:
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    IF ((rd%NEQ .EQ. 0) .OR. (p_rmatrix%NEQ .EQ. 0) .OR. &
        (p_rmatrix%NEQ .NE. rd%NEQ) ) THEN
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    IF (p_rsubnode%icycle .LT. 0) THEN
      ! Wrong cycle
      PRINT *,'Multigrid: Wrong cycle!'
      rsolverNode%iresult = 2
      RETURN
    END IF
    
    IF ((.NOT. ASSOCIATED(p_rsubnode%p_rlevelInfoHead)) .OR. &
        (.NOT. ASSOCIATED(p_rsubnode%p_rlevelInfoTail)) .OR. &
        (p_rsubnode%nlevels .LE. 0)) THEN
      PRINT *,'Multigrid: No levels attached!'
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    ireslength = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

    ! Minimum/maximum number of iterations
    nminIterations = MAX(rsolverNode%nminIterations,0)
    nmaxIterations = MAX(rsolverNode%nmaxIterations,0)
      
    ! Iteration when the residuum is printed:
    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! Status reset
    rsolverNode%iresult = 0
    rsolverNode%icurrentIteration = 0
    rsolverNode%dinitialDefect = 0.0_DP
    rsolverNode%dfinalDefect = 0.0_DP
    rsolverNode%dconvergenceRate = 0.0_DP
    rsolverNode%dasymptoticConvergenceRate = 0.0_DP
    
    ! Number of blocks in the current equation
    nblocks = rd%nblocks
    
    ! We start on the maximum level. 
    p_rcurrentLevel => p_rsubnode%p_rlevelInfoTail
    ilev = p_rsubnode%nlevels
    nlmax = p_rsubnode%nlevels
    p_rlowerLevel => p_rcurrentLevel%p_rprevLevel
    
    ! Is there only one level? Can be seen if the current level
    ! already contains a coarse grid solver.
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
    
      IF (rsolverNode%ioutputLevel .GT. 1) THEN
        PRINT *,'Multigrid: Only one level. Switching back to standard solver.'
      END IF
      CALL linsol_precondDefect(p_rcurrentLevel%p_rcoarseGridSolver,rd)
      
      ! Take the statistics from the coarse grid solver.
      rsolverNode%dinitialDefect = p_rcurrentLevel%p_rcoarseGridSolver%dinitialDefect
      rsolverNode%dfinalDefect = p_rcurrentLevel%p_rcoarseGridSolver%dfinalDefect
      rsolverNode%dconvergenceRate = p_rcurrentLevel%p_rcoarseGridSolver%dconvergenceRate
      rsolverNode%iiterations = p_rcurrentLevel%p_rcoarseGridSolver%iiterations
      rsolverNode%dasymptoticConvergenceRate = p_rcurrentLevel% &
        p_rcoarseGridSolver%dasymptoticConvergenceRate
      
    ELSE
      ! Get the norm of the initial residuum.
      ! As the initial iteration vector is zero, this is the norm
      ! of the RHS:
      dres = lsysbl_vectorNorm (rd,rsolverNode%iresNorm)
      IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                (dres .LE. 1E99_DP))) dres = 0.0_DP

      ! Initialize starting residuum
      rsolverNode%dinitialDefect = dres

      ! initialize the queue of the last residuals with RES
      Dresqueue = dres

      ! Check if out initial defect is zero. This may happen if the filtering
      ! routine filters "everything out"!
      ! In that case we can directly stop our computation.

      IF ( rsolverNode%dinitialDefect .LT. rsolverNode%drhsZero ) THEN
        ! final defect is 0, as initialised in the output variable above
        CALL lsysbl_clearVector(rd)
        rsolverNode%dfinalDefect = dres
        rsolverNode%dfinalDefect = dres
        rsolverNode%dconvergenceRate = 0.0_DP
        rsolverNode%iiterations = 0
        rsolverNode%dasymptoticConvergenceRate = 0.0_DP
        
      ELSE
        
        ! Ok, we can start with multigrid.
        ! At first check the sorting state of the vector.
        bsort = lsysbl_isVectorSorted(rd)
        
        ! To deal with sorted vectors is a bit tricky! Here a small summary
        ! how sorted vectors is dealed with, in case sorting is activated:
        !
        ! 1.) At the beginning of the algorithm, RHS, solution and temp vector
        !     on the finest level are set up as 'sorted'. The whole algorithm
        !     here deals with sorted vectors, so that matrix/vector multiplication
        !     works.
        ! 2.) At the beginning of the algorithm, all temporary vectors on the
        !     lower levels are set up as 'unsorted'.
        ! 3.) When restricting the defect,
        !     - the defect is unsorted on level l
        !     - the unsorted defect is restricted to level l-1
        !     - the restricted defect on level l-1 is sorted
        ! 4.) When prolongating a correction vector,
        !     - the correction vector on level l-1 is unsorted
        !     - the unsorted vector is prolongated to level l
        !     - the unsorted vector on level l is sorted
        ! So at the end of the MG sweep, we have the same situation as at the start
        ! of the algorithm.
        
        ! At first, reset the cycle counters of all levels
        p_rcurrentLevel => p_rsubnode%p_rlevelInfoHead
        DO WHILE(ASSOCIATED(p_rcurrentLevel%p_rnextLevel))
          IF (p_rsubnode%icycle .EQ. 0) THEN
            p_rcurrentLevel%ncyclesRemaining = 2
          ELSE
            p_rcurrentLevel%ncyclesRemaining = p_rsubnode%icycle
          END IF  
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
        END DO
        
        ! Afterwards, p_rcurrentLevel points to the maximum level again.
        ! There, we set ncyclesRemaining to 1.
        p_rcurrentLevel%ncyclesRemaining = 1
        
        ! Print out the initial residuum

        IF (rsolverNode%ioutputLevel .GE. 2) THEN
          !WRITE (MTERM,'(A,I7,A,D25.16)') 
          PRINT *,'Multigrid: Iteration ',0,',  !!RES!! = ',rsolverNode%dinitialDefect
        END IF
        
        ! Copy the initial RHS to the RHS vector on the maximum level.
        CALL lsysbl_copyVector(rd,p_rcurrentLevel%rrhsVector)
        
        ! Replace the solution vector on the finest level by rd.
        ! Afterwards, rd and the solution vector on the finest level
        ! share the same memory location, so we use rd as iteration
        ! vector directly.
        ! We set bisCopy to true to indicate that this vector is
        ! actually does not belong to us. This is just to be sure that
        ! an accidently performed releaseVector does not release a handle
        ! (although this should not happen).
        p_rcurrentLevel%rsolutionVector = rd
        p_rcurrentLevel%rsolutionVector%bisCopy = .TRUE.
        
        ! Clear the initial solution vector.
        CALL lsysbl_clearVector (p_rcurrentLevel%rsolutionVector)
        
        ! Start multigrid iteration; perform at most nmaxiterations iterations.
        DO ite = 1, nmaxiterations
        
          rsolverNode%icurrentIteration = ite
          
          ! Initialize cycle counters for all levels.
          p_rcurrentLevel => p_rsubnode%p_rlevelInfoHead
          DO WHILE(ASSOCIATED(p_rcurrentLevel%p_rnextLevel))
            p_rcurrentLevel%ncycles = p_rcurrentLevel%ncyclesRemaining
            p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          END DO
          ! Don't forget the highest level.
          p_rcurrentLevel%ncycles = p_rcurrentLevel%ncyclesRemaining
        
          ! p_rcurrentLevel now points to the meximum level.
          ilev = p_rsubnode%nlevels
          p_rlowerLevel => p_rcurrentLevel%p_rprevLevel
          
          ! Build the defect...
          CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rtempVector)
          IF (ite .NE. 1) THEN   ! initial solution vector is zero!
            CALL lsysbl_blockMatVec (&
                 p_rcurrentLevel%rsystemMatrix, &
                 p_rcurrentLevel%rsolutionVector,&
                 p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
          END IF
          IF (bfilter) THEN
            ! Apply the filter chain to the vector
            CALL filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                             p_RfilterChain)
          END IF
          
          cycleloop: DO  ! Loop for the cycles
          
            ! On the maximum level we already built out defect vector. If we are    
            ! on a lower level than NLMAX, perform smoothing+restriction down to the
            ! coarse level. We identify the coarse level by checking if
            ! the current level has a coarse grid solver.
            
            DO WHILE (.NOT. ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver))
            
              ! Perform the pre-smoothing with the current solution vector
              IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
                CALL linsol_smoothCorrection (p_rcurrentLevel%p_rpreSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)
              END IF
            
              ! Build the defect vector
              CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                      p_rcurrentLevel%rtempVector)
              CALL lsysbl_blockMatVec (&
                  p_rcurrentLevel%rsystemMatrix, &
                  p_rcurrentLevel%rsolutionVector,&
                  p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
              IF (bfilter) THEN
                ! Apply the filter chain to the vector
                CALL filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                                p_RfilterChain)
              END IF
              
              ! Restriction of the defect. The restricted defect is placed
              ! in the right hand side vector of the lower level.
              ! The projection parameters from level ilev to level
              ! ilev-1 is configured in the rprojection structure of the
              ! current level.
              !
              ! Make sure the vector is unsorted before the restriction!
              ! We can use the 'global' temporary array p_rprjTempVector
              ! as temporary memory during the resorting process.
              IF (bsort) THEN
                CALL lsysbl_sortVectorInSitu (p_rcurrentLevel%rtempVector,&
                      p_rsubnode%rprjTempVector,.FALSE.)
              END IF
              
              ! When restricting to the coarse grid, we directy restrict into
              ! the solution vector. It's used there as RHS and replaced in-situ
              ! by the solution by the coarse grid solver. So don't need the RHS vector
              ! on the coarse grid and save one vector-copy.
              ! 
              ! Otherwise, we restrict to the RHS on the lower level and continue
              ! the smoothing process there.
              IF (ASSOCIATED(p_rlowerLevel%p_rprevLevel)) THEN
                ! We don't project to the coarse grid
                CALL mlprj_performRestriction (p_rcurrentLevel%rinterlevelProjection,&
                      p_rlowerLevel%rrhsVector, &
                      p_rcurrentLevel%rtempVector, &
                      p_rsubnode%rprjTempVector)

                ! Apply the filter chain (e.g. to implement boundary conditions)
                ! on just calculated right hand side
                ! (which is a defect vector against the 0-vector).
                IF (bfilter) THEN
                  ! Apply the filter chain to the vector.
                  ! We are in 'unsorted' state; applying the filter here is
                  ! supposed to be a bit faster...
                  CALL filter_applyFilterChainVec (p_rlowerLevel%rrhsVector, &
                                                  p_RfilterChain)
                END IF

                IF (bsort) THEN
                  ! Resort the RHS on the lower level.
                  CALL lsysbl_sortVectorInSitu (p_rlowerLevel%rrhsVector,&
                        p_rsubnode%rprjTempVector,.TRUE.)

                  ! Temp-vector and solution vector there are yet uninitialised,
                  ! therefore we can mark them as sorted without calling the
                  ! resorting routine.
                  p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                  ABS(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy) 
                  
                  p_rlowerLevel%rsolutionVector%RvectorBlock(1:nblocks)% &
                    isortStrategy = ABS(p_rlowerLevel%rsolutionVector% &
                    RvectorBlock(1:nblocks)%isortStrategy)
                END IF

                ! Choose zero as initial vector on lower level. 
                CALL lsysbl_clearVector (p_rlowerLevel%rsolutionVector)
                
              ELSE
              
                ! THe vector is to be restricted to the coarse grid.
                CALL mlprj_performRestriction (p_rcurrentLevel%rinterlevelProjection,&
                      p_rlowerLevel%rsolutionVector, &
                      p_rcurrentLevel%rtempVector, &
                      p_rsubnode%rprjTempVector)

                ! Apply the filter chain (e.g. to implement boundary conditions)
                ! on just calculated right hand side
                ! (which is a defect vector against the 0-vector).
                IF (bfilter) THEN
                  ! Apply the filter chain to the vector.
                  ! We are in 'unsorted' state; applying the filter here is
                  ! supposed to be a bit faster...
                  CALL filter_applyFilterChainVec (p_rlowerLevel%rsolutionVector, &
                                                  p_RfilterChain)
                END IF

                IF (bsort) THEN
                  ! Resort the RHS on the lower level.
                  CALL lsysbl_sortVectorInSitu (p_rlowerLevel%rsolutionVector,&
                        p_rsubnode%rprjTempVector,.TRUE.)

                  ! Temp-vector and RHS can be ignored on the coarse grid.
                END IF

              END IF              
            
              ! Go down one level
              ilev = ilev - 1
              p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
              p_rlowerLevel => p_rcurrentLevel%p_rprevLevel
              
              ! If we are not on the lowest level, repeat the smoothing of 
              ! the solution/restriction of the new defect in the next loop 
              ! pass...
            END DO   ! ilev > minimum level
            
            ! Now we reached the coarse grid.
            ! Apply the filter chain (e.g. to implement boundary conditions)
            ! on the just calculated right hand side,
            ! which is currently located in rsolutionVector.
            IF (bfilter) THEN
              ! Apply the filter chain to the vector
              CALL filter_applyFilterChainVec (p_rcurrentLevel%rsolutionVector, &
                                               p_RfilterChain)
            END IF
            
            ! Solve the system on lowest level by preconditioning
            ! of the RHS=defect vector.
            CALL linsol_precondDefect(p_rcurrentLevel%p_rcoarseGridSolver,&
                                      p_rcurrentLevel%rsolutionVector)
            
            ! Now we have the solution vector on the lowest level - we have to go
            ! upwards now... but probably not to NLMAX! That depends on the cycle.
            !
            DO WHILE(ASSOCIATED(p_rcurrentLevel%p_rnextLevel))
              ilev = ilev + 1
              p_rlowerLevel => p_rcurrentLevel
              p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
              
              ! Prolongate the solution vector from the coarser level
              ! to the temp vector on the finer level.
              !
              ! Make sure the vector is unsorted before the prolongation!
              ! We can use the 'global' temporary array p_rprjTempVector
              ! as temporary memory during the resorting process.
              IF (bsort) THEN
                CALL lsysbl_sortVectorInSitu (p_rlowerLevel%rsolutionVector,&
                      p_rsubnode%rprjTempVector,.FALSE.)
                
                ! Temp-vector and RHS vector there unused now,
                ! therefore we can mark them as not sorted without calling the
                ! resorting routine.
                ! This prepares these vectors for the next sweep, when an unsorted
                ! vector comes 'from above'.
                p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                  -ABS(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)% &
                  isortStrategy) 
                p_rlowerLevel%rrhsVector%RvectorBlock(1:nblocks)% &
                  isortStrategy = -ABS(p_rlowerLevel%rrhsVector% &
                  RvectorBlock(1:nblocks)%isortStrategy) 
              END IF
              CALL mlprj_performProlongation (p_rcurrentLevel%rinterlevelProjection,&
                    p_rlowerLevel%rsolutionVector, &
                    p_rcurrentLevel%rtempVector, &
                    p_rsubnode%rprjTempVector)

              IF (bfilter) THEN
                ! Apply the filter chain to the vector.
                CALL filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                                 p_RfilterChain)
              END IF

              IF (bsort) THEN
                ! Resort the temp vector on the current level so that it fits
                ! to the RHS, solution vector and matrix again
                CALL lsysbl_sortVectorInSitu (p_rcurrentLevel%rtempVector,&
                      p_rsubnode%rprjTempVector,.TRUE.)
              END IF
              
              ! Step length control. Get the optimal damping parameter for the
              ! defect correction.
              CALL cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrix,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rsubnode%rprjTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              
              ! Perform the coarse grid correction by adding the coarse grid
              ! solution (with the calculated step-length parameter) to
              ! the current solution.
              
              CALL lsysbl_vectorLinearComb (p_rcurrentLevel%rtempVector,&
                                            p_rcurrentLevel%rsolutionVector,&
                                            dstep,1.0_DP)
                                            
              ! Perform the post-smoothing with the current solution vector
              IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
                CALL linsol_smoothCorrection (p_rcurrentLevel%p_rpostSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)
              END IF

              ! Update the iteration counter(s) for realising the MG-cycle(s).
              ! Then either repeat this loop to perform the next prolongation or   
              ! repeat the cycleloop to do perform a next MG sweep on the current      
              ! level.                                                        
              !                                                               
              ! Here icycle defines how the cycle-counters are updated.           
              ! For a W-cycle the cycle counter is resetted to 2 if the sweep is 
              ! fulfilled on the current level, for F-cycle it's set to 1 to not
              ! perform more that 1 cycle on the current level anymore.       
              
              p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncyclesRemaining-1
              IF (p_rcurrentLevel%ncyclesRemaining .LE. 0) THEN
                IF (p_rsubnode%icycle .EQ. 0) THEN
                  p_rcurrentLevel%ncyclesRemaining = 1
                ELSE
                  p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles
                  CYCLE cycleloop
                END IF
              END IF
              
            END DO ! ilev < nlmax
            
            ! We finally reached the maximum level. Quit the cycle loop
            ! to go on with the next iteration.
            EXIT cycleloop
          
          END DO cycleloop
          
          ! We have (hopefully) successfully performed one MG-sweep, starting
          ! and ending on the finest level. As we are now on the finest level
          ! again, we can update our defect vector to test the current
          ! residuum...
          !
          ! Calculate the residuum and its norm.
          CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rtempVector)
          CALL lsysbl_blockMatVec (&
                p_rcurrentLevel%rsystemMatrix, &
                p_rcurrentLevel%rsolutionVector,&
                p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
          IF (bfilter) THEN
            ! Apply the filter chain to the vector
            CALL filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                             p_RfilterChain)
          END IF
          dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,rsolverNode%iresNorm)
          IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                    (dres .LE. 1E99_DP))) dres = 0.0_DP
          
          ! Shift the queue with the last residuals and add the new
          ! residual to it. Check length of ireslength to be larger than
          ! 0 as some compilers might produce Floating exceptions
          ! otherwise! (stupid pgf95)
          IF (ireslength .GT. 0) &
            dresqueue(1:ireslength) = EOSHIFT(dresqueue(1:ireslength),1,dres)

          rsolverNode%dfinalDefect = dres
          
          ! Test if the iteration is diverged
          IF (linsol_testDivergence(rsolverNode,dres)) THEN
            PRINT *,'Multigrid: Solution diverging!'
            rsolverNode%iresult = 1
            EXIT
          END IF
       
          ! At least perform nminIterations iterations
          IF (ite .GE. nminIterations) THEN
          
            ! Check if the iteration converged
            IF (linsol_testConvergence(rsolverNode,dres)) EXIT
            
          END IF

          ! print out the current residuum

          IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
              (MOD(ite,niteResOutput).EQ.0)) THEN
            !WRITE (MTERM,'(A,I7,A,D25.16)') 
            PRINT *,'Multigrid: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
          END IF
          
        END DO  ! ite
        
        ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
        ! completed

        IF (ite .GT. rsolverNode%nmaxIterations) &
          ite = rsolverNode%nmaxIterations

        ! Finish - either with an error or if converged.
        ! Print the last residuum.

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
            (rsolverNode%iresult .GE. 0)) THEN
          !WRITE (MTERM,'(A,I7,A,D25.16)') 
          PRINT *,'Multigrid: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
        END IF
        
      END IF

      rsolverNode%iiterations = ite
      
      ! Finally, we look at some statistics:
      !
      ! Don't calculate anything if the final residuum is out of bounds -
      ! would result in NaN's,...
        
      IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
        ! Calculate asymptotic convergence rate
      
        IF (dresqueue(1) .GE. 1E-70_DP) THEN
          I = MAX(1,MIN(rsolverNode%iiterations,ireslength-1))
          rsolverNode%dasymptoticConvergenceRate = &
            (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(I,DP))
        END IF

        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting multigrid convergence rate stays zero.
        ! In the other case the multigrid convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :

        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/REAL(rsolverNode%iiterations,DP))
        END IF
      END IF
    
    END IF

    ! As the solution vector on the finest level shared its memory with rd,
    ! we just calculated the new correction vector!
      
    IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A)') 
        PRINT *,'Multigrid statistics:'
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A,I5)')     'Iterations              : ',
    !*            IPARAM(OITE)
        PRINT *,'Iterations              : ',rsolverNode%iiterations
        !WRITE (MTERM,'(A,D24.12)') '!!INITIAL RES!!         : ',
    !*            DPARAM(ODEFINI)
        PRINT *,'!!INITIAL RES!!         : ',rsolverNode%dinitialDefect
        !WRITE (MTERM,'(A,D24.12)') '!!RES!!                 : ',
    !*            DPARAM(ODEFFIN)
        PRINT *,'!!RES!!                 : ',rsolverNode%dfinalDefect
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
    !*              rparam%dfinalDefect / rparam%dinitialDefect
          PRINT *,'!!RES!!/!!INITIAL RES!! : ',&
                  rsolverNode%dfinalDefect / rsolverNode%dinitialDefect
        ELSE
          !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
    !*              0D0
          PRINT*,'!!RES!!/!!INITIAL RES!! : ',0.0_DP
        END IF
        !WRITE (MTERM,'(A)') ''
        PRINT *
        !WRITE (MTERM,'(A,D24.12)') 'Rate of convergence     : ',
    !*            DPARAM(ORHO)
        PRINT *,'Rate of convergence     : ',rsolverNode%dconvergenceRate

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'Multigrid: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'Multigrid: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
        PRINT *,&
              'MG: Iterations/Rate of convergence: ',&
              rsolverNode%iiterations,' /',rsolverNode%dconvergenceRate
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
END MODULE
