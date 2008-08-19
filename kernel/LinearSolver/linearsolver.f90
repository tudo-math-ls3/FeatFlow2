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
!#                                      *Not on the coarse grid!'
!#     b) linsol_addMultigridLevel    - Add a multigrid level. First the coarse grid
!#                                      must be added, then level 2, then level 3 etc.
!#
!#     Steps a), b) must be repeated for all levels.
!#
!#     Btw., levels can be removed from the solver with linsol_removeMultigridLevel.
!#     The level stucture can be obtained with linsol_getMultigridLevel.
!#
!# 3.) linsol_matricesCompatible      - Check if the system matrices are compatible
!#                                      to the solver; this can be omitted if one knows
!#                                      the solver and its properties.
!#     linsol_setMatrices             - Attach the system matrices of all levels to the
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
!# 1.) linsol_addMultigridLevel
!#     -> Adds a new level to the Multigrid structure
!#
!# 2.) linsol_removeMultigridLevel
!#     -> Deletes a level and attached solver structures of that level from 
!#        Multigrid
!#
!# 3.) linsol_cleanMultigridLevels
!#     -> Deletes all levels and attached solver structures
!#
!# 4.) linsol_getMultigridLevel
!#     -> Returns the level structure of a specific level
!#
!# 5.) linsol_getMultigridLevelCount
!#     -> Returns number of currently attached levels
!#
!# 6.) linsol_convertToSmoother
!#     -> Converts a solver node into a smoother for multigrid smoothing.
!#
!# 7.) linsol_alterSolver
!#     -> Allows direct modification of an existing solver. Passes a command
!#        block to all subsolvers in a solver tree.
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
!#
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
!#      Well, all solvers are build up in the same way with the same
!#      interfaces - that's the reason why every solver can be used as
!#      a preconditioner in another one. If you want to create your
!#      own solver and plug it into here, use the following guideline
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
!#
!#  8.) Which solvers are available at all at the moment?
!#
!#      Well, search for the linsol_initXXXX routines :-)
!#      Currently, we have:
!#
!#       1.) linsol_initDefCorr
!#           -> standard defect correction loop
!#
!#       2.) linsol_initJacobi
!#           -> standard Jacobi preconditioner
!#           -> see [http://www.netlib.org/linalg/html_templates/Templates.html],
!#                  [http://mathworld.wolfram.com/JacobiMethod.html]
!#
!#       3.) linsol_initJinWeiTam
!#           -> Jin-Wei-Tam-preconditioner
!#           -> see [X.-Q. Jin, Y.-M. Wei, H.-S. Tam;
!#                   Preconditioning techniques for symmetric $M$-matrices;
!#                   2005; Calcolo 42, p. 105-113; 
!#                   DOI: 10.1007/s10092-005-0100-6]
!#
!#       5.) linsol_initSOR
!#           -> SOR preconditioner / GS preconditioner (=SOR(1.0))
!#           -> see [http://www.netlib.org/linalg/html_templates/Templates.html],
!#                  [http://mathworld.wolfram.com/SuccessiveOverrelaxationMethod.html]
!#
!#       5.) linsol_initSSOR
!#           -> SSOR preconditioner
!#           -> see [http://www.netlib.org/linalg/html_templates/Templates.html],
!#                  [http://mathworld.wolfram.com/SuccessiveOverrelaxationMethod.html]
!#
!#       6.) linsol_initVANCA
!#           -> VANCA preconditioner; multiple versions
!#
!#       7.) linsol_initUMFPACK4
!#           -> UMFPACK preconditioner
!#
!#       8.) linsol_initMILUs1x1
!#           -> (M)ILU-preconditioner for 1x1-matrices from SPLIB
!#           -> see [David Hysom and A. Pothen; Level-based Incomplete LU 
!#                   factorization: Graph Model and Algorithms; 
!#                   Tech Report UCRL-JC-150789; Lawrence Livermore National Labs;
!#                   Nov 2002; http://www.cs.odu.edu/~pothen/papers.html]
!#
!#       9.) linsol_initBiCGStab
!#           -> BiCGStab preconditioner
!#           -> see [van der Vorst, H.A.; BiCGStab: A Fast and Smoothly Converging
!#                   Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems;
!#                   SIAM J. Sci. Stat. Comput. 1992, Vol. 13, No. 2, pp. 631-644]
!#
!#      10.) linsol_initCG
!#           -> Conjugate Gradient preconditioner
!#
!#      11.) linsol_initGMRES
!#           -> (flexible) Generalized Minimal-Residual preconditioner
!#           -> see [Y. Saad; A flexible inner outer preconditioned GMRES algorithm;
!#                   SIAM J. Sci. Comput. 1993, Vol. 14, No. 2, pp. 461-469]
!#                  [http://www-users.cs.umn.edu/~saad/]
!#                  [http://www.netlib.org/linalg/html_templates/Templates.html]
!#
!#      12.) linsol_initMultigrid
!#           -> Multigrid preconditioner
!#
!#      13.) linsol_initMultigrid2
!#           -> Multigrid preconditioner; alternative implementation using
!#              an array to store level information rather than a linear list.
!#
!# 9.) What is this linsol_matricesCompatible?
!#
!#     This allows to check a set of system matrices against a solver node.
!#     If the matrices are not compatible (e.g. because their properties like
!#     matrix format, being transposed,... do not fit) an error code is
!#     returned -- as this set of matrices would lead the solver to an
!#     immediate program stop if they are factorised with linsol_initStructure!
!#     It may be used for special solvers like VANCA to check, if the solver
!#     can handle what is set up in the discretisation.
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
  USE genoutput
  USE iluk
  
  USE matrixio
  
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
  
  ! SOR/GS iteration $x_1 = x_0 + (D+\omega L)^{-1}(b-Ax_0)$
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
  
  ! Jin-Wei-Tam iteration
  INTEGER, PARAMETER :: LINSOL_ALG_JINWEITAM     = 12
  
  ! Multigrid iteration
  INTEGER, PARAMETER :: LINSOL_ALG_MULTIGRID2    = 13
  
  ! EMS iteration
  INTEGER, PARAMETER :: LINSOL_ALG_EMS           = 17

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

! *****************************************************************************

!<constantblock description="Flags for matrix compatibility check.">

  ! Matrices are compatible to a solver node.
  INTEGER, PARAMETER :: LINSOL_COMP_OK               = 0
  
  ! Matrices generally not compatible, no special reason.
  INTEGER, PARAMETER :: LINSOL_COMP_ERRGENERAL       = 1
  
  ! One of the submatrices is solved in transposed structure,
  ! but one of the subsolvers cannot handle that.
  INTEGER, PARAMETER :: LINSOL_COMP_ERRTRANSPOSED    = 2

  ! Some of the submatrices are not saved transposed although
  ! they have to be. (Usual error for special VANCA solvers.)
  INTEGER, PARAMETER :: LINSOL_COMP_ERRNOTTRANSPOSED = 3

  ! One of the submatrices is a block matrix,
  ! but one of the subsolvers can only handle scalar matrices.
  INTEGER, PARAMETER :: LINSOL_COMP_ERRNOTSCALAR     = 4
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Identifiers for stopping criterium istoppingCriterion.">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. the iteration stops when both,
  !    the relative AND the absolute stopping criterium holds
  INTEGER, PARAMETER :: LINSOL_STOP_STANDARD     = 0

  ! Use 'minimum' stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use one of them, i.e. the iteration stops when the
  !    either the relative OR the absolute stopping criterium holds
  INTEGER, PARAMETER :: LINSOL_STOP_ONEOF        = 1
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Bitfield identifiers for the ability of a linear solver">

  ! Solver can handle scalar systems
  INTEGER(I32), PARAMETER :: LINSOL_ABIL_SCALAR       = 2**0

  ! Solver can handle block systems
  INTEGER(I32), PARAMETER :: LINSOL_ABIL_BLOCK        = 2**1

  ! Solver can handle multiple levels
  INTEGER(I32), PARAMETER :: LINSOL_ABIL_MULTILEVEL   = 2**2
  
  ! Solver allows checking the defect during the iteration.
  ! Solvers not capable of this perform only a fixed number of solution
  ! steps (e.g. UMFPACK performs always one step).
  INTEGER(I32), PARAMETER :: LINSOL_ABIL_CHECKDEF     = 2**3
  
  ! Solver is a direct solver (e.g. UMFPACK, ILU).
  ! Otherwise the solver is of iterative nature and might perform
  ! multiple steps to solve the problem.
  INTEGER(I32), PARAMETER :: LINSOL_ABIL_DIRECT       = 2**4
  
  ! Solver might use subsolvers (preconditioners, smoothers,...)
  INTEGER(I32), PARAMETER :: LINSOL_ABIL_USESUBSOLVER = 2**5
  
  ! Solver supports filtering
  INTEGER(I32), PARAMETER :: LINSOL_ABIL_USEFILTER    = 2**6

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

  ! Matrix-structure has changed between initStructure and initData
  INTEGER, PARAMETER :: LINSOL_ERR_MATRIXHASCHANGED = 3
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Variants of the VANCA solver">

  ! General VANCA solver
  INTEGER, PARAMETER :: LINSOL_VANCA_GENERAL           = 0   
  
  ! General VANCA solver. Specialised 'direct' version, i.e. when 
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup. 
  INTEGER, PARAMETER :: LINSOL_VANCA_GENERALDIRECT     = 1

  ! Simple VANCA, 2D Navier-Stokes problem, general discretisation
  INTEGER, PARAMETER :: LINSOL_VANCA_2DNAVST           = 2

  ! Simple VANCA, 2D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when 
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup. 
  INTEGER, PARAMETER :: LINSOL_VANCA_2DNAVSTDIRECT     = 3

  ! Full VANCA, 2D Navier-Stokes problem, general discretisation
  INTEGER, PARAMETER :: LINSOL_VANCA_2DFNAVST          = 4

  ! Full VANCA, 2D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when 
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup. 
  INTEGER, PARAMETER :: LINSOL_VANCA_2DFNAVSTDIRECT    = 5

  ! Full VANCA, 2D Navier-Stokes optimal control problem, general discretisation.
  INTEGER, PARAMETER :: LINSOL_VANCA_2DFNAVSTOC        = 20

  ! Full VANCA, 2D Navier-Stokes optimal control problem, general discretisation.
  ! Specialised 'direct' version, i.e. when 
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup. 
  INTEGER, PARAMETER :: LINSOL_VANCA_2DFNAVSTOCDIRECT  = 21

  ! Diagonal VANCA, 2D Navier-Stokes optimal control problem, general discretisation.
  INTEGER, PARAMETER :: LINSOL_VANCA_2DFNAVSTOCDIAG    = 22

  ! Diagonal VANCA, 2D Navier-Stokes optimal control problem, general discretisation.
  ! Specialised 'direct' version, i.e. when 
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup. 
  INTEGER, PARAMETER :: LINSOL_VANCA_2DFNAVSTOCDIAGDIR = 23

  ! Simple VANCA, 3D Navier-Stokes problem, general discretisation
  INTEGER, PARAMETER :: LINSOL_VANCA_3DNAVST           = 30

  ! Simple VANCA, 3D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when 
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup. 
  INTEGER, PARAMETER :: LINSOL_VANCA_3DNAVSTDIRECT     = 31

  ! Full VANCA, 3D Navier-Stokes problem, general discretisation
  INTEGER, PARAMETER :: LINSOL_VANCA_3DFNAVST          = 32

  ! Full VANCA, 3D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when 
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup. 
  INTEGER, PARAMETER :: LINSOL_VANCA_3DFNAVSTDIRECT    = 33

!</constantblock>

! *****************************************************************************

!<constantblock description="Variants of the BiCGStab solver">

  ! BiCGStab with Left-Preconditioning (or without preconditioning)
  INTEGER, PARAMETER :: LINSOL_BICGSTAB_LEFT_PRECOND           = 0
  
  ! BiCGStab with Right-Preconditioning
  INTEGER, PARAMETER :: LINSOL_BICGSTAB_RIGHT_PRECOND          = 1
  
!</constantblock>

! *****************************************************************************

!<constantblock description="default values for the GMRES(m) solver">

  ! One or two Gram-Schmidt calls per GMRES iteration
  LOGICAL, PARAMETER :: LINSOL_GMRES_DEF_TWICE_GS              = .FALSE.
    
!</constantblock>

! *****************************************************************************

!<constantblock description="Possible commands for linsol_slterSolver">

  ! Dummy command, do-nothing.
  INTEGER, PARAMETER :: LINSOL_ALTER_NOTHING                   = 0
  
  ! Change VANCA subtype.
  ! In the configuration block, there must be specified:
  ! Iconfig(1) = identifier for VANCA subtype to be changed.
  ! Iconfig(2) = destination subtype of VANCA solver
  INTEGER, PARAMETER :: LINSOL_ALTER_CHANGEVANCA               = 1
    
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
    REAL(DP)                        :: ddivRel = 1E6_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY disables the absolute divergence check.
    ! standard = SYS_INFINITY
    REAL(DP)                        :: ddivAbs = SYS_INFINITY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! RHS-vector is treated as zero if max(defect) < drhsZero
    REAL(DP)                        :: drhsZero = 1E-90_DP

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
    ! =-1: no output, =0: no output except for warning messages, 
    ! =1: basic output, =2, extended output
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

    ! Pointer to a structure for the Multigrid solver; NULL() if not set
    TYPE (t_linsolSubnodeMultigrid2), POINTER     :: p_rsubnodeMultigrid2  => NULL()

    ! Pointer to a structure for the ILU0 1x1 solver; NULL() if not set
    TYPE (t_linsolSubnodeILU01x1), POINTER        :: p_rsubnodeILU01x1     => NULL()
    
    ! Pointer to a structure for (M)ILUs 1x1 solver; NULL() if not set
    TYPE (t_linsolSubnodeMILUs1x1), POINTER       :: p_rsubnodeMILUs1x1    => NULL()

    ! Pointer to a structure for SSOR; NULL() if not set
    TYPE (t_linsolSubnodeSSOR), POINTER           :: p_rsubnodeSSOR        => NULL()
    
    ! Pointer to a structure for Jin-Wei-Tam; NULL() if not set
    TYPE (t_linsolSubNodeJinWeiTam), POINTER      :: p_rsubnodeJinWeiTam   => NULL()
    
    ! Pointer to a structure for the CG solver; NULL() if not set
    TYPE (t_linsolSubnodeCG), POINTER             :: p_rsubnodeCG          => NULL()
    
    ! Pointer to a structure for the GMRES(m) solver; NULL() if not set
    TYPE (t_linsolSubnodeGMRES), POINTER          :: p_rsubnodeGMRES       => NULL()
        
    ! Pointer to a structure for the EMS solver; NULL() if not set
    TYPE (t_linsolSubnodeEMS), POINTER            :: p_rsubnodeEMS         => NULL()

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
  
  ! This structure realises the subnode for the SSOR solver.
  
  TYPE t_linsolSubnodeSSOR
  
    ! Scaling of the preconditioned solution.
    ! =FALSE is the old FEAT style.
    ! =TRUE gives the implementation as suggested in the literature.
    LOGICAL :: bscale
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the defect correction solver.
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
  
  ! This structure realises the subnode for the EMS solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_linsolSubnodeEMS
  
    ! Temporary vector to use during the solution process
    TYPE(t_vectorBlock) :: rtempVector

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
  
    ! For general 2D/3D Navier-Stokes problem
    TYPE(t_vanca)       :: rvanca
  
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
  
    ! Type of preconditioning. One of the LINSOL_BICGSTAB_xxxx constants.
    INTEGER :: cprecondType = LINSOL_BICGSTAB_LEFT_PRECOND
  
    ! Temporary vectors to use during the solution process
    TYPE(t_vectorBlock), DIMENSION(6) :: RtempVectors
    
    ! Another temporary vector for right- and symmetrical preconditioning
    TYPE(t_vectorBlock) :: rprecondTemp
    
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
  
  ! This structure realises the subnode for the CG solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_linsolSubnodeCG
  
    ! Temporary vectors to use during the solution process
    TYPE(t_vectorBlock), DIMENSION(4) :: RtempVectors

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
  
    ! Matrix output for debug.
    ! If this is set to a value <> 0, the numerical factorisation routine
    ! writes the matrix to a text file before factorising it.
    ! The text file gets the name 'matrixN.txt' with N=imatrixDebugOutput.
    ! This is for debugging purposes and should be used with care,
    ! as the text files grow rather quickly with the dimension!
    INTEGER :: imatrixDebugOutput = 0
  
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
  
  ! This structure realises the subnode for the Jin-Wei-Tam solver.
  
  TYPE t_linsolSubNodeJinWeiTam
  
    ! Sum of all Matrix elements
    REAL(DP) :: dmatrixSum
    
    ! some sort of relax factor
    REAL(DP) :: drelax
  
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the flexible GMRES(m) solver.
  
  TYPE t_linsolSubnodeGMRES
  
    ! Maximum Restarts
    INTEGER :: imaxRestarts

    ! Dimension of Krylov Subspace
    INTEGER :: ikrylovDim
    
    ! Maybe we should apply Gram-Schmidt twice?
    LOGICAL :: btwiceGS
    
    ! Scale factor for pseudo-residuals
    REAL(DP) :: dpseudoResScale
    
    ! Some temporary 1D/2D arrays
    REAL(DP), DIMENSION(:),   POINTER :: Dc, Ds, Dq
    REAL(DP), DIMENSION(:,:), POINTER :: Dh

    ! The handles of the arrays
    INTEGER :: hDc, hDs, hDq, hDh
    
    ! Some temporary vectors
    TYPE(t_vectorBlock), DIMENSION(:), POINTER :: p_rv                => NULL()
    TYPE(t_vectorBlock), DIMENSION(:), POINTER :: p_rz                => NULL()
    TYPE(t_vectorBlock) :: rx
  
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
  
  ! This saves information about ILU(k) matrices as defined in SPLIB.
  
  TYPE t_linsolSubnodeMILUs1x1

    ! Fill-in level for decomposition
    INTEGER :: ifill
    
    ! Relaxation factor for (M)ILU(s)
    REAL(DP) :: drelax
    
    ! Scaling factor of the matrix; usually = 1.0
    REAL(DP) :: dscaleFactor

    ! A structure to store the (M)ILU-decomposition
    ! the decomposition is stored in a modified sparse
    ! row format (MSR)
    TYPE(t_MILUdecomp) :: rMILUdecomp

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
    
    ! t_matrixBlock structure for the calculation of the optimal correction.
    ! Normally, this matrix has not to be set. If not defined, rsystemMatrix is
    ! used. If defined, the calculation of the optimal defect (guided by
    ! the parameter rcoarseGridCorrection in the multigrid subnode)
    ! is done with this matrix instead of the system matrix.
    TYPE(t_matrixBlock)                :: rsystemMatrixOptCorrection
    
    ! A temporary vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_vectorBlock)                 :: rtempVector
    
    ! A temporary vector that is used for the calculation of the 
    ! coarse grid correction. This vector shares its data with the
    ! rprjTempVector vector in the multigrid solver structure 
    ! t_linsolSubnodeMultigrid, but has the shape of the RHS vector.
    ! The structure is initialised in initStructure and cleaned up
    ! in doneStructure.
    TYPE(t_vectorBlock)                 :: rcgcorrTempVector

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
    
    ! STATUS/INTERNAL: initial residuum when a solution is restricted to this
    ! level. Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    REAL(DP)                       :: dinitResCycle = 0.0_DP
    
    ! STATUS/INTERNAL: Number of current cycle on that level. 
    ! Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    INTEGER                        :: icycleCount   = 0
    
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
  
    ! INPUT PARAMETER: Cycle identifier. 
    !  0=F-cycle, 
    !  1=V-cycle, 
    !  2=W-cycle.
    INTEGER                       :: icycle                   = 0
    
    ! INPUT PARAMETER: Adaptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all 
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsRelCycle * |initial res on that level|.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    REAL(DP)                      :: depsRelCycle             = 1E99_DP
    
    ! INPUT PARAMETER: If adaptive cycles are activated by depsRelCycle < 1E99_DP,
    ! this configures the maximum number of cycles that are performed.
    ! The value =-1 deactivates this upper bouns, so coarse grid cycles are
    ! repeated until the convergence criterion given by depsRelCycle is reached.
    INTEGER                       :: nmaxAdaptiveCycles       = -1
    
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

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_linsolSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  TYPE t_linsolMGLevelInfo2
  
    ! A level identifier. Can be set to identify the global number of the
    ! level, this structure refers to. Only for output purposes.
    INTEGER                        :: ilevel                 = 0
    
    ! t_matrixBlock structure that holds the system matrix.
    TYPE(t_matrixBlock)                :: rsystemMatrix
    
    ! t_matrixBlock structure for the calculation of the optimal correction.
    ! Normally, this matrix has not to be set. If not defined, rsystemMatrix is
    ! used. If defined, the calculation of the optimal defect (guided by
    ! the parameter rcoarseGridCorrection in the multigrid subnode)
    ! is done with this matrix instead of the system matrix.
    TYPE(t_matrixBlock)                :: rsystemMatrixOptCorrection
    
    ! A temporary vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_vectorBlock)                 :: rtempVector
    
    ! A temporary vector that is used for the calculation of the 
    ! coarse grid correction. This vector shares its data with the
    ! rprjTempVector vector in the multigrid solver structure 
    ! t_linsolSubnodeMultigrid, but has the shape of the RHS vector.
    ! The structure is initialised in initStructure and cleaned up
    ! in doneStructure.
    TYPE(t_vectorBlock)                 :: rcgcorrTempVector

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
    
    ! STATUS/INTERNAL: initial residuum when a solution is restricted to this
    ! level. Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    REAL(DP)                       :: dinitResCycle = 0.0_DP
    
    ! STATUS/INTERNAL: Number of current cycle on that level. 
    ! Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    INTEGER                        :: icycleCount   = 0
    
  END TYPE
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  TYPE t_linsolSubnodeMultigrid2
  
    ! INPUT PARAMETER: Cycle identifier. 
    !  0=F-cycle, 
    !  1=V-cycle, 
    !  2=W-cycle.
    INTEGER                       :: icycle                   = 0
    
    ! INPUT PARAMETER: Adaptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all 
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsRelCycle * |initial res on that level|.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    REAL(DP)                      :: depsRelCycle             = 1E99_DP
    
    ! INPUT PARAMETER: If adaptive cycles are activated by depsRelCycle < 1E99_DP,
    ! this configures the maximum number of cycles that are performed.
    ! The value =-1 deactivates this upper bouns, so coarse grid cycles are
    ! repeated until the convergence criterion given by depsRelCycle is reached.
    INTEGER                       :: nmaxAdaptiveCycles       = -1
    
    ! INPUT PARAMETER: Coarse grid correction structure for step length control. 
    ! Defines the algorithm for computing the optimal correction as well as the
    ! minimum and maximum step length ALPHAMIN/ALPHAMAX.
    ! The standard setting/initialisation is suitable for conforming elements.
    TYPE(t_coarseGridCorrection)  :: rcoarseGridCorrection
    
    ! A pointer to the level info structures for all the levels in multigrid
    TYPE(t_linsolMGLevelInfo2), DIMENSION(:), POINTER :: p_RlevelInfo
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:),POINTER :: p_RfilterChain     => NULL()
  
    ! A temp vector for the prolongation/restriction.
    ! Memory for this vector is allocated in initStructure and released
    ! in doneStructure.
    TYPE(t_vectorScalar) :: rprjTempVector
    
  END TYPE
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>

  ! This structure realises a configuration block that can be passed to the
  ! function linsol_alterSolver to allow in-computation change of
  ! a solver or its subsolvers.
  TYPE t_linsol_alterSolverConfig
  
    ! A command flag of type LINSOL_ALTER_xxxx.
    INTEGER :: ccommand = LINSOL_ALTER_NOTHING
  
    ! An integer precision configuration block. The meaning of this block
    ! depends on ccommand.
    INTEGER, DIMENSION(16) :: Iconfig
    
    ! A double precision configuration block. The meaning of this block
    ! depends on ccommand.
    REAL(DP), DIMENSION(16) :: Dconfig
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

  ! Copy the matrix structure on the finest level to the rsolverNode
  ! structure. This corresponds to the system we want to solve.
  CALL lsysbl_duplicateMatrix (Rmatrices(UBOUND(Rmatrices,1)), &
      rsolverNode%rsystemMatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
  
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
  CASE (LINSOL_ALG_CG)
    CALL linsol_setMatrixCG (rsolverNode, Rmatrices)
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_setMatrixBiCGStab (rsolverNode, Rmatrices)
  CASE (LINSOL_ALG_GMRES)
    CALL linsol_setMatrixGMRES (rsolverNode, Rmatrices)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_setMatrixMultigrid (rsolverNode, Rmatrices)
  CASE (LINSOL_ALG_MULTIGRID2)
    CALL linsol_setMatrixMultigrid2 (rsolverNode, Rmatrices)
  CASE DEFAULT
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matricesCompatible (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  
  ! This subroutine checks the matrices in Rmatrices if they are compatible
  ! with the solver node rsolverNode. Only if the matrices are compatiable,
  ! they are allowed to be attached to the solver with linsol_setMatrices --
  ! otherwise one of the subsolvers will halt the program because not being
  ! able to handle one of the matrices.
  !
  ! The last element in Rmatrices defines the system matrix of the
  ! linear system that is to solve. If there's a solver involved
  ! that needs matrices on multiple grids (e.g. multigrid), Rmatrices
  ! must contain the matrices of all levels. In this case, Rmatrices(1)
  ! defines the matrix on the coarsest level to solve, Rmatrices(2)
  ! the matrix of the first refinement etc. up to the level where
  ! to solve the problem.
  !
  ! ccompatible is set to LINSOL_COMP_OK if the matrices are compatible, 
  ! i.e. that the solver can handle them. If the matrices are not compatible,
  ! ccompatible is set to a LINSOL_COMP_xxxx-flag that identifies what's
  ! wrong with the matrices.
  !
  ! With the optional parameter CcompatibleDetail, the caller can get more
  ! detailed information about the compatibility of the matrices.
  ! For every level i, CcompatibleDetail(i) returns a LINSOL_COMP_xxxx
  ! flag that tells whether the matrix on that level is compatible
  ! or not.
  ! Note that it is possible to create cases where this routine ALWAYS
  ! fails, e.g. if two incompatible VANCA solvers are used as pre- and
  ! postsmoother, resp., on one leve in multigrid! But at least if pre- and
  ! postsmoother on each MG level are identical, this routine will
  ! successfully determine (and return in CcompatibleDetail) if 
  ! everything is ok.
  
!</description>
  
!<input>
  ! Array of system matrices on all levels of the discretisation.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN) :: Rmatrices

  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode
!</input>

!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible
  
  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

  ! Depending on the solver type, call the corresponding check-routine
  ! or return the corresponding matrix compatibility flag directly.
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (LINSOL_ALG_VANCA)
    ! Ask VANCA if the matrices are ok.
    CALL linsol_matCompatVANCA (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_JACOBI)
    ! Ask Jacobi if the matrices are ok.
    CALL linsol_matCompatJacobi (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_JINWEITAM)
    ! Ask JWT if the matrices are ok.
    CALL linsol_matCompatJinWeiTam (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_SOR)
    ! Ask SOR if the matrices are ok.
    CALL linsol_matCompatSOR (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_SSOR)
    ! Ask SSOR if the matrices are ok.
    CALL linsol_matCompatSSOR (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_DEFCORR)
    ! Ask the defect correction and its subsolvers if the matrices are ok.
    CALL linsol_matCompatDefCorr (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_UMFPACK4)
    ! Ask VANCA if the matrices are ok.
    CALL linsol_matCompatUMFPACK4 (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
  
  CASE (LINSOL_ALG_CG)
    ! Ask CG and its subsolvers if the matrices are ok.
    CALL linsol_matCompatCG (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_BICGSTAB)
    ! Ask BiCGStab and its subsolvers if the matrices are ok.
    CALL linsol_matCompatBiCGStab (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_GMRES)
    ! Ask GMRES(m) and its subsolvers if the matrices are ok.
    CALL linsol_matCompatGMRES (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_MULTIGRID)
    ! Ask Multigrid and its subsolvers if the matrices are ok.
    CALL linsol_matCompatMultigrid (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
  CASE (LINSOL_ALG_MULTIGRID2)
    ! Ask Multigrid and its subsolvers if the matrices are ok.
    CALL linsol_matCompatMultigrid2 (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)

  CASE (LINSOL_ALG_EMS)
    ! Ask EMS if the matrices are ok.
    CALL linsol_matCompatEMS (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)

  CASE DEFAULT
    ! Nothing special. Let's assume that the matrices are ok.
    ccompatible = LINSOL_COMP_OK
    CcompatibleDetail(:) = ccompatible
    
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
    CASE (LINSOL_ALG_CG)
      CALL linsol_initStructureCG (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_initStructureBiCGStab (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_GMRES)
      CALL linsol_initStructureGMRES (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_initStructureMultigrid (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID2)
      CALL linsol_initStructureMultigrid2 (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_EMS)
      CALL linsol_initStructureEMS (rsolverNode,ierror,isubgroup)
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
    CASE (LINSOL_ALG_CG)
      CALL linsol_initDataCG (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_initDataBiCGStab (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_GMRES)
      CALL linsol_initDataGMRES (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_JINWEITAM)
      CALL linsol_initDataJinWeiTam (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_initDataMultigrid (rsolverNode,ierror,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID2)
      CALL linsol_initDataMultigrid2 (rsolverNode,ierror,isubgroup)
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
      CALL linsol_doneDataVANCA (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_doneDataDefCorr (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_doneDataUMFPACK4 (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MILUS1X1)
      CALL linsol_doneDataMILUs1x1 (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_CG)
      CALL linsol_doneDataCG (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_doneDataBiCGStab (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_GMRES)
      CALL linsol_doneDataGMRES (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_doneDataMultigrid (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID2)
      CALL linsol_doneDataMultigrid2 (rsolverNode,isubgroup)
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
    CASE (LINSOL_ALG_CG)
      CALL linsol_doneStructureCG (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_doneStructureBiCGStab (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_GMRES)
      CALL linsol_doneStructureGMRES (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_doneStructureMultigrid (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_MULTIGRID2)
      CALL linsol_doneStructureMultigrid2 (rsolverNode,isubgroup)
    CASE (LINSOL_ALG_EMS)
      CALL linsol_doneStructureEMS (rsolverNode,isubgroup)
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
    CASE (LINSOL_ALG_SSOR)
      CALL linsol_doneSSOR (p_rsolverNode)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_doneDefCorr (p_rsolverNode)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_doneUMFPACK4 (p_rsolverNode)
    CASE (LINSOL_ALG_CG)
      CALL linsol_doneCG (p_rsolverNode)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_doneBiCGStab (p_rsolverNode)
    CASE (LINSOL_ALG_GMRES)
      CALL linsol_doneGMRES (p_rsolverNode)
    CASE (LINSOL_ALG_MILUS1X1)
      CALL linsol_doneMILUs1x1 (p_rsolverNode)
    CASE (LINSOL_ALG_JINWEITAM)
      CALL linsol_doneJinWeiTam (p_rsolverNode)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_doneMultigrid (p_rsolverNode)
    CASE (LINSOL_ALG_MULTIGRID2)
      CALL linsol_doneMultigrid2 (p_rsolverNode)
    CASE (LINSOL_ALG_EMS)
      CALL linsol_doneEMS (p_rsolverNode)
    CASE DEFAULT
    END SELECT
    
    ! Clean up the associated matrix structure.
    ! Of course, the memory of the matrix is not released from memory, because
    ! if there's a matrix attached, it belongs to the application, not to the
    ! solver!
    CALL lsysbl_releaseMatrix(p_rsolverNode%rsystemMatrix)
    
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
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_alterSolver (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of an existing solver
  ! configuration. It can be used to explicitely alter the configuration
  ! of a solver or one of the subsolvers or to fetch information from
  ! it. This is especially useful if rsolverNode represents a very complex
  ! solver configuration and the caller wants to modify a very special
  ! part of one of the subsolvers in it (like changing a solvers subtype).
  !
  ! ralterConfig is a configuration block that specifies a 'command' and
  ! a couple of integer/double precision variables. This command block
  ! is passed to all solvers in the whole rsolverNode solver tree.
  ! Each subsolver then decides on it's own what to do with the configuration,
  ! depending on the command ralterConfig%ccommand.
  ! 
  ! Example: If ralterConfig%ccommand=LINSOL_ALTER_CHANGEVANCA, all VANCA
  !  solvers and smoothers in the solver will react to the configuration
  !  block and change their type, depending on ralterConfig%Iconfig.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode

  ! A command/configuration block that is passed to all solvers in the
  ! solver tree identified by rsolverNode.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Call the alterConfig-routine of the actual solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_VANCA)
      CALL linsol_alterVANCA (rsolverNode, ralterConfig)
    CASE (LINSOL_ALG_DEFCORR)
      CALL linsol_alterDefCorr (rsolverNode, ralterConfig)
    CASE (LINSOL_ALG_UMFPACK4)
      ! No alter routine for UMFPACK
    CASE (LINSOL_ALG_MILUS1X1)
      ! No structure routine for (M)ILU(s)
    CASE (LINSOL_ALG_CG)
      CALL linsol_alterCG (rsolverNode, ralterConfig)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_alterBiCGStab (rsolverNode, ralterConfig)
    CASE (LINSOL_ALG_GMRES)
      CALL linsol_alterGMRES (rsolverNode, ralterConfig)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_alterMultigrid (rsolverNode, ralterConfig)
    CASE (LINSOL_ALG_MULTIGRID2)
      CALL linsol_alterMultigrid2 (rsolverNode, ralterConfig)
    END SELECT
  
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
    
    SELECT CASE (rsolverNode%istoppingCriterion)
    
    CASE (LINSOL_STOP_ONEOF)
      ! Iteration stops if either the absolute or the relative criterium holds.
      loutput = .FALSE.
      
      ! Absolute convergence criterion? Check the norm directly.
      IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
        IF (.NOT. (dvecNorm .GT. rsolverNode%depsAbs)) THEN
          loutput = .TRUE.
          RETURN
        END IF
      END IF
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm. 
      IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
        IF (.NOT. &
            (dvecNorm .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect)) THEN
          loutput = .TRUE.
          RETURN
        END IF
      END IF
    
    CASE DEFAULT
      ! Standard stopping criterion.
      ! Iteration stops if both the absolute and the relative criterium holds.
      loutput = .TRUE.
      
      ! Absolute convergence criterion? Check the norm directly.
      IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
        IF (dvecNorm .GT. rsolverNode%depsAbs) THEN
          loutput = .FALSE.
          RETURN
        END IF
      END IF
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm. 
      IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
        IF (dvecNorm .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect) THEN
          loutput = .FALSE.
          RETURN
        END IF
      END IF
    END SELECT
    
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
    CASE (LINSOL_ALG_SOR)
      CALL linsol_precSOR (rsolverNode,rd)
    CASE (LINSOL_ALG_SSOR)
      CALL linsol_precSSOR (rsolverNode,rd)
    CASE (LINSOL_ALG_UMFPACK4)
      CALL linsol_precUMFPACK4 (rsolverNode,rd)
    CASE (LINSOL_ALG_MILUS1x1)
      CALL linsol_precMILUS1x1 (rsolverNode,rd)
    CASE (LINSOL_ALG_CG)
      CALL linsol_precCG (rsolverNode,rd)
    CASE (LINSOL_ALG_BICGSTAB)
      CALL linsol_precBiCGStab (rsolverNode,rd)
    CASE (LINSOL_ALG_GMRES)
      CALL linsol_precGMRES (rsolverNode,rd)
    CASE (LINSOL_ALG_MULTIGRID)
      CALL linsol_precMultigrid (rsolverNode,rd)
    CASE (LINSOL_ALG_MULTIGRID2)
      CALL linsol_precMultigrid2 (rsolverNode,rd)
    CASE (LINSOL_ALG_JINWEITAM)
      CALL linsol_precJinWeiTam (rsolverNode,rd)
    CASE (LINSOL_ALG_EMS)
      CALL linsol_precEMS (rsolverNode,rd)
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
  !    $$ x_{n+1}  =  x_n  +  (b-Ax_n) $$
  ! with $x_0:=0$. 
  ! It's possible to include a damping parameter to this operation by 
  ! changing rsolverNode%domega to a value $\not =1$. In this case, the
  ! defect correction iteration changes to the Richardson iteration
  !
  !    $$ x_{n+1}  =  x_n  +  \omega(b-Ax_n) $$
  !
  ! By specifying an additional preconditioner, it's possible to perform
  ! the preconditioned defect correction iteration
  !
  !    $$ x_{n+1}  =  x_n  +  \omega P^{-1} (b-Ax_n) $$
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
  
  RECURSIVE SUBROUTINE linsol_alterDefCorr (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the Defect correction solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration 
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there's a preconditioner attached. If yes, pass the command
    ! structure to that one.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_alterSolver(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,&
          ralterConfig)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDefCorr (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver 
  ! from the heap.
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

    ! Release memory if still associated
    CALL linsol_doneDataDefCorr (rsolverNode, rsolverNode%isolverSubgroup)
    CALL linsol_doneStructureDefCorr (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeDefCorr)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatDefCorr (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL linsol_matricesCompatible (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,&
          Rmatrices,ccompatible,CcompatibleDetail)
    END IF
    
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
          rsolverNode%p_rsubnodeDefCorr%rtempVector,.FALSE.,.FALSE.,&
          rsolverNode%cdefaultDataType)

    CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeDefCorr%rtempVector2,.FALSE.,.FALSE.,&
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
  
  RECURSIVE SUBROUTINE linsol_doneDataDefCorr (rsolverNode, isolverSubgroup)
  
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
  INTEGER :: ireslength,ite,i,j,niteAsymptoticCVR
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

    ireslength = 32 
    niteAsymptoticCVR = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

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
        CALL output_line ('DefCorr: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
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
          IF (rsolverNode%ioutputLevel .LT. 2) THEN
            DO i=MAX(1,SIZE(Dresqueue)-ite-1+1),SIZE(Dresqueue)
              j = ITE-MAX(1,SIZE(Dresqueue)-ite+1)+i
              CALL output_line ('DefCorr: Iteration '// &
                  TRIM(sys_siL(j,10))//',  !!RES!! = '//&
                  TRIM(sys_sdEL(Dresqueue(i),15)) )
            END DO
          END IF
          CALL output_line ('DefCorr: Solution diverging!')
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
          CALL output_line ('DefCorr: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF

      END DO

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) THEN
        ! Warning if we didn't reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        IF ((rsolverNode%ioutputLevel .GE. 0) .AND. &
            (rsolverNode%nmaxIterations .GT. rsolverNode%nminIterations)) THEN
          CALL output_line ('DefCorr: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        END IF

        ite = rsolverNode%nmaxIterations
      END IF

      ! Finish - either with an error or if converged.
      ! Print the last residuum.

      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        CALL output_line ('DefCorr: Iteration '// &
            TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
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
    
      IF (niteAsymptoticCVR .NE. 0) THEN
        I = MAX(33-ite,33-niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(33-I,DP))
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
        CALL output_lbrk()
        CALL output_line ('DefCorr statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'DefCorr: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE

  
! *****************************************************************************
! Routines for the EMS iteration
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initEMS (p_rsolverNode)
  
!<description>
  ! Creates a t_linsolNode solver structure for the EMS iteration.
!</description>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
!</output>

!</subroutine>
  
    ! Create a default solver structure
    CALL linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_EMS
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK    + &
                                LINSOL_ABIL_DIRECT
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeEMS)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneEMS (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the EMS solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of EMS which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    !CALL linsol_doneDataEMS (rsolverNode, rsolverNode%isolverSubgroup)
    CALL linsol_doneStructureEMS (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeEMS)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatEMS (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
  END SUBROUTINE
  
! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureEMS (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the EMS solver
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
    
    ! Cancel here, if we don't belong to the subgroup to be initialised
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
    
    ! Intialisation. In our case: allocate temporary vector for our data
    ! by using the associated matrix as template.
    CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeEMS%rtempVector,.FALSE.,.FALSE.,&
          rsolverNode%cdefaultDataType)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureEMS (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the EMS solver
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
    
    ! Cancel here, if we don't belong to the subgroup to be initialised
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

    ! Ok, we are the one to be released. We have a temp-vector to be released!
    IF (rsolverNode%p_rsubnodeEMS%rtempVector%NEQ .NE. 0) THEN
      CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeEMS%rtempVector)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precEMS (rsolverNode,rd)
  
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

  ! Damping parameter
  REAL(DP) :: domega, dalpha, dbeta

  ! The system matrix
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  
  ! A pointer to our temporary vector
  TYPE(t_vectorBlock), POINTER :: p_rt
  
  ! The local subnode
  TYPE(t_linsolSubnodeEMS), POINTER :: p_rsubnode
  
  INTEGER :: i, n

    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeEMS
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    IF ((rd%NEQ .EQ. 0) .OR. (p_rmatrix%NEQ .EQ. 0) .OR. &
        (p_rmatrix%NEQ .NE. rd%NEQ) ) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF
    
    ! Get our temporary vector
    p_rt => p_rsubnode%rtempVector
    
!    IF (p_rmatrix%imatrixSpec .EQ. LSYSBS_MSPEC_SADDLEPOINT) THEN
!    
!      n = p_rmatrix%ndiagBlocks
!      
!      CALL lsyssc_clearVector(p_rt%RvectorBlock(n))
!    
!      ! Go through all diagonal blocks
!      DO i = 1, p_rmatrix%ndiagBlocks-1
!        ! t_i := B_i * d_n
!        CALL lsyssc_scalarMatVec(p_rmatrix%RmatrixBlock(i,n),&
!          rd%RvectorBlock(n), p_rt%RvectorBlock(i), 1.0_DP, 0.0_DP)
!        
!        ! t_n := t_n - D_i * t_i
!        CALL lsyssc_scalarMatVec(p_rmatrix%RmatrixBlock(n,i),&
!          p_rt%RvectorBlock(i), p_rt%RvectorBlock(n), -1.0_DP, 1.0_DP)
!          
!        ! t_i := A_i * d_i
!        CALL lsyssc_scalarMatVec(p_rmatrix%RmatrixBlock(i,i),&
!          rd%RvectorBlock(i), p_rt%RvectorBlock(i), 1.0_DP, 0.0_DP)
!        !CALL lsyssc_clearVector(p_rt%RvectorBlock(i))
!          
!      END DO
!      
!      !CALL lsysbl_blockMatVec(p_rmatrix, rd, p_rt, 1.0_DP, 1.0_DP)
!
!    ELSE
    
      ! Calculate t := A*d
      CALL lsysbl_blockMatVec(p_rmatrix, rd, p_rt, 1.0_DP, 0.0_DP)
    
!    END IF
      
    ! Calculate alpha := < d, t >
    dalpha = lsysbl_scalarProduct(rd, p_rt)
      
    ! Calculate beta := < t, t >
    dbeta = lsysbl_scalarProduct(p_rt, p_rt)
    
    ! Check if alpha and beta are valid
    IF (dalpha .EQ. 0.0_DP) THEN
      CALL output_line('EMS: < d, A*d > = 0 !!!')
      rsolverNode%iresult = 1
    END IF
    IF (dbeta .EQ. 0.0_DP) THEN
      CALL output_line('EMS: < A*d, A*d > = 0 !!!')
      rsolverNode%iresult = 1
    END IF
    
    ! Okay, calculate omega := alpha / beta
    domega = dalpha / dbeta
    
    ! And scale defect by omega
    CALL lsysbl_scaleVector(rd, domega)
    
    ! That's it
  
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
  
  RECURSIVE SUBROUTINE linsol_matCompatJacobi (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
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
    REAL(DP) :: dlocOmega,dscale
    REAL(SP) :: flocOmega,fscale
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
          CALL lsyssc_getbase_Kdiagonal(p_rmatrix,p_Kdiag)
        ELSE
          CALL lsyssc_getbase_Kld(p_rmatrix,p_Kdiag)
        END IF
      
        ! Get the omega parameter for the matrix. Don't forget the scaling 
        ! factor in the matrix!
        dlocomega = rsolverNode%domega * p_rmatrix%dScaleFactor

        ! Now, which data format do we have? Single or double?
        SELECT CASE (p_rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Get the matrix data arrays
          CALL lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            IF (p_rmatrix%dscaleFactor .EQ. 1.0_DP) THEN
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = dlocOmega * p_Dvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
              END DO
            ELSE
              dscale = 1.0_DP/p_rmatrix%dscaleFactor
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = dscale * dlocOmega * p_Dvector(ieq) &
                                 / p_Dmatrix(p_Kdiag(ieq))
              END DO
            END IF
            
          CASE (ST_SINGLE)
            ! Get the data array
            CALL lsysbl_getbase_single (rd,p_Fvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            IF (p_rmatrix%dscaleFactor .EQ. 1.0_DP) THEN
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = dlocOmega * p_Fvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
              END DO
            ELSE
              dscale = 1.0_DP/p_rmatrix%dscaleFactor
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = dscale * dlocOmega * p_Fvector(ieq) &
                                 / p_Dmatrix(p_Kdiag(ieq))
              END DO
            END IF

          CASE DEFAULT
            PRINT *,'Jacobi: Unsupported vector format.'
            CALL sys_halt()
          END SELECT
          
        CASE (ST_SINGLE)

          ! Get the matrix data arrays
          CALL lsyssc_getbase_single (p_rmatrix,p_Fmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            IF (p_rmatrix%dscaleFactor .EQ. 1.0_DP) THEN
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = dlocOmega * p_Dvector(ieq) / p_Fmatrix(p_Kdiag(ieq)) 
              END DO
            ELSE
              dscale = 1.0_DP/p_rmatrix%dscaleFactor
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = dscale * dlocOmega * p_Dvector(ieq) &
                                 / p_Fmatrix(p_Kdiag(ieq)) 
              END DO
            END IF
            
          CASE (ST_SINGLE)
            ! Get the data array
            CALL lsysbl_getbase_single (rd,p_Fvector)
            
            ! Multiplication with Omega can be speeded up as we use
            ! sigle-precision only.
            flocOmega = REAL(dlocOmega,SP)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            IF (p_rmatrix%dscaleFactor .EQ. 1.0_DP) THEN
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = flocOmega * p_Fvector(ieq) / p_Fmatrix(p_Kdiag(ieq))
              END DO
            ELSE
              fscale = 1.0_DP/p_rmatrix%dscaleFactor
              DO ieq = 1,rd%NEQ
                p_Dvector(ieq) = fscale * flocOmega * p_Fvector(ieq) &
                                 / p_Fmatrix(p_Kdiag(ieq))
              END DO
            END IF

          CASE DEFAULT
            PRINT *,'Jacobi: Unsupported vector format.'
            CALL sys_halt()
          END SELECT

        CASE DEFAULT
          PRINT *,'Jacobi: Unsupported matrix format.'
          CALL sys_halt()
        END SELECT
      
      CASE DEFAULT
        PRINT *,'Jacobi: Unsupported matrix format.'
        CALL sys_halt()
      END SELECT
      
    END DO
  
  END SUBROUTINE

! *****************************************************************************
! Routines for the Jin-Wei-Tam solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initJinWeiTam (p_rsolverNode, drelax, domega)
  
!<description>
  ! Creates a t_linsolNode solver structure for the Jin-Wei-Tam solver. The node
  ! can be used to directly solve a problem or to be attached as a solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>

!<input>
  ! OPTIONAL: Relax parameter. Is saved to rsolverNode%p_rsubnodeJinWeiTam%drelax 
  !           if specified. 0.0=Standard Jacobi, 1.0=full JinWeiTam
  REAL(DP), OPTIONAL :: drelax
!</input>

!<input>
  ! OPTIONAL: Damping parameter. Is saved to rsolverNode%domega if specified.
  REAL(DP), OPTIONAL :: domega
!</input>
 

!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
!</output>

    ! Create a default solver structure
    
    CALL linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_JINWEITAM 
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
    
    ! No subnode for Jin-Wei-Tam. Only save domega to the structure.
    IF (PRESENT(domega)) THEN
      p_rsolverNode%domega = domega
    END IF
    
    ! Allocate the subnode for Jin-Wei-Tam
    ALLOCATE(p_rsolverNode%p_rsubnodeJinWeiTam)
    
    ! If a relax parameter has been passed, save it,
    ! otherwise save 1.0_DP as relax parameter
    IF (PRESENT(drelax)) THEN
      p_rsolverNode%p_rsubnodeJinWeiTam%drelax = drelax
    ELSE
      p_rsolverNode%p_rsubnodeJinWeiTam%drelax = 1.0_DP
    END IF
    
    ! Initialise the matrix sum to 0.0_DP
    p_rsolverNode%p_rsubnodeJinWeiTam%dmatrixSum = 0.0_DP
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatJinWeiTam (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! We cannot handle the matrix if...
    !
    ! ... it's not scalar
    IF (Rmatrices(UBOUND(Rmatrices,1))%ndiagBlocks .NE. 1) THEN
      ccompatible = LINSOL_COMP_ERRNOTSCALAR
    END IF

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneJinWeiTam (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the Jin-Wei-Tam solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of Jin-Wei-Tam which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeJinWeiTam)
  
  END SUBROUTINE

!<subroutine>

  ! ***************************************************************************
  
  RECURSIVE SUBROUTINE linsol_initDataJinWeiTam (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! This routine initialises the data for the Jin-Wei-Tam solver.
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Jin-Wei-Tam solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET       :: rsolverNode
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
    INTEGER :: iblockrow, iblockcol, ientry
    REAL(DP) :: dsum
    REAL(DP), DIMENSION(:), POINTER :: p_Dmatrix
    REAL(SP), DIMENSION(:), POINTER :: p_Fmatrix

    ! a pointer to the currently selected scalar matrix
    TYPE (t_matrixScalar), POINTER :: p_rmatrix

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup
    
    ! format the total sum to 0
    dSum = 0.0_DP

    ! go through all block rows of our matrix
    DO iblockrow=1, rsolverNode%rsystemMatrix%ndiagBlocks
    
      ! go through all blocks in our block row
      DO iblockcol=1, rsolverNode%rsystemMatrix%ndiagBlocks
        
        ! store our current scalar matrix
        p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblockrow,iblockcol)

        ! skip the current block if it is empty
        IF(p_rmatrix%NA .EQ. 0) THEN
          CYCLE
        END IF
        
        ! Check the scale factor
        !IF(p_rmatrix%dscaleFactor .NE. 1.0_DP) THEN
        !  PRINT *,'Jin-Wei-Tam: Unsupported scale factor'
        !  ierror = LINSOL_ERR_INITERROR
        !  CALL sys_halt()
        !END IF

        ! check the data format
        SELECT CASE(p_rmatrix%cdataType)
        CASE(ST_DOUBLE)
          ! now get the matrix' data
          CALL lsyssc_getbase_double (p_rmatrix, p_Dmatrix)
        
          ! go through all data entries of the matrix and sum them up
          DO ientry=1, p_rmatrix%NA
            dSum = dSum + p_Dmatrix(ientry)
          END DO
          
        CASE(ST_SINGLE)
          ! now get the matrix' data
          CALL lsyssc_getbase_single(p_rmatrix, p_Fmatrix)

          ! go through all data entries of the matrix and sum them up
          DO ientry=1, p_rmatrix%NA
            dSum = dSum + p_Fmatrix(ientry)
          END DO
          
        CASE DEFAULT
          PRINT *,'Jin-Wei-Tam: Unsupported matrix format'
          ierror = LINSOL_ERR_INITERROR
        
        END SELECT
          
      END DO
      
    END DO
    
    ! Make sure that the sum of all matrix elements is not 0
    IF(dSum .EQ. 0.0_DP) THEN
      PRINT *,'Jin-Wei-Tam: Sum of all matrix elements is 0.0!'
      ierror = LINSOL_ERR_INITERROR
      CALL sys_halt()
    END IF
    
    ! everything went fine, we can store the sum now
    rsolverNode%p_rsubnodeJinWeiTam%dmatrixSum = dSum
    
  END SUBROUTINE
  
  ! ***************************************************************************

  RECURSIVE SUBROUTINE linsol_precJinWeiTam (rsolverNode,rd)
  
!<description>
  ! Applies the Jin-Wei-Tam preconditioner $P$ to the defect, where
  !  $$ P :=  (\alpha * [1, ..., 1]^T * [1, ..., 1]) + D^{-1} $$
  ! where $D$ is the main diagonal of $A$ and
  !  $$ \alpha := \sum_{1 \leq i,j \leq n} A_{ij} $$
  !
  ! Instead of performing the matrix-vector-multiplication 
  !  $$ d_{new} := Pd $$
  ! which would cost O(n^2) arithmetic operations, we use another approach:
  !  $$ d_{new} := [x, x, ..., x]^T + D^{-1}d $$
  ! where
  !  $$ x := (\sum_{1 \leq i \leq n}d_i) * \alpha^{-1} $$
  ! Using this method, we only use O(n) operations to calculate the new
  ! preconditioned defect vector $d_{new}$
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initData routine must have been called to 
  ! prepare the solver for solving the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Jin-Wei-Tam solver
  TYPE(t_linsolNode), INTENT(INOUT),TARGET  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER iblock
    INTEGER(PREC_VECIDX) :: ieq
    
    TYPE (t_matrixScalar), POINTER :: p_rmatrix
    TYPE (t_vectorScalar), POINTER :: p_rvector
    INTEGER (PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiag
    REAL(DP) :: dlocOmega
    REAL(SP) :: flocOmega
    REAL(DP), DIMENSION(:), POINTER :: p_Dvector, p_Dmatrix
    REAL(SP), DIMENSION(:), POINTER :: p_Fvector, p_Fmatrix
    REAL(DP) :: dvecSum
    
    ! initialize vector sum to 0
    dvecSum = 0.0_DP
    
    ! first, we need to calculate the sum of all vector elements
    ! $$ dvecSum := \sum_{1 \leq i \leq n}d_i $$
    ! For this, we need to go through all blocks of our vector
    DO iblock = 1,rd%nblocks
      ! Get the vector block
      p_rvector => rd%RvectorBlock(iblock)

      ! check the vector's precision
      SELECT CASE (p_rvector%cdataType)
      CASE(ST_DOUBLE)
        ! get the vector's data array
        CALL storage_getbase_double(p_rvector%h_Ddata, p_Dvector)
      
        ! go through all vector elements
        DO ieq=1, rd%NEQ
          dvecSum = dvecSum + p_Dvector(ieq)
        END DO
        
      CASE(ST_SINGLE)
        ! get the vector's data array
        CALL storage_getbase_single(p_rvector%h_Ddata, p_Fvector)
      
        ! go through all vector elements
        DO ieq=1, rd%NEQ
          dvecSum = dvecSum + p_Fvector(ieq)
        END DO

      CASE DEFAULT      
        PRINT *,'Jin-Wei-Tam: Unsupported Vector format'
        CALL sys_halt()
        
      END SELECT
      
    END DO
    
    ! now, divide the vector sum through the matrix sum
    dvecSum = dvecSum / rsolverNode%p_rsubnodeJinWeiTam%dmatrixSum
    
    ! and apply our relax parameter (usually 1.0_DP)
    dvecSum = dvecSum * rsolverNode%p_rsubnodeJinWeiTam%drelax
    
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
          CALL lsyssc_getbase_Kdiagonal(p_rmatrix,p_Kdiag)
        ELSE
          CALL lsyssc_getbase_Kld(p_rmatrix,p_Kdiag)
        END IF
      
        ! Get the omega parameter for the matrix. Don't forget the scaling 
        ! factor in the matrix!
        dlocomega = rsolverNode%domega * p_rmatrix%dScaleFactor

        ! Now, which data format do we have? Single or double?
        SELECT CASE (p_rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Get the matrix data arrays
          CALL lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            DO ieq = 1,rd%NEQ
              p_Dvector(ieq) = dlocOmega * &
                  (dvecSum + (p_Dvector(ieq) / p_Dmatrix(p_Kdiag(ieq))))
            END DO
            
          CASE (ST_SINGLE)
            ! Get the data array
            CALL lsysbl_getbase_single (rd,p_Fvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            DO ieq = 1,rd%NEQ
              p_Dvector(ieq) = dlocOmega * &
                  (dvecSum + (p_Fvector(ieq) / p_Dmatrix(p_Kdiag(ieq))))
            END DO

          CASE DEFAULT
            PRINT *,'Jin-Wei-Tam: Unsupported vector format.'
            CALL sys_halt()
          END SELECT
          
        CASE (ST_SINGLE)

          ! Get the matrix data arrays
          CALL lsyssc_getbase_single (p_rmatrix,p_Fmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! and multiply all entries with the inverse of the diagonal 
            ! of the matrix.
            DO ieq = 1,rd%NEQ
              p_Dvector(ieq) = dlocOmega * &
                  (dvecSum + (p_Dvector(ieq) / p_Fmatrix(p_Kdiag(ieq))))
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
              p_Dvector(ieq) = flocOmega * &
                  (dvecSum + (p_Fvector(ieq) / p_Fmatrix(p_Kdiag(ieq))))
            END DO

          CASE DEFAULT
            PRINT *,'Jin-Wei-Tam: Unsupported vector format.'
            CALL sys_halt()
          END SELECT

        CASE DEFAULT
          PRINT *,'Jin-Wei-Tam: Unsupported matrix format.'
          CALL sys_halt()
        END SELECT
      
      CASE DEFAULT
        PRINT *,'Jin-Wei-Tam: Unsupported matrix format.'
        CALL sys_halt()
      END SELECT
      
    END DO

  END SUBROUTINE
  
! *****************************************************************************
! Routines for the SOR solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initSOR (p_rsolverNode, domega)
  
!<description>
  ! Creates a t_linsolNode solver structure for the SOR solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  !
  ! This SOR solver has no done-routine as there is no dynamic information
  ! allocated. The given damping parameter domega is saved to the solver
  ! structure rsolverNode%domega.
  !
  ! For domega = 1.0_DP this SOR solver is identical to the Gauss-Seidel
  ! algorithm.
!</description>
  
!<input>
  ! OPTIONAL: Damping parameter. Is saved to rsolverNode%domega if specified.
  REAL(DP), INTENT(IN), OPTIONAL :: domega

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
    p_rsolverNode%calgorithm = LINSOL_ALG_SOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
    
    ! Save domega to the structure.
    IF (PRESENT(domega)) THEN
      p_rsolverNode%domega = domega
    END IF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatSOR (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precSOR (rsolverNode,rd)
  
!<description>
  ! Applies the SOR preconditioner $D \in A$ to the defect 
  ! vector rd and solves 
  !  $$ (D + \omega L) d_{new} = d.$$
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the SOR solver
  TYPE(t_linsolNode), INTENT(INOUT),TARGET  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: iblock
    
    TYPE (t_matrixScalar), POINTER :: p_rmatrix
    REAL(DP), DIMENSION(:), POINTER :: p_Dvector, p_Dmatrix
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
   
    ! Loop through all diagonal blocks. Each block corresponds to one
    ! diagonal block in the matrix.
    DO iblock = 1,rd%nblocks
      ! Get the matrix
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
      
      ! Some small checks...
      IF (p_rmatrix%NEQ .EQ. 0) THEN
        PRINT *,'SOR: No diagonal submatrix for component ',iblock
        CALL sys_halt()
      END IF
      
      IF (IAND(p_rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
          .NE. 0) THEN
        PRINT *,'SOR: Transposed submatrices not supported.'
        CALL sys_halt()
      END IF

      ! Now we have to make some decisions. At first, which matrix
      ! structure do we have?
      SELECT CASE (p_rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)

        CALL lsyssc_getbase_Kdiagonal(p_rmatrix,p_Kdiagonal)
        CALL lsyssc_getbase_Kcol(p_rmatrix,p_Kcol)
        CALL lsyssc_getbase_Kld(p_rmatrix,p_Kld)
      
        ! Now, which data format do we have? Single or double?
        SELECT CASE (p_rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Get the matrix data arrays
          CALL lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! Call the SOR subroutine (see below), do the work.
            CALL performSOR9dbledble_ID119 (p_Dmatrix,p_Kcol,p_Kld,&
                                             p_Kdiagonal,rsolverNode%domega,&
                                             p_Dvector,p_rmatrix%dscaleFactor)
            
          CASE DEFAULT
            PRINT *,'SOR: Unsupported vector format.'
            CALL sys_halt()
          END SELECT
          
        CASE DEFAULT
          PRINT *,'SOR: Unsupported matrix format.'
          CALL sys_halt()
        END SELECT
      
      CASE DEFAULT
        PRINT *,'SOR: Unsupported matrix format.'
        CALL sys_halt()
      END SELECT
      
    END DO
    
  CONTAINS
  
    !--------------------------------------------------------------------------
    ! Auxiliary routine: SOR
    ! Matrix format 9, double precision matrix, double precision vector
    
    SUBROUTINE performSOR9dbledble_ID119 (DA,Kcol,Kld,Kdiagonal,domega,Dx,dscale)
    
    ! input: Matrix array
    REAL(DP), DIMENSION(:), INTENT(IN) :: DA
    
    ! input: column structure
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol
    
    ! input: row structure
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld
    
    ! input: position of diagonal entries in the matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kdiagonal
    
    ! input: Relaxation parameter; standard value is 1.2.
    REAL(DP), INTENT(IN) :: domega
    
    ! Scaling factor of the matrix; usually = 1.0
    REAL(DP), INTENT(IN) :: dscale
    
    ! input: vector to be preconditioned.
    ! output: preconditioned vector
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dx
    
      ! local variables
      INTEGER(PREC_VECIDX) :: NEQ,ieq
      INTEGER(PREC_MATIDX) :: Idiag,ICOL
      REAL(DP) :: daux,dsc
      
      NEQ = SIZE(Dx)

      ! We perform the following preconditioning step:
      !
      !    (D+\omega L) d_{new} = d
      !
      ! and replace the vector Dx in-situ.
      !
      ! We're doing it this way:
      !
      !          (D+\omega L) d = d
      !
      !                          i-1
      ! <=> a_ii d_i  +  \omega \sum a_ij d_j  =  d_i          (i=1,...,n)
      !                          j=1
      !
      !                                       i-1
      ! <=> d_i :=  1/aii * ( d_i  -  \omega \sum a_ij d_j )   (i=1,...,n)
      !                                       j=1
      !
      ! The scaling factor in introduced into that formula by multiplying
      ! all a_ij by dscale -- which reduces to dividing d_i by dscale
      ! in every step...
      
      dsc = 1.0_DP/dscale  
      
      DO ieq=1,NEQ
        daux = 0.0_DP
        ! Loop through the L; this is given by all
        ! entries below Kdiagonal.
        Idiag = Kdiagonal(ieq)
        DO ICOL=Kld(ieq),Kdiagonal(ieq)-1
          daux=daux+DA(ICOL)*Dx(Kcol(ICOL))
        END DO
        Dx(ieq) = (dsc*Dx(ieq) - daux*domega) / DA(Idiag)
      END DO
 
    END SUBROUTINE
  
  END SUBROUTINE

! *****************************************************************************
! Routines for the SSOR solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initSSOR (p_rsolverNode, domega, bscale)
  
!<description>
  ! Creates a t_linsolNode solver structure for the SSOR solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! OPTIONAL: Damping parameter. Is saved to rsolverNode%domega if specified.
  REAL(DP), INTENT(IN), OPTIONAL :: domega

  ! OPTIONAL: If set to TRUE, the solution is scaled by 1/(domega*(2-omega))
  ! which gives the SSOR preconditioner in the literature. If not existent
  ! or set to FALSE, no scaling is performed; this is the original 
  ! implementation of SSOR in FEAT.
  LOGICAL, INTENT(IN), OPTIONAL :: bscale
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
    p_rsolverNode%calgorithm = LINSOL_ALG_SSOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
    
    ! Save domega to the structure.
    IF (PRESENT(domega)) THEN
      p_rsolverNode%domega = domega
    END IF
  
    ! Allocate the subnode for SSOR.
    ALLOCATE(p_rsolverNode%p_rsubnodeSSOR)

    ! Save whether the solution should be scaled or not.
    p_rsolverNode%p_rsubnodeSSOR%bscale = .FALSE.
    IF (PRESENT(bscale)) p_rsolverNode%p_rsubnodeSSOR%bscale = bscale
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatSSOR (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible
  
  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneSSOR (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the SSOR solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of SSOR which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  DEALLOCATE(rsolverNode%p_rsubnodeSSOR)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precSSOR (rsolverNode,rd)
  
!<description>
  ! Applies the SSOR preconditioner $D \in A$ to the defect 
  ! vector rd and solves 
  !  $$ \frac{1}{\omega(2-\omega)} (D+\omega L)D^{-1}(D+\omega R) d_{new} = d.$$
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the SSOR solver
  TYPE(t_linsolNode), INTENT(INOUT),TARGET  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: iblock
    LOGICAL :: bscale
    
    TYPE (t_matrixScalar), POINTER :: p_rmatrix
    REAL(DP), DIMENSION(:), POINTER :: p_Dvector, p_Dmatrix
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
    
    ! Get bscale from the solver structure
    bscale = rsolverNode%p_rsubnodeSSOR%bscale
    
    ! Loop through all diagonal blocks. Each block corresponds to one
    ! diagonal block in the matrix.
    DO iblock = 1,rd%nblocks
      ! Get the matrix
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
      
      ! Some small checks...
      IF (p_rmatrix%NEQ .EQ. 0) THEN
        PRINT *,'SSOR: No diagonal submatrix for component ',iblock
        CALL sys_halt()
      END IF
      
      IF (IAND(p_rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
          .NE. 0) THEN
        PRINT *,'SSOR: Transposed submatrices not supported.'
        CALL sys_halt()
      END IF

      ! Now we have to make some decisions. At first, which matrix
      ! structure do we have?
      SELECT CASE (p_rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)

        CALL lsyssc_getbase_Kdiagonal(p_rmatrix,p_Kdiagonal)
        CALL lsyssc_getbase_Kcol(p_rmatrix,p_Kcol)
        CALL lsyssc_getbase_Kld(p_rmatrix,p_Kld)
      
        ! Now, which data format do we have? Single or double?
        SELECT CASE (p_rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Get the matrix data arrays
          CALL lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          SELECT CASE (rd%cdataType)
          CASE (ST_DOUBLE)
            ! Get the data array
            CALL lsysbl_getbase_double (rd,p_Dvector)
            
            ! Call the SSOR subroutine (see below), do the work.
            CALL performSSOR9dbledble_ID119 (p_Dmatrix,p_Kcol,p_Kld,&
                                             p_Kdiagonal,rsolverNode%domega,&
                                             bscale,p_Dvector,&
                                             p_rmatrix%dscaleFactor)
            
          CASE DEFAULT
            PRINT *,'SSOR: Unsupported vector format.'
            CALL sys_halt()
          END SELECT
          
        CASE DEFAULT
          PRINT *,'SSOR: Unsupported matrix format.'
          CALL sys_halt()
        END SELECT
      
      CASE DEFAULT
        PRINT *,'SSOR: Unsupported matrix format.'
        CALL sys_halt()
      END SELECT
      
    END DO
    
  CONTAINS
  
    !--------------------------------------------------------------------------
    ! Auxiliary routine: SSOR
    ! Matrix format 9, double precision matrix, double precision vector
    
    SUBROUTINE performSSOR9dbledble_ID119 (DA,Kcol,Kld,Kdiagonal,domega,bscale,&
                                           Dx,dscale)
    
    ! input: Matrix array
    REAL(DP), DIMENSION(:), INTENT(IN) :: DA
    
    ! input: column structure
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol
    
    ! input: row structure
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld
    
    ! input: position of diagonal entries in the matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kdiagonal
    
    ! input: Relaxation parameter; standard value is 1.2.
    REAL(DP), INTENT(IN) :: domega

    ! input: Scaling factor of the matrix; usually = 1.0
    REAL(DP), INTENT(IN) :: dscale
    
    ! input: Whether the solution should be scaled as suggested in 
    ! the literature
    LOGICAL, INTENT(IN) :: bscale
    
    ! input: vector to be preconditioned.
    ! output: preconditioned vector
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dx
    
      ! local variables
      INTEGER(PREC_VECIDX) :: NEQ,ieq
      INTEGER(PREC_MATIDX) :: Idiag,ICOL
      REAL(DP) :: daux,dsc
      
      NEQ = SIZE(Dx)

      ! We perform the following preconditioning step:
      !
      !    (D+\omega L) D^{-1} (D+\omega R) d_{new} = d
      !
      ! and replace the vector Dx in-situ.
      !
      ! At first we solve:
      !
      !          (D+\omega L) d = d
      !
      !                          i-1
      ! <=> a_ii d_i  +  \omega \sum a_ij d_j  =  d_i          (i=1,...,n)
      !                          j=1
      !
      !                                       i-1
      ! <=> d_i :=  1/aii * ( d_i  -  \omega \sum a_ij d_j )   (i=1,...,n)
      !                                       j=1
      !  
      ! The scaling factor in introduced into that formula by multiplying
      ! all a_ij by dscale -- which reduces to dividing d_i by dscale
      ! in every step...
      
      dsc = 1.0_DP/dscale  
      
      DO ieq=1,NEQ
        daux = 0.0_DP
        ! Loop through the L; this is given by all
        ! entries below Kdiagonal.
        Idiag = Kdiagonal(ieq)
        DO ICOL=Kld(ieq),Kdiagonal(ieq)-1
          daux=daux+DA(ICOL)*Dx(Kcol(ICOL))
        END DO
        Dx(ieq) = (dsc*Dx(ieq) - daux*domega) / DA(Idiag)
      END DO
      
      ! Next step is to solve
      !
      !          D^{-1} (D+\omega R) d  =  d
      !
      ! <=>             (D+\omega R) d  =  D d
      !
      !                           n
      ! <=> a_ii d_i  +  \omega \sum a_ij d_j  =  a_ii d_i          (i=n,...,1)
      !                         j=i+1
      !
      !                                n
      ! <=> d_i  +  (\omega / a_ii)  \sum a_ij d_j  =  d_i          (i=n,...,1)
      !                              j=i+1
      !
      !                                        n
      ! <=> d_i :=  d_i  -  (\omega / a_ii)  \sum a_ij d_j          (i=n,...,1)
      !                                      j=i+1
      !
      ! Note that the case i=n is trivial, the sum is empty; so we can
      ! start the loop with i=n-1.
      !
      ! The scaling factor cancels out in this equation, so nothing has
      ! to be done here concerning that.

      DO ieq = NEQ-1,1,-1
        daux=0.0_DP
        ! Loop through the L; this is given by all
        ! entries above Kdiagonal.
        Idiag = Kdiagonal(ieq)
        DO ICOL=Idiag+1,Kld(ieq+1)-1
          daux=daux+DA(ICOL)*Dx(Kcol(ICOL))
        END DO
        Dx(ieq) = Dx(ieq)-daux*domega / DA(Idiag)
      END DO
      
      ! The literature suggests the formula
      !
      !   1/(omega (2-omega))  (D+\omega L) D^{-1} (D+\omega R) d_{new} = d
      !
      ! See e.g.: 
      !  [Barrett et al., "Templates for the Solution of Linear systems:
      !   Building Blocks for Iterative Methods", p. 42,
      !   http://www.netlib.org/linalg/html_templates/Templates.html]        
      !
      ! We have calculated only
      !                        (D+\omega L) D^{-1} (D+\omega R) d_{new} = d
      !
      ! This gives interestingly better results when using SSOR as 
      ! preconditioner in CG or BiCGStab and corresponds to the old FEAT
      ! implementation.
      !
      ! By setting bscale=TRUE, we activate the formula suggested in the
      ! literature, which yields another scaling of our vector d by
      ! the scalar omega*(2-omega):
      
      IF (bscale) THEN
        daux = domega * (2.0_DP-domega)
        IF (daux .NE. 1.0_DP) THEN
          CALL lalg_scaleVectorDble (Dx,daux)
        END IF
      END IF
 
    END SUBROUTINE
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the VANCA CC2D/CC3D solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initVANCA (p_rsolverNode,domega,csubtypeVANCA)
  
!<description>
  ! Creates a t_linsolNode solver structure for the VANCA solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  !
  ! VANCA is somehow a special type of solver. There is one general VANCA
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
  ! OPTIONAL: Damping parameter. Is saved to rsolverNode\%domega if specified.
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
  
  RECURSIVE SUBROUTINE linsol_alterVANCA (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the VANCA solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration 
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check the command.
    IF (ralterConfig%ccommand .EQ. LINSOL_ALTER_CHANGEVANCA) THEN
    
      ! That's a command to all VANCA solvers. Are we meant?
      ! Check Iconfig(1) which specifies the VANCA subsolver group
      ! to be changed.
      IF (ralterConfig%Iconfig(1) .EQ. rsolverNode%p_rsubnodeVANCA%rvanca%csubtype) THEN
      
        ! Oops, we are meant. Set the VANCA variant as specified in Iconfig(2).
        rsolverNode%p_rsubnodeVANCA%rvanca%csubtype = ralterConfig%Iconfig(2)
      
      END IF
    
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatVANCA (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN),TARGET   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    INTEGER :: iblock,jblock
    TYPE(t_matrixBlock), POINTER :: p_rmat

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! VANCA is a bit tricky. Loop through all scalar submatrices and
    ! check them. The VANCA subtype decides on whether it can use that or not!
    
    p_rmat => Rmatrices(UBOUND(Rmatrices,1))
    ! VANCA subtype?
    SELECT CASE (rsolverNode%p_rsubnodeVANCA%csubtypeVANCA)
    CASE (LINSOL_VANCA_GENERAL,LINSOL_VANCA_GENERALDIRECT)
      
      ! Check all sub-matrices
      DO jblock = 1,p_rmat%ndiagblocks
        DO iblock = 1,p_rmat%ndiagblocks
          IF (p_rmat%RmatrixBlock(iblock,jblock)%NEQ .NE. 0) THEN
            IF (IAND(p_rmat%RmatrixBlock(iblock,jblock)%imatrixSpec,&
                LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
              ! We can't handle transposed matrices.
              ccompatible = LINSOL_COMP_ERRTRANSPOSED
            END IF
          END IF
        END DO
      END DO

    CASE (LINSOL_VANCA_2DNAVST  , LINSOL_VANCA_2DNAVSTDIRECT  ,&
          LINSOL_VANCA_2DFNAVST , LINSOL_VANCA_2DFNAVSTDIRECT )
      ! Blocks (3,1) and (3,2) must be virtually transposed
      IF ((IAND(p_rmat%RmatrixBlock(3,1)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) .OR. &
          (IAND(p_rmat%RmatrixBlock(3,2)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .EQ. 0)) THEN
        ccompatible = LINSOL_COMP_ERRNOTTRANSPOSED
      END IF

    CASE (LINSOL_VANCA_3DNAVST  , LINSOL_VANCA_3DNAVSTDIRECT  ,&
          LINSOL_VANCA_3DFNAVST , LINSOL_VANCA_3DFNAVSTDIRECT )
      ! Blocks (4,1), (4,2) and (4,3) must be virtually transposed
      IF ((IAND(p_rmat%RmatrixBlock(4,1)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) .OR. &
          (IAND(p_rmat%RmatrixBlock(4,2)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) .OR. &
          (IAND(p_rmat%RmatrixBlock(4,3)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .EQ. 0)) THEN
        ccompatible = LINSOL_COMP_ERRNOTTRANSPOSED
      END IF
      
    END SELECT
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
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
    CALL sys_halt()
  END IF
  
  ! If isubgroup does not coincide with isolverSubgroup from the solver
  ! structure, skip the rest here.
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

  ! Allocate a temporary vector
  CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
        rsolverNode%p_rsubnodeVANCA%rtempVector,.FALSE.,.FALSE.,&
        rsolverNode%cdefaultDataType)
        
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
      CALL sys_halt()
    END IF
    
    ! Which VANCA solver do we actually have?
    SELECT CASE (rsolverNode%p_rsubnodeVANCA%csubtypeVANCA)
    CASE (LINSOL_VANCA_GENERAL,LINSOL_VANCA_GENERALDIRECT)
      ! General VANCA for everything
      CALL vanca_initConformal(rsolverNode%rsystemMatrix,&
                              rsolverNode%p_rsubnodeVANCA%rvanca,&
                              VANCAPC_GENERAL,VANCATP_STANDARD)
                               
    CASE (LINSOL_VANCA_2DNAVST  , LINSOL_VANCA_2DNAVSTDIRECT  )
      ! Diagonal-type VANCA for Navier-Stokes
      CALL vanca_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANCA%rvanca,&
                                VANCAPC_2DNAVIERSTOKES,VANCATP_DIAGONAL)
                                
    CASE (LINSOL_VANCA_2DFNAVST , LINSOL_VANCA_2DFNAVSTDIRECT )
      ! Full VANCA for Navier-Stokes
      CALL vanca_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANCA%rvanca,&
                                VANCAPC_2DNAVIERSTOKES,VANCATP_FULL)

    CASE (LINSOL_VANCA_2DFNAVSTOC,LINSOL_VANCA_2DFNAVSTOCDIRECT )
      ! Full VANCA for Navier-Stokes optimal control
      CALL vanca_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANCA%rvanca,&
                                VANCAPC_2DNAVIERSTOKESOPTC,VANCATP_FULL)

    CASE (LINSOL_VANCA_2DFNAVSTOCDIAG,LINSOL_VANCA_2DFNAVSTOCDIAGDIR )
      ! Full VANCA for Navier-Stokes optimal control
      CALL vanca_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANCA%rvanca,&
                                VANCAPC_2DNAVIERSTOKESOPTC,VANCATP_DIAGOPTC)

    CASE (LINSOL_VANCA_3DNAVST  , LINSOL_VANCA_3DNAVSTDIRECT  )
      ! Diagonal-type VANCA for Navier-Stokes
      CALL vanca_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANCA%rvanca,&
                                VANCAPC_3DNAVIERSTOKES,VANCATP_DIAGONAL)
                                
    CASE (LINSOL_VANCA_3DFNAVST , LINSOL_VANCA_3DFNAVSTDIRECT )
      ! Full VANCA for Navier-Stokes
      CALL vanca_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANCA%rvanca,&
                                VANCAPC_3DNAVIERSTOKES,VANCATP_FULL)
                                
    END SELECT
      
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataVANCA (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases temporary memory of the VANCA solver allocated in 
  ! linsol_initDataVANCA
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
    
    ! Release VANCA.
    CALL vanca_doneConformal (rsolverNode%p_rsubnodeVANCA%rvanca)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureVANCA (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases temporary memory of the VANCA solver allocated in
  ! linsol_initStrutureVANCA.
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
  CALL linsol_doneDataVANCA (rsolverNode, isubgroup)
  
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
    
    ! Execute VANCA
    CALL vanca_conformal (rsolverNode%p_rsubnodeVANCA%rvanca, &
        p_rvector, rd, domega)

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
  
  RECURSIVE SUBROUTINE linsol_matCompatUMFPACK4 (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible

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
  INTEGER(I32) :: idupFlag

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
    CALL sys_halt()
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
    
    ! We can now modify p_rmatrix without harming the original matrix,
    ! so set the dup-flag for the copy command to rtempMatrix below
    ! to LSYSSC_DUP_SHARE.
    idupFlag = LSYSSC_DUP_SHARE
  
  ELSE
    p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock (1,1)

    ! As we work with the original matrix, we set idupFlag
    ! to LSYSSC_DUP_COPY. This prevents the original matrix from being
    ! destroyed in the copy-modify-part below.
    idupFlag = LSYSSC_DUP_COPY
  END IF

  SELECT CASE (p_rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
    ! Format 9 is exactly the UMFPACK matrix.
    ! Make a copy of the matrix so we can change its content for UMF4SYM.
    CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                 idupFlag,idupFlag)
  CASE (LSYSSC_MATRIX7)
    ! For format 7, we have to modify the matrix slightly.
    ! Make a copy of the matrix:
    CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                 idupFlag,idupFlag)
    ! Resort the entries to put the diagonal entry to the correct position.
    ! This means: Convert the structure-7 matrix to a structure-9 matrix:
    CALL lsyssc_convertMatrix (rtempMatrix,LSYSSC_MATRIX9)
  END SELECT
  
  IF (p_rmatrix%cdataType .NE. ST_DOUBLE) THEN
    PRINT *,'UMFPACK can only handle double precision matrices!'
    CALL sys_halt()
  END IF

  ! Modify Kcol/Kld of the matrix. Subtract 1 to get the 0-based.
  CALL lsyssc_addIndex (rtempMatrix%h_Kcol,-1_I32)
  CALL lsyssc_addIndex (rtempMatrix%h_Kld,-1_I32)
  
  ! Get the data arrays.
  CALL lsyssc_getbase_Kcol (rtempMatrix,p_Kcol)
  CALL lsyssc_getbase_Kld (rtempMatrix,p_Kld)
  CALL lsyssc_getbase_double (rtempMatrix,p_DA)
  
  ! Fill the matrix content by 1.0. That way, UMFPACK will treat
  ! all entries as nonzero.
  CALL lalg_setVectorDble (p_Da,1.0_DP)
  
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
        rsolverNode%p_rsubnodeUMFPACK4%rtempVector,.FALSE.,.FALSE.,&
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
    CALL sys_halt()
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
    CALL sys_halt()
  END IF

  SELECT CASE (p_rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
    ! Format 9 is exactly the UMFPACK matrix.
    IF (p_rmatrix%dScaleFactor .EQ. 1.0_DP) THEN
      ! Make a copy of the matrix structure, but use the same matrix entries.
      CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                    LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
    ELSE
      ! Make a copy of the whole matrix. We will resolve the scaling factor later,
      ! which changes the entries!
      CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                    LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    END IF
  CASE (LSYSSC_MATRIX7)
    ! For format 7, we have to modify the matrix slightly.
    ! Make a copy of the whole matrix:
    CALL lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                 LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    ! Resort the entries to put the diagonal entry to the correct position.
    ! This means: Convert the structure-7 matrix to a structure-9 matrix:
    CALL lsyssc_convertMatrix (p_rmatrix,LSYSSC_MATRIX9)
  END SELECT
  
  IF (rtempMatrix%dscaleFactor .NE. 1.0_DP) THEN
    ! The matrix entries have been duplicated above in this case, so we are
    ! allowed to change the entries of rtempMatrix without changing the
    ! original entries. So, 'un'-scale the matrix.
    CALL lsyssc_scaleMatrix (rtempMatrix,rtempMatrix%dscaleFactor)
    rtempMatrix%dscaleFactor = 1.0_DP
  END IF

  ! If the debug flag is set, write out the matrix to a text file.
  IF (rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput .GT. 0) THEN
    IF (rsolverNode%ioutputLevel .GE. 0) THEN
      CALL output_line ('Writing matrix to a text file.',&
          OU_CLASS_WARNING,OU_MODE_STD,'linsol_initDataUMFPACK4')
    END IF
    CALL lsyssc_getbase_double (rtempMatrix,p_DA)
    WHERE (abs(p_Da) .LT. 1.0E-12_DP) p_Da = 0.0_DP
    CALL matio_writeMatrixHR (p_rmatrix, 'matrix',.TRUE., 0, 'matrix'//&
        TRIM(sys_siL(rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput,10))//&
        '.txt', '(E15.5)')
  END IF
  
  ! DEBUG!!!
  !CALL sys_halt()

  ! Modify Kcol/Kld of the matrix. Subtract 1 to get them 0-based.
  CALL lsyssc_addIndex (rtempMatrix%h_Kcol,-1_I32)
  CALL lsyssc_addIndex (rtempMatrix%h_Kld,-1_I32)
  
  ! Get the data arrays.
  CALL lsyssc_getbase_Kcol (rtempMatrix,p_Kcol)
  CALL lsyssc_getbase_Kld (rtempMatrix,p_Kld)
  CALL lsyssc_getbase_double (rtempMatrix,p_DA)
  
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
  CASE (-11)
    ! no memory
    ierror = LINSOL_ERR_MATRIXHASCHANGED
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
  INTEGER :: KSYS
  REAL(DP) :: dres
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
      CALL sys_halt()
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

    IF (rsolverNode%ioutputLevel .GE. 2) THEN
      dres = lsysbl_vectorNorm (rd,rsolverNode%iresNorm)
      IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                (dres .LE. 1E99_DP))) dres = 0.0_DP

      CALL output_line ('UMFPACK: !!Initial RES!! = '//TRIM(sys_sdEL(dres,15)) )
    END IF

    ! Solve the system
    ! Solve the system. Note that UMFPACK expects the matrix in
    ! CSR format, which is transposed to our matrix format 9 --
    ! So we solve the transposed system:
    KSYS = 1
    CALL UMF4SOL(KSYS,p_Dx,p_Db,rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
                 rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,Dinfo)
                
    IF (rsolverNode%ioutputLevel .GE. 1) THEN
      dres = lsysbl_vectorNorm (rd,rsolverNode%iresNorm)
      IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                (dres .LE. 1E99_DP))) dres = 0.0_DP

      CALL output_line ('UMFPACK: !!solution!!    = '//TRIM(sys_sdEL(dres,15))//&
          ' Status = '//TRIM(sys_siL(INT(Dinfo(1)),10)) )
    END IF
                
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
  
  RECURSIVE SUBROUTINE linsol_matCompatMILUs1x1 (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! But we cannot handle it if it's not scalar!
    IF (Rmatrices(UBOUND(Rmatrices,1))%ndiagBlocks .NE. 1) &
      ccompatible = LINSOL_COMP_ERRNOTSCALAR
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
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
  INTEGER(I32) :: isubgroup,ierr,maxstr
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_matrixScalar), POINTER :: p_rmatrixSc
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld, p_Kcol
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  INTEGER(I32) :: ifill
  REAL(DP) :: drelax
  
  TYPE(t_MILUdecomp), POINTER :: rMILUDecomp
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there's no matrix assigned
    
    IF (rsolverNode%rsystemMatrix%NEQ .EQ. 0) THEN
      PRINT *,'Error: No matrix associated!'
      CALL sys_halt()
    END IF
    
    rMILUDecomp => rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp
    
    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
    
    p_rmatrix => rsolverNode%rsystemMatrix

    ! We only support scalar 1x1 matrices in structure 7 and 9.
    IF (p_rmatrix%ndiagBlocks .NE. 1) THEN
      PRINT *,'(M)ILU(s) supports only 1x1 matrices!'
      CALL sys_halt()
    END IF

    p_rmatrixSc => p_rmatrix%RmatrixBlock(1,1)

    ! We only support scalar 1x1 matrices in structure 7 and 9.
    IF (p_rmatrixSc%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'(M)ILU(s) supports only double precision matrices!'
      CALL sys_halt()
    END IF
    
    IF ((p_rmatrixSc%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND.&
        (p_rmatrixSc%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
      PRINT *,'(M)ILU(s) supports only structure 7 and 9 matrices!'
      CALL sys_halt()
    END IF
    
    ! Get the matrix description
    CALL lsyssc_getbase_Kld (p_rmatrixSc,p_Kld)
    CALL lsyssc_getbase_Kcol (p_rmatrixSc,p_Kcol)
    CALL lsyssc_getbase_double (p_rmatrixSc,p_DA)

    maxstr = 0

    ! Now calculate the decomposition. Probably allocate more memory if it's
    ! not enough.
    ! Calculate the (M)ILUs matrix using SPLIB.
    ! Get the parameters from the solver node; they were set during the
    ! initialisatin of the solver.
    ifill = rsolverNode%p_rsubnodeMILUs1x1%ifill
    drelax = rsolverNode%p_rsubnodeMILUs1x1%drelax
    
    CALL iluk_ilu(p_rmatrixSc%NEQ,p_DA,p_Kcol,p_Kld, &
              ifill,drelax,&
              ierr,&
              rMILUDecomp)
              
   
    ! Error?
    SELECT CASE (ierr)
    CASE (:-1)
      PRINT *,'Warning: (M)ILU(s) decomposition singular!'
    CASE (0)
      ! everything ok
    CASE DEFAULT
      ierror = LINSOL_ERR_INITERROR
      RETURN
    END SELECT
    
    ! now here we can reallocate the jlu array in our structure
    IF (rMILUDecomp%h_jlu .NE. ST_NOHANDLE) THEN
      ! If less than the half of the memory ILU wanted to have is used,
      ! we reallocate the memory. It does not make sense to have that much
      ! waste!
      
      ! nzlu is the number of bytes needed for the jlu array
      ! reallocate if it uses more than twice the amount, ILU wanted      
      IF (rMILUDecomp%nzlu .LT. rMILUDecomp%isize/2) THEN
        call storage_realloc('linsol_initDataMILUs1x1', rMILUDecomp%nzlu, &
             rMILUDecomp%h_jlu, ST_NEWBLOCK_ZERO, .true.)
      END IF
      
      ! Save 1/scaling factor of the matrix -- to support scaled matrices
      ! when preconditioning.
      rsolverNode%p_rsubnodeMILUs1x1%dscaleFactor = 1.0_DP/p_rmatrixSc%dScaleFactor
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
    ! 
    IF(rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_lu .ne. ST_NOHANDLE) THEN
      call iluk_freeDecomp(rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp)
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
    REAL(DP), DIMENSION(:), POINTER :: p_Dd, p_lu
    INTEGER(PREC_MATIDX) :: lu,jlu,ilup
    INTEGER(I32), DIMENSION(:), POINTER :: p_jlu,p_ilup

    IF (rd%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'(M)ILU(s) only supports double precision vectors!'
      CALL sys_halt()
    END IF
    
    ! Get the data array of rd
    CALL lsysbl_getbase_double (rd,p_Dd)
    
    ! Get MILUs information from the parameter block
    
    lu = rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_lu 
    jlu = rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_jlu 
    ilup = rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_ilup
    
    ! When the scaling factor is not = 1, scale the vector before
    ! preconditioning. This emulates: d = (cA)^-1 d = A^-1 (x/c)!
    ! (The value saved in the structure is 1/c!)
    CALL lsysbl_scaleVector(rd,rsolverNode%p_rsubnodeMILUs1x1%dscaleFactor)

    ! Solve the system. Call SPLIB, this overwrites the defect vector
    ! with the preconditioned one.
    !CALL lusolt (INT(SIZE(p_Dd),I32),p_Dd, p_Iwork(lu:), &
    !             p_Iwork(jlu:), p_Iwork(ilup:))

    ! Without the following pointers, the INTEL compiler would create temporary
    ! arrays which may lead to a SEGFAULT because of a full stack!

    call storage_getbase_int(jlu, p_jlu)
    call storage_getbase_int(ilup, p_ilup)
    call storage_getbase_double(lu, p_lu)

    CALL iluk_lusolt (INT(SIZE(p_Dd),I32),p_Dd,p_lu,p_jlu,p_ilup)
                 
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the CG solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initCG (p_rsolverNode,p_rpreconditioner,p_Rfilter)
  
!<description>
  ! Creates a t_linsolNode solver structure for the CG solver. The node
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
    p_rsolverNode%calgorithm = LINSOL_ALG_CG
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK    + &
                                LINSOL_ABIL_CHECKDEF + &
                                LINSOL_ABIL_USESUBSOLVER + &
                                LINSOL_ABIL_USEFILTER
    
    ! Allocate the subnode for CG.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeCG)
    
    ! Attach the preconditioner if given. 
    
    IF (PRESENT(p_rpreconditioner)) THEN 
      p_rsolverNode%p_rsubnodeCG%p_rpreconditioner => p_rpreconditioner
    END IF

    ! Attach the filter if given. 
    
    IF (PRESENT(p_Rfilter)) THEN
      p_rsolverNode%p_rsubnodeCG%p_RfilterChain => p_Rfilter
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_alterCG (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the CG solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration 
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there's a preconditioner attached. If yes, pass the command
    ! structure to that one.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_alterSolver(rsolverNode%p_rsubnodeCG%p_rpreconditioner,&
          ralterConfig)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatCG (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix. This solver can usually handle 
    ! everything.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible
    
    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_matricesCompatible ( &
          rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
          Rmatrices,ccompatible,CcompatibleDetail)
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixCG (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of CG
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
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_setMatrices (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                              Rmatrices)
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureCG (rsolverNode, ierror,isolverSubgroup)
  
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
    TYPE(t_linsolSubnodeCG), POINTER :: p_rsubnode
    
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
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_initStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                                isubgroup)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be initialised
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

    ! CG needs 3 temporary vectors + 1 for preconditioning. 
    ! Allocate that here! Use the default data type prescribed in the solver 
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeCG
    DO i=1,4
      CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%RtempVectors(i),.FALSE.,.FALSE.,&
          rsolverNode%cdefaultDataType)
    END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataCG (rsolverNode, ierror,isolverSubgroup)
  
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
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_initData (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                            isubgroup,ierror)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataCG (rsolverNode, isolverSubgroup)
  
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
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_doneData (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                            isubgroup)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureCG (rsolverNode, isolverSubgroup)
  
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
    TYPE(t_linsolSubnodeCG), POINTER :: p_rsubnode
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_doneStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                                isubgroup)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be released
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeCG
    IF (p_rsubnode%RtempVectors(1)%NEQ .NE. 0) THEN
      DO i=4,1,-1
        CALL lsysbl_releaseVector (p_rsubnode%RtempVectors(i))
      END DO
    END IF
      
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneCG (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the CG solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it's also released from the heap by calling 
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the CG solver.
  TYPE(t_linsolNode), POINTER         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there's a preconditioner attached. If yes, release it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL linsol_releaseSolver(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    END IF
    
    ! Release memory if still associated
    CALL linsol_doneDataCG (rsolverNode, rsolverNode%isolverSubgroup)
    CALL linsol_doneStructureCG (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the CG subnode
    DEALLOCATE(rsolverNode%p_rsubnodeCG)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precCG (rsolverNode,rd)
  
!<description>
  ! Applies CG preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: dalpha,dbeta,dgamma,dgammaOld,dres,dfr
  INTEGER :: ireslength,ite,i,j,niteAsymptoticCVR

  ! The queue saves the current residual and the two previous residuals.
  REAL(DP), DIMENSION(32) :: Dresqueue
  
  ! The system matrix
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  INTEGER :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  LOGICAL bprec,bfilter
  
  ! Our structure
  TYPE(t_linsolSubnodeCG), POINTER :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  TYPE(t_vectorBlock), POINTER :: p_DR,p_DP,p_DD,p_rx
  TYPE(t_linsolNode), POINTER :: p_rprecSubnode
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Get some information
    p_rsubnode => rsolverNode%p_rsubnodeCG
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

    ireslength = 32 !MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))
    niteAsymptoticCVR = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    bfilter = ASSOCIATED(rsolverNode%p_rsubnodeCG%p_RfilterChain)
    
    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! Set pointers to the temporary vectors
    ! residuum vector
    p_DR   => p_rsubnode%RtempVectors(1)
    ! preconditioned residuum vector
    p_DP   => p_rsubnode%RtempVectors(2)
    ! direction vector
    p_DD   => p_rsubnode%RtempVectors(3)
    ! solution vector
    p_rx   => p_rsubnode%RtempVectors(4)
    
    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    CALL lsysbl_assignDiscretIndirect (rd,p_DR)
    CALL lsysbl_assignDiscretIndirect (rd,p_DP)
    CALL lsysbl_assignDiscretIndirect (rd,p_DD)
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
      
    CALL lsysbl_clearVector(p_DR)
    CALL lsysbl_clearVector(p_DP)
    CALL lsysbl_clearVector(p_DD)
    
    ! Initialization

    dalpha = 1.0_DP
    dbeta  = 1.0_DP
    dgamma = 1.0_DP
    dgammaOld = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.

    CALL lsysbl_copyVector(rd,p_DR)
    IF (bfilter) THEN
      ! Apply the filter chain to the vector
      CALL filter_applyFilterChainVec (p_DR, p_RfilterChain)
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
        CALL output_line ('CG: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      ! Copy the residuum vector p_DR to the preconditioned one.
      CALL lsysbl_copyVector(p_DR,p_DP)
    
      IF (bprec) THEN
        ! Perform preconditioning with the assigned preconditioning
        ! solver structure.
        CALL linsol_precondDefect (p_rprecSubnode,p_DP)
      END IF

      ! Calculate gamma
      dgamma = lsysbl_scalarProduct (p_DR, p_DP)
      
      ! Perform at most nmaxIterations loops to get a new vector
      
      ! Copy the preconditioned residual vector to the direction vector
      CALL lsysbl_copyVector (p_DP, p_DD)
      ! And multiply it with -1.0
      CALL lsysbl_scaleVector (p_DD, -1.0_DP)

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 1: Calculate alpha
        !
        !            < g_k , h_k >           gamma
        ! alpha := ----------------- = -----------------
        !          < d_k , A * d_k >   < d_k , A * d_k >
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! We will now abuse the p_DP vector to save the result of
        ! the matrix-vector-multiplication A * d_k.
        ! Using this trick, we can save one temporary vector.
        CALL lsysbl_blockMatVec (p_rmatrix,p_DD,p_DP, 1.0_DP,0.0_DP)

        ! Calculate alpha for the CG iteration
        dalpha = dgamma / lsysbl_scalarProduct (p_DD, p_DP)

        IF (dalpha .EQ. 0.0_DP) THEN
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line('CG: Convergence failed, ALPHA=0!')
          END IF
          rsolverNode%iresult = -2
          EXIT
        END IF
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 2: Calculate x_{k+1}
        !
        ! x_{k+1} := x_k - alpha * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        CALL lsysbl_vectorLinearComb (p_DD, p_rx, -1.0_DP * dalpha, 1.0_DP)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 3: Calculate g_{k+1}
        !
        ! g_{k+1} := g_k + alpha * A * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Since we have abused p_DP for saving A * d_k, we can use it now
        ! for the linear combination of g_{k+1}
        CALL lsysbl_vectorLinearComb (p_DP, p_DR, dalpha, 1.0_DP)

        IF (bfilter) THEN
          ! Apply the filter chain to the new defect vector
          CALL filter_applyFilterChainVec (p_DR, p_RfilterChain)
        END IF

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 4: Calculate ||g_{k+1}|| and write some output
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
          IF (rsolverNode%ioutputLevel .LT. 2) THEN
            DO i=MAX(1,SIZE(Dresqueue)-ite-1+1),SIZE(Dresqueue)
              j = ITE-MAX(1,SIZE(Dresqueue)-ite+1)+i
              CALL output_line ('CG: Iteration '// &
                  TRIM(sys_siL(j,10))//',  !!RES!! = '//&
                  TRIM(sys_sdEL(Dresqueue(i),15)) )
            END DO
          END IF
          CALL output_line('CG: Solution diverging!')
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
          CALL output_line ('CG: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 5: Calculate h_{k+1}
        !
        ! In other words: Copy the residual vector and apply the
        ! preconditioner, if given.
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Copy the residuum vector p_DR to the preconditioned one.
        CALL lsysbl_copyVector(p_DR,p_DP)
    
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL linsol_precondDefect (p_rprecSubnode,p_DP)
        END IF
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 6: Calculate gamma
        !
        ! gamma' := gamma
        ! gamma  := < g_{k+1}, h_{k+1} >
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        dgammaOld = dgamma
        dgamma = lsysbl_scalarProduct (p_DR, p_DP)
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 7: Calculate beta
        !
        !          gamma
        ! beta := --------
        !          gamma'
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        dbeta = dgamma / dgammaOld

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 8: Calculate d_{k+1}
        !
        ! d_{k+1} := beta * d_k - h_{k+1}
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        CALL lsysbl_vectorLinearComb (p_DP, p_DD, -1.0_DP, dbeta)
        
        ! That's it - next iteration!
      END DO

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) THEN
        ! Warning if we didn't reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        IF ((rsolverNode%ioutputLevel .GE. 0) .AND. &
            (rsolverNode%nmaxIterations .GT. rsolverNode%nminIterations)) THEN
          CALL output_line ('CG: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        END IF

        ite = rsolverNode%nmaxIterations
      END IF

      ! Finish - either with an error or if converged.
      ! Print the last residuum.


      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        CALL output_line ('CG: Iteration '// &
            TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
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
    
      IF (niteAsymptoticCVR .NE. 0) THEN
        I = MAX(33-ite,33-niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(33-I,DP))
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
        CALL output_lbrk()
        CALL output_line ('CG statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'CG: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF

    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE

! *****************************************************************************
! Routines for the BiCGStab solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_Rfilter,&
                                 cprecondType)
  
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
  
  ! OPTIONAL: An integer specifying what type of preconditioning we want to
  ! use. Can be one of the LINSOL_BICGSTAB_* constants (see contant block).
  ! If not given, Left-preconditioning is used.
  INTEGER, OPTIONAL :: cprecondType
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
    
    ! Set brightPrecond
    
    IF (PRESENT(cprecondType)) THEN
      ! We explicitly check if cprecondType holds a valid constant to avoid
      ! errors during further execution of the algorithm.
      IF (cprecondType .EQ. LINSOL_BICGSTAB_RIGHT_PRECOND) THEN
        p_rsolverNode%p_rsubnodeBiCGStab%cprecondType = LINSOL_BICGSTAB_RIGHT_PRECOND
      ELSE
        p_rsolverNode%p_rsubnodeBiCGStab%cprecondType = LINSOL_BICGSTAB_LEFT_PRECOND
      END IF
    ELSE
      ! If no preconditioner type has been explicitly requested, we use
      ! left-preconditioning by default.
      p_rsolverNode%p_rsubnodeBiCGStab%cprecondType = LINSOL_BICGSTAB_LEFT_PRECOND
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_alterBiCGStab (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the BiCGStab solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration 
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there's a preconditioner attached. If yes, pass the command
    ! structure to that one.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_alterSolver(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner,&
          ralterConfig)
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatBiCGStab (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix. This solver can usually handle 
    ! everything.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible

    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL linsol_matricesCompatible ( &
          rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
          Rmatrices,ccompatible,CcompatibleDetail)
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
          p_rsubnode%RtempVectors(i),.FALSE.,.FALSE.,rsolverNode%cdefaultDataType)
    END DO
    
    ! If we want to use right-preconditioning, we need one more
    ! temporary vector.
    IF (p_rsubnode%cprecondType .EQ. LINSOL_BICGSTAB_RIGHT_PRECOND) THEN
      CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%rprecondTemp,.FALSE.,.FALSE.,rsolverNode%cdefaultDataType)
    END IF      
  
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
    
    ! Release temporary vector for right-preconditioning
    IF (p_rsubnode%cprecondType .EQ. LINSOL_BICGSTAB_RIGHT_PRECOND) THEN
      IF (p_rsubnode%rprecondTemp%NEQ .NE. 0) THEN
        CALL lsysbl_releaseVector (p_rsubnode%rprecondTemp)
      END IF
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
!</description>

!<inputoutput>
  ! The t_linsolNode structure of the BiCGStab solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>

!</subroutine>

    ! Now we check what type of preconditioning we have to use.
    ! If our preconditioning type is LINSOL_BICGSTAB_LEFT_PRECOND, we try to solve
    ! the system 
    !                       P^{-1} A * x = P^{1} * b
    ! where P is the preconditioner matrix.
    !
    ! If our preconditioning type is LINSOL_BICGSTAB_RIGHT_PRECOND, we try to solve
    ! the system(s)
    !                        AP^{-1} * y = b
    !                              P * x = y
    ! where P is the preconditioner matrix.
    !
    ! Depending on which constant has been passed to linsol_initBiCGStab, we call
    ! the left or right variant of the BiCGStab solver here.
    IF (rsolverNode%p_rsubnodeBiCGStab%cprecondType .EQ. &
        LINSOL_BICGSTAB_RIGHT_PRECOND) THEN
      CALL linsol_precBiCGStabRight(rsolverNode,rd)
    ELSE
      CALL linsol_precBiCGStabLeft(rsolverNode,rd)
    END IF
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE linsol_precBiCGStabLeft (rsolverNode,rd)
  
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
  INTEGER :: ireslength,ite,i,j,niteAsymptoticCVR

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

    ireslength = 32 !MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))
    niteAsymptoticCVR = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

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
        CALL output_line ('BiCGStab: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      CALL lsysbl_copyVector(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        drho1 = lsysbl_scalarProduct (p_DR0,p_DR) 

        IF (drho0*domega0 .EQ. 0.0_DP) THEN
          ! Should not happen
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line ('BiCGStab: Iteration prematurely stopped! '//&
                 'Correction vector is zero!')
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
        
        IF (ABS(dalpha) .EQ. 0.0_DP) THEN
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line ('BiCGStab: Convergence failed, ALPHA=0!')
          END IF
          rsolverNode%iresult = -2
          EXIT
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
              CALL output_line ('BiCGStab: Convergence failed: omega=0!')
            END IF
            rsolverNode%iresult = -2
            EXIT
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
          IF (rsolverNode%ioutputLevel .LT. 2) THEN
            DO i=MAX(1,SIZE(Dresqueue)-ite-1+1),SIZE(Dresqueue)
              j = ITE-MAX(1,SIZE(Dresqueue)-ite+1)+i
              CALL output_line ('BiCGStab: Iteration '// &
                  TRIM(sys_siL(j,10))//',  !!RES!! = '//&
                  TRIM(sys_sdEL(Dresqueue(i),15)) )
            END DO
          END IF
          CALL output_line ('BiCGStab: Solution diverging!')
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
          CALL output_line ('BiCGStab: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF

      END DO

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) THEN
        ! Warning if we didn't reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        IF ((rsolverNode%ioutputLevel .GE. 0) .AND. &
            (rsolverNode%nmaxIterations .GT. rsolverNode%nminIterations)) THEN
          CALL output_line ('BiCGStab: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        END IF

        ite = rsolverNode%nmaxIterations
      END IF

      ! Finish - either with an error or if converged.
      ! Print the last residuum.


      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        CALL output_line ('BiCGStab: Iteration '// &
            TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
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
    
      IF (niteAsymptoticCVR .NE. 0) THEN
        I = MAX(33-ite,33-niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(33-I,DP))
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
        CALL output_lbrk()
        CALL output_line ('BiCGStab statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'BiCGStab: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precBiCGStabRight (rsolverNode,rd)
  
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
  INTEGER :: ireslength,ite,i,j,niteAsymptoticCVR

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
  TYPE(t_vectorBlock), POINTER :: p_DR,p_DR0,p_DP,p_DPA,p_DSA,p_DZ,p_rx
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

    ireslength = 32 !MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))
    niteAsymptoticCVR = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

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
    p_DZ   => p_rsubnode%rprecondTemp
    
    
    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    CALL lsysbl_assignDiscretIndirect (rd,p_DR )
    CALL lsysbl_assignDiscretIndirect (rd,p_DR0)
    CALL lsysbl_assignDiscretIndirect (rd,p_DP )
    CALL lsysbl_assignDiscretIndirect (rd,p_DPA)
    CALL lsysbl_assignDiscretIndirect (rd,p_DSA)
    CALL lsysbl_assignDiscretIndirect (rd,p_rx)
    CALL lsysbl_assignDiscretIndirect (rd,p_DZ)
    
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
    !IF (bprec) THEN
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      !CALL linsol_precondDefect (p_rprecSubnode,p_DR)
    !END IF
    
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
        CALL output_line ('BiCGStab: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      CALL lsysbl_copyVector(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        drho1 = lsysbl_scalarProduct (p_DR0,p_DR) 

        IF (drho0*domega0 .EQ. 0.0_DP) THEN
          ! Should not happen
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line ('BiCGStab: Iteration prematurely stopped! '//&
                 'Correction vector is zero!')
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
        
        ! We need to copy p_DP to p_DZ for preconditioning now
        CALL lsysbl_copyVector(p_DP, p_DZ)
        
        ! Call preconditioner
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL linsol_precondDefect (p_rprecSubnode,p_DZ)
        END IF

        CALL lsysbl_blockMatVec (p_rmatrix, p_DZ,p_DPA, 1.0_DP,0.0_DP)
        
        IF (bfilter) THEN
          ! Apply the filter chain to the vector
          CALL filter_applyFilterChainVec (p_DPA, p_RfilterChain)
        END IF

        dalpha = lsysbl_scalarProduct (p_DR0,p_DPA)
        
        IF (dalpha .EQ. 0.0_DP) THEN
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line ('BiCGStab: Convergence failed, ALPHA=0!')
          END IF
          rsolverNode%iresult = -2
          EXIT
        END IF
        
        dalpha = drho1/dalpha

        CALL lsysbl_vectorLinearComb (p_DPA,p_DR,-dalpha,1.0_DP)

        ! We need to copy p_DR to p_DZ for preconditioning now
        CALL lsysbl_copyVector(p_DR, p_DZ)
        
        ! Call preconditioner
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL linsol_precondDefect (p_rprecSubnode,p_DZ)
        END IF

        CALL lsysbl_blockMatVec (p_rmatrix, p_DZ,p_DSA, 1.0_DP,0.0_DP)
        
        IF (bfilter) THEN
          ! Apply the filter chain to the vector
          CALL filter_applyFilterChainVec (p_DSA, p_RfilterChain)
        END IF
        
        domega1 = lsysbl_scalarProduct (p_DSA,p_DR)
        domega2 = lsysbl_scalarProduct (p_DSA,p_DSA)
        
        IF (domega1 .EQ. 0.0_DP) THEN
          domega0 = 0.0_DP
        ELSE
          IF (domega2 .EQ. 0.0_DP) THEN
            IF (rsolverNode%ioutputLevel .GE. 2) THEN
              CALL output_line ('BiCGStab: Convergence failed: omega=0!')
            END IF
            rsolverNode%iresult = -2
            EXIT
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
          IF (rsolverNode%ioutputLevel .LT. 2) THEN
            DO i=MAX(1,SIZE(Dresqueue)-ite-1+1),SIZE(Dresqueue)
              j = ITE-MAX(1,SIZE(Dresqueue)-ite+1)+i
              CALL output_line ('BiCGStab: Iteration '// &
                  TRIM(sys_siL(j,10))//',  !!RES!! = '//&
                  TRIM(sys_sdEL(Dresqueue(i),15)) )
            END DO
          END IF
          CALL output_line ('BiCGStab: Solution diverging!')
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
          CALL output_line ('BiCGStab: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF

      END DO

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) THEN
        ! Warning if we didn't reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        IF ((rsolverNode%ioutputLevel .GE. 0) .AND. &
            (rsolverNode%nmaxIterations .GT. rsolverNode%nminIterations)) THEN
          CALL output_line ('BiCGStab: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        END IF

        ite = rsolverNode%nmaxIterations
      END IF

      ! Finish - either with an error or if converged.
      ! Print the last residuum.


      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        CALL output_line ('BiCGStab: Iteration '// &
            TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      END IF

    END IF

    rsolverNode%iiterations = ite
    
    ! Since we are using right-preconditioning here, we need to apply the
    ! preconditioner to our iteration vector here.
    IF (bprec) THEN
      CALL linsol_precondDefect (p_rprecSubnode,p_rx)
    END IF

    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    CALL lsysbl_copyVector (p_rx,rd)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
    
      ! Calculate asymptotic convergence rate
    
      IF (niteAsymptoticCVR .NE. 0) THEN
        I = MAX(33-ite,33-niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(33-I,DP))
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
        CALL output_lbrk()
        CALL output_line ('BiCGStab statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'BiCGStab: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE


! *****************************************************************************
! Routines for the GMRES(m) solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initGMRES (p_rsolverNode,ikrylovDim,p_rpreconditioner,&
                               p_Rfilter,btwiceGS,dpseudoResScale)
  
!<description>
  ! Creates a t_linsolNode solver structure for the GMRES(m) solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! Maximal dimension of the Krylov subspace for the GMRES iteration.
  ! Is qual to the number of GMRES iterations before a restart is done.
  ! Must not be smaller than 1, otherwise linsol_initStructure will fail!
  INTEGER :: ikrylovDim
  
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

  ! OPTIONAL: A boolean which specifies whether we should apply the Gram-Schmidt
  ! process once or twice per GMRES iteration.
  ! Since the Gram-Schmidt algorithm suffers under numerical instability, we
  ! have a chance to improve its stability by applying the Gram-Schmidt process
  ! twice. If btwiceGS is given and it is .TRUE., then the Gram-Schmidt process
  ! will be applied twice per GMRES iteration, otherwise only once.
  LOGICAL, OPTIONAL :: btwiceGS
  
  ! OPTIONAL: Scale factor for the pseudo-residual-norms to test against the
  ! stopping criterion in the inner GMRES loop. Since the pseudo-residual-norms
  ! may be highly inaccurate, this scale factor can be used to scale the
  ! pseudo-residual-norms before testing against the stopping criterion.
  ! It is recommended to use a scale factor > 1, a factor of 2 should be usually
  ! okay. If dpseudoresscale is not given or set to < 1, it is automatically
  ! set to 2.
  ! Setting dpseudoresscale to 0 will disable the pseudo-resirudal norm check.
  REAL(DP), OPTIONAL :: dpseudoResScale
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
    p_rsolverNode%calgorithm = LINSOL_ALG_GMRES
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK    + &
                                LINSOL_ABIL_CHECKDEF + &
                                LINSOL_ABIL_USESUBSOLVER + &
                                LINSOL_ABIL_USEFILTER
    
    ! Allocate the subnode for GMRES.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeGMRES)
    
    ! Save the dimension of the Krylov subspace. This is equal to the number
    ! of GMRES iterations before a restart is done.
    p_rsolverNode%p_rsubnodeGMRES%ikrylovDim = ikrylovDim
    
    ! Save if the Gram-Schmidt process should be applied twice or not.
    IF (PRESENT(btwiceGS)) THEN
      p_rsolverNode%p_rsubNodeGMRES%btwiceGS = btwiceGS
    ELSE
      p_rsolverNode%p_rsubNodeGMRES%btwiceGS = LINSOL_GMRES_DEF_TWICE_GS
    END IF
    
    ! Save the pseudo-residual-norm scale factor if given, or set
    ! to default scale factor if not given or invalid (i.e. < 1)
    IF (PRESENT(dpseudoResScale)) THEN
      IF (dpseudoResScale .GE. 0.0_DP) THEN
        p_rsolverNode%p_rsubNodeGMRES%dpseudoResScale = dpseudoResScale
      ELSE
        p_rsolverNode%p_rsubNodeGMRES%dpseudoResScale = 0.0_DP
      ENDIF
    ELSE
      p_rsolverNode%p_rsubNodeGMRES%dpseudoResScale = 0.0_DP
    ENDIF
    
    ! Attach the preconditioner if given. 
    IF (PRESENT(p_rpreconditioner)) THEN 
      p_rsolverNode%p_rsubNodeGMRES%p_rpreconditioner => p_rpreconditioner
    END IF

    ! Attach the filter if given. 
    IF (PRESENT(p_Rfilter)) THEN
      p_rsolverNode%p_rsubNodeGMRES%p_RfilterChain => p_Rfilter
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_alterGMRES (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the GMRES solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration 
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there's a preconditioner attached. If yes, pass the command
    ! structure to that one.
    IF (ASSOCIATED(rsolverNode%p_rsubNodeGMRES%p_rpreconditioner)) THEN
      CALL linsol_alterSolver(rsolverNode%p_rsubNodeGMRES%p_rpreconditioner,&
          ralterConfig)
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE linsol_matCompatGMRES (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix. This solver can usually handle 
    ! everything.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    IF (PRESENT(CcompatibleDetail)) &
        CcompatibleDetail (UBOUND(CcompatibleDetail,1)) = ccompatible

    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) THEN
      CALL linsol_matricesCompatible ( &
          rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
          Rmatrices,ccompatible,CcompatibleDetail)
    END IF
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixGMRES (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of GMRES(m)
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) THEN
      CALL linsol_setMatrices (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                              Rmatrices)
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureGMRES (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES solver
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
    INTEGER(I32) :: idim
    INTEGER(I32) , DIMENSION(2) :: idim2
    TYPE(t_linsolSubnodeGMRES), POINTER :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    p_rsubnode => rsolverNode%p_rsubnodeGMRES

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(p_rsubnode%p_rpreconditioner)) THEN
      CALL linsol_initStructure (p_rsubnode%p_rpreconditioner, isubgroup)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be initialised
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN
    
    ! We now need to check if the dimension of the Krylov subspace is
    ! positive. If it is not, we need to cancel the initialization here.
    IF (p_rsubnode%ikrylovDim .LE. 0) THEN
      PRINT *, "Error: Dimension of Krylov subspace for GMRES(m) is <= 0 !"
      ierror = LINSOL_ERR_INITERROR
      CALL sys_halt()
    END IF
    
    ! Get the stuff out of our solver node
    idim = p_rsubnode%ikrylovDim
    idim2(1) = idim
    idim2(2) = idim

    ! Now comes the intersting part - we need to allocate a lot of arrays
    ! and vectors for the Krylov subspace for GMRES here.
    
    ! Call our storage to allocate the 1D/2D arrays
    CALL storage_new('linsol_initStructureGMRES', 'Dh', idim2, &
        ST_DOUBLE, p_rsubnode%hDh, ST_NEWBLOCK_NOINIT)
    CALL storage_new('linsol_initStructureGMRES', 'Dc', idim, &
        ST_DOUBLE, p_rsubnode%hDs, ST_NEWBLOCK_NOINIT)
    CALL storage_new('linsol_initStructureGMRES', 'Ds', idim, &
        ST_DOUBLE, p_rsubnode%hDc, ST_NEWBLOCK_NOINIT)
    CALL storage_new('linsol_initStructureGMRES', 'Dq', idim+1_I32, &
        ST_DOUBLE,  p_rsubnode%hDq, ST_NEWBLOCK_NOINIT)
    
    ! Get the pointers
    CALL storage_getbase_double2D(p_rsubnode%hDh, p_rsubnode%Dh)
    CALL storage_getbase_double(p_rsubnode%hDc, p_rsubnode%Dc)
    CALL storage_getbase_double(p_rsubnode%hDs, p_rsubnode%Ds)
    CALL storage_getbase_double(p_rsubnode%hDq, p_rsubnode%Dq)
    
    ! Allocate space for our auxiliary vectors
    ALLOCATE(p_rsubnode%p_rv(idim+1))
    ALLOCATE(p_rsubnode%p_rz(idim))
    
    ! Create them
    DO i=1, idim+1
      CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%p_rv(i),.FALSE.,.FALSE.,rsolverNode%cdefaultDataType)
    END DO
    DO i=1, idim
      CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%p_rz(i),.FALSE.,.FALSE.,rsolverNode%cdefaultDataType)
    END DO
    
    ! Create an iteration vector x
    CALL lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
        p_rsubnode%rx,.FALSE.,.FALSE.,rsolverNode%cdefaultDataType)
    
    ! Okay, that's all!
      
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataGMRES (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
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
    IF (ASSOCIATED(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) THEN
      CALL linsol_initData (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                            isubgroup,ierror)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataGMRES (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
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
    IF (ASSOCIATED(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) THEN
      CALL linsol_doneData (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                            isubgroup)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureGMRES (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
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
    INTEGER :: isubgroup,i,idim
    TYPE(t_linsolSubnodeGMRES), POINTER :: p_rsubnode
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) THEN
      CALL linsol_doneStructure (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                                isubgroup)
    END IF
    
    ! Cancel here, if we don't belong to the subgroup to be released
    IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeGMRES
    idim = p_rsubnode%ikrylovDim
    
    ! Release auxiliary vectors
    IF (ASSOCIATED(p_rsubnode%p_rz)) THEN
      DO i=1, idim
        CALL lsysbl_releaseVector(p_rsubnode%p_rz(i))
      END DO
      DO i=1, idim+1
        CALL lsysbl_releaseVector(p_rsubnode%p_rv(i))
      END DO
    
      ! Destroy them
      DEALLOCATE(p_rsubnode%p_rz)
      DEALLOCATE(p_rsubnode%p_rv)
    ENDIF

    ! Release iteration vector
    IF (p_rsubnode%rx%NEQ > 0) THEN
      CALL lsysbl_releaseVector(p_rsubnode%rx)
    END IF
    
    ! Destroy our 1D/2D arrays
    IF (p_rsubnode%hDq .NE. ST_NOHANDLE) THEN
      CALL storage_free(p_rsubnode%hDq)
      CALL storage_free(p_rsubnode%hDs)
      CALL storage_free(p_rsubnode%hDc)
      CALL storage_free(p_rsubnode%hDh)
    END IF
    
    ! That's all!
      
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneGMRES (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the GMRES(m) solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it's also released from the heap by calling 
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the GMRES(m) solver.
  TYPE(t_linsolNode), POINTER         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there's a preconditioner attached. If yes, release it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) THEN
      CALL linsol_releaseSolver(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)
    END IF
    
    ! Release memory if still associated
    CALL linsol_doneDataGMRES (rsolverNode, rsolverNode%isolverSubgroup)
    CALL linsol_doneStructureGMRES (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the GMRES subnode
    DEALLOCATE(rsolverNode%p_rsubnodeGMRES)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precGMRES (rsolverNode,rd)
  
!<description>
  ! Applies GMRES(m) preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES solver
  TYPE(t_linsolNode), INTENT(INOUT), TARGET :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: dalpha,dbeta,dres,dfr,dtmp,dpseudores,dprnsf
  INTEGER :: ireslength,ite,i,j,k,niteAsymptoticCVR
  
  ! Here come our 1D/2D arrays
  REAL(DP), DIMENSION(:), POINTER :: p_Dc, p_Ds, p_Dq
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dh

  ! The queue saves the current residual and the two previous residuals.
  REAL(DP), DIMENSION(32) :: Dresqueue
  
  ! The system matrix
  TYPE(t_matrixBlock), POINTER :: p_rmatrix

  
  ! Minimum number of iterations, print-sequence for residuals, dimension
  ! of krylov subspace
  INTEGER :: nminIterations, niteResOutput, idim, hDq
  
  ! Whether to filter/precondition, apply Gram-Schmidt twice
  LOGICAL bprec,bfilter,btwiceGS
  
  ! Our structure
  TYPE(t_linsolSubnodeGMRES), POINTER :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  TYPE(t_vectorBlock), POINTER :: p_rx
  TYPE(t_vectorBlock), DIMENSION(:), POINTER :: p_rv, p_rz
  TYPE(t_linsolNode), POINTER :: p_rprecSubnode
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Get some information
    p_rsubnode => rsolverNode%p_rsubnodeGMRES
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

    ireslength = 32 !MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))
    niteAsymptoticCVR = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(p_rsubnode%p_rpreconditioner)
    bfilter = ASSOCIATED(p_rsubnode%p_RfilterChain)
    
    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)
    
    ! Apply Gram-Schmidt twice?
    btwiceGS = p_rsubnode%btwiceGS
    
    ! Get the dimension of Krylov subspace
    idim = p_rsubnode%ikrylovdim
    
    ! Get our pseudo-residual-norm scale factor
    dprnsf = p_rsubnode%dpseudoResScale
    
    ! Now we need to check the pseudo-residual scale factor.
    ! If it is non-positive, we set it to 0.
    IF (dprnsf .LE. 0.0_DP) THEN
      dprnsf = 0.0_DP
    END IF

    ! Set pointers to the temporary vectors
    ! defect vectors
    p_rz => p_rsubnode%p_rz
    ! basis vectors
    p_rv => p_rsubnode%p_rv
    ! solution vector
    p_rx => p_rsubnode%rx

    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    DO i=1, idim+1
      CALL lsysbl_assignDiscretIndirect (rd,p_rv(i))
    END DO
    DO i=1, idim
      CALL lsysbl_assignDiscretIndirect (rd,p_rz(i))
    END DO
    CALL lsysbl_assignDiscretIndirect (rd,p_rx)
    
    ! Set pointers to the 1D/2D arrays
    p_Dh => p_rsubnode%Dh
    p_Dc => p_rsubnode%Dc
    p_Ds => p_rsubnode%Ds
    p_Dq => p_rsubnode%Dq
    
    ! And get the handle of Dq, since we need to call
    ! storage_clear during the iteration
    hDq = p_rsubnode%hDq
    
    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    CALL lsysbl_assignDiscretIndirect (rd, p_rx)
    DO i=1,idim
      CALL lsysbl_assignDiscretIndirect (rd, p_rz(i))
    END DO
    DO i=1,idim+1
      CALL lsysbl_assignDiscretIndirect (rd, p_rv(i))
    END DO
    
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
      
    ! Copy our RHS rd to p_rv(1). As the iteration vector is 0, this
    ! is also our initial defect.
    CALL lsysbl_copyVector(rd,p_rv(1))
    IF (bfilter) THEN
      CALL filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
    END IF
    
    ! Get the norm of the residuum
    ! We need to calculate the norm twice, since we need one for
    ! the stopping criterion (selected by the user) and the
    ! euclidian norm for the internal GMRES algorithm
    dfr = lsysbl_vectorNorm (p_rv(1),rsolverNode%iresNorm)
    dres = lsysbl_vectorNorm (p_rv(1),LINALG_NORMEUCLID)
    IF (.NOT.((dfr .GE. 1E-99_DP) .AND. &
              (dfr .LE. 1E99_DP))) dfr = 0.0_DP

    ! Initialize starting residuum
    rsolverNode%dinitialDefect = dfr

    ! initialize the queue of the last residuals with RES
    Dresqueue = dfr

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    IF (rsolverNode%dinitialDefect .LT. rsolverNode%drhsZero) THEN
     
      ! final defect is 0, as initialised in the output variable above
      CALL lsysbl_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dfr
    ELSE
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        CALL output_line ('GMRES('//&
             TRIM(sys_siL(idim,10))//'): Iteration '//&
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      ite = 0
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= Outer Loop
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      DO WHILE (ite < rsolverNode%nmaxIterations)
        
        ! Step O.1:
        ! Create elementary vector e_1 scaled by the euclid norm of the
        ! defect, i.e.: q = ( ||v(1)||_2, 0, ..., 0 )
        CALL storage_clear (hDq)
        p_Dq(1) = dres
        
        ! Step O.2:
        ! Now scale the defect by the inverse of its norm
        CALL lsysbl_scaleVector (p_rv(1), 1.0_DP / dres)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! -= Inner Loop (GMRES iterations)
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        DO i = 1, idim
          ! Another iteration
          ite = ite + 1
        
          ! Step I.1:
          ! Solve P * z(i) = v(i), where P is the preconditioner matrix
          CALL lsysbl_copyVector (p_rv(i), p_rz(i))
          
          ! Apply preconditioner to z(i)
          IF (bprec) THEN
            CALL linsol_precondDefect (p_rprecSubnode, p_rz(i))
          END IF
          
          ! Apply filter chain to z(i)
          IF (bfilter) THEN
            CALL filter_applyFilterChainVec (p_rz(i), p_RfilterChain)
          END IF
          
          
          ! Step I.2:
          ! Calculate v(i+1) = A * v(i)
          CALL lsysbl_blockMatVec (p_rmatrix, p_rz(i), p_rv(i+1), 1.0_DP, 0.0_DP)
          
          
          ! Step I.3:
          ! Perfom Gram-Schmidt process (1st time)
          DO k = 1, i
            p_Dh(k,i) = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
            CALL lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -p_Dh(k,i), 1.0_DP)
          END DO
          
          ! If the user wishes, perform Gram-Schmidt one more time
          ! to improve numerical stability.
          IF (btwiceGS) THEN
            DO k = 1, i
              dtmp = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
              p_Dh(k,i) = p_Dh(k,i) + dtmp;
              CALL lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -dtmp, 1.0_DP)
            END DO
          END IF
          
          ! Step I.4:
          ! Calculate alpha = ||v(i+1)||_2
          dalpha = lsysbl_vectorNorm (p_rv(i+1), LINALG_NORMEUCLID)
          
          ! Step I.5:
          ! Scale v(i+1) by the inverse of its euclid norm
          IF (dalpha > SYS_EPSREAL) THEN
            CALL lsysbl_scaleVector (p_rv(i+1), 1.0_DP / dalpha)
          ELSE
            ! Well, let's just print a warning here...
            IF(rsolverNode%ioutputLevel .GE. 2) THEN
              CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))//&
                '): Warning: !!v(i+1)!! < EPS !')
            END IF
          END IF
          
          
          ! Step I.6:
          ! Apply Givens rotations to the i-th column of h which
          ! renders the Householder matrix an upper triangular matrix
          DO k = 1, i-1
            dtmp = p_Dh(k,i)
            p_Dh(k,i)   = p_Dc(k) * dtmp + p_Ds(k) * p_Dh(k+1,i)
            p_Dh(k+1,i) = p_Ds(k) * dtmp - p_Dc(k) * p_Dh(k+1,i)
          END DO
          
          
          ! Step I.7:
          ! Calculate Beta
          ! Beta = (h(i,i)^2 + alpha^2) ^ (1/2)
          dbeta = SQRT(p_Dh(i,i)**2 + dalpha**2)
          IF (dbeta < SYS_EPSREAL) THEN
            dbeta = SYS_EPSREAL
            IF(rsolverNode%ioutputLevel .GE. 2) THEN
              CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))//&
                '): Warning: beta < EPS !')
            END IF
          END IF
          
          
          ! Step I.8:
          ! Calculate next plane rotation
          p_Ds(i) = dalpha / dbeta
          p_Dc(i) = p_Dh(i,i) / dbeta;
          
          
          ! Step I.9
          ! Set h(i,i) = beta
          p_Dh(i,i) = dbeta
          
          
          ! Step I.10:
          ! Update q(i) and calculate q(i+1)
          p_Dq(i+1) = p_Ds(i) * p_Dq(i)
          p_Dq(i)   = p_Dc(i) * p_Dq(i)
          
          
          ! Step I.11
          ! Check pseudo-residual-norm and do all the other stuff we
          ! need to do before starting a new GMRES iteration
          
          ! Get pseudo-residual-norm          
          dpseudores = ABS(p_Dq(i+1))
          
          ! Print our current pseudo-residual-norm
          IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
            CALL output_line ('GMRES('//&
              TRIM(sys_siL(idim,10))//'): Iteration '//&
              TRIM(sys_siL(ITE,10))//',  !q(i+1)! = '//&
              TRIM(sys_sdEL(dpseudores,15)) )
          END IF

          ! The euclid norm of our current defect is implicitly given by
          ! |q(i+1)|, however, it may be inaccurate, therefore, checking if
          ! |q(i+1)| fulfills the stopping criterion would not be very wise,
          ! as it may result in a lot of unnecessary restarts.
          ! We will therefore check (|q(i+1)| * dprnsf) against the tolerances
          ! instead, that should *hopefully* be okay.
          ! If dprnsf is 0, we will check |q(i+1)| against EPS to avoid NaNs
          ! or Infs during the next GMRES iteration.
          dtmp = dpseudores * dprnsf
          
          ! Is |q(i+1)| smaller than our machine's exactness?
          ! If yes, then exit the inner GMRES loop.
          IF (dpseudores .LE. SYS_EPSREAL) THEN
            EXIT
          
          ! Check if (|q(i+1)| * dprnsf) is greater than the machine's
          ! exactness (this can only happen if dprnsf is positive).
          ! If yes, call linsol_testConvergence to test against the
          ! stopping criterions.
          ELSE IF (dtmp > SYS_EPSREAL) THEN
            
            ! If (|q(i+1)| * dprnsf) fulfills the stopping criterion
            ! then exit the inner GMRES loop
            IF (linsol_testConvergence(rsolverNode,dtmp)) EXIT
          
            ! We also want to check if our solution is diverging.
            IF (linsol_testDivergence(rsolverNode,dpseudores)) THEN
              IF(rsolverNode%ioutputLevel .GE. 2) THEN
                CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))//&
                  '): Warning: Pseudo-residuals diverging!')
              END IF
              
              ! Instead of exiting the subroutine, we just exit the inner loop
              EXIT
            END IF
          
          END IF
          
          ! Maybe we have already done enough iterations?
          IF (ite .GE. rsolverNode%nmaxIterations) EXIT
          
          ! Okay, next inner loop iteration, please
        END DO
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! -= End of Inner Loop
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        
        ! Step O.3:
        ! Get the dimension of the current basis of the krylov subspace
        i = MIN(i, idim)
        
        ! Step O.4:
        ! Solve H' * q = q , where H' is the upper triangular matrix of
        ! the upper left (i x i)-block of H.
        ! We use the BLAS Level 2 subroutine 'dtrsv' to do the work for us
        CALL dtrsv('U', 'N', 'N', i, p_Dh, idim, p_Dq, 1)
        
        ! Step O.5:
        ! Update our solution vector
        ! x = x + q(1)*z(1) + q(2)*z(2) + ... + q(i)*z(i)
        DO k = 1, i
          CALL lsysbl_vectorLinearComb(p_rz(k), p_rx, p_Dq(k), 1.0_DP)
        END DO
        
        ! Step O.6:
        ! Calculate 'real' residual
        ! v(1) = b - (A * x)
        CALL lsysbl_copyVector (rd, p_rv(1))
        CALL lsysbl_blockMatVec(p_rmatrix, p_rx, p_rv(1), -1.0_DP, 1.0_DP)

        ! Step O.7:
        ! Call filter chain if given.
        IF (bfilter) THEN
          CALL filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
        END IF
        
        ! Step O.8:
        ! Calculate euclid norm of the residual (needed for next q)
        dres = lsysbl_vectorNorm (p_rv(1), LINALG_NORMEUCLID)
        
        ! Calculate residual norm for stopping criterion
        ! TODO: try to avoid calculating the norm twice
        dfr = lsysbl_vectorNorm (p_rv(1), rsolverNode%iresNorm)
        
        ! Step O.9:
        ! Test for convergence, divergence and write some output now
        
        ! Shift the queue with the last residuals and add the new
        ! residual to it. Check length of ireslength to be larger than
        ! 0 as some compilers might produce Floating exceptions
        ! otherwise! (stupid pgf95)
        IF (ireslength .GT. 0) &
          dresqueue(1:ireslength) = EOSHIFT(dresqueue(1:ireslength),1,dfr)

        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        IF (linsol_testDivergence(rsolverNode,dfr)) THEN
          IF (rsolverNode%ioutputLevel .LT. 2) THEN
            DO i=MAX(1,SIZE(Dresqueue)-ite-1+1),SIZE(Dresqueue)
              j = ITE-MAX(1,SIZE(Dresqueue)-ite+1)+i
              CALL output_line ('GMRES: Iteration '// &
                  TRIM(sys_siL(j,10))//',  !!RES!! = '//&
                  TRIM(sys_sdEL(Dresqueue(i),15)) )
            END DO
          END IF
          CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))//&
            '): Solution diverging!')
          rsolverNode%iresult = 1
          EXIT
        END IF
     
        ! At least perform nminIterations iterations
        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          IF (linsol_testConvergence(rsolverNode,dfr)) EXIT
          
        END IF
        
        ! We need to check if the euclidian norm of the
        ! residual is maybe < EPS - if this is true,
        ! then starting another GMRES iteration would
        ! result in NaNs / Infs!!!
        IF (dres .LE. SYS_EPSREAL) EXIT

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))//&
            '): Iteration '//&
            TRIM(sys_siL(ITE,10))//',  !!RES!!  = '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF

      END DO
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= End of Outer Loop
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) THEN
        ! Warning if we didn't reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        IF ((rsolverNode%ioutputLevel .GE. 0) .AND. &
            (rsolverNode%nmaxIterations .GT. rsolverNode%nminIterations)) THEN
          CALL output_line ('GMRES: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        END IF

        ite = rsolverNode%nmaxIterations
      END IF

      ! Finish - either with an error or if converged.
      ! Print the last residuum.
      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))//&
          '): Iteration '//&
          TRIM(sys_siL(ITE,10))//',  !!RES!!  = '//&
          TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
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
      IF (niteAsymptoticCVR .NE. 0) THEN
        I = MAX(33-ite,33-niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(33-I,DP))
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
        CALL output_lbrk()
        CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))// ') statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line ('GMRES('// TRIM(sys_siL(idim,10))//&
          '): Iterations/Rate of convergence: ' //&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
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
  ! multigrid solver. The application may modify the structure or throw the
  ! pointer away, it doesn't matter. Multigrid will maintain the
  ! structure internally.
  TYPE(t_linsolMGLevelInfo), POINTER     :: p_rlevelInfo
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
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
    CALL sys_halt()
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
    
    ! Clean up the associated matrix if there is one.
    ! Of course, the memory of the matrices will not be deallocated
    ! because the matrices belong to the application, not to the solver!
    CALL lsysbl_releaseMatrix(p_rlevelInfo%rsystemMatrix)

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

!<subroutine>
  
  SUBROUTINE linsol_getMultigridLevel2 (rsolverNode,ilevel,p_rlevelInfo)
                    
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
  TYPE(t_linsolMGLevelInfo2), POINTER     :: p_rlevelInfo
!</output>  
  
!</subroutine>

    NULLIFY(p_rlevelInfo)

    ! Do we have the level?
    IF ((ilevel .GT. 0) .AND. &
        (ilevel .LE. SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo))) THEN
      ! Get it.
      p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    END IF

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
  
  SUBROUTINE linsol_doneMultigrid2Level (rlevelInfo)
                    
!<description>
  ! Cleans up the multigrid level structure rlevelInfo. All memory allocated
  ! in this structure is released.
!</description>
  
!<inputoutput>
  ! The t_levelInfo structure to clean up.
  TYPE(t_linsolMGLevelInfo2), INTENT(INOUT)     :: rlevelInfo
!</inputoutput>

!</subroutine>
  
    ! Release the temporary vectors of this level. 
    ! The vector may not exist - if the application has not called
    ! initStructure!
    IF (rlevelInfo%rrhsVector%NEQ .NE. 0) &
      CALL lsysbl_releaseVector(rlevelInfo%rrhsVector)

    IF (rlevelInfo%rtempVector%NEQ .NE. 0) &
      CALL lsysbl_releaseVector(rlevelInfo%rtempVector)

    IF (rlevelInfo%rsolutionVector%NEQ .NE. 0) &
      CALL lsysbl_releaseVector(rlevelInfo%rsolutionVector)

    ! Release sub-solvers if associated.
    !
    ! Caution: Pre- and postsmoother may be identical!
    IF (ASSOCIATED(rlevelInfo%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(rlevelInfo%p_rpreSmoother, rlevelInfo%p_rpostSmoother))) THEN
      CALL linsol_releaseSolver(rlevelInfo%p_rpostsmoother)
    END IF

    IF (ASSOCIATED(rlevelInfo%p_rpresmoother)) THEN
      CALL linsol_releaseSolver(rlevelInfo%p_rpresmoother)
    END IF

    IF (ASSOCIATED(rlevelInfo%p_rcoarseGridSolver)) THEN
      CALL linsol_releaseSolver(rlevelInfo%p_rcoarseGridSolver)
    END IF
    
    ! Clean up the associated matrix if there is one.
    ! Of course, the memory of the matrices will not be deallocated
    ! because the matrices belong to the application, not to the solver!
    CALL lsysbl_releaseMatrix(rlevelInfo%rsystemMatrix)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_cleanMultigrid2Levels (rsolverNode)
                    
!<description>
  ! This routine removes all level information from the MG solver and releases
  ! all attached solver nodes on all levels (smoothers, coarse grid solvers,
  ! ...). 
  !
  ! The level info array is not released.
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver.
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!</subroutine>

  INTEGER :: ilevel

  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
  END IF

  ! Remove all the levels
  DO ilevel = 1,SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)
    CALL linsol_doneMultigrid2Level (&
        rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel))
  END DO

  END SUBROUTINE
  
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
  
  ! Initialise the coarse grid correction structure.1
  CALL cgcor_init (p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection)

  ! Attach the filter if given. 
  
  IF (PRESENT(p_Rfilter)) THEN
    p_rsolverNode%p_rsubnodeMultigrid%p_RfilterChain => p_Rfilter
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_alterMultigrid (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the Multigrid solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration 
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel

    ! Pass the command structure to all subsolvers and smoothers
    ! on all levels.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    DO WHILE(ASSOCIATED(p_rcurrentLevel))
      IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
        CALL linsol_alterSolver (p_rcurrentLevel%p_rpreSmoother, &
                                 ralterConfig)
      END IF
      ! Pre- and postsmoother may be identical!
      ! Take care not to alter the same smoother twice!
      IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
          (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) THEN
        CALL linsol_alterSolver (p_rcurrentLevel%p_rpostSmoother, &
                                 ralterConfig)
      END IF
      IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
        CALL linsol_alterSolver (p_rcurrentLevel%p_rcoarseGridSolver, &
                                 ralterConfig)
      END IF
      
      ! Next level -- if there is one.    
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel        
    END DO
    
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
    CALL sys_halt()
  END IF
  
  ! Make sure we have the right amount of matrices
  IF (SIZE(Rmatrices) .NE. rsolverNode%p_rsubnodeMultigrid%nlevels) THEN
    PRINT *,'Error: Wrong number of matrices'
    CALL sys_halt()
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
    
    ! Write a link to the matrix into the level-structure.
    CALL lsysbl_duplicateMatrix (Rmatrices(ilevel-nlmin+1), &
        p_rcurrentLevel%rsystemMatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    
    ! Next level -- if there is one.    
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel        
    ilevel = ilevel + 1
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatMultigrid (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  ! More precisely, ccompatible returns always the status of the highest
  ! level where there is an error -- or LINSOL_COMP_OK if there is no error.
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! local variables
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
    INTEGER                       :: ilevel,nlmin,ccompatLevel
    
    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK
    
    ! Reset all compatibility flags
    IF (PRESENT(CcompatibleDetail)) CcompatibleDetail (:) = ccompatible    

    ! Make sure the solver node is configured for multigrid
    IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
        (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
      PRINT *,'Error: Multigrid structure not initialised'
      CALL sys_halt()
    END IF
    
    ! Make sure we have the right amount of matrices
    IF (SIZE(Rmatrices) .NE. rsolverNode%p_rsubnodeMultigrid%nlevels) THEN
      PRINT *,'Error: Wrong number of matrices'
      CALL sys_halt()
    END IF

    ! Check for every level, if there's a presmoother, postsmoother or
    ! coarse grid solver structure attached. If yes, it must check
    ! the matrices!
    ! For this purpose, call the general check routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    nlmin  = LBOUND(Rmatrices,1)
    ilevel = LBOUND(Rmatrices,1)
    
    ! Loop through all levels
    DO WHILE(ASSOCIATED(p_rcurrentLevel))
    
      ! Compatibility flag for that level
      ccompatLevel = LINSOL_COMP_OK
    
      ! Check presmoother
      IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      
        IF (PRESENT(CcompatibleDetail)) THEN
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel,CcompatibleDetail(nlmin:ilevel))
        ELSE
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel)
        END IF
      
        IF (ccompatLevel .NE. LINSOL_COMP_OK) THEN
          ccompatible = ccompatLevel
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          ilevel = ilevel + 1
          CYCLE
        END IF
        
      END IF
      
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
          (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) THEN
        
        IF (PRESENT(CcompatibleDetail)) THEN
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel,CcompatibleDetail(nlmin:ilevel))
        ELSE
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel)
        END IF
        
        IF (ccompatLevel .NE. LINSOL_COMP_OK) THEN
          ccompatible = ccompatLevel
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          ilevel = ilevel + 1
          CYCLE
        END IF
        
      END IF
      
      ! Check coarse grid solver
      IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
        
        IF (PRESENT(CcompatibleDetail)) THEN
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(nlmin:ilevel),ccompatLevel,CcompatibleDetail(nlmin:ilevel))
        ELSE
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(nlmin:ilevel),ccompatLevel)
        END IF
        
        IF (ccompatLevel .NE. LINSOL_COMP_OK) THEN
          ccompatible = ccompatLevel
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          ilevel = ilevel + 1
          CYCLE
        END IF
        
      END IF
      
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
    CALL sys_halt()
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
      IF (ierror .NE. LINSOL_ERR_NOERROR) RETURN
    END IF
    ! Pre- and postsmoother may be identical!
    ! Take care not to initialise the same smoother twice!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
      IF (ierror .NE. LINSOL_ERR_NOERROR) RETURN
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
      IF (ierror .NE. LINSOL_ERR_NOERROR) RETURN
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
            p_rcurrentLevel%rsolutionVector,.FALSE.,.FALSE.,&
            rsolverNode%cdefaultDataType)
      p_rtemplVect => p_rcurrentLevel%rsolutionVector
    END IF
    ! On the coarse grid, we don't need memory for temporary/RHS
    ! vectors. The RHS on the coarse grid is replaced in-situ
    ! by the solution vector.
    IF (ASSOCIATED(p_rcurrentLevel%p_rprevLevel)) THEN
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rrhsVector,.FALSE.,.FALSE.,&
            rsolverNode%cdefaultDataType)
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rtempVector,.FALSE.,.FALSE.,&
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
    !
    ! Don't do anything to the RHS/temp vectors on the coarse grid
    ! and the solution vector on the fine grid -- they don't exist!
    IF (ASSOCIATED(p_rcurrentLevel%p_rprevLevel)) THEN
      p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy = &
        -ABS(p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy)
      p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy = &
        -ABS(p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy)
    END IF
    
    IF (ASSOCIATED(p_rcurrentLevel%p_rnextLevel)) THEN
      p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy = &
        -ABS(p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy)
    END IF
    
    ! And the next level...
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  ! Do we have to allocate memory for prolongation/restriction/resorting?
  IF (imemmax .GT. 0) THEN
    CALL lsyssc_createVector (rsolverNode%p_rsubnodeMultigrid%rprjTempVector,&
                              imemmax,.FALSE.,rsolverNode%cdefaultDataType)
                              
    ! Use this vector also for the coarse grid correction.
    ! Create a block vector on every level with the structure of the 
    ! RHS that shares the memory with rprjTempVector. 
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    DO WHILE(ASSOCIATED(p_rcurrentLevel))
    
      IF (ASSOCIATED(p_rcurrentLevel%p_rprevLevel)) THEN
        CALL lsysbl_createVecFromScalar (&
            rsolverNode%p_rsubnodeMultigrid%rprjTempVector,&
            p_rcurrentLevel%rcgcorrTempVector)
      
        CALL lsysbl_enforceStructure (&
            p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rcgcorrTempVector)
      END IF
    
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
    END DO
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
    CALL sys_halt()
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rpreSmoother,isubgroup,ierror)
      IF (ierror .NE. LINSOL_ERR_NOERROR) RETURN
    END IF
    
    ! Pre- and postsmoother may be identical!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rpostSmoother,isubgroup,ierror)
      IF (ierror .NE. LINSOL_ERR_NOERROR) RETURN
    END IF
    
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup,ierror)
      IF (ierror .NE. LINSOL_ERR_NOERROR) RETURN
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
    CALL sys_halt()
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
    CALL sys_halt()
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
  
    ! Release the coarse grid solver
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    
    ! Pre- and postsmoother may be identical! If not, first release the 
    ! postsmoother.
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    
    ! Release the presmoother
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

      ! Release the temp vector for the coarse grid correction.
      ! The vector shares its memory with the 'global' projection
      ! vector, so only release it if the other vector exists.
      IF (rsolverNode%p_rsubnodeMultigrid%rprjTempVector%NEQ .GT. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rcgcorrTempVector)
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
    CALL sys_halt()
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
  
  ! Release the coarse grid correction structure
  CALL cgcor_release (rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection)
  
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
    INTEGER :: iiterations
    REAL(DP) :: dres
    !DEBUG: REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
    
    ! Cancel if nmaxIterations = number of smoothing steps is =0.
    IF (rsolverNode%nmaxIterations .LE. 0) RETURN
    
    ! Some solvers can be used without the usual preconditioner approach to gain speed.
    ! First check if we have such a solver; these won't need the additional matrix
    ! vector multiplication below!
    
    SELECT CASE (rsolverNode%calgorithm)
      ! We are in a case where we can apply an adaptive 1-step defect-correction
      ! type preconditioner as smoother. This is a very special case and applies
      ! only to some special algorithms, which work directly with a solution-
      ! and RHS vector, like special VANCA variants. 
      !
      ! The nice thing: This saves one matrix-vector multiplication and 
      ! a lot of time (about 30-40%)!!!
      !
      ! The disadvantage: We cannot use any filtering with this approach,
      ! since the filters only apply to defect vectors and not solution vectors!
      !
      ! But in most situations, this disadvantage does not have such a large
      ! impact, as the MG solver applies filtering frequently to all the defect
      ! vectors that pop up during its iteration. 
      ! Or the algorithm filters internally, so we don't have to take into
      ! account the filtering at all.
    
    CASE (LINSOL_ALG_VANCA) 
      ! The basic proceeding is like below: Apply the corresponding solver multiple
      ! times to the solution- and RHS-vector to improve the solution. Let's see
      ! which solver we have. 
      ! If the previous IF-guess was wrong, we will leave this IF-clause and
      ! fall back to use the preconditioner approach.
      
      !DEBUG: CALL lsysbl_getbase_double (rx,p_Ddata)
      !DEBUG: CALL lsysbl_getbase_double (rb,p_Ddata2)
          
      SELECT CASE (rsolverNode%p_rsubnodeVANCA%csubtypeVANCA)
      CASE (LINSOL_VANCA_GENERALDIRECT,&
            LINSOL_VANCA_2DNAVSTDIRECT,&
            LINSOL_VANCA_2DFNAVSTDIRECT,&
            LINSOL_VANCA_2DFNAVSTOCDIRECT,&
            LINSOL_VANCA_2DFNAVSTOCDIAGDIR,&
            LINSOL_VANCA_3DNAVSTDIRECT,&
            LINSOL_VANCA_3DFNAVSTDIRECT)
        ! Yes, this solver can be applied to a given solution/rhs vector directly.
        ! Call it nmaxIterations times to perform the smoothing.
        DO i=1,rsolverNode%nmaxIterations
        
          ! Probably print the residuum
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL lsysbl_copyVector(rb,rtemp)
            CALL lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
            
            dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
            IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                      
            CALL output_line ('Smoother: Step '//TRIM(sys_siL(i-1,10))//&
                ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
          END IF
        
          ! Perform nmaxIterations:   x_n+1 = x_n + C^{-1} (b-Ax_n)
          ! without explicitly calculating (b-Ax_n) like below.
          CALL vanca_conformal (rsolverNode%p_rsubnodeVANCA%rvanca, &
                                rx, rb, rsolverNode%domega)
                                
        END DO

        ! Probably print the residuum
        IF (rsolverNode%ioutputLevel .GE. 2) THEN
          CALL lsysbl_copyVector(rb,rtemp)
          CALL lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
          
          dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
          IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                    
          CALL output_line ('Smoother: Step '//TRIM(sys_siL(i-1,10))//&
              ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
        END IF

        ! That's it.
        RETURN
      
      END SELECT
      
    END SELECT
      
    ! This is a 1-step solver, we have to emulate the smoothing
    ! iterations. Perform rsolverNode%nmaxIterations steps of the
    ! form
    !     $$ x_{n+1} = x_n + P^{-1}(b-Ax_n) $$
    ! with $x_0 = 0$.
    
    bfilter = ASSOCIATED(p_RfilterChain)
    
    !DEBUG: CALL lsysbl_getbase_double (rx,p_Ddata)
    !DEBUG: CALL lsysbl_getbase_double (rtemp,p_Ddata2)
    
    ! Do we have an iterative or one-step solver given?
    ! A 1-step solver performs the following loop nmaxIterations times, while an iterative
    ! solver is called only once and performs nmaxIterations steps internally.
    ! (This is a convention. Calling an iterative solver i times with j internal steps
    ! would also be possible, but we don't implement that here.)
    IF (IAND(rsolverNode%ccapability,LINSOL_ABIL_DIRECT) .NE. 0) THEN
      iiterations = rsolverNode%nmaxIterations
    ELSE
      iiterations = 1
    END IF
    
    DO i=1,iiterations
    
      CALL lsysbl_copyVector(rb,rtemp)
      CALL lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
        IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                  
        CALL output_line ('Smoother: Step '//TRIM(sys_siL(i-1,10))//&
            ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
      END IF
      
      ! Apply the filter to this defect before preconditioning
      IF (bfilter) THEN
        CALL filter_applyFilterChainVec (rtemp, p_RfilterChain)
      END IF
      
      ! Perform preconditioning
      CALL linsol_precondDefect(rsolverNode,rtemp)
      
      ! If the preconditioner broke down, cancel the smoothing,
      ! it would destroy out solution!
      IF (rsolverNode%iresult .NE. 0) THEN
        IF (rsolverNode%ioutputLevel .GE. 0) THEN
          CALL output_line ('Smoothing canceled, preconditioner broke down!', &
                            OU_CLASS_WARNING,OU_MODE_STD,'linsol_smoothCorrection')
          EXIT
        END IF
      END IF
      
      CALL lsysbl_vectorLinearComb (rtemp,rx,1.0_DP,1.0_DP)
      
    END DO

    ! Probably print the final residuum
    IF (rsolverNode%ioutputLevel .GE. 2) THEN
      CALL lsysbl_copyVector(rb,rtemp)
      CALL lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
      
      dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
      IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                
      CALL output_line ('Smoother: Step '//TRIM(sys_siL(i-1,10))//&
          ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
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
  INTEGER :: nminIterations,nmaxIterations,ireslength,niteResOutput,niteAsymptoticCVR
  INTEGER :: ite,ilev,nlmax,i,j,nblocks
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
    ireslength = 32 !MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))
    niteAsymptoticCVR = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

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
        CALL output_line ('Multigrid: Only one level. '//&
             'Switching back to standard solver.')
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
            p_rcurrentLevel%ncycles = 2
          ELSE
            p_rcurrentLevel%ncycles = p_rsubnode%icycle
          END IF  
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
        END DO
        
        ! Afterwards, p_rcurrentLevel points to the maximum level again.
        ! There, we set ncycles to 1.
        p_rcurrentLevel%ncycles = 1
        
        ! Print out the initial residuum

        IF (rsolverNode%ioutputLevel .GE. 2) THEN
          CALL output_line ('Multigrid: Iteration '// &
              TRIM(sys_siL(0,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
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
            p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles
            p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          END DO
          ! Don't forget the highest level.
          p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles
        
          ! p_rcurrentLevel now points to the maximum level.
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
              
              ! Extended output
              IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother) .AND. &
                  (rsolverNode%ioutputLevel .GE. 3) .AND. &
                  (MOD(ite,niteResOutput).EQ.0)) THEN
                  
                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                    ' after presm.:     !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
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
                
                ! Extended output and/or adaptive cycles
                IF ((rsolverNode%p_rsubnodeMultigrid%depsRelCycle .NE. 1E99_DP) .OR.&
                    (rsolverNode%ioutputLevel .GE. 3)) THEN
                
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rrhsVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  ! In case adaptive cycles are activated, save the 'initial' residual
                  ! of that level into the level structure. Then we can later check
                  ! if we have to repeat the cycle on the coarse mesh.
                  IF (rsolverNode%p_rsubnodeMultigrid%depsRelCycle .NE. 1E99_DP) THEN
                    p_rlowerLevel%dinitResCycle = dres
                    p_rlowerLevel%icycleCount = 1
                  END IF
                         
                  ! If the output level is high enough, print that residuum norm.   
                  IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                      (MOD(ite,niteResOutput).EQ.0)) THEN
                    CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
                        ' after restrict.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                  END IF
                  
                END IF

              ELSE
              
                ! The vector is to be restricted to the coarse grid.
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

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rsolutionVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
                      ' after restrict.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
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
                !
                ! Note that this shall not be done on the coarse grid as there is
                ! no temp/rhs vector!
                IF (ASSOCIATED(p_rlowerLevel%p_rprevLevel)) THEN
                  p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                    -ABS(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)% &
                    isortStrategy) 
                  p_rlowerLevel%rrhsVector%RvectorBlock(1:nblocks)% &
                    isortStrategy = -ABS(p_rlowerLevel%rrhsVector% &
                    RvectorBlock(1:nblocks)%isortStrategy) 
                END IF
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
              ! Use either the standard system matrix for this purpose or
              ! the specific matrix for the coarse grid correction if this exists
              ! in the level info structure.
              IF (p_rcurrentLevel%rsystemMatrixOptCorrection%NEQ .NE. 0) THEN
                CALL cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrixOptCorrection,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              ELSE
                CALL cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrix,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              END IF
              
              ! Perform the coarse grid correction by adding the coarse grid
              ! solution (with the calculated step-length parameter) to
              ! the current solution.
              
              CALL lsysbl_vectorLinearComb (p_rcurrentLevel%rtempVector,&
                                            p_rcurrentLevel%rsolutionVector,&
                                            dstep,1.0_DP)
                            
              ! Extended output
              IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                  (MOD(ite,niteResOutput).EQ.0)) THEN
                  
                CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                        p_rcurrentLevel%rtempVector)
                CALL lsysbl_blockMatVec (&
                    p_rcurrentLevel%rsystemMatrix, &
                    p_rcurrentLevel%rsolutionVector,&
                    p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                    ' after c.g.corr.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
              END IF
                                            
              ! Perform the post-smoothing with the current solution vector
              IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
                CALL linsol_smoothCorrection (p_rcurrentLevel%p_rpostSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                  CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector)
                  CALL lsysbl_blockMatVec (&
                      p_rcurrentLevel%rsystemMatrix, &
                      p_rcurrentLevel%rsolutionVector,&
                      p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                  dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                      ' after postsm.:    !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                END IF
                                              
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
                  ! Cycle finished. Reset counter for next cycle.
                  p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles

                  IF ((rsolverNode%p_rsubnodeMultigrid%depsRelCycle .NE. 1E99_DP) .AND. &
                      ASSOCIATED(p_rcurrentLevel%p_rnextLevel)) THEN
                      
                    ! Adaptive cycles activated. 
                    !
                    ! We are on a level < nlmax.
                    ! At first, calculate the residuum on that level.
                    CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                            p_rcurrentLevel%rtempVector)
                    CALL lsysbl_blockMatVec (&
                        p_rcurrentLevel%rsystemMatrix, &
                        p_rcurrentLevel%rsolutionVector,&
                        p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                    dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                        rsolverNode%iresNorm)
                    IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                              (dres .LE. 1E99_DP))) dres = 0.0_DP
                              
                    ! Compare it with the initial residuum. If it's not small enough
                    ! and if we haven't reached the maximum number of cycles,
                    ! repeat the complete cycle.
                    IF ( ((rsolverNode%p_rsubnodeMultigrid%nmaxAdaptiveCycles .LE. -1) &
                          .OR. &
                          (p_rcurrentLevel%icycleCount .LT. &
                           rsolverNode%p_rsubnodeMultigrid%nmaxAdaptiveCycles)) &
                        .AND. &
                        (dres .GT. rsolverNode%p_rsubnodeMultigrid%depsRelCycle * &
                                    p_rcurrentLevel%dinitResCycle) ) THEN

                      IF (rsolverNode%ioutputLevel .GE. 3) THEN
                        CALL output_line ( &
                          TRIM( &
                          sys_siL(p_rcurrentLevel%icycleCount,10)) &
                          //'''th repetition of cycle on level '// &
                          TRIM(sys_siL(ilev,5))//'.')
                      END IF

                      p_rcurrentLevel%icycleCount = p_rcurrentLevel%icycleCount+1
                      CYCLE cycleloop
                    END IF
                    
                    ! Otherwise: The cycle(s) is/are finished; 
                    ! the END DO goes up one level.
                    
                  END IF

                END IF
              ELSE

                ! Next cycle; go down starting from the current level
                CYCLE cycleloop

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
          IF (rsolverNode%ioutputLevel .LT. 2) THEN
            DO i=MAX(1,SIZE(Dresqueue)-ite-1+1),SIZE(Dresqueue)
              j = ITE-MAX(1,SIZE(Dresqueue)-ite+1)+i
              CALL output_line ('Multigrid: Iteration '// &
                  TRIM(sys_siL(j,10))//',  !!RES!! = '//&
                  TRIM(sys_sdEL(Dresqueue(i),15)) )
            END DO
          END IF
            CALL output_line ('Multigrid: Solution diverging!')
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
            CALL output_line ('Multigrid: Iteration '// &
                TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
                TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          END IF
          
        END DO  ! ite
        
        ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
        ! completed

        IF (ite .GT. rsolverNode%nmaxIterations) THEN
          ! Warning if we didn't reach the convergence criterion.
          ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
          IF ((rsolverNode%ioutputLevel .GE. 0) .AND. &
              (rsolverNode%nmaxIterations .GT. rsolverNode%nminIterations)) THEN
            CALL output_line ('Multigrid: Accuracy warning: '//&
                'Solver did not reach the convergence criterion')
          END IF

          ite = rsolverNode%nmaxIterations
        END IF

        ! Finish - either with an error or if converged.
        ! Print the last residuum.

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
            (rsolverNode%iresult .GE. 0)) THEN
          CALL output_line ('Multigrid: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF
        
        ! Final number of iterations
        rsolverNode%iiterations = ite
      
      END IF

      ! Finally, we look at some statistics:
      !
      ! Don't calculate anything if the final residuum is out of bounds -
      ! would result in NaN's,...
        
      IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
        ! Calculate asymptotic convergence rate
      
      IF (niteAsymptoticCVR .NE. 0) THEN
        I = MAX(33-ite,33-niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(33-I,DP))
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
        CALL output_lbrk()
        CALL output_line ('Multigrid statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'MG: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF

    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initMultigrid2 (p_rsolverNode,nlevels,p_Rfilter)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the Multigrid solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  ! Before the solver can used for solving the problem, the caller must add
  ! information about all levels (matrices,...) to the solver. 
  ! This can be done by linsol_setMultigridLevel.
  
!</description>
  
!<input>

  ! Number of levels supported by this solver.
  INTEGER, INTENT(IN) :: nlevels
  
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
  p_rsolverNode%calgorithm = LINSOL_ALG_MULTIGRID2
  
  ! Initialise the ability bitfield with the ability of this solver:
  p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR     + LINSOL_ABIL_BLOCK        + &
                              LINSOL_ABIL_MULTILEVEL + LINSOL_ABIL_CHECKDEF     + &
                              LINSOL_ABIL_USESUBSOLVER + &
                              LINSOL_ABIL_USEFILTER
  
  ! Allocate the subnode for Multigrid.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(p_rsolverNode%p_rsubnodeMultigrid2)
  
  ! Initialise the coarse grid correction structure.1
  CALL cgcor_init (p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
  
  ! Allocate level information data
  ALLOCATE(p_rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(nlevels))
  
  ! Attach the filter if given. 
  IF (PRESENT(p_Rfilter)) THEN
    p_rsolverNode%p_rsubnodeMultigrid2%p_RfilterChain => p_Rfilter
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_alterMultigrid2 (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the Multigrid solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration 
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT), TARGET             :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  TYPE(t_linsol_alterSolverConfig), INTENT(INOUT)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i
    TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel

    ! Pass the command structure to all subsolvers and smoothers
    ! on all levels.
    DO i=1,SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(i)
      
      IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
        CALL linsol_alterSolver (p_rcurrentLevel%p_rpreSmoother, &
                                 ralterConfig)
      END IF
      ! Pre- and postsmoother may be identical!
      ! Take care not to alter the same smoother twice!
      IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
          (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) THEN
        CALL linsol_alterSolver (p_rcurrentLevel%p_rpostSmoother, &
                                 ralterConfig)
      END IF
      IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
        CALL linsol_alterSolver (p_rcurrentLevel%p_rcoarseGridSolver, &
                                 ralterConfig)
      END IF
      
    END DO
    
  END SUBROUTINE

! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixMultigrid2 (rsolverNode,Rmatrices)
  
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
  TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel
  INTEGER                       :: ilevel,nlmax
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
  END IF
  
  ! Make sure we have the right amount of matrices
  IF (SIZE(Rmatrices) .NE. SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)) THEN
    PRINT *,'Error: Wrong number of matrices'
    CALL sys_halt()
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. If yes, attach the matrix
  ! to it.
  ! For this purpose, call the general initialisation routine, but
  ! pass only that part of the Rmatrices array that belongs
  ! to the range of levels between the coarse grid and the current grid.
  
  nlmax = SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)
  DO ilevel = 1,nlmax
  
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rpreSmoother, &
                                Rmatrices(1:ilevel))
    END IF
    ! Pre- and postsmoother may be identical!
    ! Take care not to initialise the same smoother twice!
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rpostSmoother, &
                                Rmatrices(1:ilevel))
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rcoarseGridSolver, &
                                Rmatrices(1:ilevel))
    END IF
    
    ! Write a link to the matrix into the level-structure.
    CALL lsysbl_duplicateMatrix (Rmatrices(ilevel), &
        p_rcurrentLevel%rsystemMatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_matCompatMultigrid2 (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  TYPE(t_linsolNode), INTENT(IN)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  ! More precisely, ccompatible returns always the status of the highest
  ! level where there is an error -- or LINSOL_COMP_OK if there is no error.
  INTEGER, INTENT(OUT) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! local variables
    TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel
    INTEGER                       :: ilevel,ccompatLevel
    
    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK
    
    ! Reset all compatibility flags
    IF (PRESENT(CcompatibleDetail)) CcompatibleDetail (:) = ccompatible    

    ! Make sure the solver node is configured for multigrid
    IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
        (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
      PRINT *,'Error: Multigrid structure not initialised'
      CALL sys_halt()
    END IF
    
    ! Make sure we have the right amount of matrices
    IF (SIZE(Rmatrices) .NE. SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)) THEN
      PRINT *,'Error: Wrong number of matrices'
      CALL sys_halt()
    END IF

    ! Check for every level, if there's a presmoother, postsmoother or
    ! coarse grid solver structure attached. If yes, it must check
    ! the matrices!
    ! For this purpose, call the general check routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    
    DO ilevel = 1,SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    
      ! Compatibility flag for that level
      ccompatLevel = LINSOL_COMP_OK
    
      ! Check presmoother
      IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      
        IF (PRESENT(CcompatibleDetail)) THEN
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(1:ilevel),ccompatLevel,CcompatibleDetail(1:ilevel))
        ELSE
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(1:ilevel),ccompatLevel)
        END IF
      
        IF (ccompatLevel .NE. LINSOL_COMP_OK) THEN
          ccompatible = ccompatLevel
          CYCLE
        END IF
        
      END IF
      
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
          (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) THEN
        
        IF (PRESENT(CcompatibleDetail)) THEN
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(1:ilevel),ccompatLevel,CcompatibleDetail(1:ilevel))
        ELSE
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(1:ilevel),ccompatLevel)
        END IF
        
        IF (ccompatLevel .NE. LINSOL_COMP_OK) THEN
          ccompatible = ccompatLevel
          CYCLE
        END IF
        
      END IF
      
      ! Check coarse grid solver
      IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
        
        IF (PRESENT(CcompatibleDetail)) THEN
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(1:ilevel),ccompatLevel,CcompatibleDetail(1:ilevel))
        ELSE
          CALL linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(1:ilevel),ccompatLevel)
        END IF
        
        IF (ccompatLevel .NE. LINSOL_COMP_OK) THEN
          ccompatible = ccompatLevel
          CYCLE
        END IF
        
      END IF
      
    END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureMultigrid2 (rsolverNode,ierror,isolverSubgroup)
  
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
  TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel,p_rprevLevel
  INTEGER :: isubgroup,nlmax,ilevel
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
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
  END IF

  nlmax = SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  DO ilevel = 1,nlmax
  
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

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

  END DO
  
  ! Cancel here, if we don't belong to the subgroup to be initialised
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

  imemmax = 0

  ! Multigrid needs temporary vectors for the RHS, solution and general
  ! data. Memory for that is allocated on all levels for temporary, RHS and 
  ! solution vectors.
  DO ilevel = 1,nlmax
  
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

    p_rmatrix => p_rcurrentLevel%rsystemMatrix
    NULLIFY(p_rtemplVect)
    IF (ilevel .LT. nlmax) THEN
      ! On the maximum level we don't need additional memory for the 
      ! solution vector - as solution vector on the fine grid, 
      ! the vector which comes as parameter
      ! to the multigrid  preconditioner is used.
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rsolutionVector,.FALSE.,.FALSE.,&
            rsolverNode%cdefaultDataType)
      p_rtemplVect => p_rcurrentLevel%rsolutionVector
    END IF
    ! On the coarse grid, we don't need memory for temporary/RHS
    ! vectors. The RHS on the coarse grid is replaced in-situ
    ! by the solution vector.
    IF (ilevel .GT. 1) THEN
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rrhsVector,.FALSE.,.FALSE.,&
            rsolverNode%cdefaultDataType)
      CALL lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rtempVector,.FALSE.,.FALSE.,&
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
    IF (ilevel .GT. 1) THEN ! otherwise coarse grid
      p_rprevLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel-1)
   
      ! Calculate the memory that is necessary for prolongation/restriction -
      ! in case there is temporary memory needed.
      ! The system matrix on the fine/coarse grid specifies the discretisation.
      imemmax = MAX(imemmax,mlprj_getTempMemoryMat ( &
                    p_rcurrentLevel%rinterlevelProjection, &
                    p_rprevLevel%rsystemMatrix,&
                    p_rcurrentLevel%rsystemMatrix))
    END IF

    ! All temporary vectors are marked as 'unsorted'. We can set the
    ! sorting flag directly here without using the resorting routines
    ! as the vectors have just been created.
    !
    ! Don't do anything to the RHS/temp vectors on the coarse grid
    ! and the solution vector on the fine grid -- they don't exist!
    IF (ilevel .GT. 1) THEN
      p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy = &
        -ABS(p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy)
      p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy = &
        -ABS(p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy)
    END IF
    
    IF (ilevel .LT. NLMAX) THEN
      p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy = &
        -ABS(p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy)
    END IF
    
    ! And the next level...
  END DO
  
  ! Do we have to allocate memory for prolongation/restriction/resorting?
  IF (imemmax .GT. 0) THEN
    CALL lsyssc_createVector (rsolverNode%p_rsubnodeMultigrid2%rprjTempVector,&
                              imemmax,.FALSE.,rsolverNode%cdefaultDataType)
                              
    ! Use this vector also for the coarse grid correction.
    ! Create a block vector on every level (except for the coarse grid)
    ! with the structure of the RHS that shares the memory with rprjTempVector. 
    DO ilevel = 2,nlmax
  
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    
      CALL lsysbl_createVecFromScalar (&
          rsolverNode%p_rsubnodeMultigrid2%rprjTempVector,&
          p_rcurrentLevel%rcgcorrTempVector)
    
      CALL lsysbl_enforceStructure (&
          p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rcgcorrTempVector)

    END DO
  ELSE
    rsolverNode%p_rsubnodeMultigrid2%rprjTempVector%NEQ = 0
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initDataMultigrid2 (rsolverNode,ierror,isolverSubgroup)
  
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
  TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup,ilevel
  
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
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  DO ilevel = 1,SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)
  
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

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

  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneDataMultigrid2 (rsolverNode,isolverSubgroup)
  
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
  TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup,ilevel
  
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
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  DO ilevel = SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo),1,-1
  
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    
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
    
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructureMultigrid2 (rsolverNode,isolverSubgroup)
  
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
  TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup,ilevel,nlmax
  
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
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
  END IF
  
  nlmax = SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  DO ilevel = 1,nlmax
  
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
  
    ! Release the coarse grid solver
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    
    ! Pre- and postsmoother may be identical! If not, first release the 
    ! postsmoother.
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother) .AND. &
        (.NOT. ASSOCIATED(p_rcurrentLevel%p_rpreSmoother, &
                          p_rcurrentLevel%p_rpostSmoother))) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    
    ! Release the presmoother
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    
  END DO

  ! Cancel here, if we don't belong to the subgroup to be initialised
  IF (isubgroup .NE. rsolverNode%isolverSubgroup) RETURN

  ! Release temporary memory.
  DO ilevel = 1,nlmax
  
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
  
    IF (ilevel .GT. 1) THEN
    
      IF (p_rcurrentLevel%rtempVector%NEQ .NE. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rtempVector)
      END IF
      
      IF (p_rcurrentLevel%rrhsVector%NEQ .NE. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rrhsVector)
      END IF

      ! Release the temp vector for the coarse grid correction.
      ! The vector shares its memory with the 'global' projection
      ! vector, so only release it if the other vector exists.
      IF (rsolverNode%p_rsubnodeMultigrid2%rprjTempVector%NEQ .GT. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rcgcorrTempVector)
      END IF

    END IF

    IF (ilevel .LT. nlmax) THEN
      IF (p_rcurrentLevel%rsolutionVector%NEQ .NE. 0) THEN
        CALL lsysbl_releaseVector (p_rcurrentLevel%rsolutionVector)
      END IF
    END IF

  END DO

  ! Check if we have to release our temporary vector for prolongation/
  ! restriction.
  ! Do we have to allocate memory for prolongation/restriction?
  IF (rsolverNode%p_rsubnodeMultigrid2%rprjTempVector%NEQ .GT. 0) THEN
    CALL lsyssc_releaseVector (rsolverNode%p_rsubnodeMultigrid2%rprjTempVector)
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneMultigrid2 (rsolverNode)
  
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
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID2) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    CALL sys_halt()
  END IF

  ! Release data and structure if this has not happened yet.
  CALL linsol_doneDataMultigrid2(rsolverNode)
  CALL linsol_doneStructureMultigrid2(rsolverNode)
  
  ! Remove all the levels
  CALL linsol_cleanMultigrid2Levels (rsolverNode)
  DEALLOCATE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)
  
  ! Release the coarse grid correction structure
  CALL cgcor_release (rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
  
  ! Release MG substructure, that's it.
  DEALLOCATE(rsolverNode%p_rsubnodeMultigrid2)

  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_precMultigrid2 (rsolverNode,rd)
  
!<description>
  ! Applies Multigrid preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The structure of the levels (pre/postsmoothers, interlevel projection
  ! structures) must have been prepared by setting the level parameters
  ! before calling this routine. (The parameters of a specific level
  ! can be get by calling linsol_getMultigridLevel2).
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
  INTEGER :: nminIterations,nmaxIterations,ireslength,niteResOutput,niteAsymptoticCVR
  INTEGER :: ite,ilev,nlmax,i,j,nblocks
  REAL(DP) :: dres,dstep
  LOGICAL :: bfilter,bsort
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
  ! The queue saves the current residual and the two previous residuals.
  REAL(DP), DIMENSION(32) :: Dresqueue
  
  ! The system matrix on the current level
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  
  ! Our MG structure
  TYPE(t_linsolSubnodeMultigrid2), POINTER :: p_rsubnode
  
  ! The current level and the next lower one.
  TYPE(t_linsolMGLevelInfo2), POINTER :: p_rcurrentLevel,p_rlowerLevel
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeMultigrid2
    
    bfilter = ASSOCIATED(p_rsubnode%p_RfilterChain)
    p_RfilterChain => p_rsubnode%p_RfilterChain
    
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
    
    IF (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)) THEN
      PRINT *,'Multigrid: No levels attached!'
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Maximum level
    nlmax = SIZE(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    ireslength = 32 !MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))
    niteAsymptoticCVR = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

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
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(nlmax)
    
    ! Is there only one level? Can be seen if the current level
    ! already contains a coarse grid solver.
    IF (nlmax .EQ. 1) THEN
    
      IF (rsolverNode%ioutputLevel .GT. 1) THEN
        CALL output_line ('Multigrid: Only one level. '//&
             'Switching back to standard solver.')
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
      ! True multigrid.
      !    
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
        DO i=1,nlmax-1
          p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(i)
          IF (p_rsubnode%icycle .EQ. 0) THEN
            p_rcurrentLevel%ncycles = 2
          ELSE
            p_rcurrentLevel%ncycles = p_rsubnode%icycle
          END IF  
        END DO
        
        ! We start at the maximum level.        
        ilev = nlmax

        ! Get current and next lower level.
        p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev)
        p_rlowerLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev-1)
               
        ! On the current (max.) level we set ncycles to 1.
        p_rcurrentLevel%ncycles = 1
        
        ! Print out the initial residuum

        IF (rsolverNode%ioutputLevel .GE. 2) THEN
          CALL output_line ('Multigrid: Iteration '// &
              TRIM(sys_siL(0,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        END IF

        ! Copy the initial RHS to the RHS vector on the maximum level.
        CALL lsysbl_copyVector(rd,p_rcurrentLevel%rrhsVector)
        
        ! Replace the solution vector on the finest level by rd.
        ! Afterwards, rd and the solution vector on the finest level
        ! share the same memory location, so we use rd as iteration
        ! vector directly.
        CALL lsysbl_duplicateVector (rd,p_rcurrentLevel%rsolutionVector,&
            LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
        
        ! Clear the initial solution vector.
        CALL lsysbl_clearVector (p_rcurrentLevel%rsolutionVector)
        
        ! Start multigrid iteration; perform at most nmaxiterations iterations.
        DO ite = 1, nmaxiterations
        
          rsolverNode%icurrentIteration = ite
          
          ! Initialize cycle counters for all levels.
          DO i=1,nlmax
            p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(i)
            p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles
          END DO
        
          ! p_rcurrentLevel now points to the maximum level again!
          !
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
            
            DO WHILE (ilev .GT. 1)
            
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
              
              ! Extended output
              IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother) .AND. &
                  (rsolverNode%ioutputLevel .GE. 3) .AND. &
                  (MOD(ite,niteResOutput).EQ.0)) THEN
                  
                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                    ' after presm.:     !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
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
              IF (ilev .GT. 2) THEN
              
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
                
                ! Extended output and/or adaptive cycles
                IF ((rsolverNode%p_rsubnodeMultigrid2%depsRelCycle .NE. 1E99_DP) .OR.&
                    (rsolverNode%ioutputLevel .GE. 3)) THEN
                
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rrhsVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  ! In case adaptive cycles are activated, save the 'initial' residual
                  ! of that level into the level structure. Then we can later check
                  ! if we have to repeat the cycle on the coarse mesh.
                  IF (rsolverNode%p_rsubnodeMultigrid2%depsRelCycle .NE. 1E99_DP) THEN
                    p_rlowerLevel%dinitResCycle = dres
                    p_rlowerLevel%icycleCount = 1
                  END IF
                         
                  ! If the output level is high enough, print that residuum norm.   
                  IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                      (MOD(ite,niteResOutput).EQ.0)) THEN
                    CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
                        ' after restrict.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                  END IF
                  
                END IF

              ELSE
              
                ! The vector is to be restricted to the coarse grid.
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

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rsolutionVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
                      ' after restrict.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                END IF

              END IF              
            
              ! Go down one level
              ilev = ilev - 1
              p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev)
              
              IF (ilev .NE. 1) &
                p_rlowerLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev-1)
              
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
            DO WHILE(ilev .LT. nlmax)
            
              ilev = ilev + 1
              p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev)
              p_rlowerLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev-1)
              
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
                !
                ! Note that this shall not be done on the coarse grid as there is
                ! no temp/rhs vector!
                IF (ilev .GT. 2) THEN
                  p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                    -ABS(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)% &
                    isortStrategy) 
                  p_rlowerLevel%rrhsVector%RvectorBlock(1:nblocks)% &
                    isortStrategy = -ABS(p_rlowerLevel%rrhsVector% &
                    RvectorBlock(1:nblocks)%isortStrategy) 
                END IF
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
              ! Use either the standard system matrix for this purpose or
              ! the specific matrix for the coarse grid correction if this exists
              ! in the level info structure.
              IF (p_rcurrentLevel%rsystemMatrixOptCorrection%NEQ .NE. 0) THEN
                CALL cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrixOptCorrection,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              ELSE
                CALL cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrix,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              END IF
              
              ! Perform the coarse grid correction by adding the coarse grid
              ! solution (with the calculated step-length parameter) to
              ! the current solution.
              
              CALL lsysbl_vectorLinearComb (p_rcurrentLevel%rtempVector,&
                                            p_rcurrentLevel%rsolutionVector,&
                                            dstep,1.0_DP)
                            
              ! Extended output
              IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                  (MOD(ite,niteResOutput).EQ.0)) THEN
                  
                CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                        p_rcurrentLevel%rtempVector)
                CALL lsysbl_blockMatVec (&
                    p_rcurrentLevel%rsystemMatrix, &
                    p_rcurrentLevel%rsolutionVector,&
                    p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                    ' after c.g.corr.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
              END IF
                                            
              ! Perform the post-smoothing with the current solution vector
              IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
                CALL linsol_smoothCorrection (p_rcurrentLevel%p_rpostSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                  CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector)
                  CALL lsysbl_blockMatVec (&
                      p_rcurrentLevel%rsystemMatrix, &
                      p_rcurrentLevel%rsolutionVector,&
                      p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                  dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                      ' after postsm.:    !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                END IF
                                              
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
                  ! Cycle finished. Reset counter for next cycle.
                  p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles

                  IF ((rsolverNode%p_rsubnodeMultigrid2%depsRelCycle .NE. 1E99_DP) .AND. &
                      (ilev .LT. NLMAX)) THEN
                      
                    ! Adaptive cycles activated. 
                    !
                    ! We are on a level < nlmax.
                    ! At first, calculate the residuum on that level.
                    CALL lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                            p_rcurrentLevel%rtempVector)
                    CALL lsysbl_blockMatVec (&
                        p_rcurrentLevel%rsystemMatrix, &
                        p_rcurrentLevel%rsolutionVector,&
                        p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                    dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                        rsolverNode%iresNorm)
                    IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                              (dres .LE. 1E99_DP))) dres = 0.0_DP
                              
                    ! Compare it with the initial residuum. If it's not small enough
                    ! and if we haven't reached the maximum number of cycles,
                    ! repeat the complete cycle.
                    IF ( ((rsolverNode%p_rsubnodeMultigrid2%nmaxAdaptiveCycles .LE. -1) &
                          .OR. &
                          (p_rcurrentLevel%icycleCount .LT. &
                           rsolverNode%p_rsubnodeMultigrid2%nmaxAdaptiveCycles)) &
                        .AND. &
                        (dres .GT. rsolverNode%p_rsubnodeMultigrid2%depsRelCycle * &
                                    p_rcurrentLevel%dinitResCycle) ) THEN

                      IF (rsolverNode%ioutputLevel .GE. 3) THEN
                        CALL output_line ( &
                          TRIM( &
                          sys_siL(p_rcurrentLevel%icycleCount,10)) &
                          //'''th repetition of cycle on level '// &
                          TRIM(sys_siL(ilev,5))//'.')
                      END IF

                      p_rcurrentLevel%icycleCount = p_rcurrentLevel%icycleCount+1
                      CYCLE cycleloop
                    END IF
                    
                    ! Otherwise: The cycle(s) is/are finished; 
                    ! the END DO goes up one level.
                    
                  END IF

                END IF
              ELSE

                ! Next cycle; go down starting from the current level
                CYCLE cycleloop

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
          IF (rsolverNode%ioutputLevel .LT. 2) THEN
            DO i=MAX(1,SIZE(Dresqueue)-ite-1+1),SIZE(Dresqueue)
              j = ITE-MAX(1,SIZE(Dresqueue)-ite+1)+i
              CALL output_line ('Multigrid: Iteration '// &
                  TRIM(sys_siL(j,10))//',  !!RES!! = '//&
                  TRIM(sys_sdEL(Dresqueue(i),15)) )
            END DO
          END IF
            CALL output_line ('Multigrid: Solution diverging!')
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
            CALL output_line ('Multigrid: Iteration '// &
                TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
                TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          END IF
          
        END DO  ! ite
        
        ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
        ! completed

        IF (ite .GT. rsolverNode%nmaxIterations) THEN
          ! Warning if we didn't reach the convergence criterion.
          ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
          IF ((rsolverNode%ioutputLevel .GE. 0) .AND. &
              (rsolverNode%nmaxIterations .GT. rsolverNode%nminIterations)) THEN
            CALL output_line ('Multigrid: Accuracy warning: '//&
                'Solver did not reach the convergence criterion')
          END IF

          ite = rsolverNode%nmaxIterations
        END IF

        ! Finish - either with an error or if converged.
        ! Print the last residuum.

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
            (rsolverNode%iresult .GE. 0)) THEN
          CALL output_line ('Multigrid: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF
        
        ! Final number of iterations
        rsolverNode%iiterations = ite
      
      END IF

      ! Finally, we look at some statistics:
      !
      ! Don't calculate anything if the final residuum is out of bounds -
      ! would result in NaN's,...
        
      IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
        ! Calculate asymptotic convergence rate
      
        IF (niteAsymptoticCVR .NE. 0) THEN
          I = MAX(33-ite,33-niteAsymptoticCVR)
          rsolverNode%dasymptoticConvergenceRate = &
            (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(33-I,DP))
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
      
      ! Release the temporary copy of the solution vector
      CALL lsysbl_releaseVector (p_rcurrentLevel%rsolutionVector)
    
    END IF

    ! As the solution vector on the finest level shared its memory with rd,
    ! we just calculated the new correction vector!
      
    IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        CALL output_lbrk()
        CALL output_line ('Multigrid statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'MG: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF

    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
END MODULE
