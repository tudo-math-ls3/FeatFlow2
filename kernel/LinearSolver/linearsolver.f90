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
!#         <tex> $$Ax = b$$ </tex>
!#     with A being an arbitrary block matrix
!# 2.) Saddle-point block matrices of the form
!#     <verb>
!#           [A  B1] (x) = (f1)
!#           [B2  0] (y)   (f2)
!#     </verb>
!#     with a regular block matrix A and general block matrices B1 and B2
!# 3.) Saddle-point-like block matrices of the form
!#     <verb>
!#           [A  B1] (x) = (f1)
!#           [B2  C] (y)   (f2)
!#     </verb>
!#     with a regular block matrix A, general block matrices B1 and B2
!#     and a 'stabilisation' matrix $C \sim 0$.
!#
!# Cases 2.) + 3.) are special cases of 1.) in that sense, that specialised
!# solvers are available for these kind of systems which allow faster solving.
!#
!# To solve a problem, one has basically to call the following routines
!# 1.) linsol_initXXXX                - Initialises a solver, returns a solver
!#                                      structure identifying the solver
!# 2.) linsol_setMatrix               - Attach the system matrix to the solver
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
!# For initialising a multigrid solver, this sequence changes a bit,
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
!#                                      *Not on the coarse grid!*
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
!# For initialising a multigrid solver, this sequence changes a bit,
!# since the data of multiple levels have to be associated to
!# the solver before invoking it:
!#
!# 1.) linsol_initMultigrid2          - Initialises multigrid solver, returns a
!#                                      solver structure identifying the solver
!# 2.) a) linsol_getMultigrid2Level   - Get a pointer to the level information
!#
!#     b) On the coarse grid:
!#          linsol_initXXXX           - Initialise coarse grid solver structure
!#
!#        On finer grids up to the maximum level:
!#          linsol_initXXXX
!#        + linsol_convertToSmoother  - Initialise a solver structure that can be
!#                                      attached as smoother to one level of multigrid.
!#                                      *Not on the coarse grid!*
!#     c) Add the pointer into the structure taken from a.). Write the pointer
!#        to there as coarse grid solver or smoother.
!#
!#     Steps a)-c) must be repeated for all levels. In linsol_getMultigrid2Level,
!#     counting always starts with level 1.
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
!# <!-- ------------------------------DEPRECATED START-------------------------------
!# For initialising a multigrid solver, this sequence changes a bit,
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
!#                                      *Not on the coarse grid!*
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
!# ------------------------------DEPRECATED STOP-------------------------------- -->
!# The following routines serve as auxiliary routines for the application to
!# maintain a Multigrid solver node:
!#
!# 1.) linsol_convertToSmoother
!#     -> Converts a solver node into a smoother for multigrid smoothing.
!#
!# 2.) linsol_alterSolver
!#     -> Allows direct modification of an existing solver. Passes a command
!#        block to all subsolvers in a solver tree.
!#
!# 3.) linsol_initProjMultigrid2Level
!#     -> Initialise the multigrid projection structure of one level
!#        with the standard projection or
!#     -> Initialise the multigrid projection structure of one level
!#        with a user defined projection
!#
!# 4.) linsol_initProjMultigrid2
!#     -> Initialise the multigrid projection structure on all uninitialised levels
!#        with the standard projectino
!#
!# The following routines serve as auxiliary routines for tghe MG2-solver.
!#
!# 1.) linsol_cleanMultigrid2Levels
!#     -> Deletes all levels and attached solver structures
!#
!# 2.) linsol_getMultigrid2Level
!#     -> Returns the level structure of a specific level
!#
!# 3.) linsol_getMultigrid2LevelCount
!#     -> Returns number of currently attached levels
!#
!# 4.) linsol_doneProjMultigrid2
!#     -> Release the projection information
!#
!#
!# Implementational details / structure of the solver library \\
!# ---------------------------------------------------------- \\
!# When going through this library, a newcomer would think:
!#
!#       'Where are the solvers? I can only find preconditioners!?!'
!#       'And what should that thing with the defect vectors mean?'
!#
!# The reason is simple: 'Everything is a preconditioner for defect vectors!'
!#
!# A short explaination: Let us assume we want to solve a linear system:
!#
!#                  <tex>$$ Ax=b $$</tex>
!#
!# For this we have a linear solver <tex>$P^{-1} \approx A^{-1}$</tex> (Gauss elimination
!# or whatever), so we get our solution vector x by saying:
!#
!#                 <tex>$$ x = P^{-1} b $$</tex>
!#
!# This is the usual way how to deal with a linear system in Numrics 1,
!# but it is a little bit hard to deal with in computer science. But we can
!# deal with the linear system more easily if we perform a simple
!# modification to this equation:
!#
!# <tex>                    $$ Ax = b $$
!# $$ \Rightarrow        P^{-1}Ax = P^{-1}b             $$
!# $$ \Rightarrow               0 = P^{-1} (b-Ax)       $$
!#
!# $$ \Rightarrow         x_{n+1} = x_n  +  P^{-1} (b-Ax) $$</tex>
!#
!# So the linear solver <tex>$P^{-1}$</tex> can equivalently be used in a defect
!# correction approach as a preconditioner for the defect $(b-Ax)$!
!# This can obviously be done for every linear solver.
!# So, this allows us to make a normalisation:
!#
!# ! All linear solvers can be formulated as a preconditioner to be applied
!#   to defect vectors, we do not have to take care of solution vectors !
!#
!# (and with the problems related to that like implementing boundary
!# conditions...)
!#
!# So everything that can be found in this library is a preconditioner.
!# The xxx_solve routines are all formulated that way, that they accept
!# the defect vector $d := (b-Ax)$ and overwrite it by the preconditioned
!# defect:
!#    <tex>    $$ d := P^{-1}d $$   </tex>
!#
!# which is actually the same as solving the linear system
!#
!#    <tex>    $$ Pd_{new} = d $$   </tex>
!#
!# with the right hand side $d$ being a defect vector given from outside.
!# Solving <tex>$ Pd_{new} = d $</tex> can then be done by an arbitrary linear solver,
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
!#  <tex>  $$  x_{new}  =  x + P^{-1}(b-Ax)  =  P^{-1}b  $$   </tex>
!#
!#  FAQ - Frequently asked Questions  \\
!# ---------------------------------- \\
!# Ok, the solver library seems to be a little bit complicated. Therefore
!# here a small chapter about how to set up specific tasks for newbies:
!#
!# 1.) I would like to set up a simple BiCGStab solver with ILU(0) preconditioner.
!#     How to do that?
!#
!#     Declare two solver structures, create the ILU(0) solver in the first,
!#     BiCGStab solver in the second and hang in the ILU(0) solver as
!#     preconditioner. Pass a filter chain in case you do not need it:
!#
!#     a) No filtering of the vectors:
!#
!# <verb>
!#      TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
!#
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner)
!# </verb>
!#
!#      -> p_rsolverNode is the solver identifying your BiCGStab(ILU(0))
!#         solver.
!#
!#     b) Filter the defect vector for Dirichlet boundary conditions
!#
!# <verb>
!#      TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
!#      TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
!#
!#      NULLIFY(p_rpreconditioner)
!#      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
!#      CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,RfilterChain)
!# </verb>
!#
!#      -> p_rsolverNode is the solver identifying your BiCGStab(ILU(0))
!#         solver.
!#
!# 2.) In numerics 1 I learned the pure Richardson iteration. Where can I
!#     find it?
!#
!#     Pure Richardson iteration (<tex>$x_{n+1} = x_n + \omega (b-Ax)$</tex>) is called
!#     'defect correction' here. Set up a defect correction solver and
!#     modify the domega parameter in the solver node:
!#
!# <verb>
!#      TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
!#
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner)
!#      p_rsolverNode%domega = 0.001   ! or whatever you like
!# </verb>
!#
!# 3.) And where is the good old Jacobi iteration hidden? There is only that
!#     preconditioner!?!
!#
!#     Yes, i have only the Jacobi preconditioner. To get a Jacobi iteration,
!#     you must use it as a preconditioner in a defect correction loop:
!#
!#    <tex> $$ x_{n+1}  =  x_n  +  \omega D^{-1}  (b-Ax) $$ </tex>
!#    <!--  ------------------  ^^^^^^^^^^^^^  ------  Defect correction loop
!#                              Jacobi preconditioner                         -->
!#
!#     So,
!#
!# <verb>
!#      TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
!#
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initJacobi (p_rpreconditioner,0.7_DP)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner)
!# </verb>
!#
!#    Alternatively, you can use the domega parameter from the defect
!#    correction to specify the <tex>$\omega$</tex>, while using no damping parameter
!#    in the actual Jacobi preconditioner:
!#
!# <verb>
!#      NULLIFY(p_rpreconditioner)
!#      CALL linsol_initJacobi (p_rpreconditioner)
!#      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner)
!#      p_rsolverNode%domega = 0.7_DP
!# </verb>
!#
!#  4.) I want to plug in my own solver XY. What do I have to do?
!#
!#      Well, all solvers are build up in the same way with the same
!#      interfaces - that is the reason why every solver can be used as
!#      a preconditioner in another one. If you want to create your
!#      own solver and plug it into here, use the following guideline
!#      of what to do:
!# <verb>
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
!#         Do not forget to deallocate in linsol_doneXXXX everything that
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
!# </verb>
!#      That is it, your solver should now be able to work ;-)
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
!#      lsysbl_assignDiscrIndirect.
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
!#       3.) linsol_initSOR
!#           -> SOR preconditioner / GS preconditioner (=SOR(1.0))
!#           -> see [http://www.netlib.org/linalg/html_templates/Templates.html],
!#                  [http://mathworld.wolfram.com/SuccessiveOverrelaxationMethod.html]
!#
!#       4.) linsol_initSSOR
!#           -> SSOR preconditioner
!#           -> see [http://www.netlib.org/linalg/html_templates/Templates.html],
!#                  [http://mathworld.wolfram.com/SuccessiveOverrelaxationMethod.html]
!#
!#       5.) linsol_initVANKA
!#           -> VANKA preconditioner; multiple versions
!#
!#       6.) linsol_initUMFPACK4
!#           -> UMFPACK preconditioner
!#
!#       7.) linsol_initILU0
!#           -> ILU(0)-preconditioner
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
!#           -> (flexible) Generalised Minimal-Residual preconditioner
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
!#              Simplified interface.
!#
!#      14.) linsol_initSPSOR
!#           -> SP-SOR preconditioner
!#
!#      15.) linsol_initBlockJac
!#           -> Jacobi preconditioner using the block structure of dG discretisations
!#
!#      16.) linsol_initBlockSOR
!#           -> SOR preconditioner using the block structure of dG discretisations
!#
!#      17.) linsol_initDeflGMRES
!#           -> Deflated GMRES algorithm
!#           -> see [Erlangga, Y.; Nabben, R.; Algrbraic Multilevel Krylov Methods;
!#                   SIAM J. Sci. Comput. 2009, Vol 31, Nr. 5, pp. 3417-3437]
!#
!#      18.) linsol_initBlockUpwGS
!#           -> Block GS using upwind numbering of elements (for DG discretisations)
!#
!# 9.) What is this linsol_matricesCompatible?
!#
!#     This allows to check a set of system matrices against a solver node.
!#     If the matrices are not compatible (e.g. because their properties like
!#     matrix format, being transposed,... do not fit) an error code is
!#     returned -- as this set of matrices would lead the solver to an
!#     immediate program stop if they are factorised with linsol_initStructure!
!#     It may be used for special solvers like VANKA to check, if the solver
!#     can handle what is set up in the discretisation.
!#
!# 10.) How is the damping parameter rsolverNode%domega used?
!#
!#     Well, this depends a little bit on the algorithm. Generally can be said
!#     that most preconditioners scale the preconditioned defect by this
!#     value. For example when UMFPACK is used for preconditioning,
!#     UMFPACK returns:
!#      <tex> $$         rd := domega A^(-1) rd  $$ </tex>
!#
!#     thus one should usually use domega=1. 'Direct' solvers usually
!#     act this ways. Exceptions for this rule are algorithms that use this
!#     parameter for 'internal damping' like SOR:
!#
!#      <tex> $$         rd := (D + domega L)^{-1} rd  $$ </tex>
!#
!#     i.e. here not the defect is damped but the operator itself.
!#
!# 11.) What is the difference between "Multigrid 1" and "Multigrid 2" and
!#      which one should I use?
!#
!#      You should use "Multigrid 2". This solver has a simplified
!#      and more clear interface and has some automatism for creating
!#      appropriate multilevel projection routines. It is simply easier.
!#      The original "Multigrid 1" variant used a linear list for maintaining
!#      the levels in hope to be able to add levels more easily. We skipped
!#      this idea in a later implementation as the array implementation
!#      can handle that case with a similar ease.
!#
!# 12.) I want to use a special prolongation/restriction in the multigrid,
!#      but I cannot see how to specify it.
!#
!#      It is possible to specify the prolongation/restriction by
!#      manually initialising the projection structure. There are two
!#      possibilities for that:
!#
!# <verb>
!#      a) Complete user defined projection.
!#         Here the user has do define all projection structures as follows:
!#
!#        1) linsol_initMultigrid2
!#        2) On all levels except for the coarse level
!#             linsol_getMultigrid2Level
!#             + Create a projection structure rprojection for that level
!#             + use linsol_initProjMultigrid2Level(...,rprojection)
!#               to attach the projection structure to that level
!#        3) On all levels:
!#             linsol_getMultigrid2Level
!#             + Add pointers to the coarse grid solve / all smoothers on all levels.
!#
!#        4) linsol_setMatrices
!#        5) linsol_initStructure
!#        ...
!#
!#        So you can attach the user defined projection structure with
!#        linsol_initProjMultigrid2Level after linsol_initMultigrid2.
!#        The application has to maintain the projection structure and
!#        release its memory only after releasing the solver.
!#
!#     b) Changing the default projection.
!#        Here we first initialise the projection and change it later.
!#        1) linsol_initMultigrid2
!#
!#        2) Initialise the default projection
!#             linsol_initProjMultigrid2
!#
!#        3) On all levels except for the coarse level
!#             linsol_getMultigrid2Level(...,p_rlevelInfo)
!#             + Change the p_rlevelInfo%p_rprojection structure as needed
!#
!#        3) On all levels:
!#             linsol_getMultigrid2Level
!#             + Add pointers to the coarse grid solve / all smoothers on all levels.
!#
!#        4) linsol_setMatrices
!#        5) linsol_initStructure
!#        ...
!#
!#        Here, the projection structure is maintained by the solver and
!#        automatically released.
!# </verb>
!# </purpose>
!##############################################################################

module linearsolver

!$use omp_lib
  use fsystem
  use storage
  use spatialdiscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  use filtersupport
  use coarsegridcorrection
  use vanka
  use globalsystem
  use genoutput
  use iluk
  use spsor
  use mprimitives
  use triangulation
  use statistics
  use matrixio
  
  implicit none
  
  private
  
  public :: linsol_setMatrices
  public :: linsol_setMatrix
  public :: linsol_matricesCompatible
  public :: linsol_setOnelevelMatrixDirect
  public :: linsol_initStructure
  public :: linsol_initData
  public :: linsol_updateStructure
  public :: linsol_updateData
  public :: linsol_doneData
  public :: linsol_doneStructure
  public :: linsol_releaseSolver
  public :: linsol_alterSolver
  public :: linsol_testConvergence
  public :: linsol_testDivergence
  public :: linsol_precondDefect
  public :: linsol_solveAdaptively
  public :: linsol_convertToSmoother
  public :: linsol_initSolverGeneral
  
  public :: linsol_initDefCorr
  public :: linsol_initJacobi
  public :: linsol_initSOR
  public :: linsol_initSSOR
  public :: linsol_initVANKA
  public :: linsol_initUMFPACK4
  public :: linsol_initILU0
  public :: linsol_initMILUs1x1
  public :: linsol_initBiCGStab
  public :: linsol_initCG
  public :: linsol_initGMRES
  public :: linsol_initMultigrid
  public :: linsol_initMultigrid2
  public :: linsol_initSchur
  public :: linsol_initSPSOR
  public :: linsol_initBlockJac
  public :: linsol_initBlockSOR
  public :: linsol_initBlockUpwGS
  public :: linsol_initSpecDefl
  public :: linsol_initDeflGMRES
    
  public :: linsol_addMultigridLevel,linsol_addMultigridLevel2
  
  public :: linsol_getMultigrid2Level,linsol_getMultigrid2LevelCount
  
  public :: linsol_getDeflGMRESLevel,linsol_getDeflGMRESLevelCount
  
!  interface linsol_addMultigridLevel
!    module procedure linsol_addMultigridLevel1
!    module procedure linsol_addMultigridLevel2
!  end interface

  interface linsol_initProjMultigrid2Level
    module procedure linsol_initProjMG2LvByPrj
    module procedure linsol_initProjMG2LvByDiscr
  end interface

  interface linsol_initProjDeflGMRESLevel
    module procedure linsol_initPrjDGMLvByPrj
    module procedure linsol_initPrjDGMLvByDiscr
  end interface
  
  public :: linsol_initProjMultigrid2Level
  public :: linsol_initProjMultigrid2

  public :: linsol_initProjDeflGMRESLevel
  public :: linsol_initProjDeflGMRES

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<constants>

!<constantblock description="Algorithm identifiers">

  ! Undefined algorithm
  integer, parameter, public :: LINSOL_ALG_UNDEFINED     = 0
  
  ! Preconditioned defect correction (Richardson iteration);
  ! <tex>$x_{n+1} = x_n + \omega P^{-1} (b-Ax)$</tex>
  integer, parameter, public :: LINSOL_ALG_DEFCORR       = 1
  
  ! Jacobi iteration <tex>$x_1 = x_0 + \omega D^{-1} (b-Ax_0)$</tex>
  integer, parameter, public :: LINSOL_ALG_JACOBI        = 2
  
  ! SOR/GS iteration <tex>$x_1 = x_0 + \omega (D+\omega L)^{-1}(b-Ax_0)$</tex>
  integer, parameter, public :: LINSOL_ALG_SOR           = 4
  
  ! SSOR iteration
  integer, parameter, public :: LINSOL_ALG_SSOR          = 5
  
  ! CG iteration (preconditioned)
  integer, parameter, public :: LINSOL_ALG_CG            = 6
  
  ! BiCGStab iteration (preconditioned)
  integer, parameter, public :: LINSOL_ALG_BICGSTAB      = 7

  ! GMRES iteration (preconditioned)
  integer, parameter, public :: LINSOL_ALG_GMRES         = 8
  
  ! Multigrid iteration
  integer, parameter, public :: LINSOL_ALG_MULTIGRID     = 9

  ! UMFPACK2
  integer, parameter, public :: LINSOL_ALG_UMFPACK2      = 10

  ! UMFPACK4
  integer, parameter, public :: LINSOL_ALG_UMFPACK4      = 11
  
  ! Multigrid iteration
  integer, parameter, public :: LINSOL_ALG_MULTIGRID2    = 12

  ! ILU(0) iteration
  integer, parameter, public :: LINSOL_ALG_ILU0          = 50
  
  ! (M)ILU(s) iteration (scalar system)
  integer, parameter, public :: LINSOL_ALG_MILUS1x1      = 51
  
  ! SPAI(0) iteration (scalar system)
  integer, parameter, public :: LINSOL_ALG_SPAI01x1      = 52

  ! SPAI(k) iteration (scalar system)
  integer, parameter, public :: LINSOL_ALG_SPAIK1x1      = 53
  
  ! VANKA iteration
  integer, parameter, public :: LINSOL_ALG_VANKA         = 54
  
  ! Schur-complement iteration
  integer, parameter, public :: LINSOL_ALG_SCHUR         = 55
  
  ! SP-SOR iteration
  integer, parameter, public :: LINSOL_ALG_SPSOR         = 56
  
  ! Block Jacobi iteration
  integer, parameter, public :: LINSOL_ALG_BLOCKJAC      = 57
  
  ! Spectral deflation iteration
  integer, parameter, public :: LINSOL_ALG_SPECDEFL      = 58

  ! Deflaged MG GMRES
  integer, parameter, public :: LINSOL_ALG_DEFLGMRES     = 59

  ! Block SOR iteration
  integer, parameter, public :: LINSOL_ALG_BLOCKSOR      = 60
  
  ! Block Upwind GS iteration
  integer, parameter, public :: LINSOL_ALG_BLOCKUPWGS    = 61
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Flags for matrix compatibility check.">

  ! Matrices are compatible to a solver node.
  integer, parameter, public :: LINSOL_COMP_OK               = 0
  
  ! Matrices generally not compatible, no special reason.
  integer, parameter, public :: LINSOL_COMP_ERRGENERAL       = 1
  
  ! One of the submatrices is solved in transposed structure,
  ! but one of the subsolvers cannot handle that.
  integer, parameter, public :: LINSOL_COMP_ERRTRANSPOSED    = 2

  ! Some of the submatrices are not saved transposed although
  ! they have to be. (Usual error for special VANKA solvers.)
  integer, parameter, public :: LINSOL_COMP_ERRNOTTRANSPOSED = 3

  ! One of the submatrices is a block matrix,
  ! but one of the subsolvers can only handle scalar matrices.
  integer, parameter, public :: LINSOL_COMP_ERRNOTSCALAR     = 4
  
  ! One of the diagonal submatrices (pivot) is zero, but
  ! the solver cannot handle saddle-point matrices.
  integer, parameter, public :: LINSOL_COMP_ERRZEROPIVOT     = 5
  
  ! One of the submatrices has unsupported matrix type.
  integer, parameter, public :: LINSOL_COMP_ERRMATTYPE       = 6
  
  ! The block matrix is rectangular.
  integer, parameter, public :: LINSOL_COMP_ERRMATRECT       = 7
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Identifiers for stopping criterium istoppingCriterion.">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. the iteration stops when both,
  !    the relative AND the absolute stopping criterium holds
  integer, parameter, public :: LINSOL_STOP_STANDARD     = 0

  ! Use 'minimum' stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use one of them, i.e. the iteration stops when the
  !    either the relative OR the absolute stopping criterium holds
  integer, parameter, public :: LINSOL_STOP_ONEOF        = 1
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Bitfield identifiers for the ability of a linear solver">

  ! Solver can handle scalar systems
  integer(I32), parameter, public :: LINSOL_ABIL_SCALAR       = 2**0

  ! Solver can handle block systems
  integer(I32), parameter, public :: LINSOL_ABIL_BLOCK        = 2**1

  ! Solver can handle multiple levels
  integer(I32), parameter, public :: LINSOL_ABIL_MULTILEVEL   = 2**2
  
  ! Solver allows checking the defect during the iteration.
  ! Solvers not capable of this perform only a fixed number of solution
  ! steps (e.g. UMFPACK performs always one step).
  integer(I32), parameter, public :: LINSOL_ABIL_CHECKDEF     = 2**3
  
  ! Solver is a direct solver (e.g. UMFPACK, ILU).
  ! Otherwise the solver is of iterative nature and might perform
  ! multiple steps to solve the problem.
  integer(I32), parameter, public :: LINSOL_ABIL_DIRECT       = 2**4
  
  ! Solver might use subsolvers (preconditioners, smoothers,...)
  integer(I32), parameter, public :: LINSOL_ABIL_USESUBSOLVER = 2**5
  
  ! Solver supports filtering
  integer(I32), parameter, public :: LINSOL_ABIL_USEFILTER    = 2**6
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Error constants returned by initialisation routines">

  ! Initialisation routine went fine
  integer, parameter, public :: LINSOL_ERR_NOERROR       = 0

  ! Warning: Singular matrix
  integer, parameter, public :: LINSOL_ERR_SINGULAR      = -1

  ! Error during the initialisation: Not enough memory
  integer, parameter, public :: LINSOL_ERR_NOMEMORY      = 1

  ! General error during the initialisation
  integer, parameter, public :: LINSOL_ERR_INITERROR     = 2

  ! Matrix-structure has changed between initStructure and initData
  integer, parameter, public :: LINSOL_ERR_MATRIXHASCHANGED = 3
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Variants of the VANKA solver">

  ! General VANKA solver
  integer, parameter, public :: LINSOL_VANKA_GENERAL           = 0
  
  ! General VANKA solver. Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_GENERALDIRECT     = 1

  ! Simple VANKA, 2D Navier-Stokes problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_2DNAVST           = 2

  ! Simple VANKA, 2D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_2DNAVSTDIRECT     = 3

  ! Full VANKA, 2D Navier-Stokes problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_2DFNAVST          = 4

  ! Full VANKA, 2D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTDIRECT    = 5

  ! Simple VANKA, 2D Navier-Stokes problem, general discretisation,
  ! Solution-based variant.
  integer, parameter, public :: LINSOL_VANKA_2DNAVSTSB         = 6

  ! Simple VANKA, 2D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  ! Solution-based variant.
  integer, parameter, public :: LINSOL_VANKA_2DNAVSTDIRECTSB   = 7

  ! Full VANKA, 2D Navier-Stokes optimal control problem, general discretisation.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOC        = 20

  ! Full VANKA, 2D Navier-Stokes optimal control problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOCDIRECT  = 21

  ! Diagonal VANKA, 2D Navier-Stokes optimal control problem, general discretisation.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOCDIAG    = 22

  ! Diagonal VANKA, 2D Navier-Stokes optimal control problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOCDIAGDIR = 23

  ! Diagonal VANKA, 2D Navier-Stokes optimal control problem, general discretisation,
  ! new implementation.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOCDIAG2   = 24

  ! Diagonal VANKA, 2D Navier-Stokes optimal control problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOCDIAGDIR2 = 25

  ! Full VANKA, 2D Navier-Stokes optimal control problem, general discretisation,
  ! new implementation.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOCFULL2   = 26

  ! Full VANKA, 2D Navier-Stokes optimal control problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_2DFNAVSTOCFULLDIR2 = 27

  ! Simple VANKA, 3D Navier-Stokes problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_3DNAVST           = 30

  ! Simple VANKA, 3D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_3DNAVSTDIRECT     = 31

  ! Full VANKA, 3D Navier-Stokes problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_3DFNAVST          = 32

  ! Full VANKA, 3D Navier-Stokes problem, general discretisation.
  ! Specialised 'direct' version, i.e. when
  ! used as a smoother in multigrid, this bypasses the usual defect
  ! correction approach to give an additional speedup.
  integer, parameter, public :: LINSOL_VANKA_3DFNAVSTDIRECT    = 33

  ! -------------- NEW IMPLEMENTATION --------------

  ! Simple VANKA, 2D Navier-Stokes problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_NAVST2D_DIAG = 101

  ! Full VANKA, 2D Navier-Stokes problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_NAVST2D_FULL = 102

  ! Simple VANKA, 2D Boussinesq problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_BOUSS2D_DIAG = 121

  ! Full VANKA, 2D Boussinesq problem, general discretisation
  integer, parameter, public :: LINSOL_VANKA_BOUSS2D_FULL = 122

!</constantblock>

! *****************************************************************************

!<constantblock description="Variants of the SP-SOR solver">

  ! SP-SOR solver for 2D Navier-Stokes systems
  integer, parameter, public :: LINSOL_SPSOR_NAVST2D        = 1

  ! 'diagonal' SP-SOR solver for 2D Navier-Stokes systems
  integer, parameter, public :: LINSOL_SPSOR_NAVST2D_DIAG   = 2

!</constantblock>

! *****************************************************************************

!<constantblock description="default values for the GMRES(m) solver">

  ! One or two Gram-Schmidt calls per GMRES iteration
  logical, parameter, public :: LINSOL_GMRES_DEF_TWICE_GS = .false.
    
!</constantblock>

! *****************************************************************************

!<constantblock description="types of the Schur-complement solver">

  ! Diagonal type Schur-complement solver
  integer, parameter, public :: LINSOL_SCHUR_TYPE_DIAG = 0
  
  ! Lower-triangular type Schur-complement solver
  integer, parameter, public :: LINSOL_SCHUR_TYPE_LTRI = 1
  
  ! Upper-triangular type Schur-complement solver (not implemented)
  !integer, parameter, public :: LINSOL_SCHUR_TYPE_UTRI = 2
  
  ! Full type Schur-complement solver (not implemented)
  !integer, parameter, public :: LINSOL_SCHUR_TYPE_FULL = 3
    
!</constantblock>

! *****************************************************************************

!<constantblock description="Possible commands for linsol_alterSolver">

  ! Dummy command, do-nothing.
  integer, parameter, public :: LINSOL_ALTER_NOTHING                   = 0
  
  ! Change VANKA subtype.
  ! In the configuration block, there must be specified:
  ! Iconfig(1) = identifier for VANKA subtype to be changed.
  ! Iconfig(2) = destination subtype of VANKA solver
  integer, parameter, public :: LINSOL_ALTER_CHANGEVANKA               = 1
    
!</constantblock>

! *****************************************************************************

!<constantblock description="Private constants">

  ! Length of the residual queue in iterative solvers.
  ! This defines how many residuals are saved and can be used for
  ! output, e.g. for the asymptotic convergence rate and for printing
  ! the residual queue in case of errors.
  integer, parameter :: LINSOL_NRESQUEUELENGTH                 = 32
  
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
  ! by a solver depends on the solver (e.g. UMFPACK will not respect
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
  
  type t_linsolNode
    
    ! OUTPUT: Result
    ! The result of the solution process.
    ! =0: success. =1: iteration broke down, diverging, =2: error in the parameters,
    ! <0: algorithm-specific error
    integer                    :: iresult = 0
    
    ! OUTPUT: Number of performed iterations, if the solver
    ! is of iterative nature.
    ! Is to 1 by the solver if not used (indicating at least 1 performed
    ! iteration, which is always the case).
    integer                    :: iiterations = 0
    
    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK:
    ! Norm of initial residuum
    real(DP)                        :: dinitialDefect = 0.0_DP

    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK:
    ! Norm of final residuum
    real(DP)                        :: dfinalDefect = 0.0_DP

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK:
    ! Convergence rate
    real(DP)                        :: dconvergenceRate = 0.0_DP

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK:
    ! Asymptotic convergence rate
    real(DP)                        :: dasymptoticConvergenceRate = 0.0_DP

    ! OUTPUT PARAMETER:
    ! Total time for solver
    real(DP)                        :: dtimeTotal = 0.0_DP

    ! OUTPUT PARAMETER FOR SOLVERS THAT SUPPORT FILTERING:
    ! Total time for filtering (not yet implemented)
    !real(DP)                        :: dtimeFiltering = 0.0_DP
    
    ! INPUT PARAMETER:
    ! Damping parameter for preconditioner. The t_linsolNode structure
    ! represents a preconditioner operator of the form:
    !
    !              <tex>$d --> omega * P^{-1} * d$</tex>
    !
    ! The actual result of the solver algorithm is scaled by domega
    ! after the solving process.
    real(DP)                        :: domega  = 1.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP)                        :: depsRel = 1E-5_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSABS.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP)                        :: depsAbs = 1E-5_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Difference in the residual stopping criterion. Stop iteration if
    ! !!defect_new!! - !!defect_old!! < depsDiff * !!defect_old!!.
    ! This stopping criterion holds additionally to depsAbs/depsRel
    ! in an OR sense -- i.e. the iteration is stopped if
    ! !!defect_new!! - !!defect_old!! < depsDiff * !!defect_old!!
    ! holds OR if the stopping criterion given by depsAbs/depsRel
    ! holds!
    ! =0: ignore this stopping criterion
    real(DP)                        :: depsDiff = 0.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Relative divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY_DP disables the relative divergence check.
    ! standard = 1E3
    real(DP)                        :: ddivRel = 1E6_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Absolute divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY_DP disables the absolute divergence check.
    ! standard = SYS_INFINITY_DP
    real(DP)                        :: ddivAbs = SYS_INFINITY_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! RHS-vector is treated as zero if max(defect) < drhsZero
    real(DP)                        :: drhsZero = 1E-90_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Type of stopping criterion to use. One of the
    ! LINSOL_STOP_xxxx constants.
    integer                    :: istoppingCriterion = LINSOL_STOP_STANDARD

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Minimum number of iterations top perform
    integer                    :: nminIterations = 1

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Maximum number of iterations top perform
    integer                    :: nmaxIterations = 50
    
    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK:
    ! Perform residual checks
    ! YES: check residuals (default),
    ! NO: do not check residuals, simply perform as many iterations as
    ! configured by nminIterations. Is typically set to NO for smoothing
    ! with a solver.
    integer                    :: iresCheck = YES

    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK:
    ! Type of norm to use in the residual checking (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    integer                    :: iresNorm = 2
    
    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Number of iterations that should be used to calculate
    ! the asymptotic convergtence rate, if the solver supports
    ! to generate an asymptotic convergence rate of the last
    ! couple of steps. =0: do not determine asymptotic convergence rate
    integer                    :: niteAsymptoticCVR = 0
    
    ! INPUT PARAMETER: Output level
    ! This determines the output level of the solver.
    ! =-1: no output, =0: no output except for warning messages,
    ! =1: basic output, =2, extended output
    integer                    :: ioutputLevel = 0

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK:
    ! Number of iterations to perform before printing out the
    ! norm of the residual to screen.
    ! =1: Print residual in every iteration
    integer                    :: niteResOutput = 1
    
    ! INPUT PARAMETER: Default data type.
    ! This parameter prescribes the default data type that is used
    ! for temporary verctors inside the solvers. Solvers that need
    ! temporary vectors allocate them using this data type.
    ! The type identifier (either ST_SINGLE or ST_DOUBLE (default))
    ! should match the data type of the RHS/solution vectors,
    ! otherwise there might be performance loss!
    integer                    :: cdefaultDataType = ST_DOUBLE

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
    integer                     :: isolverSubgroup = 0
    
    ! INPUT PARAMETER: t_matrixBlock structure that holds
    ! information about the linear system where to solve
    ! (i.e. usually the finest level).
    ! For multilevel-capable algorithms, information about the other
    ! levels are saved in local solver-specific structures.
    type(t_matrixBlock)                :: rsystemMatrix

    ! READ ONLY: Algorithm identifier.
    ! One of the LINSOL_ALG_xxxx constants.
    ! Depending on the value, a solver-specific structure may be
    ! assigned to this structure below.
    integer                    :: calgorithm = LINSOL_ALG_UNDEFINED
    
    ! READ ONLY: Solver ability tag.
    ! Bitfield. A combination of LINSOL_ABIL_xxxx
    ! flags that specify the ability of the actual solver (e.g.
    ! whether it can handle block-matrices or only scalar matrices,...).
    ! Calling a solver that is not able to handle a specific problem
    ! will cause an error by the framework.
    integer(I32)                    :: ccapability = 0
    
    ! STATUS FOR ITERATIVE SOLVERS: Current iteration
    integer                    :: icurrentIteration

    ! STATUS FOR ITERATIVE SOLVERS: Last defect
    real(DP) :: dlastDefect = 0.0_DP

    ! Pointer to a structure for the VANKA solver; NULL() if not set
    type (t_linsolSubnodeVANKA), pointer          :: p_rsubnodeVANKA       => null()

    ! Pointer to a structure for the SP-SOR solver; NULL() if not set
    type (t_linsolSubnodeSPSOR), pointer          :: p_rsubnodeSPSOR       => null()
    
    ! Pointer to a structure for the Defect correction solver; NULL() if not set
    type (t_linsolSubnodeDefCorr), pointer        :: p_rsubnodeDefCorr     => null()

    ! Pointer to a structure for the BiCGStab solver; NULL() if not set
    type (t_linsolSubnodeBiCGStab), pointer       :: p_rsubnodeBiCGStab    => null()

    ! Pointer to a structure for the UMFPACK4 solver; NULL() if not set
    type (t_linsolSubnodeUMFPACK4), pointer       :: p_rsubnodeUMFPACK4    => null()

    ! Pointer to a structure for the Multigrid solver; NULL() if not set
    type (t_linsolSubnodeMultigrid), pointer      :: p_rsubnodeMultigrid   => null()

    ! Pointer to a structure for the Multigrid solver; NULL() if not set
    type (t_linsolSubnodeMultigrid2), pointer     :: p_rsubnodeMultigrid2  => null()

    ! Pointer to a structure for the ILU0 solver; NULL() if not set
    type (t_linsolSubnodeILU0), pointer           :: p_rsubnodeILU0        => null()
    
    ! Pointer to a structure for (M)ILUs 1x1 solver; NULL() if not set
    type (t_linsolSubnodeMILUs1x1), pointer       :: p_rsubnodeMILUs1x1    => null()

    ! Pointer to a structure for SOR; NULL() if not set
    type (t_linsolSubnodeSOR), pointer            :: p_rsubnodeSOR         => null()

    ! Pointer to a structure for SSOR; NULL() if not set
    type (t_linsolSubnodeSSOR), pointer           :: p_rsubnodeSSOR        => null()
    
    ! Pointer to a structure for the CG solver; NULL() if not set
    type (t_linsolSubnodeCG), pointer             :: p_rsubnodeCG          => null()
    
    ! Pointer to a structure for the GMRES(m) solver; NULL() if not set
    type (t_linsolSubnodeGMRES), pointer          :: p_rsubnodeGMRES       => null()
    
    ! Pointer to a structure for the Schur-complement solver; NULL() if not set
    type (t_linsolSubnodeSchur), pointer          :: p_rsubnodeSchur       => null()

    ! Pointer to a structure for the spectral deflation preconditioner; NULL() if not set
    type (t_linsolSubnodeSpecDefl), pointer       :: p_rsubnodeSpecDefl    => null()

    ! Pointer to a structure for the deflation GMRES preconditioner; NULL() if not set
    type (t_linsolSubnodeDeflGMRES), pointer      :: p_rsubnodeDeflGMRES    => null()
    
    ! Pointer to a structure for the block upwind GS preconditioner; NULL() if not set
    type (t_linsolSubnodeBlockUpwGS), pointer     :: p_rsubnodeBlockUpwGS   => null()

  end type
  
  public :: t_linsolNode

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the SOR solver.
  type t_linsolSubnodeSOR
  
    ! Relaxation parameter for SOR. Must be in range (0,2).
    real(DP) :: drelax
  
  end type
  
  public :: t_linsolSubnodeSOR

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the SSOR solver.
  
  type t_linsolSubnodeSSOR
  
    ! Relaxation parameter for SSOR. Must be in range (0,2).
    real(DP) :: drelax
  
    ! Scaling of the preconditioned solution.
    ! =FALSE is the old FEAT style.
    ! =TRUE gives the implementation as suggested in the literature.
    logical :: bscale
    
  end type
  
  public :: t_linsolSubnodeSSOR

!</typeblock>
  
! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the defect correction solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the
  ! preconditioning.
  
  type t_linsolSubnodeDefCorr
  
    ! Temporary vector to use during the solution process
    type(t_vectorBlock) :: rtempVector

    ! 2nd Temporary vector to use during the solution process
    type(t_vectorBlock) :: rtempVector2

    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    type(t_linsolNode), pointer       :: p_rpreconditioner            => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:), pointer      :: p_RfilterChain => null()
  
  end type
  
  public :: t_linsolSubnodeDefCorr

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the VANKA solver.
  
  type t_linsolSubnodeVANKA
  
    ! Type of VANKA subsolver. One of the LINSOL_VANKA_xxxx constants.
    ! Specifies either thge general VANKA solver or a specialised version
    ! for improved speed.
    integer             :: csubtypeVANKA
  
    ! For general 2D/3D Navier-Stokes problem
    type(t_vanka)       :: rvanka
  
    ! Temporary vector to use during the solution process
    type(t_vectorBlock) :: rtempVector

  end type
  
  public :: t_linsolSubnodeVANKA
  
!</typeblock>


! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the SP-SOR solver.
  type t_linsolSubnodeSPSOR
  
    ! Type of SP-SOR subsolver. One of the LINSOL_SPSOR_xxxx constants.
    integer             :: csubtype
  
    ! The SP-SOR data structure
    type(t_spsor)       :: rdata

  end type
  
  public :: t_linsolSubnodeSPSOR
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the BiCGStab solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the
  ! preconditioning.
  
  type t_linsolSubnodeBiCGStab
  
    ! Temporary vectors to use during the solution process
    type(t_vectorBlock), dimension(6) :: RtempVectors
    
    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    type(t_linsolNode), pointer       :: p_rpreconditioner            => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:), pointer      :: p_RfilterChain => null()
  
  end type
  
  public :: t_linsolSubnodeBiCGStab
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the CG solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the
  ! preconditioning.
  
  type t_linsolSubnodeCG
  
    ! Temporary vectors to use during the solution process
    type(t_vectorBlock), dimension(4) :: RtempVectors

    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    type(t_linsolNode), pointer       :: p_rpreconditioner            => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:), pointer      :: p_RfilterChain => null()
  
  end type

  public :: t_linsolSubnodeCG
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the UMFPACK4 solver.
  
  type t_linsolSubnodeUMFPACK4
  
    ! Matrix output for debug.
    ! If this is set to a value <> 0, the numerical factorisation routine
    ! writes the matrix to a text file before factorising it.
    ! The text file gets the name 'matrixN.txt' with N=imatrixDebugOutput.
    ! This is for debugging purposes and should be used with care,
    ! as the text files grow rather quickly with the dimension!
    !
    ! =0: no output.
    ! =1: Write the matrix into a text file in human readable form
    !     using matio_writeMatrixHR.
    ! =2: Write the matrix into a text file in Maple format
    !     using matio_writeMatrixMaple.
    integer :: imatrixDebugOutput = 0
    
    ! Name of the matrix in the matrix file when using imatrixDebugOutput<>0.
    ! If "" is used here, a matrix name is automatically chosen.
    ! This string is also used as filename for the output file.
    character(LEN=SYS_NAMELEN) :: smatrixName = ""
  
    ! Accuracy string for text output of matrix entries
    ! if imatrixDebugOutput <> 0.
    character(LEN=SYS_NAMELEN) :: saccuracy = "(E15.5)"
  
    ! Control structure for UMFPACK4; contains parameter for the solver
    real(DP), dimension(20) :: Dcontrol

    ! Handle for symbolic factorisation.
    ! This is not a FEAT-Handle!
    integer :: isymbolic = 0
!    integer(I32) :: isymbolic = 0

    ! Handle for numeric factorisation
    ! This is not a FEAT-Handle!
    integer :: inumeric = 0
!    integer(I32) :: inumeric = 0
    
    ! Handle to a temporary vector for storing the solution
    type(t_vectorBlock) :: rtempVector

  end type

  public :: t_linsolSubnodeUMFPACK4
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the flexible GMRES(m) solver.
  
  type t_linsolSubnodeGMRES
  
    ! Dimension of Krylov Subspace
    integer :: ikrylovDim
    
    ! Maybe we should apply Gram-Schmidt twice?
    logical :: btwiceGS
    
    ! Scale factor for pseudo-residuals
    real(DP) :: dpseudoResScale
    
    ! Some temporary 1D/2D arrays
    real(DP), dimension(:),   pointer :: Dc, Ds, Dq
    real(DP), dimension(:,:), pointer :: Dh

    ! The handles of the arrays
    integer :: hDc, hDs, hDq, hDh
    
    ! Some temporary vectors
    type(t_vectorBlock), dimension(:), pointer :: p_rv                => null()
    type(t_vectorBlock), dimension(:), pointer :: p_rz                => null()
    type(t_vectorBlock) :: rx
  
    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    type(t_linsolNode), pointer       :: p_rpreconditioner            => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:), pointer      :: p_RfilterChain => null()

  end type

  public :: t_linsolSubnodeGMRES
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the scalar ILU(0) solver.
  
  type t_linsolSubnodeILU0
  
    ! A pointer to the ILU0-decomposition of the main matrix.
    type(t_matrixBlock) :: rmatLU
    
  end type
  
  public :: t_linsolSubnodeILU0

!</typeblock>
  
! *****************************************************************************

!<typeblock>
  
  ! This saves information about ILU(k) matrices as defined in SPLIB.
  
  type t_linsolSubnodeMILUs1x1

    ! Fill-in level for decomposition
    integer :: ifill
    
    ! Relaxation factor for (M)ILU(s)
    real(DP) :: drelax
    
    ! Scaling factor of the matrix; usually = 1.0
    real(DP) :: dscaleFactor

    ! A structure to store the (M)ILU-decomposition
    ! the decomposition is stored in a modified sparse
    ! row format (MSR)
    type(t_MILUdecomp) :: rMILUdecomp

  end type

  public :: t_linsolSubnodeMILUs1x1
  
!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure saves the information for the Schur-Complement solver.
  
  type t_linsolSubnodeSchur
  
    ! Type of preconditioner. One of the LINSOL_SCHUR_TYPE_XXXX constants.
    integer :: ctype = LINSOL_SCHUR_TYPE_DIAG
  
    ! An array of block matrices that saves the sub-matrices A.
    type(t_matrixBlock), dimension(:), pointer :: p_RmatrixA => null()
    
    ! An array of block discretisations for the sub-matrices A.
    type(t_blockDiscretisation), dimension(:), pointer :: p_RdiscrA => null()
  
    ! A block matrix that saves the sub-matrix B.
    type(t_matrixBlock) :: rmatrixB
  
    ! A block matrix that saves the sub-matrix D.
    type(t_matrixBlock) :: rmatrixD
    
    ! An array of block matrices that saves the sub-matrices S.
    type(t_matrixBlock), dimension(:), pointer :: p_RmatrixS => null()
    
    ! A pointer to a solver structure for the sub-matrix A.
    type(t_linsolNode), pointer :: p_rsolverA => null()
    
    ! A pointer to a solver structure for the Schur-complement matrix S.
    type(t_linsolNode), pointer :: p_rsolverS => null()
  
  end type
  
  public :: t_linsolSubnodeSchur

!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure collects timing information for the multigrid solver.

  type t_linsolMGTiming
    
    ! Total time for prolongation
    real(DP)                           :: d_tmProlongation

    ! Total time for restriction
    real(DP)                           :: d_tmRestriction
  
    ! Total time for defect calculation
    real(DP)                           :: d_tmDefCalc
    
    ! Total time for smoothing
    real(DP)                           :: d_tmSmooth
    
    ! Total time for coarse grid solving
    real(DP)                           :: d_tmCoarseGridSolve

    ! Total time for coarse grid corrections
    real(DP)                           :: d_tmCorrection
  
  end type

  public :: t_linsolMGTiming

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_linsolSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  type t_linsolMGLevelInfo
  
    ! A level identifier. Can be set to identify the global number of the
    ! level, this structure refers to. Only for output purposes.
    integer                        :: ilevel                 = 0
    
    ! Pointer to previous (lower) level or NULL() if the node is the
    ! head of the list
    type(t_linsolMGLevelInfo), pointer  :: p_rprevLevel           => null()

    ! Pointer to next (higher) level or NULL() if the node is the
    ! tail of the list
    type(t_linsolMGLevelInfo), pointer  :: p_rnextLevel           => null()
  
    ! t_matrixBlock structure that holds the system matrix.
    type(t_matrixBlock)                :: rsystemMatrix
    
    ! t_matrixBlock structure for the calculation of the optimal correction.
    ! Normally, this matrix has not to be set. If not defined, rsystemMatrix is
    ! used. If defined, the calculation of the optimal defect (guided by
    ! the parameter rcoarseGridCorrection in the multigrid subnode)
    ! is done with this matrix instead of the system matrix.
    type(t_matrixBlock)                :: rsystemMatrixOptCorrection
    
    ! A temporary vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_vectorBlock)                 :: rtempVector
    
    ! A temporary vector that is used for the calculation of the
    ! coarse grid correction. This vector shares its data with the
    ! rprjTempVector vector in the multigrid solver structure
    ! t_linsolSubnodeMultigrid, but has the shape of the RHS vector.
    ! The structure is initialised in initStructure and cleaned up
    ! in doneStructure.
    type(t_vectorBlock)                 :: rcgcorrTempVector

    ! A RHS vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_vectorBlock)                 :: rrhsVector

    ! A solution vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_vectorBlock)                 :: rsolutionVector

    ! A pointer to the solver node for the presmoothing or NULL(),
    ! if no presmoother is used.
    type(t_linsolNode), pointer         :: p_rpresmoother         => null()

    ! A pointer to the solver node for the postsmoothing or NULL(),
    ! if no presmoother is used.
    type(t_linsolNode), pointer         :: p_rpostsmoother        => null()

    ! A pointer to the solver node for the coarse grid solver if
    ! the level corresponding to this structure represents the coarse
    ! grid. NULL() otherwise.
    type(t_linsolNode), pointer         :: p_rcoarseGridSolver    => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:),pointer      :: p_RfilterChain => null()
    
    ! Specifies whether multigrid should set up the rprojection depending on
    ! the system matrix that is passed using the linsol_setMatrices routine or
    ! use a user-specified projection structure.
    ! If set to .true., multigrid will set up the projection structure by itself.
    ! If ser to .false., multigrid will use the projection structure that the
    ! user specified when we called the linsol_addMultiGridLevel routine.
    logical                             :: bautoProjection = .false.
    
    ! An interlevel projection structure that configures the transfer of
    ! defect/solution vectors from one level to another.
    ! For the coarse grid, this structure is ignored.
    ! For a finer grid (e.g. level 4), this defines the grid transfer
    ! between the current and the lower level (so here between level 3 and
    ! level 4).
    type(t_interlevelProjectionBlock)   :: rprojection
    
    ! STATUS/INTERNAL: MG cycle information.
    ! Number of cycles to perform on this level.
    integer                        :: ncycles

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of remaining cycles to perform on this level.
    integer                        :: ncyclesRemaining
    
    ! STATUS/INTERNAL: initial residuum when a solution is restricted to this
    ! level. Only used if adaptive cycles are activated by setting
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    real(DP)                       :: dinitResCycle = 0.0_DP
    
    ! STATUS/INTERNAL: Number of current cycle on that level.
    ! Only used if adaptive cycles are activated by setting
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    integer                        :: icycleCount   = 0
    
  end type

  public :: t_linsolMGLevelInfo
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  type t_linsolSubnodeMultigrid
  
    ! INPUT PARAMETER: Cycle identifier.
    !  0=F-cycle,
    !  1=V-cycle,
    !  2=W-cycle.
    integer                       :: icycle                   = 0
    
    ! INPUT PARAMETER: Adaptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsRelCycle * |initial res on that level|.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    real(DP)                      :: depsRelCycle             = 1E99_DP
    
    ! INPUT PARAMETER: If adaptive cycles are activated by depsRelCycle < 1E99_DP,
    ! this configures the maximum number of cycles that are performed.
    ! The value =-1 deactivates this upper bouns, so coarse grid cycles are
    ! repeated until the convergence criterion given by depsRelCycle is reached.
    integer                       :: nmaxAdaptiveCycles       = -1
    
    ! INPUT PARAMETER: Number of levels in the linked list of multigrid levels
    integer                       :: nlevels                  = 0
    
    ! INPUT PARAMETER: Coarse grid correction structure for step length control.
    ! Defines the algorithm for computing the optimal correction as well as the
    ! minimum and maximum step length ALPHAMIN/ALPHAMAX.
    ! The standard setting/initialisation is suitable for conforming elements.
    type(t_coarseGridCorrection)  :: rcoarseGridCorrection
    
    ! Pointer to the head of the linked list of levels; corresponds
    ! to the lowest level.
    type(t_linsolMGLevelInfo), pointer :: p_rlevelInfoHead         => null()

    ! Pointer to the tail of the linked list of levels; corresponds
    ! to the highest level, i.e. the level where the system should be
    ! solved.
    type(t_linsolMGLevelInfo), pointer :: p_rlevelInfoTail         => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:),pointer :: p_RfilterChain     => null()
  
    ! A temp vector for the prolongation/restriction.
    ! Memory for this vector is allocated in initStructure and released
    ! in doneStructure.
    type(t_vectorScalar) :: rprjTempVector
    
  end type

  public :: t_linsolSubnodeMultigrid
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_linsolSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  type t_linsolMG2LevelInfo
  
    ! A level identifier. Can be set to identify the global number of the
    ! level, this structure refers to. Only for output purposes.
    integer                        :: ilevel                 = 0
    
    !<!-- ----------------------------------------------------------------------
    ! Specification of the subsolver components.
    ! Must be initialised by the user prior to calling multigrid!
    ! --------------------------------------------------------------------- -->
    
    ! A pointer to the solver node for the presmoothing or NULL(),
    ! if no presmoother is used.
    type(t_linsolNode), pointer         :: p_rpresmoother         => null()

    ! A pointer to the solver node for the postsmoothing or NULL(),
    ! if no presmoother is used.
    type(t_linsolNode), pointer         :: p_rpostsmoother        => null()

    ! A pointer to the solver node for the coarse grid solver if
    ! the level corresponding to this structure represents the coarse
    ! grid. NULL() otherwise.
    type(t_linsolNode), pointer         :: p_rcoarseGridSolver    => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:),pointer      :: p_RfilterChain => null()
    
    !<!-- ----------------------------------------------------------------------
    ! OPTIONAL parameters.
    ! May be initialised by the user in case of special situations
    ! --------------------------------------------------------------------- -->
    
    ! OPTIONAL: An interlevel projection structure that configures the transport
    ! of defect/solution vectors from one level to another.
    ! For the coarse grid, this structure is ignored.
    ! For a finer grid (e.g. level 4), this defines the grid transfer
    ! between the current and the lower level (so here between level 3 and
    ! level 4).
    !
    ! If the user needs a special projection strategy, he/she can either
    ! set a pointer to a user defined projection strategy here or
    ! modify the strategy after doing linsol_initStructure().
    type(t_interlevelProjectionBlock), pointer   :: p_rprojection => null()
    
    ! <!-- ----------------------------------------------------------------------
    ! INTERNAL VARIABLES.
    ! These variables are internally maintained by the MG solver and do not have
    ! to be modified in any sense.
    ! ---------------------------------------------------------------------- -->
    
    ! Specifies if we use the standard interlevel projection strategy.
    ! If this is set to FALSE, the user set the p_rprojection to
    ! a user defined projection strategy, so we do not have to touch it.
    logical :: bstandardProjection = .false.
    
    ! t_matrixBlock structure that holds the system matrix.
    type(t_matrixBlock)                :: rsystemMatrix
    
    ! t_matrixBlock structure for the calculation of the optimal correction.
    ! Normally, this matrix has not to be set. If not defined, rsystemMatrix is
    ! used. If defined, the calculation of the optimal defect (guided by
    ! the parameter rcoarseGridCorrection in the multigrid subnode)
    ! is done with this matrix instead of the system matrix.
    type(t_matrixBlock)                :: rsystemMatrixOptCorrection
    
    ! A temporary vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_vectorBlock)                 :: rtempVector
    
    ! A temporary vector that is used for the calculation of the
    ! coarse grid correction. This vector shares its data with the
    ! rprjTempVector vector in the multigrid solver structure
    ! t_linsolSubnodeMultigrid, but has the shape of the RHS vector.
    ! The structure is initialised in initStructure and cleaned up
    ! in doneStructure.
    type(t_vectorBlock)                 :: rcgcorrTempVector

    ! A RHS vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_vectorBlock)                 :: rrhsVector

    ! A solution vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_vectorBlock)                 :: rsolutionVector

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of cycles to perform on this level.
    integer                        :: ncycles

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of remaining cycles to perform on this level.
    integer                        :: ncyclesRemaining
    
    ! STATUS/INTERNAL: initial residuum when a solution is restricted to this
    ! level. Only used if adaptive cycles are activated by setting
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    real(DP)                       :: dinitResCycle = 0.0_DP
    
    ! STATUS/INTERNAL: Number of current cycle on that level.
    ! Only used if adaptive cycles are activated by setting
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    integer                        :: icycleCount   = 0
    
  end type

  public :: t_linsolMG2LevelInfo
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  type t_linsolSubnodeMultigrid2
  
    ! INPUT PARAMETER: Cycle identifier.
    !  0=F-cycle,
    !  1=V-cycle,
    !  2=W-cycle.
    integer                       :: icycle                   = 0
    
    ! INPUT PARAMETER: Adaptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsRelCycle * |initial res on that level|.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    real(DP)                      :: depsRelCycle             = 1E99_DP
    
    ! INPUT PARAMETER: If adaptive cycles are activated by depsRelCycle < 1E99_DP,
    ! this configures the maximum number of cycles that are performed.
    ! The value =-1 deactivates this upper bouns, so coarse grid cycles are
    ! repeated until the convergence criterion given by depsRelCycle is reached.
    integer                       :: nmaxAdaptiveCycles       = -1
    
    ! INPUT PARAMETER: Coarse grid correction structure for step length control.
    ! Defines the algorithm for computing the optimal correction as well as the
    ! minimum and maximum step length ALPHAMIN/ALPHAMAX.
    ! The standard setting/initialisation is suitable for conforming elements.
    type(t_coarseGridCorrection)  :: rcoarseGridCorrection
    
    ! Number of levels level in MG.
    integer :: nlevels = 1
    
    ! A pointer to the level info structures for all the levels in multigrid
    type(t_linsolMG2LevelInfo), dimension(:), pointer :: p_RlevelInfo
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:),pointer :: p_RfilterChain     => null()
  
    ! A temp vector for the prolongation/restriction.
    ! Memory for this vector is allocated in initStructure and released
    ! in doneStructure.
    type(t_vectorScalar) :: rprjTempVector
    
  end type
  
  public :: t_linsolSubnodeMultigrid2
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>

  ! This structure realises a configuration block that can be passed to the
  ! function linsol_alterSolver to allow in-computation change of
  ! a solver or its subsolvers.
  type t_linsol_alterSolverConfig
  
    ! A command flag of type LINSOL_ALTER_xxxx.
    integer :: ccommand = LINSOL_ALTER_NOTHING
  
    ! An integer precision configuration block. The meaning of this block
    ! depends on ccommand.
    integer, dimension(16) :: Iconfig
    
    ! A double precision configuration block. The meaning of this block
    ! depends on ccommand.
    real(DP), dimension(16) :: Dconfig
  end type
  
  public :: t_linsol_alterSolverConfig

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the spectral deflation preconditioner.
  type t_linsolSubnodeSpecDefl
  
    ! Relaxation parameter.
    real(DP) :: drelax = 1.0_DP
    
    ! Guess for the maximum eigenvalue
    real(DP) :: dmaxEigenvalue = 0.0_DP
    
    ! Pointer to a preconditioner. => null if not used.
    type(t_linsolNode), pointer :: p_rpreconditioner => null()
    
    ! Temporary vector
    type(t_vectorBlock) :: rvectorTemp
  
    ! Pointer to a filter chain (i.e. an array of t_filterChain
    ! structures) if filtering should be applied to the vector during the
    ! iteration. =>null if no filter is to be applied.
    type(t_filterChain), dimension(:), pointer :: p_Rfilter => null()

  end type
  
  public :: t_linsolSubnodeSpecDefl

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_linsolSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  type t_linsolDeflGMRESLevelInfo
  
    ! A level identifier. Can be set to identify the global number of the
    ! level, this structure refers to. Only for output purposes.
    integer                        :: ilevel                 = 0
    
    !<!-- ----------------------------------------------------------------------
    ! Deflated MG GMRES parameters
    ! --------------------------------------------------------------------- -->
    
    ! Guess for the maximum eigenvalue
    real(DP) :: dmaxEigenvalue = 1.0_DP

    ! Number of solver iterations applied on this level.
    ! Applies to all levels except for the maximum level.
    ! =-1: niterations equals to ikrylowDim, the dimension of the Krylow space.
    integer :: niterations = -1
    
    !<!-- ----------------------------------------------------------------------
    ! OPTIONAL parameters.
    ! May be initialised by the user in case of special situations
    ! --------------------------------------------------------------------- -->
    
    ! OPTIONAL: An interlevel projection structure that configures the transport
    ! of defect/solution vectors from one level to another.
    ! For the coarse grid, this structure is ignored.
    ! For a finer grid (e.g. level 4), this defines the grid transfer
    ! between the current and the lower level (so here between level 3 and
    ! level 4).
    !
    ! If the user needs a special projection strategy, he/she can either
    ! set a pointer to a user defined projection strategy here or
    ! modify the strategy after doing linsol_initStructure().
    type(t_interlevelProjectionBlock), pointer   :: p_rprojection => null()
    
    ! <!-- ----------------------------------------------------------------------
    ! INTERNAL VARIABLES.
    ! These variables are internally maintained by the MG solver and do not have
    ! to be modified in any sense.
    ! ---------------------------------------------------------------------- -->
    
    ! Specifies if we use the standard interlevel projection strategy.
    ! If this is set to FALSE, the user set the p_rprojection to
    ! a user defined projection strategy, so we do not have to touch it.
    logical :: bstandardProjection = .false.
    
    ! t_matrixBlock structure that holds the system matrix on this level.
    type(t_matrixBlock)                :: rsystemMatrix

    !<!-- ----------------------------------------------------------------------
    ! General GMRES parameters
    ! --------------------------------------------------------------------- -->

    ! Dimension of Krylow Subspace
    integer :: ikrylowDim = 5
    
    ! Some temporary 1D/2D arrays
    real(DP), dimension(:),   pointer :: Dc, Ds, Dq
    real(DP), dimension(:,:), pointer :: Dh

    ! The handles of the arrays
    integer :: hDc, hDs, hDq, hDh
    
    ! Some temporary vectors
    type(t_vectorBlock), dimension(:), pointer :: p_rv => null()
    type(t_vectorBlock), dimension(:), pointer :: p_rz => null()
    type(t_vectorBlock) :: rx

    type(t_vectorBlock) :: rtempVector
    type(t_vectorBlock) :: rsolutionVector
    type(t_vectorBlock) :: rrhsVector
    type(t_vectorBlock) :: rdefTempVector
    type(t_vectorBlock) :: ractualRHS
  
    ! A pointer to the solver node for preconditioning or NULL(),
    ! if no preconditioner is used.
    type(t_linsolNode), pointer :: p_rpreconditioner => null()

  end type

  public :: t_linsolDeflGMRESLevelInfo
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  type t_linsolSubnodeDeflGMRES
  
    ! Number of levels level.
    integer :: nlevels = 1
    
    ! Maybe we should apply Gram-Schmidt twice?
    logical :: btwiceGS = LINSOL_GMRES_DEF_TWICE_GS
    
    ! Scale factor for pseudo-residuals
    real(DP) :: dpseudoResScale = 0.0_DP
    
    ! Relaxation parameter.
    real(DP) :: drelax = 1.0_DP
    
    ! Method how to determine the maximum eigenvalue.
    ! =0: Use a fixed maximum eigenvalue defined by dmaxEigenvalue
    ! on each level.
    integer :: cmaxEigenvalMethod = 0
  
    ! A pointer to the level info structures for all the levels in multigrid
    type(t_linsolDeflGMRESLevelInfo), dimension(:), pointer :: p_RlevelInfo => null()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    type(t_filterChain), dimension(:),pointer :: p_RfilterChain => null()
  
    ! A temp vector for the prolongation/restriction.
    ! Memory for this vector is allocated in initStructure and released
    ! in doneStructure.
    type(t_vectorScalar) :: rprjTempVector
    
    ! Defines whether a preconditioner is used as left or right
    ! preconditioner.
    ! =false: The preconditoner is uzsed as left preconditioner. This is
    !    the default setting.
    ! =true: The preconditioner is used as right preconditioner.
    logical :: brightprec = .false.

  end type
  
  public :: t_linsolSubnodeDeflGMRES
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the SOR solver.
  type t_linsolSubnodeBlockUpwGS
  
    ! Velocity vector
    real(DP), dimension(2) :: Dvel
    
    ! Pointer to the normals of the edges in the triangulation
    real(dp), dimension(:,:), pointer :: p_Dnormals
    
    ! Number of elements in the triangulation
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! List of elements in sorted in downwind direction
    integer, dimension(:), pointer :: p_Iqueue
    
    ! List of elements upwind to one element
    integer, dimension(:,:), pointer :: p_IupwElems
    
    
  end type
  
  public :: t_linsolSubnodeBlockUpwGS

!</typeblock>

!</types>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

contains

  ! ***************************************************************************

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_setMatrices (rsolverNode,Rmatrices)
  
!<description>
  
  ! Initialises the system matrices in the solver node rsolverNode and
  ! in all nodes attached to it (preconditioners,...). For this purpose,
  ! the corresponding initialisation routine of each solver is called.
  !
  ! The last element in Rmatrices defines the system matrix of the
  ! linear system that is to solve. If there is a solver involved
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
  type(t_matrixBlock), dimension(:), intent(in) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)             :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Copy the matrix structure on the finest level to the rsolverNode
    ! structure. This corresponds to the system we want to solve.
    call lsysbl_duplicateMatrix (Rmatrices(ubound(Rmatrices,1)), &
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
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_DEFCORR)
      call linsol_setMatrixDefCorr (rsolverNode,Rmatrices)
    case (LINSOL_ALG_UMFPACK4)
      ! UMFPACK needs no matrix initialisation routine.
    case (LINSOL_ALG_VANKA)
      ! VANKA needs no matrix initialisation routine.
    case (LINSOL_ALG_SPSOR)
      ! SPSOR needs no matrix initialisation routine.
    case (LINSOL_ALG_CG)
      call linsol_setMatrixCG (rsolverNode, Rmatrices)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_setMatrixBiCGStab (rsolverNode, Rmatrices)
    case (LINSOL_ALG_GMRES)
      call linsol_setMatrixGMRES (rsolverNode, Rmatrices)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_setMatrixMultigrid (rsolverNode, Rmatrices)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_setMatrixMultigrid2 (rsolverNode, Rmatrices)
    case (LINSOL_ALG_SCHUR)
      call linsol_setMatrixSchur (rsolverNode, Rmatrices)
    case (LINSOL_ALG_SPECDEFL)
      call linsol_setMatrixSpecDefl (rsolverNode, Rmatrices)
    case (LINSOL_ALG_DEFLGMRES)
      call linsol_setMatrixDeflGMRES (rsolverNode, Rmatrices)
    end select

  end subroutine


  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_setMatrix (rsolverNode, rmatrix)
  
!<description>
  
  ! Initialises the system matrix in the solver node rsolverNode and
  ! in all nodes attached to it (preconditioners,...).
  ! This subroutine is a single-level variant of linsol_setMatrices.
  
!</description>
  
!<input>
  ! The system matrix.
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! local variables
  type(t_matrixBlock), dimension(1) :: Rmatrices

    ! temporarily duplicate the matrix into an array of length 1
    call lsysbl_duplicateMatrix(rmatrix, Rmatrices(1),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    
    ! call the array version of this routine
    call linsol_setMatrices(rsolverNode, Rmatrices)

    ! and release the temporary copy
    call lsysbl_releaseMatrix(Rmatrices(1))

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_matricesCompatible (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  
  ! This subroutine checks the matrices in Rmatrices if they are compatible
  ! with the solver node rsolverNode. Only if the matrices are compatiable,
  ! they are allowed to be attached to the solver with linsol_setMatrices --
  ! otherwise one of the subsolvers will halt the program because not being
  ! able to handle one of the matrices.
  !
  ! The last element in Rmatrices defines the system matrix of the
  ! linear system that is to solve. If there is a solver involved
  ! that needs matrices on multiple grids (e.g. multigrid), Rmatrices
  ! must contain the matrices of all levels. In this case, Rmatrices(1)
  ! defines the matrix on the coarsest level to solve, Rmatrices(2)
  ! the matrix of the first refinement etc. up to the level where
  ! to solve the problem.
  !
  ! ccompatible is set to LINSOL_COMP_OK if the matrices are compatible,
  ! i.e. that the solver can handle them. If the matrices are not compatible,
  ! ccompatible is set to a LINSOL_COMP_xxxx-flag that identifies what is
  ! wrong with the matrices.
  !
  ! With the optional parameter CcompatibleDetail, the caller can get more
  ! detailed information about the compatibility of the matrices.
  ! For every level i, CcompatibleDetail(i) returns a LINSOL_COMP_xxxx
  ! flag that tells whether the matrix on that level is compatible
  ! or not.
  ! Note that it is possible to create cases where this routine ALWAYS
  ! fails, e.g. if two incompatible VANKA solvers are used as pre- and
  ! postsmoother, resp., on one leve in multigrid! But at least if pre- and
  ! postsmoother on each MG level are identical, this routine will
  ! successfully determine (and return in CcompatibleDetail) if
  ! everything is ok.
  
!</description>
  
!<input>
  ! Array of system matrices on all levels of the discretisation.
  type(t_matrixBlock), dimension(:), intent(in) :: Rmatrices

  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode
!</input>

!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible
  
  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Depending on the solver type, call the corresponding check-routine
    ! or return the corresponding matrix compatibility flag directly.
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      ! Ask VANKA if the matrices are ok.
      call linsol_matCompatVANKA (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)

    case (LINSOL_ALG_SPSOR)
      ! Ask SP-SOR if the matrices are ok.
      call linsol_matCompatSPSOR (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_JACOBI)
      ! Ask Jacobi if the matrices are ok.
      call linsol_matCompatJacobi (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_SOR)
      ! Ask SOR if the matrices are ok.
      call linsol_matCompatSOR (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_SSOR)
      ! Ask SSOR if the matrices are ok.
      call linsol_matCompatSSOR (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_DEFCORR)
      ! Ask the defect correction and its subsolvers if the matrices are ok.
      call linsol_matCompatDefCorr (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_UMFPACK4)
      ! Ask VANKA if the matrices are ok.
      call linsol_matCompatUMFPACK4 (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
    
    case (LINSOL_ALG_CG)
      ! Ask CG and its subsolvers if the matrices are ok.
      call linsol_matCompatCG (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_BICGSTAB)
      ! Ask BiCGStab and its subsolvers if the matrices are ok.
      call linsol_matCompatBiCGStab (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_GMRES)
      ! Ask GMRES(m) and its subsolvers if the matrices are ok.
      call linsol_matCompatGMRES (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_MULTIGRID)
      ! Ask Multigrid and its subsolvers if the matrices are ok.
      call linsol_matCompatMultigrid (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_MULTIGRID2)
      ! Ask Multigrid and its subsolvers if the matrices are ok.
      call linsol_matCompatMultigrid2 (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)
      
    case (LINSOL_ALG_SCHUR)
      ! Ask Schur and its subsolvers if the matrices are ok.
      call linsol_matCompatSchur (rsolverNode,Rmatrices,ccompatible,CcompatibleDetail)

    case default
      ! Nothing special. Let us assume that the matrices are ok.
      ccompatible = LINSOL_COMP_OK
      CcompatibleDetail(:) = ccompatible
      
    end select

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_setOnelevelMatrixDirect (rsolverNode,rmatrix)
  
!<description>
  ! This is a special-case subroutine for the case that there is only
  ! one level to solve (e.g. when solving directly with UMFPACK).
  ! rmatrix is exactly one matrix of the level where to solve.
!</description>
  
!<input>
  ! System matrices of the discretisation.
  type(t_matrixBlock), intent(in),target :: rmatrix
!</input>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode
!</inputoutput>

!</subroutine>

  ! t_matrixBlock structure array that receives the matrix pointer;
  ! it has length 1 here in this special situation.
  type(t_matrixBlock), dimension(1) :: Rmatrices
    
    ! Copy the matrix pointer.
    Rmatrices(1) = rmatrix
    
    ! Call the standard initialisation routine
    call linsol_setMatrices (rsolverNode,Rmatrices)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initStructure (rsolverNode, ierror, isolverSubgroup)
  
!<description>
  
  ! Initialises the problem structure in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  ! The initialisation of the problem structure allows the solver component
  ! to perform some 'precalculation', e.g. the UMFPACK4 or ILU solver can
  ! perform a symbolical factorisation. The problem structure usually does
  ! not change during a simulation, except when the grid moves e.g..
  
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode
!</inputoutput>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are initialised.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the structure-init routine of the specific solver
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      call linsol_initStructureVANKA (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_SPSOR)
      call linsol_initStructureSPSOR (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_DEFCORR)
      call linsol_initStructureDefCorr (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_UMFPACK4)
      call linsol_initStructureUMFPACK4 (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_ILU0)
      call linsol_initStructureILU0 (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_MILUS1X1)
      ! No structure routine for (M)ILU(s)
    case (LINSOL_ALG_CG)
      call linsol_initStructureCG (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_initStructureBiCGStab (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_GMRES)
      call linsol_initStructureGMRES (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_initStructureMultigrid (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_initStructureMultigrid2 (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_SCHUR)
      call linsol_initStructureSchur (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_SPECDEFL)
      call linsol_initStructureSpecDefl (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_DEFLGMRES)
      call linsol_initStructureDeflGMRES (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_BLOCKUPWGS)
      call linsol_initStructureBlockUpwGS (rsolverNode,ierror,isubgroup)
    end select
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initData (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Initialises the problem data in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  ! The initialisation of the problem structure allows the solver component
  ! to perform some 'precalculation', e.g. the UMFPACK4 or ILU solver can
  ! perform a numerical factorisation. The problem structure usually does
  ! not change during a simulation, except when the grid moves e.g.
!</description>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are initialised.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the data-init routine of the specific solver
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      call linsol_initDataVANKA (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_SPSOR)
      call linsol_initDataSPSOR (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_DEFCORR)
      call linsol_initDataDefCorr (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_UMFPACK4)
      call linsol_initDataUMFPACK4 (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_ILU0)
      call linsol_initDataILU0 (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_MILUS1X1)
      call linsol_initDataMILUs1x1 (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_CG)
      call linsol_initDataCG (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_initDataBiCGStab (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_GMRES)
      call linsol_initDataGMRES (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_initDataMultigrid (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_initDataMultigrid2 (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_SCHUR)
      call linsol_initDataSchur (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_SPECDEFL)
      call linsol_initDataSpecDefl (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_DEFLGMRES)
      call linsol_initDataDeflGMRES (rsolverNode,ierror,isubgroup)
    case (LINSOL_ALG_BLOCKUPWGS)
      call linsol_initDataBlockUpwGS (rsolverNode,ierror,isubgroup)
    end select
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_updateStructure (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  
  ! Reinitialises the problem structure in the solver node rsolverNode
  
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! For now, call the structure-done- and structure-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    call linsol_doneStructure (rsolverNode,isubgroup)
    call linsol_initStructure (rsolverNode,ierror,isubgroup)
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_updateData (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Reinitialises the problem data in the solver node rsolverNode.
!</description>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<inputoutput>
  ! The solver node containing the solver confuguration
  type(t_linsolNode), intent(inout)                     :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>
  
!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! For now, call the data-done- and data-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    call linsol_doneData (rsolverNode, isubgroup)
    call linsol_initData (rsolverNode, isubgroup,ierror)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_doneData (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Releases the problem data in the solver node rsolverNode.
  
!</description>
  
!<inputoutput>
  
  ! The solver node containing the solver confuguration
  type(t_linsolNode), intent(inout)                     :: rsolverNode
  
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, handle solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the data-done routine of the specific solver
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      call linsol_doneDataVANKA (rsolverNode,isubgroup)
    case (LINSOL_ALG_SPSOR)
      !call linsol_doneDataSPSOR (rsolverNode,isubgroup)
    case (LINSOL_ALG_DEFCORR)
      call linsol_doneDataDefCorr (rsolverNode,isubgroup)
    case (LINSOL_ALG_UMFPACK4)
      call linsol_doneDataUMFPACK4 (rsolverNode,isubgroup)
    case (LINSOL_ALG_ILU0)
      !  No data routine for ILU(0)
    case (LINSOL_ALG_MILUS1X1)
      call linsol_doneDataMILUs1x1 (rsolverNode,isubgroup)
    case (LINSOL_ALG_CG)
      call linsol_doneDataCG (rsolverNode,isubgroup)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_doneDataBiCGStab (rsolverNode,isubgroup)
    case (LINSOL_ALG_GMRES)
      call linsol_doneDataGMRES (rsolverNode,isubgroup)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_doneDataMultigrid (rsolverNode,isubgroup)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_doneDataMultigrid2 (rsolverNode,isubgroup)
    case (LINSOL_ALG_SCHUR)
      call linsol_doneDataSchur (rsolverNode,isubgroup)
    case (LINSOL_ALG_SPECDEFL)
      call linsol_doneDataSpecDefl (rsolverNode,isubgroup)
    case (LINSOL_ALG_DEFLGMRES)
      call linsol_doneDataDeflGMRES (rsolverNode,isubgroup)
    case (LINSOL_ALG_BLOCKUPWGS)
      call linsol_doneDataBlockUpwGS (rsolverNode,isubgroup)
    end select

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_doneStructure (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases the problem structure in the solver node rsolverNode.
  ! This is done by calling the appropriate routine of the
  ! actual solver.
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>
  
!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, handle solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the data-done routine of the specific solver
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      call linsol_doneStructureVANKA (rsolverNode,isubgroup)
    case (LINSOL_ALG_SPSOR)
      call linsol_doneStructureSPSOR (rsolverNode,isubgroup)
    case (LINSOL_ALG_DEFCORR)
      call linsol_doneStructureDefCorr (rsolverNode,isubgroup)
    case (LINSOL_ALG_UMFPACK4)
      call linsol_doneStructureUMFPACK4 (rsolverNode,isubgroup)
    case (LINSOL_ALG_ILU0)
      call linsol_doneStructureILU0 (rsolverNode,isubgroup)
    case (LINSOL_ALG_MILUS1X1)
      ! No structure routine for (M)ILU(s)
    case (LINSOL_ALG_CG)
      call linsol_doneStructureCG (rsolverNode,isubgroup)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_doneStructureBiCGStab (rsolverNode,isubgroup)
    case (LINSOL_ALG_GMRES)
      call linsol_doneStructureGMRES (rsolverNode,isubgroup)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_doneStructureMultigrid (rsolverNode,isubgroup)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_doneStructureMultigrid2 (rsolverNode,isubgroup)
    case (LINSOL_ALG_SCHUR)
      call linsol_doneStructureSchur (rsolverNode,isubgroup)
    case (LINSOL_ALG_SPECDEFL)
      call linsol_doneStructureSpecDefl (rsolverNode,isubgroup)
    case (LINSOL_ALG_DEFLGMRES)
      call linsol_doneStructureDeflGMRES (rsolverNode,isubgroup)
    case (LINSOL_ALG_BLOCKUPWGS)
      call linsol_doneStructureBlockUpwGS (rsolverNode,isubgroup)
    end select

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_releaseSolver (p_rsolverNode,bkeepSolverNode)
  
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
  logical, intent(in), optional :: bkeepSolverNode
!</input>

!<inputoutput>
  ! The solver node which is to be released. If the node contains attached
  ! subsolvers (preconditioners, smoothers,...) they are also released.
  ! On return, rsolverNode is NULL().
  type(t_linsolNode), pointer     :: p_rsolverNode
!</inputoutput>
  
!</subroutine>

    if (.not. associated(p_rsolverNode)) then
      
      ! Print a warning message and return
      call output_line ('Solver note not assigned!', &
                        OU_CLASS_WARNING, OU_MODE_STD, 'linsol_releaseSolver')
      return
      
    end if

    ! Depending on the solver type, call the corresponding done-routine
    select case(p_rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      call linsol_doneVANKA (p_rsolverNode)
    case (LINSOL_ALG_SPSOR)
      call linsol_doneSPSOR (p_rsolverNode)
    case (LINSOL_ALG_SOR)
      call linsol_doneSOR (p_rsolverNode)
    case (LINSOL_ALG_BLOCKSOR)
      call linsol_doneBlockSOR (p_rsolverNode)
    case (LINSOL_ALG_SSOR)
      call linsol_doneSSOR (p_rsolverNode)
    case (LINSOL_ALG_DEFCORR)
      call linsol_doneDefCorr (p_rsolverNode)
    case (LINSOL_ALG_UMFPACK4)
      call linsol_doneUMFPACK4 (p_rsolverNode)
    case (LINSOL_ALG_CG)
      call linsol_doneCG (p_rsolverNode)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_doneBiCGStab (p_rsolverNode)
    case (LINSOL_ALG_GMRES)
      call linsol_doneGMRES (p_rsolverNode)
    case (LINSOL_ALG_ILU0)
      call linsol_doneILU0 (p_rsolverNode)
    case (LINSOL_ALG_MILUS1X1)
      call linsol_doneMILUs1x1 (p_rsolverNode)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_doneMultigrid (p_rsolverNode)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_doneMultigrid2 (p_rsolverNode)
    case (LINSOL_ALG_SCHUR)
      call linsol_doneSchur (p_rsolverNode)
    case (LINSOL_ALG_SPECDEFL)
      call linsol_doneSpecDefl (p_rsolverNode)
    case (LINSOL_ALG_DEFLGMRES)
      call linsol_doneDeflGMRES (p_rsolverNode)
    end select
    
    ! Clean up the associated matrix structure.
    ! Of course, the memory of the matrix is not released from memory, because
    ! if there is a matrix attached, it belongs to the application, not to the
    ! solver!
    call lsysbl_releaseMatrix(p_rsolverNode%rsystemMatrix)
    
    ! Finally release the node itself (if we are allowed to).
    ! Deallocate the structure (if we are allowed to), finish.
    if (.not. present(bkeepSolverNode)) then
      deallocate(p_rsolverNode)
    else
      if (.not. bkeepSolverNode) then
        deallocate(p_rsolverNode)
      end if
    end if
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_alterSolver (rsolverNode, ralterConfig)
  
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
  ! Each subsolver then decides on it is own what to do with the configuration,
  ! depending on the command ralterConfig%ccommand.
  !
  ! Example: If ralterConfig%ccommand=LINSOL_ALTER_CHANGEVANKA, all VANKA
  !  solvers and smoothers in the solver will react to the configuration
  !  block and change their type, depending on ralterConfig%Iconfig.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode

  ! A command/configuration block that is passed to all solvers in the
  ! solver tree identified by rsolverNode.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Call the alterConfig-routine of the actual solver
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      call linsol_alterVANKA (rsolverNode, ralterConfig)
    case (LINSOL_ALG_SPSOR)
      ! No alter routine for SP-SOR
    case (LINSOL_ALG_DEFCORR)
      call linsol_alterDefCorr (rsolverNode, ralterConfig)
    case (LINSOL_ALG_UMFPACK4)
      ! No alter routine for UMFPACK
    case (LINSOL_ALG_ILU0)
      ! No alter routine for ILU(0)
    case (LINSOL_ALG_MILUS1X1)
      ! No alter routine for (M)ILU(s)
    case (LINSOL_ALG_CG)
      call linsol_alterCG (rsolverNode, ralterConfig)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_alterBiCGStab (rsolverNode, ralterConfig)
    case (LINSOL_ALG_GMRES)
      call linsol_alterGMRES (rsolverNode, ralterConfig)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_alterMultigrid (rsolverNode, ralterConfig)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_alterMultigrid2 (rsolverNode, ralterConfig)
    end select
  
  end subroutine
  
  ! ***************************************************************************

!<function>
  
  logical function linsol_testConvergence (rsolverNode, dvecNorm, rdef) result(loutput)
  
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
  type(t_linsolNode), intent(in) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks convergence depending on dvecNorm.
  type(t_vectorBlock), intent(in), optional :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector.
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using dvecNorm.
  real(DP), intent(inout) :: dvecNorm
!</inputoutput>
  
!</function>

    ! Calculate the norm of the vector or take the one given
    ! as parameter
    if (present(rdef)) then
      dvecNorm = lsysbl_vectorNorm (rdef,rsolverNode%iresNorm)
    end if
    
    ! Relative difference in the residuals small enough?
    if (rsolverNode%depsDiff .ne. 0.0_DP) then
      if (abs(dvecNorm-rsolverNode%dlastDefect) .lt. &
              rsolverNode%depsDiff*rsolverNode%dlastDefect) then
        loutput = .true.
        return
      end if
    end if

    select case (rsolverNode%istoppingCriterion)
    
    case (LINSOL_STOP_ONEOF)
      ! Iteration stops if either the absolute or the relative criterium holds.
      loutput = .false.
      
      ! Absolute convergence criterion? Check the norm directly.
      if (rsolverNode%depsAbs .ne. 0.0_DP) then
        if (.not. (dvecNorm .gt. rsolverNode%depsAbs)) then
          loutput = .true.
          return
        end if
      end if
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm.
      if (rsolverNode%depsRel .ne. 0.0_DP) then
        if (.not. &
            (dvecNorm .gt. rsolverNode%depsRel * rsolverNode%dinitialDefect)) then
          loutput = .true.
          return
        end if
      end if
    
    case default
      ! Standard stopping criterion.
      ! Iteration stops if both the absolute and the relative criterium holds.
      loutput = .true.
      
      ! Absolute convergence criterion? Check the norm directly.
      if (rsolverNode%depsAbs .ne. 0.0_DP) then
        if (dvecNorm .gt. rsolverNode%depsAbs) then
          loutput = .false.
          return
        end if
      end if
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm.
      if (rsolverNode%depsRel .ne. 0.0_DP) then
        if (dvecNorm .gt. rsolverNode%depsRel * rsolverNode%dinitialDefect) then
          loutput = .false.
          return
        end if
      end if
    end select
    
  end function
  
  ! ***************************************************************************

!<function>
  
  logical function linsol_testDivergence (rsolverNode, dvecNorm, rdef) result(loutput)
  
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
  type(t_linsolNode), intent(in) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks divergence depending on dvecNorm.
  type(t_vectorBlock), intent(in), optional :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector.
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and divergence is tested using dvecNorm.
  real(DP), intent(inout) :: dvecNorm
!</inputoutput>

!</function>

    ! Calculate the norm of the vector if not given
    ! as parameter
    if (present(rdef)) then
      dvecNorm = lsysbl_vectorNorm (rdef,rsolverNode%iresNorm)
    end if
    
    loutput = .false.
    
    ! Absolute divergence criterion? Check the norm directly.
    if (rsolverNode%ddivAbs .ne. SYS_INFINITY_DP) then
     
      ! use NOT here - gives a better handling of special cases like NaN!
      if ( .not. (dvecNorm .le. rsolverNode%ddivAbs)) then
        loutput = .true.
      end if
      
    end if
    
    ! Relative divergence criterion? Multiply with initial residuum
    ! and check the norm.
    if (rsolverNode%depsRel .ne. SYS_INFINITY_DP) then
      if ( .not. (dvecNorm .le. rsolverNode%dinitialDefect*rsolverNode%ddivRel) ) then
        loutput = .true.
      end if
    end if
  
  end function
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_precondDefect (rsolverNode,rd)
  
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
  type(t_linsolNode), intent(inout)                :: rsolverNode
  
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)               :: rd
!</inputoutput>
  
!</subroutine>

  ! a timer
  type(t_timer) :: rtimer

    ! The only condition to this routine is that matrix and vector are compatible!
    call lsysbl_isMatrixCompatible(rd,rsolverNode%rsystemMatrix,.false.)

    ! clear and start timer
    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Select the solver as configured in rsolverNode and let it perform
    ! the actual preconditioning task.
    select case(rsolverNode%calgorithm)
    case (LINSOL_ALG_VANKA)
      call linsol_precVANKA (rsolverNode,rd)
    case (LINSOL_ALG_SPSOR)
      call linsol_precSPSOR (rsolverNode,rd)
    case (LINSOL_ALG_DEFCORR)
      call linsol_precDefCorr (rsolverNode,rd)
    case (LINSOL_ALG_JACOBI)
      call linsol_precJacobi (rsolverNode,rd)
    case (LINSOL_ALG_SOR)
      call linsol_precSOR (rsolverNode,rd)
    case (LINSOL_ALG_SSOR)
      call linsol_precSSOR (rsolverNode,rd)
    case (LINSOL_ALG_UMFPACK4)
      call linsol_precUMFPACK4 (rsolverNode,rd)
    case (LINSOL_ALG_ILU0)
      call linsol_precILU0 (rsolverNode,rd)
    case (LINSOL_ALG_MILUS1x1)
      call linsol_precMILUS1x1 (rsolverNode,rd)
    case (LINSOL_ALG_CG)
      call linsol_precCG (rsolverNode,rd)
    case (LINSOL_ALG_BICGSTAB)
      call linsol_precBiCGStab (rsolverNode,rd)
    case (LINSOL_ALG_GMRES)
      call linsol_precGMRES (rsolverNode,rd)
    case (LINSOL_ALG_MULTIGRID)
      call linsol_precMultigrid (rsolverNode,rd)
    case (LINSOL_ALG_MULTIGRID2)
      call linsol_precMultigrid2 (rsolverNode,rd)
    case (LINSOL_ALG_SCHUR)
      call linsol_precSchur (rsolverNode,rd)
    case (LINSOL_ALG_BLOCKJAC)
      call linsol_precBlockJac (rsolverNode,rd)
    case (LINSOL_ALG_BLOCKSOR)
      call linsol_precBlockSOR (rsolverNode,rd)
    case (LINSOL_ALG_SPECDEFL)
      call linsol_precSpecDefl (rsolverNode,rd)
    case (LINSOL_ALG_DEFLGMRES)
      call linsol_precDeflGMRES (rsolverNode,rd)
    case (LINSOL_ALG_BLOCKUPWGS)
      call linsol_precDefBlockUpwGS (rsolverNode,rd)
    end select

    ! stop timer and update solver time
    call stat_stopTimer(rtimer)
    rsolverNode%dtimeTotal = rsolverNode%dtimeTotal + rtimer%delapsedReal

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_solveAdaptively (rsolverNode,rx,rb,rtemp)
  
!<description>
  ! This routine starts the solution process to solve a linear system of the
  ! form <tex>$Ax=b$</tex>. The solver configuration must be given by rsolverNode.
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
  type(t_vectorBlock), intent(in), target            :: rb
!</input>
  
!<inputoutput>
  ! The solver node containing the solver configuration
  type(t_linsolNode), intent(inout)                  :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  type(t_vectorBlock), intent(inout)                 :: rx
  
  ! A temporary vector of the same size and structure as rx.
  type(t_vectorBlock), intent(inout)                 :: rtemp
!</inputoutput>
  
!</subroutine>

  ! a local timer
  type(t_timer) :: rtimer

  ! current solver time
  real(DP) ::dcurTime

    ! Fetch current solver time; this might be > 0 if the solver is called more than once.
    dcurTime = rsolverNode%dtimeTotal

    ! clear and start timer
    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Method-specific remarks:
    ! The linear system <tex>$ Ax=b $</tex> is reformulated into a one-step defect-correction
    ! approach
    ! <tex>   $$ x  =  x_0  +  A^{-1}  ( b - A x_0 ) $$  </tex>
    ! The standard solver P configured in rsovlverNode above is then used to
    ! solve
    !     <tex> $$ Ay = b-Ax_0 $$ </tex>
    ! and the solution $x$ is then calculated by $x=x+y$.
    
    ! Calculate the defect:
    ! To build <tex>$ (b-Ax) $</tex>, copy the RHS to the temporary vector
    ! and make a matrix/vector multiplication.
    call lsysbl_copyVector (rb,rtemp)
    call lsysbl_blockMatVec (rsolverNode%rsystemMatrix, rx, rtemp, -1.0_DP, 1.0_DP)
    
    ! Call linsol_precondDefect to solve the subproblem <tex>$ Ay = b-Ax $</tex>.
    ! This overwrites rtemp with the correction vector.
    call linsol_precondDefect (rsolverNode,rtemp)

    ! Add the correction vector to the solution vector and release the memory.
    ! In case we have one... If the initial vector was assumed as zero, we do not
    ! have a correction vector, the result is directly in rx.
    
    ! Correct the solution vector: x = x + y
    call lsysbl_vectorLinearComb (rtemp, rx, 1.0_DP, 1.0_DP)

    ! stop timer and update solver time
    call stat_stopTimer(rtimer)
    rsolverNode%dtimeTotal = dcurTime + rtimer%delapsedReal

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_convertToSmoother (rsolverNode,nsmoothingSteps,domega)
  
!<description>
  ! Converts a t_linsolNode to a smoother structure. A smoother is a solver
  ! that performs a fixed number of iterations without respecting any
  ! residuum. nsmoothingSteps is the number of steps the smoother should
  ! perform.
!</description>
  
!<input>
  ! Number of steps the smoother should perform
  integer, intent(in)          :: nsmoothingSteps
  
  ! OPTIONAL: Damping parameter.
  ! The parameter rsolverNode%domega is set to this value in order
  ! to set up the damping.
  real(DP), intent(in), optional :: domega
!</input>
  
!<inputoutput>
  ! Solver node which should be configured as smoother.
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!</subroutine>

    rsolverNode%depsRel = 0.0_DP
    rsolverNode%depsAbs = 0.0_DP
    rsolverNode%nminIterations = nsmoothingSteps
    rsolverNode%nmaxIterations = nsmoothingSteps
    rsolverNode%iresCheck      = NO
    
    if (present(domega)) rsolverNode%domega = domega

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_initSolverGeneral (p_rsolverNode)
  
!<description>
  ! Creates a new solver node p_rsolverNode on the heap with default
  ! settings and returns a pointer to it.
!</description>
  
!<output>
  ! A pointer to a new solver node on the heap.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Allocate the node - the default initialisation will set most of the
    ! parameters in the structure.
    allocate(p_rsolverNode)

  end subroutine
  
! *****************************************************************************
! Routines for the Defect Correction iteration
! *****************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,Rfilter)
  
!<description>
  ! Creates a t_linsolNode solver structure for the defect correction iteration.
  ! The node can be used to directly solve a problem or to be attached
  ! as solver or preconditioner to another solver structure. The node can be
  ! deleted by linsol_releaseSolver.
  !
  ! The defect correction performs nmaxIterations iterations of the type
  !    <tex> $$ x_{n+1}  =  x_n  +  (b-Ax_n) $$ </tex>
  ! with $x_0:=0$.
  ! It is possible to include a damping parameter to this operation by
  ! changing rsolverNode%domega to a value <tex>$ \not =1 $</tex>. In this case, the
  ! defect correction iteration changes to the Richardson iteration
  !
  ! <tex>  $$ x_{n+1}  =  x_n  +  \omega(b-Ax_n) $$  </tex>
  !
  ! By specifying an additional preconditioner, it is possible to perform
  ! the preconditioned defect correction iteration
  !
  ! <tex>  $$ x_{n+1}  =  x_n  +  \omega P^{-1} (b-Ax_n) $$  </tex>
  !
  ! By specifying an additional filter, the defect vector is filtered before
  ! each preconditioner step.
!</description>
  
!<input>
  ! OPTIONAL: A pointer to the solver structure of a solver that should be
  ! used for preconditioning. If not given or set to NULL(), no preconditioning
  ! will be used.
  type(t_linsolNode), pointer, optional   :: p_rpreconditioner
  
  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not present, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
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
    allocate(p_rsolverNode%p_rsubnodeDefCorr)
    
    ! Attach the preconditioner if given.
    if (present(p_rpreconditioner)) then
      p_rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner => p_rpreconditioner
    end if

    ! Attach the filter if given.
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubnodeDefCorr%p_RfilterChain => Rfilter
    end if

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_alterDefCorr (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the Defect correction solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there is a preconditioner attached. If yes, pass the command
    ! structure to that one.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_alterSolver(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,&
          ralterConfig)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDefCorr (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of UMFPACK4 which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Check if there is a preconditioner attached. If yes, release it.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_releaseSolver(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    end if

    ! Release memory if still associated
    call linsol_doneDataDefCorr (rsolverNode, rsolverNode%isolverSubgroup)
    call linsol_doneStructureDefCorr (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the subnode structure
    deallocate(rsolverNode%p_rsubnodeDefCorr)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_matCompatDefCorr (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_matricesCompatible (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,&
          Rmatrices,ccompatible,CcompatibleDetail)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_setMatrixDefCorr (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_setMatrices (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                              Rmatrices)
    end if

  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initStructureDefCorr (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_initStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                                isubgroup,ierror)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Initialisation. In our case: allocate temporary vectors for our data
    ! by using the associated matrix as template.
    ! That vectors are used in the defect correction so save the intermediate
    ! 'solution' vector.
    call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeDefCorr%rtempVector,.false.,.false.,&
          rsolverNode%cdefaultDataType)

    call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeDefCorr%rtempVector2,.false.,.false.,&
          rsolverNode%cdefaultDataType)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initDataDefCorr (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_initData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                            isubgroup,ierror)
    end if
    
    ! Nothing else to do here.
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataDefCorr (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_doneData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                            isubgroup)
    end if
    
    ! Nothing else to do here.
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureDefCorr (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call linsol_doneStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                                isubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Ok, we are the one to be released. We have temp-vectors to be released!
    if (rsolverNode%p_rsubnodeDefCorr%rtempVector%NEQ .ne. 0) then
      call lsysbl_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector)
      call lsysbl_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector2)
    end if
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_precDefCorr (rsolverNode,rd)
  
!<description>
  ! Applies the Defect Correction preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: ite,i,j,niteAsymptoticCVR
  real(DP) :: dres,dfr
  type(t_vectorBlock), pointer :: p_rx,p_rdef
  type(t_linsolNode), pointer :: p_rprecSubnode
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
  ! Damping parameter
  real(DP) :: domega

  ! The queue saves the current residual and the two previous residuals.
  real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
  
  ! The system matrix
  type(t_matrixBlock), pointer :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  integer :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  logical bprec,bfilter
  
  ! The local subnode
  type(t_linsolSubnodeDefCorr), pointer :: p_rsubnode

    ! Status reset
    rsolverNode%iresult = 0
    rsolverNode%icurrentIteration = 0
    rsolverNode%dinitialDefect = 0.0_DP
    rsolverNode%dlastDefect = 0.0_DP
    rsolverNode%dfinalDefect = 0.0_DP
    rsolverNode%dconvergenceRate = 0.0_DP
    rsolverNode%dasymptoticConvergenceRate = 0.0_DP
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeDefCorr
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
    nminIterations = max(rsolverNode%nminIterations,0)
    
    ! Damping parameter
    domega = rsolverNode%domega
      
    ! Use preconditioning? Filtering?
    bprec = associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    bfilter = associated(rsolverNode%p_rsubnodeDefCorr%p_RfilterChain)
    
    if (bprec) then
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    end if
    if (bfilter) then
      p_RfilterChain => p_rsubnode%p_RfilterChain
    end if

    ! Iteration when the residuum is printed:
    niteResOutput = max(1,rsolverNode%niteResOutput)

    ! We want to perform <= nmaxIterations of the form
    !
    ! <tex>   $$ x_{n+1}  =  x_n  +  \omega P^{-1} (b-Ax) $$  </tex>
    !
    ! At first, set up our temporary vectors, which holds the current
    ! 'solution' and the current 'defect'. We already allocated them
    ! during the initialisation phase - now we want to use them!
    p_rx   => rsolverNode%p_rsubnodeDefCorr%rtempVector
    p_rdef => rsolverNode%p_rsubnodeDefCorr%rtempVector2
    
    ! The vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    call lsysbl_assignDiscrIndirect (rd,p_rx)
    call lsysbl_assignDiscrIndirect (rd,p_rdef)
  
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    call lsysbl_clearVector (p_rx)
  
    ! Copy our RHS rd to p_rdef. As the iteration vector is 0, this
    ! is also our initial defect.
    call lsysbl_copyVector(rd,p_rdef)
    if (bfilter) then
      ! Apply the filter chain to the vector
      call filter_applyFilterChainVec (p_rdef, p_RfilterChain)
    end if
    
    ! Get the norm of the residuum
    dres = lsysbl_vectorNorm (p_rdef,rsolverNode%iresNorm)
    !if (.not.((dres .ge. 1E-99_DP) .and. &
    !          (dres .le. 1E99_DP))) dres = 0.0_DP
    if (dres .lt. 1E-99_DP) dres = 0.0_DP

    ! Initialise starting residuum
    rsolverNode%dinitialDefect = dres
    rsolverNode%dlastDefect = 0.0_DP

    ! initialise the queue of the last residuals with RES
    Dresqueue = dres

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.
    if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
     
      ! final defect is 0, as initialised in the output variable above
      call lsysbl_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    else

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('DefCorr: Iteration '// &
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

      ! Perform at most nmaxIterations loops to get a new vector
      do ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call linsol_precondDefect (p_rprecSubnode,p_rdef)
        end if
        
        ! In p_rdef, we now have the current residuum <tex>$P^{-1} (b-Ax)$</tex>.
        ! Add it (damped) to the current iterate p_x to get
        ! <tex>  $$ x  :=  x  +  \omega P^{-1} (b-Ax) $$  </tex>
        call lsysbl_vectorLinearComb (p_rdef ,p_rx,domega,1.0_DP)
        
        ! Calculate the residuum for the next step : (b-Ax)
        call lsysbl_copyVector (rd,p_rdef)
        call lsysbl_blockMatVec (p_rmatrix, p_rx,p_rdef, -1.0_DP,1.0_DP)
        if (bfilter) then
          ! Apply the filter chain to the vector
          call filter_applyFilterChainVec (p_rdef, p_RfilterChain)
        end if
        
        ! Get the norm of the new (final?) residuum
        dfr = lsysbl_vectorNorm (p_rdef,rsolverNode%iresNorm)
     
        ! Shift the queue with the last residuals and add the new
        ! residual to it.
        Dresqueue = eoshift(Dresqueue,1,dfr)

        rsolverNode%dlastDefect = rsolverNode%dfinalDefect
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        if (linsol_testDivergence(rsolverNode,dfr)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('DefCorr: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
          call output_line ('DefCorr: Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
     
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (linsol_testConvergence(rsolverNode,dfr)) exit
          
        end if

        ! print out the current residual
        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('DefCorr: Iteration '// &
              trim(sys_siL(ite,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if

      end do

      ! Set ite to rsolverNode%nmaxIterations to prevent printing
      ! of "rsolverNode%nmaxIterations+1" if the loop was completed.

      if (ite .gt. rsolverNode%nmaxIterations) then
        ! Warning if we did not reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        if ((rsolverNode%ioutputLevel .ge. 0) .and. &
            (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
          call output_line ('DefCorr: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        end if

        ite = rsolverNode%nmaxIterations
      end if

      ! Finish - either with an error or if converged.
      ! Print the last residuum.
      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ite .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('DefCorr: Iteration '// &
            trim(sys_siL(ite,10))//',  !!RES!! = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if

    end if

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    call lsysbl_copyVector (p_rx,rd)
      
    ! Do not calculate anything if the final residuum is out of bounds -
    ! would result in NaN`s,...
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
    
      ! Calculate asymptotic convergence rate
      if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
        i = min(rsolverNode%iiterations,niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / &
            dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
      end if

      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :
      if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/real(rsolverNode%iiterations,DP))
      else
        rsolverNode%dconvergenceRate = 0.0_DP
      end if
      
      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('DefCorr statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'DefCorr: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if
      
    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    end if
  
  end subroutine

! *****************************************************************************
! Routines for the Jacobi solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initJacobi (p_rsolverNode, domega)
  
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
  !
  ! DEPRECATED: Do not use this parameter! If you really need to use a damped
  ! Jacobi preconditioner for whatever reason, directly set
  ! p_rsolverNode%domega to the desired value.
  ! This parameter will be removed in future from this routine!
  real(DP), optional :: domega
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_JACOBI
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
    
    ! No subnode for Jacobi. Only save domega to the structure.
    if (present(domega)) then
    
      ! DEPRECATED!
      call output_line ('domega parameter is deprecated!', &
                        OU_CLASS_WARNING,OU_MODE_STD, 'linsol_initJacobi')
                        
      p_rsolverNode%domega = domega
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatJacobi (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  integer :: i,n
  type(t_matrixScalar), pointer :: p_rblock => null()

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Get the index of the top-level matrix
    n = ubound(Rmatrices,1)
    
    ! Let us make sure the matrix is not rectangular.
    if(Rmatrices(n)%nblocksPerRow .ne. Rmatrices(n)%nblocksPerCol) then
      
      ! Rectangular block matrix
      ccompatible = LINSOL_COMP_ERRMATRECT
      
    else
    
      ! Okay, let us loop through all main diagonal blocks.
      do i = 1, Rmatrices(n)%nblocksPerRow
      
        ! Get the diagonal block
        p_rblock => Rmatrices(n)%RmatrixBlock(i,i)
        
        ! Is the matrix empty?
        if(p_rblock%NEQ .eq. 0) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! Is the scaling factor zero?
        if(p_rblock%dscaleFactor .eq. 0.0_DP) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! What about the matrix type?
        if((p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
           (p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
          ! Matrix type not supported
          ccompatible = LINSOL_COMP_ERRMATTYPE
          exit
        end if
        
        ! This block is okay, continue with next one
      
      end do
    
    end if
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precJacobi (rsolverNode,rd)
  
!<description>
  ! Applies the Jacobi preconditioner <tex>$ D \in A $<tex> to the defect
  ! vector rd and solves <tex>$ Dd_{new} = d $<tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Jacobi solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: ieq, iblock
    
    type(t_matrixScalar), pointer :: p_rmatrix
    integer, dimension(:), pointer :: p_Kdiag
    real(DP) :: dlocOmega
    real(SP) :: flocOmega
    real(DP), dimension(:), pointer :: p_Dvector, p_Dmatrix
    real(SP), dimension(:), pointer :: p_Fvector, p_Fmatrix
    
      ! Status reset
      rsolverNode%iresult = 0

      ! Loop through all blocks. Each block corresponds to one
      ! diagonal block in the matrix.
      do iblock = 1, rd%nblocks
        
        ! Get the matrix
        p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
        
        ! Now we have to make some decisions. At first, which matrix
        ! structure do we have?
        select case (p_rmatrix%cmatrixFormat)
        case (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
          ! In case of structure 9, Kdiagonal points to the diagonal elements.
          ! In case of structure 7, Kld points to the diagonal elements
          if (p_rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9) then
            call lsyssc_getbase_Kdiagonal(p_rmatrix,p_Kdiag)
          else
            call lsyssc_getbase_Kld(p_rmatrix,p_Kdiag)
          end if
        
          ! Get the omega parameter for the matrix. Do not forget the scaling
          ! factor in the matrix!
          dlocomega = rsolverNode%domega / p_rmatrix%dScaleFactor

          ! Now, which data format do we have? Single or double?
          select case (p_rmatrix%cdataType)
          case (ST_DOUBLE)
            ! Get the matrix data arrays
            call lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
            
            ! Take care of the accuracy of the vector
            select case (rd%RvectorBlock(iblock)%cdataType)
            case (ST_DOUBLE)
              ! Get the data array
              call lsyssc_getbase_double (rd%RvectorBlock(iblock),p_Dvector)
              
              ! and multiply all entries with the inverse of the diagonal
              ! of the matrix.
              if(dlocOmega .eq. 1.0_DP) then
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = p_Dvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
                end do
              else
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = dlocOmega * p_Dvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
                end do
              end if
              
            case (ST_SINGLE)
              ! Get the data array
              call lsyssc_getbase_single (rd%RvectorBlock(iblock),p_Fvector)
              
              ! and multiply all entries with the inverse of the diagonal
              ! of the matrix.
              if(dlocOmega .eq. 1.0_DP) then
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = p_Fvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
                end do
              else
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = dlocOmega * p_Fvector(ieq) / p_Dmatrix(p_Kdiag(ieq))
                end do
              end if

            case default
              call output_line ('Unsupported vector precision.', &
                                OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precJacobi')
              call sys_halt()
            end select
            
          case (ST_SINGLE)

            ! Get the matrix data arrays
            call lsyssc_getbase_single (p_rmatrix,p_Fmatrix)
            
            ! Take care of the accuracy of the vector
            select case (rd%RvectorBlock(iblock)%cdataType)
            case (ST_DOUBLE)
              ! Get the data array
              call lsyssc_getbase_double (rd%RvectorBlock(iblock),p_Dvector)
              
              ! and multiply all entries with the inverse of the diagonal
              ! of the matrix.
              if(dlocOmega .eq. 1.0_DP) then
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = p_Dvector(ieq) / p_Fmatrix(p_Kdiag(ieq))
                end do
              else
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = dlocOmega * p_Dvector(ieq) / p_Fmatrix(p_Kdiag(ieq))
                end do
              end if
              
            case (ST_SINGLE)
              ! Get the data array
              call lsyssc_getbase_single (rd%RvectorBlock(iblock),p_Fvector)
              
              ! Multiplication with Omega can be speeded up as we use
              ! sigle-precision only.
              flocOmega = real(dlocOmega,SP)
              
              ! and multiply all entries with the inverse of the diagonal
              ! of the matrix.
              if(flocOmega .eq. 1.0_SP) then
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = p_Fvector(ieq) / p_Fmatrix(p_Kdiag(ieq))
                end do
              else
                do ieq = 1, rd%RvectorBlock(iblock)%NEQ
                  p_Dvector(ieq) = flocOmega * p_Fvector(ieq) / p_Fmatrix(p_Kdiag(ieq))
                end do
              end if

            case default
              call output_line ('Unsupported vector precision.', &
                                OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precJacobi')
              call sys_halt()
            end select

          case default
            call output_line ('Unsupported matrix precision.', &
                              OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precJacobi')
            call sys_halt()
          end select
        
        case default
          call output_line ('Unsupported matrix format.', &
                            OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precJacobi')
          call sys_halt()
        end select
        
      end do
  
  end subroutine

! *****************************************************************************
! Routines for the SOR solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initSOR (p_rsolverNode, drelax)
  
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
  ! For drelax = 1.0_DP this SOR solver is identical to the Gauss-Seidel
  ! algorithm.
!</description>
  
!<input>
  ! OPTIONAL: Relaxation parameter. If not given, 1.0 is used.
  real(DP), intent(in), optional :: drelax

!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_SOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT

    ! Allocate the subnode for SOR.
    allocate(p_rsolverNode%p_rsubnodeSOR)
    
    ! Save the relaxation parameter.
    p_rsolverNode%p_rsubnodeSOR%drelax = 1.0_DP
    if (present(drelax)) then
      p_rsolverNode%p_rsubnodeSOR%drelax = drelax
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatSOR (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  integer :: i,n
  type(t_matrixScalar), pointer :: p_rblock => null()

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Get the index of the top-level matrix
    n = ubound(Rmatrices,1)
    
    ! Let us make sure the matrix is not rectangular.
    if(Rmatrices(n)%nblocksPerRow .ne. Rmatrices(n)%nblocksPerCol) then
      
      ! Rectangular block matrix
      ccompatible = LINSOL_COMP_ERRMATRECT
      
    else
    
      ! Okay, let us loop through all main diagonal blocks.
      do i = 1, Rmatrices(n)%nblocksPerRow
      
        ! Get the diagonal block
        p_rblock => Rmatrices(n)%RmatrixBlock(i,i)
        
        ! Is the matrix empty?
        if(p_rblock%NEQ .eq. 0) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! Is the scaling factor zero?
        if(p_rblock%dscaleFactor .eq. 0.0_DP) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! What about the matrix type?
        if(p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX9) then
          ! Matrix type not supported
          ccompatible = LINSOL_COMP_ERRMATTYPE
          exit
        end if
        
        ! Is the matrix transposed?
        if (iand(p_rblock%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
          ! Matrix mus not be transposed.
          ccompatible = LINSOL_COMP_ERRTRANSPOSED
          exit
        end if
        
        ! This block is okay, continue with next one
      
      end do
    
    end if
        
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneSOR (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the SOR solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of SOR which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    deallocate(rsolverNode%p_rsubnodeSOR)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precSOR (rsolverNode,rd)
  
!<description>
  ! Applies the SOR preconditioner onto the defect vector rd by solving
  ! <tex> $$ x = omega * relax * (D + relax * L)^{-1} * d $$ </tex>
  ! rd will be overwritten by the preconditioned defect vector x.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the SOR solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode

  ! On entry: The defect vector d to be preconditioned.
  ! On exit: The preconditioned defect vector x.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: iblock,jblock
  
  type (t_matrixScalar), pointer :: p_rmatrix
  real(DP), dimension(:), pointer :: p_Dvector, p_Dmatrix
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kdiagonal
  real(DP) :: domega, drelax
  
    ! Status reset
    rsolverNode%iresult = 0

    ! Okay, get omega and the relaxation parameter
    domega = rsolverNode%domega
    drelax = rsolverNode%p_rsubnodeSOR%drelax
   
    ! Loop through all block rows.
    do iblock = 1, rd%nblocks
    
      ! Loop through all blocks in the lower triangular block matrix.
      do jblock = 1, iblock-1
      
        ! Get the matrix
        p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,jblock)
        
        ! Is it empty?
        if(p_rmatrix%NEQ .eq. 0) cycle
        
        ! Otherwise update the defect.
        call lsyssc_scalarMatVec(p_rmatrix, rd%RvectorBlock(jblock), &
                    rd%RvectorBlock(iblock), drelax, 1.0_DP)
        
      end do
    
      ! Get the main diagonal matrix
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
      
      ! Now we have to make some decisions. At first, which matrix
      ! structure do we have?
      select case (p_rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9)

        call lsyssc_getbase_Kdiagonal(p_rmatrix,p_Kdiagonal)
        call lsyssc_getbase_Kcol(p_rmatrix,p_Kcol)
        call lsyssc_getbase_Kld(p_rmatrix,p_Kld)
      
        ! Now, which data format do we have? Single or double?
        select case (p_rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Get the matrix data arrays
          call lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          select case (rd%RvectorBlock(iblock)%cdataType)
          case (ST_DOUBLE)
            ! Get the data array
            call lsyssc_getbase_double (rd%RvectorBlock(iblock),p_Dvector)
            
            ! Call the SOR subroutine (see below), do the work.
            call performSOR9_dd (p_Dmatrix,p_Kcol,p_Kld,p_Kdiagonal,&
                domega,drelax,p_Dvector,p_rmatrix%dscaleFactor)
            
          case default
            call output_line ('Unsupported vector precision.', &
                              OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precSOR')
            call sys_halt()
          end select
          
        case default
          call output_line ('Unsupported matrix precision.', &
                            OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precSOR')
          call sys_halt()
        end select
      
      case default
        call output_line ('Unsupported matrix format.', &
                          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precSOR')
        call sys_halt()
      end select
      
    end do
    
  contains
  
    !--------------------------------------------------------------------------
    ! Auxiliary routine: SOR
    ! Matrix format 9, double precision matrix, double precision vector
    
    subroutine performSOR9_dd (DA,Kcol,Kld,Kdiagonal,domega,drelax,Dx,dscale)
    
    ! input: Matrix array
    real(DP), dimension(:), intent(in) :: DA
    
    ! input: column structure
    integer, dimension(:), intent(in) :: Kcol
    
    ! input: row structure
    integer, dimension(:), intent(in) :: Kld
    
    ! input: position of diagonal entries in the matrix
    integer, dimension(:), intent(in) :: Kdiagonal
    
    ! input: Damping parameter; standard value is 1.0.
    real(DP), intent(in) :: domega
    
    ! input: Relaxation parameter;
    real(DP), intent(in) :: drelax
    
    ! Scaling factor of the matrix; usually = 1.0
    real(DP), intent(in) :: dscale
    
    ! input: vector to be preconditioned.
    ! output: preconditioned vector
    real(DP), dimension(:), intent(inout) :: Dx
    
    ! local variables
    integer :: i,j,k
    real(DP) :: daux,dalpha
      
      ! We want to perform the following preconditioning step:
      !
      !    1/relax * (D + relax * L) * x = d
      !
      ! and replace the vector d by x in-situ.
      !
      ! We are doing it this way:
      !
      !   1/relax * (D + relax * L) * x = d
      !
      ! <==>
      !                         i-1
      ! a_ii / relax * x_i  +  \sum a_ij x_j  =  d_i      (i=1,...,n)
      !                         j=1
      ! <==>
      !                                 i-1
      ! x_i :=  relax/a_ii * ( d_i  -  \sum a_ij x_j )    (i=1,...,n)
      !                                 j=1
      !
      ! Additionally, we have to scale the resulting vector by the damping
      ! factor omega. Also do not forget that the a_ij must be scaled by the
      ! matrix` scaling factor.
      
      ! Take care of the damping parameter omega and the scaling factor of
      ! the matrix. We will precalculate an auxiliary scalar alpha that will
      ! be multiplied by d_i in each step - this will do the job.
      dalpha = domega / dscale
      
      ! Loop through all matrix rows
      do i = 1, size(Kld)-1
        
        ! Get the index of the main diagonal entry a_ii
        k = Kdiagonal(i)
        
        ! Loop through the lower triangular row part of A
        daux = 0.0_DP
        do j = Kld(i), k-1
          daux = daux + DA(j)*Dx(Kcol(j))
        end do
        
        ! Calculate x_i
        Dx(i) = drelax*(dalpha*Dx(i) - daux) / DA(k)
        
        ! Remark:
        ! When setting Dx(i) in the line above, daux must NOT be
        ! scaled by domega as one might possibly suggest!
        
      end do
 
    end subroutine
  
  end subroutine

! *****************************************************************************
! Routines for the SSOR solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initSSOR (p_rsolverNode, drelax, bscale)
  
!<description>
  ! Creates a t_linsolNode solver structure for the SSOR solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! OPTIONAL: Relaxation parameter. If not given, 1.0 is used.
  real(DP), intent(in), optional :: drelax

  ! OPTIONAL: If set to TRUE, the solution is scaled by 1/(drelax*(2-drelax))
  ! which gives the SSOR preconditioner in the literature. If not existent
  ! or set to FALSE, no scaling is performed; this is the original
  ! implementation of SSOR in FEAT.
  logical, intent(in), optional :: bscale
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_SSOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
  
    ! Allocate the subnode for SSOR.
    allocate(p_rsolverNode%p_rsubnodeSSOR)
    
    ! Save relaxation parameter
    p_rsolverNode%p_rsubnodeSSOR%drelax = 1.0_DP
    if (present(drelax)) p_rsolverNode%p_rsubnodeSSOR%drelax = drelax
    
    ! Save whether the solution should be scaled or not.
    p_rsolverNode%p_rsubnodeSSOR%bscale = .false.
    if (present(bscale)) p_rsolverNode%p_rsubnodeSSOR%bscale = bscale
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatSSOR (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible
  
  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  integer :: i,n
  type(t_matrixScalar), pointer :: p_rblock => null()

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Get the index of the top-level matrix
    n = ubound(Rmatrices,1)
    
    ! Let us make sure the matrix is not rectangular.
    if(Rmatrices(n)%nblocksPerRow .ne. Rmatrices(n)%nblocksPerCol) then
      
      ! Rectangular block matrix
      ccompatible = LINSOL_COMP_ERRMATRECT
      
    else
    
      ! Okay, let us loop through all main diagonal blocks.
      do i = 1, Rmatrices(n)%nblocksPerRow
      
        ! Get the diagonal block
        p_rblock => Rmatrices(n)%RmatrixBlock(i,i)
        
        ! Is the matrix empty?
        if(p_rblock%NEQ .eq. 0) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! Is the scaling factor zero?
        if(p_rblock%dscaleFactor .eq. 0.0_DP) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! What about the matrix type?
        if(p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX9) then
          ! Matrix type not supported
          ccompatible = LINSOL_COMP_ERRMATTYPE
          exit
        end if
        
        ! Is the matrix transposed?
        if (iand(p_rblock%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
          ! Matrix mus not be transposed.
          ccompatible = LINSOL_COMP_ERRTRANSPOSED
          exit
        end if
        
        ! This block is okay, continue with next one
      
      end do
    
    end if
        
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneSSOR (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the SSOR solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of SSOR which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    deallocate(rsolverNode%p_rsubnodeSSOR)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precSSOR (rsolverNode,rd)
  
!<description>
  ! Applies the SSOR preconditioner onto the defect vector rd by solving
  ! <tex> $$ (D + relax * L) * D^{-1} * (D + relax * U) x = d.$$ </tex>
  ! rd will be overwritten by the preconditioned defect vector x.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the SSOR solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode

  ! On entry: The defect vector d to be preconditioned.
  ! On exit: The preconditioned defect vector x.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: iblock
  logical :: bscale
  
  type (t_matrixScalar), pointer :: p_rmatrix
  real(DP), dimension(:), pointer :: p_Dvector, p_Dmatrix
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kdiagonal
  real(DP) :: domega, drelax
    
    ! Status reset
    rsolverNode%iresult = 0

    ! Get bscale from the solver structure
    bscale = rsolverNode%p_rsubnodeSSOR%bscale
    
    ! And get the damping and relaxation parameters
    domega = rsolverNode%domega
    drelax = rsolverNode%p_rsubnodeSSOR%drelax
    
    ! Loop through all diagonal blocks.
    do iblock = 1,rd%nblocks
      ! Get the matrix
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
      
      ! Now we have to make some decisions. At first, which matrix
      ! structure do we have?
      select case (p_rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9)

        call lsyssc_getbase_Kdiagonal(p_rmatrix,p_Kdiagonal)
        call lsyssc_getbase_Kcol(p_rmatrix,p_Kcol)
        call lsyssc_getbase_Kld(p_rmatrix,p_Kld)
      
        ! Now, which data format do we have? Single or double?
        select case (p_rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Get the matrix data arrays
          call lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          select case (rd%RvectorBlock(iblock)%cdataType)
          case (ST_DOUBLE)
            ! Get the data array
            call lsyssc_getbase_double (rd%RvectorBlock(iblock),p_Dvector)
            
            ! Call the SSOR subroutine (see below), do the work.
            call performSSOR9_dd (p_Dmatrix, p_Kcol, p_Kld, &
                p_Kdiagonal, domega, drelax, bscale, p_Dvector, &
                p_rmatrix%dscaleFactor)
            
          case default
            call output_line ('Unsupported vector precision.', &
                              OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precSSOR')
            call sys_halt()
          end select
          
        case default
          call output_line ('Unsupported matrix precision.', &
                            OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precSSOR')
          call sys_halt()
        end select
      
      case default
        call output_line ('Unsupported matrix format.', &
                          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precSSOR')
        call sys_halt()
      end select
      
    end do
    
  contains
  
    !--------------------------------------------------------------------------
    ! Auxiliary routine: SSOR
    ! Matrix format 9, double precision matrix, double precision vector
    
    subroutine performSSOR9_dd (DA,Kcol,Kld,Kdiagonal,domega,&
                                drelax,bscale,Dx,dscale)
    
    ! input: Matrix array
    real(DP), dimension(:), intent(in) :: DA
    
    ! input: column structure
    integer, dimension(:), intent(in) :: Kcol
    
    ! input: row structure
    integer, dimension(:), intent(in) :: Kld
    
    ! input: position of diagonal entries in the matrix
    integer, dimension(:), intent(in) :: Kdiagonal
    
    ! input: Damping parameter; standard value is 1.0
    real(DP), intent(in) :: domega
    
    ! input: Relaxation parameter
    real(DP), intent(in) :: drelax

    ! input: Scaling factor of the matrix; usually = 1.0
    real(DP), intent(in) :: dscale
    
    ! input: Whether the solution should be scaled as suggested in
    ! the literature
    logical, intent(in) :: bscale
    
    ! input: vector to be preconditioned.
    ! output: preconditioned vector
    real(DP), dimension(:), intent(inout) :: Dx
    
    ! local variables
    integer :: i,j,k,n
    real(DP) :: daux,dalpha,dtheta
      
      ! Get the size of the system
      n = size(Kld)-1

      ! We want to perform the following preconditioning step:
      !
      ! <tex>             $$ (D + relax * L) * D^{-1} * (D + relax * U) * x = d $$
      !
      ! $$ \Leftrightarrow   (D + relax * L) * (I + relax * D^{-1} * U) * x = d $$ </tex>
      !
      ! Now the literature comes up with the brilliant idea to scale the
      ! result x of the equation above by the factor 1/(relax*(2-relax)),
      ! see e.g.:
      !  [Barrett et al., `Templates for the Solution of Linear systems:
      !   Building Blocks for Iterative Methods`, p. 42,
      !   http://www.netlib.org/linalg/html_templates/Templates.html]
      !
      ! In the default case, this implementation does not perform this
      ! scaling, however the SSOR preconditioner can be set to this behaviour
      ! by setting the bscale parameter to .true..
      ! Interestingly, not scaling the result gives better results when
      ! using SSOR as a preconditioner for BiCGStab and it also corresponds
      ! to the SSOR implementation used in FEAT 1.x.
      !
      ! For now we are going to ignore the optional scaling factor (and also
      ! the damping factor omega) as we will take care of this later.
      !
      ! So in the first step, we are going to solve the following system
      ! to get an auxiliary vector y:
      !
      ! (D + relax * L) * y = d
      !
      ! <==>
      !                       i-1
      ! a_ii * y_i + relax * \sum a_ij * y_j = d_i         (i=1,...,n)
      !                       j=1
      ! <==>
      !                                 i-1
      ! y_i := 1/a_ii * (d_i - relax * \sum a_ij * y_j)    (i=1,...,n)
      !                                 j=1
      !
      ! Remark:
      ! Actually we are going to overwrite the defect vector d by our
      ! auxiliary vector y in-situ.
      
      ! Take care of the scaling factor of the matrix:
      dalpha = 1.0_DP / dscale
      
      ! Forward insertion
      do i = 1, n
        
        ! Get the index of the main diagonal entry a_ii
        k = Kdiagonal(i)

        ! Loop through the lower triangular row part of A
        daux = 0.0_DP
        do j = Kld(i), k-1
          daux = daux + DA(j)*Dx(Kcol(j))
        end do

        ! Calculate y_i
        Dx(i) = (dalpha*Dx(i) - drelax*daux) / DA(k)
        
      end do
      
      ! What we now have to solve is the following system to get another
      ! (auxiliary) vector z:
      !
      ! <tex> $$ (I + relax * D^{-1} * U) * z = y $$ </tex>
      !
      ! <==>
      !                       n
      ! z_i + relax/a_ii * (\sum a_ij * z_j) = y_i      (i=n,...,1)
      !                     j=i+1
      ! <==>
      !                              n
      ! z_i := y_i - relax/a_ii * (\sum a_ij * z_j)     (i=n,...,1)
      !                            j=i+1
      !
      ! Once again, we are going to overwrite y by z in-situ.
      !
      ! Now the vector z is almost what we wanted - except for the
      ! damping factor omega and the optional scaling factor
      ! 1/(relax*(2-relax)), which we now need to take care of.
      !
      ! Fortunately, we do not need to worry about the scaling factor
      ! of the matrix here, as is appears once for the a_ij inside the
      ! sum and once in 1/a_ii - so the scaling factor of the matrix
      ! simply cancels out!
      !
      ! We are now going to introduce a factor called theta, which is
      ! equal to the damping factor omega if bscale is .false., otherwise
      ! it is equal to omega/(relax*(2-relax)) as suggested in the literature.
      !
      ! So finally we want to calculate:
      !
      ! x := theta * z
      !
      ! Now inserting this into our upper definition of the z_i we get:
      !
      ! x_i := theta * z_i                                            (i=n,...,1)
      !
      ! <==>
      !                                       n
      ! x_i := theta * (y_i - relax/a_ii * (\sum a_ij * z_j))         (i=n,...,1)
      !                                     j=i+1
      ! <==>
      !                                              n
      ! x_i := theta * y_i - theta * relax/a_ii * (\sum a_ij * z_j))  (i=n,...,1)
      !                                            j=i+1
      ! <==>
      !                                      n
      ! x_i := theta * y_i - relax/a_ii * (\sum a_ij * theta * z_j))  (i=n,...,1)
      !                                    j=i+1
      ! <==>
      !                                      n
      ! x_i := theta * y_i - relax/a_ii * (\sum a_ij * x_j))          (i=n,...,1)
      !                                    j=i+1
      !
      ! Now that is it - x is our final preconditioned defect.
      
      
      ! Okay, so let us calculate theta:
      if(bscale) then
        ! Scaled as suggested in literature.
        dtheta = domega / (drelax * (2.0_DP - drelax))
      else
        dtheta = domega
      end if

      ! Loop through all matrix rows in reverse order (backward insertion).
      do i = n, 1, -1

        ! Get the index of the main diagonal entry a_ii
        k = Kdiagonal(i)
        
        ! Loop through the upper triangular row part of A
        daux = 0.0_DP
        do j = Kld(i+1)-1, k+1, -1
          daux = daux + DA(j)*Dx(Kcol(j))
        end do
        
        ! And calculate the final x_i
        Dx(i) = dtheta*Dx(i) - drelax*(daux / DA(k))
      
      end do
      
    end subroutine
  
  end subroutine
  
! *****************************************************************************
! Routines for the VANKA CC2D/CC3D solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initVANKA (p_rsolverNode,domega,csubtypeVANKA)
  
!<description>
  ! Creates a t_linsolNode solver structure for the VANKA solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  !
  ! VANKA is somehow a special type of solver. There is one general VANKA
  ! solver, which applies to many situations, but which is quite slow.
  ! For higher performance, system-specific VANKA subsolvers must be created
  ! with hardcoded treatment of matrices/vectors!
  !
  ! For this purpose, there is the csubtypeVANKA flag. If this flag is present,
  ! it indicates a special VANKA variant for a special type of situation.
  ! the application must ensure then that the correct matrix/vector structure
  ! is used when solving the system with this special subtype, as each subtype
  ! implements a hardwired processing of matrix/vector data!
!</description>
  
!<input>
  ! OPTIONAL: Damping parameter. Is saved to rsolverNode\%domega if specified.
  real(DP), optional :: domega
  
  ! OPTIONAL: VANKA subtype.
  ! If not present, the standard VANKA solver is used.
  ! If present, this is one of the LINSOL_VANKA_xxxx flags that indicate a
  ! special VANKA variant for higher performance.
  integer, intent(in), optional       :: csubtypeVANKA
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_VANKA
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
    
    ! Allocate the subnode for VANKA.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeVANKA)
    
    ! Initialise the VANKA subtype.
    p_rsolverNode%p_rsubnodeVANKA%csubtypeVANKA = LINSOL_VANKA_GENERAL
    if (present(csubtypeVANKA)) then
      p_rsolverNode%p_rsubnodeVANKA%csubtypeVANKA = csubtypeVANKA
    end if
    
    if (present(domega)) p_rsolverNode%domega = domega
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_alterVANKA (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the VANKA solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check the command.
    if (ralterConfig%ccommand .eq. LINSOL_ALTER_CHANGEVANKA) then
    
      ! That is a command to all VANKA solvers. Are we meant?
      ! Check Iconfig(1) which specifies the VANKA subsolver group
      ! to be changed.
      if (ralterConfig%Iconfig(1) .eq. rsolverNode%p_rsubnodeVANKA%rvanka%csubtype) then
      
        ! Oops, we are meant. Set the VANKA variant as specified in Iconfig(2).
        rsolverNode%p_rsubnodeVANKA%rvanka%csubtype = ralterConfig%Iconfig(2)
      
      end if
    
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatVANKA (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in),target   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  integer :: iblock,jblock
  type(t_matrixBlock), pointer :: p_rmat

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK
    
    ! VANKA is a bit tricky. Loop through all scalar submatrices and
    ! check them. The VANKA subtype decides on whether it can use that or not!
    
    p_rmat => Rmatrices(ubound(Rmatrices,1))
    ! VANKA subtype?
    select case (rsolverNode%p_rsubnodeVANKA%csubtypeVANKA)
    case (LINSOL_VANKA_GENERAL,LINSOL_VANKA_GENERALDIRECT)
      
      ! Check all sub-matrices
      do jblock = 1,p_rmat%nblocksPerRow
        do iblock = 1,p_rmat%nblocksPerCol
          if (p_rmat%RmatrixBlock(iblock,jblock)%NEQ .ne. 0) then
            if (iand(p_rmat%RmatrixBlock(iblock,jblock)%imatrixSpec,&
                LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
              ! We ca not handle transposed matrices.
              ccompatible = LINSOL_COMP_ERRTRANSPOSED
            end if
          end if
        end do
      end do

    case (LINSOL_VANKA_2DNAVST  , LINSOL_VANKA_2DNAVSTDIRECT  ,&
          LINSOL_VANKA_2DNAVSTSB, LINSOL_VANKA_2DNAVSTDIRECTSB,&
          LINSOL_VANKA_2DFNAVST , LINSOL_VANKA_2DFNAVSTDIRECT )
      ! Blocks (3,1) and (3,2) must be virtually transposed
      if ((iand(p_rmat%RmatrixBlock(3,1)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .eq. 0) .or. &
          (iand(p_rmat%RmatrixBlock(3,2)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .eq. 0)) then
        ccompatible = LINSOL_COMP_ERRNOTTRANSPOSED
      end if

    case (LINSOL_VANKA_3DNAVST  , LINSOL_VANKA_3DNAVSTDIRECT  ,&
          LINSOL_VANKA_3DFNAVST , LINSOL_VANKA_3DFNAVSTDIRECT )
      ! Blocks (4,1), (4,2) and (4,3) must be virtually transposed
      if ((iand(p_rmat%RmatrixBlock(4,1)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .eq. 0) .or. &
          (iand(p_rmat%RmatrixBlock(4,2)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .eq. 0) .or. &
          (iand(p_rmat%RmatrixBlock(4,3)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .eq. 0)) then
        ccompatible = LINSOL_COMP_ERRNOTTRANSPOSED
      end if
    
    !  ---------------- NEW IMPLEMENTATION ----------------

    case (LINSOL_VANKA_NAVST2D_DIAG, LINSOL_VANKA_NAVST2D_FULL)
      ! Blocks (3,1) and (3,2) must not be virtually transposed
      if ((iand(p_rmat%RmatrixBlock(3,1)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
          (iand(p_rmat%RmatrixBlock(3,2)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
        ccompatible = LINSOL_COMP_ERRTRANSPOSED
      end if

    case (LINSOL_VANKA_BOUSS2D_DIAG, LINSOL_VANKA_BOUSS2D_FULL)
      ! Blocks (3,1) and (3,2) must not be virtually transposed
      if ((iand(p_rmat%RmatrixBlock(3,1)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
          (iand(p_rmat%RmatrixBlock(3,2)%imatrixSpec, &
            LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
        ccompatible = LINSOL_COMP_ERRTRANSPOSED
      end if
      
    end select
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initStructureVANKA (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Memory allocation for the VANKA solver.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANKA solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there is no matrix assigned
    if (rsolverNode%rsystemMatrix%NEQ .eq. 0) then
      call output_line ('No matrix associated!', &
                        OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initStructureVANKA')
      call sys_halt()
    end if
    
    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Allocate a temporary vector
    call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeVANKA%rtempVector,.false.,.false.,&
          rsolverNode%cdefaultDataType)
        
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initDataVANKA (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs final preparation of the VANKA solver subtype.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANKA solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there is no matrix assigned
    if (rsolverNode%rsystemMatrix%NEQ .eq. 0) then
      call output_line ('No matrix associated!', &
                        OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataVANKA')
      call sys_halt()
    end if
    
    ! Which VANKA solver do we actually have?
    select case (rsolverNode%p_rsubnodeVANKA%csubtypeVANKA)
    case (LINSOL_VANKA_GENERAL,LINSOL_VANKA_GENERALDIRECT)
      ! General VANKA for everything
      call vanka_initConformal(rsolverNode%rsystemMatrix,&
                              rsolverNode%p_rsubnodeVANKA%rvanka,&
                              VANKAPC_GENERAL,VANKATP_STANDARD)
                               
    case (LINSOL_VANKA_2DNAVST  , LINSOL_VANKA_2DNAVSTDIRECT  )
      ! Diagonal-type VANKA for Navier-Stokes
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_2DNAVIERSTOKES,VANKATP_DIAGONAL)

    case (LINSOL_VANKA_2DNAVSTSB, LINSOL_VANKA_2DNAVSTDIRECTSB)
      ! Diagonal-type VANKA for Navier-Stokes
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_2DNAVIERSTOKES,VANKATP_DIAGONAL_SOLBASED)
                                
    case (LINSOL_VANKA_2DFNAVST , LINSOL_VANKA_2DFNAVSTDIRECT )
      ! Full VANKA for Navier-Stokes
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_2DNAVIERSTOKES,VANKATP_FULL)

    case (LINSOL_VANKA_2DFNAVSTOC,LINSOL_VANKA_2DFNAVSTOCDIRECT )
      ! Full VANKA for Navier-Stokes optimal control
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_2DNAVIERSTOKESOPTC,VANKATP_FULL)

    case (LINSOL_VANKA_2DFNAVSTOCDIAG,LINSOL_VANKA_2DFNAVSTOCDIAGDIR )
      ! Full VANKA for Navier-Stokes optimal control
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_2DNAVIERSTOKESOPTC,VANKATP_DIAGOPTC)

    case (LINSOL_VANKA_3DNAVST  , LINSOL_VANKA_3DNAVSTDIRECT  )
      ! Diagonal-type VANKA for Navier-Stokes
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_3DNAVIERSTOKES,VANKATP_DIAGONAL)
                                
    case (LINSOL_VANKA_3DFNAVST , LINSOL_VANKA_3DFNAVSTDIRECT )
      ! Full VANKA for Navier-Stokes
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_3DNAVIERSTOKES,VANKATP_FULL)
    
    ! ------------------------- NEW IMPLEMENTATIONS -------------------------
                               
    case (LINSOL_VANKA_NAVST2D_DIAG)
      ! Diagonal-type VANKA for Navier-Stokes
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_NAVIERSTOKES2D,VANKATP_NAVST2D_DIAG)

    case (LINSOL_VANKA_NAVST2D_FULL)
      ! Full VANKA for Navier-Stokes
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_NAVIERSTOKES2D,VANKATP_NAVST2D_FULL)

    case (LINSOL_VANKA_BOUSS2D_DIAG)
      ! Diagonal-type VANKA for Boussinesq
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_BOUSSINESQ2D,VANKATP_BOUSS2D_DIAG)
                               
    case (LINSOL_VANKA_BOUSS2D_FULL)
      ! Full VANKA for Boussinesq
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_BOUSSINESQ2D,VANKATP_BOUSS2D_FULL)

    case (LINSOL_VANKA_2DFNAVSTOCDIAG2,LINSOL_VANKA_2DFNAVSTOCDIAGDIR2 )
      ! Diagonal VANKA for Navier-Stokes optimal control
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_NAVIERSTOKESOPTC2D,VANKATP_NAVSTOPTC2D_DIAG)

    case (LINSOL_VANKA_2DFNAVSTOCFULL2,LINSOL_VANKA_2DFNAVSTOCFULLDIR2 )
      ! Full VANKA for Navier-Stokes optimal control
      call vanka_initConformal (rsolverNode%rsystemMatrix,&
                                rsolverNode%p_rsubnodeVANKA%rvanka,&
                                VANKAPC_NAVIERSTOKESOPTC2D,VANKATP_NAVSTOPTC2D_FULL)

    end select
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneDataVANKA (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases temporary memory of the VANKA solver allocated in
  ! linsol_initDataVANKA
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANKA solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release VANKA.
    call vanka_doneConformal (rsolverNode%p_rsubnodeVANKA%rvanka)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneStructureVANKA (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases temporary memory of the VANKA solver allocated in
  ! linsol_initStrutureVANKA.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANKA solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release the temp vector and associated data if associated
    if (rsolverNode%p_rsubnodeVANKA%rtempVector%NEQ .ne. 0) then
      call lsysbl_releaseVector(rsolverNode%p_rsubnodeVANKA%rtempVector)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneVANKA (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the VANKA solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of VANKA which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  integer :: isubgroup
  
    ! The done routine always releases everything.
    isubgroup = rsolverNode%isolverSubgroup

    ! Release temporary data.
    call linsol_doneDataVANKA (rsolverNode, isubgroup)

    deallocate(rsolverNode%p_rsubnodeVANKA)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precVANKA (rsolverNode,rd)
  
!<description>
  ! Applies the VANKA preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANKA solver
  type(t_linsolNode), intent(inout), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  type (t_vectorBlock), pointer :: p_rvector
  type (t_matrixBlock), pointer :: p_rmatrix
  real(DP) :: domega

    ! Status reset
    rsolverNode%iresult = 0

    ! Getch some information
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Damping parameter
    domega = rsolverNode%domega
    
    ! Get our temporary vector.
    p_rvector => rsolverNode%p_rsubnodeVANKA%rtempVector

    ! The vector shares the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vector.
    call lsysbl_assignDiscrIndirect (rd,p_rvector)
  
    ! Clear our solution vector
    call lsysbl_clearVector (p_rvector)
    
    ! Execute VANKA
    call vanka_conformal (rsolverNode%p_rsubnodeVANKA%rvanka, &
        p_rvector, rd, domega)

    ! Copy the solution vector to rd - it is our preconditioned defect now.
    call lsysbl_copyVector (p_rvector,rd)
  
  end subroutine

  
! *****************************************************************************
! Routines for the SP-SOR solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initSPSOR (p_rsolverNode, csubtype, domega)
  
!<description>
  ! Creates a t_linsolNode solver structure for the SP-SOR solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! SP-SOR subtype. One of the LINSOL_SPSOR_XXXX constants
  integer, intent(in) :: csubtype

  ! OPTIONAL: Damping parameter. Is saved to rsolverNode\%domega if specified.
  real(DP), optional :: domega
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_SPSOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_BLOCK + LINSOL_ABIL_DIRECT
    
    ! Allocate the subnode for SP-SOR.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeSPSOR)
    
    ! Store the subtype.
    p_rsolverNode%p_rsubnodeSPSOR%csubtype = csubtype
    
    if (present(domega)) p_rsolverNode%domega = domega
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatSPSOR (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in):: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), target, intent(in) :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

!  integer :: iblock,jblock
!  type(t_matrixBlock), pointer :: p_rmat

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK
    
!    ! VANKA is a bit tricky. Loop through all scalar submatrices and
!    ! check them. The VANKA subtype decides on whether it can use that or not!
!
!    p_rmat => Rmatrices(ubound(Rmatrices,1))
!
!    ! VANKA subtype?
!    select case (rsolverNode%p_rsubnodeVANKA%csubtypeVANKA)
!    case (LINSOL_VANKA_GENERAL,LINSOL_VANKA_GENERALDIRECT)
!
!      ! Check all sub-matrices
!      do jblock = 1,p_rmat%nblocksPerRow
!        do iblock = 1,p_rmat%nblocksPerCol
!          if (p_rmat%RmatrixBlock(iblock,jblock)%NEQ .ne. 0) then
!            if (iand(p_rmat%RmatrixBlock(iblock,jblock)%imatrixSpec,&
!                LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
!              ! We can not handle transposed matrices.
!              ccompatible = LINSOL_COMP_ERRTRANSPOSED
!            end if
!          end if
!        end do
!      end do
!
!    case (LINSOL_VANKA_2DNAVST  , LINSOL_VANKA_2DNAVSTDIRECT  ,&
!          LINSOL_VANKA_2DNAVSTSB, LINSOL_VANKA_2DNAVSTDIRECTSB,&
!          LINSOL_VANKA_2DFNAVST , LINSOL_VANKA_2DFNAVSTDIRECT )
!      ! Blocks (3,1) and (3,2) must be virtually transposed
!      if ((iand(p_rmat%RmatrixBlock(3,1)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .eq. 0) .or. &
!          (iand(p_rmat%RmatrixBlock(3,2)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .eq. 0)) then
!        ccompatible = LINSOL_COMP_ERRNOTTRANSPOSED
!      end if
!
!    case (LINSOL_VANKA_3DNAVST  , LINSOL_VANKA_3DNAVSTDIRECT  ,&
!          LINSOL_VANKA_3DFNAVST , LINSOL_VANKA_3DFNAVSTDIRECT )
!      ! Blocks (4,1), (4,2) and (4,3) must be virtually transposed
!      if ((iand(p_rmat%RmatrixBlock(4,1)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .eq. 0) .or. &
!          (iand(p_rmat%RmatrixBlock(4,2)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .eq. 0) .or. &
!          (iand(p_rmat%RmatrixBlock(4,3)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .eq. 0)) then
!        ccompatible = LINSOL_COMP_ERRNOTTRANSPOSED
!      end if
!
!    !  ---------------- NEW IMPLEMENTATION ----------------
!
!    case (LINSOL_VANKA_NAVST2D_DIAG, LINSOL_VANKA_NAVST2D_FULL, &
!          LINSOL_VANKA_NAVST2D_PDOF, &
!          LINSOL_VANKA_NAVST2D_SPSOR, LINSOL_VANKA_NAVST2D_SPSSOR)
!      ! Blocks (3,1) and (3,2) must not be virtually transposed
!      if ((iand(p_rmat%RmatrixBlock(3,1)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
!          (iand(p_rmat%RmatrixBlock(3,2)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
!        ccompatible = LINSOL_COMP_ERRTRANSPOSED
!      end if
!
!    case (LINSOL_VANKA_BOUSS2D_DIAG, LINSOL_VANKA_BOUSS2D_FULL) !, &
!          !LINSOL_VANKA_BOUSS2D_SPSOR)
!      ! Blocks (3,1) and (3,2) must not be virtually transposed
!      if ((iand(p_rmat%RmatrixBlock(3,1)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
!          (iand(p_rmat%RmatrixBlock(3,2)%imatrixSpec, &
!            LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
!        ccompatible = LINSOL_COMP_ERRTRANSPOSED
!      end if
!
!    end select
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initStructureSPSOR (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Memory allocation for the SP-SOR solver.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the SP-SOR solver
  type(t_linsolNode), target, intent(inout) :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there is no matrix assigned
    if (rsolverNode%rsystemMatrix%NEQ .eq. 0) then
      call output_line ('No matrix associated!', &
                        OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initStructureSPSOR')
      call sys_halt()
    end if
    
    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! And set up the SP-SOR structure
    select case(rsolverNode%p_rsubnodeSPSOR%csubtype)
    case (LINSOL_SPSOR_NAVST2D)
      ! SP-SOR for 2D Navier-Stokes systems
      call spsor_initNavSt2D(rsolverNode%p_rsubnodeSPSOR%rdata, &
                             rsolverNode%rsystemMatrix, 0_I32)

    case (LINSOL_SPSOR_NAVST2D_DIAG)
      ! diagonal SP-SOR for 2D Navier-Stokes systems
      call spsor_initNavSt2D(rsolverNode%p_rsubnodeSPSOR%rdata, &
                             rsolverNode%rsystemMatrix, SPSOR_FLAG_DIAG)
    
    case default
      call output_line ('Invalid SP-SOR subtype!', &
                        OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initStructureSPSOR')
      call sys_halt()
    end select
        
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initDataSPSOR (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs final preparation of the SP-SOR solver subtype.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the SP-SOR solver
  type(t_linsolNode), target, intent(inout) :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there is no matrix assigned
    if (rsolverNode%rsystemMatrix%NEQ .eq. 0) then
      call output_line ('No matrix associated!', &
                        OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataVANKA')
      call sys_halt()
    end if
    
    ! Initialise the SP-SOR data
    call spsor_initData(rsolverNode%p_rsubnodeSPSOR%rdata)
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneStructureSPSOR (rsolverNode, isolverSubgroup)
  
!<description>
  ! Releases temporary memory of the SP-SOR solver allocated in
  ! linsol_initStrutureSPSOR.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the SP-SOR solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release the SP-SOR data structure
    call spsor_done(rsolverNode%p_rsubnodeSPSOR%rdata)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneSPSOR (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the SP-SOR solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of SP-SOR which is to be cleaned up.
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  integer :: isubgroup
  
    ! The done routine always releases everything.
    isubgroup = rsolverNode%isolverSubgroup

    ! Release temporary data.
    call linsol_doneStructureSPSOR (rsolverNode, isubgroup)

    deallocate(rsolverNode%p_rsubnodeSPSOR)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precSPSOR (rsolverNode,rd)
  
!<description>
  ! Applies the SP-SOR preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANKA solver
  type(t_linsolNode), intent(inout), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

    ! Status reset
    rsolverNode%iresult = 0

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (rsolverNode%rsystemMatrix%NEQ .eq. 0) .or. &
        (rsolverNode%rsystemMatrix%NEQ .ne. rd%NEQ) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
      
    end if

    ! Call precond-routine of SP-SOR
    call spsor_precond(rsolverNode%p_rsubnodeSPSOR%rdata, rd, &
                       rsolverNode%domega)
  
  end subroutine
  
! *****************************************************************************
! Routines for the UMFPACK4 solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initUMFPACK4 (rsolverNode)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the UMFPACK4 solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  
!</description>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: rsolverNode
   
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(rsolverNode)
    
    ! Initialise the type of the solver
    rsolverNode%calgorithm = LINSOL_ALG_UMFPACK4
    
    ! Initialise the ability bitfield with the ability of this solver:
    rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_DIRECT
    
    ! Allocate the subnode for UMFPACK4.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(rsolverNode%p_rsubnodeUMFPACK4)
    
    ! Initialise the UMFPACK4 control structure:
    call UMF4DEF(rsolverNode%p_rsubnodeUMFPACK4%Dcontrol)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatUMFPACK4 (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initStructureUMFPACK4 (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Performs a symbolic factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  type(t_matrixScalar), pointer :: p_rmatrix
  type(t_matrixScalar) :: rtempMatrix
!  integer(I32), dimension(:), pointer :: p_Kld32
!  integer(I32), dimension(:), pointer :: p_Kcol32
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_DA
  type(t_matrixBlock), target :: rmatrixLocal
  integer :: idupFlag

    ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
    ! of a call to the solver routines.
    real(DP), dimension(90) :: Dinfo
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there is no matrix assigned
    if (rsolverNode%rsystemMatrix%NEQ .eq. 0) then
      call output_line ('No matrix associated!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initStructureUMFPACK4')
      call sys_halt()
    end if
    
    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Check out that we can handle the matrix.
    if ((rsolverNode%rsystemMatrix%nblocksPerRow .ne. 1) .or. &
        (rsolverNode%rsystemMatrix%nblocksPerCol .ne. 1)) then
      ! We have to create a global matrix first!
      call glsys_assembleGlobal (rsolverNode%rsystemMatrix,rmatrixLocal, &
                                 .true.,.true.)
      p_rmatrix => rmatrixLocal%RmatrixBlock(1,1)
      
      ! We can now modify p_rmatrix without harming the original matrix,
      ! so set the dup-flag for the copy command to rtempMatrix below
      ! to LSYSSC_DUP_SHARE.
      idupFlag = LSYSSC_DUP_SHARE
    
    else
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock (1,1)

      ! As we work with the original matrix, we set idupFlag
      ! to LSYSSC_DUP_COPY. This prevents the original matrix from being
      ! destroyed in the copy-modify-part below.
      idupFlag = LSYSSC_DUP_COPY
    end if

    select case (p_rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      ! Format 9 is exactly the UMFPACK matrix.
      ! Make a copy of the matrix so we can change its content for UMF4SYM.
      call lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                   idupFlag,idupFlag)
    case (LSYSSC_MATRIX7)
      ! For format 7, we have to modify the matrix slightly.
      ! Make a copy of the matrix:
      call lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                   idupFlag,idupFlag)
      ! Resort the entries to put the diagonal entry to the correct position.
      ! This means: Convert the structure-7 matrix to a structure-9 matrix:
      call lsyssc_convertMatrix (rtempMatrix,LSYSSC_MATRIX9)
    end select
    
    if (p_rmatrix%cdataType .ne. ST_DOUBLE) then
      call output_line ('UMFPACK can only handle double precision matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initStructureUMFPACK4')
      call sys_halt()
    end if

    ! Modify Kcol/Kld of the matrix. Subtract 1 to get the 0-based.
    call lsyssc_addIndex (rtempMatrix%h_Kcol,-1)
    call lsyssc_addIndex (rtempMatrix%h_Kld,-1)
    
    ! Get the data arrays.
    call lsyssc_getbase_double (rtempMatrix,p_DA)

    ! Fill the matrix content by 1.0. That way, UMFPACK will treat
    ! all entries as nonzero.
    call lalg_setVectorDble (p_Da,1.0_DP)
    
    ! Currently, we support only 32 bit integers in UMFPACK;
    ! if necessary, convert!
!    if (KIND(0) .ne. KIND(0_I32)) then
!      call storage_setdatatype(rtempMatrix%h_Kcol,ST_INT32)
!      call storage_setdatatype(rtempMatrix%h_Kld,ST_INT32)
!
!      call storage_getbase_int32 (rtempMatrix%h_Kcol,p_Kcol32)
!      call storage_getbase_int32 (rtempMatrix%h_Kld,p_Kld32)
!
!      ! Perform a symbolic factorization...
!      call UMF4SYM(rtempMatrix%NEQ,rtempMatrix%NEQ,p_Kld32,p_Kcol32,p_Da, &
!                  rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
!                  rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
!                  Dinfo)
!    else
      call storage_getbase_int (rtempMatrix%h_Kcol,p_Kcol)
      call storage_getbase_int (rtempMatrix%h_Kld,p_Kld)
      
      ! Perform a symbolic factorization...
      call UMF4SYM(rtempMatrix%NEQ,rtempMatrix%NEQ,p_Kld,p_Kcol,p_Da, &
                  rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
                  rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
                  Dinfo)
!    end if
                 
    ! Check Dinfo(1) if there is an error
    select case (int(Dinfo(1)))
    case (0)
      ! ok.
    case (1)
      ! Singular matrix
      ierror = LINSOL_ERR_SINGULAR
    case (-1)
      ! no memory
      ierror = LINSOL_ERR_NOMEMORY
    case default
      ! do not know what went wrong
      ierror = LINSOL_ERR_INITERROR
    end select

    ! Throw away the temporary matrix/matrices
    call lsyssc_releaseMatrix (rtempMatrix)
    if ((rsolverNode%rsystemMatrix%nblocksPerCol .ne. 1) .or. &
        (rsolverNode%rsystemMatrix%nblocksPerRow .ne. 1)) then
      call lsysbl_releaseMatrix (rmatrixLocal)
    end if
    
    ! Allocate a temporary vector
    call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeUMFPACK4%rtempVector,.false.,.false.,&
          rsolverNode%cdefaultDataType)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initDataUMFPACK4 (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs a numeric factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  type(t_matrixScalar), pointer :: p_rmatrix
  type(t_matrixBlock), target :: rmatrixLocal
  type(t_matrixScalar) :: rtempMatrix
  real(DP), dimension(:), pointer :: p_DA
  character(LEN=SYS_STRLEN) :: smatrixname

!  integer(I32), dimension(:), pointer :: p_Kld32
!  integer(I32), dimension(:), pointer :: p_Kcol32
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol

  ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
  ! of a call to the solver routines.
  real(DP), dimension(90) :: Dinfo
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there is no matrix assigned
    
    if (rsolverNode%rsystemMatrix%NEQ .eq. 0) then
      call output_line ('No matrix associated!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataUMFPACK4')
      call sys_halt()
    end if
    
    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Check out that we can handle the matrix.
    if ((rsolverNode%rsystemMatrix%nblocksPerCol .ne. 1) .or. &
        (rsolverNode%rsystemMatrix%nblocksPerRow .ne. 1)) then
      ! We have to create a global matrix first!
      call glsys_assembleGlobal (rsolverNode%rsystemMatrix,rmatrixLocal, &
                                 .true.,.true.)
      p_rmatrix => rmatrixLocal%RmatrixBlock(1,1)
    else
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock (1,1)
    end if

    if (p_rmatrix%cdataType .ne. ST_DOUBLE) then
      call output_line ('UMFPACK can only handle double precision matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataUMFPACK4')
      call sys_halt()
    end if

    select case (p_rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      ! Format 9 is exactly the UMFPACK matrix.
      if (p_rmatrix%dScaleFactor .eq. 1.0_DP) then
        ! Make a copy of the matrix structure, but use the same matrix entries.
        call lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                      LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
      else
        ! Make a copy of the whole matrix. We will resolve the scaling factor later,
        ! which changes the entries!
        call lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                      LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
      end if
    case (LSYSSC_MATRIX7)
      ! For format 7, we have to modify the matrix slightly.
      ! Make a copy of the whole matrix:
      call lsyssc_duplicateMatrix (p_rmatrix,rtempMatrix,&
                                   LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
      ! Resort the entries to put the diagonal entry to the correct position.
      ! This means: Convert the structure-7 matrix to a structure-9 matrix:
      call lsyssc_convertMatrix (p_rmatrix,LSYSSC_MATRIX9)
    end select
    
    if (rtempMatrix%dscaleFactor .ne. 1.0_DP) then
      ! The matrix entries have been duplicated above in this case, so we are
      ! allowed to change the entries of rtempMatrix without changing the
      ! original entries. So, 'un'-scale the matrix.
      call lsyssc_scaleMatrix (rtempMatrix,rtempMatrix%dscaleFactor)
      rtempMatrix%dscaleFactor = 1.0_DP
    end if

    ! If the debug flag is set, write out the matrix to a text file.
    if (rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput .gt. 0) then
      if (rsolverNode%ioutputLevel .ge. 0) then
        call output_line ('Writing matrix to a text file.',&
            OU_CLASS_WARNING,OU_MODE_STD,'linsol_initDataUMFPACK4')
      end if
      
      call lsyssc_getbase_double (rtempMatrix,p_DA)
      where (abs(p_Da) .lt. 1.0E-12_DP) p_Da = 0.0_DP
      
      ! Matrix name
      smatrixname = rsolverNode%p_rsubnodeUMFPACK4%smatrixName

      ! Create a matrix name if there is none.
      if (smatrixname .eq. "") then
        smatrixname = "matrix"//&
          trim(sys_siL(rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput,10))//&
          ".txt"
      end if
      
      select case (rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput)
      case (2)
        call matio_writeMatrixMaple (p_rmatrix, &
            trim(rsolverNode%p_rsubnodeUMFPACK4%smatrixName),&
            0, smatrixname, &
            trim(rsolverNode%p_rsubnodeUMFPACK4%saccuracy))
      case default
        call matio_writeMatrixHR (p_rmatrix, &
            trim(rsolverNode%p_rsubnodeUMFPACK4%smatrixName),&
            .true., 0, smatrixname, &
            trim(rsolverNode%p_rsubnodeUMFPACK4%saccuracy))
      end select
    end if
    
    ! DEBUG!!!
    !CALL sys_halt()

    ! Modify Kcol/Kld of the matrix. Subtract 1 to get them 0-based.
    call lsyssc_addIndex (rtempMatrix%h_Kcol,-1)
    call lsyssc_addIndex (rtempMatrix%h_Kld,-1)
    
    ! Get the data arrays.
    call lsyssc_getbase_double (rtempMatrix,p_DA)
    
    ! Currently, we support only 32 bit integers in UMFPACK;
    ! if necessary, convert!
!    if (KIND(0) .ne. KIND(0_I32)) then
!      call storage_setdatatype(rtempMatrix%h_Kcol,ST_INT32)
!      call storage_setdatatype(rtempMatrix%h_Kld,ST_INT32)
!
!      call storage_getbase_int32 (rtempMatrix%h_Kcol,p_Kcol32)
!      call storage_getbase_int32 (rtempMatrix%h_Kld,p_Kld32)
!
!      ! Perform a numeric factorization...
!      call UMF4NUM(p_Kld32,p_Kcol32,p_Da, &
!              rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
!              rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
!              rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
!              Dinfo)
!    else
      call storage_getbase_int (rtempMatrix%h_Kcol,p_Kcol)
      call storage_getbase_int (rtempMatrix%h_Kld,p_Kld)
      
      ! Perform a numeric factorization...
      call UMF4NUM(p_Kld,p_Kcol,p_Da, &
              rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
              rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
              rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
              Dinfo)
!    end if
            
    ! Check Dinfo(1) if there is an error
    select case (int(Dinfo(1)))
    case (0)
      ! ok.
    case (1)
      ! Singular matrix
      ierror = LINSOL_ERR_SINGULAR
    case (-1)
      ! no memory
      ierror = LINSOL_ERR_NOMEMORY
    case (-11)
      ! no memory
      ierror = LINSOL_ERR_MATRIXHASCHANGED
    case default
      ! do not know what went wrong
      ierror = LINSOL_ERR_INITERROR
    end select

    ! Throw away the temporary matrix/matrices
    call lsyssc_releaseMatrix (rtempMatrix)
    if ((rsolverNode%rsystemMatrix%nblocksPerCol .ne. 1) .or. &
        (rsolverNode%rsystemMatrix%nblocksPerRow .ne. 1)) then
      call lsysbl_releaseMatrix (rmatrixLocal)
    end if
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneDataUMFPACK4 (rsolverNode,isolverSubgroup)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release the numerical factorisation if associated
    if (rsolverNode%p_rsubnodeUMFPACK4%inumeric .ne. 0) then
      call UMF4FNUM(rsolverNode%p_rsubnodeUMFPACK4%inumeric)
      rsolverNode%p_rsubnodeUMFPACK4%inumeric = 0
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneStructureUMFPACK4 (rsolverNode, isolverSUbgroup)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release the symbolical factorisation if associated
    if (rsolverNode%p_rsubnodeUMFPACK4%isymbolic .ne. 0) then
      call UMF4FSYM(rsolverNode%p_rsubnodeUMFPACK4%isymbolic)
      rsolverNode%p_rsubnodeUMFPACK4%isymbolic = 0
    end if
    
    ! Release the temp vector if associated
    if (rsolverNode%p_rsubnodeUMFPACK4%rtempVector%NEQ .ne. 0) then
      call lsysbl_releaseVector(rsolverNode%p_rsubnodeUMFPACK4%rtempVector)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneUMFPACK4 (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of UMFPACK4 which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  integer :: isubgroup
  
    ! The done routine always releases everything.
    isubgroup = rsolverNode%isolverSubgroup

    ! Release symbolical and numerical factorisation if still associated...
    call linsol_doneDataUMFPACK4 (rsolverNode, isubgroup)
    call linsol_doneStructureUMFPACK4 (rsolverNode, isubgroup)
    
    deallocate(rsolverNode%p_rsubnodeUMFPACK4)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_precUMFPACK4 (rsolverNode,rd)
  
!<description>
  ! Applies UMFPACK preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>
 
  ! local variables
  integer :: KSYS
  real(DP) :: dres
  real(DP), dimension(:), pointer :: p_Dx,p_Db
  type(t_vectorBlock), pointer :: p_rb
  type(t_linsolSubnodeUMFPACK4), pointer :: p_rsubnode

  ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
  ! of a call to the solver routines.
  real(DP), dimension(90) :: Dinfo
  
    ! Check that our RHS db is double precision - UMFPACK supports only
    ! this.
    if (rd%cdataType .ne. ST_DOUBLE) then
      call output_line ('UMFPACK only supports double precision vectors!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precUMFPACK4')
      call sys_halt()
    end if

    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeUMFPACK4
    p_rb   => p_rsubnode%rtempVector
    
    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    call lsysbl_assignDiscrIndirect (rd,p_rb)
    
    ! Copy the RHS rd to the temp vector; it will be overwritten
    ! by the solution vector
    call lsysbl_copyVector (rd,p_rb)

    ! Get the RHS and solution vector data
    call lsysbl_getbase_double(rd,p_Dx)
    call lsysbl_getbase_double(p_rb,p_Db)

    if (rsolverNode%ioutputLevel .ge. 2) then
      dres = lsysbl_vectorNorm (rd,rsolverNode%iresNorm)
      if (.not.((dres .ge. 1E-99_DP) .and. &
                (dres .le. 1E99_DP))) dres = 0.0_DP

      call output_line ('UMFPACK: !!Initial RES!! = '//trim(sys_sdEL(dres,15)) )
    end if

    ! Solve the system
    ! Solve the system. Note that UMFPACK expects the matrix in
    ! CSR format, which is transposed to our matrix format 9 --
    ! So we solve the transposed system:
    KSYS = 1
    call UMF4SOL(KSYS,p_Dx,p_Db,rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
                 rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,Dinfo)
                
    if (rsolverNode%ioutputLevel .ge. 1) then
      dres = lsysbl_vectorNorm (rd,rsolverNode%iresNorm)
      if (.not.((dres .ge. 1E-99_DP) .and. &
                (dres .le. 1E99_DP))) dres = 0.0_DP

      call output_line ('UMFPACK: !!solution!!    = '//trim(sys_sdEL(dres,15))//&
          ' Status = '//trim(sys_siL(int(Dinfo(1)),10)) )
    end if
                
    ! Check the solver status
    select case (int(Dinfo(1)))
    case (0)
      ! Everything ok.
      rsolverNode%iiterations = 1
      rsolverNode%dconvergenceRate = 0.0_DP
      
      ! Just for algorithms that check the residuum, set the
      ! initial/final residuum values to valid values.
      rsolverNode%dinitialDefect = 1.0_DP
      rsolverNode%dlastDefect = 1.0_DP
      rsolverNode%dfinalDefect = 0.0_DP
      
      ! Scale defect by omega
      call lsysbl_scaleVector(rd, rsolverNode%domega)
      
    case default
      ! We had an error. Do not know which one.
      rsolverNode%iresult = -1
    end select
  
  end subroutine
  
! *****************************************************************************
! Routines for the ILU(0) solver
! *****************************************************************************
! This is a 1-step solver usually used for preconditioning: x_1 = C^-1 x_0.

!<subroutine>
  
  subroutine linsol_initILU0 (rsolverNode)
  
!<description>
  ! Creates a t_linsolNode solver structure for the ILU(0) solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer :: rsolverNode
!</output>

!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(rsolverNode)
    
    ! Normally this solver iterates only once (used e.g. in a preconditioner).
    ! This is the default value. The application can set it higher to
    ! apply it multiple times.
    rsolverNode%nminIterations = 1
    rsolverNode%nmaxIterations = 1
    
    ! Initialise the type of the solver
    rsolverNode%calgorithm = LINSOL_ALG_ILU0
    
    ! Initialise the ability bitfield with the ability of this solver.
    rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK &
                            + LINSOL_ABIL_DIRECT
    
    ! Allocate the subnode for ILU(0).
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(rsolverNode%p_rsubnodeILU0)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatILU0 (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initStructureILU0 (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Okay, are we talking about a scalar matrix here?
    if((rsolverNode%rsystemMatrix%nblocksPerRow .eq. 1) .and. &
       (rsolverNode%rsystemMatrix%nblocksPerCol .eq. 1)) then
      
      ! The system matrix is a scalar matrix - that's fine.
      ! In this case we just need to make a copy of the matrix
      ! structure.
      call lsysbl_duplicateMatrix(rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeILU0%rmatLU, &
          LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
    
    else
    
      ! The system matrix is a non-scalar block matrix.
      ! In this case we'll assemble a global scalar matrix of it.
      call glsys_assembleGlobal(rsolverNode%rsystemMatrix, &
         rsolverNode%p_rsubnodeILU0%rmatLU, .true., .false.)
      
    end if
  
  end subroutine
    
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initDataILU0 (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs a numeric factorisation on the assigned matrix, i.e.
  ! computes the ILU0 decomposition.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiag, p_Iaux
  real(DP), dimension(:), pointer :: p_DLU
  real(DP) :: dmult
  integer :: i,j,k,m,n,idx
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Okay, are we talking about a scalar matrix here?
    if((rsolverNode%rsystemMatrix%nblocksPerRow .eq. 1) .and. &
       (rsolverNode%rsystemMatrix%nblocksPerCol .eq. 1)) then
      
      ! The system matrix is a scalar matrix - that's fine.
      ! So let's copy the matrix content to our LU matrix.
      call lsysbl_duplicateMatrix(rsolverNode%rsystemMatrix, &
          rsolverNode%p_rsubnodeILU0%rmatLU, &
          !LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)
          LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
          
        ! Note: Although we already share the structure of the system matrix,
        ! there is a bug in lsysbl_duplicateMatrix which throws away the
        ! old matrix structure if we specify LSYSSC_DUP_IGNORE for the
        ! structure :(
    
    else
    
      ! The system matrix is a non-scalar block matrix.
      ! In this case we'll assemble a global scalar matrix of it.
      call glsys_assembleGlobal (rsolverNode%rsystemMatrix, &
         rsolverNode%p_rsubnodeILU0%rmatLU, .false., .true.)
    
    end if
    
    ! Now let's fetch the matrix arrays
    call lsyssc_getbase_Kld(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_Kld)
    call lsyssc_getbase_Kcol(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_Kcol)
    call lsyssc_getbase_Kdiagonal(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_Kdiag)
    call lsyssc_getbase_double(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_DLU)
    
    ! ---------- FACTORISATION CODE FOLLOWS ----------
    n = rsolverNode%p_rsubnodeILU0%rmatLU%NEQ
    
    ! allocate an auxiliary array
    allocate(p_Iaux(n))
    
    ! clear the auxiliary array
    call lalg_clearVectorInt(p_Iaux,n)
    
    ! Loop over all rows of the matrix
    do k = 1, n
    
      ! Set up the auxiliary array for this row
      do j = p_Kld(k), p_Kld(k+1)-1
        p_Iaux(p_Kcol(j)) = j
      end do ! j
      
      ! Loop over row k of L
      do j = p_Kld(k), p_Kdiag(k)-1
      
        ! Get the column index
        idx = p_Kcol(j)
        
        ! Calculate L_kj := L_kj / U_jj
        p_DLU(j) = p_DLU(j) * p_DLU(p_Kdiag(idx))
        dmult = p_DLU(j)
        
        do i = p_Kdiag(idx), p_Kld(idx+1)-1
          m = p_Iaux(p_Kcol(i))
          if(m .gt. 0) &
            p_DLU(m) = p_DLU(m) - dmult*p_DLU(i)
        end do ! i
        
      end do ! j
      
      ! Reset auxiliary array
      do j = p_Kld(k), p_Kld(k+1)-1
        p_Iaux(p_Kcol(j)) = 0
      end do ! j
      
      ! Invert main diagonal entry (if possible)
      j = p_Kdiag(k)
      if(abs(p_DLU(j)) .lt. SYS_EPSREAL_DP) then
        
        ! Oops... zero pivot
        ierror = LINSOL_ERR_SINGULAR
        
        ! Exit the main loop
        exit
        
      else
        p_DLU(j) = 1.0_DP / p_DLU(j)
      end if
      
    end do ! k
    
    ! release auxiliary array
    deallocate(p_Iaux)
    
    ! That's it
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneStructureILU0 (rsolverNode,isolverSubgroup)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release the ILU matrix.
    if(rsolverNode%p_rsubnodeILU0%rmatLU%NEQ .gt. 0) &
      call lsysbl_releaseMatrix(rsolverNode%p_rsubnodeILU0%rmatLU)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneILU0 (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of which is to be cleaned up.
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  integer :: isubgroup
  
    ! The done-routine always releases everything.
    isubgroup = rsolverNode%isolverSubgroup
    
    ! Release symbolical and numerical factorisation if still associated...
    call linsol_doneStructureILU0 (rsolverNode, isubgroup)
    
    deallocate(rsolverNode%p_rsubnodeILU0)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_precILU0 (rsolverNode,rd)
  
!<description>
  ! Applies ILU(0) preconditioner <tex>$ LU \approx A $<tex> to the defect
  ! vector rd and solves <tex>$ LUd_{new} = d $<tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the ILU(0) solver
  type(t_linsolNode), intent(inout) :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout) :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Dx, p_DLU
  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiag
  integer :: i,j,k
  real(DP) :: dt

    ! Status reset
    rsolverNode%iresult = 0

    ! Get the data array of the vector
    call lsysbl_getbase_double(rd, p_Dx)
    
    ! Fetch the data arrays of the LU decompositon
    call lsyssc_getbase_Kld(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_Kld)
    call lsyssc_getbase_Kcol(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_Kcol)
    call lsyssc_getbase_Kdiagonal(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_Kdiag)
    call lsyssc_getbase_double(rsolverNode%p_rsubnodeILU0%rmatLU%RmatrixBlock(1,1), p_DLU)
    
    ! Forward insertion
    do i = 2, rd%NEQ
    
      dt = 0.0_DP
      do j = p_Kld(i), p_Kdiag(i)-1
        dt = dt + p_DLU(j)*p_Dx(p_Kcol(j))
      end do ! j
      
      p_Dx(i) = p_Dx(i) - dt
      
    end do ! i
    
    ! Backward insertion
    do i = rd%NEQ, 1, -1
      
      k = p_Kdiag(i)
      
      dt = 0.0_DP
      do j = p_Kld(i+1)-1, k+1, -1
        dt = dt + p_DLU(j)*p_Dx(p_Kcol(j))
      end do ! j
      
      p_Dx(i) = (p_Dx(i) - dt) * p_DLU(k)
      
    end do ! i
    
    ! When the scaling factor is not = 1, scale the vector.
    call lsysbl_scaleVector(rd,rsolverNode%domega)
    
    ! That's it

  end subroutine

! *****************************************************************************
! Routines for the (M)ILUs 1x1 solver
! *****************************************************************************
! This is a 1-step solver usually used for preconditioning: x_1 = C^-1 x_0.
! As we use SPLIB, we have only a numerical factorisation, no symbolical
! factorisation.

!<subroutine>
  
  subroutine linsol_initMILUs1x1 (rsolverNode,ifill,drelax)
  
!<description>
  ! Creates a t_linsolNode solver structure for the (M)ILUs 1x1 solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! Fill-in level for factorisation.
  ! 0=ILU(0), 1=ILU(1),...
  integer, intent(in) :: ifill
  
  ! Relaxation parameter.
  ! 0.0=ILU(s), 1.0=MILU(s), between 0.0 and 1.0=multiplier to use before
  ! adding discarded fill to the diagonal.
  real(DP), intent(in) :: drelax
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: rsolverNode
!</output>

!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(rsolverNode)
    
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
    allocate(rsolverNode%p_rsubnodeMILUs1x1)
    
    ! Save the (M)ILU(s) parameters to the structure
    rsolverNode%p_rsubnodeMILUs1x1%ifill = ifill
    rsolverNode%p_rsubnodeMILUs1x1%drelax = drelax
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatMILUs1x1 (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK

    ! But we cannot handle it if it is not scalar!
    if((Rmatrices(ubound(Rmatrices,1))%nblocksPerCol .ne. 1) .or. &
       (Rmatrices(ubound(Rmatrices,1))%nblocksPerRow .ne. 1)) then
      ccompatible = LINSOL_COMP_ERRNOTSCALAR
    end if
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initDataMILUs1x1 (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Performs a numeric factorisation on the assigned matrix, i.e.
  ! computes the (M)ILU(s) decomposition.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup,maxstr
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_matrixScalar), pointer :: p_rmatrixSc
  integer, dimension(:), pointer :: p_Kld, p_Kcol
  real(DP), dimension(:), pointer :: p_DA
  integer :: ifill,ierr
  real(DP) :: drelax
  
  type(t_MILUdecomp), pointer :: rMILUDecomp
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Stop if there is no matrix assigned
    if (rsolverNode%rsystemMatrix%NEQ .eq. 0) then
      call output_line ('Error: No matrix associated!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataMILUs1x1')
      call sys_halt()
    end if
    
    rMILUDecomp => rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp
    
    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    p_rmatrix => rsolverNode%rsystemMatrix

    ! We only support scalar 1x1 matrices in structure 7 and 9.
    if ((p_rmatrix%nblocksPerCol .ne. 1) .or. &
        (p_rmatrix%nblocksPerRow .ne. 1)) then
      call output_line ('(M)ILU(s) supports only 1x1 matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataMILUs1x1')
      call sys_halt()
    end if

    p_rmatrixSc => p_rmatrix%RmatrixBlock(1,1)

    ! We only support scalar 1x1 matrices in structure 7 and 9.
    if (p_rmatrixSc%cdataType .ne. ST_DOUBLE) then
      call output_line ('(M)ILU(s) supports only double precision matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataMILUs1x1')
      call sys_halt()
    end if
    
    if ((p_rmatrixSc%cmatrixFormat .ne. LSYSSC_MATRIX9) .and.&
        (p_rmatrixSc%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ('(M)ILU(s) supports only structure 7 and 9 matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataMILUs1x1')
      call sys_halt()
    end if
    
    ! Get the matrix description
    call lsyssc_getbase_Kld (p_rmatrixSc,p_Kld)
    call lsyssc_getbase_Kcol (p_rmatrixSc,p_Kcol)
    call lsyssc_getbase_double (p_rmatrixSc,p_DA)

    maxstr = 0

    ! Now calculate the decomposition. Probably allocate more memory if it is
    ! not enough.
    ! Calculate the (M)ILUs matrix using SPLIB.
    ! Get the parameters from the solver node; they were set during the
    ! initialisatin of the solver.
    ifill = rsolverNode%p_rsubnodeMILUs1x1%ifill
    drelax = rsolverNode%p_rsubnodeMILUs1x1%drelax
    
    call iluk_ilu(p_rmatrixSc%NEQ,p_DA,p_Kcol,p_Kld, &
              ifill,drelax,ierr,rMILUDecomp)
              
   
    ! Error?
    select case (ierr)
    case (:-1)
      call output_line ('(M)ILU(s) decomposition singular!', &
          OU_CLASS_WARNING, OU_MODE_STD, 'linsol_initDataMILUs1x1')
    case (0)
      ! everything ok
    case default
      ierror = LINSOL_ERR_INITERROR
      return
    end select
    
    ! now here we can reallocate the jlu array in our structure
    if (rMILUDecomp%h_jlu .ne. ST_NOHANDLE) then
      ! If less than the half of the memory ILU wanted to have is used,
      ! we reallocate the memory. It does not make sense to have that much
      ! waste!
      
      ! nzlu is the number of bytes needed for the jlu array
      ! reallocate if it uses more than twice the amount, ILU wanted
      if (rMILUDecomp%nzlu .lt. rMILUDecomp%isize/2) then
        call storage_realloc('linsol_initDataMILUs1x1', rMILUDecomp%nzlu, &
             rMILUDecomp%h_jlu, ST_NEWBLOCK_ZERO, .true.)
      end if
      
      ! Save 1/scaling factor of the matrix -- to support scaled matrices
      ! when preconditioning.
      rsolverNode%p_rsubnodeMILUs1x1%dscaleFactor = 1.0_DP/p_rmatrixSc%dScaleFactor
    end if
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneDataMILUs1x1 (rsolverNode,isolverSubgroup)
  
!<description>
  ! Releases the memory of teh numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release the ILU matrix.
    !
    if(rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_lu .ne. ST_NOHANDLE) then
      call iluk_freeDecomp(rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneMILUs1x1 (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of UMFPACK4 which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  integer :: isubgroup
  
    ! The done-routine always releases everything.
    isubgroup = rsolverNode%isolverSubgroup
    
    ! Release symbolical and numerical factorisation if still associated...
    call linsol_doneDataMILUs1x1 (rsolverNode, isubgroup)
    
    deallocate(rsolverNode%p_rsubnodeMILUs1x1)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_precMILUS1x1 (rsolverNode,rd)
  
!<description>
  ! Applies (M)ILU(s) preconditioner <tex>$ LU \approx A $<tex> to the defect
  ! vector rd and solves <tex>$ LUd_{new} = d $<tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the (M)ILU(s) solver
  type(t_linsolNode), intent(inout)         :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Dd, p_lu
  integer :: lu,jlu,ilup
  integer, dimension(:), pointer :: p_jlu,p_ilup

    if (rd%cdataType .ne. ST_DOUBLE) then
      call output_line ('(M)ILU(s) only supports double precision vectors!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precMILUS1x1')
      call sys_halt()
    end if
    
    ! Status reset
    rsolverNode%iresult = 0

    ! Get the data array of rd
    call lsysbl_getbase_double (rd,p_Dd)
    
    ! Get MILUs information from the parameter block
    lu = rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_lu
    jlu = rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_jlu
    ilup = rsolverNode%p_rsubnodeMILUs1x1%rMILUdecomp%h_ilup
    
    ! When the scaling factor is not = 1, scale the vector before
    ! preconditioning. This emulates: d = omega(cA)^-1 d = A^-1 (omega*x/c)!
    ! (The value saved in the structure is 1/c!).
    call lsysbl_scaleVector(rd,rsolverNode%domega * &
        rsolverNode%p_rsubnodeMILUs1x1%dscaleFactor)

    ! Solve the system. Call SPLIB, this overwrites the defect vector
    ! with the preconditioned one.
    !CALL lusolt (INT(SIZE(p_Dd),I32),p_Dd, p_Iwork(lu:), &
    !             p_Iwork(jlu:), p_Iwork(ilup:))

    ! Without the following pointers, the INTEL compiler would create temporary
    ! arrays which may lead to a SEGFAULT because of a full stack!

    call storage_getbase_int(jlu, p_jlu)
    call storage_getbase_int(ilup, p_ilup)
    call storage_getbase_double(lu, p_lu)

    call iluk_lusolt (size(p_Dd),p_Dd,p_lu,p_jlu,p_ilup)
                 
  end subroutine
  
! *****************************************************************************
! Routines for the CG solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initCG (p_rsolverNode,p_rpreconditioner,Rfilter)
  
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
  type(t_linsolNode), pointer, optional   :: p_rpreconditioner
  
  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not present, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
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
    allocate(p_rsolverNode%p_rsubnodeCG)
    
    ! Attach the preconditioner if given.
    if (present(p_rpreconditioner)) then
      p_rsolverNode%p_rsubnodeCG%p_rpreconditioner => p_rpreconditioner
    end if

    ! Attach the filter if given.
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubnodeCG%p_RfilterChain => Rfilter
    end if
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_alterCG (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the CG solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there is a preconditioner attached. If yes, pass the command
    ! structure to that one.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_alterSolver(rsolverNode%p_rsubnodeCG%p_rpreconditioner,&
          ralterConfig)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_matCompatCG (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix. This solver can usually handle
    ! everything.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_matricesCompatible ( &
          rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
          Rmatrices,ccompatible,CcompatibleDetail)
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_setMatrixCG (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of CG
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_setMatrices (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                              Rmatrices)
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initStructureCG (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup,i
  type(t_linsolSubnodeCG), pointer :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_initStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                                isubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! CG needs 3 temporary vectors + 1 for preconditioning.
    ! Allocate that here! Use the default data type prescribed in the solver
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeCG
    do i=1,4
      call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%RtempVectors(i),.false.,.false.,&
          rsolverNode%cdefaultDataType)
    end do
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initDataCG (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_initData (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                            isubgroup,ierror)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataCG (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_doneData (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                            isubgroup)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureCG (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup,i
  type(t_linsolSubnodeCG), pointer :: p_rsubnode
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_doneStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                                isubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be released
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeCG
    if (p_rsubnode%RtempVectors(1)%NEQ .ne. 0) then
      do i=4,1,-1
        call lsysbl_releaseVector (p_rsubnode%RtempVectors(i))
      end do
    end if
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneCG (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the CG solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it is also released from the heap by calling
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the CG solver.
  type(t_linsolNode), pointer         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there is a preconditioner attached. If yes, release it.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call linsol_releaseSolver(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    end if
    
    ! Release memory if still associated
    call linsol_doneDataCG (rsolverNode, rsolverNode%isolverSubgroup)
    call linsol_doneStructureCG (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the CG subnode
    deallocate(rsolverNode%p_rsubnodeCG)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_precCG (rsolverNode,rd)
  
!<description>
  ! Applies CG preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP) :: dalpha,dbeta,dgamma,dgammaOld,dres,dfr
  integer :: ite,i,j,niteAsymptoticCVR

  ! The queue saves the current residual and the two previous residuals.
  real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
  
  ! The system matrix
  type(t_matrixBlock), pointer :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  integer :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  logical bprec,bfilter
  
  ! Our structure
  type(t_linsolSubnodeCG), pointer :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  type(t_vectorBlock), pointer :: p_DR,p_DP,p_DD,p_rx
  type(t_linsolNode), pointer :: p_rprecSubnode
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Get some information
    p_rsubnode => rsolverNode%p_rsubnodeCG
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
    nminIterations = max(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?
    bprec = associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    bfilter = associated(rsolverNode%p_rsubnodeCG%p_RfilterChain)
    
    ! Iteration when the residuum is printed:
    niteResOutput = max(1,rsolverNode%niteResOutput)

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
    call lsysbl_assignDiscrIndirect (rd,p_DR)
    call lsysbl_assignDiscrIndirect (rd,p_DP)
    call lsysbl_assignDiscrIndirect (rd,p_DD)
    call lsysbl_assignDiscrIndirect (rd,p_rx)
    
    if (bprec) then
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    end if
    if (bfilter) then
      p_RfilterChain => p_rsubnode%p_RfilterChain
    end if
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    call lsysbl_clearVector (p_rx)
      
    ! Initialise used vectors with zero
    call lsysbl_clearVector(p_DR)
    call lsysbl_clearVector(p_DP)
    call lsysbl_clearVector(p_DD)
    
    ! Initialization
    dalpha = 1.0_DP
    dbeta  = 1.0_DP
    dgamma = 1.0_DP
    dgammaOld = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.
    call lsysbl_copyVector(rd,p_DR)
    if (bfilter) then
      ! Apply the filter chain to the vector
      call filter_applyFilterChainVec (p_DR, p_RfilterChain)
    end if
    
    
    ! Get the norm of the residuum
    dres = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
    if (.not.((dres .ge. 1E-99_DP) .and. &
              (dres .le. 1E99_DP))) dres = 0.0_DP

    ! Initialise starting residuum
    rsolverNode%dinitialDefect = dres
    rsolverNode%dlastDefect = 0.0_DP

    ! initialise the queue of the last residuals with RES
    Dresqueue = dres

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.
    if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
     
      ! final defect is 0, as initialised in the output variable above
      call lsysbl_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    else

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('CG: Iteration '// &
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

      ! Copy the residuum vector p_DR to the preconditioned one.
      call lsysbl_copyVector(p_DR,p_DP)
    
      if (bprec) then
        ! Perform preconditioning with the assigned preconditioning
        ! solver structure.
        call linsol_precondDefect (p_rprecSubnode,p_DP)
      end if

      ! Calculate gamma
      dgamma = lsysbl_scalarProduct (p_DR, p_DP)
      
      ! Perform at most nmaxIterations loops to get a new vector
      
      ! Copy the preconditioned residual vector to the direction vector
      call lsysbl_copyVector (p_DP, p_DD)
      ! And multiply it with -1.0
      call lsysbl_scaleVector (p_DD, -1.0_DP)

      do ite = 1,rsolverNode%nmaxIterations
      
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
        call lsysbl_blockMatVec (p_rmatrix,p_DD,p_DP, 1.0_DP,0.0_DP)

        ! Calculate alpha for the CG iteration
        dalpha = dgamma / lsysbl_scalarProduct (p_DD, p_DP)

        if (dalpha .eq. 0.0_DP) then
          ! We are below machine exactness - we can not do anything more...
          ! May happen with very small problems with very few unknowns!
          if (rsolverNode%ioutputLevel .ge. 2) then
            call output_line('CG: Convergence failed, ALPHA=0!')
          end if
          rsolverNode%iresult = -2
          exit
        end if
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 2: Calculate x_{k+1}
        !
        ! x_{k+1} := x_k - alpha * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call lsysbl_vectorLinearComb (p_DD, p_rx, -1.0_DP * dalpha, 1.0_DP)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 3: Calculate g_{k+1}
        !
        ! g_{k+1} := g_k + alpha * A * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Since we have abused p_DP for saving A * d_k, we can use it now
        ! for the linear combination of g_{k+1}
        call lsysbl_vectorLinearComb (p_DP, p_DR, dalpha, 1.0_DP)

        if (bfilter) then
          ! Apply the filter chain to the new defect vector
          call filter_applyFilterChainVec (p_DR, p_RfilterChain)
        end if

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 4: Calculate ||g_{k+1}|| and write some output
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Get the norm of the new (final?) residuum
        dfr = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        ! Shift the queue with the last residuals and add the new
        ! residual to it.
        Dresqueue = eoshift(Dresqueue,1,dfr)

        rsolverNode%dlastDefect = rsolverNode%dfinalDefect
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        if (linsol_testDivergence(rsolverNode,dfr)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('CG: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
          call output_line('CG: Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
     
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (linsol_testConvergence(rsolverNode,dfr)) exit
          
        end if

        ! print out the current residuum

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('CG: Iteration '// &
              trim(sys_siL(ite,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 5: Calculate h_{k+1}
        !
        ! In other words: Copy the residual vector and apply the
        ! preconditioner, if given.
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Copy the residuum vector p_DR to the preconditioned one.
        call lsysbl_copyVector(p_DR,p_DP)
    
        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call linsol_precondDefect (p_rprecSubnode,p_DP)
        end if
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 6: Calculate gamma
        !
        ! gamma` := gamma
        ! gamma  := < g_{k+1}, h_{k+1} >
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        dgammaOld = dgamma
        dgamma = lsysbl_scalarProduct (p_DR, p_DP)
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 7: Calculate beta
        !
        !          gamma
        ! beta := --------
        !          gamma`
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        dbeta = dgamma / dgammaOld

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 8: Calculate d_{k+1}
        !
        ! d_{k+1} := beta * d_k - h_{k+1}
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call lsysbl_vectorLinearComb (p_DP, p_DD, -1.0_DP, dbeta)
        
        ! That is it - next iteration!
      end do

      ! Set ite to rsolverNode%nmaxIterations to prevent printing of
      ! "rsolverNode%nmaxIterations+1" of the loop was completed
      if (ite .gt. rsolverNode%nmaxIterations) then
        ! Warning if we did not reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        if ((rsolverNode%ioutputLevel .ge. 0) .and. &
            (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
          call output_line ('CG: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        end if

        ite = rsolverNode%nmaxIterations
      end if

      ! Finish - either with an error or if converged.
      ! Print the last residuum.
      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ite .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('CG: Iteration '// &
            trim(sys_siL(ite,10))//',  !!RES!! = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if

    end if

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning. Scale the vector by the damping parameter
    ! in the solver structure.
    call lsysbl_vectorLinearComb (p_rx,rd,rsolverNode%domega,0.0_DP)
      
    ! Do not calculate anything if the final residuum is out of bounds -
    ! would result in NaN`s,...
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
    
      ! Calculate asymptotic convergence rate
      if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
        i = min(rsolverNode%iiterations,niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / &
            dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
      end if

      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :
      if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/real(rsolverNode%iiterations,DP))
      else
        rsolverNode%dconvergenceRate = 0.0_DP
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('CG statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'CG: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if

    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    end if
  
  end subroutine

! *****************************************************************************
! Routines for the BiCGStab solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,Rfilter)
  
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
  type(t_linsolNode), pointer, optional   :: p_rpreconditioner
  
  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not present, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
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
    allocate(p_rsolverNode%p_rsubnodeBiCGStab)
    
    ! Attach the preconditioner if given.
    if (present(p_rpreconditioner)) then
      p_rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner => p_rpreconditioner
    end if

    ! Attach the filter if given.
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubnodeBiCGStab%p_RfilterChain => Rfilter
    end if
      
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_alterBiCGStab (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the BiCGStab solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there is a preconditioner attached. If yes, pass the command
    ! structure to that one.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_alterSolver(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner,&
          ralterConfig)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_matCompatBiCGStab (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix. This solver can usually handle
    ! everything.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible

    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_matricesCompatible ( &
          rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
          Rmatrices,ccompatible,CcompatibleDetail)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_setMatrixBiCGStab (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of BiCGStab
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_setMatrices (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                              Rmatrices)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initStructureBiCGStab (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup,i
  type(t_linsolSubnodeBiCGStab), pointer :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_initStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                                isubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! BiCGStab needs 5 temporary vectors + 1 for preconditioning.
    ! Allocate that here! Use the default data type prescribed in the solver
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    do i=1,6
      call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%RtempVectors(i),.false.,.false.,rsolverNode%cdefaultDataType)
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initDataBiCGStab (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_initData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                            isubgroup,ierror)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_doneData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                            isubgroup)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup,i
  type(t_linsolSubnodeBiCGStab), pointer :: p_rsubnode
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_doneStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                                isubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be released
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    if (p_rsubnode%RtempVectors(1)%NEQ .ne. 0) then
      do i=6,1,-1
        call lsysbl_releaseVector (p_rsubnode%RtempVectors(i))
      end do
    end if
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneBiCGStab (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the BiCGStab solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it is also released from the heap by calling
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the BiCGStab solver.
  type(t_linsolNode), pointer         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there is a preconditioner attached. If yes, release it.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call linsol_releaseSolver(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    end if
    
    ! Release memory if still associated
    call linsol_doneDataBiCGStab (rsolverNode, rsolverNode%isolverSubgroup)
    call linsol_doneStructureBiCGStab (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the BiCGStab subnode
    deallocate(rsolverNode%p_rsubnodeBiCGStab)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine linsol_precBiCGStab (rsolverNode,rd)
  
!<description>
  ! Applies BiCGStab preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
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
  type(t_linsolNode), intent(inout), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP) :: dalpha,dbeta,domega0,domega1,domega2,dres
  real(DP) :: drho1,drho0,dfr
  integer :: ite,i,j,niteAsymptoticCVR

  ! The queue saves the current residual and the two previous residuals.
  real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
  
  ! The system matrix
  type(t_matrixBlock), pointer :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  integer :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  logical bprec,bfilter
  
  ! Our structure
  type(t_linsolSubnodeBiCGStab), pointer :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  type(t_vectorBlock), pointer :: p_DR,p_DR0,p_DP,p_DPA,p_DSA,p_rx
  type(t_linsolNode), pointer :: p_rprecSubnode
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
    nminIterations = max(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?
    bprec = associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    bfilter = associated(rsolverNode%p_rsubnodeBiCGStab%p_RfilterChain)
    
    ! Iteration when the residuum is printed:
    niteResOutput = max(1,rsolverNode%niteResOutput)

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
    call lsysbl_assignDiscrIndirect (rd,p_DR )
    call lsysbl_assignDiscrIndirect (rd,p_DR0)
    call lsysbl_assignDiscrIndirect (rd,p_DP )
    call lsysbl_assignDiscrIndirect (rd,p_DPA)
    call lsysbl_assignDiscrIndirect (rd,p_DSA)
    call lsysbl_assignDiscrIndirect (rd,p_rx)
    
    if (bprec) then
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    end if
    if (bfilter) then
      p_RfilterChain => p_rsubnode%p_RfilterChain
    end if
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    call lsysbl_clearVector (p_rx)
      
    ! Initialise used vectors with zero
    call lsysbl_clearVector(p_DP)
    call lsysbl_clearVector(p_DPA)
    
    ! Initialise the iteration vector with zero.

    ! Initialization
    drho0  = 1.0_DP
    dalpha = 1.0_DP
    domega0 = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.
    call lsysbl_copyVector(rd,p_DR)
    if (bfilter) then
      ! Apply the filter chain to the vector
      call filter_applyFilterChainVec (p_DR, p_RfilterChain)
    end if
    if (bprec) then
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      call linsol_precondDefect (p_rprecSubnode,p_DR)
    end if
    
    ! Get the norm of the residuum
    dres = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
    if (.not.((dres .ge. 1E-99_DP) .and. &
              (dres .le. 1E99_DP))) dres = 0.0_DP

    ! Initialise starting residuum
    rsolverNode%dinitialDefect = dres
    rsolverNode%dlastDefect = 0.0_DP

    ! initialise the queue of the last residuals with RES
    Dresqueue = dres

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.
    if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
     
      ! final defect is 0, as initialised in the output variable above

      call lsysbl_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    else

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('BiCGStab: Iteration '// &
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

      call lsysbl_copyVector(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector
      do ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        drho1 = lsysbl_scalarProduct (p_DR0,p_DR)

        if (drho0*domega0 .eq. 0.0_DP) then
          ! Should not happen
          if (rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('BiCGStab: Iteration prematurely stopped! '//&
                 'Correction vector is zero!')
          end if

          ! Some tuning for the output, then cancel.
          rsolverNode%iresult = -1
          rsolverNode%iiterations = ite-1
          exit
          
        end if

        dbeta=(drho1*dalpha)/(drho0*domega0)
        drho0 = drho1

        call lsysbl_vectorLinearComb (p_DR ,p_DP,1.0_DP,dbeta)
        call lsysbl_vectorLinearComb (p_DPA ,p_DP,-dbeta*domega0,1.0_DP)

        call lsysbl_blockMatVec (p_rmatrix, p_DP,p_DPA, 1.0_DP,0.0_DP)
        
        if (bfilter) then
          ! Apply the filter chain to the vector
          call filter_applyFilterChainVec (p_DPA, p_RfilterChain)
        end if
        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call linsol_precondDefect (p_rprecSubnode,p_DPA)
        end if

        dalpha = lsysbl_scalarProduct (p_DR0,p_DPA)
        
        if (abs(dalpha) .eq. 0.0_DP) then
          ! We are below machine exactness - we can not do anything more...
          ! May happen with very small problems with very few unknowns!
          if (rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('BiCGStab: Convergence failed, ALPHA=0!')
          end if
          rsolverNode%iresult = -2
          exit
        end if
        
        dalpha = drho1/dalpha

        call lsysbl_vectorLinearComb (p_DPA,p_DR,-dalpha,1.0_DP)

        call lsysbl_blockMatVec (p_rmatrix, p_DR,p_DSA, 1.0_DP,0.0_DP)
        
        if (bfilter) then
          ! Apply the filter chain to the vector
          call filter_applyFilterChainVec (p_DSA, p_RfilterChain)
        end if
        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call linsol_precondDefect (p_rprecSubnode,p_DSA)
        end if
        
        domega1 = lsysbl_scalarProduct (p_DSA,p_DR)
        domega2 = lsysbl_scalarProduct (p_DSA,p_DSA)
        
        if (domega1 .eq. 0.0_DP) then
          domega0 = 0.0_DP
        else
          if (domega2 .eq. 0.0_DP) then
            if (rsolverNode%ioutputLevel .ge. 2) then
              call output_line ('BiCGStab: Convergence failed: omega=0!')
            end if
            rsolverNode%iresult = -2
            exit
          end if
          domega0 = domega1/domega2
        end if

        call lsysbl_vectorLinearComb (p_DP ,p_rx,dalpha,1.0_DP)
        call lsysbl_vectorLinearComb (p_DR ,p_rx,domega0,1.0_DP)

        call lsysbl_vectorLinearComb (p_DSA,p_DR,-domega0,1.0_DP)

        ! Get the norm of the new (final?) residuum
        dfr = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        ! Shift the queue with the last residuals and add the new
        ! residual to it.
        Dresqueue = eoshift(Dresqueue,1,dfr)

        rsolverNode%dlastDefect = rsolverNode%dfinalDefect
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        if (linsol_testDivergence(rsolverNode,dfr)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('BiCGStab: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
          call output_line ('BiCGStab: Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
     
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (linsol_testConvergence(rsolverNode,dfr)) exit
          
        end if

        ! print out the current residuum

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('BiCGStab: Iteration '// &
              trim(sys_siL(ite,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if

      end do

      ! Set ite to rsolverNode%nmaxIterations to prevent printing of
      ! "rsolverNode%nmaxIterations+1" of the loop was completed
      if (ite .gt. rsolverNode%nmaxIterations) then
        ! Warning if we did not reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        if ((rsolverNode%ioutputLevel .ge. 0) .and. &
            (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
          call output_line ('BiCGStab: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        end if

        ite = rsolverNode%nmaxIterations
      end if

      ! Finish - either with an error or if converged.
      ! Print the last residuum.
      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ite .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('BiCGStab: Iteration '// &
            trim(sys_siL(ite,10))//',  !!RES!! = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if

    end if

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning. Scale the vector by the damping parameter
    ! in the solver structure.
    call lsysbl_vectorLinearComb (p_rx,rd,rsolverNode%domega,0.0_DP)

    ! Do not calculate anything if the final residuum is out of bounds -
    ! would result in NaN`s,...
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
    
      ! Calculate asymptotic convergence rate
      if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
        i = min(rsolverNode%iiterations,niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / &
            dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
      end if

      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :
      if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/real(rsolverNode%iiterations,DP))
      else
        rsolverNode%dconvergenceRate = 0.0_DP
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('BiCGStab statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'BiCGStab: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if
      
    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    end if
  
  end subroutine

! *****************************************************************************
! Routines for the GMRES(m) solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initGMRES (p_rsolverNode,ikrylovDim,p_rpreconditioner,&
                               Rfilter,btwiceGS,dpseudoResScale)
  
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
  integer :: ikrylovDim
  
  ! OPTIONAL: A pointer to the solver structure of a solver that should be
  ! used for preconditioning. If not given or set to NULL(), no preconditioning
  ! will be used.
  type(t_linsolNode), pointer, optional   :: p_rpreconditioner
  
  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not present, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter

  ! OPTIONAL: A boolean which specifies whether we should apply the Gram-Schmidt
  ! process once or twice per GMRES iteration.
  ! Since the Gram-Schmidt algorithm suffers under numerical instability, we
  ! have a chance to improve its stability by applying the Gram-Schmidt process
  ! twice. If btwiceGS is given and it is .TRUE., then the Gram-Schmidt process
  ! will be applied twice per GMRES iteration, otherwise only once.
  logical, optional :: btwiceGS
  
  ! OPTIONAL: Scale factor for the pseudo-residual-norms to test against the
  ! stopping criterion in the inner GMRES loop. Since the pseudo-residual-norms
  ! may be highly inaccurate, this scale factor can be used to scale the
  ! pseudo-residual-norms before testing against the stopping criterion.
  ! It is recommended to use a scale factor > 1, a factor of 2 should be usually
  ! okay. If dpseudoresscale is not given or set to < 1, it is automatically
  ! set to 2.
  ! Setting dpseudoresscale to 0 will disable the pseudo-resirudal norm check.
  real(DP), optional :: dpseudoResScale
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
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
    allocate(p_rsolverNode%p_rsubnodeGMRES)
    
    ! Save the dimension of the Krylov subspace. This is equal to the number
    ! of GMRES iterations before a restart is done.
    p_rsolverNode%p_rsubnodeGMRES%ikrylovDim = ikrylovDim
    
    ! Save if the Gram-Schmidt process should be applied twice or not.
    if (present(btwiceGS)) then
      p_rsolverNode%p_rsubNodeGMRES%btwiceGS = btwiceGS
    else
      p_rsolverNode%p_rsubNodeGMRES%btwiceGS = LINSOL_GMRES_DEF_TWICE_GS
    end if
    
    ! Save the pseudo-residual-norm scale factor if given, or set
    ! to default scale factor if not given or invalid (i.e. < 1)
    if (present(dpseudoResScale)) then
      if (dpseudoResScale .ge. 0.0_DP) then
        p_rsolverNode%p_rsubNodeGMRES%dpseudoResScale = dpseudoResScale
      else
        p_rsolverNode%p_rsubNodeGMRES%dpseudoResScale = 0.0_DP
      endif
    else
      p_rsolverNode%p_rsubNodeGMRES%dpseudoResScale = 0.0_DP
    endif
    
    ! Attach the preconditioner if given.
    if (present(p_rpreconditioner)) then
      p_rsolverNode%p_rsubNodeGMRES%p_rpreconditioner => p_rpreconditioner
    end if

    ! Attach the filter if given.
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubNodeGMRES%p_RfilterChain => Rfilter
    end if
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_alterGMRES (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the GMRES solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

    ! Check if there is a preconditioner attached. If yes, pass the command
    ! structure to that one.
    if (associated(rsolverNode%p_rsubNodeGMRES%p_rpreconditioner)) then
      call linsol_alterSolver(rsolverNode%p_rsubNodeGMRES%p_rpreconditioner,&
          ralterConfig)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  recursive subroutine linsol_matCompatGMRES (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

    ! Normally, we can handle the matrix. This solver can usually handle
    ! everything.
    ccompatible = LINSOL_COMP_OK

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible

    ! Do we have a preconditioner given? If yes, call the matrix check
    ! routine on that.
    if (associated(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) then
      call linsol_matricesCompatible ( &
          rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
          Rmatrices,ccompatible,CcompatibleDetail)
    end if
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_setMatrixGMRES (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of GMRES(m)
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    if (associated(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) then
      call linsol_setMatrices (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                              Rmatrices)
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initStructureGMRES (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup,i
  integer :: idim
  integer , dimension(2) :: idim2
  type(t_linsolSubnodeGMRES), pointer :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    p_rsubnode => rsolverNode%p_rsubnodeGMRES

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(p_rsubnode%p_rpreconditioner)) then
      call linsol_initStructure (p_rsubnode%p_rpreconditioner, isubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! We now need to check if the dimension of the Krylov subspace is
    ! positive. If it is not, we need to cancel the initialization here.
    if (p_rsubnode%ikrylovDim .le. 0) then
      call output_line ('imension of Krylov subspace for GMRES(m) is <= 0 !', &
          OU_CLASS_ERROR, OU_MODE_STD, 'mysubroutine')
      ierror = LINSOL_ERR_INITERROR
      call sys_halt()
    end if
    
    ! Get the stuff out of our solver node
    idim = p_rsubnode%ikrylovDim
    idim2(1) = idim
    idim2(2) = idim

    ! Now comes the intersting part - we need to allocate a lot of arrays
    ! and vectors for the Krylov subspace for GMRES here.
    
    ! Call our storage to allocate the 1D/2D arrays
    call storage_new('linsol_initStructureGMRES', 'Dh', idim2, &
        ST_DOUBLE, p_rsubnode%hDh, ST_NEWBLOCK_NOINIT)
    call storage_new('linsol_initStructureGMRES', 'Dc', idim, &
        ST_DOUBLE, p_rsubnode%hDs, ST_NEWBLOCK_NOINIT)
    call storage_new('linsol_initStructureGMRES', 'Ds', idim, &
        ST_DOUBLE, p_rsubnode%hDc, ST_NEWBLOCK_NOINIT)
    call storage_new('linsol_initStructureGMRES', 'Dq', idim+1, &
        ST_DOUBLE,  p_rsubnode%hDq, ST_NEWBLOCK_NOINIT)
    
    ! Get the pointers
    call storage_getbase_double2D(p_rsubnode%hDh, p_rsubnode%Dh)
    call storage_getbase_double(p_rsubnode%hDc, p_rsubnode%Dc)
    call storage_getbase_double(p_rsubnode%hDs, p_rsubnode%Ds)
    call storage_getbase_double(p_rsubnode%hDq, p_rsubnode%Dq)
    
    ! Allocate space for our auxiliary vectors
    allocate(p_rsubnode%p_rv(idim+1))
    allocate(p_rsubnode%p_rz(idim))
    
    ! Create them
    do i=1, idim+1
      call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%p_rv(i),.false.,.false.,rsolverNode%cdefaultDataType)
    end do
    do i=1, idim
      call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
          p_rsubnode%p_rz(i),.false.,.false.,rsolverNode%cdefaultDataType)
    end do
    
    ! Create an iteration vector x
    call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
        p_rsubnode%rx,.false.,.false.,rsolverNode%cdefaultDataType)
    
    ! Okay, that is all!
      
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initDataGMRES (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) then
      call linsol_initData (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                            isubgroup,ierror)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataGMRES (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) then
      call linsol_doneData (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                            isubgroup)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureGMRES (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES(m) solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup,i,idim
  type(t_linsolSubnodeGMRES), pointer :: p_rsubnode
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) then
      call linsol_doneStructure (rsolverNode%p_rsubnodeGMRES%p_rpreconditioner, &
                                isubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be released
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeGMRES
    idim = p_rsubnode%ikrylovDim
    
    ! Release auxiliary vectors
    if (associated(p_rsubnode%p_rz)) then
      do i=1, idim
        call lsysbl_releaseVector(p_rsubnode%p_rz(i))
      end do
      do i=1, idim+1
        call lsysbl_releaseVector(p_rsubnode%p_rv(i))
      end do
    
      ! Destroy them
      deallocate(p_rsubnode%p_rz)
      deallocate(p_rsubnode%p_rv)
    endif

    ! Release iteration vector
    if (p_rsubnode%rx%NEQ > 0) then
      call lsysbl_releaseVector(p_rsubnode%rx)
    end if
    
    ! Destroy our 1D/2D arrays
    if (p_rsubnode%hDq .ne. ST_NOHANDLE) then
      call storage_free(p_rsubnode%hDq)
      call storage_free(p_rsubnode%hDs)
      call storage_free(p_rsubnode%hDc)
      call storage_free(p_rsubnode%hDh)
    end if
    
    ! That is all!
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneGMRES (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the GMRES(m) solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it is also released from the heap by calling
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the GMRES(m) solver.
  type(t_linsolNode), pointer         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there is a preconditioner attached. If yes, release it.
    if (associated(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)) then
      call linsol_releaseSolver(rsolverNode%p_rsubnodeGMRES%p_rpreconditioner)
    end if
    
    ! Release memory if still associated
    call linsol_doneDataGMRES (rsolverNode, rsolverNode%isolverSubgroup)
    call linsol_doneStructureGMRES (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the GMRES subnode
    deallocate(rsolverNode%p_rsubnodeGMRES)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_precGMRES (rsolverNode,rd)
  
!<description>
  ! Applies GMRES(m) preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the GMRES solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP) :: dalpha,dbeta,dres,dfr,dtmp,dpseudores,dprnsf
  integer :: ite,i,j,k,niteAsymptoticCVR
  
  ! Here come our 1D/2D arrays
  real(DP), dimension(:), pointer :: p_Dc, p_Ds, p_Dq
  real(DP), dimension(:,:), pointer :: p_Dh

  ! The queue saves the current residual and the two previous residuals.
  real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
  
  ! The system matrix
  type(t_matrixBlock), pointer :: p_rmatrix

  
  ! Minimum number of iterations, print-sequence for residuals, dimension
  ! of krylov subspace
  integer :: nminIterations, niteResOutput, idim, hDq
  
  ! Whether to filter/precondition, apply Gram-Schmidt twice
  logical bprec,bfilter,btwiceGS
  
  ! Our structure
  type(t_linsolSubnodeGMRES), pointer :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  type(t_vectorBlock), pointer :: p_rx
  type(t_vectorBlock), dimension(:), pointer :: p_rv, p_rz
  type(t_linsolNode), pointer :: p_rprecSubnode
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Get some information
    p_rsubnode => rsolverNode%p_rsubnodeGMRES
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
    nminIterations = max(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?
    bprec = associated(p_rsubnode%p_rpreconditioner)
    bfilter = associated(p_rsubnode%p_RfilterChain)
    
    ! Iteration when the residuum is printed:
    niteResOutput = max(1,rsolverNode%niteResOutput)
    
    ! Apply Gram-Schmidt twice?
    btwiceGS = p_rsubnode%btwiceGS
    
    ! Get the dimension of Krylov subspace
    idim = p_rsubnode%ikrylovdim
    
    ! Get our pseudo-residual-norm scale factor
    dprnsf = p_rsubnode%dpseudoResScale
    
    ! Now we need to check the pseudo-residual scale factor.
    ! If it is non-positive, we set it to 0.
    if (dprnsf .le. 0.0_DP) then
      dprnsf = 0.0_DP
    end if

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
    do i=1, idim+1
      call lsysbl_assignDiscrIndirect (rd,p_rv(i))
    end do
    do i=1, idim
      call lsysbl_assignDiscrIndirect (rd,p_rz(i))
    end do
    call lsysbl_assignDiscrIndirect (rd,p_rx)
    
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
    call lsysbl_assignDiscrIndirect (rd, p_rx)
    do i=1,idim
      call lsysbl_assignDiscrIndirect (rd, p_rz(i))
    end do
    do i=1,idim+1
      call lsysbl_assignDiscrIndirect (rd, p_rv(i))
    end do
    
    if (bprec) then
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    end if
    if (bfilter) then
      p_RfilterChain => p_rsubnode%p_RfilterChain
    end if
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    call lsysbl_clearVector (p_rx)
      
    ! Copy our RHS rd to p_rv(1). As the iteration vector is 0, this
    ! is also our initial defect.
    call lsysbl_copyVector(rd,p_rv(1))
    if (bfilter) then
      call filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
    end if
    
    ! Get the norm of the residuum
    ! We need to calculate the norm twice, since we need one for
    ! the stopping criterion (selected by the user) and the
    ! euclidian norm for the internal GMRES algorithm
    dfr = lsysbl_vectorNorm (p_rv(1),rsolverNode%iresNorm)
    dres = lsysbl_vectorNorm (p_rv(1),LINALG_NORMEUCLID)
    if (.not.((dfr .ge. 1E-99_DP) .and. &
              (dfr .le. 1E99_DP))) dfr = 0.0_DP

    ! Initialise starting residuum
    rsolverNode%dinitialDefect = dfr
    rsolverNode%dlastDefect = 0.0_DP

    ! initialise the queue of the last residuals with RES
    Dresqueue = dfr

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    if (rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero) then
     
      ! final defect is 0, as initialised in the output variable above
      call lsysbl_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dfr
    else
      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('GMRES('//&
             trim(sys_siL(idim,10))//'): Iteration '//&
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

      ite = 0
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= Outer Loop
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      do while (ite < rsolverNode%nmaxIterations)
        
        ! Step O.1:
        ! Create elementary vector e_1 scaled by the euclid norm of the
        ! defect, i.e.: q = ( ||v(1)||_2, 0, ..., 0 )
        call storage_clear (hDq)
        p_Dq(1) = dres
        
        ! Step O.2:
        ! Now scale the defect by the inverse of its norm
        call lsysbl_scaleVector (p_rv(1), 1.0_DP / dres)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! -= Inner Loop (GMRES iterations)
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        do i = 1, idim
          ! Another iteration
          ite = ite + 1
        
          ! Step I.1:
          ! Solve P * z(i) = v(i), where P is the preconditioner matrix
          call lsysbl_copyVector (p_rv(i), p_rz(i))
          
          ! Apply preconditioner to z(i)
          if (bprec) then
            call linsol_precondDefect (p_rprecSubnode, p_rz(i))
          end if
          
          ! Apply filter chain to z(i)
          if (bfilter) then
            call filter_applyFilterChainVec (p_rz(i), p_RfilterChain)
          end if
          
          
          ! Step I.2:
          ! Calculate v(i+1) = A * z(i)
          call lsysbl_blockMatVec (p_rmatrix, p_rz(i), p_rv(i+1), 1.0_DP, 0.0_DP)
          
          
          ! Step I.3:
          ! Perfom Gram-Schmidt process (1st time)
          do k = 1, i
            p_Dh(k,i) = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
            call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -p_Dh(k,i), 1.0_DP)
          end do
          
          ! If the user wishes, perform Gram-Schmidt one more time
          ! to improve numerical stability.
          if (btwiceGS) then
            do k = 1, i
              dtmp = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
              p_Dh(k,i) = p_Dh(k,i) + dtmp;
              call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -dtmp, 1.0_DP)
            end do
          end if
          
          ! Step I.4:
          ! Calculate alpha = ||v(i+1)||_2
          dalpha = lsysbl_vectorNorm (p_rv(i+1), LINALG_NORMEUCLID)
          
          ! Step I.5:
          ! Scale v(i+1) by the inverse of its euclid norm
          if (dalpha > SYS_EPSREAL_DP) then
            call lsysbl_scaleVector (p_rv(i+1), 1.0_DP / dalpha)
          else
            ! Well, let us just print a warning here...
            if(rsolverNode%ioutputLevel .ge. 2) then
              call output_line ('GMRES('// trim(sys_siL(idim,10))//&
                '): Warning: !!v(i+1)!! < EPS !')
            end if
          end if
          
          
          ! Step I.6:
          ! Apply Givens rotations to the i-th column of h which
          ! renders the Householder matrix an upper triangular matrix
          do k = 1, i-1
            dtmp = p_Dh(k,i)
            p_Dh(k,i)   = p_Dc(k) * dtmp + p_Ds(k) * p_Dh(k+1,i)
            p_Dh(k+1,i) = p_Ds(k) * dtmp - p_Dc(k) * p_Dh(k+1,i)
          end do
          
          
          ! Step I.7:
          ! Calculate Beta
          ! Beta = (h(i,i)^2 + alpha^2) ^ (1/2)
          dbeta = sqrt(p_Dh(i,i)**2 + dalpha**2)
          if (dbeta < SYS_EPSREAL_DP) then
            dbeta = SYS_EPSREAL_DP
            if(rsolverNode%ioutputLevel .ge. 2) then
              call output_line ('GMRES('// trim(sys_siL(idim,10))//&
                '): Warning: beta < EPS !')
            end if
          end if
          
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
          dpseudores = abs(p_Dq(i+1))
          
          ! Print our current pseudo-residual-norm
          if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
            call output_line ('GMRES('//&
              trim(sys_siL(idim,10))//'): Iteration '//&
              trim(sys_siL(ite,10))//',  !q(i+1)! = '//&
              trim(sys_sdEL(dpseudores,15)) )
          end if

          ! The euclid norm of our current defect is implicitly given by
          ! |q(i+1)|, however, it may be inaccurate, therefore, checking if
          ! |q(i+1)| fulfills the stopping criterion would not be very wise,
          ! as it may result in a lot of unnecessary restarts.
          ! We will therefore check (|q(i+1)| * dprnsf) against the tolerances
          ! instead, that should *hopefully* be okay.
          ! If dprnsf is 0, we will check |q(i+1)| against EPS to avoid NaNs
          ! or Infs during the next GMRES iteration.
          dtmp = dpseudores * dprnsf
          
          ! Is |q(i+1)| smaller than our machine`s exactness?
          ! If yes, then exit the inner GMRES loop.
          if (dpseudores .le. SYS_EPSREAL_DP) then
            exit
          
          ! Check if (|q(i+1)| * dprnsf) is greater than the machine`s
          ! exactness (this can only happen if dprnsf is positive).
          ! If yes, call linsol_testConvergence to test against the
          ! stopping criterions.
          else if (dtmp > SYS_EPSREAL_DP) then
            
            ! If (|q(i+1)| * dprnsf) fulfills the stopping criterion
            ! then exit the inner GMRES loop
            if (linsol_testConvergence(rsolverNode,dtmp)) exit
          
            ! We also want to check if our solution is diverging.
            if (linsol_testDivergence(rsolverNode,dpseudores)) then
              if(rsolverNode%ioutputLevel .ge. 2) then
                call output_line ('GMRES('// trim(sys_siL(idim,10))//&
                  '): Warning: Pseudo-residuals diverging!')
              end if
              
              ! Instead of exiting the subroutine, we just exit the inner loop
              exit
            end if
          
          end if
          
          ! Maybe we have already done enough iterations?
          if (ite .ge. rsolverNode%nmaxIterations) exit
          
          ! Okay, next inner loop iteration, please
        end do
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! -= End of Inner Loop
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        
        ! Step O.3:
        ! Get the dimension of the current basis of the krylov subspace
        i = min(i, idim)
        
        ! Step O.4:
        ! Solve H' * q = q , where H' is the upper triangular matrix of
        ! the upper left (i x i)-block of H.
        ! We use the BLAS Level 2 subroutine 'dtrsv' to do the work for us
        call dtrsv('U', 'N', 'N', i, p_Dh, idim, p_Dq, 1)
        
        ! Step O.5:
        ! Update our solution vector
        ! x = x + q(1)*z(1) + q(2)*z(2) + ... + q(i)*z(i)
        do k = 1, i
          call lsysbl_vectorLinearComb(p_rz(k), p_rx, p_Dq(k), 1.0_DP)
        end do
        
        ! Step O.6:
        ! Calculate 'real' residual
        ! v(1) = b - (A * x)
        call lsysbl_copyVector (rd, p_rv(1))
        call lsysbl_blockMatVec(p_rmatrix, p_rx, p_rv(1), -1.0_DP, 1.0_DP)

        ! Step O.7:
        ! Call filter chain if given.
        if (bfilter) then
          call filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
        end if
        
        ! Step O.8:
        ! Calculate euclid norm of the residual (needed for next q)
        dres = lsysbl_vectorNorm (p_rv(1), LINALG_NORMEUCLID)
        
        ! Calculate residual norm for stopping criterion
        ! TODO: try to avoid calculating the norm twice
        dfr = lsysbl_vectorNorm (p_rv(1), rsolverNode%iresNorm)
        
        ! Step O.9:
        ! Test for convergence, divergence and write some output now
        
        ! Shift the queue with the last residuals and add the new
        ! residual to it.
        Dresqueue = eoshift(Dresqueue,1,dfr)

        rsolverNode%dlastDefect = rsolverNode%dfinalDefect
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        if (linsol_testDivergence(rsolverNode,dfr)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('GMRES: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
     
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (linsol_testConvergence(rsolverNode,dfr)) exit
          
        end if
        
        ! We need to check if the euclidian norm of the
        ! residual is maybe < EPS - if this is true,
        ! then starting another GMRES iteration would
        ! result in NaNs / Infs!!!
        if (dres .le. SYS_EPSREAL_DP) exit

        ! print out the current residuum

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Iteration '//&
            trim(sys_siL(ite,10))//',  !!RES!!  = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if

      end do
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= End of Outer Loop
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      ! Set ite to rsolverNode%nmaxIterations to prevent printing of
      ! "rsolverNode%nmaxIterations+1" of the loop was completed
      if (ite .gt. rsolverNode%nmaxIterations) then
        ! Warning if we did not reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        if ((rsolverNode%ioutputLevel .ge. 0) .and. &
            (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
          call output_line ('GMRES: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        end if

        ite = rsolverNode%nmaxIterations
      end if

      ! Finish - either with an error or if converged.
      ! Print the last residuum.
      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ite .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('GMRES('// trim(sys_siL(idim,10))//&
          '): Iteration '//&
          trim(sys_siL(ite,10))//',  !!RES!!  = '//&
          trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if

    end if

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning. Scale the vector by the damping parameter
    ! in the solver structure.
    call lsysbl_vectorLinearComb (p_rx,rd,rsolverNode%domega,0.0_DP)
      
    ! Do not calculate anything if the final residuum is out of bounds -
    ! would result in NaN`s,...
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
    
      ! Calculate asymptotic convergence rate
      if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
        i = min(rsolverNode%iiterations,niteAsymptoticCVR)
        rsolverNode%dasymptoticConvergenceRate = &
          (rsolverNode%dfinalDefect / &
            dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
      end if

      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :
      if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/real(rsolverNode%iiterations,DP))
      else
        rsolverNode%dconvergenceRate = 0.0_DP
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('GMRES('// trim(sys_siL(idim,10))// ') statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line ('GMRES('// trim(sys_siL(idim,10))//&
          '): Iterations/Rate of convergence: ' //&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if

    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    end if
  
  end subroutine
  
! *****************************************************************************
! Routines for the Schur-Complement solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initSchur (p_rsolverNode,p_rsolverA,p_rsolverS,&
                               RmatrixS,ctype)
  
!<description>
  ! Creates a t_linsolNode solver structure for the Schur solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to the solver structure of a solver that should be used to
  ! solve the A-submatrix system.
  type(t_linsolNode), pointer :: p_rsolverA
  
  ! A pointer to the solver structure of a solver that should be used to
  ! solve the Schur-complement system.
  type(t_linsolNode), pointer :: p_rsolverS
  
  ! An array of Schur-complement matrices.
  type(t_matrixBlock), dimension(:), target, intent(in) :: RmatrixS
  
  ! OPTIONAL: One of the LINSOL_SCHUR_TYPE_XXXX constants which specifies
  ! which type of preconditioner is to be used.
  ! If not given, LINSOL_SCHUR_TYPE_DIAG is used.
  integer, optional, intent(in) :: ctype
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_SCHUR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT
    
    ! Allocate the subnode for Schur.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeSchur)
    
    ! Attach the solvers
    p_rsolverNode%p_rsubnodeSchur%p_rsolverA => p_rsolverA
    p_rsolverNode%p_rsubnodeSchur%p_rsolverS => p_rsolverS
    
    ! Attach the Schur-complement matrix
    p_rsolverNode%p_rsubnodeSchur%p_RmatrixS => RmatrixS
    
    ! Set type if given
    if(present(ctype)) &
      p_rsolverNode%p_rsubnodeSchur%ctype = ctype
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_matCompatSchur (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>
  
  integer :: i,j

    ! Normally, we can handle the matrix. This solver can usually handle
    ! everything.
    ccompatible = LINSOL_COMP_OK
    
    i = ubound(Rmatrices,1)
    j = Rmatrices(i)%nblocksPerRow
    
    ! Let us see whether we can handle this matrix.
    ! It needs to be at least a 2x2 block matrix.
    if((Rmatrices(i)%nblocksPerRow .lt. 2) .or. &
       (Rmatrices(i)%nblocksPerCol .lt. 2)) &
      ccompatible = LINSOL_COMP_ERRGENERAL
    
    ! And it needs to be square, of course.
    if(Rmatrices(i)%nblocksPerRow .ne. Rmatrices(i)%nblocksPerRow) &
      ccompatible = LINSOL_COMP_ERRMATRECT
    
    ! And finally, we need to make sure that the block at the lower right
    ! corner has the same dimensions as our Schur-Complement matrix.
    !if((rsolverNode%p_rsubnodeSchur%p_rmatScS%NEQ .ne. &
    !    Rmatrices(i)%RmatrixBlock(j,j)%NEQ) .or. &
    !   (rsolverNode%p_rsubnodeSchur%p_rmatScS%NCOLS .ne. &
    !    Rmatrices(i)%RmatrixBlock(j,j)%NCOLS)) then
    !  ccompatible = LINSOL_COMP_ERRGENERAL

    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_setMatrixSchur (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!</subroutine>

  integer :: i,k,n

    k = ubound(Rmatrices,1)
    n = Rmatrices(k)%nblocksPerRow

    ! Now comes the tricky part: We need to decompose the input matrix
    ! into sub-matrices.
    
    ! Allocate a new array of A-submatrices and their underlying
    ! discretisations.
    if(associated(rsolverNode%p_rsubnodeSchur%p_RmatrixA)) &
      deallocate(rsolverNode%p_rsubnodeSchur%p_RmatrixA)
    if(associated(rsolverNode%p_rsubnodeSchur%p_RdiscrA)) &
      deallocate(rsolverNode%p_rsubnodeSchur%p_RdiscrA)
    allocate(rsolverNode%p_rsubnodeSchur%p_RmatrixA(1:k))
    allocate(rsolverNode%p_rsubnodeSchur%p_RdiscrA(1:k))
    
    ! Extract the A-submatrices on all levels. We need to do this on all
    ! levels since the user may want to use multigrid as a solver for
    ! the A-submatrix!
    do i = 1, k
      
      ! Derive a sub-discretisation for the A-submatrices.
      call spdiscr_deriveBlockDiscr(Rmatrices(i)%p_rblockDiscrTrial, &
          rsolverNode%p_rsubnodeSchur%p_RdiscrA(i), 1, n-1)
      
      ! Derive a sub-matrix
      call lsysbl_deriveSubmatrix(Rmatrices(i), &
          rsolverNode%p_rsubnodeSchur%p_rmatrixA(i), &
          LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE, 1, n-1, 1, n-1)
      
      ! And attach the sub-discretisation to the sub-matrix
      rsolverNode%p_rsubnodeSchur%p_rmatrixA(i)%p_rblockDiscrTrial => &
        rsolverNode%p_rsubnodeSchur%p_RdiscrA(i)
      rsolverNode%p_rsubnodeSchur%p_rmatrixA(i)%p_rblockDiscrTest => &
        rsolverNode%p_rsubnodeSchur%p_RdiscrA(i)
      
    end do

    ! Extract the B-submatrix on finest level
    call lsysbl_deriveSubmatrix(Rmatrices(k), &
        rsolverNode%p_rsubnodeSchur%rmatrixB, &
        LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE, 1, n-1, n, n)

    ! Extract the D-submatrix on finest level
    call lsysbl_deriveSubmatrix(Rmatrices(k), &
        rsolverNode%p_rsubnodeSchur%rmatrixD, &
        LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE, n, n, 1, n-1)
        
    ! Attach the A-matrices into the corresponding solver
    call linsol_setMatrices(rsolverNode%p_rsubnodeSchur%p_rsolverA, &
                            rsolverNode%p_rsubnodeSchur%p_RmatrixA)

    ! Attach the S-matrices into the corresponding solver
    call linsol_setMatrices(rsolverNode%p_rsubnodeSchur%p_rsolverS, &
                            rsolverNode%p_rsubnodeSchur%p_RmatrixS)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initStructureSchur (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the init routines of the subsolvers
    call linsol_initStructure(rsolverNode%p_rsubnodeSchur%p_rsolverA, isubgroup)
    call linsol_initStructure(rsolverNode%p_rsubnodeSchur%p_rsolverS, isubgroup)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initDataSchur (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the init routines of the subsolvers
    call linsol_initData(rsolverNode%p_rsubnodeSchur%p_rsolverA, isubgroup)
    call linsol_initData(rsolverNode%p_rsubnodeSchur%p_rsolverS, isubgroup)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataSchur (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the done routines of the subsolvers
    call linsol_doneData(rsolverNode%p_rsubnodeSchur%p_rsolverS, isubgroup)
    call linsol_doneData(rsolverNode%p_rsubnodeSchur%p_rsolverA, isubgroup)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureSchur (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the done routines of the subsolvers
    call linsol_doneStructure(rsolverNode%p_rsubnodeSchur%p_rsolverS, isubgroup)
    call linsol_doneStructure(rsolverNode%p_rsubnodeSchur%p_rsolverA, isubgroup)
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneSchur (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it is also released from the heap by calling
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the solver.
  type(t_linsolNode), pointer         :: rsolverNode
!</input>
  
!</subroutine>

  integer :: i

    ! Release memory if still associated
    call linsol_doneDataSchur (rsolverNode, rsolverNode%isolverSubgroup)
    call linsol_doneStructureSchur (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the subsolvers
    call linsol_releaseSolver(rsolverNode%p_rsubnodeSchur%p_rsolverS)
    call linsol_releaseSolver(rsolverNode%p_rsubnodeSchur%p_rsolverA)
    
    !deallocate(rsolverNode%p_rsubnodeSchur%p_rsolverS)
    !deallocate(rsolverNode%p_rsubnodeSchur%p_rsolverA)
    
    ! Release submatrices B and D on finest level
    call lsysbl_releaseMatrix(rsolverNode%p_rsubnodeSchur%rmatrixD)
    call lsysbl_releaseMatrix(rsolverNode%p_rsubnodeSchur%rmatrixB)
    
    ! Release submatrices A and their discretisations on all levels.
    do i = ubound(rsolverNode%p_rsubnodeSchur%p_RmatrixA,1), 1, -1

      ! Release sub-matrix A
      call lsysbl_releaseMatrix(rsolverNode%p_rsubnodeSchur%p_RmatrixA(i))

      ! Release sub-discretisation
      call spdiscr_releaseBlockDiscr(rsolverNode%p_rsubnodeSchur%p_RdiscrA(i))

    end do
    
    ! Release discretisation array
    deallocate(rsolverNode%p_rsubnodeSchur%p_RdiscrA)
    
    ! Release matrix array
    deallocate(rsolverNode%p_rsubnodeSchur%p_RmatrixA)
    
    ! Release the subnode itself
    deallocate(rsolverNode%p_rsubnodeSchur)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_precSchur (rsolverNode,rd)
  
!<description>
  ! Applies Schur-Complement preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! Our structure
  type(t_linsolSubnodeSchur), pointer :: p_rsubnode
  
  ! The system matrix
  type(t_matrixblock), pointer :: p_rmatrix
  
  ! Two sub-vectors
  type(t_vectorBlock) :: ru, rp
  
  ! Type of preconditioner
  integer :: ctype
  
    ! The Schur-complement preconditioner works for saddle-point systems
    ! of the following structure:
    !
    !             ( A  B ) * ( u ) = ( f )
    !             ( D  0 )   ( p )   ( g )
    !
    ! Depending on the parameter ctype, this preconditioner solves:
    !
    ! For ctype = LINSOL_SCHUR_TYPE_DIAG:
    !
    !             ( A  0 ) * ( u ) = ( f )
    !             ( 0 -S )   ( p )   ( g )
    !
    ! For ctype = LINSOL_SCHUR_TYPE_LTRI:
    !
    !             ( A  0 ) * ( u ) = ( f )
    !             ( D -S )   ( p )   ( g )
    !
    ! Where S is an approximation of the Schur-Complement of A:
    !
    !              S \approx (D * A^-1 * B)
    !
    ! For solving the sub-systems with A and S, two linear solvers are used.

    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Get some information
    p_rsubnode => rsolverNode%p_rsubnodeSchur
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if
    
    ! Get the preconditioner type
    ctype = p_rsubnode%ctype
    
    ! Derive the sub-vectors u and p for our sub-systems
    call lsysbl_deriveSubvector(rd,ru,1,rd%nblocks-1,.true.)
    call lsysbl_deriveSubvector(rd,rp,rd%nblocks,rd%nblocks,.true.)
    
    !if((ctype .eq. LINSOL_SCHUR_TYPE_DIAG) .or. &
    !   (ctype .eq. LINSOL_SCHUR_TYPE_LTRI)) then
    
      ! Call solver for submatrix A to solve:
      ! A * u = f
      call linsol_precondDefect(p_rsubnode%p_rsolverA, ru)
    
    !end if
    
    if(ctype .eq. LINSOL_SCHUR_TYPE_LTRI) then
      
      ! Update rhs for Schur-complement:
      ! g := g - D * u
      call lsysbl_blockMatVec (p_rsubnode%rmatrixD, ru, rp, -1.0_DP, 1.0_DP)
      
    end if
    
    ! As we now need to solve '-S*p = g <==> 'S*p = -g`, we will simply
    ! scale the defect by -1.
    call lsysbl_scaleVector(rp, -1.0_DP)
    
    ! Call solver for Schur-complement matrix S to solve:
    ! S * p = g
    call linsol_precondDefect(p_rsubnode%p_rsolverS, rp)
    
    ! Release the sub-vectors
    call lsysbl_releaseVector(rp)
    call lsysbl_releaseVector(ru)

  end subroutine

! *****************************************************************************
! Routines for the Multigrid solver
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_addMultigridLevel (p_rlevelInfo,rsolverNode, &
                    rprojection, p_rpresmoother,p_rpostsmoother,&
                    p_rcoarseGridSolver, iappend)
                    
!<description>
  ! This routine adds a new level to the linked list of levels in the multigrid
  ! solver identified by rsolverNode.
  ! The given coarse-grid solver and smoother-structures are attached
  ! to the level. The routine returns a pointer to the new level-info structure
  ! in p_rlevelInfo to allow the caller to modify the standard settings of
  ! that level if necessary.
  !
  ! It is allowed to use p_rpresmoother=p_rpostsmoother. A value of
  ! NULL() is also permissable, this deactivates the corresponding smoother.
  !
  ! It is allowed to use the same rinterlevelProjection on all levels,
  ! since the structure is level independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway).
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver, where the level
  ! should be added to. Must already be initialised for the multigrid
  ! solver.
  type(t_linsolNode), intent(inout)                   :: rsolverNode
!</inputoutput>
  
!<input>
  ! An interlevel projection structure that configures the projection
  ! between the solutions on a finer and a coarser grid. The structure
  ! must have been initialised with mlprj_initProjection.
  ! If not given, multigrid will create a compatible interlevel projection
  ! structure based on the discretisation of the system matrix that is passed
  ! to the solver in the linsol_setMatrices routine.
  !
  ! Note that this structure is level-independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway), so the same structure can be used
  ! to initialise all levels!
  type(t_interlevelProjectionBlock), intent(in) :: rprojection
  
  ! Optional: A pointer to the solver structure of a solver that should be
  ! used for presmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  type(t_linsolNode), pointer, optional   :: p_rpresmoother
  
  ! Optional: A pointer to the solver structure of a solver that should be
  ! used for postsmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  type(t_linsolNode), pointer, optional   :: p_rpostsmoother

  ! Optional: A pointer to the solver structure of a solver that should be
  ! used for coarse grid solving.
  ! Should only be given for the very first level.
  type(t_linsolNode), pointer, optional   :: p_rcoarseGridSolver
  
  ! Optional: Position where to put the new structure.
  ! YES or not specified: append the structure to the end of the list as new
  ! higher level.
  ! NO: Insert the structure as new lowest level. The caller has to make sure
  ! that the coarse grid solver on the previous lowest level is removed!
  ! (i.e. removed from the p_rlevelInfo node corresponding to the lowest
  ! level!)
  integer, intent(in), optional                  :: iappend
!</input>
  
!<output>
  ! The t_levelInfo structure for the new level that was added to the
  ! multigrid solver. The application may modify the structure or throw the
  ! pointer away, it does not matter. Multigrid will maintain the
  ! structure internally.
  type(t_linsolMGLevelInfo), pointer     :: p_rlevelInfo
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_addMultigridLevel')
      call sys_halt()
    end if
    
    ! Create a new level-info structure
    nullify(p_rlevelInfo)
    allocate(p_rlevelInfo)
    
    ! Attach the sub-solvers
    if(present(p_rcoarseGridSolver)) then
      if (associated(p_rcoarseGridSolver)) then
        p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
      end if
    end if

    if(present(p_rpresmoother)) then
      if (associated(p_rpresmoother)) then
        p_rlevelInfo%p_rpresmoother => p_rpresmoother
      end if
    end if

    if(present(p_rpostsmoother)) then
      if (associated(p_rpostsmoother)) then
        p_rlevelInfo%p_rpostsmoother => p_rpostsmoother
      end if
    end if
    
    ! Attach the projection structure and set bautoProjection to .false. to
    ! indicate that the projection structure was given by the caller.
    p_rlevelInfo%rprojection = rprojection
    p_rlevelInfo%bautoProjection = .false.

    ! Attach the level-info structure to the linked list.
    ! Does the list exist?
    if (.not. associated(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)) then
    
      ! List does not exist. Create a new one.
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
      rsolverNode%p_rsubnodeMultigrid%nlevels = 1
    
    else
    
      ! Do we have to insert it as new head?
      i = YES
      if (present(iappend)) i=iappend
      
      if (i .eq. NO) then
      
        ! New lowest level
        p_rlevelInfo%p_rnextLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead%p_rprevLevel => p_rlevelInfo
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
      
      else
      
        ! New highest level
        p_rlevelInfo%p_rprevLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail%p_rnextLevel => p_rlevelInfo
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
        
      end if
      
      ! Increase the number of existing levels
      rsolverNode%p_rsubnodeMultigrid%nlevels = &
        rsolverNode%p_rsubnodeMultigrid%nlevels + 1
    
    end if

  end subroutine

! *****************************************************************************

!<subroutine>
  
  subroutine linsol_addMultigridLevel2 (p_rlevelInfo, rsolverNode, &
                 p_rpresmoother, p_rpostsmoother, p_rcoarseGridSolver, iappend)
                    
!<description>
  ! This routine adds a new level to the linked list of levels in the multigrid
  ! solver identified by rsolverNode.
  ! The given coarse-grid solver and smoother-structures are attached
  ! to the level. The routine returns a pointer to the new level-info structure
  ! in p_rlevelInfo to allow the caller to modify the standard settings of
  ! that level if necessary.
  !
  ! It is allowed to use p_rpresmoother=p_rpostsmoother. A value of
  ! NULL() is also permissable, this deactivates the corresponding smoother.
  !
  ! In contrast to the previous routine, this routine does not recieve
  ! an interlevel projection structure and will advise multigrid to build
  ! a default projection based on the system matrix.
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver, where the level
  ! should be added to. Must already be initialised for the multigrid
  ! solver.
  type(t_linsolNode),intent(inout)          :: rsolverNode
!</inputoutput>
  
!<input>
  
  ! Optional: A pointer to the solver structure of a solver that should be
  ! used for presmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  type(t_linsolNode), pointer, optional   :: p_rpresmoother
  
  ! Optional: A pointer to the solver structure of a solver that should be
  ! used for postsmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  type(t_linsolNode), pointer, optional   :: p_rpostsmoother

  ! Optional: A pointer to the solver structure of a solver that should be
  ! used for coarse grid solving.
  ! Should only be given for the very first level.
  type(t_linsolNode), pointer, optional   :: p_rcoarseGridSolver
  
  ! Optional: Position where to put the new structure.
  ! YES or not specified: append the structure to the end of the list as new
  ! higher level.
  ! NO: Insert the structure as new lowest level. The caller has to make sure
  ! that the coarse grid solver on the previous lowest level is removed!
  ! (i.e. removed from the p_rlevelInfo node corresponding to the lowest
  ! level!)
  integer, intent(in), optional                  :: iappend
!</input>
  
!<output>
  ! The t_levelInfo structure for the new level that was added to the
  ! multigrid solver. The application may modify the structure or throw the
  ! pointer away, it does not matter. Multigrid will maintain the
  ! structure internally.
  type(t_linsolMGLevelInfo), pointer     :: p_rlevelInfo
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_addMultigridLevel')
      call sys_halt()
    end if
    
    ! Create a new level-info structure
    nullify(p_rlevelInfo)
    allocate(p_rlevelInfo)
    
    ! Make sure multigrid will create an interlevel projection structure.
    p_rlevelInfo%bautoProjection = .true.
    
    ! Attach the sub-solvers
    if(present(p_rcoarseGridSolver)) then
      if (associated(p_rcoarseGridSolver)) then
        p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
      end if
    end if

    if(present(p_rpresmoother)) then
      if (associated(p_rpresmoother)) then
        p_rlevelInfo%p_rpresmoother => p_rpresmoother
      end if
    end if

    if(present(p_rpostsmoother)) then
      if (associated(p_rpostsmoother)) then
        p_rlevelInfo%p_rpostsmoother => p_rpostsmoother
      end if
    end if

    ! Attach the level-info structure to the linked list.
    ! Does the list exist?
    if (.not. associated(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)) then
    
      ! List does not exist. Create a new one.
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
      rsolverNode%p_rsubnodeMultigrid%nlevels = 1
    
    else
    
      ! Do we have to insert it as new head?
      i = YES
      if (present(iappend)) i=iappend
      
      if (i .eq. NO) then
      
        ! New lowest level
        p_rlevelInfo%p_rnextLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead%p_rprevLevel => p_rlevelInfo
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
      
      else
      
        ! New highest level
        p_rlevelInfo%p_rprevLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail%p_rnextLevel => p_rlevelInfo
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
        
      end if
      
      ! Increase the number of existing levels
      rsolverNode%p_rsubnodeMultigrid%nlevels = &
        rsolverNode%p_rsubnodeMultigrid%nlevels + 1
    
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_removeMultigridLevel (rsolverNode,bcoarseLevel)
                    
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
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!<input>
  ! OPTIONAL: Defines which level to delete.
  ! .TRUE. = Delete the coarse level.
  ! .FALSE. or not defined = Delete the finest level.
  logical, intent(in), optional :: bcoarseLevel
!</input>
  
!</subroutine>

  ! local variables
  logical :: bcoarse
  
  ! The t_levelInfo structure to delete
  type(t_linsolMGLevelInfo), pointer     :: p_rlevelInfo
  
    bcoarse = .false.
    if (present(bcoarseLevel)) bcoarse = bcoarseLevel
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_removeMultigridLevel')
      call sys_halt()
    end if
    
    ! Where is the levelInfo-structure to delete?
    if (.not. bcoarse) then
      ! Remove the fine level node - if it exists.
      p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail

      if (associated(p_rlevelInfo)) then
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo%p_rprevLevel
        
        if (rsolverNode%p_rsubnodeMultigrid%nlevels .eq. 1) then
          ! We remove the last node.
          nullify(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)
        end if
        
        if (associated(p_rlevelInfo%p_rprevLevel)) then
          nullify(p_rlevelInfo%p_rprevLevel%p_rnextLevel)
        end if
      end if
    else
      ! Remove the coarse level node - if it exists.
      p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
      if (associated(p_rlevelInfo)) then
        
        rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo%p_rnextLevel
        
        if (rsolverNode%p_rsubnodeMultigrid%nlevels .eq. 1) then
          ! We remove the last node.
          nullify(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail)
        end if

        if (associated(p_rlevelInfo%p_rnextLevel)) then
          nullify(p_rlevelInfo%p_rnextLevel%p_rprevLevel)
        end if
      end if
    end if
    
    if (associated(p_rlevelInfo)) then
      ! Decrease the number of existing levels
      rsolverNode%p_rsubnodeMultigrid%nlevels = &
        rsolverNode%p_rsubnodeMultigrid%nlevels - 1

      ! Release the temporary vectors of this level. But be careful!
      !
      ! If the node is a coarse level node...
      if (.not. associated(p_rlevelInfo%p_rprevLevel)) then
        ! ...do not release RHS/temp vector, as there is none. Instead,
        ! release RHS/solution vector of the next higher level as that
        ! this one gets the new coarse grid.
        if (associated(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)) then
          ! The vector may not exist - if the application has not called
          ! initStructure!
          if (rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead%rrhsVector%NEQ .ne. 0) &
            call lsysbl_releaseVector(rsolverNode%p_rsubnodeMultigrid% &
                                      p_rlevelInfoHead%rrhsVector)
          if (rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead%rtempVector%NEQ .ne. 0) &
            call lsysbl_releaseVector(rsolverNode%p_rsubnodeMultigrid% &
                                      p_rlevelInfoHead%rtempVector)
        end if
        
        ! But there may be a solution vector in p_rlevelInfo that
        ! is to be deleted - if the level is not also the coarse grid
        if (rsolverNode%p_rsubnodeMultigrid%nlevels .ne. 0) then
          if (p_rlevelInfo%rsolutionVector%NEQ .ne. 0) &
            call lsysbl_releaseVector(p_rlevelInfo%rsolutionVector)
        end if

      end if

      ! If the node is a fine level node...
      if (.not. associated(p_rlevelInfo%p_rnextLevel)) then
        ! ...do not release solution vector, as there is none (it shares
        ! its memory with the vector rd to be preconditioned). Instead,
        ! release solution vector of the next lower level as that
        ! this one gets the new fine grid.
        if (associated(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail)) then
          if (rsolverNode%p_rsubnodeMultigrid% &
              p_rlevelInfoTail%rsolutionVector%NEQ .ne. 0) &
            call lsysbl_releaseVector(rsolverNode%p_rsubnodeMultigrid% &
                                      p_rlevelInfoTail%rsolutionVector)
        end if
        
        ! But there may be a RHS/temporary vector in p_rlevelInfo that
        ! is to be deleted - if the level is not also the coarse grid
        if (rsolverNode%p_rsubnodeMultigrid%nlevels .ne. 0) then
          if (p_rlevelInfo%rrhsVector%NEQ .ne. 0) &
            call lsysbl_releaseVector(p_rlevelInfo%rrhsVector)
          if (p_rlevelInfo%rtempVector%NEQ .ne. 0) &
            call lsysbl_releaseVector(p_rlevelInfo%rtempVector)
        end if
      end if
    
      ! Release sub-solvers if associated.
      !
      ! Caution: Pre- and postsmoother may be identical!
      if (associated(p_rlevelInfo%p_rpostSmoother) .and. &
          (.not. associated(p_rlevelInfo%p_rpreSmoother, &
                          p_rlevelInfo%p_rpostSmoother))) then
        call linsol_releaseSolver(p_rlevelInfo%p_rpostsmoother)
      end if

      if (associated(p_rlevelInfo%p_rpresmoother)) then
        call linsol_releaseSolver(p_rlevelInfo%p_rpresmoother)
      end if

      if (associated(p_rlevelInfo%p_rcoarseGridSolver)) then
        call linsol_releaseSolver(p_rlevelInfo%p_rcoarseGridSolver)
      end if
      
      ! Check whether multigrid created the interlevel projection structure.
      if(p_rlevelInfo%bautoProjection) then
        
        ! Yes, so we need to release it now.
        call mlprj_doneProjection(p_rlevelInfo%rprojection)
        
      end if
      
      ! Clean up the associated matrix if there is one.
      ! Of course, the memory of the matrices will not be deallocated
      ! because the matrices belong to the application, not to the solver!
      call lsysbl_releaseMatrix(p_rlevelInfo%rsystemMatrix)

      ! Remove the structure
      deallocate(p_rlevelInfo)
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_cleanMultigridLevels (rsolverNode)
                    
!<description>
  ! This routine removes all level information from the MG solver and releases
  ! all attached solver nodes on all levels (smoothers, coarse grid solvers,
  ! ...).
  ! It can be used to clean up all levels before rebuilding the level
  ! structure by addMultigridLevel.
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver.
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Remove all the levels
    do while (rsolverNode%p_rsubnodeMultigrid%nlevels .gt. 0)
      call linsol_removeMultigridLevel (rsolverNode)
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_getMultigridLevel (rsolverNode,ilevel,p_rlevelInfo)
                    
!<description>
  ! Searches inside of the multigrid structure for level ilevel and returns
  ! a pointer to the corresponding p_rlevelInfo structure
  ! (or NULL() if the level does not exist).
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!<input>
  ! Number of the level to fetch.
  integer, intent(in) :: ilevel
!</input>
  
!<output>
  ! A pointer to the corresponding t_levelInfo structure or NULL()
  ! if the level does not exist.
  type(t_linsolMGLevelInfo), pointer     :: p_rlevelInfo
!</output>
  
!</subroutine>

  ! local variables
  integer :: i

    ! Do we have the level?
    if ((ilevel .lt. 0) .or. &
        (ilevel .gt. rsolverNode%p_rsubnodeMultigrid%nlevels)) then
      nullify(p_rlevelInfo)
      return
    end if
    
    ! Go through the linear list to search for the level
    p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    do i=1,ilevel-1
      p_rlevelInfo => p_rlevelInfo%p_rnextLevel
    end do

  end subroutine
  
  ! ***************************************************************************

!<function>
  
  integer function linsol_getMultigridLevelCount (rsolverNode)
                    
!<description>
  ! Returns the number of levels currently attached to the multigrid solver.
!</description>
  
!<input>
  ! The solver structure of the multigrid solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</input>

!<result>
  ! The number of levels in the MG solver.
!</result>
  
!</function>

    linsol_getMultigridLevelCount = rsolverNode%p_rsubnodeMultigrid%nlevels

  end function
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initMultigrid (p_rsolverNode,Rfilter)
  
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
  
  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not present, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter
  
!</input>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer             :: p_rsolverNode
   
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
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
    allocate(p_rsolverNode%p_rsubnodeMultigrid)
    
    ! Initialise the coarse grid correction structure.1
    call cgcor_init (p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection)

    ! Attach the filter if given.
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubnodeMultigrid%p_RfilterChain => Rfilter
    end if
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_alterMultigrid (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the Multigrid solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout)                     :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel

    ! Pass the command structure to all subsolvers and smoothers
    ! on all levels.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    do while(associated(p_rcurrentLevel))
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_alterSolver (p_rcurrentLevel%p_rpreSmoother, &
                                 ralterConfig)
      end if
      ! Pre- and postsmoother may be identical!
      ! Take care not to alter the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_alterSolver (p_rcurrentLevel%p_rpostSmoother, &
                                 ralterConfig)
      end if
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_alterSolver (p_rcurrentLevel%p_rcoarseGridSolver, &
                                 ralterConfig)
      end if
      
      ! Next level -- if there is one.
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
    end do
    
  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_setMatrixMultigrid (rsolverNode,Rmatrices)
  
!<description>
  
  ! Assigns the system matrix rsystemMatrix to the Multigrid solver on
  ! all levels.
  
!</description>
  
!<input>
  
  ! An array of system matrices on all levels.
  ! Each level in multigrid is initialised separately.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices

!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!</subroutine>

  ! local variables
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel
  integer :: ilevel,nlmin
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_setMatrixMultigrid')
      call sys_halt()
    end if
    
    ! Make sure we have the right amount of matrices
    if (size(Rmatrices) .ne. rsolverNode%p_rsubnodeMultigrid%nlevels) then
      call output_line ('Wrong number of matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_setMatrixMultigrid')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached. If yes, attach the matrix
    ! to it.
    ! For this purpose, call the general initialisation routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    nlmin  = lbound(Rmatrices,1)
    ilevel = lbound(Rmatrices,1)
    do while(associated(p_rcurrentLevel))
    
      ! Check whether we need to create an interlevel projection structure,
      ! but only on the non-coarse levels, of course.
      if(p_rcurrentLevel%bautoProjection .and. (ilevel .ne. nlmin)) then
        
        ! Yes, so create a default projection structure based on the system
        ! matrix on this level.
        call mlprj_initProjectionMat(p_rcurrentLevel%rprojection, &
                                     Rmatrices(ilevel))
        
      end if
    
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_setMatrices (p_rcurrentLevel%p_rpreSmoother, &
                                  Rmatrices(nlmin:ilevel))
      end if
      
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_setMatrices (p_rcurrentLevel%p_rpostSmoother, &
                                  Rmatrices(nlmin:ilevel))
      end if
      
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_setMatrices (p_rcurrentLevel%p_rcoarseGridSolver, &
                                  Rmatrices(nlmin:ilevel))
      end if
      
      ! Write a link to the matrix into the level-structure.
      call lsysbl_duplicateMatrix (Rmatrices(ilevel-nlmin+1), &
          p_rcurrentLevel%rsystemMatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Next level -- if there is one.
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
      ilevel = ilevel + 1
      
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_matCompatMultigrid (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  ! More precisely, ccompatible returns always the status of the highest
  ! level where there is an error -- or LINSOL_COMP_OK if there is no error.
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  ! local variables
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel
  integer                       :: ilevel,nlmin,ccompatLevel
    
    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK
    
    ! Reset all compatibility flags
    if (present(CcompatibleDetail)) CcompatibleDetail (:) = ccompatible

    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_matCompatMultigrid')
      call sys_halt()
    end if
    
    ! Make sure we have the right amount of matrices
    if (size(Rmatrices) .ne. rsolverNode%p_rsubnodeMultigrid%nlevels) then
      call output_line ('Wrong number of matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_matCompatMultigrid')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached. If yes, it must check
    ! the matrices!
    ! For this purpose, call the general check routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    nlmin  = lbound(Rmatrices,1)
    ilevel = lbound(Rmatrices,1)
    
    ! Loop through all levels
    do while(associated(p_rcurrentLevel))
    
      ! Compatibility flag for that level
      ccompatLevel = LINSOL_COMP_OK
    
      ! Check presmoother
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
      
        if (present(CcompatibleDetail)) then
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel,CcompatibleDetail(nlmin:ilevel))
        else
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel)
        end if
      
        if (ccompatLevel .ne. LINSOL_COMP_OK) then
          ccompatible = ccompatLevel
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          ilevel = ilevel + 1
          cycle
        end if
        
      end if
      
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        
        if (present(CcompatibleDetail)) then
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel,CcompatibleDetail(nlmin:ilevel))
        else
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(nlmin:ilevel),ccompatLevel)
        end if
        
        if (ccompatLevel .ne. LINSOL_COMP_OK) then
          ccompatible = ccompatLevel
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          ilevel = ilevel + 1
          cycle
        end if
        
      end if
      
      ! Check coarse grid solver
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        
        if (present(CcompatibleDetail)) then
          call linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(nlmin:ilevel),ccompatLevel,CcompatibleDetail(nlmin:ilevel))
        else
          call linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(nlmin:ilevel),ccompatLevel)
        end if
        
        if (ccompatLevel .ne. LINSOL_COMP_OK) then
          ccompatible = ccompatLevel
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          ilevel = ilevel + 1
          cycle
        end if
        
      end if
      
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
      ilevel = ilevel + 1
      
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initStructureMultigrid (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup
  integer :: imemmax
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_vectorBlock), pointer :: p_rtemplVect
  
  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.

    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initStructureMultigrid')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    do while(associated(p_rcurrentLevel))
    
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_initStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_initStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_initStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
      
    end do
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    imemmax = 0

    ! Multigrid needs temporary vectors for the RHS, solution and general
    ! data. Memory for that is allocated on all levels for temporary, RHS and
    ! solution vectors.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    do while(associated(p_rcurrentLevel))
    
      p_rmatrix => p_rcurrentLevel%rsystemMatrix
      
      nullify(p_rtemplVect)
      
      if (associated(p_rcurrentLevel%p_rnextLevel)) then
      
        ! On the maximum level we do not need additional memory for the
        ! solution vector - as solution vector on the fine grid,
        ! the vector which comes as parameter
        ! to the multigrid  preconditioner is used.
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
              p_rcurrentLevel%rsolutionVector,.false.,.false.,&
              rsolverNode%cdefaultDataType)
              
        p_rtemplVect => p_rcurrentLevel%rsolutionVector
        
      end if
      ! On the coarse grid, we do not need memory for temporary/RHS
      ! vectors. The RHS on the coarse grid is replaced in-situ
      ! by the solution vector.
      if (associated(p_rcurrentLevel%p_rprevLevel)) then
      
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
              p_rcurrentLevel%rrhsVector,.false.,.false.,&
              rsolverNode%cdefaultDataType)
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
              p_rcurrentLevel%rtempVector,.false.,.false.,&
              rsolverNode%cdefaultDataType)
              
        p_rtemplVect => p_rcurrentLevel%rtempVector
        
      end if

      ! Calculate the memory that is probably necessary for resorting
      ! vectors during the iteration.
      if (associated(p_rtemplVect)) then
        imemmax = max(imemmax,p_rtemplVect%NEQ)
      end if
         
      ! Ask the prolongation/restriction routines how much memory is necessary
      ! for prolongation/restriction on that level.
      if (associated(p_rcurrentLevel%p_rprevLevel)) then ! otherwise coarse grid
        ! Calculate the memory that is necessary for prolongation/restriction -
        ! in case there is temporary memory needed.
        ! The system matrix on the fine/coarse grid specifies the discretisation.
        imemmax = max(imemmax,mlprj_getTempMemoryMat ( &
                      p_rcurrentLevel%rprojection, &
                      p_rcurrentLevel%p_rprevLevel%rsystemMatrix,&
                      p_rcurrentLevel%rsystemMatrix))
      end if

      ! All temporary vectors are marked as 'unsorted'. We can set the
      ! sorting flag directly here without using the resorting routines
      ! as the vectors have just been created.
      !
      ! Do not do anything to the RHS/temp vectors on the coarse grid
      ! and the solution vector on the fine grid -- they do not exist!
      if (associated(p_rcurrentLevel%p_rprevLevel)) then
        p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy)
        p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy)
      end if
      
      if (associated(p_rcurrentLevel%p_rnextLevel)) then
        p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy)
      end if
      
      ! And the next level...
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
    end do
    
    ! Do we have to allocate memory for prolongation/restriction/resorting?
    if (imemmax .gt. 0) then
    
      call lsyssc_createVector (rsolverNode%p_rsubnodeMultigrid%rprjTempVector,&
                                imemmax,.false.,rsolverNode%cdefaultDataType)
                                
      ! Use this vector also for the coarse grid correction.
      ! Create a block vector on every level with the structure of the
      ! RHS that shares the memory with rprjTempVector.
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
      do while(associated(p_rcurrentLevel))
      
        if (associated(p_rcurrentLevel%p_rprevLevel)) then
          call lsysbl_createVecFromScalar (&
              rsolverNode%p_rsubnodeMultigrid%rprjTempVector,&
              p_rcurrentLevel%rcgcorrTempVector)
        
          call lsysbl_enforceStructure (&
              p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rcgcorrTempVector)
        end if
      
        p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
      end do
    else
      rsolverNode%p_rsubnodeMultigrid%rprjTempVector%NEQ = 0
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initDataMultigrid (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup
  
  ! A-priori we have no error...
  ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup
    
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataMultigrid')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
    do while(associated(p_rcurrentLevel))
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_initData(p_rcurrentLevel%p_rpreSmoother,isubgroup,ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
      ! Pre- and postsmoother may be identical!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_initData(p_rcurrentLevel%p_rpostSmoother,isubgroup,ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_initData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup,ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataMultigrid (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_doneDataMultigrid')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
    do while(associated(p_rcurrentLevel))
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_doneData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
      end if
      ! Pre- and postsmoother may be identical!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_doneData(p_rcurrentLevel%p_rpostSmoother,isubgroup)
      end if
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_doneData(p_rcurrentLevel%p_rpreSmoother,isubgroup)
      end if
      p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureMultigrid (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_doneStructureMultigrid')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
    do while(associated(p_rcurrentLevel))
    
      ! Release the coarse grid solver
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_doneStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
      end if
      
      ! Pre- and postsmoother may be identical! If not, first release the
      ! postsmoother.
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_doneStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
      end if
      
      ! Release the presmoother
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_doneStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
      end if
      
      p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
    end do

    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Release temporary memory.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
    do while(associated(p_rcurrentLevel))
    
      if (associated(p_rcurrentLevel%p_rprevLevel)) then
      
        if (p_rcurrentLevel%rtempVector%NEQ .ne. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rtempVector)
        end if
        
        if (p_rcurrentLevel%rrhsVector%NEQ .ne. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rrhsVector)
        end if

        ! Release the temp vector for the coarse grid correction.
        ! The vector shares its memory with the 'global' projection
        ! vector, so only release it if the other vector exists.
        if (rsolverNode%p_rsubnodeMultigrid%rprjTempVector%NEQ .gt. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rcgcorrTempVector)
        end if

      end if

      if (associated(p_rcurrentLevel%p_rnextLevel)) then
        if (p_rcurrentLevel%rsolutionVector%NEQ .ne. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rsolutionVector)
        end if
      end if

      p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
    end do

    ! Check if we have to release our temporary vector for prolongation/
    ! restriction.
    ! Do we have to allocate memory for prolongation/restriction?
    if (rsolverNode%p_rsubnodeMultigrid%rprjTempVector%NEQ .gt. 0) then
      call lsyssc_releaseVector (rsolverNode%p_rsubnodeMultigrid%rprjTempVector)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneMultigrid (rsolverNode)
  
!<description>
  
  ! This routine releases all temporary memory for the multigrid solver from
  ! the heap. In particular, this releases the solver structures of all
  ! subsolvers (smoother, coarse grid solver).
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
  
!</description>
  
!<inputoutput>
  ! A pointer to a t_linsolNode structure of the Multigrid solver.
  type(t_linsolNode), intent(inout)                :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! local variables
  !TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_doneMultigrid')
      call sys_halt()
    end if

!  ! Check for every level, if there is a presmoother, postsmoother or
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
    call linsol_doneDataMultigrid(rsolverNode)
    call linsol_doneStructureMultigrid(rsolverNode)
    
    ! Remove all the levels
    call linsol_cleanMultigridLevels (rsolverNode)
    
    ! Release the coarse grid correction structure
    call cgcor_release (rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection)
    
    ! Release MG substructure, that is it.
    deallocate(rsolverNode%p_rsubnodeMultigrid)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_smoothCorrection (rsolverNode,&
                       rmatrix,rx,rb,rtemp,p_RfilterChain)
  
!<description>
  ! This routine performs a smoothing process on the vector rx
  ! belonging to a linear system <tex>$ Ax=b $</tex>.
  ! rsolverNode identifies a solver structure that is converted to a
  ! smoother using linsol_convertToSmoother: The number of smoothing
  ! steps is written to rsolverNode%nmaxIterations and
  ! rsolverNode%nmaxIterations so that the solver performs a definite
  ! number of iterations regardless of the residual.
  !
  ! rx is assumed to be of 'defect' type. There are two types of smoothing
  ! processes, depending on which type of preconditioner rsolverNode is:
  ! 1.) If rsolverNode is an iterative solver, linsol_smoothCorrection
  !     calls the associated solver P to compute <tex>$x=P^{-1}b$</tex> using
  !     rsolverNode%nmaxIterations iterations.
  ! 2.) If rsolverNode is a 1-step solver, linsol_smoothCorrection
  !     performs rsolverNode%nmaxIterations defect correction steps
  !     <tex> $$ x := x + P^{-1} (b-Ax) $$
  !     with $x_0=0$ </tex>.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<input>
  ! The RHS vector of the system
  type(t_vectorBlock), intent(in), target            :: rb

  ! The system matrix on the current level.
  type(t_matrixBlock), intent(in)                    :: rmatrix
  
  ! Pointer to the filter chain to use if rsolverNode is a direct
  ! solver. NULL() if no filtering is active.
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
!</input>
  
!<inputoutput>
  ! The solver node containing the solver configuration
  type(t_linsolNode), intent(inout)                  :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  type(t_vectorBlock), intent(inout)                 :: rx
  
  ! A temporary vector of the same size and structure as rx.
  type(t_vectorBlock), intent(inout)                 :: rtemp
!</inputoutput>
  
!</subroutine>

  integer :: i
  logical :: bfilter
  integer :: iiterations
  real(DP) :: dres
  !DEBUG: REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
    
    ! Cancel if nmaxIterations = number of smoothing steps is =0.
    if (rsolverNode%nmaxIterations .le. 0) return
    
    ! Some solvers can be used without the usual preconditioner approach to gain speed.
    ! First check if we have such a solver; these will not need the additional matrix
    ! vector multiplication below!
    select case (rsolverNode%calgorithm)
      ! We are in a case where we can apply an adaptive 1-step defect-correction
      ! type preconditioner as smoother. This is a very special case and applies
      ! only to some special algorithms, which work directly with a solution-
      ! and RHS vector, like special VANKA variants.
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
      ! Or the algorithm filters internally, so we do not have to take into
      ! account the filtering at all.
    
    case (LINSOL_ALG_VANKA)
      ! The basic proceeding is like below: Apply the corresponding solver multiple
      ! times to the solution- and RHS-vector to improve the solution. Let us see
      ! which solver we have.
      ! If the previous IF-guess was wrong, we will leave this IF-clause and
      ! fall back to use the preconditioner approach.
      
      !DEBUG: CALL lsysbl_getbase_double (rx,p_Ddata)
      !DEBUG: CALL lsysbl_getbase_double (rb,p_Ddata2)
      
      ! So, let us see if we have a special implementation here...
      select case (rsolverNode%p_rsubnodeVANKA%csubtypeVANKA)
      case (LINSOL_VANKA_GENERALDIRECT,&
            LINSOL_VANKA_2DNAVSTDIRECT,&
            LINSOL_VANKA_2DNAVSTDIRECTSB,&
            LINSOL_VANKA_2DFNAVSTDIRECT,&
            LINSOL_VANKA_2DFNAVSTOCDIRECT,&
            LINSOL_VANKA_2DFNAVSTOCDIAGDIR,&
            LINSOL_VANKA_2DFNAVSTOCDIAGDIR2,&
            LINSOL_VANKA_3DNAVSTDIRECT,&
            LINSOL_VANKA_3DFNAVSTDIRECT,&
            ! --------------- NEW IMPLEMENTATION ---------------
            LINSOL_VANKA_NAVST2D_DIAG,&
            LINSOL_VANKA_NAVST2D_FULL,&
            LINSOL_VANKA_BOUSS2D_DIAG,&
            LINSOL_VANKA_BOUSS2D_FULL)
            
        ! Yes, this solver can be applied to a given solution/rhs vector directly.
        ! Call it nmaxIterations times to perform the smoothing.
        do i = 1, rsolverNode%nmaxIterations
        
          ! Probably print the residuum
          if (rsolverNode%ioutputLevel .ge. 2) then
            call lsysbl_copyVector(rb,rtemp)
            call lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
            
            dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
            if (.not.((dres .ge. 1E-99_DP) .and. (dres .le. 1E99_DP))) dres = 0.0_DP
                      
            call output_line ('Smoother: Step '//trim(sys_siL(i-1,10))//&
                ' !!RES!! = '//trim(sys_sdEL(dres,15)) )
          end if
        
          ! Perform nmaxIterations:   x_n+1 = x_n + C^{-1} (b-Ax_n)
          ! without explicitly calculating (b-Ax_n) like below.
          call vanka_conformal (rsolverNode%p_rsubnodeVANKA%rvanka, &
                                rx, rb, rsolverNode%domega)
                                
        end do

        ! Probably print the residuum
        if (rsolverNode%ioutputLevel .ge. 2) then
          call lsysbl_copyVector(rb,rtemp)
          call lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
          
          dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
          if (.not.((dres .ge. 1E-99_DP) .and. (dres .le. 1E99_DP))) dres = 0.0_DP
                    
          call output_line ('Smoother: Step '//trim(sys_siL(i-1,10))//&
              ' !!RES!! = '//trim(sys_sdEL(dres,15)) )
        end if

        ! That is it.
        return
      
      end select
      
    case (LINSOL_ALG_SPSOR)
      ! SP-SOR smoother. We can use a fast implementation for this one, however,
      ! we will only use it if we do not have to print the residuals.
      if(rsolverNode%ioutputLevel .lt. 2) then
      
        ! Now make sure that we do NOT use a SP-SSOR (i.e. symmetric SP-SOR)
        ! variant, as this would result in a fatal error...
        select case(rsolverNode%p_rsubNodeSPSOR%csubtype)
        case (LINSOL_SPSOR_NAVST2D, &
              LINSOL_SPSOR_NAVST2D_DIAG)
              
           ! Okay, let the SP-SOR solver perform nmaxIterations
           call spsor_solve (rsolverNode%p_rsubnodeSPSOR%rdata, rx, rb, &
                             rsolverNode%nmaxIterations, rsolverNode%domega)
           
           ! And jump out of this routine
           return
        
        end select
      
      end if

    end select
      
    ! This is a 1-step solver, we have to emulate the smoothing
    ! iterations. Perform rsolverNode%nmaxIterations steps of the
    ! form
    !     <tex>$$ x_{n+1} = x_n + P^{-1}(b-Ax_n) $$</tex>
    ! with $x_0 = 0$.
    
    bfilter = associated(p_RfilterChain)
    
    !DEBUG: CALL lsysbl_getbase_double (rx,p_Ddata)
    !DEBUG: CALL lsysbl_getbase_double (rtemp,p_Ddata2)
    
    ! Do we have an iterative or one-step solver given?
    ! A 1-step solver performs the following loop nmaxIterations times, while an iterative
    ! solver is called only once and performs nmaxIterations steps internally.
    ! (This is a convention. Calling an iterative solver i times with j internal steps
    ! would also be possible, but we do not implement that here.)
    if (iand(rsolverNode%ccapability,LINSOL_ABIL_DIRECT) .ne. 0) then
      iiterations = rsolverNode%nmaxIterations
    else
      iiterations = 1
    end if
    
    do i=1,iiterations
    
      call lsysbl_copyVector(rb,rtemp)
      call lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
      
      if (rsolverNode%ioutputLevel .ge. 2) then
        dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
        if (.not.((dres .ge. 1E-99_DP) .and. (dres .le. 1E99_DP))) dres = 0.0_DP
                  
        call output_line ('Smoother: Step '//trim(sys_siL(i-1,10))//&
            ' !!RES!! = '//trim(sys_sdEL(dres,15)) )
      end if
      
      ! Apply the filter to this defect before preconditioning
      if (bfilter) then
        call filter_applyFilterChainVec (rtemp, p_RfilterChain)
      end if
      
      ! Perform preconditioning
      call linsol_precondDefect(rsolverNode,rtemp)
      
      ! If the preconditioner broke down, cancel the smoothing,
      ! it would destroy out solution!
      if (rsolverNode%iresult .ne. 0) then
        if (rsolverNode%ioutputLevel .ge. 0) then
          call output_line ('Smoothing canceled, preconditioner broke down!', &
                            OU_CLASS_WARNING,OU_MODE_STD,'linsol_smoothCorrection')
          exit
        end if
      end if
      
      call lsysbl_vectorLinearComb (rtemp,rx,1.0_DP,1.0_DP)
      
    end do

    ! Probably print the final residuum
    if (rsolverNode%ioutputLevel .ge. 2) then
      call lsysbl_copyVector(rb,rtemp)
      call lsysbl_blockMatVec (rmatrix, rx, rtemp, -1.0_DP, 1.0_DP)
      
      dres = lsysbl_vectorNorm (rtemp,rsolverNode%iresNorm)
      if (.not.((dres .ge. 1E-99_DP) .and. (dres .le. 1E99_DP))) dres = 0.0_DP
                
      call output_line ('Smoother: Step '//trim(sys_siL(i-1,10))//&
          ' !!RES!! = '//trim(sys_sdEL(dres,15)) )
    end if
    
    ! Apply the filter to the solution to enforce e.g. integral-mean-value = 0
    ! (property may be lost for some VANKA smoothers!)
    if (bfilter) then
      call filter_applyFilterChainVec (rtemp, p_RfilterChain)
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_precMultigrid (rsolverNode,rd)
  
!<description>
  ! Applies Multigrid preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The structure of the levels (pre/postsmoothers, interlevel projection
  ! structures) must have been prepared with linsol_addMultigridLevel
  ! before calling this routine.
  ! The matrices <tex>$ A $</tex> on all levels must be attached to the solver previously
  ! by linsol_setMatrices.
  
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: nminIterations,nmaxIterations,niteResOutput,niteAsymptoticCVR
  integer :: ite,ilev,nlmax,i,j,nblocks
  real(DP) :: dres,dstep
  logical :: bfilter,bsort
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
  ! The queue saves the current residual and the two previous residuals.
  real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
  
  ! The system matrix on the current level
  type(t_matrixBlock), pointer :: p_rmatrix
  
  ! Our MG structure
  type(t_linsolSubnodeMultigrid), pointer :: p_rsubnode
  
  ! The current level and the next lower one.
  type(t_linsolMGLevelInfo), pointer :: p_rcurrentLevel,p_rlowerLevel
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeMultigrid
    
    bfilter = associated(rsolverNode%p_rsubnodeMultigrid%p_RfilterChain)
    p_RfilterChain => rsolverNode%p_rsubnodeMultigrid%p_RfilterChain
    
    ! Get the system matrix on the finest level:
    p_rmatrix => rsolverNode%rsystemMatrix

    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    if (p_rsubnode%icycle .lt. 0) then
      ! Wrong cycle
      call output_line ('Invalid cycle!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precMultigrid')
      rsolverNode%iresult = 2
      return
    end if
    
    if ((.not. associated(p_rsubnode%p_rlevelInfoHead)) .or. &
        (.not. associated(p_rsubnode%p_rlevelInfoTail)) .or. &
        (p_rsubnode%nlevels .le. 0)) then
      call output_line ('No levels attached!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precMultigrid')
      rsolverNode%iresult = 2
      return
    end if

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))

    ! Minimum/maximum number of iterations
    nminIterations = max(rsolverNode%nminIterations,0)
    nmaxIterations = max(rsolverNode%nmaxIterations,0)
      
    ! Iteration when the residuum is printed:
    niteResOutput = max(1,rsolverNode%niteResOutput)

    ! Status reset
    rsolverNode%iresult = 0
    rsolverNode%icurrentIteration = 0
    rsolverNode%dinitialDefect = 0.0_DP
    rsolverNode%dlastDefect = 0.0_DP
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
    if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
    
      if (rsolverNode%ioutputLevel .gt. 1) then
        call output_line ('Multigrid: Only one level. '//&
             'Switching back to standard solver.')
      end if
      call linsol_precondDefect(p_rcurrentLevel%p_rcoarseGridSolver,rd)
      
      ! Take the statistics from the coarse grid solver.
      rsolverNode%dinitialDefect = p_rcurrentLevel%p_rcoarseGridSolver%dinitialDefect
      rsolverNode%dlastDefect = p_rcurrentLevel%p_rcoarseGridSolver%dlastDefect
      rsolverNode%dfinalDefect = p_rcurrentLevel%p_rcoarseGridSolver%dfinalDefect
      rsolverNode%dconvergenceRate = p_rcurrentLevel%p_rcoarseGridSolver%dconvergenceRate
      rsolverNode%iiterations = p_rcurrentLevel%p_rcoarseGridSolver%iiterations
      rsolverNode%dasymptoticConvergenceRate = p_rcurrentLevel% &
        p_rcoarseGridSolver%dasymptoticConvergenceRate
      
    else
      ! Get the norm of the initial residuum.
      ! As the initial iteration vector is zero, this is the norm
      ! of the RHS:
      dres = lsysbl_vectorNorm (rd,rsolverNode%iresNorm)
      if (.not.((dres .ge. 1E-99_DP) .and. &
                (dres .le. 1E99_DP))) dres = 0.0_DP

      ! Initialise starting residuum
      rsolverNode%dinitialDefect = dres
      rsolverNode%dlastDefect = 0.0_DP

      ! initialise the queue of the last residuals with RES
      Dresqueue = dres

      ! Check if out initial defect is zero. This may happen if the filtering
      ! routine filters "everything out"!
      ! In that case we can directly stop our computation.
      if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
      
        ! final defect is 0, as initialised in the output variable above
        call lsysbl_clearVector(rd)
        rsolverNode%dfinalDefect = dres
        rsolverNode%dconvergenceRate = 0.0_DP
        rsolverNode%iiterations = 0
        rsolverNode%dasymptoticConvergenceRate = 0.0_DP
        
      else
        
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
        do while(associated(p_rcurrentLevel%p_rnextLevel))
          if (p_rsubnode%icycle .eq. 0) then
            p_rcurrentLevel%ncycles = 2
          else
            p_rcurrentLevel%ncycles = p_rsubnode%icycle
          end if
          p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
        end do
        
        ! Afterwards, p_rcurrentLevel points to the maximum level again.
        ! There, we set ncycles to 1.
        p_rcurrentLevel%ncycles = 1
        
        ! Print out the initial residuum

        if (rsolverNode%ioutputLevel .ge. 2) then
          call output_line ('Multigrid: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        end if
        
        ! Copy the initial RHS to the RHS vector on the maximum level.
        call lsysbl_copyVector(rd,p_rcurrentLevel%rrhsVector)
        
        ! Replace the solution vector on the finest level by rd.
        ! Afterwards, rd and the solution vector on the finest level
        ! share the same memory location, so we use rd as iteration
        ! vector directly.
        ! We set bisCopy to true to indicate that this vector is
        ! actually does not belong to us. This is just to be sure that
        ! an accidently performed releaseVector does not release a handle
        ! (although this should not happen).
        p_rcurrentLevel%rsolutionVector = rd
        p_rcurrentLevel%rsolutionVector%bisCopy = .true.
        
        ! Clear the initial solution vector.
        call lsysbl_clearVector (p_rcurrentLevel%rsolutionVector)
        
        ! Start multigrid iteration; perform at most nmaxiterations iterations.
        do ite = 1, nmaxiterations
        
          rsolverNode%icurrentIteration = ite
          
          ! Initialise cycle counters for all levels.
          p_rcurrentLevel => p_rsubnode%p_rlevelInfoHead
          do while(associated(p_rcurrentLevel%p_rnextLevel))
            p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles
            p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
          end do
          ! Do not forget the highest level.
          p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles
        
          ! p_rcurrentLevel now points to the maximum level.
          ilev = p_rsubnode%nlevels
          p_rlowerLevel => p_rcurrentLevel%p_rprevLevel
          
          ! Build the defect...
          call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rtempVector)
          if (ite .ne. 1) then   ! initial solution vector is zero!
            call lsysbl_blockMatVec (&
                 p_rcurrentLevel%rsystemMatrix, &
                 p_rcurrentLevel%rsolutionVector,&
                 p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
          end if
          if (bfilter) then
            ! Apply the filter chain to the vector
            call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                             p_RfilterChain)
          end if
          
          cycleloop: do  ! Loop for the cycles
          
            ! On the maximum level we already built out defect vector. If we are
            ! on a lower level than NLMAX, perform smoothing+restriction down to the
            ! coarse level. We identify the coarse level by checking if
            ! the current level has a coarse grid solver.
            
            do while (.not. associated(p_rcurrentLevel%p_rcoarseGridSolver))
            
              ! Perform the pre-smoothing with the current solution vector
              if (associated(p_rcurrentLevel%p_rpreSmoother)) then
                call linsol_smoothCorrection (p_rcurrentLevel%p_rpreSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)
              end if
            
              ! Build the defect vector
              call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                      p_rcurrentLevel%rtempVector)
              call lsysbl_blockMatVec (&
                  p_rcurrentLevel%rsystemMatrix, &
                  p_rcurrentLevel%rsolutionVector,&
                  p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
              if (bfilter) then
                ! Apply the filter chain to the vector
                call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                                p_RfilterChain)
              end if
              
              ! Extended output
              if (associated(p_rcurrentLevel%p_rpreSmoother) .and. &
                  (rsolverNode%ioutputLevel .ge. 3) .and. &
                  (mod(ite,niteResOutput).eq.0)) then
                  
                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                if (.not.((dres .ge. 1E-99_DP) .and. &
                          (dres .le. 1E99_DP))) dres = 0.0_DP
                          
                call output_line ('Multigrid: Level '//trim(sys_siL(ilev,5))//&
                    ' after presm.:     !!RES!! = '//trim(sys_sdEL(dres,15)) )
              end if

              ! Restriction of the defect. The restricted defect is placed
              ! in the right hand side vector of the lower level.
              ! The projection parameters from level ilev to level
              ! ilev-1 is configured in the rprojection structure of the
              ! current level.
              !
              ! Make sure the vector is unsorted before the restriction!
              ! We can use the 'global' temporary array p_rprjTempVector
              ! as temporary memory during the resorting process.
              if (bsort) then
                call lsysbl_sortVectorInSitu (p_rcurrentLevel%rtempVector,&
                      p_rsubnode%rprjTempVector,.false.)
              end if
              
              ! When restricting to the coarse grid, we directy restrict into
              ! the solution vector. It is used there as RHS and replaced in-situ
              ! by the solution by the coarse grid solver. So do not need the RHS vector
              ! on the coarse grid and save one vector-copy.
              !
              ! Otherwise, we restrict to the RHS on the lower level and continue
              ! the smoothing process there.
              if (associated(p_rlowerLevel%p_rprevLevel)) then
                ! We do not project to the coarse grid
                call mlprj_performRestriction (p_rcurrentLevel%rprojection,&
                      p_rlowerLevel%rrhsVector, &
                      p_rcurrentLevel%rtempVector, &
                      p_rsubnode%rprjTempVector)

                ! Apply the filter chain (e.g. to implement boundary conditions)
                ! on just calculated right hand side
                ! (which is a defect vector against the 0-vector).
                if (bfilter) then
                  ! Apply the filter chain to the vector.
                  ! We are in 'unsorted' state; applying the filter here is
                  ! supposed to be a bit faster...
                  call filter_applyFilterChainVec (p_rlowerLevel%rrhsVector, &
                                                   p_RfilterChain)
                end if

                if (bsort) then
                  ! Resort the RHS on the lower level.
                  call lsysbl_sortVectorInSitu (p_rlowerLevel%rrhsVector,&
                        p_rsubnode%rprjTempVector,.true.)

                  ! Temp-vector and solution vector there are yet uninitialised,
                  ! therefore we can mark them as sorted without calling the
                  ! resorting routine.
                  p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                  abs(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy)
                  
                  p_rlowerLevel%rsolutionVector%RvectorBlock(1:nblocks)% &
                    isortStrategy = abs(p_rlowerLevel%rsolutionVector% &
                    RvectorBlock(1:nblocks)%isortStrategy)
                end if

                ! Choose zero as initial vector on lower level.
                call lsysbl_clearVector (p_rlowerLevel%rsolutionVector)
                
                ! Extended output and/or adaptive cycles
                if ((rsolverNode%p_rsubnodeMultigrid%depsRelCycle .ne. 1E99_DP) .or.&
                    (rsolverNode%ioutputLevel .ge. 3)) then
                
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rrhsVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  ! In case adaptive cycles are activated, save the 'initial' residual
                  ! of that level into the level structure. Then we can later check
                  ! if we have to repeat the cycle on the coarse mesh.
                  if (rsolverNode%p_rsubnodeMultigrid%depsRelCycle .ne. 1E99_DP) then
                    p_rlowerLevel%dinitResCycle = dres
                    p_rlowerLevel%icycleCount = 1
                  end if
                         
                  ! If the output level is high enough, print that residuum norm.
                  if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                      (mod(ite,niteResOutput).eq.0)) then
                    call output_line ('Multigrid: Level '//trim(sys_siL(ilev-1,5))//&
                        ' after restrict.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
                  end if
                  
                end if

              else
              
                ! The vector is to be restricted to the coarse grid.
                call mlprj_performRestriction (p_rcurrentLevel%rprojection,&
                      p_rlowerLevel%rsolutionVector, &
                      p_rcurrentLevel%rtempVector, &
                      p_rsubnode%rprjTempVector)

                ! Apply the filter chain (e.g. to implement boundary conditions)
                ! on just calculated right hand side
                ! (which is a defect vector against the 0-vector).
                if (bfilter) then
                  ! Apply the filter chain to the vector.
                  ! We are in 'unsorted' state; applying the filter here is
                  ! supposed to be a bit faster...
                  call filter_applyFilterChainVec (p_rlowerLevel%rsolutionVector, &
                                                   p_RfilterChain)
                end if

                if (bsort) then
                  ! Resort the RHS on the lower level.
                  call lsysbl_sortVectorInSitu (p_rlowerLevel%rsolutionVector,&
                        p_rsubnode%rprjTempVector,.true.)

                  ! Temp-vector and RHS can be ignored on the coarse grid.
                end if

                ! Extended output
                if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                    (mod(ite,niteResOutput).eq.0)) then
                    
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rsolutionVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  call output_line ('Multigrid: Level '//trim(sys_siL(ilev-1,5))//&
                      ' after restrict.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
                end if

              end if
            
              ! Go down one level
              ilev = ilev - 1
              p_rcurrentLevel => p_rcurrentLevel%p_rprevLevel
              p_rlowerLevel => p_rcurrentLevel%p_rprevLevel
              
              ! If we are not on the lowest level, repeat the smoothing of
              ! the solution/restriction of the new defect in the next loop
              ! pass...
            end do   ! ilev > minimum level
            
            ! Now we reached the coarse grid.
            ! Apply the filter chain (e.g. to implement boundary conditions)
            ! on the just calculated right hand side,
            ! which is currently located in rsolutionVector.
            if (bfilter) then
              ! Apply the filter chain to the vector
              call filter_applyFilterChainVec (p_rcurrentLevel%rsolutionVector, &
                                               p_RfilterChain)
            end if
            
            ! Solve the system on lowest level by preconditioning
            ! of the RHS=defect vector.
            call linsol_precondDefect(p_rcurrentLevel%p_rcoarseGridSolver,&
                                      p_rcurrentLevel%rsolutionVector)
            
            ! Now we have the solution vector on the lowest level - we have to go
            ! upwards now... but probably not to NLMAX! That depends on the cycle.
            !
            do while(associated(p_rcurrentLevel%p_rnextLevel))
              ilev = ilev + 1
              p_rlowerLevel => p_rcurrentLevel
              p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
              
              ! Prolongate the solution vector from the coarser level
              ! to the temp vector on the finer level.
              !
              ! Make sure the vector is unsorted before the prolongation!
              ! We can use the 'global' temporary array p_rprjTempVector
              ! as temporary memory during the resorting process.
              if (bsort) then
                call lsysbl_sortVectorInSitu (p_rlowerLevel%rsolutionVector,&
                      p_rsubnode%rprjTempVector,.false.)
                
                ! Temp-vector and RHS vector there unused now,
                ! therefore we can mark them as not sorted without calling the
                ! resorting routine.
                ! This prepares these vectors for the next sweep, when an unsorted
                ! vector comes 'from above'.
                !
                ! Note that this shall not be done on the coarse grid as there is
                ! no temp/rhs vector!
                if (associated(p_rlowerLevel%p_rprevLevel)) then
                  p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                    -abs(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)% &
                    isortStrategy)
                  p_rlowerLevel%rrhsVector%RvectorBlock(1:nblocks)% &
                    isortStrategy = -abs(p_rlowerLevel%rrhsVector% &
                    RvectorBlock(1:nblocks)%isortStrategy)
                end if
              end if
              call mlprj_performProlongation (p_rcurrentLevel%rprojection,&
                    p_rlowerLevel%rsolutionVector, &
                    p_rcurrentLevel%rtempVector, &
                    p_rsubnode%rprjTempVector)

              if (bfilter) then
                ! Apply the filter chain to the vector.
                call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                                 p_RfilterChain)
              end if

              if (bsort) then
                ! Resort the temp vector on the current level so that it fits
                ! to the RHS, solution vector and matrix again
                call lsysbl_sortVectorInSitu (p_rcurrentLevel%rtempVector,&
                      p_rsubnode%rprjTempVector,.true.)
              end if
              
              ! Step length control. Get the optimal damping parameter for the
              ! defect correction.
              ! Use either the standard system matrix for this purpose or
              ! the specific matrix for the coarse grid correction if this exists
              ! in the level info structure.
              if (p_rcurrentLevel%rsystemMatrixOptCorrection%NEQ .ne. 0) then
                call cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrixOptCorrection,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              else
                call cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrix,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              end if
              
              ! Perform the coarse grid correction by adding the coarse grid
              ! solution (with the calculated step-length parameter) to
              ! the current solution.
              
              call lsysbl_vectorLinearComb (p_rcurrentLevel%rtempVector,&
                                            p_rcurrentLevel%rsolutionVector,&
                                            dstep,1.0_DP)
                            
              ! Extended output
              if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                  (mod(ite,niteResOutput).eq.0)) then
                  
                call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                        p_rcurrentLevel%rtempVector)
                call lsysbl_blockMatVec (&
                    p_rcurrentLevel%rsystemMatrix, &
                    p_rcurrentLevel%rsolutionVector,&
                    p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                if (.not.((dres .ge. 1E-99_DP) .and. &
                          (dres .le. 1E99_DP))) dres = 0.0_DP
                          
                call output_line ('Multigrid: Level '//trim(sys_siL(ilev,5))//&
                    ' after c.g.corr.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
              end if
                                            
              ! Perform the post-smoothing with the current solution vector
              if (associated(p_rcurrentLevel%p_rpostSmoother)) then
                call linsol_smoothCorrection (p_rcurrentLevel%p_rpostSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)

                ! Extended output
                if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                    (mod(ite,niteResOutput).eq.0)) then
                    
                  call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector)
                  call lsysbl_blockMatVec (&
                      p_rcurrentLevel%rsystemMatrix, &
                      p_rcurrentLevel%rsolutionVector,&
                      p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                  dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  call output_line ('Multigrid: Level '//trim(sys_siL(ilev,5))//&
                      ' after postsm.:    !!RES!! = '//trim(sys_sdEL(dres,15)) )
                end if
                                              
              end if

              ! Update the iteration counter(s) for realising the MG-cycle(s).
              ! Then either repeat this loop to perform the next prolongation or
              ! repeat the cycleloop to do perform a next MG sweep on the current
              ! level.
              !
              ! Here icycle defines how the cycle-counters are updated.
              ! For a W-cycle the cycle counter is resetted to 2 if the sweep is
              ! fulfilled on the current level, for F-cycle it is set to 1 to not
              ! perform more that 1 cycle on the current level anymore.
              
              p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncyclesRemaining-1
              if (p_rcurrentLevel%ncyclesRemaining .le. 0) then
                if (p_rsubnode%icycle .eq. 0) then
                  p_rcurrentLevel%ncyclesRemaining = 1
                else
                  ! Cycle finished. Reset counter for next cycle.
                  p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles

                  if ((rsolverNode%p_rsubnodeMultigrid%depsRelCycle .ne. 1E99_DP) .and. &
                      associated(p_rcurrentLevel%p_rnextLevel)) then
                      
                    ! Adaptive cycles activated.
                    !
                    ! We are on a level < nlmax.
                    ! At first, calculate the residuum on that level.
                    call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                            p_rcurrentLevel%rtempVector)
                    call lsysbl_blockMatVec (&
                        p_rcurrentLevel%rsystemMatrix, &
                        p_rcurrentLevel%rsolutionVector,&
                        p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                    dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                        rsolverNode%iresNorm)
                    if (.not.((dres .ge. 1E-99_DP) .and. &
                              (dres .le. 1E99_DP))) dres = 0.0_DP
                              
                    ! Compare it with the initial residuum. If it is not small enough
                    ! and if we have not reached the maximum number of cycles,
                    ! repeat the complete cycle.
                    if ( ((rsolverNode%p_rsubnodeMultigrid%nmaxAdaptiveCycles .le. -1) &
                          .or. &
                          (p_rcurrentLevel%icycleCount .lt. &
                           rsolverNode%p_rsubnodeMultigrid%nmaxAdaptiveCycles)) &
                        .and. &
                        (dres .gt. rsolverNode%p_rsubnodeMultigrid%depsRelCycle * &
                                    p_rcurrentLevel%dinitResCycle) ) then

                      if (rsolverNode%ioutputLevel .ge. 3) then
                        call output_line ( &
                          trim( &
                          sys_siL(p_rcurrentLevel%icycleCount,10)) &
                          //'''th repetition of cycle on level '// &
                          trim(sys_siL(ilev,5))//'.')
                      end if

                      p_rcurrentLevel%icycleCount = p_rcurrentLevel%icycleCount+1
                      cycle cycleloop
                    end if
                    
                    ! Otherwise: The cycle(s) is/are finished;
                    ! the END DO goes up one level.
                    
                  end if

                end if
              else

                ! Next cycle; go down starting from the current level
                cycle cycleloop

              end if
              
            end do ! ilev < nlmax
            
            ! We finally reached the maximum level. Quit the cycle loop
            ! to go on with the next iteration.
            exit cycleloop
          
          end do cycleloop
          
          ! We have (hopefully) successfully performed one MG-sweep, starting
          ! and ending on the finest level. As we are now on the finest level
          ! again, we can update our defect vector to test the current
          ! residuum...
          !
          ! Calculate the residuum and its norm.
          call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rtempVector)
          call lsysbl_blockMatVec (&
                p_rcurrentLevel%rsystemMatrix, &
                p_rcurrentLevel%rsolutionVector,&
                p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
          if (bfilter) then
            ! Apply the filter chain to the vector
            call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                             p_RfilterChain)
          end if
          dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,rsolverNode%iresNorm)
          if (.not.((dres .ge. 1E-99_DP) .and. &
                    (dres .le. 1E99_DP))) dres = 0.0_DP
          
          ! Shift the queue with the last residuals and add the new
          ! residual to it.
          Dresqueue = eoshift(Dresqueue,1,dres)

          rsolverNode%dlastDefect = rsolverNode%dfinalDefect
          rsolverNode%dfinalDefect = dres
          
          ! Test if the iteration is diverged
          if (linsol_testDivergence(rsolverNode,dres)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('Multigrid: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
            call output_line ('Multigrid: Solution diverging!')
            rsolverNode%iresult = 1
            exit
          end if
       
          ! At least perform nminIterations iterations
          if (ite .ge. nminIterations) then
          
            ! Check if the iteration converged
            if (linsol_testConvergence(rsolverNode,dres)) exit
            
          end if

          ! print out the current residuum

          if ((rsolverNode%ioutputLevel .ge. 2) .and. &
              (mod(ite,niteResOutput).eq.0)) then
            call output_line ('Multigrid: Iteration '// &
                trim(sys_siL(ite,10))//',  !!RES!! = '//&
                trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          end if
          
        end do  ! ite
        
        ! Finish - either with an error or if converged.
        ! Print the last residuum.

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (ite .ge. 1) .and. (ite .le. rsolverNode%nmaxIterations) .and. &
            (rsolverNode%iresult .ge. 0)) then
          call output_line ('Multigrid: Iteration '// &
              trim(sys_siL(ite,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if
        
        ! Set ite to NIT to prevent printing of "NIT+1" of the loop was
        ! completed

        if (ite .gt. rsolverNode%nmaxIterations) then
          ! Warning if we did not reach the convergence criterion.
          ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
          if ((rsolverNode%ioutputLevel .ge. 0) .and. &
              (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
            call output_line ('Multigrid: Accuracy warning: '//&
                'Solver did not reach the convergence criterion')
          end if

          ite = rsolverNode%nmaxIterations
        end if

        ! Final number of iterations
        rsolverNode%iiterations = ite
      
      end if

      ! Finally, we look at some statistics:
      !
      ! Do not calculate anything if the final residuum is out of bounds -
      ! would result in NaN`s,...
        
      if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
        ! Calculate asymptotic convergence rate
      
        if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
          i = min(rsolverNode%iiterations,niteAsymptoticCVR)
          rsolverNode%dasymptoticConvergenceRate = &
            (rsolverNode%dfinalDefect / &
              dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
        end if

        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting multigrid convergence rate stays zero.
        ! In the other case the multigrid convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :

        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/real(rsolverNode%iiterations,DP))
        else
          rsolverNode%dconvergenceRate = 0.0_DP
        end if
      end if
    
    end if

    ! As the solution vector on the finest level shared its memory with rd,
    ! we just calculated the new correction vector!
      
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
      ! Scale the defect by the damping parameter in the solver structure.
      call lsysbl_scaleVector (rd,rsolverNode%domega)
      
      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('Multigrid statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'MG: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if

    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    end if
  
  end subroutine
  
! *****************************************************************************
! Routines for the Multigrid-2 solver
! *****************************************************************************
  
  ! ***************************************************************************

!<function>
  
  integer function linsol_getMultigrid2LevelCount (rsolverNode)
                    
!<description>
  ! Returns the number of levels currently attached to the multigrid solver.
!</description>
  
!<input>
  ! The solver structure of the multigrid solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</input>

!<result>
  ! The number of levels in the MG solver.
!</result>
  
!</function>

    linsol_getMultigrid2LevelCount = rsolverNode%p_rsubnodeMultigrid2%nlevels

  end function

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_getMultigrid2Level (rsolverNode,ilevel,p_rlevelInfo)
                    
!<description>
  ! Searches inside of the multigrid structure for level ilevel and returns
  ! a pointer to the corresponding p_rlevelInfo structure-
  ! If the level does not exist, an error is thrown.
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!<input>
  ! Number of the level to fetch.
  integer, intent(in) :: ilevel
!</input>
  
!<output>
  ! A pointer to the corresponding t_levelInfo structure.
  type(t_linsolMG2LevelInfo), pointer     :: p_rlevelInfo
!</output>
  
!</subroutine>

    nullify(p_rlevelInfo)

    ! Do we have the level?
    if ((ilevel .gt. 0) .and. &
        (ilevel .le. rsolverNode%p_rsubnodeMultigrid2%nlevels)) then
      ! Get it.
      p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    else
      call output_line('Level out of range!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'linsol_getMultigrid2Level')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneMultigrid2Level (rlevelInfo)
                    
!<description>
  ! Cleans up the multigrid level structure rlevelInfo. All memory allocated
  ! in this structure is released.
!</description>
  
!<inputoutput>
  ! The t_levelInfo structure to clean up.
  type(t_linsolMG2LevelInfo), intent(inout)     :: rlevelInfo
!</inputoutput>

!</subroutine>
  
    ! Release the temporary vectors of this level.
    ! The vector may not exist - if the application has not called
    ! initStructure!
    if (rlevelInfo%rrhsVector%NEQ .ne. 0) &
      call lsysbl_releaseVector(rlevelInfo%rrhsVector)

    if (rlevelInfo%rtempVector%NEQ .ne. 0) &
      call lsysbl_releaseVector(rlevelInfo%rtempVector)

    if (rlevelInfo%rsolutionVector%NEQ .ne. 0) &
      call lsysbl_releaseVector(rlevelInfo%rsolutionVector)

    ! Release sub-solvers if associated.
    !
    ! Caution: Pre- and postsmoother may be identical!
    if (associated(rlevelInfo%p_rpostSmoother) .and. &
        (.not. associated(rlevelInfo%p_rpreSmoother, rlevelInfo%p_rpostSmoother))) then
      call linsol_releaseSolver(rlevelInfo%p_rpostsmoother)
    end if

    if (associated(rlevelInfo%p_rpresmoother)) then
      call linsol_releaseSolver(rlevelInfo%p_rpresmoother)
    end if

    if (associated(rlevelInfo%p_rcoarseGridSolver)) then
      call linsol_releaseSolver(rlevelInfo%p_rcoarseGridSolver)
    end if
    
    ! Clean up the associated matrix if there is one.
    ! Of course, the memory of the matrices will not be deallocated
    ! because the matrices belong to the application, not to the solver!
    call lsysbl_releaseMatrix(rlevelInfo%rsystemMatrix)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_cleanMultigrid2Levels (rsolverNode)
                    
!<description>
  ! This routine removes all level information from the MG solver and releases
  ! all attached solver nodes on all levels (smoothers, coarse grid solvers,
  ! ...).
  !
  ! The level info array is not released.
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver.
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

  integer :: ilevel

    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_cleanMultigrid2Levels')
      call sys_halt()
    end if

    ! Remove all the levels
    do ilevel = 1,rsolverNode%p_rsubnodeMultigrid2%nlevels
      call linsol_doneMultigrid2Level (&
          rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel))
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initProjMG2LvByPrj (rlevelInfo,rinterlevelProjection)
  
!<description>
  ! Initialises the interlevel projection structure for one level.
!</description>
  
!<inputoutput>
  ! The t_levelInfo structure to set up.
  type(t_linsolMG2LevelInfo), intent(inout)     :: rlevelInfo
  
  ! User-specified interlevel projection structure.
  ! Multigrid will use the given projection structure for doing prolongation/
  ! restriction.
  type(t_interlevelProjectionBlock), intent(in), target  :: rinterlevelProjection
!</inputoutput>

!</subroutine>

    ! Check if the pointer to the projection structure is already assigned.
    ! If that is the case, release it.
    if (rlevelInfo%bstandardProjection) then
      ! Release the old standard projection structure.
      if (associated(rlevelInfo%p_rprojection)) then
        call mlprj_doneProjection (rlevelInfo%p_rprojection)
        deallocate (rlevelInfo%p_rprojection)
      end if
    end if
    
    ! Remember the user specified projection
    rlevelInfo%p_rprojection => rinterlevelProjection
    rlevelInfo%bstandardProjection = .false.

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initProjMG2LvByDiscr (rlevelInfo,rdiscrCoarse,rdiscrFine)
  
!<description>
  ! Initialises the interlevel projection structure for one level
  ! based on coarse/fine grid discretisation structures.
!</description>
  
!<inputoutput>
  ! The t_levelInfo structure to set up.
  type(t_linsolMG2LevelInfo), intent(inout)     :: rlevelInfo
!</inputoutput>

!</input>
  ! Discretisation structure of the coarse level.
  type(t_blockDiscretisation), intent(in), target   :: rdiscrCoarse

  ! Discretisation structure of the fine level.
  type(t_blockDiscretisation), intent(in), target   :: rdiscrFine
!</input>

!</subroutine>

    ! Check if the pointer to the projection structure is already assigned.
    ! If that is the case, release it.
    if (rlevelInfo%bstandardProjection) then
      ! Release the old standard projection structure.
      if (associated(rlevelInfo%p_rprojection)) then
        call mlprj_doneProjection (rlevelInfo%p_rprojection)
      else
        ! Allocate a new projection structure
        allocate(rlevelInfo%p_rprojection)
      end if
    else
      ! Allocate a new projection structure
      rlevelInfo%bstandardProjection = .true.
      allocate(rlevelInfo%p_rprojection)
    end if
    
    ! Initialise the projection structure with standard parameters.
    call mlprj_initProjectionDiscr (rlevelInfo%p_rprojection,rdiscrCoarse)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initProjMultigrid2 (rsolverNode)
  
!<description>
  ! Initialises the interlevel projection structure for all levels
  ! using the standard parameters for prolongation/restriction.
  !
  ! The routine can be called only after the matrices are attached to
  ! the solver with the setMatrices routine. It will automatically detect
  ! if the caller already specified user defined projection structures
  ! by a call to linsol_initProjMultigrid2LvByPrj().
  !
  ! The routine is called during initStructure but can also be called
  ! manually in advance after setMatrices() if the caller wants to
  ! alter the projection structures manually.
!</description>
  
!<inputoutput>
  ! A t_linsolNode structure. defining the solver.
  type(t_linsolNode), intent(inout)             :: rsolverNode
!</inputoutput>

!</subroutine>

    integer :: ilevel, nlmax
    type(t_blockDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine
    
    nlmax = rsolverNode%p_rsubnodeMultigrid2%nlevels
    
    ! Loop through all levels. Whereever the projection structure
    ! is missing, create it.
    do ilevel = 2,nlmax
      if (.not. associated( &
          rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)% &
          p_rprojection)) then
        
        ! Initialise the projection structure.
        p_rdiscrCoarse => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel-1)%rsystemmatrix%&
          p_rblockDiscrTrial
        p_rdiscrFine => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)%rsystemmatrix%&
          p_rblockDiscrTrial
        call linsol_initProjMG2LvByDiscr (rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel),&
            p_rdiscrCoarse,p_rdiscrFine)
        
      end if
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneProjMultigrid2 (rsolverNode)
  
!<description>
  ! Releases all automatically allocated interlevel projection structures
  ! for all levels in the MG solver.
!</description>
  
!<inputoutput>
  ! A  t_linsolNode structure. defining the solver.
  type(t_linsolNode), intent(inout)             :: rsolverNode
!</inputoutput>

!</subroutine>

    integer :: ilevel, nlmax
    
    nlmax = rsolverNode%p_rsubnodeMultigrid2%nlevels
    
    ! Loop through all levels. Whereever the projection structure
    ! is missing, create it.
    do ilevel = 2,nlmax
      if (associated(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)% &
          p_rprojection)) then
        if (rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)% &
            bstandardProjection) then
          call mlprj_doneProjection (rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)%&
              p_rprojection)
          deallocate (rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)%p_rprojection)
        else
          ! Just nullify the pointer, the structure does not belong to us.
          nullify (rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)%p_rprojection)
        end if
        
        rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)% &
            bstandardProjection = .false.
        
      end if
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initMultigrid2 (p_rsolverNode,nlevels,Rfilter)
  
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
  integer, intent(in) :: nlevels
  
  ! Optional: A pointer to a filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not given or set to NULL(), no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter
  
!</input>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer             :: p_rsolverNode
   
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
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
    allocate(p_rsolverNode%p_rsubnodeMultigrid2)
    
    ! Initialise the coarse grid correction structure.1
    call cgcor_init (p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
    
    ! Allocate level information data
    allocate(p_rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(nlevels))
    p_rsolverNode%p_rsubnodeMultigrid2%nlevels = nlevels
    
    ! Attach the filter if given.
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubnodeMultigrid2%p_RfilterChain => Rfilter
    end if
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_alterMultigrid2 (rsolverNode, ralterConfig)
  
!<description>
  ! This routine allows on-line modification of the Multigrid solver.
  ! ralterConfig%ccommand is analysed and depending on the configuration
  ! in this structure, the solver reacts.
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_linsolNode), intent(inout), target             :: rsolverNode

  ! A command/configuration block that specifies a command which is given
  ! to the solver.
  type(t_linsol_alterSolverConfig), intent(inout)       :: ralterConfig
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel

    ! Pass the command structure to all subsolvers and smoothers
    ! on all levels.
    do i=1,rsolverNode%p_rsubnodeMultigrid2%nlevels
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(i)
      
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_alterSolver (p_rcurrentLevel%p_rpreSmoother, &
                                 ralterConfig)
      end if
      ! Pre- and postsmoother may be identical!
      ! Take care not to alter the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_alterSolver (p_rcurrentLevel%p_rpostSmoother, &
                                 ralterConfig)
      end if
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_alterSolver (p_rcurrentLevel%p_rcoarseGridSolver, &
                                 ralterConfig)
      end if
      
    end do
    
  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_setMatrixMultigrid2 (rsolverNode,Rmatrices)
  
!<description>
  
  ! Assigns the system matrix rsystemMatrix to the Multigrid solver on
  ! all levels.
  
!</description>
  
!<input>
  
  ! An array of system matrices on all levels.
  ! Each level in multigrid is initialised separately.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices

!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!</subroutine>

  ! local variables
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel
  integer                       :: ilevel,nlmax
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_setMatrixMultigrid2')
      call sys_halt()
    end if
    
    ! Make sure we have the right amount of matrices
    if (size(Rmatrices) .ne. rsolverNode%p_rsubnodeMultigrid2%nlevels) then
      call output_line ('Wrong number of matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_setMatrixMultigrid2')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached. If yes, attach the matrix
    ! to it.
    ! For this purpose, call the general initialisation routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    
    nlmax = rsolverNode%p_rsubnodeMultigrid2%nlevels
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_setMatrices (p_rcurrentLevel%p_rpreSmoother, &
                                  Rmatrices(1:ilevel))
      end if
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_setMatrices (p_rcurrentLevel%p_rpostSmoother, &
                                  Rmatrices(1:ilevel))
      end if
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_setMatrices (p_rcurrentLevel%p_rcoarseGridSolver, &
                                  Rmatrices(1:ilevel))
      end if
      
      ! Write a link to the matrix into the level-structure.
      call lsysbl_duplicateMatrix (Rmatrices(ilevel), &
          p_rcurrentLevel%rsystemMatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_matCompatMultigrid2 (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  ! More precisely, ccompatible returns always the status of the highest
  ! level where there is an error -- or LINSOL_COMP_OK if there is no error.
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  ! local variables
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel
  integer                       :: ilevel,ccompatLevel
    
    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK
    
    ! Reset all compatibility flags
    if (present(CcompatibleDetail)) CcompatibleDetail (:) = ccompatible

    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_matCompatMultigrid2')
      call sys_halt()
    end if
    
    ! Make sure we have the right amount of matrices
    if (size(Rmatrices) .ne. rsolverNode%p_rsubnodeMultigrid2%nlevels) then
      call output_line ('Wrong number of matrices!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_matCompatMultigrid2')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached. If yes, it must check
    ! the matrices!
    ! For this purpose, call the general check routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    
    do ilevel = 1,rsolverNode%p_rsubnodeMultigrid2%nlevels
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    
      ! Compatibility flag for that level
      ccompatLevel = LINSOL_COMP_OK
    
      ! Check presmoother
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
      
        if (present(CcompatibleDetail)) then
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(1:ilevel),ccompatLevel,CcompatibleDetail(1:ilevel))
        else
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpreSmoother, &
              Rmatrices(1:ilevel),ccompatLevel)
        end if
      
        if (ccompatLevel .ne. LINSOL_COMP_OK) then
          ccompatible = ccompatLevel
          cycle
        end if
        
      end if
      
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        
        if (present(CcompatibleDetail)) then
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(1:ilevel),ccompatLevel,CcompatibleDetail(1:ilevel))
        else
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpostSmoother, &
              Rmatrices(1:ilevel),ccompatLevel)
        end if
        
        if (ccompatLevel .ne. LINSOL_COMP_OK) then
          ccompatible = ccompatLevel
          cycle
        end if
        
      end if
      
      ! Check coarse grid solver
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        
        if (present(CcompatibleDetail)) then
          call linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(1:ilevel),ccompatLevel,CcompatibleDetail(1:ilevel))
        else
          call linsol_matricesCompatible (p_rcurrentLevel%p_rcoarseGridSolver, &
              Rmatrices(1:ilevel),ccompatLevel)
        end if
        
        if (ccompatLevel .ne. LINSOL_COMP_OK) then
          ccompatible = ccompatLevel
          cycle
        end if
        
      end if
      
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initStructureMultigrid2 (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel,p_rprevLevel
  integer :: isubgroup,nlmax,ilevel
  integer :: imemmax
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_vectorBlock), pointer :: p_rtemplVect
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.

    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initStructureMultigrid2')
      call sys_halt()
    end if
    
    nlmax = rsolverNode%p_rsubnodeMultigrid2%nlevels

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_initStructure(p_rcurrentLevel%p_rpreSmoother,ierror,isubgroup)
      end if
      
      ! Cancel in case of an error
      if (ierror .ne. LINSOL_ERR_NOERROR) return
      
      ! Pre- and postsmoother may be identical!
      ! Take care not to initialise the same smoother twice!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_initStructure(p_rcurrentLevel%p_rpostSmoother,ierror,isubgroup)

        ! Cancel in case of an error
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_initStructure(p_rcurrentLevel%p_rcoarseGridSolver,ierror,isubgroup)

        ! Cancel in case of an error
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if

    end do
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Initialise the interlevel projection as far as it is not already initialised
    call linsol_initProjMultigrid2 (rsolverNode)

    imemmax = 0

    ! Multigrid needs temporary vectors for the RHS, solution and general
    ! data. Memory for that is allocated on all levels for temporary, RHS and
    ! solution vectors.
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

      p_rmatrix => p_rcurrentLevel%rsystemMatrix
      nullify(p_rtemplVect)
      if (ilevel .lt. nlmax) then
        ! On the maximum level we do not need additional memory for the
        ! solution vector - as solution vector on the fine grid,
        ! the vector which comes as parameter
        ! to the multigrid  preconditioner is used.
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
              p_rcurrentLevel%rsolutionVector,.false.,.false.,&
              rsolverNode%cdefaultDataType)
        p_rtemplVect => p_rcurrentLevel%rsolutionVector
      end if
      ! On the coarse grid, we do not need memory for temporary/RHS
      ! vectors. The RHS on the coarse grid is replaced in-situ
      ! by the solution vector.
      if (ilevel .gt. 1) then
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
              p_rcurrentLevel%rrhsVector,.false.,.false.,&
              rsolverNode%cdefaultDataType)
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
              p_rcurrentLevel%rtempVector,.false.,.false.,&
              rsolverNode%cdefaultDataType)
        p_rtemplVect => p_rcurrentLevel%rtempVector
      end if

      ! Calculate the memory that is probably necessary for resorting
      ! vectors during the iteration.
      if (associated(p_rtemplVect)) then
        imemmax = max(imemmax,p_rtemplVect%NEQ)
      end if
         
      ! Ask the prolongation/restriction routines how much memory is necessary
      ! for prolongation/restriction on that level.
      if (ilevel .gt. 1) then ! otherwise coarse grid
        p_rprevLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel-1)
     
        ! Calculate the memory that is necessary for prolongation/restriction -
        ! in case there is temporary memory needed.
        ! The system matrix on the fine/coarse grid specifies the discretisation.
        imemmax = max(imemmax,mlprj_getTempMemoryMat ( &
                      p_rcurrentLevel%p_rprojection, &
                      p_rprevLevel%rsystemMatrix,&
                      p_rcurrentLevel%rsystemMatrix))
      end if

      ! All temporary vectors are marked as 'unsorted'. We can set the
      ! sorting flag directly here without using the resorting routines
      ! as the vectors have just been created.
      !
      ! Do not do anything to the RHS/temp vectors on the coarse grid
      ! and the solution vector on the fine grid -- they do not exist!
      if (ilevel .gt. 1) then
        p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy)
        p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy)
      end if
      
      if (ilevel .lt. NLMAX) then
        p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy)
      end if
      
      ! And the next level...
    end do
    
    ! Do we have to allocate memory for prolongation/restriction/resorting?
    if (imemmax .gt. 0) then
      call lsyssc_createVector (rsolverNode%p_rsubnodeMultigrid2%rprjTempVector,&
                                imemmax,.false.,rsolverNode%cdefaultDataType)
                                
      ! Use this vector also for the coarse grid correction.
      ! Create a block vector on every level (except for the coarse grid)
      ! with the structure of the RHS that shares the memory with rprjTempVector.
      do ilevel = 2,nlmax
    
        p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
      
        call lsysbl_createVecFromScalar (&
            rsolverNode%p_rsubnodeMultigrid2%rprjTempVector,&
            p_rcurrentLevel%rcgcorrTempVector)
      
        call lsysbl_enforceStructure (&
            p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rcgcorrTempVector)

      end do
    else
      rsolverNode%p_rsubnodeMultigrid2%rprjTempVector%NEQ = 0
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initDataMultigrid2 (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup,ilevel
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup
    
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataMultigrid2')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    
    do ilevel = 1,rsolverNode%p_rsubnodeMultigrid2%nlevels
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)

      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_initData(p_rcurrentLevel%p_rpreSmoother,ierror,isubgroup)

        ! Cancel in case of an error
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
      ! Pre- and postsmoother may be identical!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_initData(p_rcurrentLevel%p_rpostSmoother,ierror,isubgroup)

        ! Cancel in case of an error
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_initData(p_rcurrentLevel%p_rcoarseGridSolver,ierror,isubgroup)

        ! Cancel in case of an error
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if

    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataMultigrid2 (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup,ilevel
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_doneDataMultigrid2')
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    
    do ilevel = rsolverNode%p_rsubnodeMultigrid2%nlevels,1,-1
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
      
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_doneData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
      end if
      
      ! Pre- and postsmoother may be identical!
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_doneData(p_rcurrentLevel%p_rpostSmoother,isubgroup)
      end if
      
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_doneData(p_rcurrentLevel%p_rpreSmoother,isubgroup)
      end if
      
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureMultigrid2 (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup,ilevel,nlmax
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_doneStructureMultigrid2')
      call sys_halt()
    end if
    
    nlmax = rsolverNode%p_rsubnodeMultigrid2%nlevels

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    
      ! Release the coarse grid solver
      if (associated(p_rcurrentLevel%p_rcoarseGridSolver)) then
        call linsol_doneStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
      end if
      
      ! Pre- and postsmoother may be identical! If not, first release the
      ! postsmoother.
      if (associated(p_rcurrentLevel%p_rpostSmoother) .and. &
          (.not. associated(p_rcurrentLevel%p_rpreSmoother, &
                            p_rcurrentLevel%p_rpostSmoother))) then
        call linsol_doneStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
      end if
      
      ! Release the presmoother
      if (associated(p_rcurrentLevel%p_rpreSmoother)) then
        call linsol_doneStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
      end if
      
    end do

    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Release temporary memory.
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilevel)
    
      if (ilevel .gt. 1) then
      
        if (p_rcurrentLevel%rtempVector%NEQ .ne. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rtempVector)
        end if
        
        if (p_rcurrentLevel%rrhsVector%NEQ .ne. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rrhsVector)
        end if

        ! Release the temp vector for the coarse grid correction.
        ! The vector shares its memory with the 'global' projection
        ! vector, so only release it if the other vector exists.
        if (p_rcurrentLevel%rcgcorrTempVector%NEQ .gt. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rcgcorrTempVector)
        end if

      end if

      if (ilevel .lt. nlmax) then
        if (p_rcurrentLevel%rsolutionVector%NEQ .ne. 0) then
          call lsysbl_releaseVector (p_rcurrentLevel%rsolutionVector)
        end if
      end if

    end do

    ! Check if we have to release our temporary vector for prolongation/
    ! restriction.
    ! Do we have to allocate memory for prolongation/restriction?
    if (rsolverNode%p_rsubnodeMultigrid2%rprjTempVector%NEQ .gt. 0) then
      call lsyssc_releaseVector (rsolverNode%p_rsubnodeMultigrid2%rprjTempVector)
    end if
    
    ! Release the interlevel projection structures as far as they were
    ! automatically created.
    call linsol_doneProjMultigrid2 (rsolverNode)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneMultigrid2 (rsolverNode)
  
!<description>
  
  ! This routine releases all temporary memory for the multigrid solver from
  ! the heap. In particular, this releases the solver structures of all
  ! subsolvers (smoother, coarse grid solver).
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
  
!</description>
  
!<inputoutput>
  ! A pointer to a t_linsolNode structure of the Multigrid solver.
  type(t_linsolNode), intent(inout)                :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_MULTIGRID2) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid2))) then
      call output_line ('Multigrid structure not initialised!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_doneMultigrid2')
      call sys_halt()
    end if

    ! Release data and structure if this has not happened yet.
    call linsol_doneDataMultigrid2(rsolverNode)
    call linsol_doneStructureMultigrid2(rsolverNode)
    
    ! Remove all the levels
    call linsol_cleanMultigrid2Levels (rsolverNode)
    deallocate(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)
    
    ! Release the coarse grid correction structure
    call cgcor_release (rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
    
    ! Release MG substructure, that is it.
    deallocate(rsolverNode%p_rsubnodeMultigrid2)

  end subroutine
    
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_precMultigrid2 (rsolverNode,rd)
  
!<description>
  ! Applies Multigrid preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The structure of the levels (pre/postsmoothers, interlevel projection
  ! structures) must have been prepared by setting the level parameters
  ! before calling this routine. (The parameters of a specific level
  ! can be get by calling linsol_getMultigridLevel2).
  ! The matrices <tex>$ A $</tex> on all levels must be attached to the solver previously
  ! by linsol_setMatrices.
  
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: nminIterations,nmaxIterations,niteResOutput,niteAsymptoticCVR
  integer :: ite,ilev,nlmax,i,j,nblocks
  real(DP) :: dres,dstep
  logical :: bfilter,bsort
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
  ! The queue saves the current residual and the two previous residuals.
  real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
  
  ! The system matrix on the current level
  type(t_matrixBlock), pointer :: p_rmatrix
  
  ! Our MG structure
  type(t_linsolSubnodeMultigrid2), pointer :: p_rsubnode
  
  ! The current level and the next lower one.
  type(t_linsolMG2LevelInfo), pointer :: p_rcurrentLevel,p_rlowerLevel
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeMultigrid2
    
    bfilter = associated(p_rsubnode%p_RfilterChain)
    p_RfilterChain => p_rsubnode%p_RfilterChain
    
    ! Get the system matrix on the finest level:
    p_rmatrix => rsolverNode%rsystemMatrix
    
    ! Check the parameters
    if ((rd%NEQ .eq. 0) .or. (p_rmatrix%NEQ .eq. 0) .or. &
        (p_rmatrix%NEQ .ne. rd%NEQ) ) then
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    if (p_rsubnode%icycle .lt. 0) then
      ! Wrong cycle
      call output_line ('Invalid cycle!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precMultigrid2')
      rsolverNode%iresult = 2
      return
    end if
    
    if (.not. associated(rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo)) then
      call output_line ('No levels attached!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precMultigrid2')
      rsolverNode%iresult = 2
      return
    end if

    ! Maximum level
    nlmax = rsolverNode%p_rsubnodeMultigrid2%nlevels

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate
    niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))

    ! Minimum/maximum number of iterations
    nminIterations = max(rsolverNode%nminIterations,0)
    nmaxIterations = max(rsolverNode%nmaxIterations,0)
      
    ! Iteration when the residuum is printed:
    niteResOutput = max(1,rsolverNode%niteResOutput)

    ! Status reset
    rsolverNode%iresult = 0
    rsolverNode%icurrentIteration = 0
    rsolverNode%dinitialDefect = 0.0_DP
    rsolverNode%dlastDefect = 0.0_DP
    rsolverNode%dfinalDefect = 0.0_DP
    rsolverNode%dconvergenceRate = 0.0_DP
    rsolverNode%dasymptoticConvergenceRate = 0.0_DP
    
    ! Number of blocks in the current equation
    nblocks = rd%nblocks
    
    ! We start on the maximum level.
    p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(nlmax)
    
    ! Is there only one level? Can be seen if the current level
    ! already contains a coarse grid solver.
    if (nlmax .eq. 1) then
    
      if (rsolverNode%ioutputLevel .gt. 1) then
        call output_line ('Multigrid: Only one level. '//&
             'Switching back to standard solver.')
      end if
      call linsol_precondDefect(p_rcurrentLevel%p_rcoarseGridSolver,rd)
      
      ! Take the statistics from the coarse grid solver.
      rsolverNode%dinitialDefect = p_rcurrentLevel%p_rcoarseGridSolver%dinitialDefect
      rsolverNode%dfinalDefect = p_rcurrentLevel%p_rcoarseGridSolver%dfinalDefect
      rsolverNode%dlastDefect = p_rcurrentLevel%p_rcoarseGridSolver%dlastDefect
      rsolverNode%dconvergenceRate = p_rcurrentLevel%p_rcoarseGridSolver%dconvergenceRate
      rsolverNode%iiterations = p_rcurrentLevel%p_rcoarseGridSolver%iiterations
      rsolverNode%dasymptoticConvergenceRate = p_rcurrentLevel% &
        p_rcoarseGridSolver%dasymptoticConvergenceRate
      
    else
      ! True multigrid.
      !
      ! Get the norm of the initial residuum.
      ! As the initial iteration vector is zero, this is the norm
      ! of the RHS:
      dres = lsysbl_vectorNorm (rd,rsolverNode%iresNorm)
      if (.not.((dres .ge. 1E-99_DP) .and. &
                (dres .le. 1E99_DP))) dres = 0.0_DP

      ! Initialise starting residuum
      rsolverNode%dinitialDefect = dres
      rsolverNode%dlastDefect = 0.0_DP

      ! initialise the queue of the last residuals with RES
      Dresqueue = dres

      ! Check if out initial defect is zero. This may happen if the filtering
      ! routine filters "everything out"!
      ! In that case we can directly stop our computation.

      if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
      
        ! final defect is 0, as initialised in the output variable above
        call lsysbl_clearVector(rd)
        rsolverNode%dfinalDefect = dres
        rsolverNode%dconvergenceRate = 0.0_DP
        rsolverNode%iiterations = 0
        rsolverNode%dasymptoticConvergenceRate = 0.0_DP
        
      else
        
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
        do i=1,nlmax-1
          p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(i)
          if (p_rsubnode%icycle .eq. 0) then
            p_rcurrentLevel%ncycles = 2
          else
            p_rcurrentLevel%ncycles = p_rsubnode%icycle
          end if
        end do
        
        ! We start at the maximum level.
        ilev = nlmax

        ! Get current and next lower level.
        p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev)
        p_rlowerLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev-1)
               
        ! On the current (max.) level we set ncycles to 1.
        p_rcurrentLevel%ncycles = 1
        
        ! Print out the initial residuum

        if (rsolverNode%ioutputLevel .ge. 2) then
          call output_line ('Multigrid: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        end if

        ! Copy the initial RHS to the RHS vector on the maximum level.
        call lsysbl_copyVector(rd,p_rcurrentLevel%rrhsVector)
        
        ! Replace the solution vector on the finest level by rd.
        ! Afterwards, rd and the solution vector on the finest level
        ! share the same memory location, so we use rd as iteration
        ! vector directly.
        call lsysbl_duplicateVector (rd,p_rcurrentLevel%rsolutionVector,&
            LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
        
        ! Clear the initial solution vector.
        call lsysbl_clearVector (p_rcurrentLevel%rsolutionVector)
        
        ! Start multigrid iteration; perform at most nmaxiterations iterations.
        do ite = 1, nmaxiterations
        
          rsolverNode%icurrentIteration = ite
          
          ! Initialise cycle counters for all levels.
          do i=1,nlmax
            p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(i)
            p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles
          end do
        
          ! p_rcurrentLevel now points to the maximum level again!
          !
          ! Build the defect...
          call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rtempVector)
          if (ite .ne. 1) then   ! initial solution vector is zero!
            call lsysbl_blockMatVec (&
                 p_rcurrentLevel%rsystemMatrix, &
                 p_rcurrentLevel%rsolutionVector,&
                 p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
          end if
          if (bfilter) then
            ! Apply the filter chain to the vector
            call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                             p_RfilterChain)
          end if
          
          cycleloop: do  ! Loop for the cycles
          
            ! On the maximum level we already built out defect vector. If we are
            ! on a lower level than NLMAX, perform smoothing+restriction down to the
            ! coarse level. We identify the coarse level by checking if
            ! the current level has a coarse grid solver.
            
            do while (ilev .gt. 1)
            
              ! Perform the pre-smoothing with the current solution vector
              if (associated(p_rcurrentLevel%p_rpreSmoother)) then
                call linsol_smoothCorrection (p_rcurrentLevel%p_rpreSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)
              end if
            
              ! Build the defect vector
              call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                      p_rcurrentLevel%rtempVector)
              call lsysbl_blockMatVec (&
                  p_rcurrentLevel%rsystemMatrix, &
                  p_rcurrentLevel%rsolutionVector,&
                  p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
              if (bfilter) then
                ! Apply the filter chain to the vector
                call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                                p_RfilterChain)
              end if
              
              ! Extended output
              if (associated(p_rcurrentLevel%p_rpreSmoother) .and. &
                  (rsolverNode%ioutputLevel .ge. 3) .and. &
                  (mod(ite,niteResOutput).eq.0)) then
                  
                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                if (.not.((dres .ge. 1E-99_DP) .and. &
                          (dres .le. 1E99_DP))) dres = 0.0_DP
                          
                call output_line ('Multigrid: Level '//trim(sys_siL(ilev,5))//&
                    ' after presm.:     !!RES!! = '//trim(sys_sdEL(dres,15)) )
              end if

              ! Restriction of the defect. The restricted defect is placed
              ! in the right hand side vector of the lower level.
              ! The projection parameters from level ilev to level
              ! ilev-1 is configured in the rprojection structure of the
              ! current level.
              !
              ! Make sure the vector is unsorted before the restriction!
              ! We can use the 'global' temporary array p_rprjTempVector
              ! as temporary memory during the resorting process.
              if (bsort) then
                call lsysbl_sortVectorInSitu (p_rcurrentLevel%rtempVector,&
                      p_rsubnode%rprjTempVector,.false.)
              end if
              
              ! When restricting to the coarse grid, we directy restrict into
              ! the solution vector. It is used there as RHS and replaced in-situ
              ! by the solution by the coarse grid solver. So do not need the RHS vector
              ! on the coarse grid and save one vector-copy.
              !
              ! Otherwise, we restrict to the RHS on the lower level and continue
              ! the smoothing process there.
              if (ilev .gt. 2) then
              
                ! We do not project to the coarse grid
                call mlprj_performRestriction (p_rcurrentLevel%p_rprojection,&
                      p_rlowerLevel%rrhsVector, &
                      p_rcurrentLevel%rtempVector, &
                      p_rsubnode%rprjTempVector)

                ! Apply the filter chain (e.g. to implement boundary conditions)
                ! on just calculated right hand side
                ! (which is a defect vector against the 0-vector).
                if (bfilter) then
                  ! Apply the filter chain to the vector.
                  ! We are in 'unsorted' state; applying the filter here is
                  ! supposed to be a bit faster...
                  call filter_applyFilterChainVec (p_rlowerLevel%rrhsVector, &
                                                   p_RfilterChain)
                end if

                if (bsort) then
                  ! Resort the RHS on the lower level.
                  call lsysbl_sortVectorInSitu (p_rlowerLevel%rrhsVector,&
                        p_rsubnode%rprjTempVector,.true.)

                  ! Temp-vector and solution vector there are yet uninitialised,
                  ! therefore we can mark them as sorted without calling the
                  ! resorting routine.
                  p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                  abs(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy)
                  
                  p_rlowerLevel%rsolutionVector%RvectorBlock(1:nblocks)% &
                    isortStrategy = abs(p_rlowerLevel%rsolutionVector% &
                    RvectorBlock(1:nblocks)%isortStrategy)
                end if

                ! Choose zero as initial vector on lower level.
                call lsysbl_clearVector (p_rlowerLevel%rsolutionVector)
                
                ! Extended output and/or adaptive cycles
                if ((rsolverNode%p_rsubnodeMultigrid2%depsRelCycle .ne. 1E99_DP) .or.&
                    (rsolverNode%ioutputLevel .ge. 3)) then
                
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rrhsVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  ! In case adaptive cycles are activated, save the 'initial' residual
                  ! of that level into the level structure. Then we can later check
                  ! if we have to repeat the cycle on the coarse mesh.
                  if (rsolverNode%p_rsubnodeMultigrid2%depsRelCycle .ne. 1E99_DP) then
                    p_rlowerLevel%dinitResCycle = dres
                    p_rlowerLevel%icycleCount = 1
                  end if
                         
                  ! If the output level is high enough, print that residuum norm.
                  if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                      (mod(ite,niteResOutput).eq.0)) then
                    call output_line ('Multigrid: Level '//trim(sys_siL(ilev-1,5))//&
                        ' after restrict.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
                  end if
                  
                end if

              else
              
                ! The vector is to be restricted to the coarse grid.
                call mlprj_performRestriction (p_rcurrentLevel%p_rprojection,&
                      p_rlowerLevel%rsolutionVector, &
                      p_rcurrentLevel%rtempVector, &
                      p_rsubnode%rprjTempVector)

                ! Apply the filter chain (e.g. to implement boundary conditions)
                ! on just calculated right hand side
                ! (which is a defect vector against the 0-vector).
                if (bfilter) then
                  ! Apply the filter chain to the vector.
                  ! We are in 'unsorted' state; applying the filter here is
                  ! supposed to be a bit faster...
                  call filter_applyFilterChainVec (p_rlowerLevel%rsolutionVector, &
                                                   p_RfilterChain)
                end if

                if (bsort) then
                  ! Resort the RHS on the lower level.
                  call lsysbl_sortVectorInSitu (p_rlowerLevel%rsolutionVector,&
                        p_rsubnode%rprjTempVector,.true.)

                  ! Temp-vector and RHS can be ignored on the coarse grid.
                end if

                ! Extended output
                if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                    (mod(ite,niteResOutput).eq.0)) then
                    
                  dres = lsysbl_vectorNorm (p_rlowerLevel%rsolutionVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  call output_line ('Multigrid: Level '//trim(sys_siL(ilev-1,5))//&
                      ' after restrict.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
                end if

              end if
            
              ! Go down one level
              ilev = ilev - 1
              p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev)
              
              if (ilev .ne. 1) &
                p_rlowerLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev-1)
              
              ! If we are not on the lowest level, repeat the smoothing of
              ! the solution/restriction of the new defect in the next loop
              ! pass...
            end do   ! ilev > minimum level
            
            ! Now we reached the coarse grid.
            ! Apply the filter chain (e.g. to implement boundary conditions)
            ! on the just calculated right hand side,
            ! which is currently located in rsolutionVector.
            if (bfilter) then
              ! Apply the filter chain to the vector
              call filter_applyFilterChainVec (p_rcurrentLevel%rsolutionVector, &
                                               p_RfilterChain)
            end if
            
            ! Solve the system on lowest level by preconditioning
            ! of the RHS=defect vector.
            call linsol_precondDefect(p_rcurrentLevel%p_rcoarseGridSolver,&
                                      p_rcurrentLevel%rsolutionVector)
            
            ! Now we have the solution vector on the lowest level - we have to go
            ! upwards now... but probably not to NLMAX! That depends on the cycle.
            !
            do while(ilev .lt. nlmax)
            
              ilev = ilev + 1
              p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev)
              p_rlowerLevel => rsolverNode%p_rsubnodeMultigrid2%p_RlevelInfo(ilev-1)
              
              ! Prolongate the solution vector from the coarser level
              ! to the temp vector on the finer level.
              !
              ! Make sure the vector is unsorted before the prolongation!
              ! We can use the 'global' temporary array p_rprjTempVector
              ! as temporary memory during the resorting process.
              if (bsort) then
                call lsysbl_sortVectorInSitu (p_rlowerLevel%rsolutionVector,&
                      p_rsubnode%rprjTempVector,.false.)
                
                ! Temp-vector and RHS vector there unused now,
                ! therefore we can mark them as not sorted without calling the
                ! resorting routine.
                ! This prepares these vectors for the next sweep, when an unsorted
                ! vector comes 'from above'.
                !
                ! Note that this shall not be done on the coarse grid as there is
                ! no temp/rhs vector!
                if (ilev .gt. 2) then
                  p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)%isortStrategy = &
                    -abs(p_rlowerLevel%rtempVector%RvectorBlock(1:nblocks)% &
                    isortStrategy)
                  p_rlowerLevel%rrhsVector%RvectorBlock(1:nblocks)% &
                    isortStrategy = -abs(p_rlowerLevel%rrhsVector% &
                    RvectorBlock(1:nblocks)%isortStrategy)
                end if
              end if
              call mlprj_performProlongation (p_rcurrentLevel%p_rprojection,&
                    p_rlowerLevel%rsolutionVector, &
                    p_rcurrentLevel%rtempVector, &
                    p_rsubnode%rprjTempVector)

              if (bfilter) then
                ! Apply the filter chain to the vector.
                call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                                 p_RfilterChain)
              end if

              if (bsort) then
                ! Resort the temp vector on the current level so that it fits
                ! to the RHS, solution vector and matrix again
                call lsysbl_sortVectorInSitu (p_rcurrentLevel%rtempVector,&
                      p_rsubnode%rprjTempVector,.true.)
              end if
              
              ! Step length control. Get the optimal damping parameter for the
              ! defect correction.
              ! Use either the standard system matrix for this purpose or
              ! the specific matrix for the coarse grid correction if this exists
              ! in the level info structure.
              if (p_rcurrentLevel%rsystemMatrixOptCorrection%NEQ .ne. 0) then
                call cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrixOptCorrection,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              else
                call cgcor_calcOptimalCorrection (p_rsubnode%rcoarseGridCorrection,&
                                          p_rcurrentLevel%rsystemMatrix,&
                                          p_rcurrentLevel%rsolutionVector,&
                                          p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector,&
                                          p_rcurrentLevel%rcgcorrTempVector,&
                                          p_RfilterChain,&
                                          dstep)
              end if
              
              ! Perform the coarse grid correction by adding the coarse grid
              ! solution (with the calculated step-length parameter) to
              ! the current solution.
              
              call lsysbl_vectorLinearComb (p_rcurrentLevel%rtempVector,&
                                            p_rcurrentLevel%rsolutionVector,&
                                            dstep,1.0_DP)
                            
              ! Extended output
              if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                  (mod(ite,niteResOutput).eq.0)) then
                  
                call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                        p_rcurrentLevel%rtempVector)
                call lsysbl_blockMatVec (&
                    p_rcurrentLevel%rsystemMatrix, &
                    p_rcurrentLevel%rsolutionVector,&
                    p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                    rsolverNode%iresNorm)
                if (.not.((dres .ge. 1E-99_DP) .and. &
                          (dres .le. 1E99_DP))) dres = 0.0_DP
                          
                call output_line ('Multigrid: Level '//trim(sys_siL(ilev,5))//&
                    ' after c.g.corr.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
              end if
                                            
              ! Perform the post-smoothing with the current solution vector
              if (associated(p_rcurrentLevel%p_rpostSmoother)) then
                call linsol_smoothCorrection (p_rcurrentLevel%p_rpostSmoother,&
                          p_rcurrentLevel%rsystemMatrix,&
                          p_rcurrentLevel%rsolutionVector,&
                          p_rcurrentLevel%rrhsVector,&
                          p_rcurrentLevel%rtempVector,p_RfilterChain)

                ! Extended output
                if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                    (mod(ite,niteResOutput).eq.0)) then
                    
                  call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                          p_rcurrentLevel%rtempVector)
                  call lsysbl_blockMatVec (&
                      p_rcurrentLevel%rsystemMatrix, &
                      p_rcurrentLevel%rsolutionVector,&
                      p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                  dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  call output_line ('Multigrid: Level '//trim(sys_siL(ilev,5))//&
                      ' after postsm.:    !!RES!! = '//trim(sys_sdEL(dres,15)) )
                end if
                                              
              end if

              ! Update the iteration counter(s) for realising the MG-cycle(s).
              ! Then either repeat this loop to perform the next prolongation or
              ! repeat the cycleloop to do perform a next MG sweep on the current
              ! level.
              !
              ! Here icycle defines how the cycle-counters are updated.
              ! For a W-cycle the cycle counter is resetted to 2 if the sweep is
              ! fulfilled on the current level, for F-cycle it is set to 1 to not
              ! perform more that 1 cycle on the current level anymore.
              
              p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncyclesRemaining-1
              if (p_rcurrentLevel%ncyclesRemaining .le. 0) then
                if (p_rsubnode%icycle .eq. 0) then
                  p_rcurrentLevel%ncyclesRemaining = 1
                else
                  ! Cycle finished. Reset counter for next cycle.
                  p_rcurrentLevel%ncyclesRemaining = p_rcurrentLevel%ncycles

                  if ((rsolverNode%p_rsubnodeMultigrid2%depsRelCycle .ne. 1E99_DP) .and. &
                      (ilev .lt. NLMAX)) then
                      
                    ! Adaptive cycles activated.
                    !
                    ! We are on a level < nlmax.
                    ! At first, calculate the residuum on that level.
                    call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,&
                                            p_rcurrentLevel%rtempVector)
                    call lsysbl_blockMatVec (&
                        p_rcurrentLevel%rsystemMatrix, &
                        p_rcurrentLevel%rsolutionVector,&
                        p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)

                    dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,&
                        rsolverNode%iresNorm)
                    if (.not.((dres .ge. 1E-99_DP) .and. &
                              (dres .le. 1E99_DP))) dres = 0.0_DP
                              
                    ! Compare it with the initial residuum. If it is not small enough
                    ! and if we have not reached the maximum number of cycles,
                    ! repeat the complete cycle.
                    if ( ((rsolverNode%p_rsubnodeMultigrid2%nmaxAdaptiveCycles .le. -1) &
                          .or. &
                          (p_rcurrentLevel%icycleCount .lt. &
                           rsolverNode%p_rsubnodeMultigrid2%nmaxAdaptiveCycles)) &
                        .and. &
                        (dres .gt. rsolverNode%p_rsubnodeMultigrid2%depsRelCycle * &
                                    p_rcurrentLevel%dinitResCycle) ) then

                      if (rsolverNode%ioutputLevel .ge. 3) then
                        call output_line ( &
                          trim( &
                          sys_siL(p_rcurrentLevel%icycleCount,10)) &
                          //'''th repetition of cycle on level '// &
                          trim(sys_siL(ilev,5))//'.')
                      end if

                      p_rcurrentLevel%icycleCount = p_rcurrentLevel%icycleCount+1
                      cycle cycleloop
                    end if
                    
                    ! Otherwise: The cycle(s) is/are finished;
                    ! the END DO goes up one level.
                    
                  end if

                end if
              else

                ! Next cycle; go down starting from the current level
                cycle cycleloop

              end if
              
            end do ! ilev < nlmax
            
            ! We finally reached the maximum level. Quit the cycle loop
            ! to go on with the next iteration.
            exit cycleloop
          
          end do cycleloop
          
          ! We have (hopefully) successfully performed one MG-sweep, starting
          ! and ending on the finest level. As we are now on the finest level
          ! again, we can update our defect vector to test the current
          ! residuum...
          !
          ! Calculate the residuum and its norm.
          call lsysbl_copyVector (p_rcurrentLevel%rrhsVector,p_rcurrentLevel%rtempVector)
          call lsysbl_blockMatVec (&
                p_rcurrentLevel%rsystemMatrix, &
                p_rcurrentLevel%rsolutionVector,&
                p_rcurrentLevel%rtempVector, -1.0_DP,1.0_DP)
          if (bfilter) then
            ! Apply the filter chain to the vector
            call filter_applyFilterChainVec (p_rcurrentLevel%rtempVector, &
                                             p_RfilterChain)
          end if
          dres = lsysbl_vectorNorm (p_rcurrentLevel%rtempVector,rsolverNode%iresNorm)
          if (.not.((dres .ge. 1E-99_DP) .and. &
                    (dres .le. 1E99_DP))) dres = 0.0_DP
          
          ! Shift the queue with the last residuals and add the new
          ! residual to it.
          Dresqueue = eoshift(Dresqueue,1,dres)

          rsolverNode%dlastDefect = rsolverNode%dfinalDefect
          rsolverNode%dfinalDefect = dres
          
          ! Test if the iteration is diverged
          if (linsol_testDivergence(rsolverNode,dres)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('Multigrid: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
            call output_line ('Multigrid: Solution diverging!')
            rsolverNode%iresult = 1
            exit
          end if
       
          ! At least perform nminIterations iterations
          if (ite .ge. nminIterations) then
          
            ! Check if the iteration converged
            if (linsol_testConvergence(rsolverNode,dres)) exit
            
          end if

          ! print out the current residuum

          if ((rsolverNode%ioutputLevel .ge. 2) .and. &
              (mod(ite,niteResOutput).eq.0)) then
            call output_line ('Multigrid: Iteration '// &
                trim(sys_siL(ite,10))//',  !!RES!! = '//&
                trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          end if
          
        end do  ! ite
        
        ! Finish - either with an error or if converged.
        ! Print the last residuum.

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (ite .ge. 1) .and. (ite .le. rsolverNode%nmaxIterations) .and. &
            (rsolverNode%iresult .ge. 0)) then
          call output_line ('Multigrid: Iteration '// &
              trim(sys_siL(ite,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if
        
        ! Set ite to NIT to prevent printing of "NIT+1" of the loop was
        ! completed

        if (ite .gt. rsolverNode%nmaxIterations) then
          ! Warning if we did not reach the convergence criterion.
          ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
          if ((rsolverNode%ioutputLevel .ge. 0) .and. &
              (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
            call output_line ('Multigrid: Accuracy warning: '//&
                'Solver did not reach the convergence criterion')
          end if

          ite = rsolverNode%nmaxIterations
        end if

        ! Final number of iterations
        rsolverNode%iiterations = ite
      
        ! Release the temporary copy of the solution vector.
        call lsysbl_releaseVector (p_rcurrentLevel%rsolutionVector)
      
      end if

      ! Finally, we look at some statistics:
      !
      ! Do not calculate anything if the final residuum is out of bounds -
      ! would result in NaN`s,...
        
      if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
        ! Calculate asymptotic convergence rate
      
        if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
          i = min(rsolverNode%iiterations,niteAsymptoticCVR)
          rsolverNode%dasymptoticConvergenceRate = &
            (rsolverNode%dfinalDefect / &
              dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
        end if

        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting multigrid convergence rate stays zero.
        ! In the other case the multigrid convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :

        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/real(rsolverNode%iiterations,DP))
        end if
      end if
      
    end if

    ! As the solution vector on the finest level shared its memory with rd,
    ! we just calculated the new correction vector!
      
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
      ! Scale the defect by the damping parameter in the solver structure.
      call lsysbl_scaleVector (rd,rsolverNode%domega)
      
      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('Multigrid statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'MG: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if

    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    end if
  
  end subroutine

! *****************************************************************************
! Routines for the Block-Jacobi solver (for DG discretisations)
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initBlockJac (p_rsolverNode)
  
!<description>
  ! Creates a t_linsolNode solver structure for the Block-Jacobi solver. The
  ! node can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  ! It is a special solver taking advantage of the block-structure stemming
  ! from a discontinuous Galerkin discretisation.

!</description>
  
!<input>
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_BLOCKJAC
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + &
                                LINSOL_ABIL_DIRECT

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatBlockJac (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  integer :: i,n
  type(t_matrixScalar), pointer :: p_rblock => null()

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Get the index of the top-level matrix
    n = ubound(Rmatrices,1)
    
    ! Let us make sure the matrix is not rectangular.
    if(Rmatrices(n)%nblocksPerRow .ne. Rmatrices(n)%nblocksPerCol) then
      
      ! Rectangular block matrix
      ccompatible = LINSOL_COMP_ERRMATRECT
      
    else
    
      ! Okay, let us loop through all main diagonal blocks.
      do i = 1, Rmatrices(n)%nblocksPerRow
      
        ! Get the diagonal block
        p_rblock => Rmatrices(n)%RmatrixBlock(i,i)
        
        ! Is the matrix empty?
        if(p_rblock%NEQ .eq. 0) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! Is the scaling factor zero?
        if(p_rblock%dscaleFactor .eq. 0.0_DP) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! What about the matrix type?
        if((p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
           (p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
          ! Matrix type not supported
          ccompatible = LINSOL_COMP_ERRMATTYPE
          exit
        end if
        
        ! This block is okay, continue with next one
      
      end do
    
    end if
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precBlockJac (rsolverNode,rd)
  
  use triangulation
  use spatialdiscretisation
  use element
  use dofmapping
  
!<description>
  ! Applies the Jacobi preconditioner <tex>$ D \in A $<tex> to the defect
  ! vector rd and solves <tex>$ Dd_{new} = d $<tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Jacobi solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: iblock
    
    type(t_matrixScalar), pointer :: p_rmatrix
    integer, dimension(:), pointer :: p_KCOL, p_KLD
    integer, dimension(:), pointer :: p_IelementList
    integer :: NEL, iel, igel, iFESpace, indof, i, ig, j, jg,iloc
    integer(I32) :: celement
    real(DP) :: dlocOmega
    real(DP), dimension(:), pointer :: p_Dvector, p_Dmatrix
    real(dp), dimension(:,:), allocatable :: DlocMat
    real(dp), dimension(:), allocatable :: DlocVec, DlocSol
    integer, dimension(:), allocatable :: IdofGlob
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    type(t_elementDistribution), pointer :: p_relemDist
    logical :: bsuccess
    
      ! Status reset
      rsolverNode%iresult = 0

      ! Loop through all blocks. Each block corresponds to one
      ! diagonal block in the matrix.
      do iblock = 1, rd%nblocks
        
        ! Get the matrix
        p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
        
        ! Now we have to make some decisions. At first, which matrix
        ! structure do we have?
        select case (p_rmatrix%cmatrixFormat)
        case (LSYSSC_MATRIX9)
        
          ! Get the omega parameter for the matrix. Do not forget the scaling
          ! factor in the matrix!
          dlocomega = rsolverNode%domega / p_rmatrix%dScaleFactor

          ! Now, which data format do we have? Single or double?
          select case (p_rmatrix%cdataType)
          case (ST_DOUBLE)
            ! Get pointers to the matrix data
            call lsyssc_getbase_Kcol (p_rmatrix,p_KCOL)
            call lsyssc_getbase_Kld (p_rmatrix,p_KLD)
            call lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
            
            ! Take care of the accuracy of the vector
            select case (rd%RvectorBlock(iblock)%cdataType)
            case (ST_DOUBLE)
              ! Get the data array
              call lsyssc_getbase_double (rd%RvectorBlock(iblock),p_Dvector)
              
              ! Get pointers for quicker access
              p_rspatialDiscr => p_rmatrix%p_rspatialDiscrTest
              p_rtriangulation => p_rspatialDiscr%p_rtriangulation

              ! Get number of elements
              NEL = p_rtriangulation%NEL

              ! Loop over all element distributions
              do iFESpace = 1, p_rspatialDiscr%inumFESpaces

                ! Get the element distribution from the spatial discretisation
                p_relemDist => p_rspatialDiscr%RelementDistr(iFESpace)

                ! Get type of element in this distribution
                celement =  p_relemDist%celement

                ! Get number of local DOF
                indof = elem_igetNDofLoc(celement)
                
                ! Allocate space for global DOFs
                allocate(IdofGlob(indof))

                ! Get list of all elements in this distribution
                call storage_getbase_int(p_relemDist%h_IelementList, p_IelementList)

                ! Allocate local matrix and vector
                allocate(DlocMat(indof,indof), DlocVec(indof), DlocSol(indof))

                ! Loop over all elements in this element distribution
                do iel = 1, p_relemDist%NEL

                  ! Get global element number
                  igel = p_IelementList(iel)

                  ! Map local to global DOFs
                  call dof_locGlobMapping(p_rspatialDiscr, igel, IdofGlob)

                  !!! Get entries from the global into the local matrix
                  ! Loop over all local lines
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)

                    ! Loop over all local columns
                    do j = 1, indof
                      ! Get global columns
                      jg = IdofGlob(j)

                      do iloc = p_KLD(ig), p_KLD(ig+1)-1
                        if (jg.eq.p_KCOL(iloc)) exit
                      end do

                      DlocMat(i,j) = p_Dmatrix(iloc)

                    end do
                  end do

                  !!! Now get local vector
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)
                     
                    ! Get global entry in local vector
                    DlocVec(i) = p_Dvector(ig)
                  end do

                  !!! Now solve local linear system
                  call mprim_invertMatrixDble(DlocMat, DlocVec, DlocSol, &
                                              indof, 2, bsuccess)

                  ! Test if local system could be solved
                  if (.not.bsuccess) then
                    call output_line ( &
                      'Local matricx singular.', &
                      OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockJac')
                      do i = 1,size(DlocMat,1)
                        write(*,*) DlocMat(i,:)
                      end do
                    call sys_halt()
                  end if

                  !!! Distribute local solution to global solution
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)
                     
                    ! Get local entry in global vector
                    p_Dvector(ig) = dlocOmega * DlocSol(i)
                  end do
                  
                end do ! iel
                  
                ! Deallocate local matrix and vector
                deallocate(DlocMat, DlocVec, DlocSol)
                
                ! Deallocate space for global DOFs
                deallocate(IdofGlob)
                  
              end do ! iFESpace
              
            case default
              call output_line ('Unsupported vector precision.', &
                                OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockJac')
              call sys_halt()
            end select
         
         
          case default
            call output_line ('Unsupported matrix precision.', &
                              OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockJac')
            call sys_halt()
          end select
        
        case default
          call output_line ('Unsupported matrix format.', &
                            OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockJac')
          call sys_halt()
        end select
        
      end do
  
  end subroutine

! *****************************************************************************
! Routines for the spectral deflation preconditioner
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initSpecDefl (p_rsolverNode, drelax, &
          p_rpreconditioner,Rfilter)
  
!<description>
  ! Creates a t_linsolNode solver structure for the spectral deflation
  ! preconditioner.
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
!</description>
  
!<input>
  ! Relaxation parameter for the deflation process.
  real(DP), optional :: drelax
  
  ! OPTIONAL: A pointer to the solver structure of a solver that should be
  ! used for preconditioning. If not given or set to NULL(), no preconditioning
  ! will be used.
  type(t_linsolNode), pointer, optional   :: p_rpreconditioner

  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not present, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_SPECDEFL
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_DIRECT
    
    ! Prepare the subnode
    allocate (p_rsolverNode%p_rsubnodeSpecDefl)
    p_rsolverNode%p_rsubnodeSpecDefl%drelax = drelax
    if (present(p_rpreconditioner)) then
        p_rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner => p_rpreconditioner
    end if
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubnodeSpecDefl%p_Rfilter => Rfilter
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_setMatrixSpecDefl (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Pass the matrix to the preconditioner
    if (associated(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner)) then
      call linsol_setMatrices(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner,&
          Rmatrices)
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initStructureSpecDefl (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the initStructure for the subsolver if present
    if (associated(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner)) then
      call linsol_initStructure(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner,&
          ierror,isolverSubgroup)
    end if
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Allocate memory for the temp vector
    call lsysbl_createVecBlockIndMat (rsolverNode%rsystemMatrix, &
        rsolverNode%p_rsubnodeSpecDefl%rvectorTemp,.false.,.false.,&
        rsolverNode%cdefaultDataType)
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneStructureSpecDefl (rsolverNode,isolverSubgroup)
  
!<description>
  ! Releases the memory.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in) :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the initStructure for the subsolver if present
    if (associated(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner)) then
      call linsol_doneStructure(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner,&
          isolverSubgroup)
    end if

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
    ! Release temp vector if still there
    if (rsolverNode%p_rsubnodeSpecDefl%rvectorTemp%neq .ne. 0) then
      call lsysbl_releaseVector(rsolverNode%p_rsubnodeSpecDefl%rvectorTemp)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initDataSpecDefl (rsolverNode, ierror, isolverSubgroup)
  
!<description>
  ! Performs final preparation of the spectral deflation solver.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the VANKA solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  real(DP) :: dmin,dmax

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the initData for the subsolver if present
    if (associated(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner)) then
      call linsol_initData(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner,&
          ierror,isolverSubgroup)
    end if

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Figure out an approximation for the maximum eigenvalue
    call lsyssc_calcGerschgorin (rsolverNode%rsystemMatrix%RmatrixBlock(1,1),&
            dmin,dmax)
            
    rsolverNode%p_rsubnodeSpecDefl%dmaxEigenvalue = dmax

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneDataSpecDefl (rsolverNode,isolverSubgroup)
  
!<description>
  ! Releases the memory.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! Call the doneData for the subsolver if present
    if (associated(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner)) then
      call linsol_doneData(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner,&
          isolverSubgroup)
    end if

    ! If isubgroup does not coincide with isolverSubgroup from the solver
    ! structure, skip the rest here.
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneSpecDefl (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it is also released from the heap by calling
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the solver.
  type(t_linsolNode), pointer         :: rsolverNode
!</input>
  
!</subroutine>

    ! Release memory if still associated
    call linsol_doneDataSpecDefl (rsolverNode, rsolverNode%isolverSubgroup)
    call linsol_doneStructureSpecDefl (rsolverNode, rsolverNode%isolverSubgroup)
    
    ! Release the subsolvers
    call linsol_releaseSolver(rsolverNode%p_rsubnodeSpecDefl%p_rpreconditioner)
    
    ! Release the subnode itself
    deallocate(rsolverNode%p_rsubnodeSpecDefl)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precSpecDefl (rsolverNode,rd)
  
!<description>
  ! Applies the spectral deflation preconditioner to the defect
  ! vector rd.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

    ! The spectral deflation algorithm applies the following
    ! modifications:
    !
    !    d  :=  d - (A d  -  drelax lambda_max d)
    !        =  d + (drelax lambda_max d  -  A d)
    !
    ! If a preconditioner P is specified, it is applied to the defect,
    ! which results in the fomula
    !
    !    d  :=  d - P(A d  -  drelax lambda_max d)

    ! local variables
    type (t_linsolSubnodeSpecDefl), pointer :: p_rsubnode
    type(t_vectorBlock), pointer :: p_rtempVec
    
    ! Temp vector
    p_rsubnode => rsolverNode%p_rsubnodeSpecDefl
    p_rtempVec => p_rsubnode%rvectorTemp
    
    ! Calculate (A d  -  drelax lambda_max d)
    ! in the temp vector
    call lsysbl_copyVector (rd,p_rtempVec)
    call lsysbl_blockMatVec (rsolverNode%rsystemMatrix, rd, p_rtempVec, &
        1.0_DP, -p_rsubnode%drelax*p_rsubnode%dmaxEigenvalue)

    if (associated(p_rsubnode%p_Rfilter)) then
      call filter_applyFilterChainVec (p_rtempVec, p_rsubnode%p_Rfilter)
    end if
    
    ! Preconditioning
    if (associated(p_rsubnode%p_rpreconditioner)) then
      call linsol_precondDefect (p_rsubnode%p_rpreconditioner,p_rtempVec)
    end if

    ! Correction
    call lsysbl_vectorLinearComb (p_rtempVec,rd,-1.0_DP,1.0_DP)
      
  end subroutine
  
  ! ***************************************************************************
  ! Deflated MG GMRES algorithm
  ! ***************************************************************************

!<function>
  
  integer function linsol_getDeflGMRESLevelCount (rsolverNode)
                    
!<description>
  ! Returns the number of levels currently attached to the deflated GMRES solver.
!</description>
  
!<input>
  ! The solver structure of the deflated GMRES solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</input>

!<result>
  ! The number of levels in the MG solver.
!</result>
  
!</function>

    linsol_getDeflGMRESLevelCount = rsolverNode%p_rsubnodeDeflGMRES%nlevels

  end function

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_getDeflGMRESLevel (rsolverNode,ilevel,p_rlevelInfo)
                    
!<description>
  ! Searches inside of the deflated GMRES structure for level ilevel and returns
  ! a pointer to the corresponding p_rlevelInfo structure-
  ! If the level does not exist, an error is thrown.
!</description>
  
!<inputoutput>
  ! The solver structure of the deflated GMRES solver
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!<input>
  ! Number of the level to fetch.
  integer, intent(in) :: ilevel
!</input>
  
!<output>
  ! A pointer to the corresponding t_levelInfo structure.
  type(t_linsolDeflGMRESLevelInfo), pointer     :: p_rlevelInfo
!</output>
  
!</subroutine>

    nullify(p_rlevelInfo)

    ! Do we have the level?
    if ((ilevel .gt. 0) .and. &
        (ilevel .le. rsolverNode%p_rsubnodeDeflGMRES%nlevels)) then
      ! Get it.
      p_rlevelInfo => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)
    else
      call output_line('Level out of range!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'linsol_getDeflGMRESLevel')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneDeflGMRESLevel (rlevelInfo)
                    
!<description>
  ! Cleans up the deflated GMRES level structure rlevelInfo. All memory allocated
  ! in this structure is released.
!</description>
  
!<inputoutput>
  ! The t_levelInfo structure to clean up.
  type(t_linsolDeflGMRESLevelInfo), intent(inout)     :: rlevelInfo
!</inputoutput>

!</subroutine>
  
    ! Release the temporary vectors of this level.
    ! The vector may not exist - if the application has not called
    ! initStructure!
    if (rlevelInfo%rrhsVector%NEQ .ne. 0) &
      call lsysbl_releaseVector(rlevelInfo%rrhsVector)

    if (rlevelInfo%rtempVector%NEQ .ne. 0) &
      call lsysbl_releaseVector(rlevelInfo%rtempVector)

    if (rlevelInfo%rsolutionVector%NEQ .ne. 0) &
      call lsysbl_releaseVector(rlevelInfo%rsolutionVector)

    ! Release sub-solvers if associated.
    if (associated(rlevelInfo%p_rpreconditioner)) then
      call linsol_releaseSolver(rlevelInfo%p_rpreconditioner)
    end if
    
    ! Clean up the associated matrix if there is one.
    ! Of course, the memory of the matrices will not be deallocated
    ! because the matrices belong to the application, not to the solver!
    call lsysbl_releaseMatrix(rlevelInfo%rsystemMatrix)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_cleanDeflGMRESLevels (rsolverNode)
                    
!<description>
  ! This routine removes all level information from the MG solver and releases
  ! all attached solver nodes on all levels (smoothers, coarse grid solvers,
  ! ...).
  !
  ! The level info array is not released.
!</description>
  
!<inputoutput>
  ! The solver structure of the deflated GMRES solver.
  type(t_linsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

  integer :: ilevel

    ! Make sure the solver node is configured for deflated GMRES
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DeflGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_cleanDeflGMRESLevels")
      call sys_halt()
    end if

    ! Remove all the levels
    do ilevel = 1,rsolverNode%p_rsubnodeDeflGMRES%nlevels
      call linsol_doneDeflGMRESLevel (&
          rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel))
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initPrjDGMLvByPrj (rlevelInfo,rinterlevelProjection)
  
!<description>
  ! Initialises the interlevel projection structure for one level.
!</description>
  
!<inputoutput>
  ! The t_levelInfo structure to set up.
  type(t_linsolDeflGMRESLevelInfo), intent(inout)     :: rlevelInfo
  
  ! User-specified interlevel projection structure.
  ! Deflated GMRES will use the given projection structure for doing prolongation/
  ! restriction.
  type(t_interlevelProjectionBlock), intent(in), target  :: rinterlevelProjection
!</inputoutput>

!</subroutine>

    ! Check if the pointer to the projection structure is already assigned.
    ! If that is the case, release it.
    if (rlevelInfo%bstandardProjection) then
      ! Release the old standard projection structure.
      if (associated(rlevelInfo%p_rprojection)) then
        call mlprj_doneProjection (rlevelInfo%p_rprojection)
        deallocate (rlevelInfo%p_rprojection)
      end if
    end if
    
    ! Remember the user specified projection
    rlevelInfo%p_rprojection => rinterlevelProjection
    rlevelInfo%bstandardProjection = .false.

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initPrjDGMLvByDiscr (rlevelInfo,rdiscrCoarse,rdiscrFine)
  
!<description>
  ! Initialises the interlevel projection structure for one level
  ! based on coarse/fine grid discretisation structures.
!</description>
  
!<inputoutput>
  ! The t_levelInfo structure to set up.
  type(t_linsolDeflGMRESLevelInfo), intent(inout)     :: rlevelInfo
!</inputoutput>

!</input>
  ! Discretisation structure of the coarse level.
  type(t_blockDiscretisation), intent(in), target   :: rdiscrCoarse

  ! Discretisation structure of the fine level.
  type(t_blockDiscretisation), intent(in), target   :: rdiscrFine
!</input>

!</subroutine>

    ! Check if the pointer to the projection structure is already assigned.
    ! If that is the case, release it.
    if (rlevelInfo%bstandardProjection) then
      ! Release the old standard projection structure.
      if (associated(rlevelInfo%p_rprojection)) then
        call mlprj_doneProjection (rlevelInfo%p_rprojection)
      else
        ! Allocate a new projection structure
        allocate(rlevelInfo%p_rprojection)
      end if
    else
      ! Allocate a new projection structure
      rlevelInfo%bstandardProjection = .true.
      allocate(rlevelInfo%p_rprojection)
    end if
    
    ! Initialise the projection structure with standard parameters.
    call mlprj_initProjectionDiscr (rlevelInfo%p_rprojection,rdiscrCoarse)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initProjDeflGMRES (rsolverNode)
  
!<description>
  ! Initialises the interlevel projection structure for all levels
  ! using the standard parameters for prolongation/restriction.
  !
  ! The routine can be called only after the matrices are attached to
  ! the solver with the setMatrices routine. It will automatically detect
  ! if the caller already specified user defined projection structures
  ! by a call to linsol_initProjDeflGMRESLvByPrj().
  !
  ! The routine is called during initStructure but can also be called
  ! manually in advance after setMatrices() if the caller wants to
  ! alter the projection structures manually.
!</description>
  
!<inputoutput>
  ! A t_linsolNode structure. defining the solver.
  type(t_linsolNode), intent(inout)             :: rsolverNode
!</inputoutput>

!</subroutine>

    integer :: ilevel, nlmax
    type(t_blockDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine
    
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels
    
    ! Loop through all levels. Whereever the projection structure
    ! is missing, create it.
    do ilevel = 2,nlmax
      if (.not. associated( &
          rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)% &
          p_rprojection)) then
        
        ! Initialise the projection structure.
        p_rdiscrCoarse => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel-1)%rsystemmatrix%&
          p_rblockDiscrTrial
        p_rdiscrFine => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)%rsystemmatrix%&
          p_rblockDiscrTrial
        call linsol_initPrjDGMLvByDiscr (rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel),&
            p_rdiscrCoarse,p_rdiscrFine)
        
      end if
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneProjDeflGMRES (rsolverNode)
  
!<description>
  ! Releases all automatically allocated interlevel projection structures
  ! for all levels in the MG solver.
!</description>
  
!<inputoutput>
  ! A  t_linsolNode structure. defining the solver.
  type(t_linsolNode), intent(inout)             :: rsolverNode
!</inputoutput>

!</subroutine>

    integer :: ilevel, nlmax
    
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels
    
    ! Loop through all levels. Whereever the projection structure
    ! is missing, create it.
    do ilevel = 2,nlmax
      if (associated(rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)% &
          p_rprojection)) then
        if (rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)% &
            bstandardProjection) then
          call mlprj_doneProjection (rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)%&
              p_rprojection)
          deallocate (rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)%p_rprojection)
        else
          ! Just nullify the pointer, the structure does not belong to us.
          nullify (rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)%p_rprojection)
        end if
        
        rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)% &
            bstandardProjection = .false.
        
      end if
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_initDeflGMRES (p_rsolverNode,nlevels,btwiceGS,dpseudoResScale,&
      Rfilter)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the deflates GMRES solver. The
  ! node can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  ! Before the solver can used for solving the problem, the caller must add
  ! information about all levels (matrices,...) to the solver.
  ! This can be done by linsol_setMultigridLevel.
  
!</description>
  
!<input>

  ! Number of levels supported by this solver.
  integer, intent(in) :: nlevels
  
  ! Optional: A pointer to a filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not given or set to NULL(), no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), target, optional   :: Rfilter
  
  ! OPTIONAL: A boolean which specifies whether we should apply the Gram-Schmidt
  ! process once or twice per GMRES iteration.
  ! Since the Gram-Schmidt algorithm suffers under numerical instability, we
  ! have a chance to improve its stability by applying the Gram-Schmidt process
  ! twice. If btwiceGS is given and it is .TRUE., then the Gram-Schmidt process
  ! will be applied twice per GMRES iteration, otherwise only once.
  logical, optional :: btwiceGS

  ! OPTIONAL: Scale factor for the pseudo-residual-norms to test against the
  ! stopping criterion in the inner GMRES loop. Since the pseudo-residual-norms
  ! may be highly inaccurate, this scale factor can be used to scale the
  ! pseudo-residual-norms before testing against the stopping criterion.
  ! It is recommended to use a scale factor > 1, a factor of 2 should be usually
  ! okay. If dpseudoresscale is not given or set to < 1, it is automatically
  ! set to 2.
  ! Setting dpseudoresscale to 0 will disable the pseudo-resirudal norm check.
  real(DP), optional :: dpseudoResScale
!</input>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer             :: p_rsolverNode
   
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_DEFLGMRES
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR     + LINSOL_ABIL_BLOCK        + &
                                LINSOL_ABIL_MULTILEVEL + LINSOL_ABIL_CHECKDEF     + &
                                LINSOL_ABIL_USESUBSOLVER + &
                                LINSOL_ABIL_USEFILTER
    
    ! Allocate the subnode for Multigrid.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeDeflGMRES)
    
    ! Allocate level information data
    allocate(p_rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(nlevels))
    p_rsolverNode%p_rsubnodeDeflGMRES%nlevels = nlevels
    
    ! Save if the Gram-Schmidt process should be applied twice or not.
    if (present(btwiceGS)) then
      p_rsolverNode%p_rsubNodeDeflGMRES%btwiceGS = btwiceGS
    else
      p_rsolverNode%p_rsubNodeDeflGMRES%btwiceGS = LINSOL_GMRES_DEF_TWICE_GS
    end if
    
    ! Save the pseudo-residual-norm scale factor if given, or set
    ! to default scale factor if not given or invalid (i.e. < 1)
    if (present(dpseudoResScale)) then
      if (dpseudoResScale .ge. 0.0_DP) then
        p_rsolverNode%p_rsubNodeDeflGMRES%dpseudoResScale = dpseudoResScale
      else
        p_rsolverNode%p_rsubNodeDeflGMRES%dpseudoResScale = 0.0_DP
      endif
    else
      p_rsolverNode%p_rsubNodeDeflGMRES%dpseudoResScale = 0.0_DP
    endif

    ! Attach the filter if given.
    if (present(Rfilter)) then
      p_rsolverNode%p_rsubnodeDeflGMRES%p_RfilterChain => Rfilter
    end if
  
  end subroutine
  
! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_setMatrixDeflGMRES (rsolverNode,Rmatrices)
  
!<description>
  
  ! Assigns the system matrix rsystemMatrix to the Multigrid solver on
  ! all levels.
  
!</description>
  
!<input>
  
  ! An array of system matrices on all levels.
  ! Each level in multigrid is initialised separately.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices

!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!</subroutine>

  ! local variables
  type(t_linsolDeflGMRESLevelInfo), pointer :: p_rcurrentLevel
  integer                       :: ilevel,nlmax
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DEFLGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_setMatrixDeflGMRES")
      call sys_halt()
    end if
    
    ! Make sure we have the right amount of matrices
    if (size(Rmatrices) .ne. rsolverNode%p_rsubnodeDeflGMRES%nlevels) then
      call output_line ("Wrong number of matrices!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_setMatrixDeflGMRES")
      call sys_halt()
    end if

    ! Check for every level, if there is a preconditioner structure 
    ! attached. If yes, attach the matrix to it.
    ! For this purpose, call the general initialisation routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)

      if (associated(p_rcurrentLevel%p_rpreconditioner)) then
        call linsol_setMatrices (p_rcurrentLevel%p_rpreconditioner, &
                                  Rmatrices(1:ilevel))
      end if
      
      ! Write a link to the matrix into the level-structure.
      call lsysbl_duplicateMatrix (Rmatrices(ilevel), &
          p_rcurrentLevel%rsystemMatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_matCompatDeflGMRES (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  ! More precisely, ccompatible returns always the status of the highest
  ! level where there is an error -- or LINSOL_COMP_OK if there is no error.
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  ! local variables
  type(t_linsolDeflGMRESLevelInfo), pointer :: p_rcurrentLevel
  integer                       :: ilevel,ccompatLevel
    
    ! Normally, we can handle the matrix.
    ccompatible = LINSOL_COMP_OK
    
    ! Reset all compatibility flags
    if (present(CcompatibleDetail)) CcompatibleDetail (:) = ccompatible

    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DEFLGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_matCompatDeflGMRES")
      call sys_halt()
    end if
    
    ! Make sure we have the right amount of matrices
    if (size(Rmatrices) .ne. rsolverNode%p_rsubnodeDeflGMRES%nlevels) then
      call output_line ("Wrong number of matrices!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_matCompatDeflGMRES")
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached. If yes, it must check
    ! the matrices!
    ! For this purpose, call the general check routine, but
    ! pass only that part of the Rmatrices array that belongs
    ! to the range of levels between the coarse grid and the current grid.
    
    do ilevel = 1,rsolverNode%p_rsubnodeDeflGMRES%nlevels
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)
    
      ! Compatibility flag for that level
      ccompatLevel = LINSOL_COMP_OK
    
      ! Check presmoother
      if (associated(p_rcurrentLevel%p_rpreconditioner)) then
      
        if (present(CcompatibleDetail)) then
          call linsol_matricesCompatible (p_rcurrentLevel%p_rpreconditioner, &
              Rmatrices(1:ilevel),ccompatLevel,CcompatibleDetail(1:ilevel))
        end if
      
        if (ccompatLevel .ne. LINSOL_COMP_OK) then
          ccompatible = ccompatLevel
          cycle
        end if
        
      end if
      
    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initStructureDeflGMRES (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  type(t_linsolDeflGMRESLevelInfo), pointer :: p_rcurrentLevel,p_rprevLevel
  integer :: isubgroup,nlmax,ilevel
  integer :: imemmax,i
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_vectorBlock), pointer :: p_rtemplVect
  integer :: idim
  integer , dimension(2) :: idim2

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.

    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DEFLGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_initStructureDeflGMRES")
      call sys_halt()
    end if
    
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)

      if (associated(p_rcurrentLevel%p_rpreconditioner)) then
        call linsol_initStructure(p_rcurrentLevel%p_rpreconditioner,ierror,isubgroup)
      end if
      
      ! Cancel in case of an error
      if (ierror .ne. LINSOL_ERR_NOERROR) return
      
    end do
    
    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Initialise the interlevel projection as far as it is not already initialised
    call linsol_initProjDeflGMRES (rsolverNode)

    imemmax = 0

    ! Multigrid needs temporary vectors for the RHS, solution and general
    ! data. Memory for that is allocated on all levels for temporary, RHS and
    ! solution vectors.
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)

      p_rmatrix => p_rcurrentLevel%rsystemMatrix
      nullify(p_rtemplVect)

      call lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rsolutionVector,.false.,.false.,&
            rsolverNode%cdefaultDataType)
      p_rtemplVect => p_rcurrentLevel%rsolutionVector
      
      call lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rrhsVector,.false.,.false.,&
            rsolverNode%cdefaultDataType)
            
      call lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rtempVector,.false.,.false.,&
            rsolverNode%cdefaultDataType)

      call lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%rdefTempVector,.false.,.false.,&
            rsolverNode%cdefaultDataType)

      call lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%ractualRHS,.false.,.false.,&
            rsolverNode%cdefaultDataType)

      ! Calculate the memory that is probably necessary for resorting
      ! vectors during the iteration.
      if (associated(p_rtemplVect)) then
        imemmax = max(imemmax,p_rtemplVect%NEQ)
      end if
         
      ! Ask the prolongation/restriction routines how much memory is necessary
      ! for prolongation/restriction on that level.
      if (ilevel .gt. 1) then ! otherwise coarse grid
        p_rprevLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel-1)
     
        ! Calculate the memory that is necessary for prolongation/restriction -
        ! in case there is temporary memory needed.
        ! The system matrix on the fine/coarse grid specifies the discretisation.
        imemmax = max(imemmax,mlprj_getTempMemoryMat ( &
                      p_rcurrentLevel%p_rprojection, &
                      p_rprevLevel%rsystemMatrix,&
                      p_rcurrentLevel%rsystemMatrix))
      end if

      ! All temporary vectors are marked as "unsorted". We can set the
      ! sorting flag directly here without using the resorting routines
      ! as the vectors have just been created.
      !
      ! Do not do anything to the RHS/temp vectors on the coarse grid
      ! and the solution vector on the fine grid -- they do not exist!
      if (ilevel .gt. 1) then
        p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rrhsVector%RvectorBlock%isortStrategy)

        p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rtempVector%RvectorBlock%isortStrategy)

        p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy = &
          -abs(p_rcurrentLevel%rsolutionVector%RvectorBlock%isortStrategy)
      end if
      
      ! We now need to check if the dimension of the Krylov subspace is
      ! positive. If it is not, we need to cancel the initialization here.
      if (p_rcurrentLevel%ikrylowDim .le. 0) then
        call output_line ('imension of Krylov subspace for GMRES(m) is <= 0 !', &
            OU_CLASS_ERROR, OU_MODE_STD, 'mysubroutine')
        ierror = LINSOL_ERR_INITERROR
        call sys_halt()
      end if
      
      ! Get the stuff out of our solver node
      idim = p_rcurrentLevel%ikrylowDim
      idim2(1) = idim
      idim2(2) = idim
  
      ! Now comes the intersting part - we need to allocate a lot of arrays
      ! and vectors for the Krylov subspace for GMRES here.
      
      ! Call our storage to allocate the 1D/2D arrays
      call storage_new('linsol_initStructureGMRES', 'Dh', idim2, &
          ST_DOUBLE, p_rcurrentLevel%hDh, ST_NEWBLOCK_NOINIT)
      call storage_new('linsol_initStructureGMRES', 'Dc', idim, &
          ST_DOUBLE, p_rcurrentLevel%hDs, ST_NEWBLOCK_NOINIT)
      call storage_new('linsol_initStructureGMRES', 'Ds', idim, &
          ST_DOUBLE, p_rcurrentLevel%hDc, ST_NEWBLOCK_NOINIT)
      call storage_new('linsol_initStructureGMRES', 'Dq', idim+1, &
          ST_DOUBLE,  p_rcurrentLevel%hDq, ST_NEWBLOCK_NOINIT)
      
      ! Get the pointers
      call storage_getbase_double2D(p_rcurrentLevel%hDh, p_rcurrentLevel%Dh)
      call storage_getbase_double(p_rcurrentLevel%hDc, p_rcurrentLevel%Dc)
      call storage_getbase_double(p_rcurrentLevel%hDs, p_rcurrentLevel%Ds)
      call storage_getbase_double(p_rcurrentLevel%hDq, p_rcurrentLevel%Dq)
      
      ! Allocate space for our auxiliary vectors
      allocate(p_rcurrentLevel%p_rv(idim+1))
      allocate(p_rcurrentLevel%p_rz(idim))
      
      ! Create them
      do i=1, idim+1
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%p_rv(i),.false.,.false.,rsolverNode%cdefaultDataType)
      end do
      do i=1, idim
        call lsysbl_createVecBlockIndMat (p_rmatrix, &
            p_rcurrentLevel%p_rz(i),.false.,.false.,rsolverNode%cdefaultDataType)
      end do
      
      ! Create an iteration vector x
      call lsysbl_createVecBlockIndMat (p_rmatrix, &
          p_rcurrentLevel%rx,.false.,.false.,rsolverNode%cdefaultDataType)
      
      ! And the next level...
    end do
    
    ! Do we have to allocate memory for prolongation/restriction/resorting?
    if (imemmax .gt. 0) then
      call lsyssc_createVector (rsolverNode%p_rsubnodeDeflGMRES%rprjTempVector,&
                                imemmax,.false.,rsolverNode%cdefaultDataType)
    else
      rsolverNode%p_rsubnodeDeflGMRES%rprjTempVector%NEQ = 0
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_initDataDeflGMRES (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  type(t_linsolDeflGMRESLevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup,ilevel
  
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup
    
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DEFLGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_initDataDeflGMRES")
      call sys_halt()
    end if

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    
    do ilevel = 1,rsolverNode%p_rsubnodeDeflGMRES%nlevels
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)

      if (associated(p_rcurrentLevel%p_rpreconditioner)) then
        call linsol_initData(p_rcurrentLevel%p_rpreconditioner,ierror,isubgroup)

        ! Cancel in case of an error
        if (ierror .ne. LINSOL_ERR_NOERROR) return
      end if
      
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataDeflGMRES (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  type(t_linsolDeflGMRESLevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup,ilevel
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DEFLGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_doneDataDeflGMRES")
      call sys_halt()
    end if

    ! Check for every level, if there is a preconditioner attached.
    
    do ilevel = rsolverNode%p_rsubnodeDeflGMRES%nlevels,1,-1
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)
      
      if (associated(p_rcurrentLevel%p_rpreconditioner)) then
        call linsol_doneData(p_rcurrentLevel%p_rpreconditioner,isubgroup)
      end if
      
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureDeflGMRES (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
   
!</inputoutput>
  
!<input>
  
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup

!</input>

!</subroutine>

  ! local variables
  type(t_linsolDeflGMRESLevelInfo), pointer :: p_rcurrentLevel
  integer :: isubgroup,ilevel,nlmax,idim,i
  
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DEFLGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_doneStructureDeflGMRES")
      call sys_halt()
    end if
    
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels

    ! Check for every level, if there is a presmoother, postsmoother or
    ! coarse grid solver structure attached.
    
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)
    
      ! Release the presmoother
      if (associated(p_rcurrentLevel%p_rpreconditioner)) then
        call linsol_doneStructure(p_rcurrentLevel%p_rpreconditioner,isubgroup)
      end if
      
    end do

    ! Cancel here, if we do not belong to the subgroup to be initialised
    if (isubgroup .ne. rsolverNode%isolverSubgroup) return

    ! Release temporary memory.
    do ilevel = 1,nlmax
    
      p_rcurrentLevel => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)
    
      if (p_rcurrentLevel%rtempVector%NEQ .ne. 0) then
        call lsysbl_releaseVector (p_rcurrentLevel%rtempVector)
      end if
      
      if (p_rcurrentLevel%rrhsVector%NEQ .ne. 0) then
        call lsysbl_releaseVector (p_rcurrentLevel%rrhsVector)
      end if

      if (p_rcurrentLevel%rsolutionVector%NEQ .ne. 0) then
        call lsysbl_releaseVector (p_rcurrentLevel%rsolutionVector)
      end if

      if (p_rcurrentLevel%rdefTempVector%NEQ .ne. 0) then
        call lsysbl_releaseVector (p_rcurrentLevel%rdefTempVector)
      end if

      if (p_rcurrentLevel%ractualRHS%NEQ .ne. 0) then
        call lsysbl_releaseVector (p_rcurrentLevel%ractualRHS)
      end if

      idim = p_rcurrentLevel%ikrylowDim
      
      ! Release auxiliary vectors
      if (associated(p_rcurrentLevel%p_rz)) then
        do i=1, idim
          call lsysbl_releaseVector(p_rcurrentLevel%p_rz(i))
        end do
        do i=1, idim+1
          call lsysbl_releaseVector(p_rcurrentLevel%p_rv(i))
        end do
      
        ! Destroy them
        deallocate(p_rcurrentLevel%p_rz)
        deallocate(p_rcurrentLevel%p_rv)
      endif
  
      ! Release iteration vector
      if (p_rcurrentLevel%rx%NEQ > 0) then
        call lsysbl_releaseVector(p_rcurrentLevel%rx)
      end if
      
      ! Destroy our 1D/2D arrays
      if (p_rcurrentLevel%hDq .ne. ST_NOHANDLE) then
        call storage_free(p_rcurrentLevel%hDq)
        call storage_free(p_rcurrentLevel%hDs)
        call storage_free(p_rcurrentLevel%hDc)
        call storage_free(p_rcurrentLevel%hDh)
      end if

    end do

    ! Check if we have to release our temporary vector for prolongation/
    ! restriction.
    ! Do we have to allocate memory for prolongation/restriction?
    if (rsolverNode%p_rsubnodeDeflGMRES%rprjTempVector%NEQ .gt. 0) then
      call lsyssc_releaseVector (rsolverNode%p_rsubnodeDeflGMRES%rprjTempVector)
    end if

    ! Release the interlevel projection structures as far as they were
    ! automatically created.
    call linsol_doneProjDeflGMRES (rsolverNode)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDeflGMRES (rsolverNode)
  
!<description>
  
  ! This routine releases all temporary memory for the multigrid solver from
  ! the heap. In particular, this releases the solver structures of all
  ! subsolvers (smoother, coarse grid solver).
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
  
!</description>
  
!<inputoutput>
  ! A pointer to a t_linsolNode structure of the Multigrid solver.
  type(t_linsolNode), intent(inout)                :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. LINSOL_ALG_DEFLGMRES) .or. &
        (.not. associated(rsolverNode%p_rsubnodeDeflGMRES))) then
      call output_line ("Deflated GMRES structure not initialised!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_doneDeflGMRES")
      call sys_halt()
    end if

    ! Release data and structure if this has not happened yet.
    call linsol_doneDataDeflGMRES(rsolverNode)
    call linsol_doneStructureDeflGMRES(rsolverNode)
    
    ! Remove all the levels
    call linsol_cleanDeflGMRESLevels (rsolverNode)
    deallocate(rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo)
    
    ! Release MG substructure, that is it.
    deallocate(rsolverNode%p_rsubnodeDeflGMRES)

  end subroutine
    
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_precDeflGMRES_recright (rsolverNode,ilevel,rx,rb)
  
!<description>
  ! Actual working routine for the deflated multilevel GMRES solver.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the deflaged GMRES solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
  
  ! Current level.
  integer, intent(in) :: ilevel
   
  ! Solution vector on level iilevel.
  ! Is replaced by a new approximation to the solution
  type(t_vectorBlock), intent(inout)        :: rx

  ! Right hand side vector on level ilevel
  type(t_vectorBlock), intent(inout)        :: rb
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_linsolDeflGMRESLevelInfo), pointer :: p_rlevelInfo
    type(t_linsolDeflGMRESLevelInfo), pointer :: p_rlevelInfoBelow
    type(t_linsolSubnodeDeflGMRES), pointer :: p_rsubnode
    type(t_linsolNode), pointer :: p_rprecSubnode
    integer :: nminIterations,nmaxIterations,nlmax
    
    real(DP) :: dalpha,dbeta,dres,dfr,dtmp,dpseudores,dprnsf
    integer :: ite,i,j,k,niteAsymptoticCVR
    
    ! Here come our 1D/2D arrays
    real(DP), dimension(:), pointer :: p_Dc, p_Ds, p_Dq
    real(DP), dimension(:,:), pointer :: p_Dh
  
    ! The queue saves the current residual and the two previous residuals.
    real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
    
    ! The system matrix
    type(t_matrixBlock), pointer :: p_rmatrix
    
    ! Minimum number of iterations, print-sequence for residuals, dimension
    ! of krylov subspace
    integer :: niteResOutput, idim, hDq
    
    ! Whether to filter/precondition, apply Gram-Schmidt twice
    logical bprec,bfilter,btwiceGS
    
    ! Pointers to temporary vectors - named for easier access
    type(t_vectorBlock), pointer :: p_rtempVector,p_rdefTempVector
    type(t_vectorBlock), dimension(:), pointer :: p_rv, p_rz
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! Get pointers to the current level and the GMRES parameters
    p_rsubnode => rsolverNode%p_rsubnodeDeflGMRES
    p_rlevelInfo => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)
    
    nullify(p_rlevelInfoBelow)
    if (ilevel .gt. 1) then
      ! The lower level
      p_rlevelInfoBelow => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel-1)
    end if
    
    ! On the maximum level, apply as many iterations as configured in the
    ! solver node. On the other levels, apply as many as configured in the
    ! level structure
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels

    if (ilevel .eq. nlmax) then
      nminIterations = max(rsolverNode%nminIterations,0)
      nmaxIterations = max(rsolverNode%nmaxIterations,0)

      ! Length of the queue of last residuals for the computation of
      ! the asymptotic convergence rate
      niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))
  
      ! Iteration when the residuum is printed:
      niteResOutput = max(1,rsolverNode%niteResOutput)
      
    else
      if (p_rlevelInfo%niterations .eq. -1) then
        nminIterations = max(p_rlevelInfo%ikrylowDim,0)
        nmaxIterations = max(p_rlevelInfo%ikrylowDim,0)
      else
        nminIterations = max(p_rlevelInfo%niterations,0)
        nmaxIterations = max(p_rlevelInfo%niterations,0)
      end if
      
      niteAsymptoticCVR = 0
      niteResOutput = 0
      
    end if

    ! Use preconditioning? Filtering?
    bprec = associated(p_rlevelInfo%p_rpreconditioner)
    bfilter = associated(p_rsubnode%p_RfilterChain)
    
    ! Apply Gram-Schmidt twice?
    btwiceGS = p_rsubnode%btwiceGS
    
    ! Get the dimension of Krylov subspace
    idim = p_rlevelInfo%ikrylowDim
    
    ! Get our pseudo-residual-norm scale factor
    dprnsf = p_rsubnode%dpseudoResScale
    
    ! Now we need to check the pseudo-residual scale factor.
    ! If it is non-positive, we set it to 0.
    if (dprnsf .le. 0.0_DP) then
      dprnsf = 0.0_DP
    end if

    ! Current matrix and other parameters
    p_rmatrix => p_rlevelInfo%rsystemMatrix
    
    ! Set pointers to the temporary vectors
    ! defect vectors
    p_rz => p_rlevelInfo%p_rz
    ! basis vectors
    p_rv => p_rlevelInfo%p_rv
    ! temp vector
    p_rtempVector => p_rlevelInfo%rtempVector
    
    ! temp vector for defect creation with preconditioning
    p_rdefTempVector => p_rlevelInfo%rdefTempVector

    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    do i=1, idim+1
      call lsysbl_assignDiscrIndirect (rb,p_rv(i))
    end do
    do i=1, idim
      call lsysbl_assignDiscrIndirect (rb,p_rz(i))
    end do
    
    ! Set pointers to the 1D/2D arrays
    p_Dh => p_rlevelInfo%Dh
    p_Dc => p_rlevelInfo%Dc
    p_Ds => p_rlevelInfo%Ds
    p_Dq => p_rlevelInfo%Dq
    
    ! And get the handle of Dq, since we need to call
    ! storage_clear during the iteration
    hDq = p_rlevelInfo%hDq
    
    ! All vectors share the same boundary conditions as rb!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    do i=1,idim
      call lsysbl_assignDiscrIndirect (rb, p_rz(i))
    end do
    do i=1,idim+1
      call lsysbl_assignDiscrIndirect (rb, p_rv(i))
    end do
    
    if (bprec) then
      p_rprecSubnode => p_rlevelInfo%p_rpreconditioner
    end if
    
    if (bfilter) then
      p_RfilterChain => p_rsubnode%p_RfilterChain
    end if
    
    ! Copy our RHS to p_rv(1). As the iteration vector is 0, this
    ! is also our initial defect.
    call lsysbl_copyVector(rb,p_rv(1))

    if (bfilter) then
      call filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
    end if
    
    ! Get the norm of the residuum.
    ! We need to calculate the norm twice, since we need one for
    ! the stopping criterion (selected by the user) and the
    ! euclidian norm for the internal GMRES algorithm
    dfr = lsysbl_vectorNorm (rb,rsolverNode%iresNorm)
    dres = lsysbl_vectorNorm (p_rv(1),LINALG_NORMEUCLID)
    if (.not.((dfr .ge. 1E-99_DP) .and. &
              (dfr .le. 1E99_DP))) dfr = 0.0_DP

    ! initialise the queue of the last residuals with RES
    Dresqueue = dfr

    if (ilevel .eq. nlmax) then
      ! Initialise starting residuum
      rsolverNode%dinitialDefect = dfr
      rsolverNode%dlastDefect = 0.0_DP

      ! Check if out initial defect is zero. This may happen if the filtering
      ! routine filters "everything out"!
      ! In that case we can directly stop our computation.
  
      if (rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero) then
        ! final defect is 0, as initialised in the output variable above
        call lsysbl_clearVector(rx)
        ite = 0
        rsolverNode%dfinalDefect = dfr
        return
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('GMRES('//&
             trim(sys_siL(idim,10))//'): Iteration '//&
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

    end if

    ite = 0
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! -= Outer Loop
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    do while (ite < nmaxIterations)
      
      ! Step O.1:
      ! Create elementary vector e_1 scaled by the euclid norm of the
      ! defect, i.e.: q = ( ||v(1)||_2, 0, ..., 0 )
      call storage_clear (hDq)
      p_Dq(1) = dres
      
      ! Step O.2:
      ! Now scale the defect by the inverse of its norm
      call lsysbl_scaleVector (p_rv(1), 1.0_DP / dres)
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= Inner Loop (GMRES iterations)
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      do i = 1, idim
        ! Another iteration
        ite = ite + 1
      
        ! Step I.1:
        ! Solve P * z(i) = v(i), where P is the preconditioner matrix
        call lsysbl_copyVector (p_rv(i), p_rz(i))
        
        ! Apply deflation preconditioning to z(i)
        if (associated(p_rlevelInfoBelow)) then
        
          ! Calculate z = (A P z  -  drelax lambda_max v(i))
          ! in the temp vector
          call lsysbl_copyVector (p_rv(i),p_rtempVector)
          call lsysbl_copyVector (p_rz(i),p_rdefTempVector)
          
          ! Preconditioning
          if (bprec) then
            call linsol_precondDefect (&
                p_rlevelInfo%p_rpreconditioner,p_rdefTempVector)
          end if

          call lsysbl_blockMatVec (p_rlevelInfo%rsystemMatrix, p_rdefTempVector, p_rtempVector, &
              1.0_DP, -p_rsubnode%drelax*p_rlevelInfo%dmaxEigenvalue)

          ! Apply the filter chain
          if (bfilter) then
            call filter_applyFilterChainVec (p_rtempVector, p_rsubnode%p_RfilterChain)
          end if
          
          ! Restriction to the RHS of the lower level
          call mlprj_performRestriction (p_rlevelInfo%p_rprojection,&
                p_rlevelInfoBelow%rrhsVector, &
                p_rtempVector, &
                p_rsubnode%rprjTempVector)

          ! Apply the filter chain to the new RHS
          if (bfilter) then
            call filter_applyFilterChainVec (p_rlevelInfoBelow%rrhsVector, p_rsubnode%p_RfilterChain)
          end if
          
          ! Preconditioning on the lower level
          call lsysbl_clearVector (p_rlevelInfoBelow%rsolutionVector)
          call linsol_precDeflGMRES_recright (rsolverNode,ilevel-1,&
              p_rlevelInfoBelow%rsolutionVector,p_rlevelInfoBelow%rrhsVector)

          ! Prolongation to our current level: t = I(...)
          call mlprj_performProlongation (p_rlevelInfo%p_rprojection,&
                p_rlevelInfoBelow%rsolutionVector, &
                p_rtempVector, &
                p_rsubnode%rprjTempVector)

          ! Correction: z(i) = v(i) - t
          call lsysbl_vectorLinearComb (p_rtempVector,p_rz(i),-1.0_DP,1.0_DP)
        
        end if
        
        ! Apply filter chain to z(i)
        if (bfilter) then
          call filter_applyFilterChainVec (p_rz(i), p_RfilterChain)
        end if
        
        ! Step I.2:
        ! Calculate v(i+1) = A * P * z(i)
        
        call lsysbl_copyVector (p_rz(i),p_rdefTempVector)

        if (bprec) then
          call linsol_precondDefect (&
              p_rlevelInfo%p_rpreconditioner,p_rdefTempVector)
        end if
        
        call lsysbl_blockMatVec (p_rmatrix, p_rdefTempVector, p_rv(i+1), 1.0_DP, 0.0_DP)

        if (bfilter) then
          call filter_applyFilterChainVec (p_rv(i+1), p_RfilterChain)
        end if
        
        
        ! Step I.3:
        ! Perfom Gram-Schmidt process (1st time)
        do k = 1, i
          p_Dh(k,i) = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
          call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -p_Dh(k,i), 1.0_DP)
        end do
        
        ! If the user wishes, perform Gram-Schmidt one more time
        ! to improve numerical stability.
        if (btwiceGS) then
          do k = 1, i
            dtmp = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
            p_Dh(k,i) = p_Dh(k,i) + dtmp
            call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -dtmp, 1.0_DP)
          end do
        end if
        
        ! Step I.4:
        ! Calculate alpha = ||v(i+1)||_2
        dalpha = lsysbl_vectorNorm (p_rv(i+1), LINALG_NORMEUCLID)
        
        ! Step I.5:
        ! Scale v(i+1) by the inverse of its euclid norm
        if (dalpha > SYS_EPSREAL_DP) then
          call lsysbl_scaleVector (p_rv(i+1), 1.0_DP / dalpha)
        else
          ! Well, let us just print a warning here...
          if(rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('GMRES('// trim(sys_siL(idim,10))//&
              '): Warning: !!v(i+1)!! < EPS !')
          end if
        end if
        
        
        ! Step I.6:
        ! Apply Givens rotations to the i-th column of h which
        ! renders the Householder matrix an upper triangular matrix
        do k = 1, i-1
          dtmp = p_Dh(k,i)
          p_Dh(k,i)   = p_Dc(k) * dtmp + p_Ds(k) * p_Dh(k+1,i)
          p_Dh(k+1,i) = p_Ds(k) * dtmp - p_Dc(k) * p_Dh(k+1,i)
        end do
        
        
        ! Step I.7:
        ! Calculate Beta
        ! Beta = (h(i,i)^2 + alpha^2) ^ (1/2)
        dbeta = sqrt(p_Dh(i,i)**2 + dalpha**2)
        if (dbeta < SYS_EPSREAL_DP) then
          dbeta = SYS_EPSREAL_DP
          if(rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('GMRES('// trim(sys_siL(idim,10))//&
              '): Warning: beta < EPS !')
          end if
        end if
        
        ! Step I.8:
        ! Calculate next plane rotation
        p_Ds(i) = dalpha / dbeta
        p_Dc(i) = p_Dh(i,i) / dbeta
        
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
        dpseudores = abs(p_Dq(i+1))
        
        if (ilevel .eq. nlmax) then
          ! Print our current pseudo-residual-norm
          if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
            call output_line ('GMRES('//&
              trim(sys_siL(idim,10))//'): Iteration '//&
              trim(sys_siL(ite,10))//',  !q(i+1)! = '//&
              trim(sys_sdEL(dpseudores,15)) )
          end if
        end if

        ! The euclid norm of our current defect is implicitly given by
        ! |q(i+1)|, however, it may be inaccurate, therefore, checking if
        ! |q(i+1)| fulfills the stopping criterion would not be very wise,
        ! as it may result in a lot of unnecessary restarts.
        ! We will therefore check (|q(i+1)| * dprnsf) against the tolerances
        ! instead, that should *hopefully* be okay.
        ! If dprnsf is 0, we will check |q(i+1)| against EPS to avoid NaNs
        ! or Infs during the next GMRES iteration.
        dtmp = dpseudores * dprnsf
        
        ! Is |q(i+1)| smaller than our machine`s exactness?
        ! If yes, then exit the inner GMRES loop.
        if (dpseudores .le. SYS_EPSREAL_DP) then
          exit
        
        ! Check if (|q(i+1)| * dprnsf) is greater than the machine`s
        ! exactness (this can only happen if dprnsf is positive).
        ! If yes, call linsol_testConvergence to test against the
        ! stopping criterions.
        else if (dtmp > SYS_EPSREAL_DP) then
          
          ! If (|q(i+1)| * dprnsf) fulfills the stopping criterion
          ! then exit the inner GMRES loop
          if (linsol_testConvergence(rsolverNode,dtmp)) exit
        
          ! We also want to check if our solution is diverging.
          if (linsol_testDivergence(rsolverNode,dpseudores)) then
            if(rsolverNode%ioutputLevel .ge. 2) then
              call output_line ('GMRES('// trim(sys_siL(idim,10))//&
                '): Warning: Pseudo-residuals diverging!')
            end if
            
            ! Instead of exiting the subroutine, we just exit the inner loop
            exit
          end if
        
        end if
        
        ! Maybe we have already done enough iterations?
        if (ite .ge. rsolverNode%nmaxIterations) exit
        
        ! Okay, next inner loop iteration, please
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= End of Inner Loop
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      
      ! Step O.3:
      ! Get the dimension of the current basis of the krylov subspace
      i = min(i, idim)
      
      ! Step O.4:
      ! Solve H' * q = q , where H' is the upper triangular matrix of
      ! the upper left (i x i)-block of H.
      ! We use the BLAS Level 2 subroutine 'dtrsv' to do the work for us
      call dtrsv('U', 'N', 'N', i, p_Dh, idim, p_Dq, 1)
      
      ! Step O.5:
      ! Update our solution vector
      ! x = x + q(1)*z(1) + q(2)*z(2) + ... + q(i)*z(i)
      do k = 1, i
        call lsysbl_vectorLinearComb(p_rz(k), rx, p_Dq(k), 1.0_DP)
      end do
      
      ! Step O.6:
      ! Calculate 'real' residual
      ! v(1) = b - (A * P * x)
      call lsysbl_copyVector (rb, p_rv(1))
      call lsysbl_copyVector (rx, p_rdefTempVector)
      
      ! Preconditioning
      if (bprec) then
        call linsol_precondDefect (&
            p_rlevelInfo%p_rpreconditioner,p_rdefTempVector)
      end if

      call lsysbl_blockMatVec(p_rmatrix, p_rdefTempVector, p_rv(1), -1.0_DP, 1.0_DP)

      ! Call filter chain if given.
      if (bfilter) then
        call filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
      end if
      
      ! Step O.7:
      ! Calculate euclid norm of the residual (needed for next q)
      dres = lsysbl_vectorNorm (p_rv(1), LINALG_NORMEUCLID)
      
      ! Calculate residual norm for stopping criterion
      ! TODO: try to avoid calculating the norm twice
      dfr = lsysbl_vectorNorm (p_rv(1), rsolverNode%iresNorm)
      
      ! Step O.8:
      ! Test for convergence, divergence and write some output now
      
      ! Shift the queue with the last residuals and add the new
      ! residual to it.
      Dresqueue = eoshift(Dresqueue,1,dfr)

      rsolverNode%dlastDefect = rsolverNode%dfinalDefect
      rsolverNode%dfinalDefect = dfr

      if (ilevel .eq. nlmax) then
        ! Test if the iteration is diverged
        if (linsol_testDivergence(rsolverNode,dfr)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('GMRES: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
   
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (linsol_testConvergence(rsolverNode,dfr)) exit
          
        end if
      
        ! We need to check if the euclidian norm of the
        ! residual is maybe < EPS - if this is true,
        ! then starting another GMRES iteration would
        ! result in NaNs / Infs!!!
        if (dres .le. SYS_EPSREAL_DP) exit

        ! print out the current residuum
  
        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Iteration '//&
            trim(sys_siL(ite,10))//',  !!RES!!  = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if

      else
      
        ! We need to check if the euclidian norm of the
        ! residual is maybe < EPS - if this is true,
        ! then starting another GMRES iteration would
        ! result in NaNs / Infs!!!
        if (dres .le. SYS_EPSREAL_DP) exit

      end if

    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! -= End of Outer Loop
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    ! On the highest level, calculate the actual solution by applying 
    ! the preconditioner
    if (bprec .and. (ilevel .eq. nlmax)) then
      call linsol_precondDefect (p_rlevelInfo%p_rpreconditioner,rx)
    end if

    if (ilevel .eq. nlmax) then
      ! Set ite to rsolverNode%nmaxIterations to prevent printing of
      ! "rsolverNode%nmaxIterations+1" of the loop was completed
      if (ite .gt. nmaxIterations) then
        ! Warning if we did not reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        if ((rsolverNode%ioutputLevel .ge. 0) .and. &
            (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
          call output_line ('GMRES: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        end if
  
        ite = rsolverNode%nmaxIterations
      end if

      ! Finish - either with an error or if converged.
      ! Print the last residuum.
      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ite .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('GMRES('// trim(sys_siL(idim,10))//&
          '): Iteration '//&
          trim(sys_siL(ite,10))//',  !!RES!!  = '//&
          trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if
  
      rsolverNode%iiterations = ite
    
      ! Do not calculate anything if the final residuum is out of bounds -
      ! would result in NaN`s,...
      if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
        ! Calculate asymptotic convergence rate
        if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
          i = min(rsolverNode%iiterations,niteAsymptoticCVR)
          rsolverNode%dasymptoticConvergenceRate = &
            (rsolverNode%dfinalDefect / &
              dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
        end if
  
        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting convergence rate stays zero.
        ! In the other case the convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/real(rsolverNode%iiterations,DP))
        else
          rsolverNode%dconvergenceRate = 0.0_DP
        end if
  
        if (rsolverNode%ioutputLevel .ge. 2) then
          call output_lbrk()
          call output_line ('GMRES('// trim(sys_siL(idim,10))// ') statistics:')
          call output_lbrk()
          call output_line ('Iterations              : '//&
               trim(sys_siL(rsolverNode%iiterations,10)) )
          call output_line ('!!INITIAL RES!!         : '//&
               trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
          call output_line ('!!RES!!                 : '//&
               trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
            call output_line ('!!RES!!/!!INITIAL RES!! : '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
          else
            call output_line ('!!RES!!/!!INITIAL RES!! : '//&
                 trim(sys_sdEL(0.0_DP,15)) )
          end if
          call output_lbrk ()
          call output_line ('Rate of convergence     : '//&
               trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
  
        end if
  
        if (rsolverNode%ioutputLevel .eq. 1) then
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Iterations/Rate of convergence: ' //&
                trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
                trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
        end if
  
      else
        ! DEF=Infinity; RHO=Infinity, set to 1
        rsolverNode%dconvergenceRate = 1.0_DP
        rsolverNode%dasymptoticConvergenceRate = 1.0_DP
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_precDeflGMRES_recleft (rsolverNode,ilevel,rx,rb)
  
!<description>
  ! Actual working routine for the deflated multilevel GMRES solver.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the deflaged GMRES solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
  
  ! Current level.
  integer, intent(in) :: ilevel
   
  ! Solution vector on level iilevel.
  ! Is replaced by a new approximation to the solution
  type(t_vectorBlock), intent(inout)        :: rx

  ! Right hand side vector on level ilevel
  type(t_vectorBlock), intent(inout), target :: rb
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_linsolDeflGMRESLevelInfo), pointer :: p_rlevelInfo
    type(t_linsolDeflGMRESLevelInfo), pointer :: p_rlevelInfoBelow
    type(t_linsolSubnodeDeflGMRES), pointer :: p_rsubnode
    type(t_linsolNode), pointer :: p_rprecSubnode
    integer :: nminIterations,nmaxIterations,nlmax
    
    real(DP) :: dalpha,dbeta,dres,dfr,dtmp,dpseudores,dprnsf
    integer :: ite,i,j,k,niteAsymptoticCVR
    
    ! Here come our 1D/2D arrays
    real(DP), dimension(:), pointer :: p_Dc, p_Ds, p_Dq
    real(DP), dimension(:,:), pointer :: p_Dh
  
    ! The queue saves the current residual and the two previous residuals.
    real(DP), dimension(LINSOL_NRESQUEUELENGTH) :: Dresqueue
    
    ! The system matrix
    type(t_matrixBlock), pointer :: p_rmatrix
    
    ! Minimum number of iterations, print-sequence for residuals, dimension
    ! of krylov subspace
    integer :: niteResOutput, idim, hDq
    
    ! Whether to filter/precondition, apply Gram-Schmidt twice
    logical bprec,bfilter,btwiceGS
    
    ! Pointers to temporary vectors - named for easier access
    type(t_vectorBlock), pointer :: p_rtempVector,p_rdefTempVector,p_ractualRHS
    type(t_vectorBlock), dimension(:), pointer :: p_rv, p_rz
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! Get pointers to the current level and the GMRES parameters
    p_rsubnode => rsolverNode%p_rsubnodeDeflGMRES
    p_rlevelInfo => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel)
    
    nullify(p_rlevelInfoBelow)
    if (ilevel .gt. 1) then
      ! The lower level
      p_rlevelInfoBelow => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(ilevel-1)
    end if
    
    ! On the maximum level, apply as many iterations as configured in the
    ! solver node. On the other levels, apply as many as configured in the
    ! level structure
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels

    if (ilevel .eq. nlmax) then
      nminIterations = max(rsolverNode%nminIterations,0)
      nmaxIterations = max(rsolverNode%nmaxIterations,0)

      ! Length of the queue of last residuals for the computation of
      ! the asymptotic convergence rate
      niteAsymptoticCVR = max(0,min(LINSOL_NRESQUEUELENGTH-1,rsolverNode%niteAsymptoticCVR))
  
      ! Iteration when the residuum is printed:
      niteResOutput = max(1,rsolverNode%niteResOutput)
      
    else
      if (p_rlevelInfo%niterations .eq. -1) then
        nminIterations = max(p_rlevelInfo%ikrylowDim,0)
        nmaxIterations = max(p_rlevelInfo%ikrylowDim,0)
      else
        nminIterations = max(p_rlevelInfo%niterations,0)
        nmaxIterations = max(p_rlevelInfo%niterations,0)
      end if
      
      niteAsymptoticCVR = 0
      niteResOutput = 0
      
    end if

    ! Use preconditioning? Filtering?
    bprec = associated(p_rlevelInfo%p_rpreconditioner)
    bfilter = associated(p_rsubnode%p_RfilterChain)
    
    ! Apply Gram-Schmidt twice?
    btwiceGS = p_rsubnode%btwiceGS
    
    ! Get the dimension of Krylov subspace
    idim = p_rlevelInfo%ikrylowDim
    
    ! Get our pseudo-residual-norm scale factor
    dprnsf = p_rsubnode%dpseudoResScale
    
    ! Now we need to check the pseudo-residual scale factor.
    ! If it is non-positive, we set it to 0.
    if (dprnsf .le. 0.0_DP) then
      dprnsf = 0.0_DP
    end if

    ! Current matrix and other parameters
    p_rmatrix => p_rlevelInfo%rsystemMatrix
    
    ! Set pointers to the temporary vectors
    ! defect vectors
    p_rz => p_rlevelInfo%p_rz
    ! basis vectors
    p_rv => p_rlevelInfo%p_rv
    ! temp vector
    p_rtempVector => p_rlevelInfo%rtempVector
    
    ! temp vector for defect creation with preconditioning
    p_rdefTempVector => p_rlevelInfo%rdefTempVector

    ! All vectors share the same boundary conditions as rd!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    do i=1, idim+1
      call lsysbl_assignDiscrIndirect (rb,p_rv(i))
    end do
    do i=1, idim
      call lsysbl_assignDiscrIndirect (rb,p_rz(i))
    end do
    
    ! Set pointers to the 1D/2D arrays
    p_Dh => p_rlevelInfo%Dh
    p_Dc => p_rlevelInfo%Dc
    p_Ds => p_rlevelInfo%Ds
    p_Dq => p_rlevelInfo%Dq
    
    ! And get the handle of Dq, since we need to call
    ! storage_clear during the iteration
    hDq = p_rlevelInfo%hDq
    
    ! All vectors share the same boundary conditions as rb!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    do i=1,idim
      call lsysbl_assignDiscrIndirect (rb, p_rz(i))
    end do
    do i=1,idim+1
      call lsysbl_assignDiscrIndirect (rb, p_rv(i))
    end do
    
    if (bprec) then
      p_rprecSubnode => p_rlevelInfo%p_rpreconditioner
    end if
    
    if (bfilter) then
      p_RfilterChain => p_rsubnode%p_RfilterChain
    end if
    
    ! ---------------------------------------------------------------
    ! Create the actual RHS vector
    ! ---------------------------------------------------------------
    ! We apply left preconditioning.
    ! On the maximum level, apply the precodnitioner to the RHS
    ! to create the actual RHS.
    ! On the lower levels, the RHS contains already the preconditioner,
    ! so no preconditioning has to be applied.
    
    p_ractualRHS => rb
    if (bprec) then
    
      if (ilevel .eq. nlmax) then
        ! Maximum level. Apply the preconditioner to the RHS.
        p_ractualRHS => p_rlevelInfo%ractualRHS
        
        call lsysbl_copyVector (rb,p_ractualRHS)

        call linsol_precondDefect (&
            p_rlevelInfo%p_rpreconditioner,p_ractualRHS)

      end if
      
    end if
    
    ! Copy our RHS to p_rv(1). As the iteration vector is 0, this
    ! is also our initial defect.
    call lsysbl_copyVector(p_ractualRHS,p_rv(1))

    if (bfilter) then
      call filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
    end if
    
    ! Get the norm of the residuum.
    ! We need to calculate the norm twice, since we need one for
    ! the stopping criterion (selected by the user) and the
    ! euclidian norm for the internal GMRES algorithm
    dfr = lsysbl_vectorNorm (p_ractualRHS,rsolverNode%iresNorm)
    dres = lsysbl_vectorNorm (p_rv(1),LINALG_NORMEUCLID)
    if (.not.((dfr .ge. 1E-99_DP) .and. &
              (dfr .le. 1E99_DP))) dfr = 0.0_DP

    ! initialise the queue of the last residuals with RES
    Dresqueue = dfr

    if (ilevel .eq. nlmax) then
      ! Initialise starting residuum
      rsolverNode%dinitialDefect = dfr
      rsolverNode%dlastDefect = 0.0_DP

      ! Check if out initial defect is zero. This may happen if the filtering
      ! routine filters "everything out"!
      ! In that case we can directly stop our computation.
  
      if (rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero) then
        ! final defect is 0, as initialised in the output variable above
        call lsysbl_clearVector(rx)
        ite = 0
        rsolverNode%dfinalDefect = dfr
        return
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('GMRES('//&
             trim(sys_siL(idim,10))//'): Iteration '//&
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

    end if

    ite = 0
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! -= Outer Loop
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    do while (ite < nmaxIterations)
      
      ! Step O.1:
      ! Create elementary vector e_1 scaled by the euclid norm of the
      ! defect, i.e.: q = ( ||v(1)||_2, 0, ..., 0 )
      call storage_clear (hDq)
      p_Dq(1) = dres
      
      ! Step O.2:
      ! Now scale the defect by the inverse of its norm
      call lsysbl_scaleVector (p_rv(1), 1.0_DP / dres)
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= Inner Loop (GMRES iterations)
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      do i = 1, idim
        ! Another iteration
        ite = ite + 1
      
        ! Step I.1:
        ! Solve P * z(i) = v(i), where P is the preconditioner matrix
        call lsysbl_copyVector (p_rv(i), p_rz(i))
        
        ! Apply deflation preconditioning to z(i)
        if (associated(p_rlevelInfoBelow)) then
        
          ! Calculate z = (P A z  -  drelax lambda_max v(i))
          ! in the temp vector
          call lsysbl_copyVector (p_rv(i),p_rtempVector)
          
          call lsysbl_blockMatVec (p_rlevelInfo%rsystemMatrix, p_rz(i),p_rdefTempVector, &
              1.0_DP, 0.0_DP)
          
          ! Preconditioning
          if (bprec) then
            call linsol_precondDefect (&
                p_rlevelInfo%p_rpreconditioner,p_rdefTempVector)
          end if
          
          call lsysbl_vectorLinearComb (p_rdefTempVector, p_rtempVector, &
              1.0_DP, -p_rsubnode%drelax*p_rlevelInfo%dmaxEigenvalue)

          ! Apply the filter chain
          if (bfilter) then
            call filter_applyFilterChainVec (p_rtempVector, p_rsubnode%p_RfilterChain)
          end if
          
          ! Restriction to the RHS of the lower level
          call mlprj_performRestriction (p_rlevelInfo%p_rprojection,&
                p_rlevelInfoBelow%rrhsVector, &
                p_rtempVector, &
                p_rsubnode%rprjTempVector)

          ! Apply the filter chain to the new RHS
          if (bfilter) then
            call filter_applyFilterChainVec (p_rlevelInfoBelow%rrhsVector, p_rsubnode%p_RfilterChain)
          end if
          
          ! Preconditioning on the lower level
          call lsysbl_clearVector (p_rlevelInfoBelow%rsolutionVector)
          call linsol_precDeflGMRES_recleft (rsolverNode,ilevel-1,&
              p_rlevelInfoBelow%rsolutionVector,p_rlevelInfoBelow%rrhsVector)

          ! Prolongation to our current level: t = I(...)
          call mlprj_performProlongation (p_rlevelInfo%p_rprojection,&
                p_rlevelInfoBelow%rsolutionVector, &
                p_rtempVector, &
                p_rsubnode%rprjTempVector)

          ! Correction: z(i) = v(i) - t
          call lsysbl_vectorLinearComb (p_rtempVector,p_rz(i),-1.0_DP,1.0_DP)
        
        end if
        
        ! Apply filter chain to z(i)
        if (bfilter) then
          call filter_applyFilterChainVec (p_rz(i), p_RfilterChain)
        end if
        
        
        ! Step I.2:
        ! Calculate v(i+1) = P * A * z(i)
        call lsysbl_blockMatVec (p_rmatrix, p_rz(i), p_rv(i+1), 1.0_DP, 0.0_DP)
        
        ! Preconditioning + filtering
        if (bprec) then
          call linsol_precondDefect (&
              p_rlevelInfo%p_rpreconditioner,p_rv(i+1))
        end if

        if (bfilter) then
          call filter_applyFilterChainVec (p_rv(i+1), p_RfilterChain)
        end if
        
        
        ! Step I.3:
        ! Perfom Gram-Schmidt process (1st time)
        do k = 1, i
          p_Dh(k,i) = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
          call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -p_Dh(k,i), 1.0_DP)
        end do
        
        ! If the user wishes, perform Gram-Schmidt one more time
        ! to improve numerical stability.
        if (btwiceGS) then
          do k = 1, i
            dtmp = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
            p_Dh(k,i) = p_Dh(k,i) + dtmp
            call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -dtmp, 1.0_DP)
          end do
        end if
        
        ! Step I.4:
        ! Calculate alpha = ||v(i+1)||_2
        dalpha = lsysbl_vectorNorm (p_rv(i+1), LINALG_NORMEUCLID)
        
        ! Step I.5:
        ! Scale v(i+1) by the inverse of its euclid norm
        if (dalpha > SYS_EPSREAL_DP) then
          call lsysbl_scaleVector (p_rv(i+1), 1.0_DP / dalpha)
        else
          ! Well, let us just print a warning here...
          if(rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('GMRES('// trim(sys_siL(idim,10))//&
              '): Warning: !!v(i+1)!! < EPS !')
          end if
        end if
        
        
        ! Step I.6:
        ! Apply Givens rotations to the i-th column of h which
        ! renders the Householder matrix an upper triangular matrix
        do k = 1, i-1
          dtmp = p_Dh(k,i)
          p_Dh(k,i)   = p_Dc(k) * dtmp + p_Ds(k) * p_Dh(k+1,i)
          p_Dh(k+1,i) = p_Ds(k) * dtmp - p_Dc(k) * p_Dh(k+1,i)
        end do
        
        
        ! Step I.7:
        ! Calculate Beta
        ! Beta = (h(i,i)^2 + alpha^2) ^ (1/2)
        dbeta = sqrt(p_Dh(i,i)**2 + dalpha**2)
        if (dbeta < SYS_EPSREAL_DP) then
          dbeta = SYS_EPSREAL_DP
          if(rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('GMRES('// trim(sys_siL(idim,10))//&
              '): Warning: beta < EPS !')
          end if
        end if
        
        ! Step I.8:
        ! Calculate next plane rotation
        p_Ds(i) = dalpha / dbeta
        p_Dc(i) = p_Dh(i,i) / dbeta
        
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
        dpseudores = abs(p_Dq(i+1))
        
        if (ilevel .eq. nlmax) then
          ! Print our current pseudo-residual-norm
          if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
            call output_line ('GMRES('//&
              trim(sys_siL(idim,10))//'): Iteration '//&
              trim(sys_siL(ite,10))//',  !q(i+1)! = '//&
              trim(sys_sdEL(dpseudores,15)) )
          end if
        end if

        ! The euclid norm of our current defect is implicitly given by
        ! |q(i+1)|, however, it may be inaccurate, therefore, checking if
        ! |q(i+1)| fulfills the stopping criterion would not be very wise,
        ! as it may result in a lot of unnecessary restarts.
        ! We will therefore check (|q(i+1)| * dprnsf) against the tolerances
        ! instead, that should *hopefully* be okay.
        ! If dprnsf is 0, we will check |q(i+1)| against EPS to avoid NaNs
        ! or Infs during the next GMRES iteration.
        dtmp = dpseudores * dprnsf
        
        ! Is |q(i+1)| smaller than our machine`s exactness?
        ! If yes, then exit the inner GMRES loop.
        if (dpseudores .le. SYS_EPSREAL_DP) then
          exit
        
        ! Check if (|q(i+1)| * dprnsf) is greater than the machine`s
        ! exactness (this can only happen if dprnsf is positive).
        ! If yes, call linsol_testConvergence to test against the
        ! stopping criterions.
        else if (dtmp > SYS_EPSREAL_DP) then
          
          ! If (|q(i+1)| * dprnsf) fulfills the stopping criterion
          ! then exit the inner GMRES loop
          if (linsol_testConvergence(rsolverNode,dtmp)) exit
        
          ! We also want to check if our solution is diverging.
          if (linsol_testDivergence(rsolverNode,dpseudores)) then
            if(rsolverNode%ioutputLevel .ge. 2) then
              call output_line ('GMRES('// trim(sys_siL(idim,10))//&
                '): Warning: Pseudo-residuals diverging!')
            end if
            
            ! Instead of exiting the subroutine, we just exit the inner loop
            exit
          end if
        
        end if
        
        ! Maybe we have already done enough iterations?
        if (ite .ge. rsolverNode%nmaxIterations) exit
        
        ! Okay, next inner loop iteration, please
      end do
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! -= End of Inner Loop
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      
      ! Step O.3:
      ! Get the dimension of the current basis of the krylov subspace
      i = min(i, idim)
      
      ! Step O.4:
      ! Solve H' * q = q , where H' is the upper triangular matrix of
      ! the upper left (i x i)-block of H.
      ! We use the BLAS Level 2 subroutine 'dtrsv' to do the work for us
      call dtrsv('U', 'N', 'N', i, p_Dh, idim, p_Dq, 1)
      
      ! Step O.5:
      ! Update our solution vector
      ! x = x + q(1)*z(1) + q(2)*z(2) + ... + q(i)*z(i)
      do k = 1, i
        call lsysbl_vectorLinearComb(p_rz(k), rx, p_Dq(k), 1.0_DP)
      end do
      
      ! Step O.6:
      ! Calculate 'real' residual
      ! v(1) = b - (A * x)
      call lsysbl_copyVector (p_ractualRHS, p_rv(1))
      
      call lsysbl_blockMatVec(p_rmatrix, rx, p_rdefTempVector, 1.0_DP, 0.0_DP)
      
      ! Preconditioning
      if (bprec) then
        call linsol_precondDefect (&
            p_rlevelInfo%p_rpreconditioner,p_rdefTempVector)
      end if

      call lsysbl_vectorLinearComb(p_rdefTempVector, p_rv(1), -1.0_DP, 1.0_DP)

      ! Call filter chain if given.
      if (bfilter) then
        call filter_applyFilterChainVec (p_rv(1), p_RfilterChain)
      end if
      
      ! Step O.7:
      ! Calculate euclid norm of the residual (needed for next q)
      dres = lsysbl_vectorNorm (p_rv(1), LINALG_NORMEUCLID)
      
      ! Calculate residual norm for stopping criterion
      ! TODO: try to avoid calculating the norm twice
      dfr = lsysbl_vectorNorm (p_rv(1), rsolverNode%iresNorm)
      
      ! Step O.8:
      ! Test for convergence, divergence and write some output now
      
      ! Shift the queue with the last residuals and add the new
      ! residual to it.
      Dresqueue = eoshift(Dresqueue,1,dfr)

      rsolverNode%dlastDefect = rsolverNode%dfinalDefect
      rsolverNode%dfinalDefect = dfr

      if (ilevel .eq. nlmax) then
        ! Test if the iteration is diverged
        if (linsol_testDivergence(rsolverNode,dfr)) then
          if (rsolverNode%ioutputLevel .lt. 2) then
            do i=max(1,size(Dresqueue)-ite-1+1),size(Dresqueue)
              j = ite-max(1,size(Dresqueue)-ite+1)+i
              call output_line ('GMRES: Iteration '// &
                  trim(sys_siL(j,10))//',  !!RES!! = '//&
                  trim(sys_sdEL(Dresqueue(i),15)) )
            end do
          end if
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
   
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (linsol_testConvergence(rsolverNode,dfr)) exit
          
        end if
      
        ! We need to check if the euclidian norm of the
        ! residual is maybe < EPS - if this is true,
        ! then starting another GMRES iteration would
        ! result in NaNs / Infs!!!
        if (dres .le. SYS_EPSREAL_DP) exit

        ! print out the current residuum
  
        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Iteration '//&
            trim(sys_siL(ite,10))//',  !!RES!!  = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if

      else
      
        ! We need to check if the euclidian norm of the
        ! residual is maybe < EPS - if this is true,
        ! then starting another GMRES iteration would
        ! result in NaNs / Infs!!!
        if (dres .le. SYS_EPSREAL_DP) exit

      end if

    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! -= End of Outer Loop
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if (ilevel .eq. nlmax) then
      ! Set ite to rsolverNode%nmaxIterations to prevent printing of
      ! "rsolverNode%nmaxIterations+1" of the loop was completed
      if (ite .gt. nmaxIterations) then
        ! Warning if we did not reach the convergence criterion.
        ! No warning if we had to do exactly rsolverNode%nmaxIterations steps
        if ((rsolverNode%ioutputLevel .ge. 0) .and. &
            (rsolverNode%nmaxIterations .gt. rsolverNode%nminIterations)) then
          call output_line ('GMRES: Accuracy warning: '//&
              'Solver did not reach the convergence criterion')
        end if
  
        ite = rsolverNode%nmaxIterations
      end if

      ! Finish - either with an error or if converged.
      ! Print the last residuum.
      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ite .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('GMRES('// trim(sys_siL(idim,10))//&
          '): Iteration '//&
          trim(sys_siL(ite,10))//',  !!RES!!  = '//&
          trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if
  
      rsolverNode%iiterations = ite
    
      ! Do not calculate anything if the final residuum is out of bounds -
      ! would result in NaN`s,...
      if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
        ! Calculate asymptotic convergence rate
        if ((niteAsymptoticCVR .gt. 0) .and. (rsolverNode%iiterations .gt. 0)) then
          i = min(rsolverNode%iiterations,niteAsymptoticCVR)
          rsolverNode%dasymptoticConvergenceRate = &
            (rsolverNode%dfinalDefect / &
              dresqueue(LINSOL_NRESQUEUELENGTH-i))**(1.0_DP/real(I,DP))
        end if
  
        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting convergence rate stays zero.
        ! In the other case the convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/real(rsolverNode%iiterations,DP))
        else
          rsolverNode%dconvergenceRate = 0.0_DP
        end if
  
        if (rsolverNode%ioutputLevel .ge. 2) then
          call output_lbrk()
          call output_line ('GMRES('// trim(sys_siL(idim,10))// ') statistics:')
          call output_lbrk()
          call output_line ('Iterations              : '//&
               trim(sys_siL(rsolverNode%iiterations,10)) )
          call output_line ('!!INITIAL RES!!         : '//&
               trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
          call output_line ('!!RES!!                 : '//&
               trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
            call output_line ('!!RES!!/!!INITIAL RES!! : '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
          else
            call output_line ('!!RES!!/!!INITIAL RES!! : '//&
                 trim(sys_sdEL(0.0_DP,15)) )
          end if
          call output_lbrk ()
          call output_line ('Rate of convergence     : '//&
               trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
  
        end if
  
        if (rsolverNode%ioutputLevel .eq. 1) then
          call output_line ('GMRES('// trim(sys_siL(idim,10))//&
            '): Iterations/Rate of convergence: ' //&
                trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
                trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
        end if
  
      else
        ! DEF=Infinity; RHO=Infinity, set to 1
        rsolverNode%dconvergenceRate = 1.0_DP
        rsolverNode%dasymptoticConvergenceRate = 1.0_DP
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_precDeflGMRES (rsolverNode,rd)
  
!<description>
  ! Applies Multigrid preconditioner <tex>$ P \approx A $</tex> to the defect
  ! vector rd and solves <tex>$ Pd_{new} = d $</tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The structure of the levels (pre/postsmoothers, interlevel projection
  ! structures) must have been prepared by setting the level parameters
  ! before calling this routine. (The parameters of a specific level
  ! can be get by calling linsol_getMultigridLevel2).
  ! The matrices <tex>$ A $</tex> on all levels must be attached to the solver previously
  ! by linsol_setMatrices.
  
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Multigrid solver
  type(t_linsolNode), intent(inout), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: nlmax
    
    ! Our MG structure
    type(t_linsolSubnodeDeflGMRES), pointer :: p_rsubnode
    type(t_linsolDeflGMRESLevelInfo), pointer :: p_rlevelInfo
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeDeflGMRES
    
    if (.not. associated(rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo)) then
      call output_line ("No levels attached!", &
          OU_CLASS_ERROR, OU_MODE_STD, "linsol_precDeflGMRES")
      rsolverNode%iresult = 2
      return
    end if

    ! Call GMRES on the maximum level
    nlmax = rsolverNode%p_rsubnodeDeflGMRES%nlevels
    p_rlevelInfo => rsolverNode%p_rsubnodeDeflGMRES%p_RlevelInfo(nlmax)
    
    call lsysbl_clearVector (p_rlevelInfo%rsolutionVector)
    
    if (rsolverNode%p_rsubnodeDeflGMRES%brightprec) then
      ! Right preconditioning
      call linsol_precDeflGMRES_recright (&
          rsolverNode,nlmax,p_rlevelInfo%rsolutionVector,rd)
    else
      ! Left preconditioning
      call linsol_precDeflGMRES_recleft (&
          rsolverNode,nlmax,p_rlevelInfo%rsolutionVector,rd)
    end if
        
    ! Overwrite our previous RHS by the new correction vector.
    ! This completes the preconditioning. Scale the vector by the damping parameter
    ! in the solver structure.
    call lsysbl_vectorLinearComb (p_rlevelInfo%rsolutionVector,rd,rsolverNode%domega,0.0_DP)
  
  end subroutine
  
  
  
! *****************************************************************************
! Routines for the Block-SOR solver (for DG discretisations)
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initBlockSOR (p_rsolverNode,drelax)
  
!<description>
  ! Creates a t_linsolNode solver structure for the Block-SOR solver. The
  ! node can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  ! It is a special solver taking advantage of the block-structure stemming
  ! from a discontinuous Galerkin discretisation.

!</description>
  
!<input>
  ! Relaxation parameter (=1 for block Gauﬂ Seidel)
  real(dp), intent(in), optional :: drelax
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_BLOCKSOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_DIRECT &
                                + LINSOL_ABIL_BLOCK 

    ! Allocate the subnode for SOR.
    allocate(p_rsolverNode%p_rsubnodeSOR)
    
    ! Save the relaxation parameter.
    p_rsolverNode%p_rsubnodeSOR%drelax = 1.0_DP
    if (present(drelax)) then
      p_rsolverNode%p_rsubnodeSOR%drelax = drelax
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatBlockSOR (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  integer :: i,n
  type(t_matrixScalar), pointer :: p_rblock => null()

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Get the index of the top-level matrix
    n = ubound(Rmatrices,1)
    
    ! Let us make sure the matrix is not rectangular.
    if(Rmatrices(n)%nblocksPerRow .ne. Rmatrices(n)%nblocksPerCol) then
      
      ! Rectangular block matrix
      ccompatible = LINSOL_COMP_ERRMATRECT
      
    else
    
      ! Okay, let us loop through all main diagonal blocks.
      do i = 1, Rmatrices(n)%nblocksPerRow
      
        ! Get the diagonal block
        p_rblock => Rmatrices(n)%RmatrixBlock(i,i)
        
        ! Is the matrix empty?
        if(p_rblock%NEQ .eq. 0) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! Is the scaling factor zero?
        if(p_rblock%dscaleFactor .eq. 0.0_DP) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! What about the matrix type?
        if((p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
           (p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
          ! Matrix type not supported
          ccompatible = LINSOL_COMP_ERRMATTYPE
          exit
        end if
        
        ! This block is okay, continue with next one
      
      end do
    
    end if
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precBlockSOR (rsolverNode,rd)
  
  use triangulation
  use spatialdiscretisation
  use element
  use dofmapping
  
!<description>
  ! Applies the Jacobi preconditioner <tex>$ D \in A $<tex> to the defect
  ! vector rd and solves <tex>$ Dd_{new} = d $<tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Jacobi solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: iblock, jblock
    
    type(t_matrixScalar), pointer :: p_rmatrix
    integer, dimension(:), pointer :: p_KCOL, p_KLD
    integer :: NEL, iel, indof, i, ig, j, jg,iloc
    integer(I32) :: celement
    real(DP) :: dlocOmega
    real(DP), dimension(:), pointer :: p_Dvector, p_Dmatrix
    real(dp), dimension(:,:), allocatable :: DlocMat
    real(dp), dimension(:), allocatable :: DlocVec, DlocSol
    integer, dimension(:), allocatable :: IdofGlob
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    type(t_elementDistribution), pointer :: p_relemDist
    logical :: bsuccess
    real(dp) :: drelax, domega
    integer, dimension(:), pointer :: p_IelementDistr
    
      ! Status reset
      rsolverNode%iresult = 0
      
      ! Get omega and the relaxation parameter
      domega = rsolverNode%domega
      drelax = rsolverNode%p_rsubnodeSOR%drelax

      ! Loop through all blocks. Each block corresponds to one
      ! diagonal block in the matrix.
      do iblock = 1, rd%nblocks
        
        ! Loop through all blocks in the lower triangular block matrix.
        do jblock = 1, iblock-1
      
          ! Get the matrix
          p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,jblock)
        
          ! Is it empty?
          if(p_rmatrix%NEQ .eq. 0) cycle
        
          ! Otherwise update the defect.
          call lsyssc_scalarMatVec(p_rmatrix, rd%RvectorBlock(jblock), &
                                   rd%RvectorBlock(iblock), drelax, 1.0_DP)
        
        end do
      
        ! Get the main diagonal matrix
        p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(iblock,iblock)
        
        ! Now we have to make some decisions. At first, which matrix
        ! structure do we have?
        select case (p_rmatrix%cmatrixFormat)
        case (LSYSSC_MATRIX9)
        
          ! Get the omega parameter for the matrix. Do not forget the scaling
          ! factor in the matrix!
          dlocomega = rsolverNode%domega / p_rmatrix%dScaleFactor

          ! Now, which data format do we have? Single or double?
          select case (p_rmatrix%cdataType)
          case (ST_DOUBLE)
            ! Get pointers to the matrix data
            call lsyssc_getbase_Kcol (p_rmatrix,p_KCOL)
            call lsyssc_getbase_Kld (p_rmatrix,p_KLD)
            call lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
            
            ! Take care of the accuracy of the vector
            select case (rd%RvectorBlock(iblock)%cdataType)
            case (ST_DOUBLE)
              ! Get the data array
              call lsyssc_getbase_double (rd%RvectorBlock(iblock),p_Dvector)
              
              ! Get pointers for quicker access
              p_rspatialDiscr => p_rmatrix%p_rspatialDiscrTest
              p_rtriangulation => p_rspatialDiscr%p_rtriangulation
              
              ! Get number of elements
              NEL = p_rtriangulation%NEL
              
              ! Check if we have a uniform discretisation (simplifies a lot)
              if (p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
              
                ! Get the element distribution from the spatial discretisation
                p_relemDist => p_rspatialDiscr%RelementDistr(1)

                ! Get type of element in this distribution
                celement =  p_relemDist%celement

                ! Get number of local DOF
                indof = elem_igetNDofLoc(celement)
                
                ! Allocate space for global DOFs
                allocate(IdofGlob(indof))

                ! Allocate local matrix and vector
                allocate(DlocMat(indof,indof), DlocVec(indof), DlocSol(indof))
              
                ! Loop over all elements in this element distribution
                do iel = 1, p_relemDist%NEL
                
                  ! Map local to global DOFs
                  call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)

                  !!! Get entries from the global into the local matrix
                  ! Loop over all local lines
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)

                    ! Loop over all local columns
                    do j = 1, indof
                      ! Get global columns
                      jg = IdofGlob(j)

                      do iloc = p_KLD(ig), p_KLD(ig+1)-1
                        if (jg.eq.p_KCOL(iloc)) exit
                      end do

                      DlocMat(i,j) = p_Dmatrix(iloc)

                    end do
                  end do

                  !!! Now get local vector
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)
                     
                    ! Get global entry in local vector
                    DlocVec(i) = dlocOmega * p_Dvector(ig)
                    
                    ! Substract product with lower diagonal
                    do jg = p_KLD(ig), p_KLD(ig+1)-1
                      if (p_KCOL(jg).ge.IdofGlob(1)) exit
                      DlocVec(i) = DlocVec(i) - p_Dmatrix(jg)*p_Dvector(p_KCOL(jg))
                    end do
                    
                  end do

                  !!! Now solve local linear system
                  call mprim_invertMatrixDble(DlocMat, DlocVec, DlocSol, &
                                              indof, 2, bsuccess)

                  ! Test if local system could be solved
                  if (.not.bsuccess) then
                    call output_line ( &
                      'Local matricx singular.', &
                      OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockSOR')
                      do i = 1,size(DlocMat,1)
                        write(*,*) DlocMat(i,:)
                      end do
                    call sys_halt()
                  end if

                  !!! Distribute local solution to global solution
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)
                     
                    ! Get local entry in global vector
                    p_Dvector(ig) = drelax * DlocSol(i)
                  end do
                  
                end do !iel
                
                ! Deallocate space for global DOFs
                deallocate(IdofGlob)

                ! Deallocate local matrix and vector
                deallocate(DlocMat, DlocVec, DlocSol)
              
              else ! No uniform distribution
              
                ! Gett array which shows, which element distribution one element is in
                call storage_getbase_int (p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
              
                ! Loop over all elements in this element distribution
                do iel = 1, p_relemDist%NEL
                
                  ! Get the element distribution from the spatial discretisation
                  p_relemDist => p_rspatialDiscr%RelementDistr(p_IelementDistr(iel))

                  ! Get type of element in this distribution
                  celement =  p_relemDist%celement

                  ! Get number of local DOF
                  indof = elem_igetNDofLoc(celement)
                  
                  ! Allocate space for global DOFs
                  allocate(IdofGlob(indof))

                  ! Allocate local matrix and vector
                  allocate(DlocMat(indof,indof), DlocVec(indof), DlocSol(indof))
                
                  ! Map local to global DOFs
                  call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)

                  !!! Get entries from the global into the local matrix
                  ! Loop over all local lines
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)

                    ! Loop over all local columns
                    do j = 1, indof
                      ! Get global columns
                      jg = IdofGlob(j)

                      do iloc = p_KLD(ig), p_KLD(ig+1)-1
                        if (jg.eq.p_KCOL(iloc)) exit
                      end do

                      DlocMat(i,j) = p_Dmatrix(iloc)

                    end do
                  end do

                  !!! Now get local vector
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)
                     
                    ! Get global entry in local vector
                    DlocVec(i) = dlocOmega * p_Dvector(ig)
                    
                    ! Substract product with lower diagonal
                    do jg = p_KLD(ig), p_KLD(ig+1)-1
                      if (p_KCOL(jg).ge.IdofGlob(1)) exit
                      DlocVec(i) = DlocVec(i) - p_Dmatrix(jg)*p_Dvector(p_KCOL(jg))
                    end do
                    
                  end do

                  !!! Now solve local linear system
                  call mprim_invertMatrixDble(DlocMat, DlocVec, DlocSol, &
                                              indof, 2, bsuccess)

                  ! Test if local system could be solved
                  if (.not.bsuccess) then
                    call output_line ( &
                      'Local matricx singular.', &
                      OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockSOR')
                      do i = 1,size(DlocMat,1)
                        write(*,*) DlocMat(i,:)
                      end do
                    call sys_halt()
                  end if

                  !!! Distribute local solution to global solution
                  do i = 1, indof
                    ! Get global line
                    ig = IdofGlob(i)
                     
                    ! Get local entry in global vector
                    p_Dvector(ig) = drelax * DlocSol(i)
                  end do
                  
                  ! Deallocate space for global DOFs
                  deallocate(IdofGlob)

                  ! Deallocate local matrix and vector
                  deallocate(DlocMat, DlocVec, DlocSol)
                  
                end do !iel
              
              end if ! Type of discretisation
              
            case default
              call output_line ('Unsupported vector precision.', &
                                OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockSOR')
              call sys_halt()
            end select
         
         
          case default
            call output_line ('Unsupported matrix precision.', &
                              OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockSOR')
            call sys_halt()
          end select
        
        case default
          call output_line ('Unsupported matrix format.', &
                            OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockSOR')
          call sys_halt()
        end select
        
      end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneBlockSOR (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the Block-SOR solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of SOR which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    deallocate(rsolverNode%p_rsubnodeSOR)
  
  end subroutine
  
  
! *****************************************************************************
! Routines for the Block Upwind GS solver (for DG discretisations)
! *****************************************************************************

!<subroutine>
  
  subroutine linsol_initBlockUpwGS (p_rsolverNode,Dvelocity,p_rtriangulation,&
                                    p_DNormals)
  
!<description>
  ! Creates a t_linsolNode solver structure for the Block-UpwGS solver. The
  ! node can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  ! It is a special solver taking advantage of the block-structure stemming
  ! from a discontinuous Galerkin discretisation.

!</description>
  
!<input>
  ! Velocity vector
  real(dp), dimension(2), intent(in) :: Dvelocity
  ! Pointer to the normals of the triangulation for each edge
  real(dp), dimension(:,:), pointer :: p_DNormals
  ! Corresponding triangulation
  type(t_triangulation), pointer :: p_rtriangulation
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_linsolNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call linsol_initSolverGeneral(p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = LINSOL_ALG_BLOCKUPWGS
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_DIRECT ! &
                                ! + LINSOL_ABIL_BLOCK 

    ! Allocate the subnode for BlockUpwGS
    allocate(p_rsolverNode%p_rsubnodeBlockUpwGS)
    
    ! Set velocity vector
    p_rsolverNode%p_rsubnodeBlockUpwGS%Dvel = Dvelocity
    
    ! Set normals
    p_rsolverNode%p_rsubnodeBlockUpwGS%p_Dnormals => p_Dnormals
    
    ! Set number of elements
    p_rsolverNode%p_rsubnodeBlockUpwGS%p_rtriangulation => p_rtriangulation
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_matCompatBlockUpwGS (rsolverNode,Rmatrices,&
      ccompatible,CcompatibleDetail)
  
!<description>
  ! This routine is called to check if the matrices in Rmatrices are
  ! compatible to the solver. Calls linsol_matricesCompatible for possible
  ! subsolvers to check the compatibility.
!</description>
  
!<input>
  ! The solver node which should be checked against the matrices
  type(t_linsolNode), intent(in)             :: rsolverNode

  ! An array of system matrices which is simply passed to the initialisation
  ! routine of the preconditioner.
  type(t_matrixBlock), dimension(:), intent(in)   :: Rmatrices
!</input>
  
!<output>
  ! A LINSOL_COMP_xxxx flag that tells the caller whether the matrices are
  ! compatible (which is the case if LINSOL_COMP_OK is returned).
  integer, intent(out) :: ccompatible

  ! OPTIONAL: An array of LINSOL_COMP_xxxx that tell for every level if
  ! the matrices on that level are ok or not. Must have the same size
  ! as Rmatrices!
  integer, dimension(:), intent(inout), optional :: CcompatibleDetail
!</output>
  
!</subroutine>

  integer :: i,n
  type(t_matrixScalar), pointer :: p_rblock => null()

    ! Normally we are compatible.
    ccompatible = LINSOL_COMP_OK
    
    ! Get the index of the top-level matrix
    n = ubound(Rmatrices,1)
    
    ! Let us make sure the matrix is not rectangular.
    if(Rmatrices(n)%nblocksPerRow .ne. Rmatrices(n)%nblocksPerCol) then
      
      ! Rectangular block matrix
      ccompatible = LINSOL_COMP_ERRMATRECT
      
    else
    
      ! Okay, let us loop through all main diagonal blocks.
      do i = 1, Rmatrices(n)%nblocksPerRow
      
        ! Get the diagonal block
        p_rblock => Rmatrices(n)%RmatrixBlock(i,i)
        
        ! Is the matrix empty?
        if(p_rblock%NEQ .eq. 0) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! Is the scaling factor zero?
        if(p_rblock%dscaleFactor .eq. 0.0_DP) then
          ! Zero pivot
          ccompatible = LINSOL_COMP_ERRZEROPIVOT
          exit
        end if
        
        ! What about the matrix type?
        if((p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
           (p_rblock%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
          ! Matrix type not supported
          ccompatible = LINSOL_COMP_ERRMATTYPE
          exit
        end if
        
        ! This block is okay, continue with next one
      
      end do
    
    end if
    
    ! Set the compatibility flag only for the maximum level -- this is a
    ! one-level solver acting only there!
    if (present(CcompatibleDetail)) &
        CcompatibleDetail (ubound(CcompatibleDetail,1)) = ccompatible
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine linsol_precDefBlockUpwGS (rsolverNode,rd)
  
  use triangulation
  use spatialdiscretisation
  use element
  use dofmapping
  
!<description>
  ! Applies the Jacobi preconditioner <tex>$ D \in A $<tex> to the defect
  ! vector rd and solves <tex>$ Dd_{new} = d $<tex>.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the Jacobi solver
  type(t_linsolNode), intent(inout),target  :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_vectorBlock), intent(inout)        :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    integer, dimension(:), pointer :: p_KCOL, p_KLD
    integer :: NEL, iel, indof, i, ig, j, jg, iloc, iouter, iblock, ielRef
    integer(I32) :: celement
    real(DP) :: dlocOmega
    real(DP), dimension(:), pointer :: p_Dvector, p_Dmatrix
    real(dp), dimension(:,:), allocatable :: DlocMat
    real(dp), dimension(:), allocatable :: DlocVec, DlocSol
    integer, dimension(:), allocatable :: IdofGlob, IdofGlobRef
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    type(t_elementDistribution), pointer :: p_relemDist
    logical :: bsuccess
    real(dp) :: domega
    
      ! Status reset
      rsolverNode%iresult = 0
      
      ! Get omega and the relaxation parameter
      domega = rsolverNode%domega

      ! Get the matrix
      p_rmatrix => rsolverNode%rsystemMatrix%RmatrixBlock(1,1)
      
      ! Now we have to make some decisions. At first, which matrix
      ! structure do we have?
      select case (p_rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9)
      
        ! Get the omega parameter for the matrix. Do not forget the scaling
        ! factor in the matrix!
        dlocomega = rsolverNode%domega / p_rmatrix%dScaleFactor

        ! Now, which data format do we have? Single or double?
        select case (p_rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Get pointers to the matrix data
          call lsyssc_getbase_Kcol (p_rmatrix,p_KCOL)
          call lsyssc_getbase_Kld (p_rmatrix,p_KLD)
          call lsyssc_getbase_double (p_rmatrix,p_Dmatrix)
          
          ! Take care of the accuracy of the vector
          select case (rd%RvectorBlock(1)%cdataType)
          case (ST_DOUBLE)
            ! Get the data array of the vector
            call lsyssc_getbase_double (rd%RvectorBlock(1),p_Dvector)
            
            ! Get pointers for quicker access
            p_rspatialDiscr => p_rmatrix%p_rspatialDiscrTest
            p_rtriangulation => p_rspatialDiscr%p_rtriangulation
            
            ! Get number of elements
            NEL = p_rtriangulation%NEL
            
            ! Check if we have a uniform discretisation (simplifies a lot)
            if (p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
            
              ! Get the element distribution from the spatial discretisation
              p_relemDist => p_rspatialDiscr%RelementDistr(1)

              ! Get type of element in this distribution
              celement =  p_relemDist%celement

              ! Get number of local DOF
              indof = elem_igetNDofLoc(celement)
              
              ! Allocate space for global DOFs
              allocate(IdofGlob(indof),IdofGlobRef(indof))

              ! Allocate local matrix and vector
              allocate(DlocMat(indof,indof), DlocVec(indof), DlocSol(indof))
            
              ! Loop over all elements in this element distribution
              do iouter = 1, p_relemDist%NEL

                ! Get next element in the queue (the "most upwind" one)
                iel = rsolverNode%p_rsubnodeBlockUpwGS%p_Iqueue(iouter)

                ! Map local to global DOFs
                call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)

                !!! Get entries from the global into the local matrix
                ! Loop over all local lines
                do i = 1, indof
                  ! Get global line
                  ig = IdofGlob(i)

                  ! Loop over all local columns
                  do j = 1, indof
                    ! Get global columns
                    jg = IdofGlob(j)

                    do iloc = p_KLD(ig), p_KLD(ig+1)-1
                      if (jg.eq.p_KCOL(iloc)) exit
                    end do

                    DlocMat(i,j) = p_Dmatrix(iloc)

                  end do
                end do

                !!! Now get local vector
                do i = 1, indof
                  ! Get global line
                  ig = IdofGlob(i)
                   
                  ! Get global entry in local vector
                  DlocVec(i) = dlocOmega * p_Dvector(ig)
                  
                  ! Loop over all blocks left of the blockdiagonal
                  ! That is loop over all upwind elements
                  do iblock = 1, 3
                    
                    ! Get upwind element
                    ielRef = rsolverNode%p_rsubnodeBlockUpwGS%p_IupwElems(iblock,iel)
                    
                    ! If it exists
                    if (ielRef.ne.0) then

                      ! Get global DOFs of upwind element
                      call dof_locGlobMapping(p_rspatialDiscr, ielRef, IdofGlobRef)
                      
!                      ! Loop over all local lines
!                      do i = 1, indof
                          ig = IdofGlob(i)
                        ! Loop over all (upwind) columns
                        do j = 1, indof
                          jg = IdofGlobRef(j)
                          
                          do iloc = p_KLD(ig), p_KLD(ig+1)-1
                             if (jg.eq.p_KCOL(iloc)) then
                               exit
                             endif  
!                             if ((iloc.eq.p_KLD(ig+1)-1).and.(jg.ne.p_KCOL(iloc))) then
!                               write(*,*) 'nicht gefunden, iel, ielREF', iel, ielREF
!                             endif
                          end do

!                          Dtemp(i) = Dtemp(i)+p_DA(iloc)*p_Dsol(jg)
                          DlocVec(i) = DlocVec(i) - p_Dmatrix(iloc)*p_Dvector(jg)
                      
!                        end do ! i lines in element
                           
                      end do ! j columns in RefElement
                     
                    endif ! check if upwind element exists
                  
                  end do ! iblock
                  
!                  ! Substract product with lower diagonal
!                  do jg = p_KLD(ig), p_KLD(ig+1)-1
!                    if (p_KCOL(jg).ge.IdofGlob(1)) exit
!                    DlocVec(i) = DlocVec(i) - p_Dmatrix(jg)*p_Dvector(p_KCOL(jg))
!                  end do
                  
                end do

                !!! Now solve local linear system
                call mprim_invertMatrixDble(DlocMat, DlocVec, DlocSol, &
                                            indof, 2, bsuccess)

                ! Test if local system could be solved
                if (.not.bsuccess) then
                  call output_line ( &
                    'Local matricx singular.', &
                    OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockUpwGS')
                    do i = 1,size(DlocMat,1)
                      write(*,*) DlocMat(i,:)
                    end do
                  call sys_halt()
                end if

                !!! Distribute local solution to global solution
                do i = 1, indof
                  ! Get global line
                  ig = IdofGlob(i)
                   
                  ! Get local entry in global vector
                  p_Dvector(ig) = DlocSol(i)
                end do
                
              end do !iel
              
              ! Deallocate space for global DOFs
              deallocate(IdofGlob, IdofGlobRef)

              ! Deallocate local matrix and vector
              deallocate(DlocMat, DlocVec, DlocSol)
            
            else ! No uniform distribution
            
              call output_line ('Discretisation not uniform - not supported.', &
                              OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockUpwGS')
              call sys_halt()
            
            end if ! Type of discretisation
            
          case default
            call output_line ('Unsupported vector precision.', &
                              OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockUpwGS')
            call sys_halt()
          end select
       
       
        case default
          call output_line ('Unsupported matrix precision.', &
                            OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockUpwGS')
          call sys_halt()
        end select
      
      case default
        call output_line ('Unsupported matrix format.', &
                          OU_CLASS_ERROR, OU_MODE_STD, 'linsol_precBlockUpwGS')
        call sys_halt()
      end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsol_doneBlockUpwGS (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the Block-SOR solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of SOR which is to be cleaned up.
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    deallocate(rsolverNode%p_rsubnodeBlockUpwGS)
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initStructureBlockUpwGS (rsolverNode,ierror,isolverSubgroup)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Allocate list of elements (to be sorted in downwind direction)
    allocate(rsolverNode%p_rsubNodeBlockUpwGS%p_Iqueue &
             (rsolverNode%p_rsubNodeBlockUpwGS%p_rtriangulation%NEL))
    
    ! Allocate list of upwind elements to one element
    allocate(rsolverNode%p_rsubNodeBlockUpwGS%p_IupwElems(3, &
              rsolverNode%p_rsubNodeBlockUpwGS%p_rtriangulation%NEL))
    
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine linsol_initDataBlockUpwGS (rsolverNode, ierror,isolverSubgroup)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
  
  ! local variables for sorting
  integer:: i,j
  integer :: iupwind, idownwind, NEL, NELinqueue, iv, iv1, iv2, iv3
  real(dp) :: dvx, dvy, dnx, dny, dnv

  integer, dimension(:,:), pointer :: p_IelementsAtEdge    
  integer, dimension(:,:), allocatable :: IdownwElems
  integer, dimension(:), allocatable :: Ioutdeg

  integer, dimension(:,:), pointer :: p_IupwElems
  integer, dimension(:), pointer :: p_Iqueue
  type(t_triangulation), pointer :: p_rtriangulation

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...


    ! Now sort the elements in downwind direction
    
    ! Get pointer to triangulation
    p_rtriangulation => rsolverNode%p_rsubNodeBlockUpwGS%p_rtriangulation
    
    ! Get pointer to upwind array
    p_IupwElems => rsolverNode%p_rsubNodeBlockUpwGS%p_IupwElems
    
    ! Get pointer to (soon to be sorted) element list
    p_Iqueue => rsolverNode%p_rsubNodeBlockUpwGS%p_Iqueue
    
    ! Get number of elements
    NEL = p_rtriangulation%NEL
    
    ! Allocate space for array with information which elements are
    ! located downwind and their degree
    allocate(IdownwElems(3,NEL))
    allocate(Ioutdeg(NEL))
    
    ! Get pointer to elements at edge
    call storage_getbase_int2D(p_rtriangulation%h_IelementsAtEdge,&
                               p_IelementsAtEdge)
    
    ! Initialize arrays
    Ioutdeg(:)       = 0
    IdownwElems(:,:) = 0
    p_IupwElems(:,:) = 0
    p_Iqueue(:)      = 0

    ! Get velocity vector
    dvx = rsolverNode%p_rsubNodeBlockUpwGS%Dvel(1)
    dvy = rsolverNode%p_rsubNodeBlockUpwGS%Dvel(2)
    
    ! If velocity = 0, then numbering is arbitrary, but the solver needs
    ! a numbering, so we take an arbitrary velocity
    if (dvx*dvx+dvy*dvy<10.0_dp*SYS_EPSREAL_DP) then
      dvx = 1.0_dp
      dvy = 1.0_dp
    end if
    
    ! Give each element a degree by looping over all edges and increasing the
    ! degree of the downwind element
    do j = 1, p_rtriangulation%NMT
      
      ! Get normal vector
      dnx = rsolverNode%p_rsubNodeBlockUpwGS%p_Dnormals(1,j)
      dny = rsolverNode%p_rsubNodeBlockUpwGS%p_Dnormals(2,j)
      
      ! Normal * velocity vector
      dnv = dvx*dnx + dvy*dny
      
      ! get upwind and downwind element
      if ((dnv).lt.0.0_DP) then
        iupwind   = p_IelementsAtEdge(2,j)
        idownwind = p_IelementsAtEdge(1,j)
      elseif ((dnv).gt.0.0_DP) then
        iupwind   = p_IelementsAtEdge(1,j)
        idownwind = p_IelementsAtEdge(2,j)
      else ! Elements parallel in transport direction
      endif
       
      ! Increase degrees and get upwind/downwind information
      if (abs(dnv).gt.1E-10_DP) then
        
        ! Check if both elements exist
        if((iupwind*idownwind).gt.0) then
        
          ! Increase degree of downwind element
          Ioutdeg(idownwind) = Ioutdeg(idownwind) + 1
        
          ! Write downwind information to free space
          if (IdownwElems(1,iupwind).eq.0) then
            IdownwElems(1,iupwind) = idownwind
          elseif (IdownwElems(2,iupwind).eq.0) then
            IdownwElems(2,iupwind) = idownwind
          else
            IdownwElems(3,iupwind) = idownwind
          endif
          
          ! Write upwind information to free space
          if (p_IupwElems(1,idownwind).eq.0) then
            p_IupwElems(1,idownwind) = iupwind
          elseif (p_IupwElems(2,idownwind).eq.0) then
            p_IupwElems(2,idownwind) = iupwind
          else
            p_IupwElems(3,idownwind) = iupwind
          endif
          
        endif
      
      endif
      
    end do
    
    ! Write all Boundaryelements (deg=0) into Resultqueue
    NELinqueue = 0
    do i = 1, NEL
      ! If degree of element = 0 (downwind to no element = at boundary here)
      if (Ioutdeg(i).eq.0) then
        NELinqueue = NELinqueue + 1
        p_Iqueue(NELinqueue) = i
      endif
    end do
    
    ! Now we do the Resorting, go through the Elements starting at the Boundary
    ! NELinqueue is NOT reset
    do i = 1, NEL
      
      ! Get next element in queue
      iv = p_Iqueue(i)
      
      ! If there is no next element - error
      ! Maybe later add another element with lowest order to includevortex-like situations
      if (iv.eq.0) then
        call output_line ( &
                      'Elements could not be sorted.', &
                      OU_CLASS_ERROR, OU_MODE_STD, 'linsol_initDataBlockUpwGS')
        call sys_halt()
      endif
      
      ! Get all downwind neighbours
      iv1 = IdownwElems(1,Iv)
      iv2 = IdownwElems(2,Iv)
      iv3 = IdownwElems(3,Iv)
      
      ! Now - for all downwind neighbours:
      ! 1. Decrease their degree
      ! 2. If their degree is == 0
      !    That is all upwind neighbours of the element are alreadvy in the queue
      !    Add it to the queue
      if (iv1.gt.0) then
        Ioutdeg(iv1) = Ioutdeg(iv1)-1
        if(Ioutdeg(iv1).eq.0) then
          NELinqueue = NELinqueue + 1
          p_Iqueue(NELinqueue) = iv1
        endif
      endif
      
      if (iv2.gt.0) then
        Ioutdeg(iv2) = Ioutdeg(iv2)-1
        if(Ioutdeg(iv2).eq.0) then
          NELinqueue = NELinqueue + 1
          p_Iqueue(NELinqueue) = Iv2
        endif
      endif
      
      if (iv3.gt.0) then
        Ioutdeg(iv3)=Ioutdeg(iv3)-1
        if (Ioutdeg(iv3).eq.0) then
          NELinqueue = NELinqueue + 1
          p_Iqueue(NELinqueue) = iv3
        endif
      endif
      
    end do
    
    ! Deallocate arrays
    deallocate(IdownwElems,Ioutdeg)
    
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneDataBlockUpwGS (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine linsol_doneStructureBlockUpwGS (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_linsolNode), intent(inout)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(in)                    :: isolverSubgroup
!</input>

!</subroutine>

  ! local variables
  integer :: isubgroup
    
    ! by default, initialise solver subroup 0
    isubgroup = 0
    if (present(isolversubgroup)) isubgroup = isolverSubgroup

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there is not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Deallocate list of elements (to be sorted in downwind direction)
    if (associated(rsolverNode%p_rsubNodeBlockUpwGS%p_Iqueue)) then
      deallocate(rsolverNode%p_rsubNodeBlockUpwGS%p_Iqueue)
      nullify(rsolverNode%p_rsubNodeBlockUpwGS%p_Iqueue)
    end if
    
    ! Deallocate list of upwind elements to one element
    if (associated(rsolverNode%p_rsubNodeBlockUpwGS%p_IupwElems)) then
      deallocate(rsolverNode%p_rsubNodeBlockUpwGS%p_IupwElems)
      nullify(rsolverNode%p_rsubNodeBlockUpwGS%p_IupwElems)
    end if
      
  end subroutine

end module