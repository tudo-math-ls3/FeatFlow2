!##############################################################################
!# ****************************************************************************
!# <name> spacepreconditioner </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the routines to solve the stationary core equation
!# of the problem with a nonlinear solver. The core equation aims to solve
!# the coupled KKT system of the minimisation problem
!#
!#   $$J(y)  :=  1/2 ||y-z||_{L_2}  +  \gamma/2 ||y(T)-z(T)||_{L_2}  +  \alpha_C||u||^2  ->  min! $$
!#
!# with $z$ being a given 'desired' flow field in the domain and
!# u being an unknown control.
!#
!# The discretised core equation reads at the moment:
!#
!#  $$        A_1 y   +  \eta_1 B p   +  \mu_1 M \lambda         = f_1 $$
!#  $$ \tau_1 B^T y   +  \kappa_1 I p                            = f_2 $$
!#
!#  $$        R_2 y   +  A_2 \lambda  + \eta_2   B \xi           = f_3 $$
!#  $$            \tau_2 B^T \lambda  + \kappa_2 I \xi           = f_4 $$
!#
!# with
!#
!#   $$ A_1 = \iota_1 I  +  \alpha_1 M  +  \theta_1 L  +  \gamma_1 N(y) + dnewton_1 N*(y)$$
!#   $$ A_2 = \iota_2 I  +  \alpha_2 M  +  \theta_2 L  +  \gamma_2 N(y) + dnewton_2 N*(y)$$
!#   $$ R_2 =               \mu_2    M  +                                 dr_2      N*(\lambda) $$
!#
!# and
!#
!#   $I$     = identity matrix,
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity including stabilisation, depending on the
!#             primal velocity, i.e.
!#                   $$ (y\Delta)\cdot $$
!#   $N*(y)$ = Newton matrix, depending on the primal velocity, i.e.
!#                  $$ (\Delta y)\cdot $$
!#   $N*(\lambda)$ = Newton matrix, depending on the dual velocity, i.e.
!#                  $$ (\Delta \lambda)\cdot $$
!#
!#   $\iota_i$  = 0/1     - switches the identity matrix on/off,
!#   $\alpha_i$ = 0/1     - switches the mass matrix on/off;
!#                          =0 for stationary problem,
!#   $\theta_í$           - weight for the Laplace matrix,
!#   $\gamma_i$ = 0/1     - Switches the nonlinearity on/off;
!#                          =0 for Stokes system,
!#   $dnewton_i \in R$    - Switches the Newton matrix on/off.
!#   $\eta_i$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau_i$   = 0/1     - Switches the 'B^T'-term on/off,
!#   $\mu_i$              - Weight for the 'coupling' mass matrix.
!#   $\kappa_i$ = 0/1     - Switches of the identity matrix I for the pressure
!#                          in the continuity equation
!#   $\dr_i \in R$        - Switches the 'reactive coupling mass matrix' on/off
!#
!# (y,p) is the velocity/pressure solution pair.
!# (lambda,xi) is the dual velocity/pressure solution.
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
!# 1.) cc_createSpacePreconditioner / cc_releasePreconditioner
!#     -> Creates/Releases a basic preconditioner structure for a spatial
!#        preconditioner.
!#
!# 2.) cc_precondSpaceDefect
!#     -> Executes spatial preconditioning on a given defect vector
!#
!# </purpose>
!##############################################################################

module spacepreconditioner

  use fsystem
  use storage
  use boundary
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
  use filtersupport
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use matrixmodification
  use linearsolver
  use linearsolverautoinitialise
  use matrixrestriction
  use trilinearformevaluation
  
  use collection
  use convection
    
  use user_callback
  
  use matrixio
  use vectorio
  
  use optcanalysis
  use spacematvecassembly
  
  implicit none
  
!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  integer, parameter :: CCPREC_NONE         = -1

  ! Preconditioning with inverse mass matrix (not yet implemented)
  integer, parameter :: CCPREC_INVERSEMASS   = 0

  ! Preconditioning by linear solver, solving the linearised system
  integer, parameter :: CCPREC_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  integer, parameter :: CCPREC_NEWTON        = 2
  
  ! Preconditioning by inexact/adaptive Newton iteration
  integer, parameter :: CCPREC_INEXACTNEWTON = 3

!</constantblock>

!</constants>

  
!<types>
  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it's noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It's more or less
  ! a collection if different flags.
  type t_ccPreconditionerSpecials
  
    ! Whether to use 'adaptive matrices', i.e. set up coarse grid matrices
    ! with the help of fine grid matrices. This is used for very special
    ! discretisations only (e.g. Q1~/Q0). =0: deactivate
    integer :: iadaptiveMatrices    = 0
    
    ! A configuration parameter for adaptive matrices.
    real(DP) :: dadMatThreshold     = 0.0_DP

    ! If the preconditioner is a linear solver:
    ! Type of solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Multigrid solver
    integer :: isolverType = 0
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of smoother.
    ! =0: general VANKA (slow, but independent of the discretisation and of the problem)
    ! =1: general VANKA; 'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 0, but slightly faster)
    ! =2: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =3: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 8, but faster)
    ! =4: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =5: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 10, but faster)
    integer :: ismootherType = 3
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of coarse grid solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Defect correction with diagonal VANKA preconditioning.
    ! =2: BiCGStab with diagonal VANKA preconditioning
    integer :: icoarseGridSolverType = 1
        
    ! This flag is set to .TRUE. if there are no Neumann boundary
    ! components. In that case, the pressure matrices of direct
    ! solvers must be changed.
    logical :: bneedPressureDiagonalBlock = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on all levels except for the coarse mesh.
    logical :: bneedVirtTransposedD = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on the coarse mesh.
    logical :: bneedVirtTransposedDonCoarse = .false.
    
  end type

!</typeblock>

!<typeblock>

  ! Represents the core equation on one level of the discretisation.
  ! Collects all information that are necessary to assemble the
  ! (linearised) system matrix and RHS vector.
  type t_cccoreEquationOneLevel
  
    ! The (linearised) system matrix for that specific level.
    type(t_matrixBlock), pointer :: p_rmatrix => null()
    
    ! Reference to the discretisation structure.
    type(t_blockDiscretisation), pointer :: p_rdiscrBlock => null()

    ! Reference to the static matrices on this level (Stokes, B,...)
    type(t_staticLevelInfo), pointer :: p_rstaticInfo

    ! Temporary vectors for the interpolation of a solution to a lower level.
    ! Exists only on levels NLMIN..NLMAX-1 !
    type(t_vectorBlock) :: rtempVector1
    type(t_vectorBlock) :: rtempVector2
    type(t_vectorBlock) :: rtempVector3
    
    ! Bondary conditions
    type(t_discreteBC) :: rdiscreteBC
    
    ! Fictitious boundary conditions
    type(t_discreteFBC) :: rdiscreteFBC

    ! Block matrix, which is used in the defect correction / Newton
    ! algorithm as preconditioner matrix of the correspnding underlying
    ! linear sytem. Is usually the (linearise) system matrix or
    ! a Newton matrix. This matrix is changed during the
    ! nonlinear iteration and used e.g. if a linear solver (Multigrid) is
    ! used for preconditioning.
    type(t_matrixBlock) :: rmatrixPreconditioner
    
    ! Flag that specifies if there are Neumann BC`s on this level.
    logical :: bhasNeumann
    
  end type

!</typeblock>

!<typeblock>

!<typeblock>

  ! Preconditioner structure for CCxD. This structure saves the configuration of the
  ! spatial preconditioner.
  
  type t_ccspatialPreconditioner
  
    ! Type of preconditioner.
    ! This is one of the CCPREC_xxxx flags as defined above (CCPREC_INVERSEMASS for
    ! preconditioning with inverse mass matrix, CCPREC_LINEARSOLVER for solving a linear
    ! system, CCPREC_NEWTON / CCPREC_INEXACTNEWTON for a Newton iteration,...)
    integer :: ctypePreconditioning = CCPREC_NONE
    
    ! Name of the section in the DAT file configuring this preconditioner.
    character(LEN=SYS_STRLEN) :: spreconditionerSection = ''
    
    ! Minimum discretisation level
    integer :: nlmin = 0
    
    ! Maximum discretisation level
    integer :: nlmax = 0
    
    ! A t_ccPreconditionerSpecials structure that saves information about
    ! special 'tweaks' in matrices such that everything works.
    type(t_ccPreconditionerSpecials) :: rprecSpecials
    
    ! An array of t_cccoreEquationOneLevel structures for all levels
    ! of the discretisation.
    type(t_cccoreEquationOneLevel), dimension(:), pointer :: RcoreEquation => null()

    ! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    type(t_linsolNode), pointer :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), dimension(:), pointer :: p_Rprojection
    
    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar) :: rtempVectorSc

    ! A filter chain that is used for implementing boundary conditions or other
    ! things when invoking the linear solver.
    type(t_filterChain), dimension(4) :: RfilterChain
    
  end type

!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

  !<subroutine>
  
    subroutine cc_preparePrecondMatrixAssembly (rnonlinearSpatialMatrix,&
        ilev,nlmin,nlmax,rprecSpecials)
  
  !<description>
    ! Prepares a rnonlinearCCMatrix structure for the assembly according
    ! to a preconditioner. rprecSpecials specifies a couple of preconditioner
    ! flags that configure the shape of the system matrix that the preconditioner
    ! needs.
    !
    ! cc_initNonlinMatrix must have been called prior to this routine to
    ! initialise the basic matrix. cc_preparePrecondMatrixAssembly will then
    ! add assembly-specific parameters of the preconditioner.
  !</description>

  !<input>
    ! Current assembly level.
    integer, intent(in) :: ilev
    
    ! Minimum assembly level.
    integer, intent(in) :: nlmin
    
    ! Maximum assembly level.
    integer, intent(in) :: nlmax
  
    ! Structure with assembly-specific parameters of the preconditioner.
    type(t_ccPreconditionerSpecials), intent(in) :: rprecSpecials
  !</input>
  
  !<inputoutput>
    ! Nonlinear matrix structure.
    ! Basic parameters in this structure are filled with data.
    type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
  !</inputoutput>
               
  !</subroutine>
      
      ! Parameters for adaptive matrices for Q1~ with anisotropic elements
      rnonlinearSpatialMatrix%iadaptiveMatrices = rprecSpecials%iadaptiveMatrices
      rnonlinearSpatialMatrix%dadmatthreshold = rprecSpecials%dadmatthreshold
      
      ! Depending on the level, we have to set up information about
      ! transposing B-matrices.
      if (ilev .eq. nlmin) then
        rnonlinearSpatialMatrix%bvirtualTransposedD = rprecSpecials%bneedVirtTransposedDonCoarse
      else
        rnonlinearSpatialMatrix%bvirtualTransposedD = rprecSpecials%bneedVirtTransposedD
      end if
      
    end subroutine

  ! ***************************************************************************
  ! Routines to create a nonlinear iteration structure, to save it
  ! to a collection, to rebuild it from there and to clean it up.
  ! ***************************************************************************

!<subroutine>

  subroutine cc_createSpacePreconditioner (rpreconditioner,NLMIN,NLMAX)
  
!<description>
  ! This routine creates a spational preconditioner structure. The structure is
  ! initialised to handle NLMAX-NLMIN+1 discretisation levels.
!</description>

!<input>
  ! Minimum discretisation level to be maintained
  integer, intent(IN) :: NLMIN
  
  ! Maximum discretisation level to be maintained. The maximum level coincides
  ! with the level where to solve the system.
  integer, intent(IN) :: NLMAX
!</input>

!<output>
  ! A spatial preconditioner structure to be initialised.
  type(t_ccspatialPreconditioner), intent(OUT) :: rpreconditioner
!</output>

!</subroutine>

    rpreconditioner%NLMIN = NLMIN
    rpreconditioner%NLMAX = NLMAX

    ! Initialise the matrix pointers on all levels that we have to maintain.
    allocate(rpreconditioner%RcoreEquation(NLMIN:NLMAX))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_releaseSpacePreconditioner (rpreconditioner)
  
!<description>
  ! Releases allocated memory in the spatial preconditioner structure.
!</description>

!<inputoutput>
  ! The spatial preconditioner structure that should be cleaned up.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>
    
    if (associated(rpreconditioner%RcoreEquation)) &
      deallocate(rpreconditioner%RcoreEquation)

    rpreconditioner%NLMIN = 0
    rpreconditioner%NLMAX = 0

  end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine cc_precondSpaceDefect (rpreconditioner,rnonlinearSpatialMatrix,&
        rd,rx1,rx2,rx3,bsuccess,rcollection)
  
    use linearsystemblock
    use collection
    
  !<description>
    ! Defect preconditioning routine. Based on the current iteration
    ! vector rx, this routine has to perform preconditioning on the defect
    ! vector rd.
    !
    ! Note that no boundary conditions have to be attached to rd, rx1, rx2 and rx3.
    ! The preconditioner internally applies the global boundary conditions of the
    ! problem where necessary.
  !</description>

  !<input>
    ! Configuration of the core equation on the maximum level.
    type(t_nonlinearSpatialMatrix), intent(IN)      :: rnonlinearSpatialMatrix

  !</input>

  !<inputoutput>
    ! Spatial preconditioner structure that defines all parameters how to perform
    ! preconditioning.
    type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner

    ! Defect vector b-A(rx)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(INOUT)            :: rd

    ! Ccollection structure of the application.
    type(t_collection)                            :: rcollection
    
    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Iteration vector of the 'previous' timestep. Can be undefined if there
    ! is no previous timestep.
    type(t_vectorBlock), intent(IN), target       :: rx1

    ! Iteration vector of the 'current' timestep.
    type(t_vectorBlock), intent(IN), target       :: rx2

    ! Iteration vector of the 'next' timestep. Can be undefined if there
    ! is no next timestep.
    type(t_vectorBlock), intent(IN), target       :: rx3
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode
    
    integer :: i
    logical :: bassembleNewton
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    

    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata
    type(t_vectorBlock) :: rxtemp,rdtemp
    integer, dimension(:), pointer :: p_IelementList

    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    type(t_matrixblock) :: rpremat
    type(t_linsolNode), pointer :: p_rsolver
    call lsysbl_getbase_double (rd,p_def)
    call lsysbl_getbase_double (rx2,p_vec)

      select case (rpreconditioner%ctypePreconditioning)
      case (CCPREC_NONE)
        ! No preconditioning. Do nothing.
      case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_INEXACTNEWTON)
        ! Preconditioning with a linear solver.
        !
        ! At first, assemble the preconditioner matrices on all levels
        ! and incorporate all boundary conditions.
        !
        ! Should we assemble the Newton part?
        
        bassembleNewton = .false.
        
        if ((rpreconditioner%ctypePreconditioning .eq. CCPREC_NEWTON) .or. &
            (rpreconditioner%ctypePreconditioning .eq. CCPREC_INEXACTNEWTON)) then
            
          ! Use Newton in any case.
          bassembleNewton = .true.
          
        end if
        
        ! Assemble the preconditioner matrices in rpreconditioner
        ! on all levels that the solver uses.
        call assembleLinsolMatrices (rpreconditioner,rnonlinearSpatialMatrix,&
            rcollection,bassembleNewton,rx1,rx2,rx3)
          
        ! Our 'parent' (the caller of the nonlinear solver) has prepared
        ! a preconditioner node for us (a linear solver with symbolically
        ! factorised matrices). Get this from the collection.
      
        p_rsolverNode => rpreconditioner%p_rsolverNode

        ! Re-attach the system matrices to the solver.
        ! Note that no pointers and no handles are changed, so we can savely do
        ! that without calling linsol_doneStructure/linsol_doneStructure.
        ! This simply informs the solver about possible new scaling factors
        ! in the matrices in case they have changed...
        allocate(Rmatrices(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
        
        ! Attach the system matrix
        do i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
          call lsysbl_duplicateMatrix ( &
            rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end do
          
        ! DEBUG!!!
        !CALL matio_writeBlockMatrixHR (Rmatrices(rpreconditioner%NLMIN), 'matrix',&
        !                               .TRUE., 0, 'matrixstat.txt','(E10.2)')
        
        call linsol_setMatrices(rpreconditioner%p_rsolverNode,Rmatrices(:))
        
        ! DEBUG!!!
        !DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        !  CALL storage_getbase_double (Rmatrices(i)% &
        !      RmatrixBlock(4,1)%h_Da,p_Ddata)
        !END DO
        
        ! Initialise data of the solver. This in fact performs a numeric
        ! factorisation of the matrices in UMFPACK-like solvers.
        call linsol_updateStructure (rpreconditioner%p_rsolverNode,ierror)
        call linsol_initData (p_rsolverNode, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) then
          print *,'linsol_initData failed!'
          call sys_halt()
        end if
        
        ! The solver got the matrices; clean up Rmatrices, it was only of temporary
        ! nature...
        do i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)

        ! Solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        call linsol_precondDefect (p_rsolverNode,rd)

        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
        ! for the next iteration.
        call linsol_doneData (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .eq. 0
        
        if (bsuccess) then
          ! Filter the final defect
          call filter_applyFilterChainVec (rd, rpreconditioner%RfilterChain)
        end if
        
        if (p_rsolverNode%dfinalDefect .gt. p_rsolverNode%dinitialDefect*0.99_DP) then
          ! Ignore the correction, it cannot be good enough!
          call output_line (&
            'Space-Time-Preconditioner: Warning. Solution ignored for missing accuracy.')
            
          call lsysbl_clearVector (rd)
        end if
        
      end select
      
    contains
      
      subroutine assembleLinsolMatrices (rpreconditioner,rnonlinearSpatialMatrix,rcollection,&
          bassembleNewton,rx1,rx2,rx3)

      use linearsystemblock
      use collection

      ! Assembles on every level a matrix for the linear-solver/Newton preconditioner.
      ! bnewton allows to specify whether the Newton matrix or only the standard
      ! system matrix is evaluated. The output is written to the p_rpreconditioner
      ! matrices specified in the rnonlinearIteration structure.

      ! Spatial preconditioner structure that defines all parameters how to perform
      ! preconditioning.
      type(t_ccspatialPreconditioner), intent(IN),target    :: rpreconditioner

      ! Level independent configuration of the core equation
      type(t_nonlinearSpatialMatrix), intent(IN)      :: rnonlinearSpatialMatrix

      ! Reference to a collection structure that contains all parameters of the
      ! discretisation (for nonlinearity, etc.).
      type(t_collection), intent(INOUT)                :: rcollection

      ! TRUE  = Assemble the Newton preconditioner.
      ! FALSE = Assemble the standard defect correction preconditioner
      !         (i.e. the linearised system matrix).
      logical, intent(IN) :: bassembleNewton
      
      ! Current iteration vector of the 'previous' timestep. May be undefined
      ! if there is no previous timestep.
      type(t_vectorBlock), intent(IN), target          :: rx1

      ! Current iteration vector.
      type(t_vectorBlock), intent(IN), target          :: rx2

      ! Current iteration vector of the 'next' timestep. May be undefined
      ! if there is no previous timestep.
      type(t_vectorBlock), intent(IN), target          :: rx3

      ! local variables
      real(DP) :: dnewton
      integer :: ilev
      type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
      type(t_vectorScalar), pointer :: p_rvectorTemp
      type(t_vectorBlock), pointer :: p_rvectorFine1,p_rvectorFine2,p_rvectorFine3
      type(t_vectorBlock), pointer :: p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3
      type(t_interlevelProjectionBlock), pointer :: p_rprojection
      type(t_nonlinearSpatialMatrix) :: rlocalNonlSpatialMatrix
      integer, dimension(1), parameter :: Irows = (/1/)

      ! DEBUG!!!
      type(t_cccoreEquationOneLevel), pointer :: p_rcore
    !    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !    call lsysbl_getbase_double (rd,p_def)
    !    call lsysbl_getbase_double (rx,p_vec)

        ! Get the interlevel projection structure and the temporary vector
        ! from the collection.
        ! Our 'parent' prepared there how to interpolate the solution on the
        ! fine grid to coarser grids.
        p_rvectorTemp => rpreconditioner%rtempVectorSc

        ! Initialise the matrix assembly structure rlocalNonlSpatialMatrix to describe the
        ! matrix we want to have. We have to initialise the adaptivity constants,
        ! which are not part of the standard initialisation.
        rlocalNonlSpatialMatrix = rnonlinearSpatialMatrix
        
        ! On all levels, we have to set up the nonlinear system matrix,
        ! so that the linear solver can be applied to it.
        
        nullify(p_rmatrix)

        do ilev=rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
        
          ! Get the matrix on the current level.
          ! Shift the previous matrix to the pointer of the fine grid matrix.
          p_rmatrixFine => p_rmatrix
          p_rmatrix => rpreconditioner%RcoreEquation(ilev)%rmatrixPreconditioner
        
          ! On the highest level, we use rx as solution to build the nonlinear
          ! matrix. On lower levels, we have to create a solution
          ! on that level from a fine-grid solution before we can use
          ! it to build the matrix!
          if (ilev .eq. rpreconditioner%NLMAX) then
          
            p_rvectorCoarse1 => rx1
            p_rvectorCoarse2 => rx2
            p_rvectorCoarse3 => rx3
            
          else
            ! We have to discretise a level hierarchy and are on a level < NLMAX.

            ! Get the mujltilevel projection structure that describes the
            ! projection from the finer level to the current one.
            p_rprojection => rpreconditioner%p_Rprojection(ilev+1)
            
            ! Get the temporary vector on level i. Will receive the solution
            ! vector on that level.
            p_rvectorCoarse1 => rpreconditioner%RcoreEquation(ilev)%rtempVector1
            p_rvectorCoarse2 => rpreconditioner%RcoreEquation(ilev)%rtempVector2
            p_rvectorCoarse3 => rpreconditioner%RcoreEquation(ilev)%rtempVector3
            
            ! Get the solution vector on level i+1. This is either the temporary
            ! vector on that level, or the solution vector on the maximum level.
            if (ilev .lt. rpreconditioner%NLMAX-1) then
              p_rvectorFine1 => rpreconditioner%RcoreEquation(ilev+1)%rtempVector1
              p_rvectorFine2 => rpreconditioner%RcoreEquation(ilev+1)%rtempVector2
              p_rvectorFine3 => rpreconditioner%RcoreEquation(ilev+1)%rtempVector3
            else
              p_rvectorFine1 => rx1
              p_rvectorFine2 => rx2
              p_rvectorFine3 => rx3
            end if

            ! Interpolate the solution from the finer grid to the coarser grid.
            ! The interpolation is configured in the interlevel projection
            ! structure we got from the collection.
            call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse1, &
                                             p_rvectorFine1,p_rvectorTemp)
            call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse2, &
                                             p_rvectorFine2,p_rvectorTemp)
            call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse3, &
                                             p_rvectorFine3,p_rvectorTemp)

            ! Apply the filter chain to the temp vector.
            ! This implements the boundary conditions that are attached to it.
            ! NOTE: Deactivated for standard CC2D compatibility -- and because
            ! it has to be checked whether the correct boundary conditions
            ! are attached to that vector!
            ! CALL filter_applyFilterChainVec (p_rvectorCoarse, p_RfilterChain)

          end if

          ! Set the pointers in the rlocalNonlSpatialMatrix structure according
          ! to the current level.
          p_rcore => rpreconditioner%RcoreEquation(ilev)

          ! Adjust the nonlinear matrix to our current level.
          rlocalNonlSpatialMatrix%p_rdiscretisation     => &
              p_rmatrix%p_rblockDiscrTrial

          rlocalNonlSpatialMatrix%p_rstaticInfo         => &
              rpreconditioner%RcoreEquation(ilev)%p_rstaticInfo

          call cc_preparePrecondMatrixAssembly (rlocalNonlSpatialMatrix,&
              ilev,rpreconditioner%NLMIN,&
              rpreconditioner%NLMAX,rpreconditioner%rprecSpecials)

          ! Assemble the matrix.
          ! If we are on a lower level, we can specify a 'fine-grid' matrix.
          if (ilev .eq. rpreconditioner%NLMAX) then
            call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rlocalNonlSpatialMatrix,&
                p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3)
          else
            call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rlocalNonlSpatialMatrix,&
                p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3,&
                p_rmatrixFine)
          end if

          ! Boundary conditions
          ! ---------------------------------------------------

          if (associated(p_RfilterChain)) then
            ! Apply the filter chain to the matrix.
            ! As the filter consists only of an implementation filter for
            ! boundary conditions, this implements the boundary conditions
            ! into the system matrix.
            call filter_applyFilterChainMat (p_rmatrix, p_RfilterChain)
          else
            ! Call the matrix filter for the boundary conditions to include the BC's
            ! into the matrix.
            call matfil_discreteBC (p_rmatrix)
            call matfil_discreteFBC (p_rmatrix)
          end if
            
          ! 'Nonlinear' boundary conditions like slip boundary conditions
          ! are not implemented with a filter chain into a matrix.
          ! Call the appropriate matrix filter of 'nonlinear' boundary
          ! conditions manually:
          call matfil_discreteNLSlipBC (p_rmatrix,.true.)
            
          ! DEBUG!!!
          !CALL matio_writeBlockMatrixHR (p_rmatrix, 'matrix',&
          !                              .TRUE., 0, 'matrix.txt','(E20.10)')

        end do
        
        if (rpreconditioner%rprecSpecials%bneedPressureDiagonalBlock) then
          
          ! The 3,3-matrix must exist! This is ensured by the initialisation routine.
          !
          ! We have a pure Dirichlet problem. This may give us some difficulties
          ! in the case, the preconditioner uses a direct solver (UMFPACK).
          ! In this case, we have to include a unit vector to the pressure
          ! matrix to make the problem definite!
          if (rpreconditioner%rprecSpecials%isolverType .eq. 0) then
            p_rmatrix => rpreconditioner%RcoreEquation(rpreconditioner%NLMAX)%&
                rmatrixPreconditioner
            
            ! Include a unit vector to the matrix part of the pressure in
            ! the primal equation -- as long as there is not a full identity
            ! matrix in the pressure matrix (what would be the case for
            ! the initial condition).
            if (rlocalNonlSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
              ! Switch the pressure matrix on and clear it; we don't know what is inside.
              p_rmatrix%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
              call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
              call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,4),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,5),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,6),Irows)
            end if

            ! Also in the dual equation, as the BC type coincides
            if (rlocalNonlSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
              ! Switch the pressure matrix on and clear it; we don't know what is inside.
              p_rmatrix%RmatrixBlock(6,6)%dscaleFactor = 1.0_DP
              call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(6,6))
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,1),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,2),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,3),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,4),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,5),Irows)
              call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(6,6),Irows)
            end if
            
          end if
          
          if (rpreconditioner%rprecSpecials%isolverType .eq. 1) then
          
            ! If we have a MG solver, We also check the coarse grid solver for
            ! the same thing!
            ! What we don't check is the smoother, thus we assume that smoothers
            ! are always solvers that allow the applicance of a filter chain.
            if (rpreconditioner%rprecSpecials%icoarseGridSolverType .eq. 0) then
              p_rmatrix => rpreconditioner%RcoreEquation(rpreconditioner%NLMIN)%&
                  rmatrixPreconditioner
              
              ! Include a unit vector to the matrix part of the pressure in
              ! the primal equation -- as long as there is not a full identity
              ! matrix in the pressure matrix (what would be the case for
              ! the initial condition).
              if (rlocalNonlSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
                ! Switch the pressure matrix on and clear it; we don't know what is inside.
                p_rmatrix%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
                call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
                call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,4),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,5),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,6),Irows)
              end if

              ! Also in the dual equation, as the BC type coincides
              if (rlocalNonlSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
                ! Switch the pressure matrix on and clear it; we don't know what is inside.
                p_rmatrix%RmatrixBlock(6,6)%dscaleFactor = 1.0_DP
                call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(6,6))
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,1),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,2),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,3),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,4),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,5),Irows)
                call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(6,6),Irows)
              end if
              
            end if
            
          end if
            
        end if
        
      end subroutine
      
    end subroutine

end module
