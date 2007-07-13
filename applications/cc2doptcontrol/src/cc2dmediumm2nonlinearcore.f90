!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2nonlinearcore </name>
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
!#  $$   \mu_2 M y    +  A_2 \lambda  +  \eta_2 B \xi            = f_3 $$
!#  $$            \tau_2 B^T \lambda  + \kappa_2 I \xi           = f_4 $$
!#
!# with
!#
!#   $$ A_1 = \iota_1 I  +  \alpha_1 M  +  \theta_1 L  +  \gamma_1 N(y) + dnewton_1 N*(y)$$
!#   $$ A_2 = \iota_2 I  +  \alpha_2 M  +  \theta_2 L  +  \gamma_2 N(y) + dnewton_2 N*(y)$$
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
!#   
!#   $\iota_i$  = 0/1     - switches the identity matrix on/off,
!#   $\alpha_i$ = 0/1     - switches the mass matrix on/off;
!#                          =0 for stationary problem,
!#   $\theta_�$           - weight for the Laplace matrix,
!#   $\gamma_i$ = 0/1     - Switches the nonlinearity on/off;
!#                          =0 for Stokes system,
!#   $dnewton_i \in R$    - Switches the Newton matrix on/off.
!#   $\eta_i$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau_i$   = 0/1     - Switches the 'B^T'-term on/off,
!#   $\mu_i$              - Weight for the 'coupling' mass matrix.
!#   $\kappa_i$ = 0/1     - Switches of the identity matrix I for the pressure
!#                          in the continuity equation
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
  
  USE matrixio
  USE vectorio
  
  USE cc2dmediumm2optcanalysis
  USE cc2dmediumm2matvecassembly
  
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

  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it's noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It's more or less
  ! a collection if different flags.
  TYPE t_ccFinalAssemblyInfo
    ! This flag is set to YES if the B matrices must be assembled
    ! transposedly. This may be necessary for special VANCA type smoothers
    ! if a linear solver is used as preconditioner.
    INTEGER :: iBmatricesTransposed = NO
    
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
    
    ! Temporary vector for the interpolation of a solution to a lower level.
    ! Exists only on levels NLMIN..NLMAX-1 !
    TYPE(t_vectorBlock), POINTER :: p_rtempVector => NULL()

    ! Identty matrix for the pressure
    TYPE(t_matrixScalar), POINTER :: p_rmatrixIdentityPressure => NULL()

    ! Block matrix, which is used in the defect correction / Newton
    ! algorithm as preconditioner matrix of the correspnding underlying
    ! linear sytem. Is usually the (linearise) system matrix or
    ! a Newton matrix. This matrix is changed during the
    ! nonlinear iteration and used e.g. if a linear solver (Multigrid) is
    ! used for preconditioning.
    TYPE(t_matrixBlock), POINTER :: p_rmatrixPreconditioner => NULL()
  
  END TYPE

!</typeblock>

!<typeblock>

!<typeblock>

  ! Preconditioner structure for CCxD. This structure saves the configuration of the
  ! spatial preconditioner.
  
  TYPE t_ccspatialPreconditioner
  
    ! Type of preconditioner.
    ! This is one of the CCPREC_xxxx flags as defined above (CCPREC_INVERSEMASS for
    ! preconditioning with inverse mass matrix, CCPREC_LINEARSOLVER for solving a linear
    ! system, CCPREC_NEWTON for a Newton iteration,...)
    INTEGER :: ctypePreconditioning = CCPREC_NONE
    
    ! Minimum discretisation level
    INTEGER :: NLMIN = 0
    
    ! Maximum discretisation level
    INTEGER :: NLMAX = 0
    
    ! A t_ccFinalAssemblyInfo structure that saves information about
    ! special 'tweaks' in matrices such that everything works.
    TYPE(t_ccFinalAssemblyInfo) :: rfinalAssembly
    
    ! An array of t_cccoreEquationOneLevel structures for all levels
    ! of the discretisation.
    TYPE(t_cccoreEquationOneLevel), DIMENSION(:), POINTER :: RcoreEquation => NULL()

    ! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock), POINTER :: p_rprojection
    
    ! Configuration block for the adaptive Newton preconditioner.
    ! Is only valid if ctypePreconditioning=CCPREC_NEWTONDYNAMIC!
    TYPE(t_ccDynamicNewtonControl) :: radaptiveNewton

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    TYPE(t_vectorScalar), POINTER :: p_rtempVectorSc

    ! A filter chain that is used for implementing boundary conditions or other
    ! things when invoking the linear solver.
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
  END TYPE

!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************
  ! Routines to create a nonlinear iteration structure, to save it
  ! to a collection, to rebuild it from there and to clean it up.
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_createPreconditioner (rpreconditioner,NLMIN,NLMAX)
  
!<description>
  ! This routine creates a spational preconditioner structure. The structure is
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
  ! A spatial preconditioner structure to be initialised.
  TYPE(t_ccspatialPreconditioner), INTENT(OUT) :: rpreconditioner
!</output>

!</subroutine>

    rpreconditioner%NLMIN = NLMIN
    rpreconditioner%NLMAX = NLMAX

    ! Initialise the matrix pointers on all levels that we have to maintain.
    ALLOCATE(rpreconditioner%RcoreEquation(NLMIN:NLMAX))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_releasePreconditioner (rpreconditioner)
  
!<description>
  ! Releases allocated memory in the spatial preconditioner structure.
!</description>

!<inputoutput>
  ! The spatial preconditioner structure that should be cleaned up.
  TYPE(t_ccspatialPreconditioner), INTENT(OUT) :: rpreconditioner
!</inputoutput>

!</subroutine>
    
    IF (ASSOCIATED(rpreconditioner%RcoreEquation)) &
      DEALLOCATE(rpreconditioner%RcoreEquation)

    rpreconditioner%NLMIN = 0
    rpreconditioner%NLMAX = 0

  END SUBROUTINE

  ! ***************************************************************************

  !<subroutine>

    SUBROUTINE c2d2_precondDefect (rpreconditioner,rmatrixComponents,&
        rd,rx,bsuccess,rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! Defect preconditioning routine. Based on the current iteration 
    ! vector rx, this routine has to perform
    ! preconditioning on the defect vector rd. The routine accepts a pointer 
    ! a collection structure rcollection, which allows the routine to access 
    ! information from the main application (e.g. system matrices).
  !</description>

  !<input>
    ! Configuration of the core equation on the maximum level.
    TYPE(t_ccmatrixComponents), INTENT(IN)      :: rmatrixComponents

  !</input>

  !<inputoutput>
    ! Spatial preconditioner structure that defines all parameters how to perform
    ! preconditioning.
    TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd

    ! Ccollection structure of the application.
    TYPE(t_collection)                            :: rcollection
    
    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    LOGICAL, INTENT(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    INTEGER :: ierror
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    INTEGER :: i
    LOGICAL :: bassembleNewton
    TYPE(t_matrixBlock), DIMENSION(:), ALLOCATABLE :: Rmatrices
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    call lsysbl_getbase_double (rd,p_def)
    call lsysbl_getbase_double (rx,p_vec)
!    NLMAX = collct_getvalue_int (p_rcollection,'NLMAX')

      SELECT CASE (rpreconditioner%ctypePreconditioning)
      CASE (CCPREC_NONE)
        ! No preconditioning. Do nothing.
      CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
        ! Preconditioning with a linear solver.
        !
        ! At first, assemble the preconditioner matrices on all levels
        ! and incorporate all boundary conditions.
        !
        ! Should we assemble the Newton part?
        
        bassembleNewton = .FALSE.
        
        IF (rpreconditioner%ctypePreconditioning .EQ. CCPREC_NEWTON) THEN
            
          ! Use Newton in any case.
          bassembleNewton = .TRUE.
          
        END IF
        
        ! Assemble the preconditioner matrices in rpreconditioner
        ! on all levels that the solver uses.
        CALL assembleLinsolMatrices (rpreconditioner,rmatrixComponents,&
            rcollection,bassembleNewton,rx)
          
        ! Our 'parent' (the caller of the nonlinear solver) has prepared
        ! a preconditioner node for us (a linear solver with symbolically
        ! factorised matrices). Get this from the collection.
      
        p_rsolverNode => rpreconditioner%p_rsolverNode

        ! Re-attach the system matrices to the solver.
        ! Note that no pointers and no handles are changed, so we can savely do
        ! that without calling linsol_doneStructure/linsol_doneStructure.
        ! This simply informs the solver about possible new scaling factors
        ! in the matrices in case they have changed...
        ALLOCATE(Rmatrices(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
        DO i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
          CALL lsysbl_duplicateMatrix ( &
            rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        END DO
        
        ! DEBUG!!!
        !CALL matio_writeBlockMatrixHR (Rmatrices(rnonlinearIteration%NLMAX), 'matrix',&
        !                               .TRUE., 0, 'matrixstat.txt','(E10.2)')
        
        CALL linsol_setMatrices(rpreconditioner%p_rsolverNode,Rmatrices(:))
        
        ! DEBUG!!!
        !DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        !  CALL storage_getbase_double (Rmatrices(i)% &
        !      RmatrixBlock(4,1)%h_Da,p_Ddata)
        !END DO
            
        ! The solver got the matrices; clean up Rmatrices, it was only of temporary
        ! nature...
        DO i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
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
        
        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
        ! for the next iteration.
        CALL linsol_doneData (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .EQ. 0
        
      END SELECT
      
      IF (bsuccess) THEN
        ! Filter the final defect
        p_RfilterChain => rpreconditioner%p_RfilterChain
        CALL filter_applyFilterChainVec (rd, p_RfilterChain)
      END IF
      
    CONTAINS
      
      SUBROUTINE assembleLinsolMatrices (rpreconditioner,rmatrixComponents,rcollection,&
          bassembleNewton,rx)

      USE linearsystemblock
      USE collection

      ! Assembles on every level a matrix for the linear-solver/Newton preconditioner.
      ! bnewton allows to specify whether the Newton matrix or only the standard
      ! system matrix is evaluated. The output is written to the p_rpreconditioner 
      ! matrices specified in the rnonlinearIteration structure.

      ! Spatial preconditioner structure that defines all parameters how to perform
      ! preconditioning.
      TYPE(t_ccspatialPreconditioner), INTENT(IN)    :: rpreconditioner

      ! Level independent configuration of the core equation
      TYPE(t_ccmatrixComponents), INTENT(IN)      :: rmatrixComponents

      ! Reference to a collection structure that contains all parameters of the
      ! discretisation (for nonlinearity, etc.).
      TYPE(t_collection), INTENT(INOUT)                :: rcollection

      ! TRUE  = Assemble the Newton preconditioner.
      ! FALSE = Assemble the standard defect correction preconditioner
      !         (i.e. the linearised system matrix).
      LOGICAL, INTENT(IN) :: bassembleNewton
      
      ! Current iteration vector. 
      TYPE(t_vectorBlock), INTENT(IN), TARGET          :: rx

      ! local variables
      REAL(DP) :: dnewton
      INTEGER :: ilev
      TYPE(t_matrixBlock), POINTER :: p_rmatrix,p_rmatrixFine
      TYPE(t_vectorScalar), POINTER :: p_rvectorTemp
      TYPE(t_vectorBlock), POINTER :: p_rvectorFine,p_rvectorCoarse
      TYPE(t_interlevelProjectionBlock), POINTER :: p_rprojection
      TYPE(t_ccmatrixComponents) :: rmatrixAssembly

      ! A filter chain for the linear solver
      TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
      
      ! DEBUG!!!
    !    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !    call lsysbl_getbase_double (rd,p_def)
    !    call lsysbl_getbase_double (rx,p_vec)
    !    NLMAX = collct_getvalue_int (p_rcollection,'NLMAX')

        ! Get the interlevel projection structure and the temporary vector
        ! from the collection.
        ! Our 'parent' prepared there how to interpolate the solution on the
        ! fine grid to coarser grids.
        p_rprojection => rpreconditioner%p_rprojection
        p_rvectorTemp => rpreconditioner%p_rtempVectorSc

        ! Get the filter chain. We need tghat later to filter the matrices.        
        p_RfilterChain => rpreconditioner%p_RfilterChain

        ! Initialise the matrix assembly structure rmatrixAssembly to describe the
        ! matrix we want to have. We have to initialise the adaptivity constants,
        ! which are not part of the standard initialisation.
        rmatrixAssembly = rmatrixComponents
        rmatrixAssembly%iadaptiveMatrices = &
            rpreconditioner%rfinalAssembly%iadaptiveMatrices
        rmatrixAssembly%dadmatthreshold = &
            rpreconditioner%rfinalAssembly%dadmatthreshold

        ! On all levels, we have to set up the nonlinear system matrix,
        ! so that the linear solver can be applied to it.
        
        NULLIFY(p_rmatrix)

        DO ilev=rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
        
          ! Get the matrix on the current level.
          ! Shift the previous matrix to the pointer of the fine grid matrix.
          p_rmatrixFine => p_rmatrix
          p_rmatrix => rpreconditioner%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
          ! On the highest level, we use rx as solution to build the nonlinear
          ! matrix. On lower levels, we have to create a solution
          ! on that level from a fine-grid solution before we can use
          ! it to build the matrix!
          IF (ilev .EQ. rpreconditioner%NLMAX) THEN
          
            p_rvectorCoarse => rx
            
          ELSE
            ! We have to discretise a level hierarchy and are on a level < NLMAX.
            
            ! Get the temporary vector on level i. Will receive the solution
            ! vector on that level. 
            p_rvectorCoarse => rpreconditioner%RcoreEquation(ilev)%p_rtempVector
            
            ! Get the solution vector on level i+1. This is either the temporary
            ! vector on that level, or the solution vector on the maximum level.
            IF (ilev .LT. rpreconditioner%NLMAX-1) THEN
              p_rvectorFine => rpreconditioner%RcoreEquation(ilev+1)%p_rtempVector
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

          ! Set the pointers in the rmatrixAssembly structure according
          ! to the current level.
          rmatrixAssembly%p_rdiscretisation         => &
              rpreconditioner%RcoreEquation(ilev)%p_rmatrix%p_rblockDiscretisation
          rmatrixAssembly%p_rmatrixStokes           => &
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixStokes          
          rmatrixAssembly%p_rmatrixB1             => &
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixB1              
          rmatrixAssembly%p_rmatrixB2             => &
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixB2              
          rmatrixAssembly%p_rmatrixMass           => &
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixMass            
          rmatrixAssembly%p_rmatrixIdentityPressure => &
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixIdentityPressure

          ! Assemble the matrix.
          ! If we are on a lower level, we can specify a 'fine-grid' matrix.
          IF (ilev .EQ. rpreconditioner%NLMAX) THEN
            CALL c2d2_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rmatrixAssembly,p_rvectorCoarse)
          ELSE
            CALL c2d2_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rmatrixAssembly,p_rvectorCoarse,p_rmatrixFine)
          END IF

          ! Boundary conditions
          ! ---------------------------------------------------

          IF (ASSOCIATED(p_RfilterChain)) THEN
            ! Apply the filter chain to the matrix.
            ! As the filter consists only of an implementation filter for
            ! boundary conditions, this implements the boundary conditions
            ! into the system matrix.
            CALL filter_applyFilterChainMat (p_rmatrix, p_RfilterChain)
          ELSE
            ! Call the matrix filter for the boundary conditions to include the BC's
            ! into the matrix.
            CALL matfil_discreteBC (p_rmatrix)
            CALL matfil_discreteFBC (p_rmatrix)
          END IF
            
          ! 'Nonlinear' boundary conditions like slip boundary conditions
          ! are not implemented with a filter chain into a matrix.
          ! Call the appropriate matrix filter of 'nonlinear' boundary
          ! conditions manually:
          CALL matfil_discreteNLSlipBC (p_rmatrix,.TRUE.)
            
          ! DEBUG!!!
          !CALL matio_writeBlockMatrixHR (p_rmatrix, 'matrix',&
          !                              .TRUE., 0, 'matrix.txt','(1X,E20.10)')

        END DO
        
      END SUBROUTINE
      
    END SUBROUTINE

END MODULE
