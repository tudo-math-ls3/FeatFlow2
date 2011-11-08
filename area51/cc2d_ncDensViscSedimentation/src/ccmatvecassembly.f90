!##############################################################################
!# ****************************************************************************
!# <name> cc2dmatrixassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the very basic matrix assembly routines for the
!# core equation. It's independent of any nonlinear iteration and provides
!# just one functionality: Assemble a matrix or a vector based on a given
!# set of parameters.
!#
!# The discretised core equation reads at the moment:
!#
!#  $$        A_1 y   +  \eta B p   = f_1 $$
!#  $$   \tau B^T y                 = f_2 $$
!#
!# with
!#
!#   $$ A_1 = \alpha M  +  \theta L  +  \gamma N(y) + \text{newton }N*(y)$$
!#
!# and
!#
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity $y\delta(\cdot)$ includung stabilisation,
!#   $N*(y)$ = Adjoint term $\cdot\delta(y)$ of the nonlinearity,
!#             used for the Newton matrix
!#
!#   $\alpha$ - weight in front of the mass matrix;
!#                =0 for stationary problem,
!#   $\theta$ - weight for the Laplace matrix,
!#   $\gamma$ - weight in front of the nonlinearity;
!#                =0 for Stokes system,
!#   $\eta$   - Switches the 'B'-term on/off,
!#   $\tau$   - Switches the 'B^T'-term on/off,
!#   newton   - Weight for the Newton term
!#
!# This equation can be written as a nonlinear system $A(y)(y,p) = (f1,f2)$
!# with a nonlinear matrix $A(\cdot)$. The structure t_nonlinearCCmatrix
!# contains a description of this matrix, With this description, it's possible
!# to do matrix vector multiplication or to 'evaluate' the matrix at a
!# given 'point' $y$ to get the 'linearised' matrix $A(y)$.
!#
!# The module contains the following routines:
!#
!# 1.) cc_assembleMatrix
!#     -> Assembles a matrix based on a set of input parameters, i.e.
!#        evaluates the system matrix $A(\cdot)$ at a point $y$ to create
!#        the linear matrix $A(y)$.
!#
!# 2.) cc_nonlinearMatMul
!#     -> Performs a matrix vector multiplication $d:=A(y)x+d$.
!#
!# </purpose>
!##############################################################################

module ccmatvecassembly

  use fsystem
  use storage
  use genoutput

  use linearsystemblock
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use scalarpde
  use derivatives
  use feevaluation
  use vectorio
  use cubature
  use basicgeometry
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use trilinearformevaluation
  use matrixio

  use statistics
  
  use convection
  use ccbasic
  use stdoperators
  use linearsystemscalar
  use geometry
  implicit none
  
!<constants>

!<constantblock description="Identifiers for the 'coperation' input parameter of the matrix assembly routine">

  ! Allocate memory if necessary.
  integer(I32), parameter :: CCMASM_ALLOCMEM              = 1
  
  ! Compute all matrix entries.
  integer(I32), parameter :: CCMASM_COMPUTE               = 2
  
  ! Allocate memory and compute matrix entries.
  integer(I32), parameter :: CCMASM_ALLOCANDCOMPUTE       = 3
  
  ! Bypass memory allocation for matrices.
  integer(I32), parameter :: CMASM_QUICKREFERENCES        = 4
  

  
!</constantblock>

!<constantblock description="Identifiers for the IUPWIND parameter that specifies how to set up the nonlinearity or stabilisation.">

  ! Streamline diffusion; configured by dupsam
  integer, parameter :: CCMASM_STAB_STREAMLINEDIFF    = 0

  ! 1st-order upwind; configured by dupsam
  integer, parameter :: CCMASM_STAB_UPWIND            = 1
  
  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'
  integer, parameter :: CCMASM_STAB_EDGEORIENTED      = 2

  ! Fast Edge-oriented stabilisation; configured by dupsam as 'gamma'. Preconputed matrix.
  integer, parameter :: CCMASM_STAB_FASTEDGEORIENTED  = 3

  ! Streamline diffusion; configured by dupsam. New implementation
  integer, parameter :: CCMASM_STAB_STREAMLINEDIFF2   = 4
  
  ! no stabilisation, just a bilinearform is evaluated
  integer, parameter :: CCMASM_NOSTAB_BILINEARFORM   = 5
!</constantblock>

!<constantblock description="Matrix type ID's specifying the general matrix class to set up.">

  ! Standard matrix.
  integer, parameter :: CCMASM_MTP_AUTOMATIC         = 0
  
  ! Standard matrix with decoupled velocity blocks
  integer, parameter :: CCMASM_MTP_DECOUPLED         = 1
  
  ! Extended 'full-tensor' matrix with submatrices A11, A12, A21, A22, all independent from
  ! each other.
  integer, parameter :: CCMASM_MTP_FULLTENSOR        = 2

!</constantblock>

  logical, parameter :: CCMASM_NONCNSTDENS           = .true.

!</constants>

!<types>

!<typeblock>

  ! This routine describes the nonlinear system matrix. The system matrix
  ! does actually not exist in memory -- since it's nonlinear! Therefore,
  ! this structure contains all parameters and settings which are necessary
  ! do apply(!) the matrix to a vector or to evaluate it.
  ! ('Evaluate a nonlinear matrix' means: Using a given FE-function $y$,
  ! assemble the linear matrix A(y).)
  type t_nonlinearCCMatrix
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    real(DP) :: dalpha = 0.0_DP
    
    ! THETA-parameter that controls the weight of the Stokes matrix
    ! in the core equation. =1.0 for stationary simulations.
    real(DP) :: dtheta = 0.0_DP
    
    ! GAMMA-parameter that controls the weight in front of the
    ! nonlinearity N(u). =1.0 for Navier-Stokes, =0.0 for Stokes equation.
    real(DP) :: dgamma = 0.0_DP

    ! ETA-parameter that switch the B-term on/off in the matrix.
    real(DP) :: deta = 0.0_DP
    
    ! TAU-parameter that switch the B^T-term on/off in the matrix.
    real(DP) :: dtau = 0.0_DP

    ! Weight for the Newton matrix N*(u).
    ! = 0.0 deactivates the Newton part.
    real(DP) :: dnewton = 0.0_DP

    ! STABILISATION: Parameter that defines how to set up the nonlinearity and
    ! whether to use some kind of stabilisation. One of the CCMASM_STAB_xxxx
    ! constants. Standard is CCMASM_STAB_STREAMLINEDIFF.
    integer :: iupwind = CCMASM_STAB_STREAMLINEDIFF
    
    ! STABILISATION: Viscosity parameter. Used for stabilisation schemes when
    ! a nonlinearity is set up.
    real(DP) :: dnu = 0.0_DP
    
    ! STABILISATION: Stabilisation parameter for streamline diffusion, upwind and
    ! edge-oriented stabilisation. If iupwind=CCMASM_STAB_STREAMLINEDIFF, a value of
    ! 0.0 deactivates any stabilisation.
    real(DP) :: dupsam = 0.0_DP
    
    ! STABILISATION: Specifies how the local H should be calculated for
    ! streamline diffusion.
    integer :: clocalH
    
    ! MATRIX RESTRICTION: Parameter to activate matrix restriction.
    ! Can be used to generate parts of the matrices on coarse grids where the
    ! aspect ratio of the cells is large. Only applicable for $\tilde Q_1$
    ! discretisations.
    ! Standard = 0 = deactivate matrix restriction
    integer :: iadaptiveMatrices = 0
    
    ! MATRIX RESTRICTION: Threshold parameter for adaptive matrix generation
    ! of coarse grid matrices (here: maximum aspect ratio).
    ! Only applicable if iadaptiveMatrices <> 0.
    ! Standard = 20.0
    real(DP) :: dadmatthreshold = 20.0_DP
    
    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...).
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! Pointer to a template FEM matrix that defines the structure of
    ! Laplace/Stokes/... matrices.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM => null()

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2) matrices.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateGradient => null()

    ! Pointer to Stokes matrix (=nu*Laplace).
    type(t_matrixScalar), pointer :: p_rmatrixStokes => null()

    ! Pointer to a B1-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB1 => null()

    ! Pointer to a B2-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB2 => null()

    ! Pointer to a B1^T-matrix.
    ! This pointer may point to NULL(). In this case, B1^T is created
    ! by 'virtually transposing' the B1 matrix.
    type(t_matrixScalar), pointer :: p_rmatrixD1 => null()

    ! Pointer to a B2-matrix.
    ! This pointer may point to NULL(). In this case, B2^T is created
    ! by 'virtually transposing' the B2 matrix.
    type(t_matrixScalar), pointer :: p_rmatrixD2 => null()

    ! Pointer to a Mass matrix.
    ! May point to NULL() during matrix creation.
    type(t_matrixScalar), pointer :: p_rmatrixMass => null()

    ! An object specifying the block discretisation
    ! for the (edge) stabilisation.
    type(t_blockDiscretisation), pointer :: p_rdiscretisationStabil => null()

    ! Pointer to the Jump stabilisation matrix.
    ! Only active if iupwind=CCMASM_STAB_FASTEDGEORIENTED, otherwise not associated
    type(t_matrixScalar), pointer :: p_rmatrixStabil => NULL()
    
    real(DP) :: drho1 = 1.0_DP
    
    real(DP) :: drho2 = 1.0_DP

  end type

!</typeblock>


!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleMatrix (coperation,cmatrixType,rmatrix,&
      rnonlinearCCMatrix,rproblem,rvector,rfineMatrix)

!<description>
  ! This routine assembles a global matrix, i.e. it evaluates the nonlinear
  ! system matrix $A(\cdot)$ at a point $y$. The rnonlinearCCMatrix
  ! must contain a description of the matrix.
  ! The 'coperation' parameter tells the routine what to do.
  ! The destination matrix rmatrix, which receives the evaluated matrix
  ! A(rvector), is then set up or updated.
  !
  ! The parameters rvector and rfineMatrix are optional. rvector must be
  ! specified, if the nonlinearity is activated (parameter $\gamma\not=0$ in
  ! rnonlinearCCMatrix). This vector specifies the 'solution' where the
  ! nonlinearity $u\nabla u$ is evaluated.
  ! rfineMatrix allows to specify a matrix of a 'one level refined mesh'. This
  ! is usually used when setting up preconditioners over multiple levels.
  ! Specifying such a matrix allows the routine (depending on the discretisation)
  ! to include some special stabilisation terms into the matrix rmatrix.
  !
  ! The routine does not include any boundary conditions in the matrix.
!</description>

!<input>

  ! One of the CCMASM_xxxx-constants. This parameter specifies 'what to do'.
  ! Using the CCMASM_ALLOCMEM constant, the routine will allocate memory in
  ! rmatrix for all the matrix blocks but will not compute any entries.
  ! Using the CCMASM_COMPUTE constant, the routine will compute the matrix
  ! entries, while assuming that the memory was allocated previously.
  ! Using the CCMASM_ALLOCANDCOMPUTE constant, the routine will do both,
  ! allocate memory and compute the entries.
  !
  ! If any of the CCMASM_ALLOCxxxx constants is specified here and rmatrix is
  ! already initialised, rmatrix is released and completely rebuild in memory!
  ! Therefore, CCMASM_COMPUTE should be used to update an existing matrix.
  !
  ! The constant CMASM_QUICKREFERENCES may be specified additional to one of
  ! the other constants (e.g. as 'CCMASM_ALLOCANDCOMPUTE+CMASM_QUICKREFERENCES').
  ! If specified, the routine tries to avoid memory allocation. This means e.g.
  ! that references to the original gradient (B-)matrices from rnonlinearCCMatrix
  ! are written to rmatrix; matrix entries are not copied!
  ! (This can be used e.g. for setting up a matrix for building a defect
  !  vector without copying matrix data.)
  ! In this case, the caller MUST NOT CHANGE rmatrix in any way, otherwise
  ! the original (template) matrices would be changed!
  integer(I32), intent(IN) :: coperation

  ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
  ! constants.
  ! Usually, CCMASM_MTP_AUTOMATIC is used here. This will automatically determine
  ! the 'smallest possible' matrix structure fitting the needs of the input
  ! parameters. By specifying another matrix type, the caller can explicitly
  ! take influence on the general matrix structure.
  !
  ! If the matrix already exists and only its entries are to be computed,
  ! CCMASM_MTP_AUTOMATIC should be specified here.
  integer, intent(IN) :: cmatrixType

  ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, p_rmatrixTemplateXXXX
  ! must be initialised as well as p_rdiscretisation!
  ! The new matrix is created based p_rmatrixTemplateXXXX as well as
  ! p_rdiscretisation. Memory is allocated automatically if it's missing.
  type(t_nonlinearCCMatrix), intent(IN) :: rnonlinearCCMatrix

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity.
  type(t_vectorBlock), intent(IN), optional :: rvector

  ! OPTIONAL: This parameter allows to specify a 'fine grid matrix'. This is
  ! usually done when assembling matrices on multiple levels. If specified, the
  ! routine will (if possible) try to include a level-dependent stabilisation
  ! term into the matrix (-> e.g. constant matrix restriction for nonparametric
  ! Rannacher-Turek element if cells are too anisotropic).
  type(t_matrixBlock), intent(IN), optional :: rfineMatrix
  
!</input>

!<inputoutput>

  ! The destination matrix which should be set up.
  ! If not initialised, a new matrix is created (as if CCMASM_ALLOCxxxx
  ! was specified).
  ! If initialised, the existing matrix is updated or recreated, depending on
  ! coperation.
  type(t_matrixBlock), intent(INOUT) :: rmatrix
  
  type(t_problem), intent(inout) :: rproblem
  
!</inputoutput>
  
!</subroutine>

    ! local variables
    logical :: ballocate
    
    ballocate = .false.
    if ((rmatrix%NEQ .le. 0) .or. &
        iand(coperation,CCMASM_ALLOCMEM) .ne. 0) then
      ballocate = .true.
    end if
    
    ! What should we do? Allocate memory?
    if (ballocate) then
    
      ! Release the matrix if present.
      call lsysbl_releaseMatrix (rmatrix)
    
      ! Create a complete new matrix.
      call allocMatrix (cmatrixType,rnonlinearCCMatrix,rmatrix)
    end if
   
    if (iand(coperation,CCMASM_COMPUTE) .ne. 0) then

      ! The system matrix looks like:
      !
      !    ( A11  A12  B1  )
      !    ( A21  A22  B2  )
      !    ( B1^T B2^T     )
      !
      ! Assemble the velocity submatrices
      !
      !    ( A11  A12   .  )
      !    ( A21  A22   .  )
      !    (  .    .    .  )
      
      call assembleVelocityBlocks (&
          rnonlinearCCMatrix,rmatrix,rvector,1.0_DP,rproblem)
      
      ! Assemble the gradient submatrices
      !
      !    ( .    .    B1  )
      !    ( .    .    B2  )
      !    ( B1^T B2^T .   )
      
      call assembleGradientMatrices (rnonlinearCCMatrix,rmatrix,&
        iand(coperation,CMASM_QUICKREFERENCES) .ne. 0)

      ! 2.) Initialise the weights for the B-matrices
      !
      !    (  .    .   B1  )
      !    (  .    .   B2  )
      !    ( B1^T B2^T  .  )
      
      rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearCCMatrix%deta
      rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearCCMatrix%deta
      
      rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearCCMatrix%dtau
      rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearCCMatrix%dtau
      
      ! Matrix restriction
      ! ---------------------------------------------------
      !
      ! For the construction of matrices on lower levels, call the matrix
      ! restriction. In case we have a uniform discretisation with Q1~,
      ! iadaptivematrix is <> 0 and so this will rebuild some matrix entries
      ! by a Galerkin approach using constant prolongation/restriction.
      ! This helps to stabilise the solver if there are elements in the
      ! mesh with high aspect ratio.
      if (present(rfineMatrix)) then
        call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(1,1), &
            rmatrix%RmatrixBlock(1,1), &
            rnonlinearCCMatrix%iadaptiveMatrices, rnonlinearCCMatrix%dadmatthreshold)
            
        if (.not. lsyssc_isMatrixContentShared(rmatrix%RmatrixBlock(2,2))) then
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(2,2), &
              rmatrix%RmatrixBlock(2,2), &
              rnonlinearCCMatrix%iadaptiveMatrices, rnonlinearCCMatrix%dadmatthreshold)
        end if
      end if
      
    end if
    
  contains
  
    ! -----------------------------------------------------
  
    subroutine allocMatrix (cmatrixType,rnonlinearCCMatrix,rmatrix)
    
    ! Allocates memory for the system matrix. rnonlinearCCMatrix provides information
    ! about the submatrices that are 'plugged into' rmatrix.
    ! Therefore, before this routine is called, rnonlinearCCMatrix must have been set up.

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    integer, intent(IN) :: cmatrixType

    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(IN), target :: rnonlinearCCMatrix

    ! A block matrix that receives the basic system matrix.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
      ! local variables
      logical :: bdecoupled, bfulltensor

      ! A pointer to the system matrix and the RHS/solution vectors.
      type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM

      ! A pointer to the discretisation structure with the data.
      type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
      ! Determine the shape of the matrix
      bdecoupled = cmatrixType .eq. CCMASM_MTP_DECOUPLED
      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
      
      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
        bfulltensor = rnonlinearCCMatrix%dnewton .ne. 0.0_DP
      end if
    
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rnonlinearCCMatrix%p_rdiscretisation
      
      ! Get a pointer to the template FEM matrix. If that doesn't exist,
      ! take the Stokes matrix as template.
      p_rmatrixTemplateFEM => rnonlinearCCMatrix%p_rmatrixTemplateFEM
      if (.not. associated(p_rmatrixTemplateFEM)) &
        p_rmatrixTemplateFEM => rnonlinearCCMatrix%p_rmatrixStokes
      if (.not. associated(p_rmatrixTemplateFEM)) then
        call output_line ('Cannot set up A matrices in system matrix!', &
            OU_CLASS_ERROR,OU_MODE_STD,'allocMatrix')
        call sys_halt()
      end if

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      if (associated(p_rdiscretisation)) then
        call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
      else
        ! No discretisation structure; create the matrix directly as 3x3 matrix.
        call lsysbl_createEmptyMatrix (rmatrix,NDIM2D+1)
      end if
        
      ! Let's consider the global system in detail. The standard matrix It has
      ! roughly the following shape:
      !
      !    ( A11       B1  ) = ( A11  A12  A13 )
      !    (      A22  B2  )   ( A21  A22  A23 )
      !    ( B1^T B2^T .   )   ( A31  A32  A33 )
      !
      ! All matrices may have multiplication factors in their front.
      !
      ! The structure of the matrices A11 and A22 of the global system matrix
      ! is governed by the template FEM matrix.
      ! Initialise them with the same structure, i.e. A11, A22 share (!) their
      ! structure (not the entries) with that of the template matrix.
      !
      ! For this purpose, use the "duplicate matrix" routine.
      ! The structure of the matrix is shared with the template FEM matrix.
      ! For the content, a new empty array is allocated which will later receive
      ! the entries.
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      if (.not. bdecoupled .and. .not. bfulltensor) then
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content.
        call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
      else
      
        ! Otherwise, create another copy of the template matrix.
        call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                    
      end if
      
      ! Manually change the discretisation structure of the Y-velocity
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))

      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      if (bfulltensor) then

        ! We have a matrix in the following shape:
        !
        !    ( A11  A12  B1  )
        !    ( A21  A22  B2  )
        !    ( B1^T B2^T .   )
        !
        ! Create A12 and A21.
      
        if (rmatrix%RmatrixBlock(1,2)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(1,2), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          ! Create a new matrix A21 in memory. create a new matrix
          ! using the template FEM matrix...
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(2,1), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
            
        end if
        
      end if

      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create empty space for the entries.
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,3),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,3),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
      ! Now, prepare D1 and D2.
      
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixD1, &
                                    rmatrix%RmatrixBlock(3,1),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixD2, &
                                    rmatrix%RmatrixBlock(3,2),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      ! That's it, all submatrices are basically set up.
      !
      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      call lsysbl_updateMatStrucInfo (rmatrix)
        
    end subroutine
    
    ! -----------------------------------------------------

    subroutine assembleVelocityBlocks (rnonlinearCCMatrix,rmatrix,rvector,dvectorWeight,rproblem)
        
    ! Assembles the velocity matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    
    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(IN) :: rnonlinearCCMatrix
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    type(t_problem), intent(inout) :: rproblem
    
    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    type(t_vectorBlock), optional :: rvector
    
    ! Weight for the velocity vector; standard = -1.
    real(DP), intent(IN), optional :: dvectorWeight
    
    ! local variables
    logical :: bshared
    integer :: iupwind
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_jumpStabilisation) :: rjumpStabil
    type(t_matrixScalar) :: rMass1
    type(t_matrixBlock) :: rmatrix1,rmatrix2
    real(DP) :: dvecWeight
    type(t_bilinearform) :: rform
    
    rproblem%rcollection%Dquickaccess(5) = rproblem%drho1
    rproblem%rcollection%Dquickaccess(6) = rproblem%drho2
    
    
      ! Standard value for dvectorWeight is = -1.
      dvecWeight = -1.0_DP
      if (present(dvectorWeight)) dvecWeight = dvectorWeight
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      ! Allocate memory if necessary. Normally this should not be necessary...
      ! A11:
      if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) then
        call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
      end if
    
      ! A22:
      if (.not. bshared) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
        end if
      end if

      ! A12/ A21:
      if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,2))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
        end if
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
        end if
      end if
    
      ! ---------------------------------------------------
      ! Plug in the mass matrix?
      if (rnonlinearCCMatrix%dalpha .ne. 0.0_DP) then
       
        if(CCMASM_NONCNSTDENS)then
        
        ! Allocate memory if necessary. Normally this should not be necessary...
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
        end if
        
        
        ! now set up a bilinearform
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        !    ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .false.
        rform%BconstantCoeff = .false.
        rform%Dcoefficients(1)  = rnonlinearCCMatrix%dalpha

        ! we here want to create a matrix for the term (rho/dT)*u_n
        ! and to incorporate this into our matrix rmatrix
        call bilf_createMatrixStructure(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial,&
                                        rmatrix%RmatrixBlock(1,1)%cmatrixFormat,&
                                        rMass1,&
                                        rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest,BILF_MATC_EDGEBASED)


        ! Now we can build the matrix entries.
        ! We specify the callback function coeff_Stokes for the coefficients.
        ! As long as we use constant coefficients, this routine is not used.
        ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
        ! the framework will call the callback routine to get analytical data.
        !
        ! We pass our collection structure as well to this routine,
        ! so the callback routine has access to everything what is
        ! in the collection.
        call bilf_buildMatrixScalar (rform,.true.,&
                                     rMass1,cc_Rho,&
                                     rproblem%rcollection)
                                     
        ! lsyssc_matrixLinearComb (rmatrixA,cA,rmatrixB,cB,rmatrixC,&
        !    bmemory,bsymb,bnumb,bisExactStructure)

        ! Adds constant times a matrix to another matrix
        !   rmatrixC = cA*rmatrixA + cB*rmatrixB
        
!        call lsysbl_duplicateMatrix(rmatrix,rmatrix2,LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
!        call matio_writeMatrixHR (rmatrix2%RmatrixBlock(1,1), 'rmatrix2',&
!                                  .TRUE., 0, 'rmatrix2.txt', '(E20.5)')
!
!        call matio_writeMatrixHR (rMass1, 'rMass1',&
!                                  .TRUE., 0, 'rMass1.txt', '(E20.5)')
        
        !call lsysbl_clearMatrix(rmatrix2)
        ! hier selber die Masse Matrix aufbauen mit (phi_i,phi_j) und
        ! rho als nicht konstantem Koeffizienten
        call lsyssc_matrixLinearComb (&
            rMass1,rnonlinearCCMatrix%dalpha,&
            rmatrix%RmatrixBlock(1,1),0.0_DP,&
            rmatrix%RmatrixBlock(1,1),&
            .false.,.false.,.true.,.true.)
            
!        call matio_writeMatrixHR (rform%RmatrixBlock(1,1), 'massbilf',&
!                                  .TRUE., 0, 'massbilf.txt', '(E20.5)')
            
            
        if (.not. bshared) then
          ! Allocate memory if necessary. Normally this should not be necessary...
          if (.not. lsyssc_hasMatrixContent (rmatrix2%RmatrixBlock(2,2))) then
            call lsyssc_allocEmptyMatrix (rmatrix2%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
          end if
          
          call lsyssc_matrixLinearComb (&
              rMass1,rnonlinearCCMatrix%dalpha,&
              rmatrix2%RmatrixBlock(2,2),0.0_DP,&
              rmatrix2%RmatrixBlock(2,2),&
              .false.,.false.,.true.,.true.)
        end if
            
        ! release the temporary matrix
        call lsyssc_releaseMatrix(rMass1)
        !call lsysbl_releaseMatrix(rmatrix2)
   
        else
        ! Allocate memory if necessary. Normally this should not be necessary...
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
        end if
      
!        call matio_writeMatrixHR (rmatrix%RmatrixBlock(1,1), 'rmatrix',&
!                                  .TRUE., 0, 'rmatrix.txt', '(E20.5)')
!
!        call matio_writeMatrixHR (rnonlinearCCMatrix%p_rmatrixMass, 'p_rmatrixMass',&
!                                  .TRUE., 0, 'p_rmatrixMass.txt', '(E20.5)')
        
        
      
        call lsyssc_matrixLinearComb (&
            rnonlinearCCMatrix%p_rmatrixMass       ,rnonlinearCCMatrix%dalpha,&
            rmatrix%RmatrixBlock(1,1),0.0_DP,&
            rmatrix%RmatrixBlock(1,1),&
            .false.,.false.,.true.,.true.)

!        call matio_writeMatrixHR (rmatrix%RmatrixBlock(1,1), 'mass',&
!                                  .TRUE., 0, 'mass.txt', '(E20.5)')

            
        if (.not. bshared) then
          ! Allocate memory if necessary. Normally this should not be necessary...
          if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) then
            call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
          end if

          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rmatrixMass     ,rnonlinearCCMatrix%dalpha,&
              rmatrix%RmatrixBlock(2,2),0.0_DP,&
              rmatrix%RmatrixBlock(2,2),&
              .false.,.false.,.true.,.true.)
        end if
        
        end if
        
      else
      
        ! Otherwise, initialise the basic matrix with 0.0
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
        
        if (.not. bshared) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
        end if
        
      end if
      
      ! If the submatrices A12 and A21 exist, fill them with zero.
      ! If they don't exist, we don't have to do anything.
      if (rnonlinearCCMatrix%dnewton .ne. 0.0_DP) then
        rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
      else
        rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
      end if
      
      if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
      end if
        
      ! ---------------------------------------------------
      ! Plug in the Stokes matrix?
      if (rnonlinearCCMatrix%dtheta .ne. 0.0_DP) then
        call lsyssc_matrixLinearComb (&
            rnonlinearCCMatrix%p_rmatrixStokes     ,rnonlinearCCMatrix%dtheta,&
            rmatrix%RmatrixBlock(1,1),1.0_DP,&
            rmatrix%RmatrixBlock(1,1),&
            .false.,.false.,.true.,.true.)
            
        if (.not. bshared) then
          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rmatrixStokes   ,rnonlinearCCMatrix%dtheta,&
              rmatrix%RmatrixBlock(2,2),1.0_DP,&
              rmatrix%RmatrixBlock(2,2),&
              .false.,.false.,.true.,.true.)
        end if
      end if
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if (rnonlinearCCMatrix%dgamma .ne. 0.0_DP) then
      
        if (.not. present(rvector)) then
          call output_line ('Velocity vector not present!', &
                             OU_CLASS_ERROR,OU_MODE_STD,'cc_assembleMatrix')
          stop
        end if
      
        select case (rnonlinearCCMatrix%iupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Streamline diffusion.

          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)

          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rnonlinearCCMatrix%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion2%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion2%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion2%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)


        case (CCMASM_STAB_UPWIND)
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = rnonlinearCCMatrix%dupsam

          ! Matrix weight
          rupwind%dtheta = rnonlinearCCMatrix%dgamma
          
          if (rnonlinearCCMatrix%dnewton .ne. 0.0_DP) then
            call output_line ('Warning: Upwind does not support assembly '&
                //'of the Newton matrix!',OU_CLASS_TRACE1)
          end if
          
          ! Call the upwind method to calculate the nonlinear matrix.
          call conv_upwind2d (rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rupwind, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))
                              
          if (.not. bshared) then
            ! Modify also the matrix block (2,2)
            call conv_upwind2d (rvector, rvector, &
                                dvecWeight, 0.0_DP,&
                                rupwind, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))
          end if

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.

          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          call cc_nostabBilin (rvector, rvector, &
                               rstreamlineDiffusion, CONV_MODMATRIX, &
                               rmatrix,rproblem)
          
          ! DEBUG
!          call matio_writeMatrixHR (rmatrix%RmatrixBlock(1,1), 'nonlin',&
!                                    .TRUE., 0, 'nonlin.txt', '(E20.5)')

          
          
          ! Call the SD method to calculate the nonlinearity.
!          call conv_streamlineDiffusionBlk2d (&
!                              rvector, rvector, &
!                              dvecWeight, 0.0_DP,&
!                              rstreamlineDiffusion, CONV_MODMATRIX, &
!                              rmatrix)
!
!          ! DEBUG
!          call matio_writeMatrixHR (rmatrix%RmatrixBlock(1,1), 'nonlinSD',&
!                                    .TRUE., 0, 'nonlinSD.txt', '(E20.5)')
!
!
!          ! DEBUG
!          call matio_writeMatrixHR (rmatrix%RmatrixBlock(2,2), 'nonlinSD22',&
!                                    .TRUE., 0, 'nonlinSD22.txt', '(E20.5)')

!          stop
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rnonlinearCCMatrix%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
              rdiscretisation=rnonlinearCCMatrix%p_rdiscretisationStabil%RspatialDiscr(1))

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
              rdiscretisation=rnonlinearCCMatrix%p_rdiscretisationStabil%RspatialDiscr(1))
          end if

        case (CCMASM_STAB_FASTEDGEORIENTED)
          ! Jump stabilisation with precomputed matrix.

          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%dtheta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)
        
          ! Sum up the precomputed edge stabilisation matrix.
          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rmatrixStabil     ,rnonlinearCCMatrix%dtheta,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
          
          if (.not. bshared) then
            call lsyssc_matrixLinearComb (&
                rnonlinearCCMatrix%p_rmatrixStabil   ,rnonlinearCCMatrix%dtheta,&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if
        
        ! no stabilisation, bilinearform evaluation
        case(CCMASM_NOSTAB_BILINEARFORM)
          ! copy the streamline diffusion
          ! Streamline diffusion.

          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the nonlinearity.
          ! here we use the CONV_MODMATRIX parameter to
          ! control the behaviour of the assembly function
          ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
          !                 nonlinear part is added to the matrix.
          ! rstreamlineDiffusion: configuration structure
          ! Configuration block for the streamline diffusion scheme
          
          ! rvector is a block vector with the
          ! X- and Y-velocity.
          ! p_rvelX1 => rvector%RvectorBlock(1)
          ! p_rvelY1 => rvector%RvectorBlock(2)
          call cc_nostabBilin (rvector, rvector, &
                               rstreamlineDiffusion, CONV_MODMATRIX, &
                               rmatrix,rproblem)
          !DEBUG
!          call matio_writeMatrixHR (rmatrix%RmatrixBlock(1,1), 'matrix',&
!                                    .TRUE., 0, 'matrixNOSTAB.txt', '(E20.5)')
        
          
        case default
          call output_line ('Don''t know how to set up nonlinearity!?!', &
              OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityBlocks')
          call sys_halt()
        
        
        
        end select
        
        

      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
      
        select case (rnonlinearCCMatrix%iupwind)
        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rnonlinearCCMatrix%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(1,1),&
              rdiscretisation=rnonlinearCCMatrix%p_rdiscretisationStabil%RspatialDiscr(1))

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                rjumpStabil,CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
                rdiscretisation=rnonlinearCCMatrix%p_rdiscretisationStabil%RspatialDiscr(1))
          end if

        case (CCMASM_STAB_FastEDGEORIENTED)
          ! Fast Jump stabilisation. Precomputed matrix.
          !
          ! Sum
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rnonlinearCCMatrix%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Sum up the precomputed edge stabilisation matrix.
          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rmatrixStabil     ,rnonlinearCCMatrix%dtheta,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
          
          if (.not. bshared) then
            call lsyssc_matrixLinearComb (&
                rnonlinearCCMatrix%p_rmatrixStabil   ,rnonlinearCCMatrix%dtheta,&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if

        case default
          ! No stabilisation
        
        end select
      
      end if ! gamma <> 0
      
    end subroutine
      
    ! -----------------------------------------------------
    
    subroutine assembleGradientMatrices (rnonlinearCCMatrix,rmatrix,bsharedMatrix)
    
    ! Initialises the gradient/divergence matrices with entries from
    ! the rnonlinearCCMatrix structure.
    !
    ! The routine copies references from the submatrices tormatrix,
    ! but it does not initialise any matrix weights / scaling factors.
    !
    ! If bsharedMatrix=TRUE, the matrix is created using references to the
    ! matrix building blocks in rlevelInfo, thus sharing all information
    ! with those matrices in rnonlinearCCMatrix. In this case, the caller must
    ! not change the matrix entries, because this would change the
    ! original 'template' matrices!
    ! (This can be used e.g. for setting up a matrix for building a defect
    !  vector without copying matrix data.)
    ! If bsharedMatrix=TRUE on the other hand, the matrix entries of the
    ! original template (B-) matrices are copied in memory,
    ! so the new matrix is allowed to be changed!

    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(IN) :: rnonlinearCCMatrix

    ! Block matrix where the B-matrices should be set up
    type(t_matrixBlock), intent(INOUT) :: rmatrix

    ! Whether or not the matrix entries of the source gradient-matrices
    ! should be copied in memory.
    !
    ! If set to FALSE, the routine tries to initialise rmatrix
    ! only with references to the original matrices, thus the caller must not
    ! change the entries. Nevertheless, if rmatrix is the owner of one of the
    ! submatrices, the routine will always copy the matrix entries,
    ! as otherwise memory would have to be deallocated!
    !
    ! If set to TRUE, the entries of the source matrices in rnonlinearCCMatrix are
    ! copied, so the caller can change rmatrix afterwards (e.g. to implement
    ! boundary conditions).
    logical, intent(IN) :: bsharedMatrix

      ! local variables
      integer :: idubStructure,idubContent
      
      ! Initialise a copy flag that tells the duplicateMatrix-routine whether to
      ! copy the entries or to create references.
      if (bsharedMatrix) then
      
        idubContent = LSYSSC_DUP_SHARE
        
        ! Normally we share entries -- except for if the submatrices belong to
        ! rmatrix! To avoid memory deallocation in this case, we copy
        ! the entries.
        if ((.not. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(1,3))) .or.&
            (.not. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(2,3)))) then
          idubContent = LSYSSC_DUP_COPY
        end if
        
      else
      
        idubContent = LSYSSC_DUP_COPY
        
      end if
      
      idubStructure = LSYSSC_DUP_SHARE
      
      ! Let's consider the global system in detail:
      !
      !    ( A11  A12  B1  ) = ( A11  A12  A13 )
      !    ( A21  A22  B2  )   ( A21  A22  A23 )
      !    ( B1^T B2^T 0   )   ( A31  A32  A33 )
      !
      ! We exclude the velocity submatrices here, so our system looks like:
      !
      !    (           B1 ) = (           A13 )
      !    (           B2 )   (           A23 )
      !    ( B1^T B2^T    )   ( A31  A32  A33 )

      ! The B1/B2 matrices exist up to now only in rnonlinearCCMatrix.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The B-blocks
      ! are already prepared and memory for the entries is already allocated;
      ! so we only have to copy the entries.
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,3),&
                                    idubStructure,idubContent)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,3),&
                                    idubStructure,idubContent)
      
      ! Now, prepare B1^T and B2^T. These matrices always share
      ! their data with the 'template' matrices as the data in these
      ! matrices is usually not overwritten by boundary conditions...
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixD1, &
                                    rmatrix%RmatrixBlock(3,1),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixD2, &
                                    rmatrix%RmatrixBlock(3,2),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                                    
    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_nonlinearMatMul (rnonlinearCCMatrix,rx,rd,dcx,dcd,rproblem,ry)

!<description>
  ! This routine performs a matrix-vector multiplication with a nonlinear
  ! matrix:
  !      rd := cx A(ry) rx + cd rd
  ! with the system matrix A(.) defined by the configuration in rnonlinearCCMatrix.
  ! The caller must initialise the rnonlinearCCMatrix according to how the
  ! matrix should look like.
  !
  ! The parameter ry is optional. If specified, this parameter defines where to
  ! evaluate the nonlinearity (if the system matrix $A$ contains a nonlinearity).
  ! If not specified, ry=rx is assumed.
  !
  ! The routine will not include any boundary conditions in the defect.
!</description>

  ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! The caller must provide either p_rmatrixTemplateXXXX in this structure
  ! or set the p_rmatrixTemplateXXXX as well as p_rdiscretisation to
  ! appropriate values. This is necessary for exploiting then structure
  ! of the matrix.
  type(t_nonlinearCCMatrix), intent(IN) :: rnonlinearCCMatrix

  ! This vector specifies the 'x' that is multiplied to the matrix.
  type(t_vectorBlock), intent(IN), target :: rx

  ! Multiplication factor in front of the term 'A(ry) rx'.
  real(DP), intent(IN) :: dcx

  ! Multiplication factor in front of the term 'rd'.
  real(DP), intent(IN) :: dcd

  ! OPTIONAL: Point where to evaluate the nonlinearity. If not specified,
  ! ry=rx is assumed.
  type(t_vectorBlock), intent(IN), target, optional :: ry

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  type(t_vectorBlock), intent(INOUT) :: rd
  
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_ry
    type(t_matrixBlock) :: rmatrix
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_DdataX,p_DdataD
    
    rproblem%rcollection%Dquickaccess(5) = rproblem%drho1
    rproblem%rcollection%Dquickaccess(6) = rproblem%drho2
    
    call lsysbl_getbase_double (rx,p_DdataX)
    call lsysbl_getbase_double (rd,p_DdataD)
    
    p_ry => rx
    if (present(ry)) p_ry => ry

    ! Probably weight the input vector.
    if (dcd .ne. 1.0_DP) then
      call lsysbl_scaleVector (rd,dcd)
    end if
    
    ! The system matrix looks like:
    !
    !    ( A11  A12  B1  )
    !    ( A21  A22  B2  )
    !    ( B1^T B2^T .   )
    !
    ! Create a temporary matrix that covers this structure.
    call lsysbl_createMatBlockByDiscr (rnonlinearCCMatrix%p_rdiscretisation,rmatrix)
    
    ! Put references to the Stokes- and B-matrices to Aij. assembleVelocityDefect
    ! needs this template matrix to provide the structure for the stabilisation
    ! routines! The B-matrices are needed later.
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
    if (rnonlinearCCMatrix%dnewton .ne. 0.0_DP) then
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    end if
    
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1,&
        rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2,&
        rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixD1,&
        rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixD2,&
        rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    ! In the first step, we assemble the defect that arises in the velocity
    ! components. This is characterised by the following submatrix:
    !
    !    ( A11  A12  .  )
    !    ( A21  A22  .  )
    !    ( .    .    .  )
    !
    ! assembleVelocityDefect handles exactly these submatrices.

    call assembleVelocityDefect (rnonlinearCCMatrix,rmatrix,rx,rd,p_ry,-dcx,rproblem)
    
    ! Now, we treat all the remaining blocks. Let's see what is missing:
    !
    !    ( .    .    B1  )
    !    ( .    .    B2  )
    !    ( B1^T B2^T .   )

    ! To build the appropriate defect, we first remove the velocity blocks:
    
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,2))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,2))

    ! Initialise the weights for the B/B^T matrices
    rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearCCMatrix%deta
    rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearCCMatrix%deta
    
    rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearCCMatrix%dtau
    rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearCCMatrix%dtau

    ! ------------------------------------------------
    ! Build the defect by matrix-vector multiplication
    !
    ! Note that no time step or whatever is included here; everything
    ! is initialised with the multiplication factors in the submatrices
    ! from above!
    call lsysbl_blockMatVec (rmatrix, rx, rd, dcx, 1.0_DP)
    
    ! Release the temporary matrix, we don't need it anymore.
    call lsysbl_releaseMatrix (rmatrix)

  contains

    subroutine assembleVelocityDefect (rnonlinearCCMatrix,&
        rmatrix,rvector,rdefect,rvelocityVector,dvectorWeight,rproblem)
        
    ! Assembles the velocity defect in the block matrix rmatrix at position
    ! itop..itop+1 in the velocity vector. rdefect must have been initialised
    ! with the right hand side vector.
    !
    ! With a matrix 'A' of the theoretical form
    !
    !       A := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    !
    ! and c=dvectorWeight, the routine will construct
    !
    !       rdefect = rdefect - c A rvector
    
    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(IN) :: rnonlinearCCMatrix

    ! Reference to the system matrix. Only the structure of the matrix
    ! is used to reconstruct the structure of the discretisation.
    ! The content of the matrix is not changed or used.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    type(t_problem), intent(INOUT) :: rproblem
    
    ! Solution vector.
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Weight for the velocity vector; usually = 1.0
    real(DP), intent(IN) :: dvectorWeight
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. The first two blocks in that block vector are
    ! used as velocity field.
    type(t_vectorBlock), intent(IN) :: rvelocityVector

    ! local variables
    logical :: bshared
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(T_jumpStabilisation) :: rjumpStabil
    type(t_matrixScalar) :: rMass1
    type(t_bilinearForm) :: rform
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_DdataX,p_DdataD
    
    call lsysbl_getbase_double (rvector,p_DdataX)
    call lsysbl_getbase_double (rdefect,p_DdataD)

      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2)) .or.&
                (.not. lsyssc_hasMatrixContent(rmatrix%RmatrixBlock(1,1)) .and.&
                 .not. lsyssc_hasMatrixContent(rmatrix%RmatrixBlock(2,2)))

      ! ---------------------------------------------------
      ! Subtract the mass matrix stuff?
      if (rnonlinearCCMatrix%dalpha .ne. 0.0_DP) then
      if(CCMASM_NONCNSTDENS)then
      
        ! now set up a bilinearform
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        !    ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .false.
        rform%BconstantCoeff = .false.
        rform%Dcoefficients(1)  = rnonlinearCCMatrix%dalpha

        call bilf_createMatrixStructure(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial,&
                                        rmatrix%RmatrixBlock(1,1)%cmatrixFormat,&
                                        rMass1,&
                                        rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest,BILF_MATC_EDGEBASED)

        call bilf_buildMatrixScalar (rform,.true.,&
                                     rMass1,cc_Rho,&
                                     rproblem%rcollection)

        ! DEBUG
!        call matio_writeMatrixHR (rMass1, 'Mass',&
!                                  .TRUE., 0, 'myMass.txt', '(E20.5)')
                                     
        ! Incorporate the mass matrix into the rhs vector
        ! subroutine lsyssc_scalarMatVec (rmatrix, rx, ry, cx, cy, btranspose)
        ! Performs a matrix vector multiplicationwith a given scalar matrix:
        !    $$ Dy   =   cx * rMatrix * rx   +   cy * ry $$
        ! multiply block11 by u1
        call lsyssc_scalarMatVec (rMass1,&
            rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)
            
        ! multiply block22 by u2
        call lsyssc_scalarMatVec (rMass1,&
            rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)
            
        ! release the temporary matrix
        call lsyssc_releaseMatrix(rMass1)
      
      else
        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixMass, &
            rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixMass, &
            rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)
      end if
      end if
      
      ! ---------------------------------------------------
      ! Subtract the Stokes matrix stuff?
      if (rnonlinearCCMatrix%dtheta .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixStokes, &
            rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rnonlinearCCMatrix%dtheta, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixStokes, &
            rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
            -rnonlinearCCMatrix%dtheta, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if (rnonlinearCCMatrix%dgamma .ne. 0.0_DP) then
      
        ! Type of stablilisation?
        select case (rnonlinearCCMatrix%iupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the defect of the nonlinearity.
          ! As rrhsTemp shares its entries with rdefect, the result is
          ! directly written to rdefect!
          ! As velocity field, we specify rvelocityVector here. The first two
          ! subvectors are used as velocity field.
          
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityVector, rvelocityVector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)
                  
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rnonlinearCCMatrix%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion2%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion2%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion2%dnewton = rnonlinearCCMatrix%dnewton
                              
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rvector,rdefect,rvelocityVector)
                              
        case (CCMASM_STAB_UPWIND)
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = rnonlinearCCMatrix%dupsam

          ! Matrix weight
          rupwind%dtheta = rnonlinearCCMatrix%dgamma
          
          ! Call the upwind method to calculate the nonlinear defect.
          call conv_upwind2d (rvector, rvector, &
                              dvectorWeight, 0.0_DP,&
                              rupwind, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rvector,rdefect)
                              
          if (.not. bshared) then
            call output_line ('Upwind does not support independent A11/A22!', &
                OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
            call sys_halt()
          end if

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          
          if (rnonlinearCCMatrix%dnewton .eq. 0.0_DP) then

            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there's content at all),
            ! these matrices are ignored then by the kernel.
            
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
            
          else

            ! Clear A12/A21 that receives parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they aren't.
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          end if
         
          call cc_nostabBilin (rvector, rvector, &
                               rstreamlineDiffusion, CONV_MODDEFECT, &
                               rmatrix,rproblem,rsolution=rvector,rdefect=rdefect)
         
!          call vecio_writeBlockVectorHR (rdefect, 'def', .false.,&
!                                    0, 'defnostab.txt', '(E20.5)')
!
!         stop
          ! Call the SD method to calculate the nonlinearity.
!          call conv_streamlineDiffusionBlk2d (&
!                              rvector, rvector, &
!                              dvectorWeight, 0.0_DP,&
!                              rstreamlineDiffusion, CONV_MODDEFECT, &
!                              rmatrix,rsolution=rvector,rdefect=rdefect)
!
!          call vecio_writeBlockVectorHR (rdefect, 'defSD', .false.,&
!                                    0, 'defnSD.txt', '(E20.5)')
!          stop
!
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...<
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rnonlinearCCMatrix%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
              rsolution=rvector,rdefect=rdefect,&
              rdiscretisation=rnonlinearCCMatrix%p_rdiscretisationStabil%RspatialDiscr(1))

          if (.not. bshared) then
            call output_line (&
                'Edge oriented stabilisation does not support independent A11/A22!', &
                OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
            call sys_halt()
          end if

        case (CCMASM_STAB_FASTEDGEORIENTED)
          ! Fast Jump stabilisation. Precomputed matrix.
          
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%dtheta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          if (rnonlinearCCMatrix%dnewton .eq. 0.0_DP) then

            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there's content at all),
            ! these matrices are ignored then by the kernel.
            
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
            
          else

            ! Clear A12/A21 that receives parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they aren't.
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          end if
         
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)
        
          ! Subtract the stabilisation matrix stuff.
          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixStabil, &
              rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixStabil, &
              rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)
              
        ! no stabilisation, bilinearform evaluation
        case(CCMASM_NOSTAB_BILINEARFORM)
        
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
                              
          ! Call the modified SD method to calculate the nonlinearity.
          ! rvector is a block vector with the
          ! X- and Y-velocity.
          ! p_rvelX1 => rvector%RvectorBlock(1)
          ! p_rvelY1 => rvector%RvectorBlock(2)
          !
          ! we need to make the matrix multiplication
          ! b=b-(rho*u*nabla)*u = b-N(u)*u
          ! and then return in rdefect the defect vector
          call cc_nostabBilin (rvector, rvector, &
                               rstreamlineDiffusion, CONV_MODDEFECT, &
                               rmatrix,rproblem,rsolution=rvector,rdefect=rdefect)
                               
          call vecio_writeBlockVectorHR (rdefect, 'def', .false.,&
                                    0, 'defnostab.txt', '(E20.5)')
              
              
                      

        case default
          call output_line ('Don''t know how to set up nonlinearity!?!', &
              OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
          call sys_halt()
        
        end select
      
      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
        !
        ! Type of stablilisation?
        select case (rnonlinearCCMatrix%iupwind)
        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rnonlinearCCMatrix%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
              rsolution=rvector,rdefect=rdefect,&
              rdiscretisation=rnonlinearCCMatrix%p_rdiscretisationStabil%RspatialDiscr(1))

        case (CCMASM_STAB_FASTEDGEORIENTED)
          ! Fast Jump stabilisation. Precomputed matrix.
          
          ! Subtract the stabilisation matrix stuff.
          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixStabil, &
              rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixStabil, &
              rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

        case default
          ! No stabilisation
        
        end select
      
      end if ! gamma <> 0
    
    end subroutine

  end subroutine
  
! ***************************************************************************

!<subroutine>

  subroutine cc_nostabBilin(rvecPrimary, rvecSecondary,&
                           rconfig, cdef, &
                           rmatrix, rproblem, rsolution ,rdefect)

!<description>
  ! Standard streamline diffusion method to set up the operator
  !
  ! $$ dtheta  *  (                dalpha * MASS
  !                  +              dbeta * STOKES
  !                  +             ddelta * u_1 * grad(.)
  !                  +            dnewton * (.) * grad(u_1)
  !                  +   ddeltaTransposed * grad(.)^T * u_1
  !                  +  dnewtonTransposed * grad(u_1)^T * (.) ) $$
  !
  ! in a matrix or to build a defect vector.
  ! 2D-version (X- and Y-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  !
  ! $u_2 = rsolution(.)$ defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! rmeshVelocity is an optional mesh velocity field that must be present
  ! if the ALE method should be used.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
  !
  ! The routine works with a block matrix rmatrix and allows to include
  ! extended operators like the Newton operator (Frechet-derivative of the
  ! convective part).
!</description>

!<input>

  ! Primary velocity field for the computation of $u_1$
  type(t_vectorBlock), intent(IN), target :: rvecPrimary
  
  ! Secondary velocity field for the computation of $u_1$
  type(t_vectorBlock), intent(IN), target :: rvecSecondary
  
  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamlineDiffusion), intent(IN) :: rconfig
  
  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(IN) :: cdef

  ! optional: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(IN), target, optional :: rsolution
  
!</input>

!<inputoutput>
  ! System block matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  ! The blocks A11,A12,A21 and A22 of this matrix are tackled by streamline
  ! diffusion.
  type(t_matrixBlock), intent(INOUT) :: rmatrix
  
  type(t_problem), intent(inout) :: rproblem
  
  ! optional: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! $r = r - \theta * u_1*grad(u_2)$
  type(t_vectorBlock), intent(INOUT), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    integer(I32) :: celement
    type(t_vectorScalar), pointer :: p_rvelX1,p_rvelX2,p_rvelY1,p_rvelY2
    type(t_vectorScalar), pointer :: p_rsolX,p_rsolY,p_rdefectX,p_rdefectY
    type(t_matrixScalar) :: rmatrixScalar
    real(DP), dimension(:), pointer :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2
    real(DP), dimension(:), pointer :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY
    ! a bilinearform
    ! The bilinear form specifying the operator
    type(t_bilinearForm) :: rform
    
    ! assign the vector pointer
    rproblem%rcollection%p_rvectorQuickAccess1 => rvecPrimary
    
    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
        print *,'SD: Solution/defect vector not present!'
        call sys_halt()
      end if
    end if
    
    ! Get the actual subvectors from the velocity vectors that define
    ! the X- and Y-velocity.
    p_rvelX1 => rvecPrimary%RvectorBlock(1)
    p_rvelY1 => rvecPrimary%RvectorBlock(2)
    p_rvelX2 => rvecSecondary%RvectorBlock(1)
    p_rvelY2 => rvecSecondary%RvectorBlock(2)
    
    if (present(rsolution)) then
      p_rsolX => rsolution%RvectorBlock(1)
      p_rsolY => rsolution%RvectorBlock(2)
    else
      nullify(p_rsolX)
      nullify(p_rsolY)
    end if
    
    if (present(rdefect)) then
      p_rdefectX => rdefect%RvectorBlock(1)
      p_rdefectY => rdefect%RvectorBlock(2)
    else
      nullify(p_rdefectX)
      nullify(p_rdefectY)
    end if
      
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%RmatrixBlock(1,1)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      print *,'SD: Unsupported matrix format'
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(2,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      print *,'SD: Unsupported matrix format'
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,1,2) .and. &
        (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      print *,'SD: Unsupported matrix format'
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,1,2) .and. &
        (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      print *,'SD: Unsupported matrix format'
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,2,1) .and. &
        (rmatrix%RmatrixBlock(2,1)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(2,1)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      print *,'SD: Unsupported matrix format'
      call sys_halt()
    end if

    ! If Newton must be calculated, make sure A12 and A21 exists and that
    ! all A11, A12, A21 and A22 are independent of each other!
    if ((rconfig%dnewton .ne. 0.0_DP) .or. (rconfig%dnewtonTransposed .ne. 0.0_DP)) then
      if (.not. lsysbl_isSubmatrixPresent(rmatrix,1,2) .or. &
          .not. lsysbl_isSubmatrixPresent(rmatrix,2,1)) then
        print *,'SD: For the Newton matrix, A12 and A21 must be defined!'
        call sys_halt()
      end if
      if (lsyssc_isMatrixContentShared ( &
              rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2)) .or. &
          lsyssc_isMatrixContentShared ( &
              rmatrix%RmatrixBlock(1,2),rmatrix%RmatrixBlock(2,1)) ) then
        print *,'SD: For the Newton matrix, the matrix blocks must be indepentent!'
        call sys_halt()
      end if
    end if

    celement = rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest% &
                RelementDistr(1)%celement
    if (rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%ccomplexity &
        .ne. SPDISC_UNIFORM) then
      print *,'SD: Unsupported discretisation.'
      call sys_halt()
    end if

    if ((rvecPrimary%cdataType .ne. ST_DOUBLE) .or. &
        (rvecSecondary%cdataType .ne. ST_DOUBLE)) then
      print *,'SD: Unsupported vector data type in velocity.'
      call sys_halt()
    end if
    
    if (present(rdefect)) then
      if ((rsolution%cdataType .ne. ST_DOUBLE) .or. &
          (rdefect%cdataType .ne. ST_DOUBLE)) then
        print *,'SD: Unsupported vector data type in solution/defect'
        call sys_halt()
      end if
    end if
    
    if (.not. rconfig%bconstViscosity) then
      print *,'SD: Only constant viscosity supported at the moment!'
      call sys_halt()
    end if
    
    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      print *,'SD: Viscosity parameter nu not initialised!'
      call sys_halt()
    end if
    
    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers don't like that ^^

    call lsyssc_getbase_double (p_rvelX1,p_DvelX1)
    call lsyssc_getbase_double (p_rvelY1,p_DvelY1)
    call lsyssc_getbase_double (p_rvelX2,p_DvelX2)
    call lsyssc_getbase_double (p_rvelY2,p_DvelY2)
    
    !!! DEBUG:
    !WHERE (abs(p_DvelX1) .LT. 1E-12_DP) p_DvelX1 = 0.0_DP
    !WHERE (abs(p_DvelY1) .LT. 1E-12_DP) p_DvelY1 = 0.0_DP
    !call vecio_writeArray_Dble (p_DvelX1, 'vecx1', &
    !                               0, 'vectorx1.txt', '(D10.3)')
    !call vecio_writeArray_Dble (p_DvelY1, 'vecx2', &
    !                               0, 'vectorx2.txt', '(D10.3)')
    
    if (present(rdefect)) then
      ! get the actual entries of the vector
      call lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
      call lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
      call lsyssc_getbase_double (p_rdefectX,p_DdefectX)
      call lsyssc_getbase_double (p_rdefectY,p_DdefectY)
      
      ! now set up a bilinearform
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_FUNC
      !    ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff = .false.
      rform%Dcoefficients(1)  = rconfig%ddelta
      rform%Dcoefficients(2)  = rconfig%ddelta

      ! Here create a matrix for the defect b=b-N(u)*u
      call bilf_createMatrixStructure(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial,&
                                      rmatrix%RmatrixBlock(1,1)%cmatrixFormat,&
                                      rmatrixScalar,&
                                      rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest)


      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Stokes for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      call bilf_buildMatrixScalar (rform,.true.,&
                                   rmatrixScalar,cc_ST,&
                                   rproblem%rcollection)

!      call matio_writeMatrixHR (rmatrixScalar, 'Nu',&
!                                .TRUE., 0, 'matrixNu.txt', '(E20.5)')
                                
      ! perform the actual multiplication: b=b-N(u)*u and write out
      ! the defect vectors:
      ! first the u1-defect
      call lsyssc_scalarMatVec(rmatrixScalar,p_rvelX1,p_rdefectX,-1.0_dp,1.0_dp)
      ! then the u2-defect
      call lsyssc_scalarMatVec(rmatrixScalar,p_rvelY1,p_rdefectY,-1.0_dp,1.0_dp)
      
      ! free the temporary matrix
      call lsyssc_releaseMatrix(rmatrixScalar)
      
      
      
    else
      ! we will set up a bilinearform
      
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_FUNC
      !    ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff = .false.
      rform%Dcoefficients(1)  = rconfig%ddelta
      rform%Dcoefficients(2)  = rconfig%ddelta

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Stokes for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      call bilf_buildMatrixScalar (rform,.false.,&
                                   rmatrix%RmatrixBlock(1,1),cc_ST,&
                                   rproblem%rcollection)

      !!! DEBUG:
!      call matio_writeMatrixHR (rmatrix%RmatrixBlock(1,1), 'matrix',&
!                                .TRUE., 0, 'matrixL.txt', '(D20.5)')
                    
    end if

  end subroutine

  ! ***************************************************************************
  
  subroutine cc_ST (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    ! A list of the corner vertices of all elements in progress.
    ! array [1..dimension,1..#vertices per element,1..Number of elements] of double
    ! real(DP), dimension(:,:,:), pointer           :: p_Dcoords => null()
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
    integer :: iel,icup,ive,nve,iin
    type(t_vectorScalar), pointer :: p_rvectorScalarUx
    type(t_vectorScalar), pointer :: p_rvectorScalarUy
    ! here we store the values of u and rho
    real(dp), dimension(:,:,:), allocatable :: Dvalues
    real(dp) :: dxcenter, dycenter, dradius, ddist,drho,drho1,drho2
    
    type(t_geometryObject), pointer :: p_rgeometryObject
    
    ! get a pointer to the geometry object
    p_rgeometryObject => collct_getvalue_geom (rcollection, 'mini')
    
    ! Definition of the circle
!    dxcenter = 0.5
!    dycenter = 0.5
!    dradius  = 0.15
    
    drho1    = rcollection%Dquickaccess(5)
    drho2    = rcollection%Dquickaccess(6)

    ! get the pointers form the collection
    p_rvectorScalarUx => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)
    p_rvectorScalarUy => rcollection%p_rvectorQuickAccess1%RvectorBlock(2)

    ! allocate memory
    allocate(Dvalues(2,npointsPerElement,nelements))
    ! evaluate u in the cubature points
    call fevl_evaluate_sim4 (p_rvectorScalarUx, &
                             rdomainIntSubset, DER_FUNC, Dvalues, 1)

    call fevl_evaluate_sim4 (p_rvectorScalarUy, &
                             rdomainIntSubset, DER_FUNC, Dvalues, 2)
    
    ! loop over the elements and cubature points
    ! and assign the coefficients
    do iel=1,nelements
      do icup=1,npointsPerElement
      
!        ddist = sqrt( (Dpoints(1,icup,iel) - dxcenter)**2 + (Dpoints(2,icup,iel)-dycenter)**2)
!        if(ddist .le. dradius)then
        call geom_isInGeometry (p_rgeometryObject, Dpoints(:,icup,iel), iin)
        if(iin .eq. 1)then
          Dcoefficients(1,icup,iel) = rform%Dcoefficients(1) * Dvalues(1,icup,iel) * drho2
          Dcoefficients(2,icup,iel) = rform%Dcoefficients(2) * Dvalues(2,icup,iel) * drho2
        else
          Dcoefficients(1,icup,iel) = rform%Dcoefficients(1) * Dvalues(1,icup,iel) * drho1
          Dcoefficients(2,icup,iel) = rform%Dcoefficients(2) * Dvalues(2,icup,iel) * drho1
        end if
      end do
    end do
    
    deallocate(Dvalues)
  end subroutine

  ! ***************************************************************************
  
  subroutine cc_Rho(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
        
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
    integer :: iel,icup,ive,nve,iin
    ! here we store the values of u and rho
    real(dp), dimension(:,:,:), allocatable :: Dvalues
    real(dp) :: dxcenter, dycenter, dradius, ddist,drho,drho1,drho2
    real(dp), dimension(:), allocatable :: DrhoElement
    
    type(t_geometryObject), pointer :: p_rgeometryObject
    
    ! get a pointer to the geometry object
    p_rgeometryObject => collct_getvalue_geom (rcollection, 'mini')
    
    drho1    = rcollection%Dquickaccess(5)
    drho2    = rcollection%Dquickaccess(6)

    ! Definition of the circle
!    dxcenter = 0.5
!    dycenter = 0.5
!    dradius  = 0.15

    ! loop over the elements and cubature points
    ! and assign the coefficients
    do iel=1,nelements
      do icup=1,npointsPerElement
        
!        ddist = sqrt( (Dpoints(1,icup,iel) - dxcenter)**2 + (Dpoints(2,icup,iel)-dycenter)**2)
!        if(ddist .le. dradius)then
        call geom_isInGeometry (p_rgeometryObject, Dpoints(:,icup,iel), iin)
        if(iin .eq. 1)then
          Dcoefficients(1,icup,iel) = drho2
        else
          Dcoefficients(1,icup,iel) = drho1
        end if
      end do
    end do
    
  end subroutine

  ! ***************************************************************************
  

  !****************************************************************************

end module
