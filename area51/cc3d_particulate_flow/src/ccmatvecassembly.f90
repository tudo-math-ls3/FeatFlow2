!##############################################################################
!# ****************************************************************************
!# <name> ccmatvecassembly </name>
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
  use basicgeometry
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use derivatives
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use matrixrestriction
  use linearsystemscalar
  use bilinearformevaluation
  use linearformevaluation
  use trilinearformevaluation
  use matrixio
  use statistics
  use linearsystemblock
  use linearsolver
  use linearsolverautoinitialise
  use collection
  use feevaluation
  
  use convection
    
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

    ! Pointer to a B3-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB3 => null()

    ! Pointer to a B1^T-matrix.
    ! This pointer may point to NULL(). In this case, B1^T is created
    ! by 'virtually transposing' the B1 matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB1T => null()

    ! Pointer to a B2-matrix.
    ! This pointer may point to NULL(). In this case, B2^T is created
    ! by 'virtually transposing' the B2 matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB2T => null()

    ! Pointer to a B3-matrix.
    ! This pointer may point to NULL(). In this case, B3^T is created
    ! by 'virtually transposing' the B3 matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB3T => null()

    ! Pointer to a Mass matrix.
    ! May point to NULL() during matrix creation.
    type(t_matrixScalar), pointer :: p_rmatrixMass => null()

  end type

!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleMatrix (coperation,cmatrixType,rmatrix,&
      rnonlinearCCMatrix,rvector,rfineMatrix)

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
      !    / A11  A12  A13  B1  \
      !    | A21  A22  A23  B2  |
      !    | A21  A22  A33  B3  |
      !    \ B1^T B2^T B3^T     /
      !
      ! Assemble the velocity submatrices
      !
      !    / A11  A12  A13   .  \
      !    | A21  A22  A23   .  |
      !    | A21  A22  A33   .  |
      !    \  .    .    .       /
      
      call assembleVelocityBlocks (&
          rnonlinearCCMatrix,rmatrix,rvector,1.0_DP)
      
      ! Assemble the gradient submatrices
      !
      !    /  .    .    .   B1  \
      !    |  .    .    .   B2  |
      !    |  .    .    .   B3  |
      !    \ B1^T B2^T B3^T     /
      
      call assembleGradientMatrices (rnonlinearCCMatrix,rmatrix,&
        iand(coperation,CMASM_QUICKREFERENCES) .ne. 0)

      ! 2.) Initialise the weights for the B-matrices
      !
      !    /  .    .    .   B1  \
      !    |  .    .    .   B2  |
      !    |  .    .    .   B3  |
      !    \ B1^T B2^T B3^T     /
      
      rmatrix%RmatrixBlock(1,4)%dscaleFactor = rnonlinearCCMatrix%deta
      rmatrix%RmatrixBlock(2,4)%dscaleFactor = rnonlinearCCMatrix%deta
      rmatrix%RmatrixBlock(3,4)%dscaleFactor = rnonlinearCCMatrix%deta
      
      rmatrix%RmatrixBlock(4,1)%dscaleFactor = rnonlinearCCMatrix%dtau
      rmatrix%RmatrixBlock(4,2)%dscaleFactor = rnonlinearCCMatrix%dtau
      rmatrix%RmatrixBlock(4,3)%dscaleFactor = rnonlinearCCMatrix%dtau
      
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

        if (.not. lsyssc_isMatrixContentShared(rmatrix%RmatrixBlock(3,3))) then
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(3,3), &
              rmatrix%RmatrixBlock(3,3), &
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
      type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient

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

      ! In the global system, there are three gradient matrices B1, B2 and B3.
      ! Get a pointer to the template structure for these.
      ! If there is no pointer, try to get use a pointer to one of these
      ! matrices directly.
      p_rmatrixTemplateGradient => rnonlinearCCMatrix%p_rmatrixTemplateGradient
      if (.not. associated(p_rmatrixTemplateGradient)) &
        p_rmatrixTemplateGradient => rnonlinearCCMatrix%p_rmatrixB1
      if (.not. associated(p_rmatrixTemplateGradient)) then
        call output_line ('Cannot set up B matrices in system matrix!', &
            OU_CLASS_ERROR,OU_MODE_STD,'allocMatrix')
        call sys_halt()
      end if

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      if (associated(p_rdiscretisation)) then
        call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
      else
        ! No discretisation structure; create the matrix directly as 4x4 matrix.
        call lsysbl_createEmptyMatrix (rmatrix,NDIM3D+1)
      end if
        
      ! Let's consider the global system in detail. The standard matrix It has
      ! roughly the following shape:
      !
      !    / A11            B1 \   / A11  A12  A13  A14 \
      !    |      A22       B2 | = | A21  A22  A23  A24 |
      !    |           A33  B2 |   | A31  A32  A33  A34 |
      !    \ B1^T B2^T B3^T  . /   \ A41  A42  A43  A44 /
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
          
      if ((.not. bdecoupled) .and. (.not. bfulltensor)) then
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content. Same thing for A33.
        call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
        call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      else
      
        ! Otherwise, create another copy of the template matrix.
        call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                    
        call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      end if
      
      ! Manually change the discretisation structure of the Y-/Z-velocity
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))
      call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(3,3),&
          p_rdiscretisation%RspatialDiscr(3))
                                          
      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      if (bfulltensor) then

        ! We have a matrix in the following shape:
        !
        !    / A11  A12  A13  B1  \
        !    | A21  A22  A23  B2  |
        !    | A21  A22  A33  B3  |
        !    \ B1^T B2^T B3^T     /
        !
        ! Create Aij with i /= j
      
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

        if (rmatrix%RmatrixBlock(1,3)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(1,3), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(1,3),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(2,1), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
            
        end if
        
        if (rmatrix%RmatrixBlock(2,3)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(2,3), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(2,3),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rmatrix%RmatrixBlock(3,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(3,1), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(3,1),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rmatrix%RmatrixBlock(3,2)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(3,2), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(3,2),LSYSSC_SETM_UNDEFINED)
            
        end if

      end if

      ! The B1/B2/B3 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2/B3 with those B1/B2/B3 of the
      ! block matrix, while we create empty space for the entries.
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,4),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,4),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB3, &
                                    rmatrix%RmatrixBlock(3,4),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
      ! Now, prepare B1^T, B2^T and B3^T. Are these matrices given?
      !
      ! If yes, take the given matrices. If not, create them by
      ! 'virtually transposing' B1, B2 and B3 (i.e. they share the same
      ! data as B1, B2 and B3 but hate the 'transpose'-flag set).
      
      if (associated(rnonlinearCCMatrix%p_rmatrixB1T)) then
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1T, &
                                      rmatrix%RmatrixBlock(4,1),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      else
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB1, &
                                      rmatrix%RmatrixBlock(4,1),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      if (associated(rnonlinearCCMatrix%p_rmatrixB2T)) then
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2T, &
                                      rmatrix%RmatrixBlock(4,2),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      else
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB2, &
                                      rmatrix%RmatrixBlock(4,2),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      if (associated(rnonlinearCCMatrix%p_rmatrixB3T)) then
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB3T, &
                                      rmatrix%RmatrixBlock(4,3),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      else
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB3, &
                                      rmatrix%RmatrixBlock(4,3),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      ! That's it, all submatrices are basically set up.
      !
      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      call lsysbl_updateMatStrucInfo (rmatrix)
        
    end subroutine
    
    ! -----------------------------------------------------

    subroutine assembleVelocityBlocks (rnonlinearCCMatrix,rmatrix,rvector,dvectorWeight)
        
    ! Assembles the velocity matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    
    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(IN) :: rnonlinearCCMatrix
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
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
    type(T_jumpStabilisation) :: rjumpStabil
    real(DP) :: dvecWeight
    
      ! Standard value for dvectorWeight is = -1.
      dvecWeight = -1.0_DP
      if (present(dvectorWeight)) dvecWeight = dvectorWeight
    
      ! Is A11=A22 physically?
      ! We either share both A22 and A33 with A11 or none of them -
      ! so it's sufficient to check whether A22 is shared with A11...
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      ! Allocate memory if necessary. Normally this should not be necessary...
      ! A11:
      if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) then
        call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
      end if
    
      ! A22,A33:
      if (.not. bshared) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
        end if
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(3,3))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(3,3),LSYSSC_SETM_UNDEFINED)
        end if
      end if

      ! A12,A21:
      if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,2))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
        end if
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
        end if
      end if
      ! A13,A31
      if (lsysbl_isSubmatrixPresent (rmatrix,1,3)) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,3))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,3),LSYSSC_SETM_UNDEFINED)
        end if
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(3,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(3,1),LSYSSC_SETM_UNDEFINED)
        end if
      end if
      ! A23,A32
      if (lsysbl_isSubmatrixPresent (rmatrix,2,3)) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,3))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,3),LSYSSC_SETM_UNDEFINED)
        end if
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(3,2))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(3,2),LSYSSC_SETM_UNDEFINED)
        end if
      end if
    
      ! ---------------------------------------------------
      ! Plug in the mass matrix?
      if (rnonlinearCCMatrix%dalpha .ne. 0.0_DP) then
       
        ! Allocate memory if necessary. Normally this should not be necessary...
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
        end if
      
        call lsyssc_matrixLinearComb (&
            rnonlinearCCMatrix%p_rmatrixMass       ,rnonlinearCCMatrix%dalpha,&
            rmatrix%RmatrixBlock(1,1),0.0_DP,&
            rmatrix%RmatrixBlock(1,1),&
            .false.,.false.,.true.,.true.)
            
        if (.not. bshared) then

          ! Allocate memory if necessary. Normally this should not be necessary...
          if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) then
            call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
          end if
          if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(3,3))) then
            call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(3,3),LSYSSC_SETM_UNDEFINED)
          end if

          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rmatrixMass     ,rnonlinearCCMatrix%dalpha,&
              rmatrix%RmatrixBlock(2,2),0.0_DP,&
              rmatrix%RmatrixBlock(2,2),&
              .false.,.false.,.true.,.true.)

          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rmatrixMass     ,rnonlinearCCMatrix%dalpha,&
              rmatrix%RmatrixBlock(3,3),0.0_DP,&
              rmatrix%RmatrixBlock(3,3),&
              .false.,.false.,.true.,.true.)
        end if
        
      else
      
        ! Otherwise, initialise the basic matrix with 0.0
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
        
        if (.not. bshared) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,3))
        end if
          
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
          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rmatrixStokes   ,rnonlinearCCMatrix%dtheta,&
              rmatrix%RmatrixBlock(3,3),1.0_DP,&
              rmatrix%RmatrixBlock(3,3),&
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
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%dupsam
          
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          if (rnonlinearCCMatrix%dnewton .eq. 0.0_DP) then
          
            ! If the submatrices Aij, i /= j,  exist, fill them with zero.
            ! If they don't exist, we don't have to do anything.
            if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
            end if
            if (lsysbl_isSubmatrixPresent (rmatrix,1,3)) then
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,3))
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,1))
            end if
            if (lsysbl_isSubmatrixPresent (rmatrix,2,3)) then
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,3))
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,2))
            end if
            
         else

            ! Clear Aij, i /= j, that may receive parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,3))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,1))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,3))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,2))
          
         end if
         
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk3d (rvector, rvector, &
              dvecWeight,0.0_DP,rstreamlineDiffusion,CONV_MODMATRIX,rmatrix)

        case (CCMASM_STAB_UPWIND)

          ! Upwinding in 3D has still to be implemented...
          print *, 'ERROR: upwind not implemented in 3D yet!'
          stop

!          ! Set up the upwind structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rupwind%dnu = rnonlinearCCMatrix%dnu
!
!          ! Set stabilisation parameter
!          rupwind%dupsam = rnonlinearCCMatrix%dupsam
!
!          ! Matrix weight
!          rupwind%dtheta = rnonlinearCCMatrix%dgamma
!
!          IF (rnonlinearCCMatrix%dnewton .NE. 0.0_DP) THEN
!            CALL output_line ('Warning: Upwind does not support assembly '&
!                //'of the Newton matrix!',OU_CLASS_TRACE1)
!          END IF
!
!          ! Call the upwind method to calculate the nonlinear matrix.
!          CALL conv_upwind2d (rvector, rvector, &
!                              dvecWeight, 0.0_DP,&
!                              rupwind, CONV_MODMATRIX, &
!                              rmatrix%RmatrixBlock(1,1))
!
!          IF (.NOT. bshared) THEN
!            ! Modify also the matrix block (2,2)
!            !CALL conv_upwind2d (rvector, rvector, &
!            !                    dvecWeight, 0.0_DP,&
!            !                    rupwind, CONV_MODMATRIX, &
!            !                    rmatrix%RmatrixBlock(2,2))
!          END IF

        case DEFAULT
          call output_line ('Don''t know how to set up nonlinearity!?!', &
              OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityBlocks')
          call sys_halt()
        
        end select

      end if
    
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
        if ((.not. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(1,4))) .or.&
            (.not. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(2,4))) .or.&
            (.not. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(3,4)))) then
          idubContent = LSYSSC_DUP_COPY
        end if
        
      else
      
        idubContent = LSYSSC_DUP_COPY
        
      end if
      
      idubStructure = LSYSSC_DUP_SHARE
      
      ! Let's consider the global system in detail:
      !
      !    / A11  A12  A13  B1  \   / A11  A12  A13  A14 \
      !    | A21  A22  A23  B2  | = | A21  A22  A23  A24 |
      !    | A31  A32  A33  B3  |   | A31  A32  A33  A34 |
      !    \ B1^T B2^T B3^T     /   \ A41  A42  A43  A44 /
      !
      ! We exclude the velocity submatrices here, so our system looks like:
      !
      !    /                B1  \   /                A14 \
      !    |                B2  | = |                A24 |
      !    |                B3  |   |                A34 |
      !    \ B1^T B2^T B3^T     /   \ A41  A42  A43  A44 /

      ! The B1/B2/B3 matrices exist up to now only in rnonlinearCCMatrix.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2/B3 with those B1/B2/B3 of the
      ! block matrix, while we create copies of the entries. The B-blocks
      ! are already prepared and memory for the entries is already allocated;
      ! so we only have to copy the entries.
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,4),&
                                    idubStructure,idubContent)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,4),&
                                    idubStructure,idubContent)
      
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB3, &
                                    rmatrix%RmatrixBlock(3,4),&
                                    idubStructure,idubContent)

      ! Now, prepare B1^T, B2^T and B3^T. Are these matrices given?
      !
      ! If yes, take the given matrices. If not, create them by
      ! 'virtually transposing' B1, B2 and B3 (i.e. they share the same
      ! data as B1 and B2 but hate the 'transpose'-flag set).
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      if (associated(rnonlinearCCMatrix%p_rmatrixB1T)) then
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1T, &
                                      rmatrix%RmatrixBlock(4,1),&
                                      idubStructure,idubContent)
      else
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB1, &
                                      rmatrix%RmatrixBlock(4,1),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      if (associated(rnonlinearCCMatrix%p_rmatrixB2T)) then
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2T, &
                                      rmatrix%RmatrixBlock(4,2),&
                                      idubStructure,idubContent)
      else
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB2, &
                                      rmatrix%RmatrixBlock(4,2),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      if (associated(rnonlinearCCMatrix%p_rmatrixB3T)) then
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB3T, &
                                      rmatrix%RmatrixBlock(4,3),&
                                      idubStructure,idubContent)
      else
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB3, &
                                      rmatrix%RmatrixBlock(4,3),&
                                      LSYSSC_TR_VIRTUAL)
      end if

    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_nonlinearMatMul (rnonlinearCCMatrix,rx,rd,dcx,dcd,ry)

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
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_ry
    type(t_matrixBlock) :: rmatrix
    
    p_ry => rx
    if (present(ry)) p_ry => ry
    
    ! Probably weight the input vector.
    if (dcd .ne. 1.0_DP) then
      call lsysbl_scaleVector (rd,dcd)
    end if
    
    ! The system matrix looks like:
    !
      !    / A11  A12  A13  B1  \
      !    | A21  A22  A23  B2  |
      !    | A21  A22  A33  B3  |
      !    \ B1^T B2^T B3^T     /
    !
    ! Create a temporary matrix that covers this structure.
    call lsysbl_createEmptyMatrix (rmatrix,NDIM3D+1)
    
    ! Put references to the Stokes- and B-matrices to Aij. assembleVelocityDefect
    ! needs this template matrix to provide the structure for the stabilisation
    ! routines! The B-matrices are needed later.
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
    if (rnonlinearCCMatrix%dnewton .ne. 0.0_DP) then
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    end if
    
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB1,&
        rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB2,&
        rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rmatrixB3,&
        rmatrix%RmatrixBlock(3,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB1, &
                                  rmatrix%RmatrixBlock(4,1),&
                                  LSYSSC_TR_VIRTUAL)

    call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB2, &
                                  rmatrix%RmatrixBlock(4,2),&
                                  LSYSSC_TR_VIRTUAL)

    call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rmatrixB3, &
                                  rmatrix%RmatrixBlock(4,3),&
                                  LSYSSC_TR_VIRTUAL)

    ! Update the structural information of the block matrix, as we manually
    ! changed the submatrices:
    call lsysbl_updateMatStrucInfo (rmatrix)
    
    ! In the first step, we assemble the defect that arises in the velocity
    ! components. This is characterised by the following submatrix:
    !
    !    / A11  A12  A13   .  \
    !    | A21  A22  A23   .  |
    !    | A21  A22  A33   .  |
    !    \  .    .    .       /
    !
    ! assembleVelocityDefect handles exactly these submatrices.

    call assembleVelocityDefect (rnonlinearCCMatrix,rmatrix,rx,rd,p_ry,-dcx)
    
    ! Now, we treat all the remaining blocks. Let's see what is missing:
    !
    !    /  .    .    .   B1  \
    !    |  .    .    .   B2  |
    !    |  .    .    .   B3  |
    !    \ B1^T B2^T B3^T     /

    ! To build the appropriate defect, we first remove the velocity blocks:
    
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,2))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,3))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,2))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,3))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(3,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(3,2))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(3,3))

    ! Initialise the weights for the B/B^T matrices
    rmatrix%RmatrixBlock(1,4)%dscaleFactor = rnonlinearCCMatrix%deta
    rmatrix%RmatrixBlock(2,4)%dscaleFactor = rnonlinearCCMatrix%deta
    rmatrix%RmatrixBlock(3,4)%dscaleFactor = rnonlinearCCMatrix%deta
    
    rmatrix%RmatrixBlock(4,1)%dscaleFactor = rnonlinearCCMatrix%dtau
    rmatrix%RmatrixBlock(4,2)%dscaleFactor = rnonlinearCCMatrix%dtau
    rmatrix%RmatrixBlock(4,3)%dscaleFactor = rnonlinearCCMatrix%dtau

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
        rmatrix,rvector,rdefect,rvelocityVector,dvectorWeight)
        
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
    type(T_jumpStabilisation) :: rjumpStabil
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))

      ! ---------------------------------------------------
      ! Subtract the mass matrix stuff?
      if (rnonlinearCCMatrix%dalpha .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixMass, &
            rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixMass, &
            rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixMass, &
            rvector%RvectorBlock(3), rdefect%RvectorBlock(3), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)
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

        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rmatrixStokes, &
            rvector%RvectorBlock(3), rdefect%RvectorBlock(3), &
            -rnonlinearCCMatrix%dtheta, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if (rnonlinearCCMatrix%dgamma .ne. 0.0_DP) then
      
        ! Type of stablilisation?
        select case (rnonlinearCCMatrix%iupwind)
        case (0)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%dupsam
          
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
          
          call conv_streamlineDiffusionBlk3d (&
                              rvelocityVector, rvelocityVector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)
                              
        case (1)

          ! Upwinding in 3D has still to be implemented...
          print *, 'ERROR: upwind not implemented in 3D yet!'
          stop

!          ! Set up the upwind structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rupwind%dnu = rnonlinearCCMatrix%dnu
!
!          ! Set stabilisation parameter
!          rupwind%dupsam = rnonlinearCCMatrix%dupsam
!
!          ! Matrix weight
!          rupwind%dtheta = rnonlinearCCMatrix%dgamma
!
!          ! Call the upwind method to calculate the nonlinear defect.
!          !CALL conv_upwind2d (rvector, rvector, &
!          !                    dvectorWeight, 0.0_DP,&
!          !                    rupwind, CONV_MODDEFECT, &
!          !                    rmatrix%RmatrixBlock(1,1),rvector,rdefect)
!
!          IF (.NOT. bshared) THEN
!            CALL output_line ('Upwind does not support independent A11/A22!', &
!                OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
!            CALL sys_halt()
!          END IF

        case DEFAULT
          call output_line ('Don''t know how to set up nonlinearity!?!', &
              OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
          call sys_halt()
        
        end select

      end if ! gamma <> 0
    
    end subroutine

  end subroutine

end module
