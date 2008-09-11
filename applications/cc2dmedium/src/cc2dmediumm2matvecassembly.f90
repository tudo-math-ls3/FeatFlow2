!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2matrixassembly </name>
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
!#   $$ A_1 = \alpha M  +  \theta L  +  \gamma N(y) $$
!#
!# and
!#
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity includung stabilisation,
!#   
!#   $\alpha$ = 0/1     - switches the mass matrix on/off;
!#                          =0 for stationary problem,
!#   $\theta$           - weight for the Laplace matrix,
!#   $\gamma$ = 0/1     - Switches the nonlinearity on/off;
!#                          =0 for Stokes system,
!#   $\eta$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau$   = 0/1     - Switches the 'B^T'-term on/off,
!#
!# To assemble such a matrix, the application has to follow two steps:
!# 
!# 1.) Create a structure of type t_ccmatrixComponents and set the
!#     parameters in this structure according to the matrix which should
!#     be assembled.
!#
!# 2.) Call the matrix/vector assembly routine with this structure as input.
!#
!# The matrix assembly routine will return an initialised matrix structure
!# of a system matrix based on the provided parameters.
!#
!# The module contains the following routines:
!#
!# 1.) c2d2_assembleMatrix
!#     -> Assembles a matrix based on a set of input parameters.
!#
!# 2.) c2d2_assembleDefect
!#     -> Set up a defect vector d:=b-A(x)x
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2matvecassembly

  USE fsystem
  USE storage
  USE linearsystemblock
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
  USE matrixio
  
  USE convection
    
  IMPLICIT NONE
  
!<constants>

!<constantblock description="Identifiers for the 'coperation' input parameter of the matrix assembly routine">

  ! Allocate memory if necessary.
  INTEGER(I32), PARAMETER :: CCMASM_ALLOCMEM              = 1
  
  ! Compute all matrix entries.
  INTEGER(I32), PARAMETER :: CCMASM_COMPUTE               = 2
  
  ! Allocate memory and compute matrix entries.
  INTEGER(I32), PARAMETER :: CCMASM_ALLOCANDCOMPUTE       = 3
  
  ! Bypass memory allocation for matrices.
  INTEGER(I32), PARAMETER :: CMASM_QUICKREFERENCES        = 4
  
!</constantblock>

!<constantblock description="Identifiers for the IUPWIND parameter that specifies how to set up the nonlinearity or stabilisation.">

  ! Streamline diffusion; configured by dupsam
  INTEGER, PARAMETER :: CCMASM_STAB_STREAMLINEDIFF    = 0

  ! 1st-order upwind; configured by dupsam
  INTEGER, PARAMETER :: CCMASM_STAB_UPWIND            = 1
  
  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'
  INTEGER, PARAMETER :: CCMASM_STAB_EDGEORIENTED      = 2

!</constantblock>

!<constantblock description="Matrix type ID's specifying the general matrix class to set up.">

  ! Standard matrix.
  INTEGER, PARAMETER :: CCMASM_MTP_AUTOMATIC         = 0
  
  ! Standard matrix with decoupled velocity blocks
  INTEGER, PARAMETER :: CCMASM_MTP_DECOUPLED         = 1
  
  ! Extended 'full-tensor' matrix with submatrices A11, A12, A21, A22, all independent from
  ! each other.
  INTEGER, PARAMETER :: CCMASM_MTP_FULLTENSOR        = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This structure provides a set of input parametes for the matrix assembly
  ! routine.
  TYPE t_ccmatrixComponents
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    REAL(DP) :: dalpha = 0.0_DP
    
    ! THETA-parameter that controls the weight of the Stokes matrix
    ! in the core equation. =1.0 for stationary simulations.
    REAL(DP) :: dtheta = 0.0_DP
    
    ! GAMMA-parameter that controls the weight in front of the
    ! nonlinearity N(u). =1.0 for Navier-Stokes, =0.0 for Stokes equation.
    REAL(DP) :: dgamma = 0.0_DP

    ! ETA-parameter that switch the B-term on/off in the matrix.
    REAL(DP) :: deta = 0.0_DP
    
    ! TAU-parameter that switch the B^T-term on/off in the matrix.
    REAL(DP) :: dtau = 0.0_DP

    ! Weight for the Newton matrix N*(u).
    ! = 0.0 deactivates the Newton part.
    REAL(DP) :: dnewton = 0.0_DP

    ! STABILISATION: Parameter that defines how to set up the nonlinearity and 
    ! whether to use some kind of stabilisation. One of the CCMASM_STAB_xxxx 
    ! constants. Standard is CCMASM_STAB_STREAMLINEDIFF.
    INTEGER :: iupwind = CCMASM_STAB_STREAMLINEDIFF
    
    ! STABILISATION: Viscosity parameter. Used for stabilisation schemes when 
    ! a nonlinearity is set up.
    REAL(DP) :: dnu = 0.0_DP
    
    ! STABILISATION: Stabilisation parameter for streamline diffusion, upwind and 
    ! edge-oriented stabilisation. If iupwind=CCMASM_STAB_STREAMLINEDIFF, a value of 
    ! 0.0 deactivates any stabilisation.
    REAL(DP) :: dupsam = 0.0_DP
    
    ! MATRIX RESTRICTION: Parameter to activate matrix restriction.
    ! Can be used to generate parts of the matrices on coarse grids where the
    ! aspect ratio of the cells is large. Only applicable for $\tilde Q_1$
    ! discretisations.
    ! Standard = 0 = deactivate matrix restriction
    INTEGER :: iadaptiveMatrices = 0
    
    ! MATRIX RESTRICTION: Threshold parameter for adaptive matrix generation
    ! of coarse grid matrices (here: maximum aspect ratio). 
    ! Only applicable if iadaptiveMatrices <> 0.
    ! Standard = 20.0
    REAL(DP) :: dadmatthreshold = 20.0_DP
    
    ! A pointer to a block matrix that provides preallocated memory.
    ! If not associated, new memory will be allocated by the assembly
    ! routine automatically.
    TYPE(t_matrixBlock), POINTER :: p_rpreallocatedMatrix => NULL()
  
    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...).
    ! Only used during matrix creation if p_rpreallocatedMatrix=>NULL(), 
    ! otherwise this is not used.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation => NULL()

    ! Pointer to a template FEM matrix that defines the structure of 
    ! Laplace/Stokes/... matrices. Only used during matrix creation
    ! if p_rpreallocatedMatrix=>NULL(), otherwise this is not used.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixTemplateFEM => NULL()

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2) matrices. Only used during matrix creation
    ! if p_rpreallocatedMatrix=>NULL(), otherwise this is not used.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixTemplateGradient => NULL()

    ! Pointer to Stokes matrix (=nu*Laplace). 
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes => NULL()

    ! Pointer to a B1-matrix.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB1 => NULL()

    ! Pointer to a B2-matrix.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB2 => NULL()

    ! Pointer to a Mass matrix.
    ! May point to NULL() during matrix creation.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixMass => NULL()

  END TYPE

!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_assembleMatrix (coperation,cmatrixType,rmatrix,rmatrixComponents,&
      rvector,rfineMatrix)

!<description>
  ! This routine assembles a global matrix. The caller must initialise the 
  ! rmatrixComponents according to how the matrix should look like.
  ! The 'coperation' parameter tells the routine what to do.
  ! The destination matrix rmatrix is then set up or updated.
  !
  ! The parameters rvector and rfineMatrix are optional. rvector must be
  ! specified, if the nonlinearity is activated (parameter $\gamma\not=0$ in 
  ! rmatrixComponents). This vector specifies the 'solution' where the nonlinearity 
  ! $u\nabla u$ is evaluated.
  ! rfineMatrix allows to specify a matrix of a 'one level refined mesh'. This
  ! is usually used when setting up preconditioners over multiple levels.
  ! Specifying such a matrix allows the routine (depending on the discretisation)
  ! to include some special stabilisation terms into the matrix rmatrix.
  !
  ! The routine will not include any boundary conditions in the matrix.
  ! Nevertheless, if the p_rpreallocatedMatrix parameter in rmatrixComponents
  ! is set and rmatrix is not initialised, a reference to boundary condition 
  ! structure of p_rpreallocatedMatrix is transferred to rmatrix.
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
  ! that references to the original gradient (B-)matrices from rmatrixComponents
  ! are written to rmatrix; matrix entries are not copied!
  ! (This can be used e.g. for setting up a matrix for building a defect
  !  vector without copying matrix data.)
  ! In this case, the caller MUST NOT CHANGE rmatrix in any way, otherwise
  ! the original (template) matrices would be changed!
  INTEGER(I32), INTENT(IN) :: coperation

  ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
  ! constants.
  ! Usually, CCMASM_MTP_AUTOMATIC is used here. This will automatically determine
  ! the 'smallest possible' matrix structure fitting the needs of the input
  ! parameters. By specifying another matrix type, the caller can explicitly
  ! take influence on the general matrix structure.
  !
  ! If the matrix already exists and only its entries are to be computed,
  ! CCMASM_MTP_AUTOMATIC should be specified here.
  INTEGER, INTENT(IN) :: cmatrixType

  ! A t_ccmatrixComponents structure providing all necessary 'source' information
  ! about how to set up the matrix. 
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, either
  ! p_rpreallocatedMatrix must be initialised or p_rmatrixTemplateXXXX
  ! as well as p_rdiscretisation!
  ! If the p_rpreallocatedMatrix pointer is associated and 'coperation' the 
  ! routine will set up rmatrix based on p_rpreallocatedMatrix. It will create
  ! rmatrix as a 'shared copy', i.e. in such a way that rmatrix shares its 
  ! information with p_rpreallocatedMatrix. That helps to save memory and 
  ! computational time for allocation/deallocation of memory!
  ! If p_rpreallocatedMatrix is not set, the new matrix is created based
  ! p_rmatrixTemplateXXXX as well as p_rdiscretisation.
  ! In both cases, memory is allocated automatically if it's missing.
  TYPE(t_ccmatrixComponents), INTENT(IN) :: rmatrixComponents

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity.
  TYPE(t_vectorBlock), INTENT(IN), OPTIONAL :: rvector

  ! OPTIONAL: This parameter allows to specify a 'fine grid matrix'. This is 
  ! usually done when assembling matrices on multiple levels. If specified, the
  ! routine will (if possible) try to include a level-dependent stabilisation
  ! term into the matrix (-> e.g. constant matrix restriction for nonparametric
  ! Rannacher-Turek element if cells are too anisotropic).
  TYPE(t_matrixBlock), INTENT(IN), OPTIONAL :: rfineMatrix
  
!</input>

!<inputoutput>

  ! The destination matrix which should be set up.
  ! If not initialised, a new matrix is created (as if CCMASM_ALLOCxxxx 
  ! was specified).
  ! If initialised, the existing matrix is updated or recreated, depending on
  ! coperation.
  TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
  
!</inputoutput>
  
!</subroutine>

    ! local variables
    LOGICAL :: ballocate
    
    ballocate = .FALSE.
    IF ((rmatrix%NEQ .LE. 0) .OR. &
        IAND(coperation,CCMASM_ALLOCMEM) .NE. 0) THEN
      ballocate = .TRUE.
    END IF
    
    ! What should we do? Allocate memory?
    IF (ballocate) THEN
    
      ! Release the matrix if present.
      CALL lsysbl_releaseMatrix (rmatrix)
    
      ! Allocate a new one. Do we have a template matrix to use as source?
      IF (ASSOCIATED(rmatrixComponents%p_rpreallocatedMatrix)) THEN
        ! Copy the template, share everything -- structure as well as content
        ! (if there is any).
        CALL lsysbl_duplicateMatrix (rmatrixComponents%p_rpreallocatedMatrix,&
            rmatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      ELSE
        ! Oops we have to create a complete new one. That task is a little bit
        ! more advanced... see below!
        CALL allocMatrix (cmatrixType,rmatrixComponents,rmatrix)
      END IF
    END IF   
   
    IF (IAND(coperation,CCMASM_COMPUTE) .NE. 0) THEN

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
      
      CALL assembleVelocityBlocks (&
          rmatrixComponents,rmatrix,rvector,1.0_DP)
      
      ! Assemble the gradient submatrices
      !          
      !    ( .    .    B1  ) 
      !    ( .    .    B2  ) 
      !    ( B1^T B2^T .   ) 
      
      CALL assembleGradientMatrices (rmatrixComponents,rmatrix,&
        IAND(coperation,CMASM_QUICKREFERENCES) .NE. 0)

      ! 2.) Initialise the weights for the B-matrices
      !
      !    (  .    .   B1  ) 
      !    (  .    .   B2  ) 
      !    ( B1^T B2^T  .  )
      
      rmatrix%RmatrixBlock(1,3)%dscaleFactor = rmatrixComponents%deta
      rmatrix%RmatrixBlock(2,3)%dscaleFactor = rmatrixComponents%deta
      
      rmatrix%RmatrixBlock(3,1)%dscaleFactor = rmatrixComponents%dtau
      rmatrix%RmatrixBlock(3,2)%dscaleFactor = rmatrixComponents%dtau
      
      ! Matrix restriction
      ! ---------------------------------------------------
      !
      ! For the construction of matrices on lower levels, call the matrix
      ! restriction. In case we have a uniform discretisation with Q1~,
      ! iadaptivematrix is <> 0 and so this will rebuild some matrix entries
      ! by a Galerkin approach using constant prolongation/restriction.
      ! This helps to stabilise the solver if there are elements in the
      ! mesh with high aspect ratio.
      IF (PRESENT(rfineMatrix)) THEN
        CALL mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(1,1), &
            rmatrix%RmatrixBlock(1,1), &
            rmatrixComponents%iadaptiveMatrices, rmatrixComponents%dadmatthreshold)
            
        IF (.NOT. lsyssc_isMatrixContentShared(rmatrix%RmatrixBlock(2,2))) THEN
          CALL mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(2,2), &
              rmatrix%RmatrixBlock(2,2), &
              rmatrixComponents%iadaptiveMatrices, rmatrixComponents%dadmatthreshold)
        END IF
      END IF
      
    END IF
    
  CONTAINS
  
    ! -----------------------------------------------------
  
    SUBROUTINE allocMatrix (cmatrixType,rmatrixComponents,rmatrix)
    
    ! Allocates memory for the system matrix. rmatrixComponents provides information
    ! about the submatrices that are 'plugged into' rmatrix.
    ! Therefore, before this routine is called, rmatrixComponents must have been set up.

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    INTEGER, INTENT(IN) :: cmatrixType

    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    TYPE(t_ccmatrixComponents), INTENT(IN), TARGET :: rmatrixComponents

    ! A block matrix that receives the basic system matrix.
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
    
      ! local variables
      LOGICAL :: bdecoupled, bfulltensor

      ! A pointer to the system matrix and the RHS/solution vectors.
      TYPE(t_matrixScalar), POINTER :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient

      ! A pointer to the discretisation structure with the data.
      TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
      ! Determine the shape of the matrix
      bdecoupled = cmatrixType .EQ. CCMASM_MTP_DECOUPLED
      bfulltensor = cmatrixType .EQ. CCMASM_MTP_FULLTENSOR
      
      IF (cmatrixType .EQ. CCMASM_MTP_AUTOMATIC) THEN
        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
        bfulltensor = rmatrixComponents%dnewton .NE. 0.0_DP
      END IF
    
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rmatrixComponents%p_rdiscretisation
      
      ! Get a pointer to the template FEM matrix. If that doesn't exist,
      ! take the Stokes matrix as template.
      p_rmatrixTemplateFEM => rmatrixComponents%p_rmatrixTemplateFEM
      IF (.NOT. ASSOCIATED(p_rmatrixTemplateFEM)) &
        p_rmatrixTemplateFEM => rmatrixComponents%p_rmatrixStokes
      IF (.NOT. ASSOCIATED(p_rmatrixTemplateFEM)) THEN
        PRINT *,'allocMatrix: Cannot set up A matrices in system matrix!'
        STOP
      END IF

      ! In the global system, there are two gradient matrices B1 and B2.
      ! Get a pointer to the template structure for these.
      ! If there is no pointer, try to get use a pointer to one of these
      ! matrices directly.
      p_rmatrixTemplateGradient => rmatrixComponents%p_rmatrixTemplateGradient
      IF (.NOT. ASSOCIATED(p_rmatrixTemplateGradient)) &
        p_rmatrixTemplateGradient => rmatrixComponents%p_rmatrixB1
      IF (.NOT. ASSOCIATED(p_rmatrixTemplateGradient)) THEN
        PRINT *,'allocMatrix: Cannot set up B matrix in system matrix!'
        STOP
      END IF

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      IF (ASSOCIATED(p_rdiscretisation)) THEN
        CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)    
      ELSE
        ! No discretisation structure; create the matrix directly as 3x3 matrix.
        CALL lsysbl_createEmptyMatrix (rmatrix,NDIM2D+1)
      END IF
        
      ! Let's consider the global system in detail. The standard matrix It has 
      ! roughly the following shape:
      !
      !    ( A11       B1  ) = ( A11  A12  A13 )
      !    (      A22  B2  )   ( A21  A22  A23 )
      !    ( B1^T B2^T .   )   ( A31  A32  A33 )
      !
      ! All matices may have multiplication factors in their front.
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
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      IF (.NOT. bdecoupled .AND. .NOT. bfulltensor) THEN     
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix 
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content.
        CALL lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
      ELSE
      
        ! Otherwise, create another copy of the template matrix.
        CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                    
      END IF
      
      ! Manually change the discretisation structure of the Y-velocity 
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      CALL lsyssc_assignDiscretDirectMat (rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))
                                          
      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      IF (bfulltensor) THEN

        ! We have a matrix in the following shape:
        !
        !    ( A11  A12  B1  ) 
        !    ( A21  A22  B2  ) 
        !    ( B1^T B2^T .   ) 
        !
        ! Create A12 and A21.
      
        IF (rmatrix%RmatrixBlock(1,2)%cmatrixFormat &
            .EQ. LSYSSC_MATRIXUNDEFINED) THEN
            
          CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(1,2), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
            
        END IF

        IF (rmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .EQ. LSYSSC_MATRIXUNDEFINED) THEN
            
          ! Create a new matrix A21 in memory. create a new matrix
          ! using the template FEM matrix...
          CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(2,1), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
            
        END IF
        
      END IF

      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create empty space for the entries. 
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,3),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,3),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
      ! Furthermore, put B1^T and B2^T to the block matrix.
      ! These matrices will not change during the whole computation,
      ! so we can put refereces to the original ones to the system matrix.
      CALL lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(3,1),&
                                    LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(3,2),&
                                    LSYSSC_TR_VIRTUAL)
                                    
      ! That's it, all submatrices are basically set up.
      !
      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      CALL lsysbl_updateMatStrucInfo (rmatrix)
        
    END SUBROUTINE
    
    ! -----------------------------------------------------

    SUBROUTINE assembleVelocityBlocks (rmatrixComponents,rmatrix,rvector,dvectorWeight)
        
    ! Assembles the velocity matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    TYPE(t_ccmatrixComponents), INTENT(IN) :: rmatrixComponents
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
    
    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    TYPE(t_vectorBlock), OPTIONAL :: rvector
    
    ! Weight for the velocity vector; standard = -1.
    REAL(DP), INTENT(IN), OPTIONAL :: dvectorWeight
    
    ! local variables
    LOGICAL :: bshared
    INTEGER :: iupwind
    TYPE(t_convUpwind) :: rupwind
    TYPE(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    TYPE(T_jumpStabilisation) :: rjumpStabil
    REAL(DP) :: dvecWeight
    
      ! Standard value for dvectorWeight is = -1.
      dvecWeight = -1.0_DP
      IF (PRESENT(dvectorWeight)) dvecWeight = dvectorWeight
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      ! Allocate memory if necessary. Normally this should not be necessary...
      ! A11:
      IF (.NOT. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) THEN
        CALL lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
      END IF
    
      ! A22:
      IF (.NOT. bshared) THEN
        IF (.NOT. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) THEN
          CALL lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
        END IF
      END IF

      ! A12/ A21:
      IF (lsysbl_isSubmatrixPresent (rmatrix,1,2)) THEN
        IF (.NOT. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,2))) THEN
          CALL lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
        END IF
        IF (.NOT. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,1))) THEN
          CALL lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
        END IF
      END IF
    
      ! ---------------------------------------------------
      ! Plug in the mass matrix?
      IF (rmatrixComponents%dalpha .NE. 0.0_DP) THEN
       
        ! Allocate memory if necessary. Normally this should not be necessary...
        IF (.NOT. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) THEN
          CALL lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
        END IF
      
        CALL lsyssc_matrixLinearComb (&
            rmatrixComponents%p_rmatrixMass       ,rmatrixComponents%dalpha,&
            rmatrix%RmatrixBlock(1,1),0.0_DP,&
            rmatrix%RmatrixBlock(1,1),&
            .FALSE.,.FALSE.,.TRUE.,.TRUE.)
            
        IF (.NOT. bshared) THEN

          ! Allocate memory if necessary. Normally this should not be necessary...
          IF (.NOT. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) THEN
            CALL lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
          END IF

          CALL lsyssc_matrixLinearComb (&
              rmatrixComponents%p_rmatrixMass     ,rmatrixComponents%dalpha,&
              rmatrix%RmatrixBlock(2,2),0.0_DP,&
              rmatrix%RmatrixBlock(2,2),&
              .FALSE.,.FALSE.,.TRUE.,.TRUE.)
        END IF
        
      ELSE
      
        ! Otherwise, initialise the basic matrix with 0.0
        CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
        
        IF (.NOT. bshared) THEN
          CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
        END IF
          
      END IF
      
      ! ---------------------------------------------------
      ! Plug in the Stokes matrix?
      IF (rmatrixComponents%dtheta .NE. 0.0_DP) THEN
        CALL lsyssc_matrixLinearComb (&
            rmatrixComponents%p_rmatrixStokes     ,rmatrixComponents%dtheta,&
            rmatrix%RmatrixBlock(1,1),1.0_DP,&
            rmatrix%RmatrixBlock(1,1),&
            .FALSE.,.FALSE.,.TRUE.,.TRUE.)
            
        IF (.NOT. bshared) THEN
          CALL lsyssc_matrixLinearComb (&
              rmatrixComponents%p_rmatrixStokes   ,rmatrixComponents%dtheta,&
              rmatrix%RmatrixBlock(2,2),1.0_DP,&
              rmatrix%RmatrixBlock(2,2),&
              .FALSE.,.FALSE.,.TRUE.,.TRUE.)
        END IF
      END IF
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      IF (rmatrixComponents%dgamma .NE. 0.0_DP) THEN
      
        IF (.NOT. PRESENT(rvector)) THEN
          CALL output_line ('Velocity vector not present!', &
                             OU_CLASS_ERROR,OU_MODE_STD,'c2d2_assembleMatrix')
          STOP
        END IF
      
        SELECT CASE (rmatrixComponents%iupwind)
        CASE (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rmatrixComponents%dupsam
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rmatrixComponents%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rmatrixComponents%dnewton
          
          IF (rmatrixComponents%dnewton .EQ. 0.0_DP) THEN
          
            ! If the submatrices A12 and A21 exist, fill them with zero.
            ! If they don't exist, we don't have to do anything.
            IF (lsysbl_isSubmatrixPresent (rmatrix,1,2)) THEN
              CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
              CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
            END IF
            
         ELSE

            ! Clear A12/A21 that may receive parts of the Newton matrix
            CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
         END IF
         
          ! Call the SD method to calculate the nonlinearity.
          CALL conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)

        CASE (CCMASM_STAB_UPWIND)
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = rmatrixComponents%dupsam

          ! Matrix weight
          rupwind%dtheta = rmatrixComponents%dgamma
          
          IF (rmatrixComponents%dnewton .NE. 0.0_DP) THEN
            CALL output_line ('Warning: Upwind does not support assembly '&
                //'of the Newton matrix!',OU_CLASS_TRACE1)
          END IF
          
          ! Call the upwind method to calculate the nonlinear matrix.
          CALL conv_upwind2d (rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rupwind, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1)) 
                              
          IF (.NOT. bshared) THEN
            ! Modify also the matrix block (2,2)
            CALL conv_upwind2d (rvector, rvector, &
                                dvecWeight, 0.0_DP,&
                                rupwind, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2)) 
          END IF     

        CASE (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%dtheta = rmatrixComponents%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rmatrixComponents%dnewton
          
          IF (rmatrixComponents%dnewton .EQ. 0.0_DP) THEN

            ! If the submatrices A12 and A21 exist, fill them with zero.
            ! If they don't exist, we don't have to do anything.
            IF (lsysbl_isSubmatrixPresent (rmatrix,1,2)) THEN
              CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
              CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
            END IF
            
          ELSE

            ! Clear A12/A21 that receives parts of the Newton matrix
            CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they aren't.
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          END IF
         
          ! Call the SD method to calculate the nonlinearity.
          CALL conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)          
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rmatrixComponents%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rmatrixComponents%dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          CALL conv_jumpStabilisation2d (&
                              rvector, rvector, dvecWeight, 0.0_DP,&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))   

          IF (.NOT. bshared) THEN
            CALL conv_jumpStabilisation2d (&
                                rvector, rvector, dvecWeight, 0.0_DP,&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))   
          END IF

        CASE DEFAULT
          PRINT *,'Don''t know how to set up nonlinearity!?!'
          STOP
        
        END SELECT

      ELSE
      
        ! That's the Stokes-case. Jump stabilisation is possible...
      
        SELECT CASE (rmatrixComponents%iupwind)
        CASE (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rmatrixComponents%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rmatrixComponents%dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          CALL conv_jumpStabilisation2d (&
                              rvector, rvector, dvecWeight, 0.0_DP,&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))   

          IF (.NOT. bshared) THEN
            CALL conv_jumpStabilisation2d (&
                                rvector, rvector, dvecWeight, 0.0_DP,&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))   
          END IF

        CASE DEFAULT
          ! No stabilisation
        
        END SELECT
      
      END IF ! gamma <> 0
    
    END SUBROUTINE  
      
    ! -----------------------------------------------------
    
    SUBROUTINE assembleGradientMatrices (rmatrixComponents,rmatrix,bsharedMatrix)
    
    ! Initialises the gradient/divergence matrices with entries from
    ! the rmatrixComponents structure.
    !
    ! The routine copies references from the submatrices tormatrix,
    ! but it does not initialise any matrix weights / scaling factors.
    !
    ! If bsharedMatrix=TRUE, the matrix is created using references to the
    ! matrix building blocks in rlevelInfo, thus sharing all information
    ! with those matrices in rmatrixComponents. In this case, the caller must
    ! not change the matrix entries, because this would change the
    ! original 'template' matrices!
    ! (This can be used e.g. for setting up a matrix for building a defect
    !  vector without copying matrix data.)
    ! If bsharedMatrix=TRUE on the other hand, the matrix entries of the
    ! original template (B-) matrices are copied in memory,
    ! so the new matrix is allowed to be changed!

    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    TYPE(t_ccmatrixComponents), INTENT(IN) :: rmatrixComponents

    ! Block matrix where the B-matrices should be set up
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix

    ! Whether or not the matrix entries of the source gradient-matrices 
    ! should be copied in memory. 
    !
    ! If set to FALSE, the routine tries to initialise rmatrix
    ! only with references to the original matrices, thus the caller must not
    ! change the entries. Nevertheless, if rmatrix is the owner of one of the
    ! submatrices, the routine will always copy the matrix entries,
    ! as otherwise memory would have to be deallocated!
    !
    ! If set to TRUE, the entries of the source matrices in rmatrixComponents are
    ! copied, so the caller can change rmatrix afterwards (e.g. to implement
    ! boundary conditions).
    LOGICAL, INTENT(IN) :: bsharedMatrix

      ! local variables
      INTEGER :: idubStructure,idubContent
      
      ! Initialise a copy flag that tells the duplicateMatrix-routine whether to
      ! copy the entries or to create references.
      IF (bsharedMatrix) THEN
      
        idubContent = LSYSSC_DUP_SHARE
        
        ! Normally we share entries -- except for if the submatrices belong to
        ! rmatrix! To avoid memory deallocation in this case, we copy
        ! the entries.
        IF ((.NOT. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(1,3))) .OR.&
            (.NOT. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(2,3)))) THEN
          idubContent = LSYSSC_DUP_COPY
        END IF
        
      ELSE
      
        idubContent = LSYSSC_DUP_COPY
        
      END IF
      
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

      ! The B1/B2 matrices exist up to now only in rmatrixComponents.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The B-blocks
      ! are already prepared and memory for the entries is already allocated;
      ! so we only have to copy the entries.
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,3),&
                                    idubStructure,idubContent)

      CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,3),&
                                    idubStructure,idubContent)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      ! These matrices are always 'shared'.
      CALL lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(3,1),&
                                    LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(3,2),&
                                    LSYSSC_TR_VIRTUAL)

    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_assembleDefect (rmatrixComponents,rx,rd,cx,ry)

!<description>
  ! This routine assembles the nonlinear defect
  !      rd := rd - cx A(ry) rx
  ! with the system matrix A(.) defined by the configuration in rmatrixComponents.
  ! The caller must initialise the rmatrixComponents according to how the 
  ! matrix should look like.
  !
  ! The parameter ry is optional. If specified, this parameter defines where to
  ! evaluate the nonlinearity (if the system matrix $A$ contains a nonlinearity).
  ! If not specified, ry=rx is assumed.
  ! The multiplication factor cx is optional as well. If not specified,
  ! cx=1.0 is assumed.
  !
  ! The routine will not include any boundary conditions in the defect.
!</description>

  ! A t_ccmatrixComponents structure providing all necessary 'source' information
  ! about how to set up the matrix. 
  !
  ! The caller must provide either p_rmatrixTemplateXXXX in this structure
  ! or set the p_rmatrixTemplateXXXX as well as p_rdiscretisation to
  ! appropriate values. This is necessary for exploiting then structure
  ! of the matrix.
  TYPE(t_ccmatrixComponents), INTENT(IN) :: rmatrixComponents

  ! This vector specifies the 'x' that is multiplied to the matrix.
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rx

  ! OPTIONAL: Multiplication factor in front of the term 'A(ry) rx'.
  ! If not specified, cx=1.0 is assumed.
  REAL(DP), INTENT(IN), OPTIONAL :: cx

  ! OPTIONAL: Point where to evaluate the nonlinearity. If not specified,
  ! ry=rx is assumed.
  TYPE(t_vectorBlock), INTENT(IN), TARGET, OPTIONAL :: ry

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    REAL(DP) :: dcx
    TYPE(t_vectorBlock), POINTER :: p_ry
    TYPE(t_matrixBlock) :: rmatrix
    
    dcx = 1.0_DP
    IF (PRESENT(cx)) dcx = cx
    
    p_ry => rx
    IF (PRESENT(ry)) p_ry => ry
    
    ! The system matrix looks like:
    !          
    !    ( A11  A12  B1  ) 
    !    ( A21  A22  B2  ) 
    !    ( B1^T B2^T .   ) 
    ! 
    ! Create a temporary matrix that covers this structure.
    CALL lsysbl_createEmptyMatrix (rmatrix,NDIM2D+1)
    
    ! Put references to the Stokes- and B-matrices to Aij. assembleVelocityDefect 
    ! needs this template matrix to provide the structure for the stabilisation
    ! routines! The B-matrices are needed later.
    CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
    IF (rmatrixComponents%dnewton .NE. 0.0_DP) THEN
      CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    END IF
    
    CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1,&
        rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2,&
        rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    CALL lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                  rmatrix%RmatrixBlock(3,1),&
                                  LSYSSC_TR_VIRTUAL)

    CALL lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                  rmatrix%RmatrixBlock(3,2),&
                                  LSYSSC_TR_VIRTUAL)

    ! Update the structural information of the block matrix, as we manually
    ! changed the submatrices:
    CALL lsysbl_updateMatStrucInfo (rmatrix)
    
    ! In the first step, we assemble the defect that arises in the velocity 
    ! components. This is characterised by the following submatrix:
    !
    !    ( A11  A12  .  ) 
    !    ( A21  A22  .  ) 
    !    ( .    .    .  ) 
    !
    ! assembleVelocityDefect handles exactly these submatices.

    CALL assembleVelocityDefect (rmatrixComponents,rmatrix,rx,rd,p_ry,dcx)
    
    ! Now, we treat all the remaining blocks. Let's see what is missing:
    !
    !    ( .    .    B1  ) 
    !    ( .    .    B2  ) 
    !    ( B1^T B2^T .   ) 

    ! To build the appropriate defect, we firat remove the velocity blocks:
    
    CALL lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,1))
    CALL lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,2))
    CALL lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,1))
    CALL lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,2))

    ! Initialise the weights for the B/B^T matrices
    rmatrix%RmatrixBlock(1,3)%dscaleFactor = rmatrixComponents%deta
    rmatrix%RmatrixBlock(2,3)%dscaleFactor = rmatrixComponents%deta
    
    rmatrix%RmatrixBlock(3,1)%dscaleFactor = rmatrixComponents%dtau
    rmatrix%RmatrixBlock(3,2)%dscaleFactor = rmatrixComponents%dtau

    ! ------------------------------------------------
    ! Build the defect by matrix-vector multiplication
    !
    ! Note that no time step or whatever is included here; everything
    ! is initialised with the multiplication factors in the submatrices
    ! from above!
    CALL lsysbl_blockMatVec (rmatrix, rx, rd, -dcx, 1.0_DP)
    
    ! Release the temporary matrix, we don't need it anymore.
    CALL lsysbl_releaseMatrix (rmatrix)

  CONTAINS

    SUBROUTINE assembleVelocityDefect (rmatrixComponents,&
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
    
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    TYPE(t_ccmatrixComponents), INTENT(IN) :: rmatrixComponents

    ! Reference to the system matrix. Only the structure of the matrix
    ! is used to reconstruct the structure of the discretisation.
    ! The content of the matrix is not changed or used.
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
    
    ! Solution vector.
    TYPE(t_vectorBlock), INTENT(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    TYPE(t_vectorBlock), INTENT(INOUT) :: rdefect
    
    ! Weight for the velocity vector; usually = 1.0
    REAL(DP), INTENT(IN) :: dvectorWeight
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. The first two blocks in that block vector are
    ! used as velocity field.
    TYPE(t_vectorBlock), INTENT(IN) :: rvelocityVector

    ! local variables
    LOGICAL :: bshared
    INTEGER :: iupwind
    TYPE(t_convUpwind) :: rupwind
    TYPE(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    TYPE(T_jumpStabilisation) :: rjumpStabil
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))

      ! ---------------------------------------------------
      ! Subtract the mass matrix stuff?
      IF (rmatrixComponents%dalpha .NE. 0.0_DP) THEN
        CALL lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, &
            rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rmatrixComponents%dalpha, 1.0_DP)

        CALL lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, &
            rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
            -rmatrixComponents%dalpha, 1.0_DP)
      END IF
      
      ! ---------------------------------------------------
      ! Subtract the Stokes matrix stuff?
      IF (rmatrixComponents%dtheta .NE. 0.0_DP) THEN
        CALL lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixStokes, &
            rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rmatrixComponents%dtheta, 1.0_DP)

        CALL lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixStokes, &
            rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
            -rmatrixComponents%dtheta, 1.0_DP)
      END IF
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      IF (rmatrixComponents%dgamma .NE. 0.0_DP) THEN
      
        ! Type of stablilisation?
        iupwind = rmatrixComponents%iupwind
      
        SELECT CASE (iupwind)
        CASE (0)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rmatrixComponents%dupsam
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rmatrixComponents%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rmatrixComponents%dnewton
          
          ! Call the SD method to calculate the defect of the nonlinearity.
          ! As rrhsTemp shares its entries with rdefect, the result is
          ! directly written to rdefect!
          ! As velocity field, we specify rvelocityVector here. The first two
          ! subvectors are used as velocity field.
          
          CALL conv_streamlineDiffusionBlk2d (&
                              rvelocityVector, rvelocityVector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)
                              
        CASE (1)
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = rmatrixComponents%dupsam

          ! Matrix weight
          rupwind%dtheta = rmatrixComponents%dgamma
          
          ! Call the upwind method to calculate the nonlinear defect.
          CALL conv_upwind2d (rvector, rvector, &
                              dvectorWeight, 0.0_DP,&
                              rupwind, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rvector,rdefect) 
                              
          IF (.NOT. bshared) THEN
            PRINT *,'Upwind does not support independent A11/A22!'
            STOP
          END IF     

        CASE (2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%dtheta = rmatrixComponents%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rmatrixComponents%dnewton
          
          IF (rmatrixComponents%dnewton .EQ. 0.0_DP) THEN

            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there's content at all),
            ! these matrices are ignored then by the kernel.
            
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
            
          ELSE

            ! Clear A12/A21 that receives parts of the Newton matrix
            CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            CALL lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they aren't.
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          END IF
         
          ! Call the SD method to calculate the nonlinearity.
          CALL conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)          
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rmatrixComponents%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rmatrixComponents%dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          CALL conv_jumpStabilisation2d (&
                              rvector, rvector, dvectorWeight, 0.0_DP,&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),&
                              rsolution=rvector,rdefect=rdefect)   

          IF (.NOT. bshared) THEN
            PRINT *,'Edge oriented stabilisation does not support independent A11/A22!'
            STOP
          END IF

        CASE DEFAULT
          PRINT *,'Don''t know how to set up nonlinearity!?!'
          STOP
        
        END SELECT
      
      ELSE
      
        ! That's the Stokes-case. Jump stabilisation is possible...
        !
        ! Type of stablilisation?
        iupwind = rmatrixComponents%iupwind
      
        SELECT CASE (iupwind)
        CASE (2)
          ! Jump stabilisation.
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = rmatrixComponents%dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = rmatrixComponents%dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          CALL conv_jumpStabilisation2d (&
                              rvector, rvector, dvectorWeight, 0.0_DP,&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),&
                              rsolution=rvector,rdefect=rdefect)   

          IF (.NOT. bshared) THEN
            CALL conv_jumpStabilisation2d (&
                                rvector, rvector, dvectorWeight, 0.0_DP,&
                                rjumpStabil, CONV_MODDEFECT, &
                                rmatrix%RmatrixBlock(2,2),&
                                rsolution=rvector,rdefect=rdefect)   
          END IF

        CASE DEFAULT
          ! No stabilisation
        
        END SELECT
      
      END IF ! gamma <> 0
    
    END SUBROUTINE

  END SUBROUTINE

END MODULE
