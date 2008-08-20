!##############################################################################
!# ****************************************************************************
!# <name> vanca </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the implementations of the VANCA preconditioner.
!# These are more or less auxiliary routines called by the VANCA 
!# preconditioner in the linearsolver.f90 solver library.
!#
!# The following routines can be found here:
!#
!#  1.) vanca_initConformal
!#      -> Initialise the VANCA for conformal discretisations.
!# 
!#  2.) vanca_conformal
!#      -> Perform one step of the VANCA solver conformal discretisations. Calls a
!#         specialised sub-VANCA-variant to do the work.
!#
!#  3.) vanca_doneConformal
!#      -> Clean up the VANCA for conformal discretisations.
!#
!# The following list of routines are all used internally. There is no need to
!# call them directly.
!# 
!#  4.) vanca_initGeneralVanca
!#      -> Initialise the general Vanca solver
!#
!#  5.) vanca_general
!#      -> Perform one step of the full VANCA solver for general block systems.
!#
!#  6.) vanca_doneGeneralVanca
!#      -> Clean up the general VANCA solver
!#
!#  7.) vanca_init2DSPQ1TQ0simple
!#      -> Initialise specialised VANCA solver for 2D saddle point problems
!#         with $\tilde Q_1/Q_0$ discretisation.
!#         Deprecated, ist not used.
!#
!#  8.) vanca_init2DNavierStokes 
!#      -> Initialise the VANCA solver for the problem class 
!#         '2D Navier-Stokes equation'. 
!#
!#  9.) vanca_2DNavierStokes 
!#      -> Apply the VANCA solver for the problem class 
!#         '2D Navier-Stokes equation'. 
!#
!# 10.) vanca_2DSPQ1TQ0simple
!#      -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#         problems with $\tilde Q_1/Q_0$ discretisation
!#
!# 11.) vanca_2DSPQ1TQ0simpleConf
!#      -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#         problems with $\tilde Q_1/Q_0$ discretisation.
!#         Applies VANCA only to a subset of all elements in the domain.
!#
!# 12.) vanca_2DSPQ1TQ0simpleCoupConf
!#      -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#         problems with $\tilde Q_1/Q_0$ discretisation.
!#         Applies VANCA only to a subset of all elements in the domain.
!#         This variant can handle fully coupled matrices.
!#
!# 13.) vanca_2DSPQ1TQ0fullConf
!#      -> Perform one step of the specialised 'full' VANCA solver for 2D saddle 
!#         point problems with $\tilde Q_1/Q_0$ discretisation.
!#         Applies VANCA only to a subset of all elements in the domain.
!#
!# 14.) vanca_2DSPQ1TQ0fullCoupConf
!#      -> Perform one step of the specialised 'full' VANCA solver for 2D saddle 
!#         point problems with $\tilde Q_1/Q_0$ discretisation.
!#         Applies VANCA only to a subset of all elements in the domain.
!#         This variant can handle fully coupled matrices.
!#
!# 15.) vanca_2DSPQ2QP1simple
!#      -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#         problems with $Q_2/QP_1$ discretisation. Diagonal VANCA approach.
!#
!# 16.) vanca_2DSPQ2QP1full
!#      -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#         problems with $Q_2/QP_1$ discretisation. Full VANCA approach.
!#
!# 17.) vanca_2DSPQ2QP1simpleConf
!#      -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#         problems with $Q_2/QP_1$ discretisation. Diagonal VANCA approach.
!#         Applies VANCA only to a subset of all elements in the domain.
!#
!# 18.) vanca_2DSPQ2QP1fullConf
!#      -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#         problems with $Q_2/QP_1$ discretisation. Full VANCA approach.
!#         Applies VANCA only to a subset of all elements in the domain.
!#
!# 19.) vanca_init2DNavierStokesOptC
!#      -> Initialise the VANCA solver for 2D Navier-Stokes optimal control
!#         problems. Specialised $\tilde Q1/Q0$ version, full VANCA approach.
!#
!# 20.) vanca_2DNSSOCQ1TQ0fullCoupConf
!#      -> Perform one step of the VANCA solver for 2D Navier-Stokes optimal
!#         control problems. Specialised $\tilde Q1/Q0$ version, full VANCA approach.
!#
!# 21.) vanca_init3DNavierStokes 
!#      -> Initialise the VANCA solver for the problem class 
!#         '3D Navier-Stokes equation'. 
!#
!# 22.) vanca_3DNavierStokes 
!#      -> Apply the VANCA solver for the problem class 
!#         '3D Navier-Stokes equation'. 
!#
!# 23.) vanca_3DSPQ1TQ0simple
!#      -> Perform one step of the specialised VANCA solver for 3D saddle point 
!#         problems with $\tilde Q_1/Q_0$ discretisation
!#
!# 24.) vanca_3DSPQ1TQ0simpleConf
!#      -> Perform one step of the specialised VANCA solver for 3D saddle point 
!#         problems with $\tilde Q_1/Q_0$ discretisation.
!#         Applies VANCA only to a subset of all elements in the domain.
!#
!# 25.) vanca_3DSPQ1TQ0fullConf
!#      -> Perform one step of the specialised 'full' VANCA solver for 3D saddle 
!#         point problems with $\tilde Q_1/Q_0$ discretisation.
!#         Applies VANCA only to a subset of all elements in the domain.
!#
!#  History
!# ---------
!# Originally, the following VANCA variants were implemented:
!#  - vanca_general
!#  - vanca_init2DSPQ1TQ0simple
!# vanca_general was designed to work with everything. vanca_init2DSPQ1TQ0simple was
!# designed to work with the 2D Navier Stokes problem, uniform $\tilde Q_1/Q_0$
!# discretisation. It used the 'diagonal' VANCA approach.
!#
!# Over the time, there was a 'full' VANCA variant added as well as variants that
!# were designed to work with $Q_2/QP_1$ discretisations. However, each variant
!# had to be chosen carefully, as all variants were rather specialised to a
!# special uniform discretisation.
!#
!# To work with more general configurations, a VANCA wrapper was introduced.
!# This replaces the needs for carefully choosing the correct VANCA variant
!# and can be seen as 'more general VANCA', as it's also designed to work
!# with conformal discretisations.
!#
!# For this purpose, all VANCA variants were modified to work with element
!# lists, e.g.
!#   vanca_2DSPQ1TQ0simple -> vanca_2DSPQ1TQ0simpleConf
!#   vanca_2DSPQ2QP1simple -> vanca_2DSPQ2QP1simpleConf
!#   ...
!# For the Navier-Stokes problem, the corresponding wrapper function is
!#  - vanca_initConformal / vanca_conformal
!# These will check the discretisation and call the correct
!# vanca_xxxxConf subroutines automatically.
!#
!# To support other types of problems, the user should proceed as follows:
!#  - Write a VANCA variant in the '....Conf' style that accepts an element
!#    list, which specifies a subset of the domain where to apply VANCA
!#    (like vanca_2DSPQ1TQ0simpleConf).
!#  - Add a VANCA subtype constant and a VANCA problem class constant
!#    to the VANCATP_xxxx and VANCACP_xxxx list if necessary.
!#  - Introduce new VANCA structures as necessary and add references
!#    to that to the t_vanca structure.
!#  - Modify the vanca_initConformal / vanca_conformal routines so that
!#    they initialise and call the new VANCA variant correctly.
!# </purpose>
!##############################################################################

MODULE vanca

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE genoutput

  IMPLICIT NONE

!<constants>

!<constantblock description="Identifiers for the different problem classes that can be handled by VANCA">

  ! General VANCA
  INTEGER, PARAMETER :: VANCAPC_GENERAL        = 0

  ! 2D Navier Stokes problem
  INTEGER, PARAMETER :: VANCAPC_2DNAVIERSTOKES = 1

  ! 2D Navier Stokes optimal control problem
  INTEGER, PARAMETER :: VANCAPC_2DNAVIERSTOKESOPTC = 2

  ! 3D Navier Stokes problem
  INTEGER, PARAMETER :: VANCAPC_3DNAVIERSTOKES = 3

!</constantblock>


!<constantblock description="Identifiers for the different VANCA subtypes. Which one is supported depends on the problem class.">

  ! Standard VANCA, most suitable for the corresponding situation
  INTEGER, PARAMETER :: VANCATP_STANDARD  = 0

  ! Diagonal-type VANCA
  INTEGER, PARAMETER :: VANCATP_DIAGONAL  = 0

  ! 'Full' VANCA
  INTEGER, PARAMETER :: VANCATP_FULL      = 1

  ! 'Full' VANCA for optimal control problems, primal equation processing
  INTEGER, PARAMETER :: VANCATP_FULLOPTC_PRIMAL = 2

  ! 'Full' VANCA for optimal control problems, dual equation processing
  INTEGER, PARAMETER :: VANCATP_FULLOPTC_DUAL   = 3

  ! 'Diagonal' VANCA for optimal control problems
  INTEGER, PARAMETER :: VANCATP_DIAGOPTC        = 4

  ! Diagonal-type VANCA, solution based
  INTEGER, PARAMETER :: VANCATP_DIAGONAL_SOLBASED = 5

!</constantblock>

!</constants>


!<types>
  
!<typeblock>
  
  ! A structure that accepts a pointer to the column/row/data arrays
  ! of a structure-7/structure 9 matrix. This is usually passed to
  ! the VANCA preconditioner(s) to specify the matrices to handle.
  
  
  TYPE t_matrixPointer79Vanca
    ! Is set to FALSE if the matrix does not exist/is empty.
    ! In this case, the pointers below are undefined!
    LOGICAL :: bexists
    
    ! TRUE if the matrix is saved transposed
    LOGICAL :: btransposed
    
    ! The scaling factor of the matrix; from the matrix structure.
    REAL(DP) :: dscaleFactor
  
    ! Pointer to the data - currently only double precision 
    REAL(DP), DIMENSION(:), POINTER :: p_DA
    
    ! Pointer to the column structure
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    
    ! Pointer to the row structure
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A structure for saving precalculated information for the general VANCA.
  ! This is initialised by vanca_initGeneralVanca and released by
  ! vanca_doneGeneralVanca.
  
  TYPE t_vancaGeneral
  
    ! Number of blocks in the global matrix
    INTEGER                                              :: nblocks
    
    ! Pointer to the block matrix
    TYPE(t_matrixBlock), POINTER                         :: p_rmatrix
  
    ! Pointers to t_matrixPointer79Vanca structures specifying
    ! the submatrices and their properties.
    TYPE(t_matrixPointer79Vanca), DIMENSION(:,:),POINTER :: p_Rmatrices
    
    ! Maximum number of local DOF's.
    INTEGER(PREC_DOFIDX)                              :: nmaxLocalDOFs
    
    ! Total number of local DOF's
    INTEGER                                           :: ndofsPerElement

    ! Number of local DOF's in the element distributions of all blocks.
    ! Note that this VANCA supports only uniform discretisations, so
    ! each entry corresponds to one block in the solution vector.
    INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_InDofsLocal => NULL()
    
    ! Offset indices of the blocks in the solution vector. IblockOffset(i)
    ! points to the beginning of the i'th block of the solution vector.
    INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IblockOffset => NULL()
    
    ! Temporary array that saves the DOF's that are in processing when
    ! looping over an element set.
    ! DIMENSION(nmaxLocalDOFs,VANCA_NELEMSIM,nblocks)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), POINTER :: p_IelementDOFs => NULL()

  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D-Navier-Stokes
  ! VANCA method for Navier-Stokes systems.
  
  TYPE t_vancaPointer2DNavSt
    ! Pointer to the column structure of the velocity matrix A11 and A22
    ! A11 and A22 must have the same structure.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA => NULL()
    
    ! Pointer to the row structure of the velocity matrix A11 and A22.
    ! They must have the same structure.
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA => NULL()
    
    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA => NULL()

    ! Pointer to the matrix entries of the velocity matrix A11
    REAL(DP), DIMENSION(:), POINTER             :: p_DA => NULL()

    ! Pointer to the matrix entries of the velocity matrix A22 or NULL
    ! A11=A22.
    REAL(DP), DIMENSION(:), POINTER             :: p_DA22 => NULL()

    ! Pointer to the column structure of the velocity matrix A11 and A22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12 => NULL()
    
    ! Pointer to the row structure of the velocity matrix A11 and A22
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12 => NULL()
    
    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA12 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A12 or NULL
    ! if not present
    REAL(DP), DIMENSION(:), POINTER             :: p_DA12 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A21 or NULL
    ! if not present
    REAL(DP), DIMENSION(:), POINTER             :: p_DA21 => NULL()

    ! Pointer to the column structure of the B/D-matrices.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB => NULL()
    
    ! Pointer to the row structure of the B/D-matrices
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB => NULL()
    
    ! Pointer to the entries of the B1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1 => NULL()

    ! Pointer to the entries of the B2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2 => NULL()
    
    ! Pointer to the entries of the D1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1 => NULL()

    ! Pointer to the entries of the D2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2 => NULL()

    ! Pointer to the matrix entries of the pressure identity matrix A33
    ! (if it exists).
    REAL(DP), DIMENSION(:), POINTER             :: p_DA33 => NULL()

    ! Pointer to diagonal entries of A33
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA33 => NULL()

    ! Spatial discretisation structure for X-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrU => NULL()
    
    ! Spatial discretisation structure for Y-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrV => NULL()
    
    ! Spatial discretisation structure for pressure
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrP => NULL()
    
    ! Multiplication factors for the submatrices; taken from the system matrix.
    ! (-> Not used in the current implementation! Although it's easy to include
    ! that into VANCA, some further speed analysis has to be done to make
    ! sure there's not too much speed impact when using these!)
    REAL(DP), DIMENSION(3,3) :: Dmultipliers
    
    ! A temporary vector for solution based VANCA variants.
    ! Undefined if not needed.
    TYPE(t_vectorBlock) :: rtempVector
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D-Navier-Stokes
  ! VANCA method for Navier-Stokes optimal control systems.
  
  TYPE t_vancaPointer2DNavStOptC
    ! Pointer to the column structure of the velocity matrix A11 and A22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA11 => NULL()
    
    ! Pointer to the row structure of the velocity matrix A11 and A22
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA11 => NULL()
    
    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA11 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A11
    REAL(DP), DIMENSION(:), POINTER             :: p_DA11 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A22 or NULL
    ! if not present
    REAL(DP), DIMENSION(:), POINTER             :: p_DA22 => NULL()

    ! Pointer to the column structure of the velocity matrix A11 and A22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12 => NULL()
    
    ! Pointer to the row structure of the velocity matrix A11 and A22
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12 => NULL()
    
    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA12 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A12 or NULL
    ! if not present
    REAL(DP), DIMENSION(:), POINTER             :: p_DA12 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A21 or NULL
    ! if not present
    REAL(DP), DIMENSION(:), POINTER             :: p_DA21 => NULL()


    ! Pointer to the matrix entries of the velocity matrix A44
    REAL(DP), DIMENSION(:), POINTER             :: p_DA44 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A55
    REAL(DP), DIMENSION(:), POINTER             :: p_DA55 => NULL()


    ! Pointer to the matrix entries of the pressure identity matrix A33
    ! (if it exists).
    REAL(DP), DIMENSION(:), POINTER             :: p_DA33 => NULL()

    ! Pointer to diagonal entries of A33
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA33 => NULL()

    ! Pointer to the matrix entries of the pressure identity matrix A66
    ! (if it exists).
    REAL(DP), DIMENSION(:), POINTER             :: p_DA66 => NULL()

    ! Pointer to diagonal entries of A66
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA66 => NULL()


    ! Pointer to the column structure of the matrix A45 and A54
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA45 => NULL()
    
    ! Pointer to the row structure of the matrix A45 and A54
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA45 => NULL()
    
    ! Pointer to diagonal entries in the matrix A45 and A54
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA45 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A45 or NULL
    ! if not present
    REAL(DP), DIMENSION(:), POINTER             :: p_DA45 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A54 or NULL
    ! if not present
    REAL(DP), DIMENSION(:), POINTER             :: p_DA54 => NULL()


    ! Pointer to the column structure of the mass matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolM => NULL()
    
    ! Pointer to the row structure of the mass matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldM => NULL()
    
    ! Pointer to diagonal entries in the mass matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalM => NULL()

    ! Pointer to the matrix entries of the mass matrix at position
    ! (1,4) and (2,5) in the pimal system, or NULL if not present.
    ! Has the same structure as the mass matrix.
    REAL(DP), DIMENSION(:), POINTER             :: p_DM14 => NULL()

    ! Pointer to the coupling system at position (4,1), or NULL if not present
    ! Has the same structure as the mass matrix.
    REAL(DP), DIMENSION(:), POINTER             :: p_DR41 => NULL()

    ! Pointer to the coupling system at position (5,2), or NULL if not present
    ! Has the same structure as the mass matrix.
    REAL(DP), DIMENSION(:), POINTER             :: p_DR52 => NULL()

    ! Pointer to the coupling system at position (4,2), or NULL if not present
    ! Has the same structure as the mass matrix.
    REAL(DP), DIMENSION(:), POINTER             :: p_DR42 => NULL()

    ! Pointer to the coupling system at position (5,1), or NULL if not present
    ! Has the same structure as the mass matrix.
    REAL(DP), DIMENSION(:), POINTER             :: p_DR51 => NULL()


    ! Pointer to the column structure of the B/D-matrices.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB => NULL()
    
    ! Pointer to the row structure of the B/D-matrices
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB => NULL()
    
    ! Pointer to the entries of the B1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1 => NULL()

    ! Pointer to the entries of the B2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2 => NULL()
    
    ! Pointer to the entries of the D1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1 => NULL()

    ! Pointer to the entries of the D2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2 => NULL()
    

    ! Spatial discretisation structure for X-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrU => NULL()
    
    ! Spatial discretisation structure for Y-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrV => NULL()
    
    ! Spatial discretisation structure for pressure
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrP => NULL()
    
    ! Multiplication factors for the submatrices; taken from the system matrix.
    ! (-> Not used in the current implementation! Although it's easy to include
    ! that into VANCA, some further speed analysis has to be done to make
    ! sure there's not too much speed impact when using these!)
    REAL(DP), DIMENSION(6,6) :: Dmultipliers
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 3D-Navier-Stokes
  ! VANCA method for Navier-Stokes systems.
  
  TYPE t_vancaPointer3DNavSt
    ! Pointer to the column structure of the velocity matrix A11, A22 and A33
    ! A11, A22 and A33 must have the same structure.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA => NULL()
    
    ! Pointer to the row structure of the velocity matrix A11, A22 and A33.
    ! They must have the same structure.
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA => NULL()
    
    ! Pointer to diagonal entries in the velocity matrix A11, A22 and A33
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA => NULL()

    ! Pointer to the matrix entries of the velocity matrix A11
    REAL(DP), DIMENSION(:), POINTER             :: p_DA => NULL()

    ! Pointer to the matrix entries of the velocity matrix A22 or NULL
    ! if A11=A22.
    REAL(DP), DIMENSION(:), POINTER             :: p_DA22 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A33 or NULL
    ! if A11=A33.
    REAL(DP), DIMENSION(:), POINTER             :: p_DA33 => NULL()

!    ! Pointer to the column structure of the velocity matrix A12 and A21
!    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12 => NULL()
!    
!    ! Pointer to the row structure of the velocity matrix A12 and A21
!    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12 => NULL()
!    
!    ! Pointer to diagonal entries in the velocity matrix A12 and A21
!    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA12 => NULL()
!
!    ! Pointer to the matrix entries of the velocity matrix A12 or NULL
!    ! if not present
!    REAL(DP), DIMENSION(:), POINTER             :: p_DA12 => NULL()
!
!    ! Pointer to the matrix entries of the velocity matrix A21 or NULL
!    ! if not present
!    REAL(DP), DIMENSION(:), POINTER             :: p_DA21 => NULL()

    ! Pointer to the column structure of the B/D-matrices.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB => NULL()
    
    ! Pointer to the row structure of the B/D-matrices
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB => NULL()
    
    ! Pointer to the entries of the B1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1 => NULL()

    ! Pointer to the entries of the B2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2 => NULL()
    
    ! Pointer to the entries of the B3-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB3 => NULL()

    ! Pointer to the entries of the D1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1 => NULL()

    ! Pointer to the entries of the D2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2 => NULL()

    ! Pointer to the entries of the D3-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD3 => NULL()

    ! Pointer to the matrix entries of the pressure identity matrix A44
    ! (if it exists).
    REAL(DP), DIMENSION(:), POINTER             :: p_DA44 => NULL()

    ! Pointer to diagonal entries of A33
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA44 => NULL()

    ! Spatial discretisation structure for X-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrU => NULL()
    
    ! Spatial discretisation structure for Y-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrV => NULL()
    
    ! Spatial discretisation structure for Z-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrW => NULL()

    ! Spatial discretisation structure for pressure
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrP => NULL()
    
    ! Multiplication factors for the submatrices; taken from the system matrix.
    ! (-> Not used in the current implementation! Although it's easy to include
    ! that into VANCA, some further speed analysis has to be done to make
    ! sure there's not too much speed impact when using these!)
    REAL(DP), DIMENSION(4,4) :: Dmultipliers
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A structure that represents a general VANCA configuration.
  
  TYPE t_vanca
  
    ! The problem class, the VANCA should be applient to. One of the
    ! VANCAPC_xxxx constants, e.g. VANCAPC_2DNAVIERSTOKES.
    ! This constants decides where in the t_vancaXXXX structures below
    ! the actual data about the VANCA can be found.
    INTEGER :: cproblemClass
  
    ! The subtype of VANCA that should handle the above problem class.
    ! One of the VANCATP_xxxx constants, e.g. VANCATP_DIAGONAL.
    INTEGER :: csubtype
    
    ! This flag is set in the init-routine of VANCA and indicates whether
    ! an extended version of VANCA must be called that supports scaled
    ! matrices or a diagonal matrix in the pressure block.
    ! =0: standard type. =1: extended version
    INTEGER :: csubsubtype = 0
    
    ! Configuration block with parameters for the general VANCA;
    ! only vaid if if cproblemClassVanca==VANCAPC_GENERAL.
    TYPE(t_vancaGeneral) :: rvancaGeneral
    
    ! Configuration block with parameters for the 2D Navier-Stokes VANCA;
    ! only vaid if if cproblemClassVanca==VANCAPC_2DNAVIERSTOKES.
    TYPE(t_vancaPointer2DNavSt) :: rvanca2DNavSt
    
    ! Configuration block with parameters for the 2D Navier-Stokes VANCA
    ! for optimal control problems;
    ! only vaid if if cproblemClassVanca==VANCAPC_2DNAVIERSTOKESOPTC.
    TYPE(t_vancaPointer2DNavStOptC) :: rvanca2DNavStOptC
    
    ! Configuration block with parameters for the 3D Navier-Stokes VANCA;
    ! only vaid if if cproblemClassVanca==VANCAPC_3DNAVIERSTOKES.
    TYPE(t_vancaPointer3DNavSt) :: rvanca3DNavSt

  END TYPE
  
!</typeblock>

!</types>

!<constants>
!<constantblock description="Constants defining the blocking of element sets in VANCA">

  ! Number of elements to handle simultaneously in general VANCA
  INTEGER :: VANCA_NELEMSIM   = 1000
  
!</constantblock>
!</constants>

CONTAINS

  ! ***************************************************************************
  ! General VANCA for conformal discretisations.
  !
  ! The general VANCA is configured in vanca_initConformal to a special type
  ! problem (2D Navier Stokes e.g.) and applies the most suitable VANCA
  ! algorithm to a vector when being executed. Roughtly said, it's a wrapper
  ! for all types of VANCA that are realised in this module.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_initConformal (rmatrix,rvanca,cproblemClass,csubtype)
  
!<description>
  ! Initialises the VANCA for conformal discretisations.
  ! Checks if the VANCA variant for conformal discretisations
  ! can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! VANCA does in general not support scaled matrices in rmatrix. The only
  ! exception is that by setting dscaleFactor=0.0 in one of the submatrices,
  ! a matrix can be deactivated.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: rmatrix
  
  ! The problem class that should be handled with this VANCA. One of the
  ! VANCAPC_xxxx constants, e.g. VANCAPC_2DNAVIERSTOKES
  INTEGER, INTENT(IN) :: cproblemClass
  
  ! The VANCA solver subtype that should handle the above problem class.
  ! One of the VANCATP_xxxx constants, e.g. VANCATP_DIAGONAL.
  INTEGER, INTENT(IN) :: csubtype
!</input>

!<output>
  ! t_vanca structure that saves algorithm-specific parameters.
  TYPE(t_vanca), INTENT(OUT) :: rvanca
!</output>

!</subroutine>

    SELECT CASE (cproblemClass)
    CASE (VANCAPC_GENERAL)
      ! General VANCA
      CALL vanca_initGeneralVanca (rmatrix,rvanca%rvancaGeneral)    
    
    CASE (VANCAPC_2DNAVIERSTOKES)
      ! Vanca for 2D Navier-Stokes problems
      CALL vanca_init2DNavierStokes (rmatrix,rvanca,csubtype)

    CASE (VANCAPC_2DNAVIERSTOKESOPTC)
      ! Vanca for 2D Navier-Stokes optimal control problems
      CALL vanca_init2DNavierStokesOptC (rmatrix,rvanca)

    CASE (VANCAPC_3DNAVIERSTOKES)
      ! Vanca for 3D Navier-Stokes problems
      CALL vanca_init3DNavierStokes (rmatrix,rvanca)

    END SELECT
    
    ! Initialise general data in the VANCA structure.
    rvanca%cproblemClass = cproblemClass
    rvanca%csubtype = csubtype

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_doneConformal (rvanca)

!<description>
  ! This routine cleans up a general VANCA solver. All memory allocated in
  ! rvancaGeneral is released.
!</description>
  
!<inputoutput>
  ! The VANCA structure to be cleaned up.
  TYPE(t_vanca), INTENT(INOUT)       :: rvanca
!</inputoutput>

!</subroutine>

    ! Type of VANCA?
    SELECT CASE (rvanca%cproblemClass)
    CASE (VANCAPC_2DNAVIERSTOKES)
      ! Vanca for 2D Navier-Stokes problems
      CALL vanca_done2DNavierStokes (rvanca%rvanca2DNavSt)
    CASE (VANCAPC_GENERAL)
      ! Release data of the general VANCA
      CALL vanca_doneGeneralVanca (rvanca%rvancaGeneral)
      
    ! ELSE: nothing to do.
    END SELECT
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_conformal (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the VANCA algorithm to the system $Ax=b$.
  ! This VANCA variant is rather general and can be applied to arbitrary
  ! conformal discretisations. It automatically chooses the correct
  ! VANCA 'subvariant' depending on the discretisation.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANCA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! t_vanca structure that saves algorithm-specific parameters.
  TYPE(t_vanca), INTENT(INOUT) :: rvanca

  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(INOUT)      :: rvector
!</inputoutput>

!</subroutine>

    ! Which type of VANCA should be used?
    SELECT CASE (rvanca%cproblemClass)
    CASE (VANCAPC_GENERAL)
      ! Use the general VANCA. This one handles the element distributions
      ! internally, so afterwards we can jump out of the loop without handling
      ! the other element distributions.
      
      CALL vanca_general (rvanca%rvancaGeneral, rvector, rrhs, domega)
      
    CASE (VANCAPC_2DNAVIERSTOKES)
      ! 2D Navier Stokes problem.
      CALL vanca_2DNavierStokes (rvanca%rvanca2DNavSt, rvector, rrhs, domega,&
          rvanca%csubtype,rvanca%csubsubtype)

    CASE (VANCAPC_2DNAVIERSTOKESOPTC)
      ! 2D Navier Stokes problem.
      CALL vanca_2DNavierStokesOptC (rvanca%rvanca2DNavStOptC, rvector, rrhs, domega,&
          rvanca%csubtype)

    CASE (VANCAPC_3DNAVIERSTOKES)
      ! 3D Navier Stokes problem.
      CALL vanca_3DNavierStokes (rvanca%rvanca3DNavSt, rvector, rrhs, domega,&
          rvanca%csubtype,rvanca%csubsubtype)

    CASE DEFAULT
      CALL output_line ('Unknown VANCA problem class!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_conformal')
      CALL sys_halt()  
      
    END SELECT
    
  END SUBROUTINE

  ! ***************************************************************************
  ! GENERAL VANCA! Supports (non-transposed) matrices of any kind and size.
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_initGeneralVanca (rmatrix,rvancaGeneral)

!<description>
  ! This routine initialises the general VANCA solver and allocates
  ! necessary memory for the iteration.
!</description>

!<input>
  ! The system matrix that is used during the VANCA iteration.
  ! Remark: VANCA saves a pointer to this matrix, so the structure must exist
  !  until the system is solved! (Usually this points to the system matrix in
  !  the corresponding solver structure of the underlying linear solver...)
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: rmatrix
!</input>
  
!<output>
  ! VANCA spiecific structure. Contains internal data and allocated memory.
  TYPE(t_vancaGeneral), INTENT(OUT)       :: rvancaGeneral
!</output>

!</subroutine>

    ! local variables
    LOGICAL :: bfirst
    INTEGER :: nblocks,i,j,nmaxLocalDOFs,ndofsPerElement
    TYPE(t_spatialDiscretisation), POINTER            :: p_rdiscretisation
    
    nblocks = rmatrix%ndiagBlocks
    nmaxLocalDOFs = 0
    ndofsPerElement = 0
    
    ! Allocate memory for the matrix structures
    rvancaGeneral%nblocks = nblocks
    ALLOCATE(rvancaGeneral%p_Rmatrices(nblocks,nblocks))
    
    ALLOCATE(rvancaGeneral%p_InDofsLocal(rmatrix%ndiagBlocks))
    ALLOCATE(rvancaGeneral%p_IblockOffset(rmatrix%ndiagBlocks+1))
    
    ! This type of VANCA only supports a uniform discretisation
    ! and matrix format 7 or 9. Transposed matrices are not allowed.
    !
    ! Check all matrices that this is the case.
    ! Build the Rmatrices structure array. We manually create a
    ! (nblock,nblock) array inside of a 1-dimensional array and cast a rank
    ! change of the array on call to the actual VANCA subroutine later.
    !
    ! Loop through the columns of the block matrix.
    !
    ! Offset position of the first block is = 0.
    rvancaGeneral%p_IblockOffset(1) = 0
    
    DO i=1,nblocks
    
      ! Note this block as 'not processed'
      bfirst = .TRUE.
      
      ! Loop through the rows of the current matrix column.
      DO j=1,nblocks
        IF (lsysbl_isSubmatrixPresent(rmatrix,j,i)) THEN
          ! Get a/the discretisation structure of the current block/matrix column
          p_rdiscretisation => rmatrix%RmatrixBlock(j,i)%p_rspatialDiscrTrial
          
          IF ((p_rdiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .AND. &
              (p_rdiscretisation%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
            CALL output_line (&
                'General VANCA supports only uniform and conformal discretisations!',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_initGeneralVanca')
            CALL sys_halt()
          END IF
          
          IF ((rmatrix%RmatrixBlock(j,i)%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND. &
              (rmatrix%RmatrixBlock(j,i)%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
            CALL output_line (&
                'General VANCA supports only matrix structure 7 and 9!',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_initGeneralVanca')
            CALL sys_halt()
          END IF
          
          rvancaGeneral%p_Rmatrices(j,i)%bexists = .TRUE.
          rvancaGeneral%p_Rmatrices(j,i)%dscaleFactor = &
            rmatrix%RmatrixBlock(j,i)%dscaleFactor
          rvancaGeneral%p_Rmatrices(j,i)%btransposed = &
            IAND(rmatrix%RmatrixBlock(j,i)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0
            
          CALL lsyssc_getbase_double (rmatrix%RmatrixBlock(j,i),&
                                      rvancaGeneral%p_Rmatrices(j,i)%p_DA)
          CALL lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(j,i),&
                                    rvancaGeneral%p_Rmatrices(j,i)%p_Kcol)
          CALL lsyssc_getbase_Kld (rmatrix%RmatrixBlock(j,i),&
                                   rvancaGeneral%p_Rmatrices(j,i)%p_Kld)
                                    
          IF (bfirst) THEN
            ! This block has not yet been processed.
            !
            ! Get the NEQ of the current block and save it as offset position
            ! in the global vector for the next block.
            rvancaGeneral%p_IblockOffset(i+1) = &
              rvancaGeneral%p_IblockOffset(i) + rmatrix%RmatrixBlock(j,i)%NCOLS

            ! We need some information for calculating DOF's later.
            ! Get the number of local DOF's in the current block.
            ! Note that we restrict to uniform discretisations!
            rvancaGeneral%p_InDofsLocal(i) = elem_igetNDofLoc(p_rdiscretisation% &
                                                RelementDistr(1)%celement)
            
            ! Calculate the maximum number of local DOF's
            nmaxLocalDOFs = MAX(nmaxLocalDOFs,rvancaGeneral%p_InDofsLocal(i))
            
            ! Calculate the total number of local DOF's
            ndofsPerElement = ndofsPerElement + rvancaGeneral%p_InDofsLocal(i)
            
            bfirst = .FALSE.
          
          END IF
          
        ELSE
          rvancaGeneral%p_Rmatrices(j,i)%bexists = .FALSE.
        END IF
        
      END DO
    END DO

    ! Save the max. and total number of local DOF's
    rvancaGeneral%nmaxLocalDOFs = nmaxLocalDOFs
    rvancaGeneral%ndofsPerElement = ndofsPerElement
    
    ! We know the maximum number of DOF's now. For the later loop over the 
    ! elements, allocate memory for storing the DOF's of an element set.
    ALLOCATE(rvancaGeneral%p_IelementDOFs(nmaxLocalDOFs,VANCA_NELEMSIM,nblocks))
    
    ! Remember the matrix
    rvancaGeneral%p_rmatrix => rmatrix
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_doneGeneralVanca (rvancaGeneral)

!<description>
  ! This routine cleans up a general VANCA solver. All memory allocated in
  ! rvancaGeneral is released.
!</description>
  
!<inputoutput>
  ! The general-VANCA structure to be cleaned up.
  TYPE(t_vancaGeneral), INTENT(INOUT)       :: rvancaGeneral
!</inputoutput>

!</subroutine>

    ! Release memory allocated in the init-routine
    
    IF (ASSOCIATED(rvancaGeneral%p_IelementDOFs)) &
      DEALLOCATE(rvancaGeneral%p_IelementDOFs)
    
    IF (ASSOCIATED(rvancaGeneral%p_InDofsLocal)) &
      DEALLOCATE(rvancaGeneral%p_InDofsLocal)
      
    IF (ASSOCIATED(rvancaGeneral%p_IblockOffset)) &
      DEALLOCATE(rvancaGeneral%p_IblockOffset)

    IF (ASSOCIATED(rvancaGeneral%p_Rmatrices)) &
      DEALLOCATE(rvancaGeneral%p_Rmatrices)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_general (rvancaGeneral, rvector, rrhs, domega)

!<description>
  ! This routine applies the general-VANCA algorithm to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvancaGeneral structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix A.
!</description>

!<input>
  
  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The general-VANCA structure. Must have been initialised with 
  ! vanca_initGeneralVanca before.
  TYPE(t_vancaGeneral), INTENT(INOUT)     :: rvancaGeneral

  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,j
  INTEGER(PREC_ELEMENTIDX) :: IELmax, IELset, iel, ieldistr
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementList
  REAL(DP), DIMENSION(:), POINTER                 :: p_Drhs,p_Dvector
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER     :: p_Ipermutation
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
  ! Saved matrix and the vector(s) must be compatible!
  CALL lsysbl_isMatrixCompatible(rvector,rvancaGeneral%p_rmatrix)
  CALL lsysbl_isVectorCompatible(rvector,rrhs)
    
  ! Get the data arrays of the vector/rhs
  CALL lsysbl_getbase_double (rvector,p_Dvector)
  CALL lsysbl_getbase_double (rrhs,p_Drhs)
  
  ! Get the discretisation structure that tells us which elements form
  ! element groups...
  p_rdiscretisation => rvector%RvectorBlock(1)%p_rspatialDiscr
  
  ! Loop over the element distributions/groups
  DO ieldistr = 1,p_rdiscretisation%inumFEspaces
  
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions.
    CALL storage_getbase_int (p_rdiscretisation% &
                              RelementDistr(ieldistr)%h_IelementList, &
                              p_IelementList)
      
    ! Loop over the elements - blockwise.
    DO IELset = 1, SIZE(p_IelementList), VANCA_NELEMSIM
    
      ! We always handle LINF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most VANCA_NELEMSIM
      ! elements simultaneously.
      IELmax = MIN(SIZE(p_IelementList),IELset-1+VANCA_NELEMSIM)
    
      ! Loop over the blocks in the block vector to get the DOF's everywhere.
      
      DO i=1,rvector%nblocks

        ! Calculate the global DOF's of all blocks.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our VANCA_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(rvector%RvectorBlock(i)%p_rspatialDiscr,&
                                     p_IelementList(IELset:IELmax), &
                                     rvancaGeneral%p_IelementDOFs(:,:,i))

        ! If the vector is sorted, push the DOF's through the permutation to get
        ! the actual DOF's.
        IF (rvector%RvectorBlock(i)%isortStrategy .GT. 0) THEN
        
          CALL storage_getbase_int(rvector%RvectorBlock(i)%h_IsortPermutation,&
                                  p_Ipermutation)

          DO iel=1,IELmax-IELset+1
            DO j=1,rvancaGeneral%p_InDofsLocal(i)
              ! We are not resorting the vector but determining the 'sorted'
              ! DOF's - this needs the 2nd half of the permutation.
              rvancaGeneral%p_IelementDOFs(j,iel,i) = &
                 p_Ipermutation(rvancaGeneral%p_IelementDOFs(j,iel,i) &
                +rvector%RvectorBlock(i)%NEQ)
            END DO
          END DO
        END IF
      
      END DO  
    
      ! Now, IdofsTotal contains all DOF's on each element, over all discretisations.
      !
      ! Call the actual VANCA to process the DOF's on each element.
      CALL vanca_general_double_mat79 (p_Dvector, p_Drhs, domega, &
          rvancaGeneral%p_Rmatrices,IELmax-IELset+1_PREC_ELEMENTIDX,&
          rvancaGeneral%p_IblockOffset,rvancaGeneral%nblocks,&
          rvancaGeneral%p_InDofsLocal,rvancaGeneral%ndofsPerElement,&
          rvancaGeneral%p_IelementDOFs)
                   
    END DO
    
  END DO ! ieldistr
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_general_double_mat79 (Dvector, Drhs, domega, Rmatrices,&
             nelements,IblockOffset,nblocks,InDofsLocal,ndofsPerElement,&
             IelementDofs)

!<description>
  ! This routine applies one step of the VANCA solver (local block 
  ! gauss-seidel) to a given set of solution/RHS vector. 
  ! All given matrices are assumed to be format 7/9,
  ! double precision. The given vector is assumed to be double precision.
!</description>

!<input>
  ! The (block) RHS vector, given as one large array.
  REAL(DP), DIMENSION(:), INTENT(IN)             :: Drhs
  
  ! A relaxation parameter. Standard = 1.0_DP.
  REAL(DP), INTENT(IN)                           :: domega
  
  ! A list of matrices to handle; directly specified by pointers
  ! to the substructures (data/columns/rows).
  TYPE(t_matrixPointer79Vanca), DIMENSION(:,:),&
                                INTENT(IN)       :: Rmatrices

  ! Number of elements that should be processed in this sweep.
  INTEGER(PREC_ELEMENTIDX), INTENT(IN)           :: nelements
  
  ! Number of blocks in the vectors
  INTEGER, INTENT(IN)                            :: nblocks
  
  ! Offset position of the blocks in the vector.
  ! Block i starts at position IblockOffset(i)+1 in Dvector / Drhs.
  ! IblockOffset(nblocks+1) gives the number of equations in Dvector/Drhs.
  INTEGER(PREC_DOFIDX), DIMENSION(MAX(nblocks,1)+1), INTENT(IN) :: IblockOffset
  
  ! Number of local DOF's in each block.
  INTEGER(PREC_DOFIDX), DIMENSION(nblocks), INTENT(IN)   :: InDofsLocal
  
  ! Total number of local DOF's per element
  INTEGER, INTENT(IN)                                    :: ndofsPerElement
  
  ! List of DOF's on every element for every block.
  ! DIMENSION(nmaxDOFs,nelements,nblocks)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), INTENT(IN)     :: IelementDOFs
  
!</input>

!<inputoutput>
  ! The initial (block) solution vector. Is overwritten by the new (block) 
  ! solution vector.
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dvector
!</inputoutput>

!</subroutine>

    ! One iteration of vanka smother (block gauss-seidel) on the system
    !
    !    A11 A12 ... A1n | U1    F1
    !    A21 A22 ... A2n | U2  = F2
    !     :   :  ...  :  |  :     :
    !    An1 An2 ... Ann | Un  = Fn
    !
    ! Let's first describe the method used here.
    !
    ! The above block system can be written as a general block system of 
    ! the form
    !
    !             A x = f
    !
    ! Consider a general defect correction approach for this system:
    !
    !     x_{n+1}  =  x_n  +  \omega C^{-1} (b - A x_n)
    !
    ! With some damping parameter \omega and some preconditioner C^{-1}.
    ! Normally, this iteration is used globally - and the usual linear
    ! solver that calls this routine usually realises this approach.
    ! VANCA now applies this defect correction loop *locally*, i.e.
    ! not for the full system at once.
    !
    ! This local approach is based on a geometric point of view.
    ! In general, one could imagine a set of connected cells in the 
    ! global domain \Omega where to apply the algorithm to (the LMPSC 
    ! approach), but for simplicity, consider only one cell. Again for
    ! simplicity, imagine that our FE-spaces is Q1~/Q0. 
    !
    ! We loop over each cell in the domain, one after the other, and
    ! change the DOF's in the solution vector there:
    !
    ! We fetch all the data (e.g. velocity, pressure) on that cell. On the
    ! first cell, we have only "old" velocity entries. These values
    ! are updated and then the calculation proceeds with the 2nd cell.
    !
    !        old                      new
    !     |---X---|                |---X---|
    !     |       |                |       |
    ! old X   1   X old   -->  new X   1   X new
    !     |       |                |       |
    !     |---X---|                |---X---|
    !        old                      new
    !
    ! From the second cell on, there might be "old" data and "new" 
    ! data on that cell - the old data that has not been updated and
    ! perhaps some already updated velocity data from a neighbor cell.
    ! 
    !        new     old                   new     new
    !     |---X---|---X---|             |---X---|---X---|
    !     |       |       |             |       |       |
    ! new X   1   X   2   X old --> new X   1   X   2   X new
    !     |       |new    |             |       |newer  |
    !     |---X---|---X---|             |---X---|---X---|
    !        new     old                   new     new
    !
    ! These values are updated and then the calculation proceeds
    ! with the next cell.
    ! As can be seen in the above picture, the "new" node in the
    ! middle is even going to be a "newer" node when handled again
    ! for the 2nd cell. This is meant by "Gauss-Seldel" character:
    ! Information is updated subsequently by using "old" data and
    ! "new" data from a previous calculation.
    !
    ! So how to do this in detail?
    !
    ! Again consider the problem:
    !
    !             A x = f
    !
    ! We assume, that all components in the vector x are      
    ! given - except for the unknowns current element; these unknowns
    ! are located anywhere in the vector x. The idea is to    
    ! shift "everything known" to the right hand side to obtain   
    ! a system for only these unknowns!  
    !                         
    ! We extract all the lines of the system that correspond to
    ! our 'unknown' DOF's on our element; this results in a rectangular      
    ! system of the form                                            
    !                                                               
    !    [ === A~ === ] x  = (f~)
    !
    ! So #rows(A~)=#rows(f~)=#DOF's on the element! Furthermore we set
    ! Now we make a defect-correction approach for this system:
    !
    !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
    !                                     -----------
    !                                        =d~
    !
    ! Here the 'projection' operator simply converts the small
    ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
    ! of the same size as x - what is easy using the number of
    ! the DOF's on the element.
    !
    ! The only question now will be: What is C^{-1}?
    !
    ! Well, here we have different choices. Everything depends on the
    ! matrix A~, which is unfortunately rectangular: It's a
    ! (#DOF's on the element, #DOF's in the space) matrix - so
    ! in case of a Q1~/Q0 discretisation, it's a (5,NEQ) matrix!
    !
    ! For full linear systems, one would choose C=A, which ist the
    ! theoretically best preconditioner. What we do here is simply
    ! extracting all columns of A~ that correspond to the DOF's
    ! on the current element: 
    !
    !   C:=delete columns of A~ that don't belong to DOF's on the element
    !
    ! This then leads to a square preconditioner C^{-1} - and that's the
    ! full method, because C^{-1} can be applied directly using Lapack e.g.!
    !
    !
    ! So all in all we write our algorithm in short:
    !
    ! loop over all elements in the given element set
    !   extract local matrix aa
    !   extract the local rhs ff
    !   compute local residuum: ff := ff - aa * x
    !   solve: xx := aa^(-1) ff
    !   update global solution vector: x := x + omega*xx
    ! end loop

    ! local variables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER, DIMENSION(MAX(nblocks,1)+1) :: IlocalIndex
    INTEGER(PREC_DOFIDX), DIMENSION(ndofsPerElement) :: IlocalDOF,IglobalDOF
    INTEGER :: i,j,k,iidx,iminiDOF
    INTEGER(PREC_VECIDX) :: irow,idof
    INTEGER(PREC_MATIDX) :: icol
    
    ! Memory for our local system; let's hope it's not too big :)
    REAL(DP), DIMENSION(ndofsPerElement,ndofsPerElement) :: Daa
    REAL(DP), DIMENSION(ndofsPerElement)                 :: Dff
    REAL(DP) :: dscale
    INTEGER(I32), DIMENSION(ndofsPerElement) :: Ipiv
    INTEGER(I32) :: iinfo
    
    ! Quickly check the matrices if one of them is saved transposed.
    ! VANCA does not support transposed matrices; would kill computational
    ! time, as we cannot extract columns from a structure-9 matrix with
    ! reasonable effort!
    DO i=1,SIZE(Rmatrices,1)
      DO j=1,SIZE(Rmatrices,2)
        IF (Rmatrices(i,j)%bexists .AND. Rmatrices(i,j)%btransposed) THEN
          CALL output_line (&
              'General VANCA does not support transposed matrices!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanca_general_double_mat79')
          CALL sys_halt()
        END IF 
      END DO ! j
    END DO ! i
        
    ! Build an index pointer for accessing local DOF's
    IlocalIndex(1) = 0
    DO i=2,nblocks+1
      IlocalIndex(i) = IlocalIndex(i-1)+InDOFsLocal(i-1)
    END DO
        
    ! Ok, let's start with the loop over the elements in the given
    ! element list.
    DO iel = 1,nelements
    
      ! IelementDOFs (.,iel,.) gives now for every block in the system
      ! the DOF's on this element.
      !
      ! First copy the RHS entries of f to f~.
      
      iidx = 1
      DO i=1,nblocks
        DO j=1,InDOFsLocal(i)
          iidx = IlocalIndex(i)+j
          
          ! Get the DOF on the element:
          idof = IelementDOFs(j,iel,i)
          
          ! Number of DOF relative to this block
          IlocalDOF(iidx) = idof
          
          ! Calculate corresponding global DOF:
          IglobalDOF(iidx) = idof+IblockOffset(i)
          
          Dff(iidx) = Drhs(IglobalDOF(iidx))
        END DO
      END DO
      
      ! Compute  ff := ff - A x  for the local unknowns
      ! to build the local residuum Dff = f~-A~x.
      !
      ! Simultaneously extract the local matrix into one array Daa(:,:).
      ! But first clear the matrix - maybe that there are unused subblocks!
      Daa = 0
      
      ! Loop through the rows of the block matrix
      DO i=1,SIZE(Rmatrices,1)
        
        ! Loop through the columns of the block matrix
        DO j=1,SIZE(Rmatrices,2)
        
          ! Is there a matrix saved at this position?
          IF (Rmatrices(i,j)%bexists .AND. &
              Rmatrices(i,j)%dscaleFactor .NE. 0.0_DP) THEN

            dscale = Rmatrices(i,j)%dscaleFactor
        
            ! Block column j in the global matrix corresponds to
            ! block column j in the small local matrix we have to fill with data.
            !
            !   ....      :           
            !   ....      :           a11 a12         a15
            !                         a21 a22         a25
            !        .... :   ==>             a33 a34 a35
            !        .... :                   a43 a44 a45
            !   .... ....             a51 a52 a53 a54
            !   ==== ==== =           ======= ======= ===
            !    1    2   3 corr. to     1       2     3
            !
            !     n x n-mat.          5x5-mat or so
            !
            ! From IlocalIndex, get the starting address of the j'th block
            ! in the local solution vector. In the above example, for j=2 e.g.
            ! this gives the starting address of the two global DOF's
            ! that correspond to the columns 3 and 4 in the local matrix.
            iidx = IlocalIndex(j)

            ! Loop through the DOF's that correspond to this block:
            DO irow = IlocalIndex(i)+1,IlocalIndex(i+1)
              
              ! Get the actual DOF, relative to this block.
              ! This is the row of the matrix that must be multiplied by x
              ! and be subtracted from the local RHS.
              idof = IlocalDOF(irow)
              
              ! Loop through the row of the matrix to its contribution to "b-Ax".
              DO k = Rmatrices(i,j)%p_Kld(idof) , Rmatrices(i,j)%p_Kld(idof+1)-1

                ! Get the column number in the global matrix. This gives
                ! the global DOF = number of the element in x that must be multiplied
                ! with that matrix entry.             
                icol = Rmatrices(i,j)%p_Kcol(k)+IblockOffset(j)

                ! Build the defect
                Dff(irow) = Dff(irow) &
                          - dscale * Rmatrices(i,j)%p_Da(k) * Dvector(icol)

                ! icol is the number of a DOF.
                ! Check if this DOF belongs to the DOF's we have to
                ! extract to our local system.
                ! Loop through the DOF's corresponding to column j of the block system.
                ! In the above example, this checks only the two DOF's corresponding
                ! to a?3 and a?4.
                DO iminiDOF = 1,InDOFsLocal(j)
                  IF (icol .EQ. IglobalDOF(iidx+iminiDOF)) THEN
                    ! Yes. Get the matrix entry, write it to the local matrix
                    Daa(irow,iidx+iminiDOF) = dscale*Rmatrices(i,j)%p_Da(k)
                    EXIT
                  END IF
                END DO ! idof
                
              END DO ! k
            
            END DO ! irow
                        
          END IF ! exists
        
        END DO ! j
        
      END DO ! i
      
      ! Ok, we now have our local matrix and our local defect vector d~.
      ! Apply LAPACK to the local system to solve C^{-1} d~.
      
      CALL DGETRF( ndofsPerElement, ndofsPerElement, Daa, ndofsPerElement, &
                   Ipiv, iinfo )
                  
      ! Note: It may happen that the matrix is singular!
      !  That is the case if all DOF's are Dirichlet-DOF's - for example
      !  if the element is completely inside of a rigid fictitious boundary
      !  object.
      ! What to do in this case? Nothing! Ignore the system!
      ! Why? 
      !  - The values for the 'velocity' DOF's are fixed, so it's no use
      !    to try to calculate them.
      !  - If there's a zero-block involved in a saddle-point problem,
      !    the pressure (that corresponds to the zero-block) is not
      !    connected to the velocity - it's completely free!
      !    By ignoring this system, we let the pressure as it is.
      !
      ! One can theoretically also solve a least-squares system with
      ! LAPACK. This calculate an y with ||f~ - C y|| -> min and |y|->min.
      ! But this is (because of the 0-block in A and therefore also in C
      ! and because of the unit-vectors in the other rows of C)
      ! exactly the case if the 'velocity' y matches f~ and the 'pressure'
      ! components are zero - which means nothing else than that there's
      ! no contribution in the preconditioned defect to correct the
      ! pressure entries of x.
                   
      IF (iinfo .EQ. 0) THEN
        CALL DGETRS('N', ndofsPerElement, 1, Daa, ndofsPerElement, &
                    Ipiv, Dff, ndofsPerElement, iinfo )
        IF (iinfo .EQ. 0) THEN
          
          ! Dff is in-situ replaced by the solution - the preconditioned
          ! defect. Add this to our original solution vector to complete
          ! the 'local' defect-correction.
          
          DO i=1,ndofsPerElement
            j = IglobalDOF(i)
            Dvector(j) = Dvector(j) + domega*Dff(i)
          END DO
        
        END IF
      END IF
    
    END DO ! iel

  END SUBROUTINE

! *****************************************************************************
! Problem class: VANCA variants for 2D Navier-Stokes problems
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_init2DNavierStokes (rmatrix,rvanca,csubtype)
  
!<description>
  ! Initialises the VANCA variant for 2D Navier-Stokes problems 
  ! for conformal discretisations.
  ! Checks if the "2D-Navier-Stokes" VANCA variant 
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! The substructure rvanca%rvanca2DNavSt is intitialised according
  ! to the information provided in rmatrix.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: rmatrix

  ! Desired subtype
  INTEGER, INTENT(IN) :: csubtype  
!</input>

!<inputoutput>
  ! t_vancaPointer2DSPNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vanca), INTENT(INOUT) :: rvanca
!</inputoutput>

!</subroutine>

    INTEGER :: i,j
    LOGICAL :: bextended
    TYPE(t_blockDiscretisation), POINTER :: p_rblockDiscr
    
    bextended = .FALSE.
    
    ! Matrix must be 3x3.
    IF (rmatrix%ndiagBlocks .NE. 3) THEN
      CALL output_line ('System matrix is not 3x3.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
      CALL sys_halt()
    END IF
    
    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(3,:) must be (virtually) transposed. All matrices must be double precision.
    DO i=1,3
      DO j=1,3
      
        IF (lsysbl_isSubmatrixPresent(rmatrix,i,j)) THEN
        
          IF (i .LE. 2) THEN
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .NE. 0) THEN
              CALL output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
              CALL sys_halt()
            END IF
          ELSE
            IF ((i .LE. 2) .AND. &
               (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .EQ. 0)) THEN
              CALL output_line ('B1/B2 submatrices must be virtually',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
              CALL output_line ('transposed (LSYSSC_MSPEC_TRANSPOSED)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
              CALL sys_halt()
            END IF
          END IF
          
          IF ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX7) .AND. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
            CALL output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
            CALL sys_halt()
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            CALL output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
            CALL sys_halt()
          END IF

          ! For scaled matrices, we have to use an extended sub-version of VANCA.
          IF ((rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) .AND. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 0.0_DP)) THEN
            bextended = .TRUE.  
          END IF
          
          IF ((i .eq. j) .AND. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) ) THEN
            bextended = .TRUE. 
          END IF
          
        END IF ! neq != 0
      END DO
    END DO
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    IF ((rmatrix%RmatrixBlock(1,3)%NA .NE. rmatrix%RmatrixBlock(3,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .NE. rmatrix%RmatrixBlock(3,1)%NCOLS)) THEN
      CALL output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
      CALL sys_halt()
    END IF

    IF ((rmatrix%RmatrixBlock(2,3)%NA .NE. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .NE. rmatrix%RmatrixBlock(3,2)%NCOLS)) THEN
      CALL output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
      CALL sys_halt()
    END IF      
  
    ! Fill the output structure with data of the matrices.
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanca%rvanca2DNavSt%p_DA )
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),&
        rvanca%rvanca2DNavSt%p_DB1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),&
        rvanca%rvanca2DNavSt%p_DB2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),&
        rvanca%rvanca2DNavSt%p_DD1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),&
        rvanca%rvanca2DNavSt%p_DD2)
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),&
        rvanca%rvanca2DNavSt%p_KcolB)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), &
        rvanca%rvanca2DNavSt%p_KldB )
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanca%rvanca2DNavSt%p_KcolA)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanca%rvanca2DNavSt%p_KldA )
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanca%rvanca2DNavSt%p_KdiagonalA)
    ELSE
      rvanca%rvanca2DNavSt%p_KdiagonalA => rvanca%rvanca2DNavSt%p_KldA
    END IF
    
    IF (lsysbl_isSubmatrixPresent(rmatrix,3,3)) THEN
    
      ! The matrix must be of format 7 or 9.
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
          rvanca%rvanca2DNavSt%p_DA33 )

      IF (rmatrix%RmatrixBlock(3,3)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(3,3), &
                                rvanca%rvanca2DNavSt%p_KdiagonalA33)
      ELSE
        CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,3), &
                                rvanca%rvanca2DNavSt%p_KdiagonalA33)
      END IF

      ! The presence of A(3,3) forces the extended VANCA to be used
      bextended = .TRUE.

    END IF
    
    IF (bextended) rvanca%csubsubtype = 1
    
    ! What is with A22? Is it the same as A11?
    IF (.NOT. lsyssc_isMatrixContentShared (&
        rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2)) ) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
          rvanca%rvanca2DNavSt%p_DA22 )
    END IF
    
    ! What is with A12 and A21? Do they exist? With a scale factor = 1.0?
    IF (lsysbl_isSubmatrixPresent(rmatrix,1,2) .AND. &
        (rmatrix%RmatrixBlock(1,2)%dscaleFactor .EQ. 1.0_DP)) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
          rvanca%rvanca2DNavSt%p_DA12 )
          
      CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,2),&
          rvanca%rvanca2DNavSt%p_KcolA12)
      CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,2), &
          rvanca%rvanca2DNavSt%p_KldA12 )
          
      ! Get the structure. It's assumed that A12 and A21 have the same!
      IF (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,2), &
                                rvanca%rvanca2DNavSt%p_KdiagonalA12)
      ELSE
        rvanca%rvanca2DNavSt%p_KdiagonalA12 => rvanca%rvanca2DNavSt%p_KldA12
      END IF
      
      IF (.NOT. lsysbl_isSubmatrixPresent(rmatrix,2,1)) THEN
        CALL output_line ('If A12 is given, A21 must also be given!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesO')
        CALL sys_halt()
      END IF
      
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
          rvanca%rvanca2DNavSt%p_DA21 )
    END IF

    ! Get the multiplication factors of the submatrices.
    rvanca%rvanca2DNavSt%Dmultipliers(1:3,1:3) = &
        rmatrix%RmatrixBlock(1:3,1:3)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest
    
    IF (.NOT. ASSOCIATED(p_rblockDiscr)) THEN
      CALL output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
      CALL sys_halt()
    END IF
    
    ! Get the discretisation structure of U,V and P from the block
    ! discretisation structure.
    rvanca%rvanca2DNavSt%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanca%rvanca2DNavSt%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanca%rvanca2DNavSt%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)
    
    IF (rvanca%rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces .NE. &
        rvanca%rvanca2DNavSt%p_rspatialDiscrV%inumFESpaces) THEN
      CALL output_line (&
          'Discretisation structures of X- and Y-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
      CALL sys_halt()
    END IF

    IF ((rvanca%rvanca2DNavSt%p_rspatialDiscrP%inumFESpaces .NE. 1) .AND. &
        (rvanca%rvanca2DNavSt%p_rspatialDiscrP%inumFESpaces .NE. &
          rvanca%rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces)) THEN
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      CALL output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokes')
      CALL sys_halt()
    END IF
    
    ! Solution based VANCA variants need an additional temp vector
    IF (csubtype .EQ. VANCATP_DIAGONAL_SOLBASED) then
      CALL lsysbl_createVecBlockIndMat (rmatrix,rvanca%rvanca2DNavSt%rtempVector, .false.)
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_done2DNavierStokes (rvanca2DNavSt)
  
!<description>
  ! Releases any memory allocated in rvanca2DNavSt.
!</description>

!<inputoutput>
  ! t_vanca structure that to be cleaned up.
  TYPE(t_vancaPointer2DNavSt), INTENT(INOUT) :: rvanca2DNavSt
!</inputoutput>

!</subroutine>

    IF (rvanca2DNavSt%rtempVector%NEQ .NE. 0) THEN
      CALL lsysbl_releaseVector (rvanca2DNavSt%rtempVector)
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DNavierStokes (rvanca2DNavSt, rvector, rrhs, domega, &
      csubtype, csubsubtype)
  
!<description>
  ! This routine applies the VANCA variant for 2D Navier-Stokes problems
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANCA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! The subtype of VANCA that should handle the above problem class.
  ! One of the VANCATP_xxxx constants, e.g. VANCATP_DIAGONAL.
  INTEGER :: csubtype

  ! The sub-subtype of VANCA that should handle the above problem class.
  ! =0: use standard VANCA. =1: use extended VANCA (e.g. with different 
  !     multipliers in the matrices)
  INTEGER :: csubsubtype
  
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(INOUT)         :: rvector

  ! t_vanca structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(INOUT) :: rvanca2DNavSt
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ielementdist
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementList
    TYPE(t_elementDistribution), POINTER :: p_relementDistrU
    TYPE(t_elementDistribution), POINTER :: p_relementDistrV
    TYPE(t_elementDistribution), POINTER :: p_relementDistrP 
    
    ! 2D Navier Stokes problem.

    ! Loop through the element distributions of the velocity.
    DO ielementdist = 1,rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces
    
      ! Get the corresponding element distributions of U, V and P.
      p_relementDistrU => &
          rvanca2DNavSt%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanca2DNavSt%p_rspatialDiscrV%RelementDistr(ielementdist)
      
      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      IF (rvanca2DNavSt%p_rspatialDiscrP%inumFESpaces .GT. 1) THEN
        p_relementDistrP => &
            rvanca2DNavSt%p_rspatialDiscrP%RelementDistr(ielementdist)
      ELSE
        p_relementDistrP => &
            rvanca2DNavSt%p_rspatialDiscrP%RelementDistr(1)
      END IF
      
      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y-velocity (and to that
      ! of the pressure).
      CALL storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)
      
      ! Which element combination do we have now?
      IF ((elem_getPrimaryElement(p_relementDistrU%celement) .EQ. EL_Q1T) .AND. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .EQ. EL_Q1T) .AND. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .EQ. EL_Q0)) THEN
        ! Q1~/Q1~/Q0 discretisation
        
        ! Which VANCA subtype do we have? The diagonal VANCA of the full VANCA?
        SELECT CASE (csubtype)
        CASE (VANCATP_DIAGONAL)
          ! Diagonal VANCA; check if A12 exists.
          IF (.NOT. ASSOCIATED(rvanca2DNavSt%p_DA12)) THEN
            ! Call the VANCA subsolver to apply VANCA to our current element list.
            IF (rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces .EQ. 1) THEN
              ! Uniform discretisation
              CALL vanca_2DSPQ1TQ0simple (rvanca2DNavSt, &
                  rvector, rrhs, domega)
            ELSE
              ! Conformal discretisation
              CALL vanca_2DSPQ1TQ0simpleConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
          ELSE
            ! Apply the conformal VANCA that allows different matrices
            ! in A11, A12, A21 and A22!
            CALL vanca_2DSPQ1TQ0simpleCoupConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
          END IF
          
        CASE (VANCATP_DIAGONAL_SOLBASED)
          ! Diagonal VANCA. Solution based if possible. Check if A12 exists.
          IF (.NOT. ASSOCIATED(rvanca2DNavSt%p_DA12)) THEN
            ! Call the VANCA subsolver to apply VANCA to our current element list.
            IF (rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces .EQ. 1) THEN
              ! Uniform discretisation
              CALL vanca_2DSPQ1TQ0simpleSol (rvanca2DNavSt, &
                  rvector, rrhs, domega,rvanca2DNavSt%rtempVector)
            ELSE
              ! Conformal discretisation
              CALL vanca_2DSPQ1TQ0simpleConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
          ELSE
            ! Apply the conformal VANCA that allows different matrices
            ! in A11, A12, A21 and A22!
            CALL vanca_2DSPQ1TQ0simpleCoupConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
          END IF

        CASE (VANCATP_FULL)
          ! Full VANCA; check if A12 exists.
          IF (.NOT. ASSOCIATED(rvanca2DNavSt%p_DA12)) THEN
            ! Call the VANCA subsolver to apply VANCA to our current element list.
            ! Note: Up to now, there is no 'full' variant -- has to be implemented!
            IF (rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces .EQ. 1) THEN
              ! uniform discretisation;
              ! here, use the same as for the general conformal discretisation.
              ! Could be speeded up by introducing another variant...
              CALL vanca_2DSPQ1TQ0fullConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            ELSE
              ! Conformal discretisation
              CALL vanca_2DSPQ1TQ0fullConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
          ELSE
            ! Apply the conformal VANCA that allows different matrices
            ! in A11, A12, A21 and A22!
            ! If we have multiplication factors, we even have to use an extended
            ! version of this.
            IF (csubsubtype .EQ. 0) THEN
              CALL vanca_2DSPQ1TQ0fullCoupConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            ELSE
              CALL vanca_2DNSQ1TQ0fullCoupConfExt (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
          END IF
        
        CASE DEFAULT
          CALL output_line (&
              'Unknown VANCA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNavierStokes')
          CALL sys_halt()  
        
        END SELECT
        
      ELSE IF &
        ((elem_getPrimaryElement(p_relementDistrU%celement) .EQ. EL_Q2) .AND.&
          (elem_getPrimaryElement(p_relementDistrV%celement) .EQ. EL_Q2) .AND.&
          (elem_getPrimaryElement(p_relementDistrP%celement) .EQ. EL_QP1)) THEN
        ! Q2/Q2/QP1 discretisation
        
        ! Which VANCA subtype do we have? The diagonal VANCA of the full VANCA?
        SELECT CASE (csubtype)
        CASE (VANCATP_DIAGONAL)
          ! Diagonal VANCA; check if A12 exists.
          IF (.NOT. ASSOCIATED(rvanca2DNavSt%p_DA12)) THEN
            ! Call the VANCA subsolver to apply VANCA to our current element list.
            IF (rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces .EQ. 1) THEN
              ! Uniform discretisation
              CALL vanca_2DSPQ2QP1simple (rvanca2DNavSt, &
                  rvector, rrhs, domega)
            ELSE
              ! Conformal discretisation
              CALL vanca_2DSPQ2QP1simpleConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
          ELSE
            ! Apply the conformal VANCA that allows different matrices
            ! in A11, A12, A21 and A22!
            CALL output_line (&
                'VANCA has no support for A12 and A21!',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNavierStokes')
            CALL sys_halt()
          END IF
          
        CASE (VANCATP_FULL)
          ! Full VANCA; check if A12 exists.
          IF (.NOT. ASSOCIATED(rvanca2DNavSt%p_DA12)) THEN
            ! Call the VANCA subsolver to apply VANCA to our current element list.
            IF (rvanca2DNavSt%p_rspatialDiscrU%inumFESpaces .EQ. 1) THEN
              ! Uniform discretisation
              CALL vanca_2DSPQ2QP1full (rvanca2DNavSt, &
                  rvector, rrhs, domega)
            ELSE
              ! Conformal discretisation
              CALL vanca_2DSPQ2QP1fullConf (rvanca2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
          ELSE
            ! Apply the conformal VANCA that allows different matrices
            ! in A11, A12, A21 and A22!
            CALL output_line (&
                'VANCA has no support for A12 and A21!',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNavierStokes')
            CALL sys_halt()
          END IF
          
        CASE DEFAULT
          CALL output_line (&
              'Unknown VANCA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNavierStokes')
          CALL sys_halt()  
          
        END SELECT
          
      ELSE
        CALL output_line (&
            'Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNavierStokes')
        CALL sys_halt()
        
      END IF
      
    END DO
      
  END SUBROUTINE

  ! ***************************************************************************
  ! 2D Navier-Stokes VANCA, simple diagonal version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_init2DSPQ1TQ0simple (rmatrix,rvanca)
  
!<description>
  ! Checks if the "2D-Navier-Stokes-Q1T-Q0" VANCA variant can be applied to
  ! the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>

!<output>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(OUT) :: rvanca
!</output>

!</subroutine>

    INTEGER :: i,j

    ! Matrix must be 3x3.
    IF (rmatrix%ndiagBlocks .NE. 3) THEN
      CALL output_line (&
          'System matrix is not 3x3.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_check2DSPQ1TQ0')
      CALL sys_halt()
    END IF
    
    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(3,:) must be (virtually) transposed. All matrices must be double precision.
    DO i=1,3
      DO j=1,3
      
        IF (lsysbl_isSubmatrixPresent(rmatrix,i,j)) THEN
        
          IF (i .LE. 2) THEN
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .NE. 0) THEN
              CALL output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DSPQ1TQ0simple')
              CALL sys_halt()
            END IF
          ELSE
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .EQ. 0) THEN
              CALL output_line ('B1/B2 submatrices must be virtually '//&
                  'transposed (LSYSSC_MSPEC_TRANSPOSED)',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DSPQ1TQ0simple')
              CALL sys_halt()
            END IF
          END IF
          
          IF ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX7) .AND. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
            CALL output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DSPQ1TQ0simple')
            CALL sys_halt()
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            CALL output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DSPQ1TQ0simple')
            CALL sys_halt()
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) THEN
            CALL output_line ('Scaled matrices not supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DSPQ1TQ0simple')
            CALL sys_halt()
          END IF
          
        END IF ! neq != 0
      END DO
    END DO
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    IF ((rmatrix%RmatrixBlock(1,3)%NA .NE. rmatrix%RmatrixBlock(3,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .NE. rmatrix%RmatrixBlock(3,1)%NCOLS)) THEN
      CALL output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DSPQ1TQ0simple')
      CALL sys_halt()
    END IF

    IF ((rmatrix%RmatrixBlock(2,3)%NA .NE. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .NE. rmatrix%RmatrixBlock(3,2)%NCOLS)) THEN
      CALL output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DSPQ1TQ0simple')
      CALL sys_halt()
    END IF
    
    ! Fill the output structure with data of the matrices.
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),rvanca%p_DA )
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),rvanca%p_DB1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),rvanca%p_DB2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),rvanca%p_DD1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),rvanca%p_DD2)
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),rvanca%p_KcolB)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), rvanca%p_KldB )
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),rvanca%p_KcolA)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), rvanca%p_KldA )
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                               rvanca%p_KdiagonalA)
    ELSE
      rvanca%p_KdiagonalA => rvanca%p_KldA
    END IF
    
    ! Get the multiplication factors of the submatrices
    rvanca%Dmultipliers(:,:) = rmatrix%RmatrixBlock(1:3,1:3)%dscaleFactor

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ1TQ0simple (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8
    REAL(DP) :: daux
    !REAL(DP), DIMENSION(3,3) :: Dmult
    
    ! Local arrays for informations about one element
    REAL(DP), DIMENSION(4) :: AA,BB1,BB2,DD1,DD2
    REAL(DP), DIMENSION(9) :: FF,UU
    INTEGER(PREC_VECIDX), DIMENSION(4) :: idofGlobal
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanca%Dmultipliers(:,:)
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    DO iel=1,NEL
    
      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF's    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      DO inode=1,4
      
        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        
        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANCA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let's hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the four velocity unknowns and the         
        ! pressure unknown on the current element; these five unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity      
        ! components and B being an (2*4) x 1 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |        |             |        |        |                
        !   --|---X----|--           |--------|--------|                
        !                                                               
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          
          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
          
          ! For support of scaled matrices, use the following lines; currently switched off.
          !FF(inode)       = FF(inode)      -Dmult(1,3)*p_DB1(ib)*daux
          !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,3)*p_DB2(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
          
            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            
            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            
            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,3)*p_DB1(ib)
            !BB2(inode) = Dmult(2,3)*p_DB2(ib)
            !DD1(inode) = Dmult(3,1)*p_DD1(ib)
            !DD2(inode) = Dmult(3,2)*p_DD2(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,4
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ1TQ0simpleSol (rvanca, rvector, rrhs, domega, rtempVector)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! Solution-based algorithm which does not use the defect correction
  ! approach (for FEAT1-compatibility).
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(INOUT)         :: rvector
  
  ! A temporary vector in the size and structure of rvector
  TYPE(t_vectorBlock), INTENT(INOUT)         :: rtempVector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8
    REAL(DP) :: daux
    !REAL(DP), DIMENSION(3,3) :: Dmult
    
    ! Local arrays for informations about one element
    REAL(DP), DIMENSION(4) :: AA,BB1,BB2,DD1,DD2
    REAL(DP), DIMENSION(9) :: FF,UU
    INTEGER(PREC_VECIDX), DIMENSION(4) :: idofGlobal
    
    ! WARNING: DOCUMENTATION PARTIALLY WRONG AND INCOMPLETE!
    ! Preconditioner was build from FEAT1 in a quick-and-dirty way...
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanca%Dmultipliers(:,:)
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Backup the solution vector
    CALL lsysbl_copyVector (rvector,rtempVector)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    DO iel=1,NEL
    
      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF's    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      DO inode=1,4
      
        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        
        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANCA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let's hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the four velocity unknowns and the         
        ! pressure unknown on the current element; these five unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity      
        ! components and B being an (2*4) x 1 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |        |             |        |        |                
        !   --|---X----|--           |--------|--------|                
        !                                                               
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          
          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        END DO
        
        ! The diagonal entry was also subtracted, but this is not
        ! desired in this approach. We therefore revert the subtraction
        ! of the diagonal element.
        FF(inode)       = FF(inode)       + AA(inode)*p_Dvector(idof)
        FF(inode+lofsv) = FF(inode+lofsv) + AA(inode)*p_Dvector(idof+ioffsetv)
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1

        DO ib = ib1,ib2
          ! Subtract contributions from the RHS which don't belong to our element.
          IF (p_KcolB(ib) .NE. IEL) THEN 
          
            J = p_KcolB(ib)
            daux = p_Dvector(j+ioffsetp)
            FF(inode)       = FF(inode)      -p_DB1(ib)*daux
            FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
            
            ! For support of scaled matrices, use the following lines; currently switched off.
            !FF(inode)       = FF(inode)      -Dmult(1,3)*p_DB1(ib)*daux
            !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,3)*p_DB2(ib)*daux
          
          ELSE
          
            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            
            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            
            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,3)*p_DB1(ib)
            !BB2(inode) = Dmult(2,3)*p_DB2(ib)
            !DD1(inode) = Dmult(3,1)*p_DD1(ib)
            !DD2(inode) = Dmult(3,2)*p_DD2(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            !
            ! COMMENTED OUT, we are working solution based!
            ! FF(1+lofsp) = FF(1+lofsp) &
            !             - DD1(inode)*p_Dvector(idof) &
            !             - DD2(inode)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            ! EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,4
        p_Dvector(idofGlobal(inode)) = UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) = UU(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = UU(1+lofsp)
    
    END DO ! iel
    
    ! The final vector is a mean between the old and the new vector.
    CALL lsysbl_vectorLinearComb (rtempVector,rvector,1.0_DP-domega,domega)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ1TQ0simpleConf (rvanca, rvector, rrhs, domega,IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DSPQ1TQ0simpleConf is the same as vanca_2DSPQ1TQ0simple except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANCA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
  
  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8
    REAL(DP) :: daux
    
    ! Local arrays for informations about one element
    REAL(DP), DIMENSION(4) :: AA,BB1,BB2,DD1,DD2
    REAL(DP), DIMENSION(9) :: FF,UU
    INTEGER(PREC_VECIDX), DIMENSION(4) :: idofGlobal
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF's    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      DO inode=1,4
      
        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the four velocity unknowns and the         
        ! pressure unknown on the current element; these five unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity      
        ! components and B being an (2*4) x 1 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |        |             |        |        |                
        !   --|---X----|--           |--------|--------|                
        !                                                               
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
          
            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            
            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,4
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ************************************************************************

!<subroutine>

  PURE SUBROUTINE vanca_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)
  
!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two 
  ! velocity vectors and one pressure vector. It's used as auxiliary 
  ! routine in the simple VANCA solver to calculate an update vector
  ! for velocity/pressure.
!</description>

!<input>
  ! Diagonal elements of the local system matrix.
  REAL(DP), DIMENSION(4), INTENT(IN) :: AA
  
  ! Entries in the submatrix B1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB1

  ! Entries in the submatrix B2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB2
  
  ! Entries in the submatrix D1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD1

  ! Entries in the submatrix D2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD2

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(IN) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(OUT) :: UU
!</output>

!</subroutine>

    ! local variables

    INTEGER :: inode
    REAL(DP) :: PP,dpres
    REAL(DP), DIMENSION(9) :: AI,dff
    
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA(1)                                                   BB1(1) )
    !     (        AA(2)                                            BB1(2) )
    !     (               AA(3)                                     BB1(3) )
    !     (                      AA(4)                              BB1(4) )
    !     (                             AA(1)                       BB2(1) )
    !     (                                    AA(2)                BB2(2) )
    !     (                                           AA(3)         BB2(3) )
    !     (                                                  AA(4)  BB2(4) )
    !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the 
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial! 
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity 
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    DO inode=1,4
    
      ! Quick check if everything is ok - we don't want to divide by 0.
      IF (AA(inode)*AA(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      ! AI(.) saves the diagonal matrix S^-1:
    
      AI(inode)=1E0_DP/AA(inode)
        
    END DO

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 4x1 matrix 
    ! b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = 0.0_DP
    dff = FF
      
    DO inode = 1,4
      dpres        = dpres &
                   - AI(inode)*(DD1(inode)*BB1(inode)+DD2(inode)*BB2(inode))
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI(inode)*(DD1(inode)*dff(inode)+DD2(inode)*dff(inode+lofsv))
    END DO

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF's are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    IF (dpres*dpres .LT. 1E-20_DP)  THEN
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      RETURN
    ENDIF
      
    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:
      
    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP
      
    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:
      
    DO inode=1,4
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
    END DO

  END SUBROUTINE

  ! ************************************************************************

!<subroutine>

  PURE SUBROUTINE vanca_getcorr_2DSPQ1TQ0simple2 (UU,FF,AA1,AA2,BB1,BB2,DD1,DD2,di)
  
!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two 
  ! velocity vectors and one pressure vector. It's used as auxiliary 
  ! routine in the simple VANCA solver to calculate an update vector
  ! for velocity/pressure.
  !
  ! In contrast to vanca_getcorr_2DSPQ1TQ0simple, the two diagonal blocks
  ! in the matrix may be different from each other.
!</description>

!<input>
  ! Diagonal elements of the local system matrix A11
  REAL(DP), DIMENSION(4), INTENT(IN) :: AA1

  ! Diagonal elements of the local system matrix A22
  REAL(DP), DIMENSION(4), INTENT(IN) :: AA2
  
  ! Entries in the submatrix B1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB1

  ! Entries in the submatrix B2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB2
  
  ! Entries in the submatrix D1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD1

  ! Entries in the submatrix D2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD2
  
  ! Entry in the submatrix I (usually =0).
  ! This is the matrix in the diagonal block of the pressure, which is usually
  ! zero in saddle point problems.
  REAL(DP), INTENT(IN) :: di
  
  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(IN) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(OUT) :: UU
!</output>

!</subroutine>

    ! local variables

    INTEGER :: inode
    REAL(DP) :: PP,dpres
    REAL(DP), DIMENSION(9) :: AI1,AI2,dff
    
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA1(1)                                                  BB1(1) )
    !     (        AA1(2)                                           BB1(2) )
    !     (               AA1(3)                                    BB1(3) )
    !     (                      AA1(4)                             BB1(4) )
    !     (                             AA2(1)                      BB2(1) )
    !     (                                    AA2(2)               BB2(2) )
    !     (                                           AA2(3)        BB2(3) )
    !     (                                                  AA2(4) BB2(4) )
    !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4) di1    )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T i )
    !      ( D1  D2  I  )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the 
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t I ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !       D^t u + Ip = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !      D^t S^-1 (f-Bp) + Ip = g
    !
    !      <=>  Ip -D^t S^-1 B p  =  g - D^t S^-1 f
    !              ***********       **************
    !                 =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial! 
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity 
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    DO inode=1,4
    
      ! Quick check if everything is ok - we don't want to divide by 0.
      IF (AA1(inode)*AA1(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      IF (AA2(inode)*AA2(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      ! AI(.) saves the diagonal matrix S^-1:
    
      AI1(inode)=1E0_DP/AA1(inode)
      AI2(inode)=1E0_DP/AA2(inode)
        
    END DO

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !      ( I - D^t S^-1 B ) p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 4x1 matrix 
    ! b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           = (I - D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = di
    dff = FF
      
    DO inode = 1,4
      dpres        = dpres &
                   - AI1(inode)*DD1(inode)*BB1(inode) &
                   - AI2(inode)*DD2(inode)*BB2(inode)
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI1(inode)*DD1(inode)*dff(inode) &
                   - AI2(inode)*DD2(inode)*dff(inode+lofsv)
    END DO

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF's are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    IF (dpres*dpres .LT. 1E-20_DP)  THEN
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      RETURN
    ENDIF
      
    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:
      
    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP
      
    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:
      
    DO inode=1,4
      UU(inode)       = AI1(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI2(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D Navier-Stokes VANCA, simple diagonal version for fully coupled
  ! velocity matrix.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1 )
  !    ( A21  A22  B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ1TQ0simpleCoupConf (rvanca, rvector, rrhs, domega,IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DSPQ1TQ0simpleCoupConf is the same as vanca_2DSPQ1TQ0simpleConf,
  ! but supports fully coupled velocity submatrices.
  ! The matrices A11 and A22 must have the same structure. The matrices A12
  ! and A21 must also have the same structure. The structure of A11 and A12
  ! may be different from each other.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
  
  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12,p_KdiagonalA12
    REAL(DP), DIMENSION(:), POINTER             :: p_DA,p_DA12,p_DA21,p_DA22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8
    REAL(DP) :: daux
    
    ! Local arrays for informations about one element
    REAL(DP), DIMENSION(2*4) :: AA,BB1,BB2,DD1,DD2
    REAL(DP), DIMENSION(9) :: FF,UU
    INTEGER(PREC_VECIDX), DIMENSION(4) :: idofGlobal
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    ! Structure of A11 is assumed to be the same as A22
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_DA22 => rvanca%p_DA22

    ! Structure of A12 is assumed to be the same as A21
    p_KcolA12 => rvanca%p_KcolA12
    p_KldA12 => rvanca%p_KldA12
    p_KdiagonalA12 => rvanca%p_KdiagonalA12
    p_DA12 => rvanca%p_DA12
    p_DA21 => rvanca%p_DA21
    
    ! Structure of B1 is assumed to be the same as B2
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF's    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      DO inode=1,4
      
        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+4) = p_DA22(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the four velocity unknowns and the         
        ! pressure unknown on the current element; these five unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity      
        ! components and B being an (2*4) x 1 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |        |             |        |        |                
        !   --|---X----|--           |--------|--------|                
        !                                                               
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)
        END DO
        
        ! Tackle 'offdiagonal' matrices A12 and A21
        
        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
          
            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            
            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_2DSPQ1TQ0CPsimple (UU,FF,AA,BB1,BB2,DD1,DD2)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,4
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ************************************************************************

!<subroutine>

  PURE SUBROUTINE vanca_getcorr_2DSPQ1TQ0CPsimple (UU,FF,AA,BB1,BB2,DD1,DD2)
  
!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two 
  ! velocity vectors and one pressure vector. It's used as auxiliary 
  ! routine in the simple VANCA solver to calculate an update vector
  ! for velocity/pressure for system where the velocity is fully coupled.
!</description>

!<input>
  ! Diagonal elements of the local system matrix. (Blocks A11 and A22)
  REAL(DP), DIMENSION(2*4), INTENT(IN) :: AA
  
  ! Entries in the submatrix B1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB1

  ! Entries in the submatrix B2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB2
  
  ! Entries in the submatrix D1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD1

  ! Entries in the submatrix D2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD2

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(IN) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(OUT) :: UU
!</output>

!</subroutine>

    ! local variables

    INTEGER :: inode
    REAL(DP) :: PP,dpres
    REAL(DP), DIMENSION(9) :: AI,dff
    
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA(1)                                                   BB1(1) )
    !     (        AA(2)                                            BB1(2) )
    !     (               AA(3)                                     BB1(3) )
    !     (                      AA(4)                              BB1(4) )
    !     (                             AA(5)                       BB2(1) )
    !     (                                    AA(6)                BB2(2) )
    !     (                                           AA(7)         BB2(3) )
    !     (                                                  AA(8)  BB2(4) )
    !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the 
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial! 
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity 
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    DO inode=1,8
    
      ! Quick check if everything is ok - we don't want to divide by 0.
      IF (AA(inode)*AA(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      ! AI(.) saves the diagonal matrix S^-1:
    
      AI(inode)=1E0_DP/AA(inode)
        
    END DO

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 4x1 matrix 
    ! b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = 0.0_DP
    dff = FF
      
    DO inode = 1,4
      dpres        = dpres &
                   - AI(inode)  *(DD1(inode)*BB1(inode)) &
                   - AI(inode+4)*(DD2(inode)*BB2(inode))
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI(inode)  *(DD1(inode)*dff(inode)) &
                   - AI(inode+4)*(DD2(inode)*dff(inode+lofsv))
    END DO

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF's are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    IF (dpres*dpres .LT. 1E-20_DP)  THEN
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      RETURN
    ENDIF
      
    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:
      
    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP
      
    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:
      
    DO inode=1,4
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode+4)*(dff(inode+lofsv)-BB2(inode)*PP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D Navier-Stokes VANCA, 'full' version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ1TQ0fullConf (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DSPQ2QP1fullConf is the same as vanca_2DSPQ2QP1full except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANCA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 4      ! Q1T = 4 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 1 ! QQ0 = 1 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp,j
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,k
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new     
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new     
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new       new    
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new    
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          4    P    2      
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = edge.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

      ! Loop over all U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 8 x 2*(NMT) matrix for the two velocity      
        ! components and B~ being an (2*4) x 2 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! 2*1 pressure elements on the adjacent cells, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 2 submatrices (originally, every velocity couples with
        ! the pressure DOF's on that cell, so we have       
        ! 1 column in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              EXIT
            END IF
          END DO          
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
            
            ! Get the entries in the B-matrices
            AA(inode,      2*nnvel+1) = p_DB1(ib)
            AA(inode+nnvel,2*nnvel+1) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(2*nnvel+1,inode)       = p_DD1(ib)
            AA(2*nnvel+1,inode+nnvel) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(2*nnvel+1,inode)*p_Dvector(idof) &
                        - AA(2*nnvel+1,inode+nnvel)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)  ..............                                   :::::: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(4,4)                               :::::: )
      !     (                            AA( 5, 5) ..............       :::::: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA( 8, 8)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
      
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!
        
        DO inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
      ELSE IF (ilapackInfo .LT. 0) THEN

        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DSPQ1TQ0fullConf')

      END IF
      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D Navier-Stokes VANCA, 'full' version for fully coupled systems.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1 )
  !    ( A21  A22  B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ1TQ0fullCoupConf (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DSPQ1TQ0fullCoupConf is the same as vanca_2DSPQ1TQ0fullConf,
  ! but supports fully coupled velocity submatrices.
  ! The matrices A11 and A22 must have the same structure. The matrices A12
  ! and A21 must also have the same structure. The structure of A11 and A12
  ! may be different from each other.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA,p_DA12,p_DA21,p_DA22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12,p_KdiagonalA12
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 4      ! Q1T = 4 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 1 ! QQ0 = 1 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp,j
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,k
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    ! Structure of A11 is assumed to be the same as A22
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_DA22 => rvanca%p_DA22

    ! Structure of A12 is assumed to be the same as A21
    p_KcolA12 => rvanca%p_KcolA12
    p_KldA12 => rvanca%p_KldA12
    p_KdiagonalA12 => rvanca%p_KdiagonalA12
    p_DA12 => rvanca%p_DA12
    p_DA21 => rvanca%p_DA21
    
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new     
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new     
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new       new    
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new    
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          4    P    2      
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = edge.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

      ! Loop over all U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 8 x 2*(NMT) matrix for the two velocity      
        ! components and B~ being an (2*4) x 2 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! 2*1 pressure elements on the adjacent cells, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with four square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 2 submatrices (originally, every velocity couples with
        ! the pressure DOF's on that cell, so we have       
        ! 1 column in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode,k) = p_DA(ia)
              AA (inode+nnvel,k+nnvel) = p_DA22(ia)
              EXIT
            END IF
          END DO          
        END DO

        ! Handle the 'off-diagonal' matrices A12 and A21
        
        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode,k+nnvel) = p_DA12(ia)
              AA (inode+nnvel,k) = p_DA21(ia)
              EXIT
            END IF
          END DO          
        END DO
                
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
            
            ! Get the entries in the B-matrices
            AA(inode,      2*nnvel+1) = p_DB1(ib)
            AA(inode+nnvel,2*nnvel+1) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(2*nnvel+1,inode)       = p_DD1(ib)
            AA(2*nnvel+1,inode+nnvel) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(2*nnvel+1,inode)*p_Dvector(idof) &
                        - AA(2*nnvel+1,inode+nnvel)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)  ..............    AA( 1, 5) ..............       :::::: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B1): )
      !     (    ................ AA(4,4)    ............... AA( 4, 8)  :::::: )
      !     ( AA(5,1)  ..............    AA( 5, 5) ..............       :::::: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B2): )
      !     (    ................ AA(8,4)    ............... AA( 8, 8)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
      
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!
        
        DO inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
                                 
      ELSE IF (ilapackInfo .LT. 0) THEN
      
        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DSPQ1TQ0fullCoupConf')
        
      END IF

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D Navier-Stokes VANCA, simple diagonal and full version.
  ! Supports only Q2/QP1 discretisation.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ2QP1simple (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_VERTEXIDX)   :: NVT
    INTEGER(PREC_EDGEIDX)    :: NMT
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 9      ! Q2 = 9 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 3 ! QP1 = 3 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j,isubdof
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    DO iel=1,NEL
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          8    9    6      
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Put on AA(.) the diagonal entry of matrix A -- the 1st and the
        ! 2nd block
        AA(inode,inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+nnvel,inode+nnvel) = p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity      
        ! components and B~ being an (2*9) x 12 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most  
        ! 4*3 pressure elements on the neighbour cell, so we have       
        ! 12 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        ! or
        !
        !              IEL   
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to  
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF's on that cell, so we have       
        ! 3 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of the original B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF's), 6 (if the velocity DOF is an edge with 
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on 
        ! the boundary and there is no neighbour, or if it's the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.
        
        DO ib = ib1,ib2
        
          IF (p_KcolB(ib) .EQ. IEL) THEN
            isubdof = 1
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL) THEN
            isubdof = 2
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL*2) THEN
            isubdof = 3
          ELSE
            ! Cycle the loop - the entry belongs to another 
            ! element, not to the current one
            CYCLE
          END IF
          
          J = p_KcolB(ib)
          
          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)                                                   :::::: )
      !     (          ..                                               :AA :: )
      !     (               ..                                          :(B1): )
      !     (                     AA(9,9)                               :::::: )
      !     (                            AA(10,10)                      :::::: )
      !     (                                      ..                   :AA :: )
      !     (                                           ..              :(B2): )
      !     (                                                AA(18,18)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
      
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!
        
        DO inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
      ELSE IF (ilapackInfo .LT. 0) THEN
        
        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DSPQ2QP1simple')
        
      END IF

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ2QP1simpleConf (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DSPQ2QP1simpleConf is the same as vanca_2DSPQ2QP1simple except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANCA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_VERTEXIDX)   :: NVT
    INTEGER(PREC_EDGEIDX)    :: NMT
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 9      ! Q2 = 9 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 3 ! QP1 = 3 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j,isubdof
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          8    9    6      
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Put on AA(.) the diagonal entry of matrix A -- the 1st and the
        ! 2nd block
        AA(inode,inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+nnvel,inode+nnvel) = p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity      
        ! components and B~ being an (2*9) x 12 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most  
        ! 4*3 pressure elements on the neighbour cell, so we have       
        ! 12 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        ! or
        !
        !              IEL   
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to  
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF's on that cell, so we have       
        ! 3 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of the original B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF's), 6 (if the velocity DOF is an edge with 
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on 
        ! the boundary and there is no neighbour, or if it's the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.
        
        DO ib = ib1,ib2
        
          IF (p_KcolB(ib) .EQ. IEL) THEN
            isubdof = 1
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL) THEN
            isubdof = 2
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL*2) THEN
            isubdof = 3
          ELSE
            ! Cycle the loop - the entry belongs to another 
            ! element, not to the current one
            CYCLE
          END IF
          
          J = p_KcolB(ib)
          
          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)                                                   :::::: )
      !     (          ..                                               :AA :: )
      !     (               ..                                          :(B1): )
      !     (                     AA(9,9)                               :::::: )
      !     (                            AA(10,10)                      :::::: )
      !     (                                      ..                   :AA :: )
      !     (                                           ..              :(B2): )
      !     (                                                AA(18,18)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
      
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!
        
        DO inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
                                        
      ELSE IF (ilapackInfo .LT. 0) THEN
        
        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DSPQ2QP1simpleConf')
        
      END IF

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ2QP1full (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_VERTEXIDX)   :: NVT
    INTEGER(PREC_EDGEIDX)    :: NMT
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 9      ! Q2 = 9 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 3 ! QP1 = 3 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp,j
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,isubdof,k
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    DO iel=1,NEL
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          8    9    6      
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity      
        ! components and B~ being an (2*9) x 12 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most  
        ! 4*3 pressure elements on the adjacent cells, so we have       
        ! 12 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        ! or
        !
        !              IEL   
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to  
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF's on that cell, so we have       
        ! 3 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              EXIT
            END IF
          END DO          
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF's), 6 (if the velocity DOF is an edge with 
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on 
        ! the boundary and there is no neighbour, or if it's the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.
        
        DO ib = ib1,ib2
        
          IF (p_KcolB(ib) .EQ. IEL) THEN
            isubdof = 1
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL) THEN
            isubdof = 2
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL*2) THEN
            isubdof = 3
          ELSE
            ! Cycle the loop - the entry belongs to another 
            ! element, not to the current one
            CYCLE
          END IF
          
          J = p_KcolB(ib)
          
          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)  ..............                                   :::::: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(9,9)                               :::::: )
      !     (                            AA(10,10) ..............       :::::: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA(18,18)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      IF (ilapackInfo .EQ. 0) THEN
      
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!
        
        DO inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
      ELSE IF (ilapackInfo .LT. 0) THEN
        
        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DSPQ2QP1full')
        
      END IF

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ2QP1fullConf (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DSPQ2QP1fullConf is the same as vanca_2DSPQ2QP1full except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANCA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_VERTEXIDX)   :: NVT
    INTEGER(PREC_EDGEIDX)    :: NMT
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 9      ! Q2 = 9 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 3 ! QP1 = 3 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp,j
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,isubdof,k
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          8    9    6      
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity      
        ! components and B~ being an (2*9) x 12 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most  
        ! 4*3 pressure elements on the adjacent cells, so we have       
        ! 12 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        ! or
        !
        !              IEL   
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to  
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF's on that cell, so we have       
        ! 3 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              EXIT
            END IF
          END DO          
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF's), 6 (if the velocity DOF is an edge with 
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on 
        ! the boundary and there is no neighbour, or if it's the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.
        
        DO ib = ib1,ib2
        
          IF (p_KcolB(ib) .EQ. IEL) THEN
            isubdof = 1
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL) THEN
            isubdof = 2
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL*2) THEN
            isubdof = 3
          ELSE
            ! Cycle the loop - the entry belongs to another 
            ! element, not to the current one
            CYCLE
          END IF
          
          J = p_KcolB(ib)
          
          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)  ..............                                   :::::: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(9,9)                               :::::: )
      !     (                            AA(10,10) ..............       :::::: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA(18,18)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
      
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!
        
        DO inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
      ELSE IF (ilapackInfo .LT. 0) THEN
        
        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DSPQ2QP1fullConf')
        
      END IF

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D VANCA, 'full' version for fully coupled Navier-Stokes.
  ! Extended version with support for diagonal matrices in the pressure and
  ! arbitrary scaled submatrices.
  ! Supports only Q1~/Q0.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1  )
  !    ( A21  A22  B2  )
  !    ( D1^T D2^T I1  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  !
  ! I1 is a diagonal matrix in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! don't exist, what results in a '0' block in these positions.
  ! ***************************************************************************

!<subroutine>
             
  SUBROUTINE vanca_2DNSQ1TQ0fullCoupConfExt (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Navier Stokes optimal control problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DNSSQ1TQ0fullCoupConf supports fully coupled velocity submatrices.
  ! The matrices A11, A22 must have the same structure. 
  ! The matrices A12 and A21 must have the same structure. 
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA11
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA11
    REAL(DP), DIMENSION(:), POINTER             :: p_DA11,p_DA12,p_DA21,p_DA22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    REAL(DP), DIMENSION(:), POINTER             :: p_Da33
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA33
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 4      ! Q1T = 4 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 1 ! QQ0 = 1 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel + nnpressure
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! Offsets of the 'local' solution parts in the 'local' solution vector
    INTEGER, PARAMETER :: lofsu = 0
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! Offset information in arrays.
    INTEGER(PREC_VECIDX)     :: ioffsetu,ioffsetv,ioffsetp,j
    
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,k
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    ! Structure of A11 is assumed to be the same as A22
    p_KcolA11 => rvanca%p_KcolA
    p_KldA11 => rvanca%p_KldA
    p_DA11 => rvanca%p_DA
    p_DA22 => rvanca%p_DA22
    IF (.NOT. ASSOCIATED(p_DA22)) p_DA22 => p_DA11

    ! Structure of A12 is assumed to be the same as A21.
    ! Get A12 and A21 -- except for if the multipliers are =0, then
    ! we switch them off by nullifying the pointers.
    IF (rvanca%Dmultipliers(1,2) .NE. 0.0_DP) THEN
      p_KcolA12 => rvanca%p_KcolA12
      p_KldA12 => rvanca%p_KldA12
      p_DA12 => rvanca%p_DA12
      p_DA21 => rvanca%p_DA21
    ELSE
      NULLIFY(p_KcolA12)
      NULLIFY(p_KldA12) 
      NULLIFY(p_DA12 )
      NULLIFY(p_DA21 )
    END IF
    
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Diagonal submatrices A33 and A66 (if they exist)
    IF (rvanca%Dmultipliers(3,3) .NE. 0.0_DP) THEN
      p_Da33 => rvanca%p_DA33
      p_KdiagonalA33 => rvanca%p_KdiagonalA33
    ELSE
      NULLIFY(p_Da33)
      NULLIFY(p_KdiagonalA33)
    END IF
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new     
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new     
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new       new    
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new    
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          4    P    2      
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !                                               
      !
      ! Fetch the pressure P on the current element into FF.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      ! Get the pressure
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = edge.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

      ! Loop over all U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode.
        
        ! Primal equation
        FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u)  = (f  )                                        
        !    [ B^t 0 ] (p)    (g  )                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!   
        !
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in rectangular
        ! systems of the form                                            
        !                                                               
        !    [ === A^ === B~  ] (| ) = (f1 )                                
        !    [ B~^t       I1~ ] (u )   (f2 )                                
        !                       (| )   (g  )                                   
        !                       (p )   
        !                                     
        !                                                               
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! 2*1 pressure elements on the adjacent cells, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL.
        !                                                               
        !  (d1 ) = (f1 ) -  [ === A^ === B~  ] (| )                                 
        !  (d2 )   (f2 )    [ B~^t       I1~ ] (u )                                 
        !  (dg )   (g  )                       (| )                                    
        !                                      (p ) 
        !
        ! Extract those entries in the A-, B- and M-matrices to our local
        ! matrix AA, which belong to DOF's in our current solution vector.
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA11(idof)
        ia2 = p_KldA11(idof+1)-1
        DO ia = ia1,ia2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - ( A11  .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .  A22  .   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   ) ( p  )

          J = p_KcolA11(ia)
          
          ! Primal equation:
          FF(inode+lofsu) = FF(inode+lofsu) &
                          - rvanca%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
          FF(inode+lofsv) = FF(inode+lofsv) &
                          - rvanca%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanca%Dmultipliers(1,1)
              AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanca%Dmultipliers(2,2)
              EXIT
            END IF
          END DO          
        END DO

        ! Handle the 'off-diagonal' matrices A12 and A21
        
        IF (ASSOCIATED(p_KldA12)) THEN
          ia1 = p_KldA12(idof)
          ia2 = p_KldA12(idof+1)-1
          DO ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .  A12  .  ) ( u  )
            !   ( dv  )   ( dv  )   ( A21  .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .  ) ( p  )

            J = p_KcolA12(ia)
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanca%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanca%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)
            
            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanca%Dmultipliers(1,2)
                AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanca%Dmultipliers(2,1)
                EXIT
              END IF
            END DO          
          END DO
        END IF
        
        ! Process A33 if it exists
                
        IF (ASSOCIATED(p_KdiagonalA33)) THEN

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   I1  ) ( p  )
          !
          ! IEL is the pressure DOF which we have to tackle.
          
          daux = rvanca%Dmultipliers(3,3)
          FF(1+lofsp) = FF(1+lofsp) &
                      - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
          AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

        END IF
                
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .  B1   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .  B2   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   ) ( p  )

          J = p_KcolB(ib)
          
          daux = p_Dvector(j+ioffsetp) 
          FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanca%Dmultipliers(1,3)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanca%Dmultipliers(2,3)

          ! Don't incorporate the B-matrices into AA yet; this will come later!
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
        
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   ) ( v  )
          !   ( dp  )   ( dp  )   ( D1  D2   .   ) ( p  )
          !
          ! In AA, we simultaneously set up (locally):
          !
          !   (  .   .  B1   ) 
          !   (  .   .  B2   ) 
          !   ( D1  D2   .   ) 

          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
            
            ! Get the entries in the B-matrices.
            ! Primal equation
            AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanca%Dmultipliers(1,3)
            AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanca%Dmultipliers(2,3)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanca%Dmultipliers(3,1)
            AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanca%Dmultipliers(3,2)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                        - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! For C, we use our local AA, i.e. applying C^{-1} means to
      ! solve the local system AA dd = FF for dd. The local defect dd is then
      ! added back to the global solution vector.
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
        
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y
        
        DO inode=1,nnvel
          ! Update of the primal velocity vectors
          p_Dvector(idofGlobal(inode)+ioffsetu) &
            = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
      
      ELSE IF (ilapackInfo .LT. 0) THEN
        
        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNSSQ1TQ0fullCoupConf')
        
      END IF

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

! *****************************************************************************
! Problem class: VANCA for STEADY OPTIMAL CONTROL PROBLEMS
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_init2DNavierStokesOptC (rmatrix,rvanca)
  
!<description>
  ! Initialises the VANCA variant for 2D Navier-Stokes problems 
  ! for conformal discretisations.
  ! Checks if the "2D-Navier-Stokes" VANCA variant 
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! The substructure rvanca%rvanca2DNavSt is intitialised according
  ! to the information provided in rmatrix.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: rmatrix
  !</input>

!<inputoutput>
  ! t_vancaPointer2DSPNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vanca), INTENT(INOUT) :: rvanca
!</inputoutput>

!</subroutine>

    INTEGER :: i,j
    TYPE(t_blockDiscretisation), POINTER :: p_rblockDiscr
    
    ! Matrix must be 6x6.
    IF (rmatrix%ndiagBlocks .NE. 6) THEN
      CALL output_line ('System matrix is not 6x6.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF
    
    ! A(1:2,1:3), A(4:5,4:6) must not be virtually transposed and of format 9.
    ! A(3,:),A(6,:) must be (virtually) transposed. All matrices must be double precision.
    DO i=1,6
      DO j=1,6
      
        IF (lsysbl_isSubmatrixPresent (rmatrix,i,j)) THEN
        
          IF ( ((i .GE. 1) .AND. (i .LE. 2)) .OR. &
               ((i .GE. 4) .AND. (i .LE. 5)) ) THEN
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .NE. 0) THEN
              CALL output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
              CALL sys_halt()
            END IF
          ELSE
            IF ((i .NE. j) .AND. &
                (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                 .EQ. 0)) THEN
              CALL output_line ('B1/B2 submatrices must be virtually',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
              CALL output_line ('transposed (LSYSSC_MSPEC_TRANSPOSED)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
              CALL sys_halt()
            END IF
          END IF
          
          IF ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX7) .AND. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
            CALL output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
            CALL sys_halt()
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            CALL output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
            CALL sys_halt()
          END IF

        END IF ! neq != 0
      END DO
    END DO
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    IF ((rmatrix%RmatrixBlock(1,3)%NA .NE. rmatrix%RmatrixBlock(3,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .NE. rmatrix%RmatrixBlock(3,1)%NCOLS)) THEN
      CALL output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF

    IF ((rmatrix%RmatrixBlock(2,3)%NA .NE. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .NE. rmatrix%RmatrixBlock(3,2)%NCOLS)) THEN
      CALL output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF      
  
    ! The structure of A(4,6) must be identical to A(6,4) and
    ! that of A(5,6) must be identical to A(6,5).
    IF ((rmatrix%RmatrixBlock(4,6)%NA .NE. rmatrix%RmatrixBlock(6,4)%NA) .OR. &
        (rmatrix%RmatrixBlock(4,6)%NEQ .NE. rmatrix%RmatrixBlock(6,4)%NCOLS)) THEN
      CALL output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF

    ! Primal and dual B- and D-matrices must share the same data.
    ! The matrices may be switched off by dscaleFactor=0, but they must have the
    ! same data arrays! Just for 'efficiency'...
    IF ((rmatrix%RmatrixBlock(5,6)%NA .NE. rmatrix%RmatrixBlock(6,5)%NA) .OR. &
        (rmatrix%RmatrixBlock(5,6)%NEQ .NE. rmatrix%RmatrixBlock(6,5)%NCOLS)) THEN
      CALL output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF      
    
    IF (.NOT. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(1,3),rmatrix%RmatrixBlock(4,6))) THEN
      CALL output_line ('Content of primal and dual B1-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF

    IF (.NOT. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(2,3),rmatrix%RmatrixBlock(5,6))) THEN
      CALL output_line ('Content of primal and dual B2-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF
  
    IF (.NOT. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(3,1),rmatrix%RmatrixBlock(6,4))) THEN
      CALL output_line ('Content of primal and dual D1-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF

    IF (.NOT. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(3,2),rmatrix%RmatrixBlock(6,5))) THEN
      CALL output_line ('Content of primal and dual D2-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF
  
    ! Fill the output structure with data of the matrices.
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),&
        rvanca%rvanca2DNavStOptC%p_DB1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),&
        rvanca%rvanca2DNavStOptC%p_DB2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),&
        rvanca%rvanca2DNavStOptC%p_DD1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),&
        rvanca%rvanca2DNavStOptC%p_DD2)
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),&
        rvanca%rvanca2DNavStOptC%p_KcolB)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), &
        rvanca%rvanca2DNavStOptC%p_KldB )
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanca%rvanca2DNavStOptC%p_KcolA11)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanca%rvanca2DNavStOptC%p_KldA11 )
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanca%rvanca2DNavStOptC%p_KdiagonalA11)
    ELSE
      rvanca%rvanca2DNavStOptC%p_KdiagonalA11 => rvanca%rvanca2DNavStOptC%p_KldA11
    END IF
    
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanca%rvanca2DNavStOptC%p_DA11 )

    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
        rvanca%rvanca2DNavStOptC%p_DA22 )
    
    ! What is with A12 and A21? Do they exist? With a scale factor = 1.0?
    IF (lsysbl_isSubmatrixPresent(rmatrix,1,2)) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
          rvanca%rvanca2DNavStOptC%p_DA12 )
          
      CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,2),&
          rvanca%rvanca2DNavStOptC%p_KcolA12)
      CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,2), &
          rvanca%rvanca2DNavStOptC%p_KldA12 )
          
      ! Get the structure. It's assumed that A12 and A21 have the same!
      IF (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,2), &
                                rvanca%rvanca2DNavStOptC%p_KdiagonalA12)
      ELSE
        rvanca%rvanca2DNavStOptC%p_KdiagonalA12 => rvanca%rvanca2DNavStOptC%p_KldA12
      END IF
      
      IF (.NOT. lsysbl_isSubmatrixPresent(rmatrix,2,1)) THEN
        CALL output_line ('If A12 is given, A21 must also be given!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
        CALL sys_halt()
      END IF
      
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
          rvanca%rvanca2DNavStOptC%p_DA21 )
    END IF
    
    ! Now we come to the dual equation -> A(4:6,4:6).
    ! We assume that A44 and A55 has the same structure as A11 and A22.
    
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,4),&
        rvanca%rvanca2DNavStOptC%p_DA44 )

    ! What is with A55? Is it the same as A44?
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(5,5),&
        rvanca%rvanca2DNavStOptC%p_DA55 )
  
    ! What is with A12 and A21? Do they exist? With a scale factor != 0.0?
    IF (lsysbl_isSubmatrixPresent(rmatrix,4,5)) THEN

      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,5),&
          rvanca%rvanca2DNavStOptC%p_DA45 )
          
      CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(4,5),&
          rvanca%rvanca2DNavStOptC%p_KcolA45)
      CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,5), &
          rvanca%rvanca2DNavStOptC%p_KldA45 )
          
      ! Get the structure. It's assumed that A12 and A21 have the same!
      IF (rmatrix%RmatrixBlock(4,5)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(4,5), &
                                rvanca%rvanca2DNavStOptC%p_KdiagonalA45)
      ELSE
        rvanca%rvanca2DNavStOptC%p_KdiagonalA45 => rvanca%rvanca2DNavStOptC%p_KldA45
      END IF
      
      IF (.NOT. lsysbl_isSubmatrixPresent(rmatrix,5,4)) THEN
        CALL output_line ('If A45 is given, A54 must also be given!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
        CALL sys_halt()
      END IF
      
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(5,4),&
          rvanca%rvanca2DNavStOptC%p_DA54 )
    END IF    
    
    ! Is there an identity matrix present at A(3,3) and/or A(6,6)?
    IF (lsysbl_isSubmatrixPresent(rmatrix,3,3)) THEN
      ! The matrix must be of format 9.
      IF (rmatrix%RmatrixBlock(3,3)%dscaleFactor .NE. 0.0_DP) THEN
        CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
            rvanca%rvanca2DNavStOptC%p_DA33 )

        IF (rmatrix%RmatrixBlock(3,3)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
          CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(3,3), &
                                  rvanca%rvanca2DNavStOptC%p_KdiagonalA33)
        ELSE
          CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,3), &
                                  rvanca%rvanca2DNavStOptC%p_KdiagonalA33)
        END IF

      END IF
    END IF

    IF (lsysbl_isSubmatrixPresent(rmatrix,6,6)) THEN
      IF (rmatrix%RmatrixBlock(6,6)%dscaleFactor .NE. 0.0_DP) THEN
        CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(6,6),&
            rvanca%rvanca2DNavStOptC%p_DA66 )

        IF (rmatrix%RmatrixBlock(6,6)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
          CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(6,6), &
                                  rvanca%rvanca2DNavStOptC%p_KdiagonalA66)
        ELSE
          CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(6,6), &
                                  rvanca%rvanca2DNavStOptC%p_KdiagonalA66)
        END IF
      END IF
    END IF
    
    ! Get the mass matrix/matrices -- if they are present.
    ! It's assumed that all mass matrices are the same except for their
    ! multiplication factors!
    IF (lsysbl_isSubmatrixPresent (rmatrix,1,4)) THEN
      IF (rmatrix%RmatrixBlock(1,4)%cmatrixFormat .EQ. LSYSSC_MATRIXD) THEN
        CALL output_line ('Lumped mass matrices not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
        CALL sys_halt()
      END IF
      
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,4),&
          rvanca%rvanca2DNavStOptC%p_DM14 )
          
      CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,4),&
          rvanca%rvanca2DNavStOptC%p_KcolM)
      CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,4), &
          rvanca%rvanca2DNavStOptC%p_KldM )
          
      IF (rmatrix%RmatrixBlock(1,4)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,4), &
                                rvanca%rvanca2DNavStOptC%p_KdiagonalM)
      ELSE
        rvanca%rvanca2DNavStOptC%p_KdiagonalM => rvanca%rvanca2DNavStOptC%p_KldM
      END IF
      
    END IF

    ! Get the coupling matrix from the primal to the dual system. 
    ! This may be the mass matrix or a completely decoupled matrix.
    ! In all cases, the submatrices have the same structure as the mass
    ! matrix, so if the structure of the mass matrix is not yet set,
    ! we have to fetch it!
    IF (lsysbl_isSubmatrixPresent (rmatrix,4,1)) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,1),&
          rvanca%rvanca2DNavStOptC%p_DR41 )

      IF (.NOT. ASSOCIATED(rvanca%rvanca2DNavStOptC%p_KcolM)) THEN
        CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(4,1),&
            rvanca%rvanca2DNavStOptC%p_KcolM)
        CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,1), &
            rvanca%rvanca2DNavStOptC%p_KldM )
      END IF
    END IF

    IF (lsysbl_isSubmatrixPresent (rmatrix,4,2)) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,2),&
          rvanca%rvanca2DNavStOptC%p_DR42 )

      IF (.NOT. ASSOCIATED(rvanca%rvanca2DNavStOptC%p_KcolM)) THEN
        CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(4,2),&
            rvanca%rvanca2DNavStOptC%p_KcolM)
        CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,2), &
            rvanca%rvanca2DNavStOptC%p_KldM )
      END IF
    END IF

    IF (lsysbl_isSubmatrixPresent (rmatrix,5,1)) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(5,1),&
          rvanca%rvanca2DNavStOptC%p_DR51 )

      IF (.NOT. ASSOCIATED(rvanca%rvanca2DNavStOptC%p_KcolM)) THEN
        CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(5,1),&
            rvanca%rvanca2DNavStOptC%p_KcolM)
        CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(5,1), &
            rvanca%rvanca2DNavStOptC%p_KldM )
      END IF
    END IF

    IF (lsysbl_isSubmatrixPresent (rmatrix,5,2)) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(5,2),&
          rvanca%rvanca2DNavStOptC%p_DR52 )

      IF (.NOT. ASSOCIATED(rvanca%rvanca2DNavStOptC%p_KcolM)) THEN
        CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(5,2),&
            rvanca%rvanca2DNavStOptC%p_KcolM)
        CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(5,2), &
            rvanca%rvanca2DNavStOptC%p_KldM )
      END IF
    END IF

    ! Get the multiplication factors of the submatrices.
    ! (-> for a later implementation; currently, the multipliers are not used!)
    rvanca%rvanca2DNavStOptC%Dmultipliers(1:6,1:6) = &
        rmatrix%RmatrixBlock(1:6,1:6)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest
    
    IF (.NOT. ASSOCIATED(p_rblockDiscr)) THEN
      CALL output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF
    
    ! Get the discretisation structure of U,V and P from the block
    ! discretisation structure.
    ! We assume that the discretisation of the dual equations are the same
    ! as for the primal equations!
    rvanca%rvanca2DNavStOptC%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanca%rvanca2DNavStOptC%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanca%rvanca2DNavStOptC%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)
    
    IF (rvanca%rvanca2DNavStOptC%p_rspatialDiscrU%inumFESpaces .NE. &
        rvanca%rvanca2DNavStOptC%p_rspatialDiscrV%inumFESpaces) THEN
      CALL output_line (&
          'Discretisation structures of X- and Y-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF

    IF ((rvanca%rvanca2DNavStOptC%p_rspatialDiscrP%inumFESpaces .NE. 1) .AND. &
        (rvanca%rvanca2DNavStOptC%p_rspatialDiscrP%inumFESpaces .NE. &
          rvanca%rvanca2DNavStOptC%p_rspatialDiscrU%inumFESpaces)) THEN
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      CALL output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesOptC')
      CALL sys_halt()
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DNavierStokesOptC (rvanca2DNavStOptC, rvector, rrhs, domega, csubtype)
  
!<description>
  ! This routine applies the VANCA variant for 2D Navier-Stokes 
  ! optimal control problems to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANCA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector

  ! The subtype of VANCA that should handle the above problem class.
  ! One of the VANCATP_xxxx constants, e.g. VANCATP_DIAGONAL.
  INTEGER :: csubtype
  
!</input>

!<inputoutput>
  ! t_vanca structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavStOptC), INTENT(INOUT) :: rvanca2DNavStOptC
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ielementdist
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementList
    TYPE(t_elementDistribution), POINTER :: p_relementDistrU
    TYPE(t_elementDistribution), POINTER :: p_relementDistrV
    TYPE(t_elementDistribution), POINTER :: p_relementDistrP 
    
    ! 2D Navier Stokes problem.

    ! Loop through the element distributions of the velocity.
    DO ielementdist = 1,rvanca2DNavStOptC%p_rspatialDiscrU%inumFESpaces
    
      ! Get the corresponding element distributions of U, V and P.
      p_relementDistrU => &
          rvanca2DNavStOptC%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanca2DNavStOptC%p_rspatialDiscrV%RelementDistr(ielementdist)
      
      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      IF (rvanca2DNavStOptC%p_rspatialDiscrP%inumFESpaces .GT. 1) THEN
        p_relementDistrP => &
            rvanca2DNavStOptC%p_rspatialDiscrP%RelementDistr(ielementdist)
      ELSE
        p_relementDistrP => &
            rvanca2DNavStOptC%p_rspatialDiscrP%RelementDistr(1)
      END IF
      
      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y-velocity (and to that
      ! of the pressure).
      CALL storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)
      
      ! Which element combination do we have now?
      IF ((elem_getPrimaryElement(p_relementDistrU%celement) .EQ. EL_Q1T) .AND. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .EQ. EL_Q1T) .AND. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .EQ. EL_Q0)) THEN
        
        ! Q1~/Q1~/Q0 discretisation
        
        SELECT CASE (csubtype)
        
        CASE (VANCATP_FULLOPTC_PRIMAL)

          ! Apply the conformal VANCA that allows different matrices
          ! in A11, A12, A21 and A22!
          CALL vanca_2DNSSOCQ1TQ0fullCoupCfFB (rvanca2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList,1)
              
        CASE (VANCATP_FULLOPTC_DUAL)
                
          ! Apply the conformal VANCA that allows different matrices
          ! in A11, A12, A21 and A22!
          CALL vanca_2DNSSOCQ1TQ0fullCoupCfFB (rvanca2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList,2)

        CASE (VANCATP_DIAGOPTC)
        
          ! Apply the conformal diagonal VANCA that allows different 
          ! matrices in A11, A22, A33 and A44!
          CALL vanca_2DNSSOCQ1TQ0diagCoupConf (rvanca2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList)
              
        CASE DEFAULT
        
          ! Apply the conformal VANCA that allows different matrices
          ! in A11, A12, A21 and A22!
          CALL vanca_2DNSSOCQ1TQ0fullCoupConf (rvanca2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList)

        END SELECT          

      ELSE

        CALL output_line ('Unsupported discretisation.',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNavierStokesOptC')
        CALL sys_halt()
        
      END IF
      
    END DO
      
  END SUBROUTINE

  ! ***************************************************************************
  ! 2D VANCA, 'diagonal' version for fully coupled Navier-Stokes systems with
  ! primal and dual equations.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11       B1               )
  !    (      A22  B2               )
  !    ( D1^T D2^T I1               )
  !    (               A44       B1 )
  !    (                    A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! The mass matrices must all be the same except for their scaling factors!
  ! Here, a and b are arbitrary multiplication factors for the mass matrices.
  !
  ! I1 and I2 are two diagonal matrices in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! don't exist, what results in a '0' block in these positions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DNSSOCQ1TQ0diagCoupConf (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Navier Stokes optimal control problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DNSSOCQ1TQ0fullCoupConf supports fully coupled velocity submatrices.
  ! The matrices A11, A22, A44 and A55 must have the same structure. 
  ! The matrices A12 and A21 must have the same structure. 
  ! The matrices A45 and A54 must have the same structure. 
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavStOptC), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    REAL(DP) :: dmult11,dmult22,dmult44,dmult55
    REAL(DP) :: dmultb1,dmultb2,dmultb3,dmultb4
    REAL(DP) :: dmultd1,dmultd2,dmultd3,dmultd4
    REAL(DP) :: dmult33,dmult66
    REAL(DP) :: di1,di2
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA11,p_Da22,p_Da44,p_Da55
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1,p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DB3,p_DB4
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1,p_DD2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD3,p_DD4
    REAL(DP), DIMENSION(:), POINTER             :: p_DA33,p_DA66
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA33,p_KdiagonalA66
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetu,ioffsetv,ioffsetp
    INTEGER(PREC_VECIDX)     :: ioffsetl1,ioffsetl2,ioffsetxi
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8
    REAL(DP) :: daux1,daux2,daux3,daux4
    
    ! Local arrays for informations about one element -- for primal and dual space.
    REAL(DP), DIMENSION(4) :: AA11,AA22,BB1,BB2,DD1,DD2
    REAL(DP), DIMENSION(9) :: FFp,UUp
    REAL(DP), DIMENSION(4) :: AA33,AA44,BB3,BB4,DD3,DD4
    REAL(DP), DIMENSION(9) :: FFd,UUd
    INTEGER(PREC_VECIDX), DIMENSION(4) :: idofGlobal
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA11
    p_KldA => rvanca%p_KldA11
    p_KdiagonalA => rvanca%p_KdiagonalA11
    p_DA11 => rvanca%p_DA11
    p_DA22 => rvanca%p_DA22
    p_Da44 => rvanca%p_Da44
    p_Da55 => rvanca%p_Da55
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DB3 => rvanca%p_DB1
    p_DB4 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    p_DD3 => rvanca%p_DD1
    p_DD4 => rvanca%p_DD2
    p_DA33 => rvanca%p_DA33
    p_DA66 => rvanca%p_DA66
    p_KdiagonalA33 => rvanca%p_KdiagonalA33
    p_KdiagonalA66 => rvanca%p_KdiagonalA66
    
    dmult11 = rvanca%Dmultipliers(1,1)
    dmult22 = rvanca%Dmultipliers(2,2)
    dmult44 = rvanca%Dmultipliers(4,4)
    dmult55 = rvanca%Dmultipliers(5,5)
    
    dmultb1 = rvanca%Dmultipliers(1,3)
    dmultb2 = rvanca%Dmultipliers(2,3)
    dmultb3 = rvanca%Dmultipliers(4,6)
    dmultb4 = rvanca%Dmultipliers(5,6)
    
    dmultd1 = rvanca%Dmultipliers(3,1)
    dmultd2 = rvanca%Dmultipliers(3,2)
    dmultd3 = rvanca%Dmultipliers(6,4)
    dmultd4 = rvanca%Dmultipliers(6,5)
    
    dmult33 = rvanca%Dmultipliers(3,3)
    dmult66 = rvanca%Dmultipliers(6,6)
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the solution components
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetl1 = ioffsetp+rvector%RvectorBlock(3)%NEQ
    ioffsetl2 = ioffsetl1+rvector%RvectorBlock(4)%NEQ
    ioffsetxi = ioffsetl2+rvector%RvectorBlock(5)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF's    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP
    
      FFp(1+lofsp) = p_Drhs(iel+ioffsetp)
      FFd(1+lofsp) = p_Drhs(iel+ioffsetxi)

      ! Loop over all 4 U-nodes of that element.
      DO inode=1,4
      
        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA11(inode) = dmult11*p_DA11(p_KdiagonalA(idof))
        AA22(inode) = dmult22*p_DA22(p_KdiagonalA(idof))
        AA33(inode) = dmult44*p_Da44(p_KdiagonalA(idof))
        AA44(inode) = dmult55*p_Da55(p_KdiagonalA(idof))
        
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FFp(inode)       = p_Drhs(idof+ioffsetu)
        FFp(inode+lofsv) = p_Drhs(idof+ioffsetv)

        FFd(inode)       = p_Drhs(idof+ioffsetl1)
        FFd(inode+lofsv) = p_Drhs(idof+ioffsetl2)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the four velocity unknowns and the         
        ! pressure unknown on the current element; these five unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity      
        ! components and B being an (2*4) x 1 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |        |             |        |        |                
        !   --|---X----|--           |--------|--------|                
        !                                                               
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux1 = dmult11*p_DA11(ia)
          daux2 = dmult22*p_DA22(ia)
          daux3 = dmult44*p_Da44(ia)
          daux4 = dmult55*p_Da55(ia)

          FFp(inode)       = FFp(inode)      -daux1*p_Dvector(J+ioffsetu)
          FFp(inode+lofsv) = FFp(inode+lofsv)-daux2*p_Dvector(J+ioffsetv)

          FFd(inode)       = FFd(inode)      -daux3*p_Dvector(J+ioffsetl1)
          FFd(inode+lofsv) = FFd(inode+lofsv)-daux4*p_Dvector(J+ioffsetl2)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux1 = p_Dvector(j+ioffsetp)
          daux2 = p_Dvector(j+ioffsetxi)

          FFp(inode)       = FFp(inode)      -dmultb1*p_DB1(ib)*daux1
          FFp(inode+lofsv) = FFp(inode+lofsv)-dmultb2*p_DB2(ib)*daux1

          FFd(inode)       = FFd(inode)      -dmultb3*p_DB3(ib)*daux2
          FFd(inode+lofsv) = FFd(inode+lofsv)-dmultb4*p_DB4(ib)*daux2
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
          
            ! Get the entries in the B-matrices
            BB1(inode) = dmultb1*p_DB1(ib)
            BB2(inode) = dmultb2*p_DB2(ib)

            BB3(inode) = dmultb3*p_DB3(ib)
            BB4(inode) = dmultb4*p_DB4(ib)
            
            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = dmultd1*p_DD1(ib)
            DD2(inode) = dmultd2*p_DD2(ib)

            DD3(inode) = dmultd3*p_DD3(ib)
            DD4(inode) = dmultd4*p_DD4(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FFp(1+lofsp) = FFp(1+lofsp) &
                          - DD1(inode)*p_Dvector(idof+ioffsetu) &
                          - DD2(inode)*p_Dvector(idof+ioffsetv)
          
            FFd(1+lofsp) = FFd(1+lofsp) &
                         - DD3(inode)*p_Dvector(idof+ioffsetl1) &
                         - DD4(inode)*p_Dvector(idof+ioffsetl2)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
      
      ! If we have blocks at A33 or A66, get the corresponding elements.
      ! We need the IF-commands here as the arrays p_DaXX/p_KdiagonalXX may 
      ! be undefined.

      di1 = 0.0_DP
      di2 = 0.0_DP

      IF (dmult33 .NE. 0.0_DP) THEN
        di1 = dmult33*p_DA33(p_KdiagonalA33(iel))
      END IF

      IF (dmult66 .NE. 0.0_DP) THEN
        di2 = dmult66*p_DA66(p_KdiagonalA66(iel))
      END IF
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_2DSPQ1TQ0simple2 (UUp,FFp,AA11,AA22,BB1,BB2,DD1,DD2,di1)
      CALL vanca_getcorr_2DSPQ1TQ0simple2 (UUd,FFd,AA33,AA33,BB3,BB4,DD3,DD4,di2)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,4
        p_Dvector(idofGlobal(inode)+ioffsetu) &
          = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * UUp(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UUp(inode+lofsv)

        p_Dvector(idofGlobal(inode)+ioffsetl1) &
          = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * UUd(inode)
        p_Dvector(idofGlobal(inode)+ioffsetl2) &
          = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * UUd(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UUp(1+lofsp)

      p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + domega * UUd(1+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D VANCA, 'full' version for fully coupled Navier-Stokes systems with
  ! primal and dual equations.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1  aM           )
  !    ( A21  A22  B2       aM      )
  !    ( D1^T D2^T I1               )
  !    ( bR            A44  A45  B1 )
  !    (      bR       A54  A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! The mass matrices must all be the same except for their scaling factors!
  ! Here, a and b are arbitrary multiplication factors for the mass matrices.
  !
  ! I1 and I2 are two diagonal matrices in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! don't exist, what results in a '0' block in these positions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DNSSOCQ1TQ0fullCoupConf (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Navier Stokes optimal control problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DNSSOCQ1TQ0fullCoupConf supports fully coupled velocity submatrices.
  ! The matrices A11, A22, A44 and A55 must have the same structure. 
  ! The matrices A12 and A21 must have the same structure. 
  ! The matrices A45 and A54 must have the same structure. 
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavStOptC), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA11
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA11
    REAL(DP), DIMENSION(:), POINTER             :: p_DA11,p_DA12,p_DA21,p_DA22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA45
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA45
    REAL(DP), DIMENSION(:), POINTER             :: p_DA44,p_DA45,p_DA54,p_DA55
    REAL(DP), DIMENSION(:), POINTER             :: p_DR41,p_DR52,p_DR51,p_DR42
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolM
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldM
    REAL(DP), DIMENSION(:), POINTER             :: p_DM
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    REAL(DP), DIMENSION(:), POINTER             :: p_Da33,p_Da66
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA33,p_KdiagonalA66
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 4      ! Q1T = 4 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 1 ! QQ0 = 1 DOF's per pressure
    INTEGER, PARAMETER :: nndualvel = 4      ! Q1T = 4 DOF's per dual velocity
    INTEGER, PARAMETER :: nndualpressure = 1 ! QQ0 = 1 DOF's per dual pressure
    INTEGER, PARAMETER :: nnprimal = 2*nnvel+nnpressure ! Q1~/Q1~/Q0 = 4+4+1 = 9 DOF's per element
    INTEGER, PARAMETER :: nnld = 2*nnprimal
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! Offsets of the 'local' solution parts in the 'local' solution vector
    INTEGER, PARAMETER :: lofsu = 0
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    INTEGER, PARAMETER :: lofsl1 = 2*nnvel+1
    INTEGER, PARAMETER :: lofsl2 = 2*nnvel+1+nnvel
    INTEGER, PARAMETER :: lofsxi = 2*nnvel+1+2*nnvel
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! Offset information in arrays.
    ! Primal variables
    INTEGER(PREC_VECIDX)     :: ioffsetu,ioffsetv,ioffsetp,j
    
    ! Dual variables
    INTEGER(PREC_VECIDX)     :: ioffsetl1,ioffsetl2,ioffsetxi
    
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,k
    REAL(DP) :: daux,daux2
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    ! Structure of A11 is assumed to be the same as A22
    p_KcolA11 => rvanca%p_KcolA11
    p_KldA11 => rvanca%p_KldA11
    p_DA11 => rvanca%p_DA11
    p_DA22 => rvanca%p_DA22

    ! Structure of A12 is assumed to be the same as A21.
    ! Get A12 and A21 -- except for if the multipliers are =0, then
    ! we switch them off by nullifying the pointers.
    IF (rvanca%Dmultipliers(1,2) .NE. 0.0_DP) THEN
      p_KcolA12 => rvanca%p_KcolA12
      p_KldA12 => rvanca%p_KldA12
      p_DA12 => rvanca%p_DA12
      p_DA21 => rvanca%p_DA21
    ELSE
      NULLIFY(p_KcolA12)
      NULLIFY(p_KldA12) 
      NULLIFY(p_DA12 )
      NULLIFY(p_DA21 )
    END IF
    
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Structure of A44 is assumed to be the same as A55, A11 and A22
    p_DA44 => rvanca%p_DA44
    p_DA55 => rvanca%p_DA55

    ! Structure of A45 is assumed to be the same as A54
    p_KcolA45 => rvanca%p_KcolA45
    p_KldA45 => rvanca%p_KldA45
    p_DA45 => rvanca%p_DA45
    p_DA54 => rvanca%p_DA54
    
    ! Mass matrix - if it's given, otherwise the pointers will be set to NULL
    ! because of the initialisation of the structure!
    p_KcolM => rvanca%p_KcolM
    p_KldM => rvanca%p_KldM
    p_DM => rvanca%p_DM14
    
    ! Coupling matrix in the dual equation at position (4:5,1:2). For a standard
    ! system, there is A(4,1) = A(5,2) = M and A(5,1) = A(4,2) = 0.
    ! For a Newton system, this block is completely decoupled!
    p_DR41 => rvanca%p_DR41
    p_DR42 => rvanca%p_DR42
    p_DR51 => rvanca%p_DR51
    p_DR52 => rvanca%p_DR52
    
    ! Diagonal submatrices A33 and A66 (if they exist)
    IF (rvanca%Dmultipliers(3,3) .NE. 0.0_DP) THEN
      p_Da33 => rvanca%p_DA33
      p_KdiagonalA33 => rvanca%p_KdiagonalA33
    ELSE
      NULLIFY(p_Da33)
      NULLIFY(p_KdiagonalA33)
    END IF
    
    IF (rvanca%Dmultipliers(6,6) .NE. 0.0_DP) THEN
      p_Da66 => rvanca%p_DA66
      p_KdiagonalA66 => rvanca%p_KdiagonalA66
    ELSE
      NULLIFY(p_Da66)
      NULLIFY(p_KdiagonalA66)
    END IF

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    ! Get the offsets of lambda1, lambda2 and xi, so the offsets
    ! of the dual solution vectors.
    ioffsetl1 = ioffsetp+rvector%RvectorBlock(3)%NEQ
    ioffsetl2 = ioffsetl1+rvector%RvectorBlock(4)%NEQ
    ioffsetxi = ioffsetl2+rvector%RvectorBlock(5)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new     
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new     
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new       new    
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new    
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          4    P    2      
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !                                               
      !
      ! Fetch the pressure P on the current element into FF.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      ! Get the primal pressure
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      
      ! Get the dual pressure
      FF(1+lofsxi) = p_Drhs(iel+ioffsetxi)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = edge-NVT.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

      ! Loop over all U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode.
        
        ! Primal equation
        FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! dual equation
        FF(inode+lofsl1) = p_Drhs(idof+ioffsetl1)
        FF(inode+lofsl2) = p_Drhs(idof+ioffsetl2)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B  M     ] (u)  = (f  )                                        
        !    [ B^t 0        ] (p)    (g  )                                        
        !    [ M      A   B ] (l)    (fl )                                        
        !    [        B^t 0 ] (xi)   (fxi)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!   
        !
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in rectangular
        ! systems of the form                                            
        !                                                               
        !    [ === A^ === B~  === M^ ====   ] (| ) = (f1 )                                
        !    [ B~^t       I1~               ] (u )   (f2 )                                
        !    [ === M^ ===     === A^ === B~ ] (| )   (g  )                                   
        !    [                B~^t       I2~] (p )   (fl1)
        !                                     (| )   (fl2)
        !                                     (l )   (flg)
        !                                     (| )
        !                                     (xi)
        !                                     
        !                                                               
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! 2*1 pressure elements on the adjacent cells, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL.
        !                                                               
        !  (d1 ) = (f1 ) -  [ === A^ === B~  === M^ ====   ] (| )                                 
        !  (d2 )   (f2 )    [ B~^t       I1~               ] (u )                                 
        !  (dg )   (g  )    [ === M^ ===     === A^ === B~ ] (| )                                    
        !  (dl1)   (fl1)    [                B~^t       I2~] (p ) 
        !  (dl2)   (fl2)                                     (| ) 
        !  (dlg)   (flg)                                     (l ) 
        !                                                    (| )
        !                                                    (xi)
        !
        ! Extract those entries in the A-, B- and M-matrices to our local
        ! matrix AA, which belong to DOF's in our current solution vector.
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA11(idof)
        ia2 = p_KldA11(idof+1)-1
        DO ia = ia1,ia2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - ( A11  .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .  A22  .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .  A44  .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .  A55  .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

          J = p_KcolA11(ia)
          
          ! Primal equation:
          FF(inode+lofsu) = FF(inode+lofsu) &
                          - rvanca%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
          FF(inode+lofsv) = FF(inode+lofsv) &
                          - rvanca%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! dual equation
          FF(inode+lofsl1) = FF(inode+lofsl1) &
                           - rvanca%Dmultipliers(4,4)*p_DA44(ia)*p_Dvector(J+ioffsetl1)
          FF(inode+lofsl2) = FF(inode+lofsl2) &
                           - rvanca%Dmultipliers(5,5)*p_DA55(ia)*p_Dvector(J+ioffsetl2)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanca%Dmultipliers(1,1)
              AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanca%Dmultipliers(2,2)
              
              AA (inode+lofsl1,k+lofsl1) = p_DA44(ia)*rvanca%Dmultipliers(4,4)
              AA (inode+lofsl2,k+lofsl2) = p_DA55(ia)*rvanca%Dmultipliers(5,5)
              EXIT
            END IF
          END DO          
        END DO

        ! Handle the 'off-diagonal' matrices A12 and A21
        
        IF (ASSOCIATED(p_KldA12)) THEN
          ia1 = p_KldA12(idof)
          ia2 = p_KldA12(idof+1)-1
          DO ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .  A12  .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   ( A21  .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA12(ia)
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanca%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanca%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)
            
            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanca%Dmultipliers(1,2)
                AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanca%Dmultipliers(2,1)
                EXIT
              END IF
            END DO          
          END DO
        END IF
        
        ! Handle the 'off-diagonal' matrices A45 and A54 if they exist
        IF (ASSOCIATED(p_KldA45)) THEN
          ia1 = p_KldA45(idof)
          ia2 = p_KldA45(idof+1)-1
          DO ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .  A45  .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .  A54  .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA45(ia)
            FF(inode+lofsl1) = FF(inode+lofsl1) &
                             - rvanca%Dmultipliers(4,5)*p_DA45(ia)*p_Dvector(J+ioffsetl2)
            FF(inode+lofsl2) = FF(inode+lofsl2) &
                             - rvanca%Dmultipliers(5,4)*p_DA54(ia)*p_Dvector(J+ioffsetl1)
            
            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsl1,k+lofsl2) = p_DA45(ia)*rvanca%Dmultipliers(4,5)
                AA (inode+lofsl2,k+lofsl1) = p_DA54(ia)*rvanca%Dmultipliers(5,4)
                EXIT
              END IF
            END DO          
          END DO
        END IF
                
        ! Process A33 if it exists
                
        IF (ASSOCIATED(p_KdiagonalA33)) THEN

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   I1  .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
          !
          ! IEL is the pressure DOF which we have to tackle.
          
          daux = rvanca%Dmultipliers(3,3)
          FF(1+lofsp) = FF(1+lofsp) &
                      - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
          AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

        END IF
                
        ! Process A66 if it exists
                
        IF (ASSOCIATED(p_KdiagonalA66)) THEN

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   I2 ) ( xi )
          !
          ! IEL is the pressure DOF which we have to tackle.
          
          daux = rvanca%Dmultipliers(6,6)
          FF(1+lofsxi) = FF(1+lofsxi) &
                       - daux*p_DA66(p_KdiagonalA66(IEL))*p_Dvector(IEL+ioffsetxi)
          AA(1+lofsxi,1+lofsxi) = daux*p_DA66(p_KdiagonalA66(IEL))

        END IF
                
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .  B1   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .  B2   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .  B1  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .  B2  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

          J = p_KcolB(ib)
          
          ! primal equation
          daux = p_Dvector(j+ioffsetp) 
          FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanca%Dmultipliers(1,3)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanca%Dmultipliers(2,3)

          ! dual equation
          daux2 = p_Dvector(j+ioffsetxi)
          FF(inode+lofsl1) = FF(inode+lofsl1)-p_DB1(ib)*daux2 * rvanca%Dmultipliers(4,6)
          FF(inode+lofsl2) = FF(inode+lofsl2)-p_DB2(ib)*daux2 * rvanca%Dmultipliers(5,6)
          
          ! Don't incorporate the B-matrices into AA yet; this will come later!
        END DO
        
        ! The mass matrix defect.
        IF (ASSOCIATED(p_DM)) THEN
          ! We assume: multiplier of A(1,4) = multiplier of A(2,5)
          daux = rvanca%Dmultipliers(1,4)
          
          ! We assume: multiplier of A(4,1) = multiplier of A(5,2)
          daux2 = rvanca%Dmultipliers(4,1)
          
          ia1 = p_KldM(idof)
          ia2 = p_KldM(idof+1)-1
          DO ia = ia1,ia2

            J = p_KcolM(ia)

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .  aM   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .  aM   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            FF(inode+lofsu) = FF(inode+lofsu)-daux*p_DM(ia)*p_Dvector(J+ioffsetl1)
            FF(inode+lofsv) = FF(inode+lofsv)-daux*p_DM(ia)*p_Dvector(J+ioffsetl2)

            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsu,k+lofsl1) = daux*p_DM(ia)
                AA (inode+lofsv,k+lofsl2) = daux*p_DM(ia)
                EXIT
              END IF
            END DO          
          END DO
        END IF

        ! The defect in the coupling matrix from the primal to the dual system
        IF (ASSOCIATED(p_DR41)) THEN
          ! Get the multipliers
          daux = rvanca%Dmultipliers(4,1)
          daux2 = rvanca%Dmultipliers(5,2)
          
          ia1 = p_KldM(idof)
          ia2 = p_KldM(idof+1)-1
          DO ia = ia1,ia2

            J = p_KcolM(ia)

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   ( bR   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .  bR   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            FF(inode+lofsl1) = FF(inode+lofsl1)-daux*p_DR41(ia)*p_Dvector(J+ioffsetu)
            FF(inode+lofsl2) = FF(inode+lofsl2)-daux2*p_DR52(ia)*p_Dvector(J+ioffsetv)
            
            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsl1,k+lofsu) = daux*p_DR41(ia)
                AA (inode+lofsl2,k+lofsv) = daux2*p_DR52(ia)
                EXIT
              END IF
            END DO          
          END DO
        END IF

        IF (ASSOCIATED(p_DR51)) THEN
          ! Get the multipliers
          daux = rvanca%Dmultipliers(5,1)
          daux2 = rvanca%Dmultipliers(4,2)
          
          ia1 = p_KldM(idof)
          ia2 = p_KldM(idof+1)-1
          DO ia = ia1,ia2

            J = p_KcolM(ia)

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .  bR   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   ( bR   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            FF(inode+lofsl1) = FF(inode+lofsl1)-daux2*p_DR42(ia)*p_Dvector(J+ioffsetv)
            FF(inode+lofsl2) = FF(inode+lofsl2)-daux*p_DR51(ia)*p_Dvector(J+ioffsetu)
            
            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsl2,k+lofsu) = daux*p_DR51(ia)
                AA (inode+lofsl1,k+lofsv) = daux2*p_DR42(ia)
                EXIT
              END IF
            END DO          
          END DO
        END IF
        
        ! THe next loop will determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
        
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   ( D1  D2   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .  D1  D2   .  ) ( xi )
          !
          ! In AA, we simultaneously set up (locally):
          !
          !   (  .   .  B1   .   .   .  ) 
          !   (  .   .  B2   .   .   .  ) 
          !   ( D1  D2   .   .   .   .  ) 
          !   (  .   .   .   .   .  B1  ) 
          !   (  .   .   .   .   .  B2  ) 
          !   (  .   .   .  D1  D2   .  ) 

          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
            
            ! Get the entries in the B-matrices.
            ! Primal equation
            AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanca%Dmultipliers(1,3)
            AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanca%Dmultipliers(2,3)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanca%Dmultipliers(3,1)
            AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanca%Dmultipliers(3,2)

            ! The same for the dual equation
            AA(inode+lofsl1,1+lofsxi) = p_DB1(ib) * rvanca%Dmultipliers(4,6)
            AA(inode+lofsl2,1+lofsxi) = p_DB2(ib) * rvanca%Dmultipliers(5,6)

            AA(1+lofsxi,inode+lofsl1) = p_DD1(ib) * rvanca%Dmultipliers(6,4)
            AA(1+lofsxi,inode+lofsl2) = p_DD2(ib) * rvanca%Dmultipliers(6,5)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                        - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)
          
            ! The same for the dual pressure
            FF(1+lofsxi) = FF(1+lofsxi) &
                         - AA(1+lofsxi,inode+lofsl1)*p_Dvector(idof+ioffsetl1) &
                         - AA(1+lofsxi,inode+lofsl2)*p_Dvector(idof+ioffsetl2)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! For C, we use our local AA, i.e. applying C^{-1} means to
      ! solve the local system AA dd = FF for dd. The local defect dd is then
      ! added back to the global solution vector.
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
        
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y
        
        DO inode=1,nnvel
          ! Update of the primal velocity vectors
          p_Dvector(idofGlobal(inode)+ioffsetu) &
            = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)

          ! Update of the dual velocity vectors
          p_Dvector(idofGlobal(inode)+ioffsetl1) &
            = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * FF(inode+lofsl1)
          p_Dvector(idofGlobal(inode)+ioffsetl2) &
            = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * FF(inode+lofsl2)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)

        p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + &
                                  domega * FF(1+lofsxi)
      
      ELSE IF (ilapackInfo .LT. 0) THEN
        
        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNSSOCQ1TQ0fullCoupConf')
        
      END IF

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D VANCA, 'full' version for fully coupled Navier-Stokes systems with
  ! primal and dual equations. Forward-Backward-variant that handles only either
  ! the primal or the dual system.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1  aM           )
  !    ( A21  A22  B2       aM      )
  !    ( D1^T D2^T I1               )
  !    ( bM            A44  A45  B1 )
  !    (      bM       A54  A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! The mass matrices must all be the same except for their scaling factors!
  ! Here, a and b are arbitrary multiplication factors for the mass matrices.
  !
  ! I1 and I2 are two diagonal matrices in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! don't exist, what results in a '0' block in these positions.
  !
  ! This variant decouples the primal system from the dual one. Depending on
  ! the parameter csystemType in the VANCA structure, only parts of the system
  ! are processed:
  !
  ! csystemType = 1: Apply VANCA to the system
  !
  !    ( A11  A12  B1  aM           )
  !    ( A21  A22  B2       aM      )
  !    ( D1^T D2^T I1               )
  !    (               I            )
  !    (                    I       )
  !    (                         I  )
  !
  ! csystemType = 2: Apply VANCA to the system
  !
  !    ( I                          )
  !    (      I                     )
  !    (           I                )
  !    ( bM            A44  A45  B1 )
  !    (      bM       A54  A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! This variant is typically used in a nonstationary environment as Gauss-Seidel
  ! preconditioner which is at first applied to the primal equation before being
  ! applied to the dual equation.
  !
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DNSSOCQ1TQ0fullCoupCfFB (rvanca, rvector, rrhs, domega, IelementList,&
      csystemPart)
  
!<description>
  ! This routine applies the specialised forward-backward full local system 
  ! VANCA algorithm for 2D Navier Stokes optimal control problems with 
  ! Q1~/Q0 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_2DNSSOCQ1TQ0fullCoupCfFB supports fully coupled velocity submatrices.
  ! The matrices A11, A22, A44 and A55 must have the same structure. 
  ! The matrices A12 and A21 must have the same structure. 
  ! The matrices A45 and A54 must have the same structure. 
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vancaPointer2DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DNavStOptC), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList

  ! Identifier for the part of the equation, where VANCA should be
  ! applied.
  ! =0: apply VANCA to the whole primal-dual-coupled system
  ! =1: apply VANCA only to the primal system, take the dual one as constant
  ! =2: apply VANCA only to the dual system, take the primal one as constant
  INTEGER, INTENT(IN) :: csystemPart
  
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA11
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA11
    REAL(DP), DIMENSION(:), POINTER             :: p_DA11,p_DA12,p_DA21,p_DA22
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA45
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA45
    REAL(DP), DIMENSION(:), POINTER             :: p_DA44,p_DA45,p_DA54,p_DA55
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA12
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA12
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolM
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldM
    REAL(DP), DIMENSION(:), POINTER             :: p_DM
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    REAL(DP), DIMENSION(:), POINTER             :: p_Da33,p_Da66
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA33,p_KdiagonalA66
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 4      ! Q1T = 4 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 1 ! QQ0 = 1 DOF's per pressure
    INTEGER, PARAMETER :: nndualvel = 4      ! Q1T = 4 DOF's per dual velocity
    INTEGER, PARAMETER :: nndualpressure = 1 ! QQ0 = 1 DOF's per dual pressure
    INTEGER, PARAMETER :: nnprimal = 2*nnvel+nnpressure ! Q1~/Q1~/Q0 = 4+4+1 = 9 DOF's per element
    INTEGER, PARAMETER :: nnld = 2*nnprimal
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! Offsets of the 'local' solution parts in the 'local' solution vector
    INTEGER, PARAMETER :: lofsu = 0
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    INTEGER, PARAMETER :: lofsl1 = 2*nnvel+1
    INTEGER, PARAMETER :: lofsl2 = 2*nnvel+1+nnvel
    INTEGER, PARAMETER :: lofsxi = 2*nnvel+1+2*nnvel
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! Offset information in arrays.
    ! Primal variables
    INTEGER(PREC_VECIDX)     :: ioffsetu,ioffsetv,ioffsetp,j
    
    ! Dual variables
    INTEGER(PREC_VECIDX)     :: ioffsetl1,ioffsetl2,ioffsetxi
    
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,k
    REAL(DP) :: daux,daux2
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    ! Structure of A11 is assumed to be the same as A22
    p_KcolA11 => rvanca%p_KcolA11
    p_KldA11 => rvanca%p_KldA11
    p_DA11 => rvanca%p_DA11
    p_DA22 => rvanca%p_DA22

    ! Structure of A12 is assumed to be the same as A21.
    ! Get A12 and A21 -- except for if the multipliers are =0, then
    ! we switch them off by nullifying the pointers.
    IF (rvanca%Dmultipliers(1,2) .NE. 0.0_DP) THEN
      p_KcolA12 => rvanca%p_KcolA12
      p_KldA12 => rvanca%p_KldA12
      p_DA12 => rvanca%p_DA12
      p_DA21 => rvanca%p_DA21
    ELSE
      NULLIFY(p_KcolA12)
      NULLIFY(p_KldA12) 
      NULLIFY(p_DA12 )
      NULLIFY(p_DA21 )
    END IF
    
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Structure of A44 is assumed to be the same as A55, A11 and A22
    p_DA44 => rvanca%p_DA44
    p_DA55 => rvanca%p_DA55

    ! Structure of A45 is assumed to be the same as A54
    p_KcolA45 => rvanca%p_KcolA45
    p_KldA45 => rvanca%p_KldA45
    p_DA45 => rvanca%p_DA45
    p_DA54 => rvanca%p_DA54
    
    ! Mass matrix - if it's given, otherwise the pointers will be set to NULL
    ! because of the initialisation of the structure!
    p_KcolM => rvanca%p_KcolM
    p_KldM => rvanca%p_KldM
    p_DM => rvanca%p_DM14
    
    ! Diagonal submatrices A33 and A66 (if they exist)
    IF (rvanca%Dmultipliers(3,3) .NE. 0.0_DP) THEN
      p_Da33 => rvanca%p_DA33
      p_KdiagonalA33 => rvanca%p_KdiagonalA33
    ELSE
      NULLIFY(p_Da33)
      NULLIFY(p_KdiagonalA33)
    END IF
    
    IF (rvanca%Dmultipliers(6,6) .NE. 0.0_DP) THEN
      p_Da66 => rvanca%p_DA66
      p_KdiagonalA66 => rvanca%p_KdiagonalA66
    ELSE
      NULLIFY(p_Da66)
      NULLIFY(p_KdiagonalA66)
    END IF

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    ! Get the offsets of lambda1, lambda2 and xi, so the offsets
    ! of the dual solution vectors.
    ioffsetl1 = ioffsetp+rvector%RvectorBlock(3)%NEQ
    ioffsetl2 = ioffsetl1+rvector%RvectorBlock(4)%NEQ
    ioffsetxi = ioffsetl2+rvector%RvectorBlock(5)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new     
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new     
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new       new    
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new    
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    SELECT CASE (csystemPart)
    ! -------------------------------------------------------------------------
    CASE (0)
      ! VANCA applied to full 18x18 system

      DO ielidx=1,SIZE(IelementList)
      
        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)
      
        ! Clear the 'local system matrix'.
        AA(:,:) = 0.0_DP
        
        ! We now have the element
        !                                               
        ! +---------+                       +----3----+
        ! |         |                       |         |
        ! |   IEL   |   with DOF's          4    P    2      
        ! |         |                       |    Q0   |
        ! +---------+                       +----1----+
        !                                               
        !
        ! Fetch the pressure P on the current element into FF.
        ! The numbers of the DOF's coincide with the definition
        ! in dofmapping.f90!
      
        ! Get the primal pressure
        FF(1+lofsp) = p_Drhs(iel+ioffsetp)
        
        ! Get the dual pressure
        FF(1+lofsxi) = p_Drhs(iel+ioffsetxi)
        
        ! Get the velocity DOF's on the current element.
        ! We assume: DOF 1..4 = edge.
        ! That's the same implementation as in dofmapping.f90!
        IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

        ! Loop over all U-nodes of that element.
        DO inode=1,nnvel
        
          ! Get the DOF we have to tackle:
          idof = IdofGlobal(inode)
          
          ! Set FF initially to the value of the right hand
          ! side vector that belongs to our current DOF corresponding
          ! to inode.
          
          ! Primal equation
          FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
          FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

          ! dual equation
          FF(inode+lofsl1) = p_Drhs(idof+ioffsetl1)
          FF(inode+lofsl2) = p_Drhs(idof+ioffsetl2)
          
          ! What do we have at this point?                           
          ! FF     : "local" RHS vector belonging to the DOF's on the
          !          current element                                 
          ! AA     : Diagonal entries of A belonging to these DOF's  
          !                                                          
          ! And at the moment:                                       
          ! idof      : number of current DOF on element IEL            
          ! inode     : "local" number of DOF on element IEL, i.e.      
          !              number of the edge         
          !                     
          ! Now comes the crucial point with the "update": How to         
          ! subsequently update the vertex values, such that the whole    
          ! thing still converges to the solution, even if a node         
          ! is updated more than once? Here, we use a typical             
          ! matrix-decomposition approach:                                
          !                                                               
          ! Again consider the problem:                                   
          !                                                               
          !    [ A   B  M     ] (u)  = (f  )                                        
          !    [ B^t 0        ] (p)    (g  )                                        
          !    [ M      A   B ] (l)    (fl )                                        
          !    [        B^t 0 ] (xi)   (fxi)                                        
          !                                                               
          ! We assume, that all components in the vector (u,p) are        
          ! given - except for the velocity and pressure unknowns 
          ! on the current element; these 21 unknowns  
          ! are located anywhere in the (u,p) vector. The idea is to      
          ! shift "everything known" to the right hand side to obtain     
          ! a system for only these unknowns!   
          !
          ! Extracting all the lines of the system that correspond to     
          ! DOF's on our single element IEL results in rectangular
          ! systems of the form                                            
          !                                                               
          !    [ === A^ === B~  === M^ ====   ] (| ) = (f1 )                                
          !    [ B~^t       I1~               ] (u )   (f2 )                                
          !    [ === M^ ===     === A^ === B~ ] (| )   (g  )                                   
          !    [                B~^t       I2~] (p )   (fl1)
          !                                     (| )   (fl2)
          !                                     (l )   (flg)
          !                                     (| )
          !                                     (xi)
          !                                     
          !                                                               
          ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
          ! 2*1 pressure elements on the adjacent cells, so we have       
          ! 2 columns in the B-matrix.                                    
          !                                                               
          !        IEL                              IEL                   
          !     |--------|             |--------|--------|                
          !     |        |             |        |        |                
          !     |   P    |      or     |   Q    X   P    |                
          !     |   X    |             |        |        |                
          !   --|--------|--           |--------|--------|                
          !
          !
          ! Now, throw all summands to the RHS vector to build a local
          ! 'defect' on our single element IEL.
          !                                                               
          !  (d1 ) = (f1 ) -  [ === A^ === B~  === M^ ====   ] (| )                                 
          !  (d2 )   (f2 )    [ B~^t       I1~               ] (u )                                 
          !  (dg )   (g  )    [ === M^ ===     === A^ === B~ ] (| )                                    
          !  (dl1)   (fl1)    [                B~^t       I2~] (p ) 
          !  (dl2)   (fl2)                                     (| ) 
          !  (dlg)   (flg)                                     (l ) 
          !                                                    (| )
          !                                                    (xi)
          !
          ! Extract those entries in the A-, B- and M-matrices to our local
          ! matrix AA, which belong to DOF's in our current solution vector.
          !
          ! At first build: fi = fi-Aui
          
          ia1 = p_KldA11(idof)
          ia2 = p_KldA11(idof+1)-1
          DO ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - ( A11  .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .  A22  .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .  A44  .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .  A55  .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA11(ia)
            
            ! Primal equation:
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanca%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanca%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

            ! dual equation
            FF(inode+lofsl1) = FF(inode+lofsl1) &
                            - rvanca%Dmultipliers(4,4)*p_DA44(ia)*p_Dvector(J+ioffsetl1)
            FF(inode+lofsl2) = FF(inode+lofsl2) &
                            - rvanca%Dmultipliers(5,5)*p_DA55(ia)*p_Dvector(J+ioffsetl2)
            
            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanca%Dmultipliers(1,1)
                AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanca%Dmultipliers(2,2)
                
                AA (inode+lofsl1,k+lofsl1) = p_DA44(ia)*rvanca%Dmultipliers(4,4)
                AA (inode+lofsl2,k+lofsl2) = p_DA55(ia)*rvanca%Dmultipliers(5,5)
                EXIT
              END IF
            END DO          
          END DO

          ! Handle the 'off-diagonal' matrices A12 and A21
          
          IF (ASSOCIATED(p_KldA12)) THEN
            ia1 = p_KldA12(idof)
            ia2 = p_KldA12(idof+1)-1
            DO ia = ia1,ia2
              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .  A12  .   .   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   ( A21  .   .   .   .   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              J = p_KcolA12(ia)
              FF(inode+lofsu) = FF(inode+lofsu) &
                              - rvanca%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
              FF(inode+lofsv) = FF(inode+lofsv) &
                              - rvanca%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)
              
              ! Whereever we find a DOF that couples to another DOF on the 
              ! same element, we put that to both A-blocks of our local matrix.
              DO k=1,nnvel
                IF (j .EQ. IdofGlobal(k)) THEN
                  AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanca%Dmultipliers(1,2)
                  AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanca%Dmultipliers(2,1)
                  EXIT
                END IF
              END DO          
            END DO
          END IF
          
          ! Handle the 'off-diagonal' matrices A45 and A54 if they exist
          IF (ASSOCIATED(p_KldA45)) THEN
            ia1 = p_KldA45(idof)
            ia2 = p_KldA45(idof+1)-1
            DO ia = ia1,ia2
              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !   ( dl1 )   ( dl1 )   (  .   .   .   .  A45  .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .   .   .  A54  .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              J = p_KcolA45(ia)
              FF(inode+lofsl1) = FF(inode+lofsl1) &
                              - rvanca%Dmultipliers(4,5)*p_DA45(ia)*p_Dvector(J+ioffsetl2)
              FF(inode+lofsl2) = FF(inode+lofsl2) &
                              - rvanca%Dmultipliers(5,4)*p_DA54(ia)*p_Dvector(J+ioffsetl1)
              
              ! Whereever we find a DOF that couples to another DOF on the 
              ! same element, we put that to both A-blocks of our local matrix.
              DO k=1,nnvel
                IF (j .EQ. IdofGlobal(k)) THEN
                  AA (inode+lofsl1,k+lofsl2) = p_DA45(ia)*rvanca%Dmultipliers(4,5)
                  AA (inode+lofsl2,k+lofsl1) = p_DA54(ia)*rvanca%Dmultipliers(5,4)
                  EXIT
                END IF
              END DO          
            END DO
          END IF
                  
          ! Process A33 if it exists
                  
          IF (ASSOCIATED(p_KdiagonalA33)) THEN

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   I1  .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.
            
            daux = rvanca%Dmultipliers(3,3)
            FF(1+lofsp) = FF(1+lofsp) &
                        - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
            AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

          END IF
                  
          ! Process A66 if it exists
                  
          IF (ASSOCIATED(p_KdiagonalA66)) THEN

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   I2 ) ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.
            
            daux = rvanca%Dmultipliers(6,6)
            FF(1+lofsxi) = FF(1+lofsxi) &
                        - daux*p_DA66(p_KdiagonalA66(IEL))*p_Dvector(IEL+ioffsetxi)
            AA(1+lofsxi,1+lofsxi) = daux*p_DA66(p_KdiagonalA66(IEL))

          END IF
                  
          ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
          
          ib1=p_KldB(idof)
          ib2=p_KldB(idof+1)-1
          DO ib = ib1,ib2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .  B1   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .  B2   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .  B1  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .  B2  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolB(ib)
            
            ! primal equation
            daux = p_Dvector(j+ioffsetp) 
            FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanca%Dmultipliers(1,3)
            FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanca%Dmultipliers(2,3)

            ! dual equation
            daux2 = p_Dvector(j+ioffsetxi)
            FF(inode+lofsl1) = FF(inode+lofsl1)-p_DB1(ib)*daux2 * rvanca%Dmultipliers(4,6)
            FF(inode+lofsl2) = FF(inode+lofsl2)-p_DB2(ib)*daux2 * rvanca%Dmultipliers(5,6)
            
            ! Don't incorporate the B-matrices into AA yet; this will come later!
          END DO
          
          ! The mass matrix defect.
          IF (ASSOCIATED(p_KldM)) THEN
            ! We assume: multiplier of A(1,4) = multiplier of A(2,5)
            daux = rvanca%Dmultipliers(1,4)
            
            ! We assume: multiplier of A(4,1) = multiplier of A(5,2)
            daux2 = rvanca%Dmultipliers(4,1)
            
            ia1 = p_KldM(idof)
            ia2 = p_KldM(idof+1)-1
            DO ia = ia1,ia2

              J = p_KcolM(ia)

              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .   .   .  aM   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   (  .   .   .   .  aM   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !   ( dl1 )   ( dl1 )   ( bM   .   .   .   .   .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .  bM   .   .   .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              FF(inode+lofsu) = FF(inode+lofsu)-daux*p_DM(ia)*p_Dvector(J+ioffsetl1)
              FF(inode+lofsv) = FF(inode+lofsv)-daux*p_DM(ia)*p_Dvector(J+ioffsetl2)

              FF(inode+lofsl1) = FF(inode+lofsl1)-daux2*p_DM(ia)*p_Dvector(J+ioffsetu)
              FF(inode+lofsl2) = FF(inode+lofsl2)-daux2*p_DM(ia)*p_Dvector(J+ioffsetv)
              
              ! Whereever we find a DOF that couples to another DOF on the 
              ! same element, we put that to both A-blocks of our local matrix.
              DO k=1,nnvel
                IF (j .EQ. IdofGlobal(k)) THEN
                  AA (inode+lofsu,k+lofsl1) = daux*p_DM(ia)
                  AA (inode+lofsv,k+lofsl2) = daux*p_DM(ia)
                  
                  AA (inode+lofsl1,k+lofsu) = daux2*p_DM(ia)
                  AA (inode+lofsl2,k+lofsv) = daux2*p_DM(ia)
                  ! AA (k+lofsl1,inode+lofsu) = daux2*p_DM(ia)
                  ! AA (k+lofsl2,inode+lofsv) = daux2*p_DM(ia)
                  EXIT
                END IF
              END DO          
            END DO
          END IF
          
          ! Ok, up to now, all loops are clean and vectoriseable. Now the only
          ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
          ! We have to find in the B-matrices the column that corresponds
          ! to our element and pressure DOF IEL - which makes it necessary
          ! to compare the column numbers in KcolB with IEL.
          ! Remember: The column numbers in B correspond to the pressure-DOF's
          ! and so to element numbers. 
          !
          ! Btw: Each row of B has at most two entries:
          !
          !      IEL                              IEL
          !   |--------|             |--------|--------|
          !   |        |             |        |        |
          !   |   P1   |      or     |   P2   X   P1   |
          !   |        |             |        |        |
          ! --|---X----|--           |--------|--------|
          !
          ! Either two (if the velocity DOF is an edge with two neighbouring
          ! elements) or one (if the velocity DOF is at an edge on the boundary
          ! and there is no neighbour).
          DO ib = ib1,ib2
          
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   ( D1  D2   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .  D1  D2   .  ) ( xi )
            !
            ! In AA, we simultaneously set up (locally):
            !
            !   (  .   .  B1   .   .   .  ) 
            !   (  .   .  B2   .   .   .  ) 
            !   ( D1  D2   .   .   .   .  ) 
            !   (  .   .   .   .   .  B1  ) 
            !   (  .   .   .   .   .  B2  ) 
            !   (  .   .   .  D1  D2   .  ) 

            IF (p_KcolB(ib) .EQ. IEL) THEN
            
              J = p_KcolB(ib)
              
              ! Get the entries in the B-matrices.
              ! Primal equation
              AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanca%Dmultipliers(1,3)
              AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanca%Dmultipliers(2,3)

              ! The same way, get DD1 and DD2.
              ! Note that DDi has exacty the same matrix structrure as BBi and is noted
              ! as 'transposed matrix' only because of the transposed-flag.
              ! So we can use "ib" as index here to access the entry of DDi:
              AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanca%Dmultipliers(3,1)
              AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanca%Dmultipliers(3,2)

              ! The same for the dual equation
              AA(inode+lofsl1,1+lofsxi) = p_DB1(ib) * rvanca%Dmultipliers(4,6)
              AA(inode+lofsl2,1+lofsxi) = p_DB2(ib) * rvanca%Dmultipliers(5,6)

              AA(1+lofsxi,inode+lofsl1) = p_DD1(ib) * rvanca%Dmultipliers(6,4)
              AA(1+lofsxi,inode+lofsl2) = p_DD2(ib) * rvanca%Dmultipliers(6,5)

              ! Build the pressure entry in the local defect vector:
              !   f_i = (f_i-Aui) - D_i pi
              ! or more precisely (as D is roughly B^T):
              !   f_i = (f_i-Aui) - (B^T)_i pi
              FF(1+lofsp) = FF(1+lofsp) &
                          - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                          - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)
            
              ! The same for the dual pressure
              FF(1+lofsxi) = FF(1+lofsxi) &
                          - AA(1+lofsxi,inode+lofsl1)*p_Dvector(idof+ioffsetl1) &
                          - AA(1+lofsxi,inode+lofsl2)*p_Dvector(idof+ioffsetl2)
            
              ! Quit the loop - the other possible entry belongs to another 
              ! element, not to the current one
              EXIT
            END IF
          END DO ! ib
          
        END DO ! inode
      
        ! Now we make a defect-correction approach for this system:
        !
        !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
        !                                     -----------
        !                                        =d~
        !
        ! Here the 'projection' operator simply converts the small
        ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
        ! of the same size as x - what is easy using the number of
        ! the DOF's on the element.
        !
        ! For C, we use our local AA, i.e. applying C^{-1} means to
        ! solve the local system AA dd = FF for dd. The local defect dd is then
        ! added back to the global solution vector.
        
        CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
        
        IF (ilapackInfo .EQ. 0) THEN
          
          ! Ok, we got the update vector in FF. Incorporate this now into our
          ! solution vector with the update formula
          !
          !  x_{n+1} = x_n + domega * y
          
          DO inode=1,nnvel
            ! Update of the primal velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetu) &
              = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
            p_Dvector(idofGlobal(inode)+ioffsetv) &
              = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)

            ! Update of the dual velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetl1) &
              = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * FF(inode+lofsl1)
            p_Dvector(idofGlobal(inode)+ioffsetl2) &
              = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * FF(inode+lofsl2)
          END DO
          
          p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                    domega * FF(1+lofsp)

          p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + &
                                    domega * FF(1+lofsxi)
        
        ELSE IF (ilapackInfo .LT. 0) THEN
          
          CALL output_line (&
              'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
              OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNSSOCQ1TQ0fullCoupCfFB')
          
        END IF

        ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
        ! coarse grid with all boundaries = Dirichlet.
        ! In this case, nothing must be changed in the vector!
      
      END DO ! iel

    ! -------------------------------------------------------------------------      
    CASE (1)
      ! VANCA applied to 9x9 primal system

      DO ielidx=1,SIZE(IelementList)
      
        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)
      
        ! Clear the 'primal local system matrix'.
        DO ib=1,nnprimal
          DO ia=1,nnprimal
            AA(ia,ib) = 0.0_DP
          END DO
        END DO
        
        ! We now have the element
        !                                               
        ! +---------+                       +----3----+
        ! |         |                       |         |
        ! |   IEL   |   with DOF's          4    P    2      
        ! |         |                       |    Q0   |
        ! +---------+                       +----1----+
        !                                               
        !
        ! Fetch the pressure P on the current element into FF.
        ! The numbers of the DOF's coincide with the definition
        ! in dofmapping.f90!
      
        ! Get the primal pressure
        FF(1+lofsp) = p_Drhs(iel+ioffsetp)
        
        ! Get the velocity DOF's on the current element.
        ! We assume: DOF 1..4 = edge.
        ! That's the same implementation as in dofmapping.f90!
        IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

        ! Loop over all U-nodes of that element.
        DO inode=1,nnvel
        
          ! Get the DOF we have to tackle:
          idof = IdofGlobal(inode)
          
          ! Set FF initially to the value of the right hand
          ! side vector that belongs to our current DOF corresponding
          ! to inode.
          
          ! Primal equation
          FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
          FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

          ! What do we have at this point?                           
          ! FF     : "local" RHS vector belonging to the DOF's on the
          !          current element                                 
          ! AA     : Diagonal entries of A belonging to these DOF's  
          !                                                          
          ! And at the moment:                                       
          ! idof      : number of current DOF on element IEL            
          ! inode     : "local" number of DOF on element IEL, i.e.      
          !              number of the edge         
          !                     

          ! Extract those entries in the A-, B- and M-matrices to our local
          ! matrix AA, which belong to DOF's in our current solution vector.
          !
          ! At first build: fi = fi-Aui
          
          ia1 = p_KldA11(idof)
          ia2 = p_KldA11(idof+1)-1
          DO ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - ( A11  .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .  A22  .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )

            J = p_KcolA11(ia)
            
            ! Primal equation:
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanca%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanca%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanca%Dmultipliers(1,1)
                AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanca%Dmultipliers(2,2)
                EXIT
              END IF
            END DO          
          END DO

          ! Handle the 'off-diagonal' matrices A12 and A21
          
          IF (ASSOCIATED(p_KldA12)) THEN
            ia1 = p_KldA12(idof)
            ia2 = p_KldA12(idof+1)-1
            DO ia = ia1,ia2
              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .  A12  .   .   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   ( A21  .   .   .   .   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !                                                   ( l1 )
              !                                                   ( l2 )
              !                                                   ( xi )

              J = p_KcolA12(ia)
              FF(inode+lofsu) = FF(inode+lofsu) &
                              - rvanca%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
              FF(inode+lofsv) = FF(inode+lofsv) &
                              - rvanca%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)
              
              ! Whereever we find a DOF that couples to another DOF on the 
              ! same element, we put that to both A-blocks of our local matrix.
              DO k=1,nnvel
                IF (j .EQ. IdofGlobal(k)) THEN
                  AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanca%Dmultipliers(1,2)
                  AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanca%Dmultipliers(2,1)
                  EXIT
                END IF
              END DO          
            END DO
          END IF
          
          ! Process A33 if it exists
                  
          IF (ASSOCIATED(p_KdiagonalA33)) THEN

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   I1  .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.
            
            daux = rvanca%Dmultipliers(3,3)
            FF(1+lofsp) = FF(1+lofsp) &
                        - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
            AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

          END IF
                  
          ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
          
          ib1=p_KldB(idof)
          ib2=p_KldB(idof+1)-1
          DO ib = ib1,ib2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .  B1   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .  B2   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )

            J = p_KcolB(ib)
            
            ! primal equation
            daux = p_Dvector(j+ioffsetp) 
            FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanca%Dmultipliers(1,3)
            FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanca%Dmultipliers(2,3)

            ! Don't incorporate the B-matrices into AA yet; this will come later!
          END DO
          
          ! The mass matrix defect.
          IF (ASSOCIATED(p_KldM)) THEN
            ! We assume: multiplier of A(1,4) = multiplier of A(2,5)
            daux = rvanca%Dmultipliers(1,4)
            
            ia1 = p_KldM(idof)
            ia2 = p_KldM(idof+1)-1
            DO ia = ia1,ia2

              J = p_KcolM(ia)

              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .   .   .  aM   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   (  .   .   .   .  aM   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !                                                   ( l1 )
              !                                                   ( l2 )
              !                                                   ( xi )

              FF(inode+lofsu) = FF(inode+lofsu)-daux*p_DM(ia)*p_Dvector(J+ioffsetl1)
              FF(inode+lofsv) = FF(inode+lofsv)-daux*p_DM(ia)*p_Dvector(J+ioffsetl2)

              ! Whereever we find a DOF that couples to another DOF on the 
              ! same element, we put that to both A-blocks of our local matrix.
              DO k=1,nnvel
                IF (j .EQ. IdofGlobal(k)) THEN
                  AA (inode+lofsu,k+lofsl1) = daux*p_DM(ia)
                  AA (inode+lofsv,k+lofsl2) = daux*p_DM(ia)
                  EXIT
                END IF
              END DO          
            END DO
          END IF
          
          ! Ok, up to now, all loops are clean and vectoriseable. Now the only
          ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
          ! We have to find in the B-matrices the column that corresponds
          ! to our element and pressure DOF IEL - which makes it necessary
          ! to compare the column numbers in KcolB with IEL.
          ! Remember: The column numbers in B correspond to the pressure-DOF's
          ! and so to element numbers. 
          !
          ! Btw: Each row of B has at most two entries:
          !
          !      IEL                              IEL
          !   |--------|             |--------|--------|
          !   |        |             |        |        |
          !   |   P1   |      or     |   P2   X   P1   |
          !   |        |             |        |        |
          ! --|---X----|--           |--------|--------|
          !
          ! Either two (if the velocity DOF is an edge with two neighbouring
          ! elements) or one (if the velocity DOF is at an edge on the boundary
          ! and there is no neighbour).
          DO ib = ib1,ib2
          
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   ( D1  D2   .   .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )
            !
            ! In AA, we simultaneously set up (locally):
            !
            !   (  .   .  B1   .   .   .  ) 
            !   (  .   .  B2   .   .   .  ) 
            !   ( D1  D2   .   .   .   .  ) 
            !   (  .   .   .   .   .   .  ) 
            !   (  .   .   .   .   .   .  ) 
            !   (  .   .   .   .   .   .  ) 

            IF (p_KcolB(ib) .EQ. IEL) THEN
            
              J = p_KcolB(ib)
              
              ! Get the entries in the B-matrices.
              ! Primal equation
              AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanca%Dmultipliers(1,3)
              AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanca%Dmultipliers(2,3)

              ! The same way, get DD1 and DD2.
              ! Note that DDi has exacty the same matrix structrure as BBi and is noted
              ! as 'transposed matrix' only because of the transposed-flag.
              ! So we can use "ib" as index here to access the entry of DDi:
              AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanca%Dmultipliers(3,1)
              AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanca%Dmultipliers(3,2)

              ! Build the pressure entry in the local defect vector:
              !   f_i = (f_i-Aui) - D_i pi
              ! or more precisely (as D is roughly B^T):
              !   f_i = (f_i-Aui) - (B^T)_i pi
              FF(1+lofsp) = FF(1+lofsp) &
                          - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                          - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)
            
              ! Quit the loop - the other possible entry belongs to another 
              ! element, not to the current one
              EXIT
            END IF
          END DO ! ib
          
        END DO ! inode
      
        ! Now we make a defect-correction approach for this system:
        !
        !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
        !                                     -----------
        !                                        =d~
        !
        ! Here the 'projection' operator simply converts the small
        ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
        ! of the same size as x - what is easy using the number of
        ! the DOF's on the element.
        !
        ! For C, we use our local AA, i.e. applying C^{-1} means to
        ! solve the local system AA dd = FF for dd. The local defect dd is then
        ! added back to the global solution vector.
        !
        ! We only set up the primal equation, so we only have to solve a
        ! 9x9 subsystem of the complete 18x18 system. Lapack can handle this
        ! using the 'leading dimension'...
        
        CALL DGESV (nnprimal, 1, AA, nnld, Ipiv, FF(1:nnprimal), nnprimal, ilapackInfo)
        
        IF (ilapackInfo .EQ. 0) THEN
          
          ! Ok, we got the update vector in FF. Incorporate this now into our
          ! solution vector with the update formula
          !
          !  x_{n+1} = x_n + domega * y
          
          DO inode=1,nnvel
            ! Update of the primal velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetu) &
              = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
            p_Dvector(idofGlobal(inode)+ioffsetv) &
              = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)

          END DO
          
          p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                    domega * FF(1+lofsp)

        ELSE IF (ilapackInfo .LT. 0) THEN
          
          CALL output_line (&
              'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
              OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNSSOCQ1TQ0fullCoupCfFB')
          
        END IF

        ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
        ! coarse grid with all boundaries = Dirichlet.
        ! In this case, nothing must be changed in the vector!
      
      END DO ! iel
  
    ! -------------------------------------------------------------------------
    CASE (2)
      ! VANCA applied to 9x9 dual system

      DO ielidx=1,SIZE(IelementList)
      
        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)
      
        ! Clear the 'dual local system matrix'.
        DO ib=nnprimal+1,2*nnprimal
          DO ia=nnprimal+1,2*nnprimal
            AA(ia,ib) = 0.0_DP
          END DO
        END DO
        
        ! We now have the element
        !                                               
        ! +---------+                       +----3----+
        ! |         |                       |         |
        ! |   IEL   |   with DOF's          4    P    2      
        ! |         |                       |    Q0   |
        ! +---------+                       +----1----+
        !                                               
        !
        ! Fetch the pressure P on the current element into FF.
        ! The numbers of the DOF's coincide with the definition
        ! in dofmapping.f90!
      
        ! Get the dual pressure
        FF(1+lofsxi) = p_Drhs(iel+ioffsetxi)
        
        ! Get the velocity DOF's on the current element.
        ! We assume: DOF 1..4 = edge.
        ! That's the same implementation as in dofmapping.f90!
        IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

        ! Loop over all U-nodes of that element.
        DO inode=1,nnvel
        
          ! Get the DOF we have to tackle:
          idof = IdofGlobal(inode)
          
          ! Set FF initially to the value of the right hand
          ! side vector that belongs to our current DOF corresponding
          ! to inode.
          
          ! dual equation
          FF(inode+lofsl1) = p_Drhs(idof+ioffsetl1)
          FF(inode+lofsl2) = p_Drhs(idof+ioffsetl2)
          
          ! What do we have at this point?                           
          ! FF     : "local" RHS vector belonging to the DOF's on the
          !          current element                                 
          ! AA     : Diagonal entries of A belonging to these DOF's  
          !                                                          
          ! And at the moment:                                       
          ! idof      : number of current DOF on element IEL            
          ! inode     : "local" number of DOF on element IEL, i.e.      
          !              number of the edge         
          !                     
          ! At first build: fi = fi-Aui
          
          ia1 = p_KldA11(idof)
          ia2 = p_KldA11(idof+1)-1
          DO ia = ia1,ia2
            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .  A44  .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .  A55  .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA11(ia)
            
            ! dual equation
            FF(inode+lofsl1) = FF(inode+lofsl1) &
                            - rvanca%Dmultipliers(4,4)*p_DA44(ia)*p_Dvector(J+ioffsetl1)
            FF(inode+lofsl2) = FF(inode+lofsl2) &
                            - rvanca%Dmultipliers(5,5)*p_DA55(ia)*p_Dvector(J+ioffsetl2)
            
            ! Whereever we find a DOF that couples to another DOF on the 
            ! same element, we put that to both A-blocks of our local matrix.
            DO k=1,nnvel
              IF (j .EQ. IdofGlobal(k)) THEN
                AA (inode+lofsl1,k+lofsl1) = p_DA44(ia)*rvanca%Dmultipliers(4,4)
                AA (inode+lofsl2,k+lofsl2) = p_DA55(ia)*rvanca%Dmultipliers(5,5)
                EXIT
              END IF
            END DO          
          END DO


          ! Handle the 'off-diagonal' matrices A45 and A54 if they exist
          IF (ASSOCIATED(p_KldA45)) THEN
            ia1 = p_KldA45(idof)
            ia2 = p_KldA45(idof+1)-1
            DO ia = ia1,ia2
              ! Calculate:
              !
              !                                                   ( u  )
              !                                                   ( v  )
              !                                                   ( p  )
              !   ( dl1 ) = ( dl1 )   (  .   .   .   .  A45  .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .   .   .  A54  .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              J = p_KcolA45(ia)
              FF(inode+lofsl1) = FF(inode+lofsl1) &
                              - rvanca%Dmultipliers(4,5)*p_DA45(ia)*p_Dvector(J+ioffsetl2)
              FF(inode+lofsl2) = FF(inode+lofsl2) &
                              - rvanca%Dmultipliers(5,4)*p_DA54(ia)*p_Dvector(J+ioffsetl1)
              
              ! Whereever we find a DOF that couples to another DOF on the 
              ! same element, we put that to both A-blocks of our local matrix.
              DO k=1,nnvel
                IF (j .EQ. IdofGlobal(k)) THEN
                  AA (inode+lofsl1,k+lofsl2) = p_DA45(ia)*rvanca%Dmultipliers(4,5)
                  AA (inode+lofsl2,k+lofsl1) = p_DA54(ia)*rvanca%Dmultipliers(5,4)
                  EXIT
                END IF
              END DO          
            END DO
          END IF
                  
          ! Process A66 if it exists
                  
          IF (ASSOCIATED(p_KdiagonalA66)) THEN

            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   I2 ) ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.
            
            daux = rvanca%Dmultipliers(6,6)
            FF(1+lofsxi) = FF(1+lofsxi) &
                        - daux*p_DA66(p_KdiagonalA66(IEL))*p_Dvector(IEL+ioffsetxi)
            AA(1+lofsxi,1+lofsxi) = daux*p_DA66(p_KdiagonalA66(IEL))

          END IF
                  
          ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
          
          ib1=p_KldB(idof)
          ib2=p_KldB(idof+1)-1
          DO ib = ib1,ib2
            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .   .   .  B1  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .  B2  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolB(ib)
            
            ! dual equation
            daux2 = p_Dvector(j+ioffsetxi)
            FF(inode+lofsl1) = FF(inode+lofsl1)-p_DB1(ib)*daux2 * rvanca%Dmultipliers(4,6)
            FF(inode+lofsl2) = FF(inode+lofsl2)-p_DB2(ib)*daux2 * rvanca%Dmultipliers(5,6)
            
            ! Don't incorporate the B-matrices into AA yet; this will come later!
          END DO
          
          ! The mass matrix defect.
          IF (ASSOCIATED(p_KldM)) THEN
            ! We assume: multiplier of A(4,1) = multiplier of A(5,2)
            daux2 = rvanca%Dmultipliers(4,1)
            
            ia1 = p_KldM(idof)
            ia2 = p_KldM(idof+1)-1
            DO ia = ia1,ia2

              J = p_KcolM(ia)

              ! Calculate:
              !
              !                                                   ( u  )
              !                                                   ( v  )
              !                                                   ( p  )
              !   ( dl1 ) = ( dl1 )   ( bM   .   .   .   .   .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .  bM   .   .   .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              FF(inode+lofsl1) = FF(inode+lofsl1)-daux2*p_DM(ia)*p_Dvector(J+ioffsetu)
              FF(inode+lofsl2) = FF(inode+lofsl2)-daux2*p_DM(ia)*p_Dvector(J+ioffsetv)
              
              ! Whereever we find a DOF that couples to another DOF on the 
              ! same element, we put that to both A-blocks of our local matrix.
              DO k=1,nnvel
                IF (j .EQ. IdofGlobal(k)) THEN
                  AA (inode+lofsl1,k+lofsu) = daux2*p_DM(ia)
                  AA (inode+lofsl2,k+lofsv) = daux2*p_DM(ia)
                  EXIT
                END IF
              END DO          
            END DO
          END IF
          
          ! Ok, up to now, all loops are clean and vectoriseable. Now the only
          ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
          ! We have to find in the B-matrices the column that corresponds
          ! to our element and pressure DOF IEL - which makes it necessary
          ! to compare the column numbers in KcolB with IEL.
          ! Remember: The column numbers in B correspond to the pressure-DOF's
          ! and so to element numbers. 
          !
          ! Btw: Each row of B has at most two entries:
          !
          !      IEL                              IEL
          !   |--------|             |--------|--------|
          !   |        |             |        |        |
          !   |   P1   |      or     |   P2   X   P1   |
          !   |        |             |        |        |
          ! --|---X----|--           |--------|--------|
          !
          ! Either two (if the velocity DOF is an edge with two neighbouring
          ! elements) or one (if the velocity DOF is at an edge on the boundary
          ! and there is no neighbour).
          DO ib = ib1,ib2
          
            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .  D1  D2   .  ) ( xi )
            !
            ! In AA, we simultaneously set up (locally):
            !
            !   (  .   .   .   .   .   .  ) 
            !   (  .   .   .   .   .   .  ) 
            !   (  .   .   .   .   .   .  ) 
            !   (  .   .   .   .   .  B1  ) 
            !   (  .   .   .   .   .  B2  ) 
            !   (  .   .   .  D1  D2   .  ) 

            IF (p_KcolB(ib) .EQ. IEL) THEN
            
              J = p_KcolB(ib)
              
              ! Get the entries in the B-matrices.
              ! Primal equation
              AA(inode+lofsl1,1+lofsxi) = p_DB1(ib) * rvanca%Dmultipliers(4,6)
              AA(inode+lofsl2,1+lofsxi) = p_DB2(ib) * rvanca%Dmultipliers(5,6)

              ! The same way, get DD1 and DD2.
              ! Note that DDi has exacty the same matrix structrure as BBi and is noted
              ! as 'transposed matrix' only because of the transposed-flag.
              ! So we can use "ib" as index here to access the entry of DDi:
              AA(1+lofsxi,inode+lofsl1) = p_DD1(ib) * rvanca%Dmultipliers(6,4)
              AA(1+lofsxi,inode+lofsl2) = p_DD2(ib) * rvanca%Dmultipliers(6,5)

              ! Build the pressure entry in the local defect vector:
              !   f_i = (f_i-Aui) - D_i pi
              ! or more precisely (as D is roughly B^T):
              !   f_i = (f_i-Aui) - (B^T)_i pi
              FF(1+lofsxi) = FF(1+lofsxi) &
                          - AA(1+lofsxi,inode+lofsl1)*p_Dvector(idof+ioffsetl1) &
                          - AA(1+lofsxi,inode+lofsl2)*p_Dvector(idof+ioffsetl2)
            
              ! Quit the loop - the other possible entry belongs to another 
              ! element, not to the current one
              EXIT
            END IF
          END DO ! ib
          
        END DO ! inode
      
        ! Now we make a defect-correction approach for this system:
        !
        !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
        !                                     -----------
        !                                        =d~
        !
        ! Here the 'projection' operator simply converts the small
        ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
        ! of the same size as x - what is easy using the number of
        ! the DOF's on the element.
        !
        ! For C, we use our local AA, i.e. applying C^{-1} means to
        ! solve the local system AA dd = FF for dd. The local defect dd is then
        ! added back to the global solution vector.
        !
        ! We only set up the primal equation, so we only have to solve a
        ! 9x9 subsystem of the complete 18x18 system. Lapack can handle this
        ! using the 'leading dimension'...
        !
        ! Note that we explicitely using the fact that AA(10+i*18,j) = AA(10,i+j),
        ! so LAPACK will really work on the submatrix A(10:18,10:18) !
        
        CALL DGESV (nnprimal, 1, AA(nnprimal+1,nnprimal+1), nnld, &
            Ipiv, FF(nnprimal+1), nnprimal, ilapackInfo)
        
        IF (ilapackInfo .EQ. 0) THEN
          
          ! Ok, we got the update vector in FF. Incorporate this now into our
          ! solution vector with the update formula
          !
          !  x_{n+1} = x_n + domega * y
          
          DO inode=1,nnvel
            ! Update of the dual velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetl1) &
              = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * FF(inode+lofsl1)
            p_Dvector(idofGlobal(inode)+ioffsetl2) &
              = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * FF(inode+lofsl2)
          END DO
          
          p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + &
                                    domega * FF(1+lofsxi)
        
        ELSE IF (ilapackInfo .LT. 0) THEN
          
          CALL output_line (&
              'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
              OU_CLASS_ERROR,OU_MODE_STD,'vanca_2DNSSOCQ1TQ0fullCoupCfFB')
          
        END IF

        ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
        ! coarse grid with all boundaries = Dirichlet.
        ! In this case, nothing must be changed in the vector!
      
      END DO ! iel

    END SELECT

  END SUBROUTINE

! *****************************************************************************
! Problem class: VANCA variants for 3D Navier-Stokes problems
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_init3DNavierStokes (rmatrix,rvanca)
  
!<description>
  ! Initialises the VANCA variant for 3D Navier-Stokes problems 
  ! for conformal discretisations. Checks if the "3D-Navier-Stokes" VANCA variant 
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! The substructure rvanca%rvanca3DNavSt is intitialised according
  ! to the information provided in rmatrix.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: rmatrix
  !</input>

!<inputoutput>
  ! t_vancaPointer3DSPNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vanca), INTENT(INOUT) :: rvanca
!</inputoutput>

!</subroutine>

    INTEGER :: i,j
    LOGICAL :: bextended
    TYPE(t_blockDiscretisation), POINTER :: p_rblockDiscr
    
    bextended = .FALSE.
    
    ! Matrix must be 4x4.
    IF (rmatrix%ndiagBlocks .NE. 4) THEN
      CALL output_line ('System matrix is not 4x4.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
      CALL sys_halt()
    END IF
    
    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(4,:) must be (virtually) transposed. All matrices must be double precision.
    DO i=1,4
      DO j=1,4
      
        IF (lsysbl_isSubmatrixPresent(rmatrix,i,j)) THEN
        
          IF (i .LE. 3) THEN
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .NE. 0) THEN
              CALL output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
              CALL sys_halt()
            END IF
          ELSE
            IF ((i .LE. 3) .AND. &
               (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .EQ. 0)) THEN
              CALL output_line ('B1/B2/B3 submatrices must be virtually',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
              CALL output_line ('transposed (LSYSSC_MSPEC_TRANSPOSED)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
              CALL sys_halt()
            END IF
          END IF
          
          IF ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX7) .AND. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
            CALL output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
            CALL sys_halt()
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            CALL output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
            CALL sys_halt()
          END IF

          ! For scaled matrices, we have to use an extended sub-version of VANCA.
          IF ((rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) .AND. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 0.0_DP)) THEN
            bextended = .TRUE.  
          END IF
          
          IF ((i .eq. j) .AND. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) ) THEN
            bextended = .TRUE. 
          END IF
          
        END IF ! neq != 0
      END DO
    END DO
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    IF ((rmatrix%RmatrixBlock(1,4)%NA .NE. rmatrix%RmatrixBlock(4,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,4)%NEQ .NE. rmatrix%RmatrixBlock(4,1)%NCOLS)) THEN
      CALL output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
      CALL sys_halt()
    END IF

    IF ((rmatrix%RmatrixBlock(2,4)%NA .NE. rmatrix%RmatrixBlock(4,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,4)%NEQ .NE. rmatrix%RmatrixBlock(4,2)%NCOLS)) THEN
      CALL output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
      CALL sys_halt()
    END IF      
  
    IF ((rmatrix%RmatrixBlock(3,4)%NA .NE. rmatrix%RmatrixBlock(4,3)%NA) .OR. &
        (rmatrix%RmatrixBlock(3,4)%NEQ .NE. rmatrix%RmatrixBlock(4,3)%NCOLS)) THEN
      CALL output_line ('Structure of B3 and B3^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
      CALL sys_halt()
    END IF      

    ! Fill the output structure with data of the matrices.
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanca%rvanca3DNavSt%p_DA )
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,4),&
        rvanca%rvanca3DNavSt%p_DB1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,4),&
        rvanca%rvanca3DNavSt%p_DB2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,4),&
        rvanca%rvanca3DNavSt%p_DB3)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,1),&
        rvanca%rvanca3DNavSt%p_DD1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,2),&
        rvanca%rvanca3DNavSt%p_DD2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,3),&
        rvanca%rvanca3DNavSt%p_DD3)
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,4),&
        rvanca%rvanca3DNavSt%p_KcolB)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,4), &
        rvanca%rvanca3DNavSt%p_KldB )
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanca%rvanca3DNavSt%p_KcolA)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanca%rvanca3DNavSt%p_KldA )
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanca%rvanca3DNavSt%p_KdiagonalA)
    ELSE
      rvanca%rvanca3DNavSt%p_KdiagonalA => rvanca%rvanca3DNavSt%p_KldA
    END IF
    
    IF (lsysbl_isSubmatrixPresent(rmatrix,4,4)) THEN
    
      ! The matrix must be of format 7 or 9.
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,4),&
          rvanca%rvanca3DNavSt%p_DA44 )

      IF (rmatrix%RmatrixBlock(4,4)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(4,4), &
                                rvanca%rvanca3DNavSt%p_KdiagonalA44)
      ELSE
        CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,4), &
                                rvanca%rvanca3DNavSt%p_KdiagonalA44)
      END IF

      ! The presence of A(4,4) forces the extended VANCA to be used
      bextended = .TRUE.

    END IF
    
    IF (bextended) rvanca%csubsubtype = 1
    
    ! What is with A22? Is it the same as A11?
    IF (.NOT. lsyssc_isMatrixContentShared (&
        rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2)) ) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
          rvanca%rvanca3DNavSt%p_DA22 )
    END IF
    
    ! What is with A33? Is it the same as A11?
    IF (.NOT. lsyssc_isMatrixContentShared (&
        rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(3,3)) ) THEN
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
          rvanca%rvanca3DNavSt%p_DA33 )
    END IF

    ! What is with A12 and A21? Do they exist? With a scale factor = 1.0?
!    IF (lsysbl_isSubmatrixPresent(rmatrix,1,2) .AND. &
!        (rmatrix%RmatrixBlock(1,2)%dscaleFactor .EQ. 1.0_DP)) THEN
!      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
!          rvanca%rvanca2DNavSt%p_DA12 )
!          
!      CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,2),&
!          rvanca%rvanca2DNavSt%p_KcolA12)
!      CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,2), &
!          rvanca%rvanca2DNavSt%p_KldA12 )
!          
!      ! Get the structure. It's assumed that A12 and A21 have the same!
!      IF (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
!        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,2), &
!                                rvanca%rvanca2DNavSt%p_KdiagonalA12)
!      ELSE
!        rvanca%rvanca2DNavSt%p_KdiagonalA12 => rvanca%rvanca2DNavSt%p_KldA12
!      END IF
!      
!      IF (.NOT. lsysbl_isSubmatrixPresent(rmatrix,2,1)) THEN
!        CALL output_line ('If A12 is given, A21 must also be given!',&
!            OU_CLASS_ERROR,OU_MODE_STD,'vanca_init2DNavierStokesO')
!        CALL sys_halt()
!      END IF
!      
!      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
!          rvanca%rvanca2DNavSt%p_DA21 )
!    END IF

    ! Get the multiplication factors of the submatrices.
    rvanca%rvanca3DNavSt%Dmultipliers(1:4,1:4) = &
        rmatrix%RmatrixBlock(1:4,1:4)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest
    
    IF (.NOT. ASSOCIATED(p_rblockDiscr)) THEN
      CALL output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
      CALL sys_halt()
    END IF
    
    ! Get the discretisation structure of U,V,W and P from the block
    ! discretisation structure.
    rvanca%rvanca3DNavSt%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanca%rvanca3DNavSt%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanca%rvanca3DNavSt%p_rspatialDiscrW => p_rblockDiscr%RspatialDiscr(3)
    rvanca%rvanca3DNavSt%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(4)
    
    IF ((rvanca%rvanca3DNavSt%p_rspatialDiscrU%inumFESpaces .NE. &
         rvanca%rvanca3DNavSt%p_rspatialDiscrV%inumFESpaces) .OR. &
        (rvanca%rvanca3DNavSt%p_rspatialDiscrU%inumFESpaces .NE. &
         rvanca%rvanca3DNavSt%p_rspatialDiscrW%inumFESpaces)) THEN
      CALL output_line (&
          'Discretisation structures of X-, Y- and Z-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
      CALL sys_halt()
    END IF

    IF ((rvanca%rvanca3DNavSt%p_rspatialDiscrP%inumFESpaces .NE. 1) .AND. &
        (rvanca%rvanca3DNavSt%p_rspatialDiscrP%inumFESpaces .NE. &
         rvanca%rvanca3DNavSt%p_rspatialDiscrU%inumFESpaces)) THEN
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      CALL output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DNavierStokes')
      CALL sys_halt()
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_3DNavierStokes (rvanca3DNavSt, rvector, rrhs, domega, &
      csubtype, csubsubtype)
  
!<description>
  ! This routine applies the VANCA variant for 3D Navier-Stokes problems
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANCA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector

  ! The subtype of VANCA that should handle the above problem class.
  ! One of the VANCATP_xxxx constants, e.g. VANCATP_DIAGONAL.
  INTEGER :: csubtype

  ! The sub-subtype of VANCA that should handle the above problem class.
  ! =0: use standard VANCA. =1: use extended VANCA (e.g. with different 
  !     multipliers in the matrices)
  INTEGER :: csubsubtype
  
!</input>

!<inputoutput>
  ! t_vanca structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer3DNavSt), INTENT(INOUT) :: rvanca3DNavSt
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ielementdist
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementList
    TYPE(t_elementDistribution), POINTER :: p_relementDistrU
    TYPE(t_elementDistribution), POINTER :: p_relementDistrV
    TYPE(t_elementDistribution), POINTER :: p_relementDistrW
    TYPE(t_elementDistribution), POINTER :: p_relementDistrP 
    
    ! 3D Navier Stokes problem.

    ! Loop through the element distributions of the velocity.
    DO ielementdist = 1,rvanca3DNavSt%p_rspatialDiscrU%inumFESpaces
    
      ! Get the corresponding element distributions of U, V, W and P.
      p_relementDistrU => &
          rvanca3DNavSt%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanca3DNavSt%p_rspatialDiscrV%RelementDistr(ielementdist)
      p_relementDistrW => &
          rvanca3DNavSt%p_rspatialDiscrW%RelementDistr(ielementdist)
      
      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      IF (rvanca3DNavSt%p_rspatialDiscrP%inumFESpaces .GT. 1) THEN
        p_relementDistrP => &
            rvanca3DNavSt%p_rspatialDiscrP%RelementDistr(ielementdist)
      ELSE
        p_relementDistrP => &
            rvanca3DNavSt%p_rspatialDiscrP%RelementDistr(1)
      END IF
      
      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y- and Z-velocity (and to that
      ! of the pressure).
      CALL storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)
      
      ! Which element combination do we have now?
      IF ((elem_getPrimaryElement(p_relementDistrU%celement) .EQ. EL_Q1T_3D) .AND. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .EQ. EL_Q1T_3D) .AND. &
          (elem_getPrimaryElement(p_relementDistrW%celement) .EQ. EL_Q1T_3D) .AND. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .EQ. EL_Q0_3D)) THEN
        ! Q1~/Q1~/Q1~/Q0 discretisation
        
        ! Which VANCA subtype do we have? The diagonal VANCA of the full VANCA?
        SELECT CASE (csubtype)
        CASE (VANCATP_DIAGONAL)
!          ! Diagonal VANCA; check if A12 exists.
!          IF (.NOT. ASSOCIATED(rvanca3DNavSt%p_DA12)) THEN
            ! Call the VANCA subsolver to apply VANCA to our current element list.
            IF (rvanca3DNavSt%p_rspatialDiscrU%inumFESpaces .EQ. 1) THEN
              ! Uniform discretisation
              CALL vanca_3DSPQ1TQ0simple (rvanca3DNavSt, &
                  rvector, rrhs, domega)
            ELSE
              ! Conformal discretisation
              CALL vanca_3DSPQ1TQ0simpleConf (rvanca3DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
!          ELSE
!            ! Apply the conformal VANCA that allows different matrices
!            ! in A11, A12, A21 and A22!
!            CALL vanca_2DSPQ1TQ0simpleCoupConf (rvanca3DNavSt, &
!                  rvector, rrhs, domega,p_IelementList)
!          END IF
          
        CASE (VANCATP_FULL)
!          ! Full VANCA; check if A12 exists.
!          IF (.NOT. ASSOCIATED(rvanca3DNavSt%p_DA12)) THEN
            ! Call the VANCA subsolver to apply VANCA to our current element list.
            ! Note: Up to now, there is no 'full' variant -- has to be implemented!
            IF (rvanca3DNavSt%p_rspatialDiscrU%inumFESpaces .EQ. 1) THEN
              ! uniform discretisation;
              ! here, use the same as for the general conformal discretisation.
              ! Could be speeded up by introducing another variant...
              CALL vanca_3DSPQ1TQ0fullConf (rvanca3DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            ELSE
              ! Conformal discretisation
              CALL vanca_3DSPQ1TQ0fullConf (rvanca3DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            END IF
!          ELSE
!            ! Apply the conformal VANCA that allows different matrices
!            ! in A11, A12, A21 and A22!
!            ! If we have multiplication factors, we even have to use an extended
!            ! version of this.
!            IF (csubsubtype .EQ. 0) THEN
!              CALL vanca_2DSPQ1TQ0fullCoupConf (rvanca3DNavSt, &
!                  rvector, rrhs, domega,p_IelementList)
!            ELSE
!              CALL vanca_2DNSQ1TQ0fullCoupConfExt (rvanca3DNavSt, &
!                  rvector, rrhs, domega,p_IelementList)
!            END IF
!          END IF
        
        CASE DEFAULT
          CALL output_line (&
              'Unknown VANCA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanca_3DNavierStokes')
          CALL sys_halt()  
        
        END SELECT
        
      ELSE
        CALL output_line (&
            'Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_3DNavierStokes')
        CALL sys_halt()
        
      END IF
      
    END DO
      
  END SUBROUTINE
  
  ! ***************************************************************************
  ! 3D Navier-Stokes VANCA, simple diagonal version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    / A              B1 \
  !    |      A         B2 |
  !    |           A    B3 |
  !    \ D1^T D2^T D3^T 0  /
  !
  ! with D1/D2/D3 having the same structure as B1/B2/D3 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1, B2 and B3 are the same matrices as D1, D2 and D3. The only
  ! difference: Some rows in B1/B2/B3 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_init3DSPQ1TQ0simple (rmatrix,rvanca)
  
!<description>
  ! Checks if the "3D-Navier-Stokes-Q1T-Q0" VANCA variant can be applied to
  ! the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>

!<output>
  ! t_vancaPointer3DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer3DNavSt), INTENT(OUT) :: rvanca
!</output>

!</subroutine>

    INTEGER :: i,j

    ! Matrix must be 4x4.
    IF (rmatrix%ndiagBlocks .NE. 4) THEN
      CALL output_line (&
          'System matrix is not 4x4.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_check3DSPQ1TQ0')
      CALL sys_halt()
    END IF
    
    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(4,:) must be (virtually) transposed. All matrices must be double precision.
    DO i=1,4
      DO j=1,4
      
        IF (lsysbl_isSubmatrixPresent(rmatrix,i,j)) THEN
        
          IF (i .LE. 3) THEN
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .NE. 0) THEN
              CALL output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
              CALL sys_halt()
            END IF
          ELSE
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .EQ. 0) THEN
              CALL output_line ('B1/B2/B3 submatrices must be virtually '//&
                  'transposed (LSYSSC_MSPEC_TRANSPOSED)',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
              CALL sys_halt()
            END IF
          END IF
          
          IF ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX7) .AND. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
            CALL output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
            CALL sys_halt()
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            CALL output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
            CALL sys_halt()
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) THEN
            CALL output_line ('Scaled matrices not supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
            CALL sys_halt()
          END IF
          
        END IF ! neq != 0
      END DO
    END DO
    
    ! The structure of A(1,4) must be identical to A(4,1) and
    ! that of A(2,4) must be identical to A(4,2).
    IF ((rmatrix%RmatrixBlock(1,4)%NA .NE. rmatrix%RmatrixBlock(4,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,4)%NEQ .NE. rmatrix%RmatrixBlock(4,1)%NCOLS)) THEN
      CALL output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
      CALL sys_halt()
    END IF

    IF ((rmatrix%RmatrixBlock(2,4)%NA .NE. rmatrix%RmatrixBlock(4,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,4)%NEQ .NE. rmatrix%RmatrixBlock(4,2)%NCOLS)) THEN
      CALL output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
      CALL sys_halt()
    END IF
    
    IF ((rmatrix%RmatrixBlock(3,4)%NA .NE. rmatrix%RmatrixBlock(4,3)%NA) .OR. &
        (rmatrix%RmatrixBlock(3,4)%NEQ .NE. rmatrix%RmatrixBlock(4,3)%NCOLS)) THEN
      CALL output_line ('Structure of B3 and B3^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanca_init3DSPQ1TQ0simple')
      CALL sys_halt()
    END IF

    ! Fill the output structure with data of the matrices.
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),rvanca%p_DA )
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,4),rvanca%p_DB1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,4),rvanca%p_DB2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,4),rvanca%p_DB3)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,1),rvanca%p_DD1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,2),rvanca%p_DD2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(4,3),rvanca%p_DD3)
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,4),rvanca%p_KcolB)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,4), rvanca%p_KldB )
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),rvanca%p_KcolA)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), rvanca%p_KldA )
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                               rvanca%p_KdiagonalA)
    ELSE
      rvanca%p_KdiagonalA => rvanca%p_KldA
    END IF
    
    ! Get the multiplication factors of the submatrices
    rvanca%Dmultipliers(:,:) = rmatrix%RmatrixBlock(1:4,1:4)%dscaleFactor

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_3DSPQ1TQ0simple (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 3D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vancaPointer3DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer3DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DB3
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD3
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_FACEIDX), DIMENSION(:,:), POINTER :: p_IfacesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetw,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 6
    INTEGER, PARAMETER :: lofsw = 12
    INTEGER, PARAMETER :: lofsp = 18
    REAL(DP) :: daux
    !REAL(DP), DIMENSION(6,6) :: Dmult
    
    ! Local arrays for informations about one element
    REAL(DP), DIMENSION(6) :: AA,BB1,BB2,BB3,DD1,DD2,DD3
    REAL(DP), DIMENSION(19) :: FF,UU
    INTEGER(PREC_VECIDX), DIMENSION(6) :: idofGlobal
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DB3 => rvanca%p_DB3
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    p_DD3 => rvanca%p_DD3
    
    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanca%Dmultipliers(:,:)
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd, 3rd and 4th component of the solution
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetw = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetp = ioffsetw+rvector%RvectorBlock(3)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================
    !
    ! If you need a describtion with fancy ASCII arts abou the algorithm
    ! that is implemented below, scroll up to the 2D version.

    ! So we start with a loop over all elements
    DO iel=1,NEL
    
      ! Fetch the pressure P on the current element into FFP
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 6 U-nodes of that element.
      DO inode=1,6
      
        ! Set idof to the DOF that belongs to our face node:
        idof = p_IfacesAtElement(inode,iel)

        ! Write the number of the face node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        
        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANCA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let's hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        FF(inode+lofsw) = p_Drhs(idof+ioffsetw)
        
        ! At first build: fi = fi-Aui
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)
          
          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          !FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
          FF(inode+lofsw) = FF(inode+lofsw)-p_DB3(ib)*daux
          
          ! For support of scaled matrices, use the following lines; currently switched off.
          !FF(inode)       = FF(inode)      -Dmult(1,4)*p_DB1(ib)*daux
          !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,4)*p_DB2(ib)*daux
          !FF(inode+lofsw) = FF(inode+lofsw)-Dmult(3,4)*p_DB3(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local BX and DX.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        ! Either two (if the velocity DOF is an face with two neighbouring
        ! elements) or one (if the velocity DOF is at an face on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
          
            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            BB3(inode) = p_DB3(ib)
            
            ! The same way, get DD1, DD2 and DD3.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            DD3(inode) = p_DD3(ib)
            
            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,4)*p_DB1(ib)
            !BB2(inode) = Dmult(2,4)*p_DB2(ib)
            !BB3(inode) = Dmult(3,4)*p_DB3(ib)
            !DD1(inode) = Dmult(4,1)*p_DD1(ib)
            !DD2(inode) = Dmult(4,2)*p_DD2(ib)
            !DD3(inode) = Dmult(4,3)*p_DD3(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv) &
                        - DD3(inode)*p_Dvector(idof+ioffsetw)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_3DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,BB3,DD1,DD2,DD3)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,6
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
        p_Dvector(idofGlobal(inode)+ioffsetw) &
          = p_Dvector(idofGlobal(inode)+ioffsetw) + domega * UU(inode+lofsw)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_3DSPQ1TQ0simpleConf (rvanca, rvector, rrhs, domega,IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 3D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_3DSPQ1TQ0simpleConf is the same as vanca_3DSPQ1TQ0simple except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANCA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vancaPointer3DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer3DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
  
  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DB3
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD3
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_FACEIDX), DIMENSION(:,:), POINTER :: p_IfacesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetw,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 6
    INTEGER, PARAMETER :: lofsw = 12
    INTEGER, PARAMETER :: lofsp = 18
    REAL(DP) :: daux
    
    ! Local arrays for informations about one element
    REAL(DP), DIMENSION(6) :: AA,BB1,BB2,BB3,DD1,DD2,DD3
    REAL(DP), DIMENSION(19) :: FF,UU
    INTEGER(PREC_VECIDX), DIMENSION(6) :: idofGlobal
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DB3 => rvanca%p_DB3
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    p_DD3 => rvanca%p_DD3
    
    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanca%Dmultipliers(:,:)
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd, 3rd and 4th component of the solution
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetw = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetp = ioffsetw+rvector%RvectorBlock(3)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================
    !
    ! If you need a describtion with fancy ASCII arts abou the algorithm
    ! that is implemented below, scroll up to the 2D version.

    ! So we start with a loop over all elements
    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Fetch the pressure P on the current element into FFP
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 6 U-nodes of that element.
      DO inode=1,6
      
        ! Set idof to the DOF that belongs to our face node:
        idof = p_IfacesAtElement(inode,iel)

        ! Write the number of the face node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        
        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANCA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let's hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        FF(inode+lofsw) = p_Drhs(idof+ioffsetw)
        
        ! At first build: fi = fi-Aui
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)
          
          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          !FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
          FF(inode+lofsw) = FF(inode+lofsw)-p_DB3(ib)*daux
          
          ! For support of scaled matrices, use the following lines; currently switched off.
          !FF(inode)       = FF(inode)      -Dmult(1,4)*p_DB1(ib)*daux
          !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,4)*p_DB2(ib)*daux
          !FF(inode+lofsw) = FF(inode+lofsw)-Dmult(3,4)*p_DB3(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local BX and DX.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        ! Either two (if the velocity DOF is an face with two neighbouring
        ! elements) or one (if the velocity DOF is at an face on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
          
            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            BB3(inode) = p_DB3(ib)
            
            ! The same way, get DD1, DD2 and DD3.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            DD3(inode) = p_DD3(ib)
            
            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,4)*p_DB1(ib)
            !BB2(inode) = Dmult(2,4)*p_DB2(ib)
            !BB3(inode) = Dmult(3,4)*p_DB3(ib)
            !DD1(inode) = Dmult(4,1)*p_DD1(ib)
            !DD2(inode) = Dmult(4,2)*p_DD2(ib)
            !DD3(inode) = Dmult(4,3)*p_DD3(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv) &
                        - DD3(inode)*p_Dvector(idof+ioffsetw)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
          
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_3DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,BB3,DD1,DD2,DD3)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,6
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
        p_Dvector(idofGlobal(inode)+ioffsetw) &
          = p_Dvector(idofGlobal(inode)+ioffsetw) + domega * UU(inode+lofsw)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ************************************************************************

!<subroutine>

  PURE SUBROUTINE vanca_getcorr_3DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,BB3,DD1,DD2,DD3)
  
!<description>
  ! This routine solves a 19x19 Jacobi-type Schur complement system for three 
  ! velocity vectors and one pressure vector. It's used as auxiliary 
  ! routine in the simple VANCA solver to calculate an update vector
  ! for velocity/pressure.
!</description>

!<input>
  ! Diagonal elements of the local system matrix.
  REAL(DP), DIMENSION(6), INTENT(IN) :: AA
  
  ! Entries in the submatrix B1.
  REAL(DP), DIMENSION(6), INTENT(IN) :: BB1

  ! Entries in the submatrix B2.
  REAL(DP), DIMENSION(6), INTENT(IN) :: BB2

  ! Entries in the submatrix B3.
  REAL(DP), DIMENSION(6), INTENT(IN) :: BB3
  
  ! Entries in the submatrix D1.
  REAL(DP), DIMENSION(6), INTENT(IN) :: DD1

  ! Entries in the submatrix D2.
  REAL(DP), DIMENSION(6), INTENT(IN) :: DD2

  ! Entries in the submatrix D3.
  REAL(DP), DIMENSION(6), INTENT(IN) :: DD3

  ! Local RHS vector; FF(1..6)=X-velocity, FF(7..12)=Y-velocity,
  ! FF(13..18)=Y-velocity, FF(19)=pressure.
  REAL(DP), DIMENSION(19), INTENT(IN) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..6)=X-velocity, UU(7..12)=Y-velocity,
  ! UU(13..18)=Y-velocity, UU(19)=pressure.
  REAL(DP), DIMENSION(19), INTENT(OUT) :: UU
!</output>

!</subroutine>

    ! local variables

    INTEGER :: inode
    REAL(DP) :: PP,dpres
    REAL(DP), DIMENSION(19) :: AI,dff
    
    INTEGER, PARAMETER :: lofsv = 6
    INTEGER, PARAMETER :: lofsw = 12
    INTEGER, PARAMETER :: lofsp = 18

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the 
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial! 
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity 
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    DO inode=1,6
    
      ! Quick check if everything is ok - we don't want to divide by 0.
      IF (AA(inode)*AA(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      ! AI(.) saves the diagonal matrix S^-1:
      AI(inode)=1.0_DP/AA(inode)
        
    END DO

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 6x1 matrix 
    ! b) S^-1 is a diagonal matrix, given by the 6 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = 0.0_DP
    dff = FF
      
    DO inode = 1,6
      dpres        = dpres - AI(inode)*( &
                     DD1(inode)*BB1(inode)+ &
                     DD2(inode)*BB2(inode)+ &
                     DD3(inode)*BB3(inode))
      dff(1+lofsp) = dff(1+lofsp) - AI(inode)*( &
                     DD1(inode)*dff(inode)+ &
                     DD2(inode)*dff(inode+lofsv)+ &
                     DD3(inode)*dff(inode+lofsw))
    END DO

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF's are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    IF (dpres*dpres .LT. 1E-20_DP)  THEN
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      RETURN
    ENDIF
      
    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:
      
    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP
      
    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:
      
    DO inode=1,6
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
      UU(inode+lofsw) = AI(inode)*(dff(inode+lofsw)-BB3(inode)*PP)
    END DO

  END SUBROUTINE

  ! ************************************************************************

!<subroutine>

  PURE SUBROUTINE vanca_getcorr_3DSPQ1TQ0simple2 (UU,FF,AA1,AA2,AA3,&
                                          BB1,BB2,BB3,DD1,DD2,DD3,di)
  
!<description>
  ! This routine solves a 19x19 Jacobi-type Schur complement system for three 
  ! velocity vectors and one pressure vector. It's used as auxiliary 
  ! routine in the simple VANCA solver to calculate an update vector
  ! for velocity/pressure.
  !
  ! In contrast to vanca_getcorr_3DSPQ1TQ0simple, the three diagonal blocks
  ! in the matrix may be different from each other.
!</description>

!<input>
  ! Diagonal elements of the local system matrix A11
  REAL(DP), DIMENSION(6), INTENT(IN) :: AA1

  ! Diagonal elements of the local system matrix A22
  REAL(DP), DIMENSION(6), INTENT(IN) :: AA2
  
  ! Diagonal elements of the local system matrix A33
  REAL(DP), DIMENSION(6), INTENT(IN) :: AA3

  ! Entries in the submatrix B1.
  REAL(DP), DIMENSION(6), INTENT(IN) :: BB1

  ! Entries in the submatrix B2.
  REAL(DP), DIMENSION(6), INTENT(IN) :: BB2
  
  ! Entries in the submatrix B3.
  REAL(DP), DIMENSION(6), INTENT(IN) :: BB3

  ! Entries in the submatrix D1.
  REAL(DP), DIMENSION(6), INTENT(IN) :: DD1

  ! Entries in the submatrix D2.
  REAL(DP), DIMENSION(6), INTENT(IN) :: DD2

  ! Entries in the submatrix D3.
  REAL(DP), DIMENSION(6), INTENT(IN) :: DD3
  
  ! Entry in the submatrix I (usually =0).
  ! This is the matrix in the diagonal block of the pressure, which is usually
  ! zero in saddle point problems.
  REAL(DP), INTENT(IN) :: di
  
  ! Local RHS vector; FF(1..6)=X-velocity, FF(7..12)=Y-velocity,
  ! FF(13..18)=Y-velocity, FF(19)=pressure.
  REAL(DP), DIMENSION(19), INTENT(IN) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..6)=X-velocity, UU(7..12)=Y-velocity,
  ! UU(13..18)=Y-velocity, UU(19)=pressure.
  REAL(DP), DIMENSION(19), INTENT(OUT) :: UU
!</output>

!</subroutine>

    ! local variables

    INTEGER :: inode
    REAL(DP) :: PP,dpres
    REAL(DP), DIMENSION(19) :: AI1,AI2,AI3,dff
    
    INTEGER, PARAMETER :: lofsv = 6
    INTEGER, PARAMETER :: lofsw = 12
    INTEGER, PARAMETER :: lofsp = 18

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T i )
    !      ( D1  D2  I  )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the 
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t I ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !       D^t u + Ip = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !      D^t S^-1 (f-Bp) + Ip = g
    !
    !      <=>  Ip -D^t S^-1 B p  =  g - D^t S^-1 f
    !              ***********       **************
    !                 =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial! 
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity 
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    DO inode=1,6
    
      ! Quick check if everything is ok - we don't want to divide by 0.
      IF (AA1(inode)*AA1(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      IF (AA2(inode)*AA2(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      IF (AA3(inode)*AA3(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      ! AI(.) saves the diagonal matrix S^-1:
      AI1(inode)=1.0_DP/AA1(inode)
      AI2(inode)=1.0_DP/AA2(inode)
      AI3(inode)=1.0_DP/AA3(inode)
        
    END DO

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !      ( I - D^t S^-1 B ) p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 6x1 matrix 
    ! b) S^-1 is a diagonal matrix, given by the 6 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           = (I - D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = di
    dff = FF
      
    DO inode = 1,6
      dpres        = dpres &
                   - AI1(inode)*DD1(inode)*BB1(inode) &
                   - AI2(inode)*DD2(inode)*BB2(inode) &
                   - AI3(inode)*DD3(inode)*BB3(inode)
                   
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI1(inode)*DD1(inode)*dff(inode) &
                   - AI2(inode)*DD2(inode)*dff(inode+lofsv) &
                   - AI3(inode)*DD3(inode)*dff(inode+lofsw)
    END DO

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF's are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    IF (dpres*dpres .LT. 1E-20_DP)  THEN
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      RETURN
    ENDIF
      
    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:
      
    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP
      
    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:
      
    DO inode=1,6
      UU(inode)       = AI1(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI2(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
      UU(inode+lofsw) = AI3(inode)*(dff(inode+lofsw)-BB3(inode)*PP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  ! 3D Navier-Stokes VANCA, 'full' version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    / A              B1 \
  !    |      A         B2 |
  !    |           A    B3 |
  !    \ D1^T D2^T D3^T 0  /
  !
  ! with D1/D2/D3 having the same structure as B1/B2/B3 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1, B2 and B3 are the same matrices as D1, D2 and D3. The only
  ! difference: Some rows in B1/B2/B3 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_3DSPQ1TQ0fullConf (rvanca, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 3D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanca_3DSPQ2QP1fullConf is the same as vanca_3DSPQ2QP1full except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANCA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vancaPointer3DNavSt structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer3DNavSt), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! A list of element numbers where VANCA should be applied to.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:)     :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DB3
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD3
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_FACEIDX), DIMENSION(:,:), POINTER :: p_IfacesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 6      ! Q1T = 6 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 1 ! QQ0 = 1 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 3*nnvel+nnpressure
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetw,ioffsetp,j
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,k
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsw = 2*nnvel
    INTEGER, PARAMETER :: lofsp = 3*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DB3 => rvanca%p_DB3
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    p_DD3 => rvanca%p_DD3
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetw = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetp = ioffsetw+rvector%RvectorBlock(3)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! So we start with a loop over all elements in the list
    DO ielidx=1,SIZE(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..6 = face.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:6) = p_IfacesAtElement(1:6,iel)

      ! Loop over all U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        FF(inode+lofsw) = p_Drhs(idof+ioffsetw)
        
        ! At first build: fi = fi-Aui
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to all A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode,k) = daux
              AA (inode+lofsv,k+lofsv) = daux
              AA (inode+lofsw,k+lofsw) = daux
              EXIT
            END IF
          END DO          
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
          FF(inode+lofsw) = FF(inode+lofsw)-p_DB3(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        ! Either two (if the velocity DOF is an face with two neighbouring
        ! elements) or one (if the velocity DOF is at an face on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
            
            ! Get the entries in the B-matrices
            AA(inode,      lofsp+1) = p_DB1(ib)
            AA(inode+lofsv,lofsp+1) = p_DB2(ib)
            AA(inode+lofsw,lofsp+1) = p_DB3(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(lofsp+1,inode)       = p_DD1(ib)
            AA(lofsp+1,inode+lofsv) = p_DD2(ib)
            AA(lofsp+1,inode+lofsw) = p_DD3(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(lofsp+1,inode)*p_Dvector(idof) &
                        - AA(lofsp+1,inode+lofsv)*p_Dvector(idof+ioffsetv) &
                        - AA(lofsp+1,inode+lofsw)*p_Dvector(idof+ioffsetw)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)  ..............                                   :::::: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(4,4)                               :::::: )
      !     (                            AA( 5, 5) ..............       :::::: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA( 8, 8)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      
      IF (ilapackInfo .EQ. 0) THEN
      
        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!
        
        DO inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
          p_Dvector(idofGlobal(inode)+ioffsetw) &
            = p_Dvector(idofGlobal(inode)+ioffsetw) + domega * FF(inode+lofsw)
        END DO
        
        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
      ELSE IF (ilapackInfo .LT. 0) THEN

        CALL output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanca_3DSPQ1TQ0fullConf')

      END IF
      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!
    
    END DO ! iel

  END SUBROUTINE

END MODULE 
