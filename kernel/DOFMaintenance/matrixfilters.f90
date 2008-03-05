!##############################################################################
!# ****************************************************************************
!# <name> matrixfilters </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of filters that can be applied to a scalar
!# or block matrix. They can be used e.g. to modify a matrix during a nonlinear
!# solution process to prepare it for a linear solver (like implementing
!# boundary conditions).
!#
!# Typical matrix filters are those working together with discrete
!# boundary conditions. For Dirichlet boundary conditions e.g. some lines
!# of the global matrix are usually replaced by unit vectors. The corresponding
!# matrix modification routine is a kind of 'filter' for a matrix and
!# can be found in this module.
!#
!# The following routines can be found here:
!#
!#  1.) matfil_discreteBC
!#      -> Apply the 'discrete boundary conditions for matrices' filter
!#         onto a given (block) matrix.
!#         This e.g. replaces lines of a matrix corresponding to Dirichlet DOF's
!#         by unit vectors for all scalar submatrices where configured.
!#
!#  2.) matfil_discreteFBC
!#      -> Apply the 'discrete fictitious boundary conditions for matríces'
!#         filter onto a given (block) matrix.
!#         This e.g. replaces lines of a matrix corresponding to Dirichlet DOF's
!#         by unit vectors for all scalar submatrices where configured.
!#
!#  3.) matfil_discreteNLSlipBC
!#      -> Implement slip boundary conditions into a matrix.
!#         Slip boundary conditions are not implemented by matfil_discreteBC!
!#         They must be implemented separately, as they are a special type
!#         of boundary conditions.
!#
!#  4.) matfil_normaliseToL20
!#      -> Adds a row to a matrix containing 'ones'. This implements the
!#         explicit condition "int_{\Omega} v = 0" into a matrix.
!#
!# Auxiliary routines:
!#
!#  1.) matfil_imposeDirichletBC
!#      -> Imposes Dirichlet BC's into a scalar matrix
!#
!#  2.) matfil_imposeNLSlipBC
!#      -> Imposes nonlinear slip boundary conditions into a scalar matrix
!#
!#  3.) matfil_imposeDirichletFBC
!#      -> Imposes Difichlet BC's for fictitious boundary components
!#         into a scalar matrix
!#
!# </purpose>
!##############################################################################

MODULE matrixfilters

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE discretebc
  USE discretefbc
  USE dofmapping
  USE matrixmodification
  USE genoutput
  
  IMPLICIT NONE

CONTAINS

! *****************************************************************************
! Scalar matrix filters
! *****************************************************************************

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_imposeDirichletBC (rmatrix,boffDiag,rdbcStructure)
  
!<description>
  ! Implements discrete Dirichlet BC's into a scalar matrix.
  ! This is normally done by replacing some lines of the matrix rmatrix
  ! (those belonging to Dirichlet nodes) by unit vectors or by zero-vectors
  ! (depending on whether the matrix is a 'diagonal' matrix or an
  ! 'off-diagonal' matrix in a larger block-system).
!</description>

!<input>
  
  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteBCDirichlet), INTENT(IN), TARGET  :: rdbcStructure
  
  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it's assumed that the matrix is not
  ! a main, guiding system matrix, but an 'off-diagonal' matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  LOGICAL :: boffDiag

!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrix
  
!</inputoutput>
  
!</subroutine>
    
  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_idx
  INTEGER, PARAMETER :: NBLOCKSIZE = 1000
  INTEGER(PREC_VECIDX), DIMENSION(1000) :: Idofs
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  INTEGER i,ilenleft

  ! If nDOF=0, there are no DOF's the current boundary condition segment,
  ! so we don't have to do anything. Maybe the case if the user selected
  ! a boundary region that is just too small.
  IF (rdbcStructure%nDOF .EQ. 0) RETURN

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.
  
  CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF (.NOT.ASSOCIATED(p_idx)) THEN
    CALL output_line ('DBC not configured',&
        OU_CLASS_ERROR,OU_MODE_STD,'matfil_imposeDirichletBC')
    CALL sys_halt()
  END IF
  
  ! Only handle nDOF DOF's, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  !
  ! Is the matrix sorted?
  IF (rmatrix%isortStrategy .LE. 0) THEN
    ! Use mmod_replaceLinesByUnit/mmod_replaceLinesByZero to replace the 
    ! corresponding rows in the matrix by unit vectors. 
    ! For more complicated FE spaces, this might have
    ! to be modified in the future...
    !
    ! We use mmod_replaceLinesByUnit for 'diagonal' matrix blocks (main
    ! system matrices) and mmod_replaceLinesByZero for 'off-diagonal'
    ! matrix blocks (in case of larger block systems)
    
    IF (boffDiag) THEN
      CALL mmod_replaceLinesByZero (rmatrix,p_idx(1:rdbcStructure%nDOF))
    ELSE
      CALL mmod_replaceLinesByUnit (rmatrix,p_idx(1:rdbcStructure%nDOF))
    END IF
  ELSE
    ! Ok, matrix is sorted, so we have to filter all the DOF's through the
    ! permutation before using them for implementing boundary conditions.
    ! We do this in blocks with 1000 DOF's each to prevent the stack
    ! from being destroyed!
    !
    ! Get the permutation from the matrix - or more precisely, the
    ! back-permutation, as we need this one for the loop below.
    CALL storage_getbase_int (rmatrix%h_IsortPermutation,p_Iperm)
    p_Iperm => p_Iperm(rmatrix%NEQ+1:)
    
    ! How many entries to handle in the first block?
    ilenleft = MIN(rdbcStructure%nDOF,NBLOCKSIZE)
    
    ! Loop through the DOF-blocks
    DO i=0,rdbcStructure%nDOF / NBLOCKSIZE
      ! Filter the DOF's through the permutation
      Idofs(1:ilenleft) = p_Iperm(p_idx(1+i*NBLOCKSIZE:i*NBLOCKSIZE+ilenleft))
      
      ! And implement the BC's with mmod_replaceLinesByUnit/
      ! mmod_replaceLinesByZero, depending on whether the matrix is a
      ! 'main' system matrix or an 'off-diagonal' system matrix in a larger
      ! block system.
      IF (boffDiag) THEN
        CALL mmod_replaceLinesByZero (rmatrix,Idofs(1:ilenleft))
      ELSE
        CALL mmod_replaceLinesByUnit (rmatrix,Idofs(1:ilenleft))
      END IF

      ! How many DOF's are left?
      ilenleft = MIN(rdbcStructure%nDOF-(i+1)*NBLOCKSIZE,NBLOCKSIZE)
    END DO
  
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_imposeNLSlipBC (rmatrix,boffDiag,bforprec,rslipBCStructure)
  
!<description>
  ! Implements discrete Slip BC's into a scalar matrix.
  ! Slip BC's are treated like Dirichlet BC's.
  ! This is normally done by replacing some lines of the matrix rmatrix
  ! (those belonging to Dirichlet nodes) by unit vectors or by zero-vectors
  ! (depending on whether the matrix is a 'diagonal' matrix or an
  ! 'off-diagonal' matrix in a larger block-system).
!</description>

!<input>
  ! Prepare matrix for preconditioning.
  ! =false: Slip-nodes are handled as Dirichlet. Rows in the
  !         matrix are replaced by unit vectors.
  ! =true : Standard. Prepare matrix for preconditioning.
  !         Only the off-diagonals of rows corresponding to 
  !         the slip-nodes are cleared.
  !         The matrix therefore is prepared to act as 
  !         Jacobi-preconditioner for all slip-nodes.
  LOGICAL, INTENT(IN) :: bforprec
  
  ! The t_discreteBCSlip that describes the discrete Slip BC's
  TYPE(t_discreteBCSlip), INTENT(IN), TARGET  :: rslipBCStructure
  
  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it's assumed that the matrix is not
  ! a main, guiding system matrix, but an 'off-diagonal' matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  LOGICAL, INTENT(IN) :: boffDiag
!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrix
  
!</inputoutput>
  
!</subroutine>
    
  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_idx
  INTEGER, PARAMETER :: NBLOCKSIZE = 1000
  INTEGER(PREC_VECIDX), DIMENSION(1000) :: Idofs
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  INTEGER i,ilenleft
  
  ! If nDOF=0, there are no DOF's the current boundary condition segment,
  ! so we don't have to do anything. Maybe the case if the user selected
  ! a boundary region that is just too small.
  IF (rslipBCStructure%nDOF .EQ. 0) RETURN

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.
  
  CALL storage_getbase_int(rslipBCStructure%h_IslipDOFs,p_idx)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF (.NOT.ASSOCIATED(p_idx)) THEN
    CALL output_line ('DBC not configured',&
        OU_CLASS_ERROR,OU_MODE_STD,'matfil_imposeNLSlipBC')
    CALL sys_halt()
  END IF
  
  ! Only handle nDOF DOF's, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  !
  ! Is the matrix sorted?
  IF (rmatrix%isortStrategy .LE. 0) THEN
    ! Use mmod_clearOffdiags/mmod_replaceLinesByZero clear the
    ! offdiagonals of the system matrix.
    ! For more complicated FE spaces, this might have
    ! to be modified in the future...
    !
    ! We use mmod_clearOffdiags for 'diagonal' matrix blocks (main
    ! system matrices) and mmod_replaceLinesByZero for 'off-diagonal'
    ! matrix blocks (in case of larger block systems)
    
    IF (boffDiag) THEN
      CALL mmod_replaceLinesByZero (rmatrix,p_idx(1:rslipBCStructure%nDOF))
    ELSE
      IF (bforprec) THEN
        CALL mmod_clearOffdiags (rmatrix,p_idx(1:rslipBCStructure%nDOF))
      ELSE
        CALL mmod_replaceLinesByUnit (rmatrix,p_idx(1:rslipBCStructure%nDOF))
      END IF
    END IF
  ELSE
    ! Ok, matrix is sorted, so we have to filter all the DOF's through the
    ! permutation before using them for implementing boundary conditions.
    ! We do this in blocks with 1000 DOF's each to prevent the stack
    ! from being destroyed!
    !
    ! Get the permutation from the matrix - or more precisely, the
    ! back-permutation, as we need this one for the loop below.
    CALL storage_getbase_int (rmatrix%h_IsortPermutation,p_Iperm)
    p_Iperm => p_Iperm(rmatrix%NEQ+1:)
    
    ! How many entries to handle in the first block?
    ilenleft = MIN(rslipBCStructure%nDOF,NBLOCKSIZE)
    
    ! Loop through the DOF-blocks
    DO i=0,rslipBCStructure%nDOF / NBLOCKSIZE
      ! Filter the DOF's through the permutation
      Idofs(1:ilenleft) = p_Iperm(p_idx(1+i*NBLOCKSIZE:i*NBLOCKSIZE+ilenleft))
      
      ! And implement the BC's with mmod_clearOffdiags/
      ! mmod_replaceLinesByZero, depending on whether the matrix is a
      ! 'main' system matrix or an 'off-diagonal' system matrix in a larger
      ! block system.
      IF (boffDiag) THEN
        CALL mmod_replaceLinesByZero (rmatrix,Idofs(1:ilenleft))
      ELSE
        IF (bforprec) THEN
          CALL mmod_clearOffdiags (rmatrix,Idofs(1:ilenleft))
        ELSE
          CALL mmod_replaceLinesByUnit (rmatrix,Idofs(1:ilenleft))
        END IF
      END IF
      
      ! How many DOF's are left?
      ilenleft = MIN(rslipBCStructure%nDOF-(i+1)*NBLOCKSIZE,NBLOCKSIZE)
    END DO
  
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_imposeDirichletFBC (rmatrix,boffDiag,rdbcStructure)
  
!<description>
  ! Implements discrete Dirichlet BC's of fictitious boundary components
  ! into a scalar matrix.
  ! This is normally done by replacing some lines of the matrix rmatrix
  ! (those belonging to Dirichlet nodes) by unit vectors or by zero-vectors
  ! (depending on whether the matrix is a 'diagonal' matrix or an
  ! 'off-diagonal' matrix in a larger block-system).
!</description>

!<input>
  
  ! The t_discreteFBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteFBCDirichlet), INTENT(IN), TARGET  :: rdbcStructure
  
  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it's assumed that the matrix is not
  ! a main, guiding system matrix, but an 'off-diagonal' matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  LOGICAL :: boffDiag

!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrix
  
!</inputoutput>
  
!</subroutine>
    
  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_idx
  INTEGER, PARAMETER :: NBLOCKSIZE = 1000
  INTEGER(PREC_VECIDX), DIMENSION(1000) :: Idofs
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  INTEGER i,ilenleft

  ! If nDOF=0, there are no DOF's the current boundary condition segment,
  ! so we don't have to do anything. Maybe the case if the user selected
  ! a boundary region that is just too small.
  IF (rdbcStructure%nDOF .EQ. 0) RETURN

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.
  
  CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF (.NOT.ASSOCIATED(p_idx)) THEN
    CALL output_line ('DBC not configured',&
        OU_CLASS_ERROR,OU_MODE_STD,'matfil_imposeDirichletFBC')
    CALL sys_halt()
  END IF
  
  ! Only handle nDOF DOF's, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  !
  ! Is the matrix sorted?
  IF (rmatrix%isortStrategy .LE. 0) THEN
    ! Use mmod_replaceLinesByUnit/mmod_replaceLinesByZero to replace the 
    ! corresponding rows in the matrix by unit vectors. 
    ! For more complicated FE spaces, this might have
    ! to be modified in the future...
    !
    ! We use mmod_replaceLinesByUnit for 'diagonal' matrix blocks (main
    ! system matrices) and mmod_replaceLinesByZero for 'off-diagonal'
    ! matrix blocks (in case of larger block systems)
    
    IF (boffDiag) THEN
      CALL mmod_replaceLinesByZero (rmatrix,p_idx(1:rdbcStructure%nDOF))
    ELSE
      CALL mmod_replaceLinesByUnit (rmatrix,p_idx(1:rdbcStructure%nDOF))
    END IF
  ELSE
    ! Ok, matrix is sorted, so we have to filter all the DOF's through the
    ! permutation before using them for implementing boundary conditions.
    ! We do this in blocks with 1000 DOF's each to prevent the stack
    ! from being destroyed!
    !
    ! Get the permutation from the matrix - or more precisely, the
    ! back-permutation, as we need this one for the loop below.
    CALL storage_getbase_int (rmatrix%h_IsortPermutation,p_Iperm)
    p_Iperm => p_Iperm(rmatrix%NEQ+1:)
    
    ! How many entries to handle in the first block?
    ilenleft = MIN(rdbcStructure%nDOF,NBLOCKSIZE)
    
    ! Loop through the DOF-blocks
    DO i=0,rdbcStructure%nDOF / NBLOCKSIZE
      ! Filter the DOF's through the permutation
      Idofs(1:ilenleft) = p_Iperm(p_idx(1+i*NBLOCKSIZE:i*NBLOCKSIZE+ilenleft))
      
      ! And implement the BC's with mmod_replaceLinesByUnit/
      ! mmod_replaceLinesByZero, depending on whether the matrix is a
      ! 'main' system matrix or an 'off-diagonal' system matrix in a larger
      ! block system.
      IF (boffDiag) THEN
        CALL mmod_replaceLinesByZero (rmatrix,Idofs(1:ilenleft))
      ELSE
        CALL mmod_replaceLinesByUnit (rmatrix,Idofs(1:ilenleft))
      END IF
    END DO
  
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_normaliseToL20 (rmatrix,istartColumn,iendColumn)
  
!<description>
  ! Modifies a scalar matrix to add an additional equation which implements
  ! $int_{\Omega) v = 0$ for all $v$. This adds another row to the matrix
  ! which sums up all vector entries when being applied to a vector.
  ! Note that this produces a rectangular matrix with NROW+1 rows in comparison
  ! to the original one!
  !
  ! This filter can be used e.g. to impose the condition $int_{\Omega) v = 0$
  ! to a matrix which is passed to an external linear solver (like UMFPACK)
  ! which cannot use filtering to cope with indefiniteness!
  !
  ! WARNING: UNTESTED!!!
!</description>

!<input>
  ! OPTIONAL: Start column in the matrix.
  ! This parameter can specify where to start the vector sum.
  ! If not specified, 1 is assumed.
  INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: istartColumn

  ! OPTIONAL: End column in the matrix.
  ! This parameter can specify where to end the vector sum.
  ! If not specified, rmatrix%NCOLS is assumed.
  INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: iendColumn
!</input>

!<inputoutput>
  ! The matrix which is to be modified.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>
  
!</subroutine>
    
    ! Local variables
    TYPE(t_matrixScalar) :: rmatrixTemp
    INTEGER(PREC_VECIDX) :: irowLen,istart,iend,i
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    
    ! Matrix must not be transposed
    IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
      CALL output_line (&
          'Virtually transposed matrices not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'matfil_normaliseToL20')
      CALL sys_halt()
    END IF

    ! Only double precision supported.
    IF (rmatrix%cdataType .NE. ST_DOUBLE) THEN
      CALL output_line (&
          'Only double precision matrices supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'matfil_normaliseToL20')
      CALL sys_halt()
    END IF
    
    ! If structure and/or content is shared, duplicate the matrix to
    ! make the entries belong to rmatrix.
    IF (lsyssc_isMatrixStructureShared(rmatrix) .AND. &
        lsyssc_isMatrixContentShared(rmatrix)) THEN
      rmatrixTemp = rmatrix
      CALL lsyssc_duplicateMatrix (rmatrixTemp,rmatrix,LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    ELSE IF (lsyssc_isMatrixStructureShared(rmatrix)) THEN
      rmatrixTemp = rmatrix
      CALL lsyssc_duplicateMatrix (rmatrixTemp,rmatrix,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
    ELSE IF (lsyssc_isMatrixContentShared(rmatrix)) THEN
      rmatrixTemp = rmatrix
      CALL lsyssc_duplicateMatrix (rmatrixTemp,rmatrix,LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
    END IF
    
    ! Length of the row
    istart = 1
    iend = rmatrix%NCOLS
    
    IF (PRESENT(istartColumn)) &
      istart = MAX(istart,MIN(istartColumn,rmatrix%NCOLS))

    IF (PRESENT(iendColumn)) &
      iend = MAX(istart,MIN(iendColumn,rmatrix%NCOLS))
      
    irowLen = iend-istart+1
    
    ! Add another row to the matrix
    CALL lsyssc_resizeMatrixDirect (rmatrix, rmatrix%NEQ+1_PREC_VECIDX, &
        rmatrix%NCOLS, rmatrix%NA+irowLen, .FALSE.)
        
    rmatrix%NEQ = rmatrix%NEQ + 1
    rmatrix%NA = rmatrix%NA+irowLen
    
    ! Add a row containing "1". TOgether with a RHS "0",
    ! this implements "sum_j a_ij v_j = 0".
    CALL lsyssc_getbase_double (rmatrix,p_Ddata)
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX1)
      istart = istart + rmatrix%NEQ * rmatrix%NCOLS
      iend = iend + rmatrix%NEQ * rmatrix%NCOLS
      DO i=istart,iend
        p_Ddata(i) = 1.0_DP
      END DO
      
    CASE (LSYSSC_MATRIX7)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      
      ! New Kld entry
      p_Kld(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ) + irowLen
      
      ! Insert the new row
      DO i=p_Kld(rmatrix%NEQ),p_Kld(rmatrix%NEQ+1)-1
        p_Ddata(i) = 1.0_DP
        p_Kcol(i) = istart + (p_Kld(rmatrix%NEQ)-i)
      END DO
    
    CASE (LSYSSC_MATRIX9)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kdiagonal)
      
      ! New Kld entry
      p_Kld(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ) + irowLen
      
      ! Insert the new row
      DO i=p_Kld(rmatrix%NEQ),p_Kld(rmatrix%NEQ+1)-1
        p_Ddata(i) = 1.0_DP
        p_Kcol(i) = istart + (p_Kld(rmatrix%NEQ)-i)
      END DO
      
      ! New Kdiagonal entry
      p_Kdiagonal(rmatrix%NEQ) = p_Kld(rmatrix%NEQ) + rmatrix%NCOLS-istart
    
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_imposeFeastMirrorBC (rmatrix,boffDiag,rfmbcStructure)
  
!<description>
  ! Implements discrete Feast mirror BC's into a scalar matrix.
  ! The FEAST mirror boundary condition is basically a Neumann boundary
  ! condition which is used for domain decomposition.
  ! One assumes that there is an additional 'virtual' layer of cells added to
  ! a boundary edge. This leads to a slight matrix modification for all
  ! DOF's on that boundary edge. 
  ! Example: For a 5-point stencil with $Q_1$, boundary DOF's get matrix
  ! weights "2, 1, -1/2, -1/2" (center, left, top, bottom), while inner 
  ! points get matrix weights "4, -1, -1, -1, -1" (center and all surroundings).
  ! To make bondary DOF's behave like inner DOF's, the entries in 
  ! the matrices belonging to such an edge have to be doubled,
  ! leading to "4, -1, -1".
  ! So this filter loops through the matrix and doubles all matrix entries
  ! that belong to DOF's on FEAST mirror boundary edges.
!</description>

!<input>
  
  ! The t_discreteBCfeastMirror that describes the discrete FEAST mirror BC's
  TYPE(t_discreteBCfeastMirror), INTENT(IN), TARGET  :: rfmbcStructure
  
  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it's assumed that the matrix is not
  ! a main, guiding system matrix, but an 'off-diagonal' matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  LOGICAL :: boffDiag

!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrix
  
!</inputoutput>
  
!</subroutine>
    
  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_ImirrorDOFs,p_ImirrorDOFsClosed
  INTEGER :: i,j
  REAL(DP), DIMENSION(:), POINTER :: p_Da
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol,p_Iperm,p_IpermInverse
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  INTEGER(PREC_MATIDX) :: ia
  INTEGER(PREC_DOFIDX) :: idof
  REAL(DP) :: dmirrorWeight

  ! Offdiagonal matrices are not processed by this routine up to now.
  IF (boffDiag) RETURN

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF ((rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND. &
      (rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
    CALL output_line ('Only matrix format 7 and 9 supported!',&
        OU_CLASS_ERROR,OU_MODE_STD,'matfil_imposeFeastMirrorBC')
    CALL sys_halt()
  END IF
  
  IF (rmatrix%cdataType .NE. ST_DOUBLE) THEN
    CALL output_line ('Matrix must be double precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'matfil_imposeFeastMirrorBC')
    CALL sys_halt()
  END IF
  
  IF (rfmbcStructure%icomponent .EQ. 0) THEN
    CALL output_line ('FMBC not configured!',&
        OU_CLASS_ERROR,OU_MODE_STD,'matfil_imposeFeastMirrorBC')
    CALL sys_halt()
  END IF
  
  IF (rfmbcStructure%h_ImirrorDOFs .EQ. ST_NOHANDLE) THEN
    ! No data inside of this structure.
    ! May happen if the region is not large enough to cover at least one DOF.
    RETURN
  END IF
  
  ! Get the matrix data
  CALL lsyssc_getbase_double (rmatrix,p_Da)
  CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
  CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
  
  ! Get the weight of the entries.
  ! =2 on finest level, =1.5 on level NLMAX-1,...
  !dmirrorWeight = 1.0_DP+REAL(4**rfmbcStructure%icoarseningLevel,DP)
  dmirrorWeight = 1.0_DP+1.0_DP*REAL(2**rfmbcStructure%icoarseningLevel,DP)
  
  ! Get pointers to the list of DOF's that belong to that region and have
  ! to be tackled.
  ! p_ImirrorDOFs is a list of all DOF's in the region.
  ! p_ImirrorDOFsClosed is a list of all DOF's in the closure of the region.
  ! For every DOF in the region, it's neighbours have to be found in the
  ! clusure. If that's the case, the corresponding matrix entry has to be doubled.
  
  CALL storage_getbase_int(rfmbcStructure%h_ImirrorDOFs,p_ImirrorDOFs)
  CALL storage_getbase_int(rfmbcStructure%h_ImirrorDOFsClosed,p_ImirrorDOFsClosed)

  ! The matrix column corresponds to the DOF. For every DOF decide on
  ! whether it's on the FEAST mirror boundary component or not.
  ! If yes, double the matrix entry.
  
  ! Is the matrix sorted?
  IF (rmatrix%isortStrategy .LE. 0) THEN
    
    ! Loop through the DOF's. Each DOF gives us a matrix row to change.
    DO i=1,SIZE(p_ImirrorDOFs)
    
      ! Loop through the matrix row. All DOF's in that matrix row that
      ! belong to the closed region have to be changed.
      DO ia=p_Kld(p_ImirrorDOFs(i)),p_Kld(p_ImirrorDOFs(i)+1)-1
        ! Get the DOF.
        idof = p_Kcol(ia)
        
        ! Search the DOF in our list. Ok, that's an n^2 algorithm.
        ! It could easily replaced by an n log n algorithm using binary
        ! search since the list of DOF's is sorted!
        ! Probably in a later implementation...
        DO j=1,SIZE(p_ImirrorDOFsClosed)
          IF (p_ImirrorDOFsClosed(j) .EQ. idof) THEN
            p_Da(ia) = dmirrorWeight * p_Da(ia)
            EXIT
          END IF
        END DO
        
      END DO
      
    END DO
    
  ELSE
  
    ! Ok, matrix is sorted, so we have to filter all the DOF's through the
    ! permutation before using them for implementing boundary conditions.
    !
    ! Get the permutation/inverse permutation from the matrix to renumber the columns into
    ! the actual DOF numbers.
    CALL storage_getbase_int (rmatrix%h_IsortPermutation,p_Iperm)
    p_IpermInverse => p_Iperm(1:rmatrix%NEQ)
    p_Iperm => p_Iperm(rmatrix%NEQ+1:)
    
    ! Loop through the DOF's. Each DOF gives us a matrix row to change.
    DO i=1,SIZE(p_ImirrorDOFs)
    
      ! Loop through the matrix row. All DOF's in that matrix row that
      ! belong to our region have to be changed.
      DO ia=p_Kld(p_Iperm(p_ImirrorDOFs(i))),&
            p_Kld(p_Iperm(p_ImirrorDOFs(i))+1)-1
        ! Get the DOF.
        idof = p_IpermInverse(p_Kcol(ia))
        
        ! Search the DOF in our list. Ok, that's an n^2 algorithm.
        ! It could easily replaced by an n log n algorithm since the list
        ! of DOF's is sorted!
        DO j=1,SIZE(p_ImirrorDOFsClosed)
          IF (p_ImirrorDOFsClosed(j) .EQ. idof) THEN
            p_Da(ia) = dmirrorWeight * p_Da(ia)
            EXIT
          END IF
        END DO
        
      END DO
      
    END DO

  END IF
  
  END SUBROUTINE
  
! ***************************************************************************
! Block matrix filters
! ***************************************************************************

  ! *************************************************************************
  ! Implementation of discrete boundary conditions into block matrices
  ! *************************************************************************

!<subroutine>

  SUBROUTINE matfil_discreteBC (rmatrix,rdiscreteBC)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to 
  ! block matrix' filter.
  ! The matrix is modified either to rdiscreteBC (if specified) or
  ! to the default boundary conditions associated to the matrix 
  ! (if rdiscreteBC is not specified).
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the matrix.
  ! If not specified, the default boundary conditions associated to the
  ! matrix rmatrix are imposed to the matrix.
  TYPE(t_discreteBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteBC
!</input>
  
!<inputoutput>
  ! The block matrix where the boundary conditions should be imposed.
  TYPE(t_matrixBlock), INTENT(INOUT),TARGET :: rmatrix
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,jblock,i,inumEntries !,icp
    LOGICAL :: boffdiagSubmatrix
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC
    
    ! Imposing boundary conditions normally changes the whole matrix!
    ! Take the given BC structure or
    ! grab the boundary condition entry list from the matrix. This
    ! is a list of all discretised boundary conditions in the system.
    IF (PRESENT(rdiscreteBC)) THEN
      p_RdiscreteBC => rdiscreteBC%p_RdiscBCList
      inumEntries = rdiscreteBC%inumEntriesUsed
    ELSE
      IF (.NOT. ASSOCIATED(rmatrix%p_rdiscreteBC)) THEN
        ! There are no BC's available, so we cannot do anything!
        RETURN
      END IF
      p_RdiscreteBC => rmatrix%p_rdiscreteBC%p_RdiscBCList  
      inumEntries = rmatrix%p_rdiscreteBC%inumEntriesUsed
    END IF
    
    IF (.NOT. ASSOCIATED(p_RdiscreteBC)) RETURN
    
    ! Is the matrix a 'primal' matrix or is it a submatrix of another block matrix?
    boffdiagSubmatrix = rmatrix%imatrixSpec .EQ. LSYSBS_MSPEC_OFFDIAGSUBMATRIX
    
    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBC)
    DO i=1, inumEntries
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteBC(i)%itype)
      CASE (DISCBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_RdiscreteBC(i)%rdirichletBCs%icomponent
        
        ! Loop through this matrix row and implement the boundary conditions
        ! into the scalar submatrices.
        ! For now, this implements unit vectors into the diagonal matrices
        ! and zero-vectors into the offdiagonal matrices.
        ! Only exception: If the matrix is a submatrix of another matrix
        ! and not on the diagonal of its parent, we must replace the rows
        ! by zero vectors!
        DO jblock = 1,rmatrix%ndiagBlocks
          IF (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) THEN
            CALL matfil_imposeDirichletBC (&
                        rmatrix%RmatrixBlock(iblock,jblock), &
                        (iblock .NE. jblock) .OR. boffdiagSubmatrix,&
                        p_RdiscreteBC(i)%rdirichletBCs)
          END IF
        END DO
        
      CASE (DISCBC_TPPRESSUREDROP)  
        ! Nothing to do; pressure drop BC's are implemented only into the RHS.

      CASE (DISCBC_TPSLIP)  
        ! Slip boundary conditions are treated like Dirichlet for all
        ! velocity components.
        ! This is a separate filter must be called manually.
        ! Therefore, there's nothing to do here.
        
      CASE (DISCBC_TPFEASTMIRROR)  
        ! FEAST mirror boundary conditions.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_RdiscreteBC(i)%rfeastMirrorBCs%icomponent

        ! Loop through this matrix row and implement the boundary conditions
        ! into the scalar submatrices.
        ! For now, this implements unit vectors into the diagonal matrices
        ! and zero-vectors into the offdiagonal matrices.
        ! Only exception: If the matrix is a submatrix of another matrix
        ! and not on the diagonal of its parent, we must replace the rows
        ! by zero vectors!
        DO jblock = 1,rmatrix%ndiagBlocks
          IF (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) THEN
            CALL matfil_imposeFeastMirrorBC (&
                        rmatrix%RmatrixBlock(iblock,jblock), &
                        (iblock .NE. jblock) .OR. boffdiagSubmatrix,&
                        p_RdiscreteBC(i)%rfeastMirrorBCs)
          END IF
        END DO

      CASE DEFAULT
        CALL output_line (&
            'Unknown boundary condition'//sys_siL(p_RdiscreteBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'matfil_discreteBC')
        CALL sys_halt()
        
      END SELECT
    END DO
  
  END SUBROUTINE

  ! *************************************************************************

!<subroutine>

  SUBROUTINE matfil_discreteNLSlipBC (rmatrix,bforprec,rdiscreteBC)

!<description>
  ! Imposes nonlinear slip boundary conditions into a given matrix.
  ! The matrix is modified either to rdiscreteBC (if specified) or
  ! to the default boundary conditions associated to the matrix 
  ! (if rdiscreteBC is not specified).
  !
  ! Note that slip boundary conditions are implemented by so called
  ! 'nonlinear filters' (c.f. vectorfilters.f). Similarly to the vector
  ! filters, slip boundary conditions are not automatically implemented
  ! by matfil_discreteBC. To implement them, this routine must be called
  ! separately from matfil_discreteBC.
!</description>

!<input>
  ! Prepare matrix for preconditioning.
  ! =false: Slip-nodes are handled as Dirichlet. Rows in the
  !         matrix are replaced by unit vectors.
  ! =true : Standard. Prepare matrix for preconditioning.
  !         Only the off-diagonals of rows corresponding to 
  !         the slip-nodes are cleared.
  !         The matrix therefore is prepared to act as 
  !         Jacobi-preconditioner for all slip-nodes.
  LOGICAL, INTENT(IN) :: bforprec

  ! OPTIONAL: boundary conditions to impose into the matrix.
  ! If not specified, the default boundary conditions associated to the
  ! matrix rmatrix are imposed to the matrix.
  TYPE(t_discreteBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteBC
!</input>
  
!<inputoutput>
  ! The block matrix where the boundary conditions should be imposed.
  TYPE(t_matrixBlock), INTENT(INOUT),TARGET :: rmatrix
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,jblock,i,icp
    LOGICAL :: boffdiagSubmatrix
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC
    
    ! Imposing boundary conditions normally changes the whole matrix!
    ! Grab the boundary condition entry list from the matrix. This
    ! is a list of all discretised boundary conditions in the system.
    p_RdiscreteBC => rmatrix%p_rdiscreteBC%p_RdiscBCList  
    
    IF (.NOT. ASSOCIATED(p_RdiscreteBC)) RETURN
    
    ! Is the matrix a 'primal' matrix or is it a submatrix of another block matrix?
    boffdiagSubmatrix = rmatrix%imatrixSpec .EQ. LSYSBS_MSPEC_OFFDIAGSUBMATRIX

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBC)
    DO i=1, rmatrix%p_rdiscreteBC%inumEntriesUsed
    
      ! Only implement slip boundary conditions.
      IF (p_RdiscreteBC(i)%itype .EQ. DISCBC_TPSLIP) THEN
        ! Slip boundary conditions are treated like Dirichlet for all
        ! velocity components.

        ! Loop through all affected components to implement the BC's.
        DO icp = 1,p_RdiscreteBC(i)%rslipBCs%ncomponents
          iblock = p_RdiscreteBC(i)%rslipBCs%Icomponents(icp)

          DO jblock = 1,rmatrix%ndiagBlocks
            IF (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) THEN
              CALL matfil_imposeNLSlipBC (&
                          rmatrix%RmatrixBlock(iblock,jblock), &
                          (iblock .NE. jblock) .OR. boffdiagSubmatrix,bforprec,&
                          p_RdiscreteBC(i)%rslipBCs)
            END IF
          END DO
          
        END DO
      END IF
    
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into 
  ! block matrices
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_discreteFBC (rmatrix,rdiscreteFBC)

!<description>
  ! This routine realises the 'impose discrete fictitious boundary 
  ! conditions to block matrix' filter.
  ! The matrix is modified either to rdiscreteFBC (if specified) or
  ! to the default boundary conditions associated to the matrix 
  ! (if rdiscreteFBC is not specified).
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the matrix.
  ! If not specified, the default boundary conditions associated to the
  ! matrix rmatrix are imposed to the matrix.
  TYPE(t_discreteFBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteFBC
!</input>
  
!<inputoutput>
  ! The block matrix where the boundary conditions should be imposed.
  TYPE(t_matrixBlock), INTENT(INOUT),TARGET :: rmatrix
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,jblock,i,j
    LOGICAL :: boffdiagSubmatrix
    TYPE(t_discreteFBCEntry), DIMENSION(:), POINTER :: p_RdiscreteFBC
    
    ! Imposing boundary conditions normally changes the whole matrix!
    ! Grab the boundary condition entry list from the matrix. This
    ! is a list of all discretised boundary conditions in the system.
    IF (.NOT. ASSOCIATED(rmatrix%p_rdiscreteBCfict)) RETURN

    p_RdiscreteFBC => rmatrix%p_rdiscreteBCfict%p_RdiscFBCList  
    
    IF (.NOT. ASSOCIATED(p_RdiscreteFBC)) RETURN
    
    ! Is the matrix a 'primal' matrix or is it a submatrix of another block matrix?
    boffdiagSubmatrix = rmatrix%imatrixSpec .EQ. LSYSBS_MSPEC_OFFDIAGSUBMATRIX

    ! Now loop through all entries in this list:
    DO i=1,SIZE(p_RdiscreteFBC)
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteFBC(i)%itype)
      CASE (DISCFBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCFBC_TPDIRICHLET)
        ! Dirichlet boundary conditions.
        ! 
        ! Loop through all blocks where to impose the BC's:
        DO j=1,p_RdiscreteFBC(i)%rdirichletFBCs%ncomponents
        
          iblock = p_RdiscreteFBC(i)%rdirichletFBCs%Icomponents(j)
        
          ! Loop through this matrix row and implement the boundary conditions
          ! into the scalar submatrices.
          ! For now, this implements unit vectors into the diagonal matrices
          ! and zero-vectors into the offdiagonal matrices.
          ! Only exception: If the matrix is a submatrix of another matrix
          ! and not on the diagonal of its parent, we must replace the rows
          ! by zero vectors!
          DO jblock = 1,rmatrix%ndiagBlocks
            IF (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) THEN
              CALL matfil_imposeDirichletFBC (&
                          rmatrix%RmatrixBlock(iblock,jblock), &
                          (iblock .NE. jblock) .OR. boffdiagSubmatrix,&
                          p_RdiscreteFBC(i)%rdirichletFBCs)
            END IF
          END DO
          
        END DO
        
      CASE DEFAULT
        CALL output_line (&
            'Unknown boundary condition'//sys_siL(p_RdiscreteFBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'matfil_discreteFBC')
        CALL sys_halt()
        
      END SELECT
    END DO
  
  END SUBROUTINE

END MODULE
