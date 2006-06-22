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
!# 1.) matfil_discreteBC
!#     -> Apply the 'discrete boundary conditions for solution vectors' filter
!#        onto a given (block) matrix.
!#        This e.g. replaces lines of a matrix corresponding to Dirichlet DOF's
!#        by unit vectors for all scalar submatrices where configured.
!#
!# </purpose>
!##############################################################################

MODULE matrixfilters

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE discretebc
  USE dofmapping
  USE matrixmodification
  
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

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.
  
  CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF (.NOT.ASSOCIATED(p_idx)) THEN
    PRINT *,'Error: DBC not configured'
    STOP
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
! Block matrix filters
! ***************************************************************************

  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block matrices,
  ! ***************************************************************************

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

  INTEGER :: iblock,jblock,i
  TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC
  
  ! Imposing boundary conditions normally changes the whole matrix!
  ! Grab the boundary condition entry list from the matrix. This
  ! is a list of all discretised boundary conditions in the system.
  p_RdiscreteBC => rmatrix%p_rdiscreteBC%p_RdiscBCList  
  
  IF (.NOT. ASSOCIATED(p_RdiscreteBC)) RETURN
  
  ! Now loop through all entries in this list:
  DO i=1,SIZE(p_RdiscreteBC)
  
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
      DO jblock = 1,rmatrix%ndiagBlocks
        CALL matfil_imposeDirichletBC (&
                    rmatrix%RmatrixBlock(iblock,jblock), &
                    iblock .NE. jblock,p_RdiscreteBC(i)%rdirichletBCs)
      END DO
      
    CASE DEFAULT
      PRINT *,'matfil_discreteBC: unknown boundary condition: ',&
              p_RdiscreteBC(i)%itype
      STOP
      
    END SELECT
  END DO
  
  END SUBROUTINE

END MODULE
