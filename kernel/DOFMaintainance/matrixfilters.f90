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
!# of the matrix are usually replaced by unit vectors. The corresponding
!# matrix modification routine is a kind of 'filter' for a matrix and
!# can be found in this module.
!#
!# The following routines can be found here:
!#
!# 1.) matfil_discreteBCSca
!#     -> Apply the 'discrete boundary conditions for solution vectors' filter
!#        onto a given (scalar) matrix.
!#        This e.g. replaces lines of a matrix corresponding to Dirichlet DOF's
!#        by unit vectors. 
!#
!# 2.) matfil_discreteBC
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
  ! Implementation of discrete boundary conditions into scalar matrix,
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_discreteBCSca (rmatrix,RdiscreteBC)

!<description>
  
  ! This routine serves as a wrapper for implementing discrete boundary
  ! conditions into a (scalar) matrix. Depending on the type of 
  ! boundary  conditions, the correct 'imposing' routine will be called that 
  ! imposes the actual boundary conditions.
  !
  ! RdiscreteBC is an optional argument describing the discrete boundary
  ! conditions. If not given, the boudnary conditions that are associated
  ! to the vector rx are imposed into rx.
  
!</description>
  
!<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_matrixScalar), INTENT(INOUT),TARGET :: rmatrix
  
  ! OPTIONAL: The boundary conditions that are to be imposed into the vector.
  ! If not given, the discrete boundary conditions associated to the vector
  ! rx (in rx%p_RdiscreteBC) are imposed to rx.
  TYPE(t_discreteBC), DIMENSION(:), OPTIONAL, INTENT(IN), TARGET :: RdiscreteBC
  
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: ibc, ibctype
  TYPE(t_discreteBC), DIMENSION(:), POINTER :: p_RdiscreteBC
  
  ! Which BC to impose?
  IF (PRESENT(RdiscreteBC)) THEN
    p_RdiscreteBC => RdiscreteBC
  ELSE
    p_RdiscreteBC => rmatrix%p_RdiscreteBC
  END IF
  
  ! Maybe that there are no BC to be imposed - e.g. in pure Neumann problems!
  IF (.NOT. ASSOCIATED(p_RdiscreteBC)) RETURN
  
  ! Loop over the BC's that are to be imposed
  
  DO ibc = 1,SIZE(RdiscreteBC)
    
    ! Choose the right boundary condition implementation routine
    ! and call it for the matrix.
  
    ibctype = p_RdiscreteBC(ibc)%itype

    SELECT CASE(ibctype) 
    CASE(DISCBC_TPDIRICHLET)
      CALL matfil_imposeDirichletBC (rmatrix, p_RdiscreteBC(ibc)%rdirichletBCs)
    END SELECT 
    
  END DO
      
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_imposeDirichletBC (rmatrix,rdbcStructure)
  
!<description>
  ! Implements discrete Dirichlet BC's into a scalar matrix.
  ! This is normally done by replacing some lines of the matrix rmatrix
  ! (those belonging to Dirichlet nodes) by unit vectors.
!</description>

!<input>
  
  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteBCDirichlet), INTENT(IN), TARGET  :: rdbcStructure
  
!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrix
  
!</inputoutput>
  
!</subroutine>
    
  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_idx

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
  ! Use mmod_replaceLinesByUnit to replace the corresponding rows in the
  ! matrix by unit vectors. For more complicated FE spaces, this might have
  ! to be modified in the future...
  
  CALL mmod_replaceLinesByUnit (rmatrix,p_idx(1:rdbcStructure%nDOF))
  
  END SUBROUTINE
  
! ***************************************************************************
! Block matrix filters
! ***************************************************************************

  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block matrices,
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matfil_discreteBC (rmatrix)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to solution'
  ! filter. This filter imposes the discrete boundary conditions which are
  ! associated to the vector rx (with rx%p_discreteBC) to this (block) vector.
!</description>
  
!<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_matrixBlock), INTENT(INOUT),TARGET :: rmatrix
  
!</inputoutput>

!</subroutine>

  INTEGER :: iblock,jblock
  
  ! Loop over the blocks - at least over all allocated ones...
  
  DO jblock = 1,rmatrix%ndiagBlocks
    DO iblock = 1,rmatrix%ndiagBlocks

      ! Impose the discrete BC into the scalar subvector.
      IF (rmatrix%RmatrixBlock(iblock,jblock)%NA .NE. 0) THEN
        CALL matfil_discreteBCSca (rmatrix%RmatrixBlock(iblock,jblock))
      END IF
      
    END DO
  END DO
    
  END SUBROUTINE

END MODULE
