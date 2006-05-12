!##############################################################################
!# ****************************************************************************
!# <name> VectorFilters </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a set of filters to impose discrete boundary
!# conditions into vectors. The structures that define discrete BC's are
!# defined in the module DiscreteBC.
!# Discrete BC's can generally be imposed into solution vectors as well as
!# into defect vectors. The routines
!# - vecfil_imposeDiscreteBC and
!# - vecfil_imposeDiscreteDefectBC
!# serve as a wrapper for the t_discreteBC structure, which call the right
!# 'imposing' routine, depending on which boundary conditions the
!# t_discreteBC structure describes.
!#
!# The following types of discrete boundary conditions are available:
!# - Dirichlet boundary conditions are typically represented as a
!#   list of DOF's where a value must be prescribed, plus the value that
!#   must be described in that special DOF.
!#   
!# </purpose>
!##############################################################################

MODULE vectorfilters

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE discretebc
  USE dofmapping
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeDiscreteBC (rx)

  !<description>
  
  ! This routine serves as a wrapper for implementing discrete boundary
  ! conditions into a 'solution' vector. Depending on the type of boundary 
  ! conditions, the correct 'imposing' routine will be called that imposes 
  ! the actual boundary conditions.
  
  !</description>
  
  !<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
  
  !</inputoutput>

!</subroutine>

  INTEGER :: iblock, ibc, ibctype
  
  ! Loop over the blocks
  
  DO iblock = 1,LSYSBL_MAXBLOCKS

    ! Is there a vector in block iblock of rx?
    
    IF (rx%RvectorBlock(iblock)%h_Ddata .NE. ST_NOHANDLE) THEN
    
      ! Loop over the BC's that are to be imposed
      
      IF (ASSOCIATED(rx%RvectorBlock(iblock)%p_RdiscreteBC)) THEN
      
        DO ibc = 1,SIZE(rx%RvectorBlock(iblock)%p_RdiscreteBC)
          
          ! Choose the right boundary condition implementation routine
          ! and call it for the vector.
        
          ibctype = rx%RvectorBlock(iblock)%p_RdiscreteBC(ibc)%itype

          SELECT CASE(ibctype) 
          CASE(DISCBC_TPDIRICHLET)
            CALL vecfil_imposeDirichletBC (rx%RvectorBlock(iblock), &
                rx%RvectorBlock(iblock)%p_RdiscreteBC(ibc)%rdirichletBCs)
          END SELECT 
          
        END DO
        
      END IF
    
    END IF
    
  END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeDiscreteDefectBC (rx)

  !<description>
  
  ! This routine implements discrete boundary conditions into a block
  ! defect vector. Depending on the type of boundary 
  ! conditions, the correct 'imposing' routine will be called that imposes 
  ! the actual boundary conditions into each block.
  
  !</description>
  
  !<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
  
  !</inputoutput>

!</subroutine>

  INTEGER :: iblock, ibc, ibctype
  
  ! Loop over the blocks
  
  DO iblock = 1,LSYSBL_MAXBLOCKS

    ! Is there a vector in block iblock of rx?
    
    IF (rx%RvectorBlock(iblock)%h_Ddata .NE. ST_NOHANDLE) THEN
    
      ! Loop over the BC's that are to be imposed
      
      IF (ASSOCIATED(rx%RvectorBlock(iblock)%p_RdiscreteBC)) THEN
      
        DO ibc = 1,SIZE(rx%RvectorBlock(iblock)%p_RdiscreteBC)
          
          ! Choose the right boundary condition implementation routine
          ! and call it for the vector.
        
          ibctype = rx%RvectorBlock(iblock)%p_RdiscreteBC(ibc)%itype

          SELECT CASE(ibctype) 
          CASE(DISCBC_TPDIRICHLET)
            CALL vecfil_imposeDirichletDefectBC (rx%RvectorBlock(iblock), &
                rx%RvectorBlock(iblock)%p_RdiscreteBC(ibc)%rdirichletBCs)
          END SELECT 
          
        END DO
        
      END IF
    
    END IF
    
  END DO
    
  END SUBROUTINE

! *****************************************************************************
! Discrete Dirichlet boundary conditions
! *****************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeDirichletBC (rx,rdbcStructure)
  
  !<description>
  
  ! Implements discrete Dirichlet BC's into a scalar vector.

  !<input>
  
  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteBCDirichlet), INTENT(IN), TARGET  :: rdbcStructure
  
  !</input>

  !<inputoutput>

  ! The scalar vector where the boundary conditions should be imposed.
  TYPE(t_vectorScalar), INTENT(INOUT), TARGET :: rx
  
  !</inputoutput>
  
!</subroutine>
    
  ! local variables
  ! INTEGER(PREC_DOFIDX) :: i
  REAL(DP), DIMENSION(:), POINTER    :: p_vec
  INTEGER(I32), DIMENSION(:), POINTER :: p_idx
  REAL(DP), DIMENSION(:), POINTER    :: p_val

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.
  
  CALL storage_getbase_double (rx%h_Ddata, p_vec)  
  CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
  CALL storage_getbase_double(rdbcStructure%h_IdirichletValues,p_val)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF ((.NOT.ASSOCIATED(p_idx)).OR.(.NOT.ASSOCIATED(p_val))) THEN
    PRINT *,'Error: DBC not configured'
    STOP
  END IF
  
  ! Only handle nDOF DOF's, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  
  p_vec(p_idx(1:rdbcStructure%nDOF)) = p_val(1:rdbcStructure%nDOF)

!  DO i=1,rdbcStructure%nDOF
!    p_vec(p_idx(i)) = p_val(i)
!  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeDirichletDefectBC (rx,rdbcStructure)
  
  !<description>
  
  ! Implements discrete Dirichlet BC's into a defect block vector.

  !<input>
  
  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteBCDirichlet), INTENT(IN),TARGET  :: rdbcStructure
  
  !</input>

  !<inputoutput>

  ! The scalar vector where the boundary conditions should be imposed.
  TYPE(t_vectorScalar), INTENT(INOUT),TARGET :: rx
  
  !</inputoutput>
  
!</subroutine>

  ! local variables
  ! INTEGER(PREC_DOFIDX) :: i
  REAL(DP), DIMENSION(:), POINTER :: p_vec
  INTEGER(I32), DIMENSION(:), POINTER :: p_idx
  
  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.
  
  CALL storage_getbase_double (rx%h_Ddata, p_vec)  
  CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
  
  ! Impose the BC-DOF's directly - more precisely, into the
  ! components of the subvector that is indexed by icomponent.

  IF (.NOT.ASSOCIATED(p_idx)) THEN
    PRINT *,'Error: DBC not configured'
    STOP
  END IF
  
  ! Only handle nDOF DOF's, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  
  p_vec(p_idx(1:rdbcStructure%nDOF)) = 0.0_DP

!  DO i=1,rdbcStructure%nDOF
!    p_vec(p_idx(i)) =0 
!  END DO
  
  END SUBROUTINE
  
END MODULE
