!##############################################################################
!# ****************************************************************************
!# <name> vectorfilters </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of filters that can be applied to a scalar
!# or block vector. These can be used e.g. during a solution process to impose
!# different quantities into solution and/or defect vectors.
!#
!# The discrete boundary conditions realised in the module 'bcassembly' are
!# one type of filter. While being created in the module 'bcassembly', 
!# this module now provides the functionality to impose discrete boundary 
!# conditions into a vector.
!# Also other filters can be found here, e.g. normalisation ov vectors, etc.
!#
!# Filters can even be collected to a complete 'filter chain' and applied
!# 'en block' onto a vector. For this purpose, there exists a higher-level
!# module 'filtersupport', which realises such a filter chain.
!#
!# The following routines can be found here:
!#
!# 3.) vecfil_normaliseToL20Sca
!#     -> Normalise a scalar vector to be in the space $L^2_0$.
!#
!# 4.) vecfil_discreteBCsol
!#     -> Apply the 'discrete boundary conditions for solution vectors' filter
!#        onto a given (block) solution vector. 
!#
!# 5.) vecfil_discreteBCrhs
!#     -> Apply the 'discrete boundary conditions for RHS vectors' filter
!#        onto a given (block) vector. 
!#
!# 6.) vecfil_discreteBCdef
!#     -> Apply the 'discrete boundary conditions for defect vectors' filter 
!#        onto a given (block) vector. 
!#
!# 7.) vecfil_subvectorToL20
!#     -> Normalise a subvector of a block vector to be in the space $L^2_0$.
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

! *****************************************************************************
! Scalar vector filters
! *****************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeDirichletBC (rx,rdbcStructure)
  
!<description>
  ! Implements discrete Dirichlet BC's into a scalar vector.
!</description>

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
    INTEGER(PREC_DOFIDX) :: i
    REAL(DP), DIMENSION(:), POINTER    :: p_vec
    INTEGER(I32), DIMENSION(:), POINTER :: p_idx
    REAL(DP), DIMENSION(:), POINTER    :: p_val
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm

    ! If nDOF=0, there are no DOF's the current boundary condition segment,
    ! so we don't have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    IF (rdbcStructure%nDOF .EQ. 0) RETURN

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.
    
    CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    CALL storage_getbase_double(rdbcStructure%h_DdirichletValues,p_val)

    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvector that is indexed by icomponent.
    
    IF ((.NOT.ASSOCIATED(p_idx)).OR.(.NOT.ASSOCIATED(p_val))) THEN
      PRINT *,'Error: DBC not configured'
      STOP
    END IF
    
    CALL lsyssc_getbase_double (rx, p_vec)  
    
    IF (.NOT.ASSOCIATED(p_vec)) THEN
      PRINT *,'Error: No vector'
      STOP
    END IF

    ! Only handle nDOF DOF's, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted?
    IF (rx%isortStrategy .LE. 0) THEN
      ! No. Implement directly.
      DO i=1,rdbcStructure%nDOF
        p_vec(p_idx(i)) = p_val(i)
      END DO
    ELSE
      ! Ups, vector sorted. At first get the permutation how its sorted
      ! - or more precisely, the back-permutation, as we need this one for 
      ! the loop below.
      CALL storage_getbase_int (rx%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)
      
      ! And 'filter' each DOF during the boundary value implementation!
      DO i=1,rdbcStructure%nDOF
        p_vec(p_Iperm(p_idx(i))) = p_val(i)
      END DO
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeDirichletDefectBC (rx,rdbcStructure)
  
!<description>
  ! Implements discrete Dirichlet BC's into a scalar defect vector.
!</description>

!<input>
  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteBCDirichlet), INTENT(IN),TARGET  :: rdbcStructure
!</input>

!<inputoutput>
  ! The scalar vector where the boundary conditions should be imposed.
  TYPE(t_vectorScalar), INTENT(INOUT), TARGET :: rx
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER(PREC_DOFIDX) :: i
    REAL(DP), DIMENSION(:), POINTER :: p_vec
    INTEGER(I32), DIMENSION(:), POINTER :: p_idx
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
    
    ! If nDOF=0, there are no DOF's the current boundary condition segment,
    ! so we don't have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    IF (rdbcStructure%nDOF .EQ. 0) RETURN

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.
    
    CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    
    IF (.NOT.ASSOCIATED(p_idx)) THEN
      PRINT *,'Error: DBC not configured'
      STOP
    END IF
    
    CALL lsyssc_getbase_double (rx, p_vec)  
    
    IF (.NOT.ASSOCIATED(p_vec)) THEN
      PRINT *,'Error: No vector'
      STOP
    END IF

    ! Impose the BC-DOF's directly - more precisely, into the
    ! components of the subvector that is indexed by icomponent.
    !
    ! Only handle nDOF DOF's, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted?
    IF (rx%isortStrategy .LE. 0) THEN
      ! No. Implement directly.
      DO i=1,rdbcStructure%nDOF
        p_vec(p_idx(i)) = 0.0_DP
      END DO
    ELSE
      ! Ups, vector sorted. At first get the permutation how its sorted -
      ! or more precisely, the back-permutation, as we need this one for 
      ! the loop below.
      CALL storage_getbase_int (rx%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)
      
      ! And 'filter' each DOF during the boundary value implementation!
      DO i=1,rdbcStructure%nDOF
        p_vec(p_Iperm(p_idx(i))) = 0.0_DP
      END DO
    END IF
  
  END SUBROUTINE

  ! ***************************************************************************
  ! Other scalar filters
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE vecfil_normaliseToL20Sca (rx)

!<description>
  ! This routine normalises a scalar vector to bring it into the space $L^2_0$.
!</description>
  
!<inputoutput>
  ! The vector which is to be normalised.
  TYPE(t_vectorScalar), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    REAL(DP), DIMENSION(:), POINTER        :: p_DelementArea,p_Ddata
    REAL(DP) :: dpintegral,c
    INTEGER(PREC_ELEMENTIDX) :: iel,nel
    
    ! Get the discretisation...
    IF (.NOT. ASSOCIATED(rx%p_rspatialDiscretisation)) RETURN

    p_rdiscretisation => rx%p_rspatialDiscretisation
    
    ! ... and check it. If we have a uniform discretisation with P_0/Q_0,
    ! we have the easy case, that the integral of the function rx is
    ! representing is area*rx. Otherwise we have to calculate the integral
    ! which is somehow more costly...
    
    IF ((p_rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) .AND. &
        ((p_rdiscretisation%RelementDistribution(1)%itrialElement .EQ. EL_P0) .OR. &
        (p_rdiscretisation%RelementDistribution(1)%itrialElement .EQ. EL_Q0))) THEN
      
      ! Ok, easy case. Get from the triangulation the AREA-array for calculating
      ! a simple integral of rx:
      
      CALL storage_getbase_double (p_rdiscretisation%p_rtriangulation%h_DelementArea, &
                                  p_DelementArea)
                                   
      ! Get the vector data of rx
      CALL lsyssc_getbase_double (rx,p_Ddata)
      
      nel = SIZE(p_DelementArea)-1
                                   
      ! Build the integral
      !   int_Omega p dx
      ! This is approximated by
      !   dpintegral = SUM_Elements P(Element)*Volume(Element)
      
      dpintegral=0D0
      DO iel=1,nel 
        dpintegral = dpintegral + p_Ddata(iel)*p_DelementArea(iel)
      END DO
      
      ! Divide dpintegral by the volume of the domain; this gives the integral
      ! mean value of the pressure:
        
      C = dpintegral / p_DelementArea(nel+1)
      
      ! Subtract the integral mean value C of the pressure from all
      ! pressure components. Afterwards, P has integral mean value = 0.

      DO iel=1,nel
        p_Ddata(iel) = p_Ddata(iel) - C
      END DO
         
    ELSE
      PRINT *,'Normalisation of non-P0 vectors not implemented!'
      STOP
    END IF

  END SUBROUTINE

! ***************************************************************************
! Block vector filters
! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposePressureDropBC (rx,rpdbcStructure)
  
!<description>
  ! Implements discrete pressure drop BC's into a block vector.
!</description>

!<input>
  ! The t_discreteBCpressureDrop that describes the discrete pressure 
  ! drop BC's
  TYPE(t_discreteBCpressureDrop), INTENT(IN), TARGET  :: rpdbcStructure
!</input>

!<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorblock), INTENT(INOUT), TARGET :: rx
  
!</inputoutput>
  
!</subroutine>
    
    ! local variables
    INTEGER :: j,icp
    INTEGER(PREC_DOFIDX) :: i
    REAL(DP), DIMENSION(:), POINTER    :: p_vec
    INTEGER(I32), DIMENSION(:), POINTER :: p_idx
    REAL(DP), DIMENSION(:,:), POINTER    :: p_val
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm

    ! If nDOF=0, there are no DOF's the current boundary condition segment,
    ! so we don't have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    IF (rpdbcStructure%nDOF .EQ. 0) RETURN

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.
    
    CALL storage_getbase_int(rpdbcStructure%h_IpressureDropDOFs,p_idx)
    CALL storage_getbase_double2d(rpdbcStructure%h_Dmodifier,p_val)

    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvectors that is indexed by Icomponent(1..NDIM2D).
    
    IF ((.NOT.ASSOCIATED(p_idx)).OR.(.NOT.ASSOCIATED(p_val))) THEN
      PRINT *,'Error: pressure drop BC not configured'
      STOP
    END IF
    
    ! Currently, only 2D is supported. 
    ! First handle the X-velocity (j=1), then the Y-velocity (j=2)
    DO j=1,NDIM2D
    
      ! Get the subvector
      icp = rpdbcStructure%Icomponents(j)
      CALL lsyssc_getbase_double (rx%RvectorBlock(icp), p_vec)
    
      IF (.NOT.ASSOCIATED(p_vec)) THEN
        PRINT *,'Error: No vector'
        STOP
      END IF

      ! Only handle nDOF DOF's, not the complete array!
      ! Probably, the array is longer (e.g. has the length of the vector), but
      ! contains only some entries...

      ! Is the vector sorted?
      IF (rx%RvectorBlock(j)%isortStrategy .LE. 0) THEN
        ! No. Implement directly.
        DO i=1,rpdbcStructure%nDOF
          p_vec(p_idx(i)) = p_vec(p_idx(i)) + p_val(j,i)
        END DO
      ELSE
        ! Oops, vector sorted. At first get the permutation how its sorted
        ! - or more precisely, the back-permutation, as we need this one for 
        ! the loop below.
        CALL storage_getbase_int (rx%RvectorBlock(j)%h_IsortPermutation,p_Iperm)
        p_Iperm => p_Iperm(rx%RvectorBlock(j)%NEQ+1:)

        ! And 'filter' each DOF during the boundary value implementation!
        DO i=1,rpdbcStructure%nDOF
          p_vec(p_Iperm(p_idx(i))) = p_vec(p_Iperm(p_idx(i))) + p_val(j,i)
        END DO
      END IF
      
    END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block solution vectors
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteBCsol (rx,dtimeWeight,rdiscreteBC)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to solution'
  ! filter. This filter imposes the discrete boundary conditions rdiscreteBC
  ! (if specified) or (if rdiscreteBC is not specified) the boundary conditions
  ! which are  associated to the vector rx (with rx%p_discreteBC) to this 
  ! (block) vector.
!</description>
  
!<input>
  ! OPTIONAL: Time-step weight. This weight is multiplied to time-dependent
  ! boundary conditions before these are added to the vector rx.
  ! The parameter can be omitted in stationary simulations or when filtering
  ! vectors during the solution process of a linear system.
  REAL(DP), INTENT(IN), OPTIONAL :: dtimeWeight

  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  TYPE(t_discreteBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,i
    REAL(DP) :: dtweight
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC

    IF (.NOT. PRESENT(rdiscreteBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      p_RdiscreteBC => rx%p_rdiscreteBC%p_RdiscBCList  
    ELSE
      p_RdiscreteBC => rdiscreteBC%p_RdiscBCList
    END IF
    
    IF (.NOT. ASSOCIATED(p_RdiscreteBC)) RETURN
    
    ! If the time-weight is not specified, 1.0 is assumed.
    IF (PRESENT(dtimeWeight)) THEN
      dtweight = dtimeWeight
    ELSE
      dtweight = 1.0_DP
    END IF
    ! Note: Time-step weight not used by any filter up to now!
    ! Perhaps in a later implementation it's needed anywhere...
    
    ! Now loop through all entries in this list:
    DO i=1,SIZE(p_RdiscreteBC)
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteBC(i)%itype)
      CASE (DISCBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_RdiscreteBC(i)%rdirichletBCs%icomponent
        
        ! Implement the Dirichlet boundary conditions into that component
        ! of the vector.
        CALL vecfil_imposeDirichletBC (rx%RvectorBlock(iblock),&
                                      p_RdiscreteBC(i)%rdirichletBCs)
        
      CASE (DISCBC_TPPRESSUREDROP)  
        ! Nothing to do; pressure drop BC's are implemented only into the RHS.
        
      CASE DEFAULT
        PRINT *,'vecfil_discreteBCsol: unknown boundary condition: ',&
                p_RdiscreteBC(i)%itype
        STOP
        
      END SELECT
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block RHS vectors
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteBCrhs (rx,dtimeWeight,rdiscreteBC)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to RHS'
  ! filter. This filter imposes the discrete boundary conditions rdiscreteBC
  ! (if specified) or (if rdiscreteBC is not specified) the boundary conditions
  ! which are  associated to the vector rx (with rx%p_discreteBC) to this 
  ! (block) vector.
!</description>
  
!<input>
  ! OPTIONAL: Time-step weight. This weight is multiplied to time-dependent
  ! boundary conditions before these are added to the vector rx.
  ! The parameter can be omitted in stationary simulations.
  REAL(DP), INTENT(IN), OPTIONAL :: dtimeWeight

  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  TYPE(t_discreteBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,i
    REAL(DP) :: dtweight
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC

    IF (.NOT. PRESENT(rdiscreteBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      p_RdiscreteBC => rx%p_rdiscreteBC%p_RdiscBCList  
    ELSE
      p_RdiscreteBC => rdiscreteBC%p_RdiscBCList
    END IF
    
    IF (.NOT. ASSOCIATED(p_RdiscreteBC)) RETURN
    
    ! If the time-weight is not specified, 1.0 is assumed.
    IF (PRESENT(dtimeWeight)) THEN
      dtweight = dtimeWeight
    ELSE
      dtweight = 1.0_DP
    END IF
    
    ! Now loop through all entries in this list:
    DO i=1,SIZE(p_RdiscreteBC)
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteBC(i)%itype)
      CASE (DISCBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_RdiscreteBC(i)%rdirichletBCs%icomponent
        
        ! Implement the Dirichlet boundary conditions into that component
        ! of the vector.
        CALL vecfil_imposeDirichletBC (rx%RvectorBlock(iblock),&
                                      p_RdiscreteBC(i)%rdirichletBCs)
      
      CASE (DISCBC_TPPRESSUREDROP)  
        CALL vecfil_imposePressureDropBC (rx,p_RdiscreteBC(i)%rpressureDropBCs)
        
      CASE DEFAULT
        PRINT *,'vecfil_discreteBCrhs: unknown boundary condition: ',&
                p_RdiscreteBC(i)%itype
        STOP
        
      END SELECT
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteBCdef (rx,dtimeWeight,rdiscreteBC)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to defect'
  ! filter. This filter imposes the discrete boundary conditions rdiscreteBC
  ! (if specified) or (if rdiscreteBC is not specified) the boundary conditions
  ! which are  associated to the defect vector rx (with rx%p_discreteBC) to 
  ! this (block) vector.
!</description>
  
!<input>
  ! OPTIONAL: Time-step weight. This weight is multiplied to time-dependent
  ! boundary conditions before these are added to the vector rx.
  ! The parameter can be omitted in stationary simulations or when filtering
  ! vectors during the solution process of a linear system.
  REAL(DP), INTENT(IN), OPTIONAL :: dtimeWeight

  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  TYPE(t_discreteBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,i
    REAL(DP) :: dtweight
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC

    IF (.NOT. PRESENT(rdiscreteBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      p_RdiscreteBC => rx%p_rdiscreteBC%p_RdiscBCList  
    ELSE
      p_RdiscreteBC => rdiscreteBC%p_RdiscBCList
    END IF
    
    IF (.NOT. ASSOCIATED(p_RdiscreteBC)) RETURN
    
    ! If the time-weight is not specified, 1.0 is assumed.
    IF (PRESENT(dtimeWeight)) THEN
      dtweight = dtimeWeight
    ELSE
      dtweight = 1.0_DP
    END IF
    ! Note: Time-step weight not used by any filter up to now!
    ! Perhaps in a later implementation it's needed anywhere...
    
    ! Now loop through all entries in this list:
    DO i=1,SIZE(p_RdiscreteBC)
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteBC(i)%itype)
      CASE (DISCBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_RdiscreteBC(i)%rdirichletBCs%icomponent
        
        ! Implement the Dirichlet boundary conditions into that component
        ! of the vector.
        CALL vecfil_imposeDirichletDefectBC (rx%RvectorBlock(iblock),&
                                            p_RdiscreteBC(i)%rdirichletBCs)
        
      CASE (DISCBC_TPPRESSUREDROP)  
        ! Nothing to do; pressure drop BC's are implemented only into the RHS.
        
      CASE DEFAULT
        PRINT *,'vecfil_discreteBCdef: unknown boundary condition: ',&
                p_RdiscreteBC(i)%itype
        STOP
        
      END SELECT
    END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************
  ! Other block filters
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_subvectorToL20 (rx,isubvector)

!<description>
  ! This routine realises the 'subvector to $L^2_0$ filter.
  ! The subvector isubvector of the block vector rx is normalised
  ! with vecfil_normaliseScalarToL20 to bring it to the space $L^2_0$.
!</description>
  
!<inputoutput>

  ! The block vector which is partially to be normalised.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
  
  ! The number of the subvector in rx which is to be normalised.
  INTEGER, INTENT(IN) :: isubvector
  
!</inputoutput>

!</subroutine>

    IF ((isubvector .LE. 0) .OR. (isubvector .GT. SIZE(rx%RvectorBlock)) .OR. &
        (rx%RvectorBlock(isubvector)%h_Ddata .EQ. ST_NOHANDLE) .OR. &
        (rx%RvectorBlock(isubvector)%NEQ .LE. 0)) RETURN
        
    ! Normalise the subvector isubvector
    CALL vecfil_normaliseToL20Sca (rx%RvectorBlock(isubvector))

  END SUBROUTINE

END MODULE
