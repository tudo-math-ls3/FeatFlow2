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
!# 'Linear' and 'nonlinear' vector filters
!# ---------------------------------------
!# There are two elemental types of vector filters realised here:
!# 'Linear' filters and 'nonlinear' filters:
!#
!# a) 'Linear' filters
!#
!#   These filters realise simple vector filters like 
!#   - implementation of Dirichlet boundary conditions
!#   - filter a vector to be in $L^2_0$
!#   or similar.
!#   The name comes from the ability to be able to be used as a filter
!#   when solving a *linear system*. Some iterative solvers like BiCGStab
!#   or Multigrid allow a vector to be filtered during the iteration. 
!#   All filters marked as 'linear filters' can be collected to a
!#   'filter chain' that is applied to a vector in such an iteration.
!#   For this purpose, there exists a higher-level module 'filtersupport',
!#   which realises such a filter chain.
!#   Of course, they can be applied to a vector anytime manually, e.g. for
!#   implementing Dirichlet boundary conditions 'by hand'.
!#
!# b) 'Nonlinear' filters
!#
!#   All filters that can not be collected in a filter chain to be
!#   applied during the solution of a linear system are called
!#   'nonlinear filters'. These filters are usually called inside of
!#   a nonlinear iteration to perform a special filtering to a vector
!#   which is probably not possible to formulate in the calling convention
!#   of a linear filter - e.g. if additional auxiliary vectors are used
!#   or if a filter consists of a predictor-corrector step or similar.
!#   Example:
!#   - implementation of Slip boundary conditions.
!#   - implementation of Pressure Drop boundary conditions.
!#
!# Note that filters are allowed consist a linear and a nonlinear part!
!# In such a case, the 'nonlinear' filter part is usually called during 
!# the nonlinear iteration, while the 'linear' part is used during the
!# solution of a linear system. An example for this may be the
!# predictor-corrector implementation of Slip boundary conditions (not
!# realised here), where the nonlinear iteration treats the BC while
!# in the solution process of the linear system, all respective nodes
!# are handled as Dirichlet.
!# 
!# The following routines can be found here:
!#
!#  1.) vecfil_normaliseToL20Sca
!#      -> Linear filter
!#      -> Normalise a scalar vector to be in the space $L^2_0$.
!#
!#  2.) vecfil_discreteBCsol
!#      -> Linear filter
!#      -> Apply the 'discrete boundary conditions for solution vectors' filter
!#         onto a given (block) solution vector.
!#
!#  3.) vecfil_discreteBCrhs
!#      -> Linear filter
!#      -> Apply the 'discrete boundary conditions for RHS vectors' filter
!#         onto a given (block) vector. 
!#
!#  4.) vecfil_discreteBCdef
!#      -> Linear filter
!#      -> Apply the 'discrete boundary conditions for defect vectors' filter 
!#         onto a given (block) vector. 
!#
!#  5.) vecfil_discreteFBCsol
!#      -> Linear filter
!#      -> Apply the 'discrete fictitious boundary conditions for solution vectors'
!#         filter onto a given (block) solution vector. 
!#
!#  6.) vecfil_discreteFBCrhs
!#      -> Linear filter
!#      -> Apply the 'discrete fictitious boundary conditions for RHS vectors' 
!#         filter onto a given (block) vector. 
!#
!#  7.) vecfil_discretFBCdef
!#      -> Linear filter
!#      -> Apply the 'discrete fictitious boundary conditions for defect vectors' 
!#         filter onto a given (block) vector. 
!#
!#  8.) vecfil_subvectorToL20
!#      -> Linear filter
!#      -> Normalise a subvector of a block vector to be in the space $L^2_0$.
!#
!#  9.) vecfil_discreteNLSlipBCdef
!#      -> Nonlinear filter
!#      -> Implements discrete nonlinear slip BC's into a scalar defect vector.
!#
!# 10.) vecfil_discreteNLPDropBCrhs
!#      -> Nonlinear filter
!#      -> Implements discrete pressure drop BC's into a RHS vector.
!#
!# Auxiliary routines, usually not called by the main program:
!#
!#  1.) vecfil_imposeDirichletBC
!#      -> Implements discrete Dirichlet BC's into a scalar vector.
!#
!#  2.) vecfil_imposeDirichletDefectBC
!#      -> Implements discrete Dirichlet BC's into a scalar defect vector.
!#
!#  3.)  vecfil_imposeDirichletFBC (rx,icomponent,rdbcStructure)
!#      -> Implements discrete Dirichlet fictitious boundary conditions into a 
!#      -> scalar vector.
!#
!#  4.) vecfil_normaliseToL20Sca (rx)
!#      -> Normalises a scalar vector to bring it into the space $L^2_0$.
!#
!#  5.) vecfil_imposePressureDropBC (rx,dtimeweight,rpdbcStructure)
!#      -> Implements discrete pressure drop BC's into a block vector.
!#
!#  6.) vecfil_imposeNLSlipDefectBC 
!#      -> Implements discrete nonlinear slip BC's into a scalar defect vector
!#         as configured in the slip BC structure.
!#
!# </purpose>
!##############################################################################

MODULE vectorfilters

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE discretebc
  USE discretefbc
  USE dofmapping
  USE genoutput
  
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
    
    IF (rdbcStructure%h_DdirichletValues .EQ. ST_NOHANDLE) THEN
      CALL output_line('Dirichlet BC''s not correctly configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      CALL output_line('Are the BC''s only configured for defect values?!?',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      CALL sys_halt()
    END IF
    
    IF (rdbcStructure%h_IdirichletDOFs .EQ. ST_NOHANDLE) THEN
      CALL output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      CALL sys_halt()
    END IF
    
    CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    CALL storage_getbase_double(rdbcStructure%h_DdirichletValues,p_val)

    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvector that is indexed by icomponent.
    
    IF ((.NOT.ASSOCIATED(p_idx)).OR.(.NOT.ASSOCIATED(p_val))) THEN
      CALL output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      CALL sys_halt()
    END IF
    
    CALL lsyssc_getbase_double (rx, p_vec)  
    
    IF (.NOT.ASSOCIATED(p_vec)) THEN
      CALL output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      CALL sys_halt()
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
    
    IF (rdbcStructure%h_IdirichletDOFs .EQ. ST_NOHANDLE) THEN
      CALL output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectBC')
      CALL sys_halt()
    END IF
    
    CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    
    IF (.NOT.ASSOCIATED(p_idx)) THEN
      CALL output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectBC')
      CALL sys_halt()
    END IF
    
    CALL lsyssc_getbase_double (rx, p_vec)  
    
    IF (.NOT.ASSOCIATED(p_vec)) THEN
      CALL output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectBC')
      CALL sys_halt()
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

!<subroutine>

  SUBROUTINE vecfil_imposeDirichletFBC (rx,icomponent,rdbcStructure)
  
!<description>
  ! Implements discrete Dirichlet fictitious boundary conditions into a 
  ! scalar vector.
!</description>

!<input>
  ! The t_discreteFBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteFBCDirichlet), INTENT(IN), TARGET  :: rdbcStructure
  
  ! Index of the solution component in rdbcStructure\%Icomponent
  INTEGER, INTENT(IN) :: icomponent
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
    REAL(DP), DIMENSION(:,:), POINTER    :: p_val
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm

    ! If nDOF=0, there are no DOF's the current boundary condition segment,
    ! so we don't have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    IF (rdbcStructure%nDOF .EQ. 0) RETURN

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.
    
    IF ((rdbcStructure%h_IdirichletDOFs .EQ. ST_NOHANDLE) .OR. &
        (rdbcStructure%h_DdirichletValues .EQ. ST_NOHANDLE)) THEN
      CALL output_line('Dirichlet BC''s not correctly configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      CALL output_line('Are the BC''s only configured for defect values?!?',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      CALL sys_halt()
    END IF
    
    CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    CALL storage_getbase_double2d(rdbcStructure%h_DdirichletValues,p_val)

    IF ((.NOT.ASSOCIATED(p_idx)).OR.(.NOT.ASSOCIATED(p_val))) THEN
      CALL output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      CALL sys_halt()
    END IF
    
    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvector that is indexed by icomponent.
    
    CALL lsyssc_getbase_double (rx, p_vec)  
    
    IF (.NOT.ASSOCIATED(p_vec)) THEN
      CALL output_line('Error: No vector',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      CALL sys_halt()
    END IF

    ! Only handle nDOF DOF's, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted?
    IF (rx%isortStrategy .LE. 0) THEN
      ! No. Implement directly.
      DO i=1,rdbcStructure%nDOF
        p_vec(p_idx(i)) = p_val(icomponent,i)
      END DO
    ELSE
      ! Ups, vector sorted. At first get the permutation how its sorted
      ! - or more precisely, the back-permutation, as we need this one for 
      ! the loop below.
      CALL storage_getbase_int (rx%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)
      
      ! And 'filter' each DOF during the boundary value implementation!
      DO i=1,rdbcStructure%nDOF
        p_vec(p_Iperm(p_idx(i))) = p_val(icomponent,i)
      END DO
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeDirichletDefectFBC (rx,rdbcStructure)
  
!<description>
  ! Implements discrete Dirichlet fictitious boundary conditions into 
  ! a scalar defect vector.
!</description>

!<input>
  ! The t_discreteFBCDirichlet that describes the discrete Dirichlet BC's
  TYPE(t_discreteFBCDirichlet), INTENT(IN),TARGET  :: rdbcStructure
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
    
    IF (rdbcStructure%h_IdirichletDOFs .EQ. ST_NOHANDLE) THEN
      CALL output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectFBC')
      CALL sys_halt()
    END IF

    CALL storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    
    IF (.NOT.ASSOCIATED(p_idx)) THEN
      CALL output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectFBC')
      CALL sys_halt()
    END IF
    
    CALL lsyssc_getbase_double (rx, p_vec)  
    
    IF (.NOT.ASSOCIATED(p_vec)) THEN
      CALL output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectFBC')
      CALL sys_halt()
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
    INTEGER(PREC_ELEMENTIDX) :: iel,nel,itrial,itrialspace
    
    ! Get the discretisation...
    IF (.NOT. ASSOCIATED(rx%p_rspatialDiscretisation)) RETURN

    p_rdiscretisation => rx%p_rspatialDiscretisation
    
    ! ... and check it. If we have a uniform discretisation with P_0/Q_0,
    ! we have the easy case, that the integral of the function rx is
    ! representing is area*rx. Otherwise we have to calculate the integral
    ! which is somehow more costly...
    

    IF (p_rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
    
      itrialspace = elem_getPrimaryElement(&
          p_rdiscretisation%RelementDistribution(1)%itrialElement)
  
      SELECT CASE (itrialspace)
      CASE (EL_P0, EL_Q0) 
    
        ! Ok, easy case. Get from the triangulation the AREA-array for calculating
        ! a simple integral of rx:
        
        CALL storage_getbase_double (p_rdiscretisation%p_rtriangulation%h_DelementVolume, &
                                     p_DelementArea)
                                      
        ! Get the vector data of rx
        CALL lsyssc_getbase_double (rx,p_Ddata)
        
        nel = SIZE(p_DelementArea)-1
                                      
        ! Build the integral
        !   int_Omega p dx
        ! This is approximated by
        !   dpintegral = SUM_Elements P(Element)*Volume(Element)
        
        dpintegral=0.0_DP
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
        
      CASE (EL_QP1)
    
        ! Ok, quadrilateral P1 element. Get from the triangulation the AREA-array for 
        ! calculating a simple integral of rx:
        
        CALL storage_getbase_double (p_rdiscretisation%p_rtriangulation%h_DelementVolume, &
                                    p_DelementArea)
                                      
        ! Get the vector data of rx
        CALL lsyssc_getbase_double (rx,p_Ddata)
        
        nel = SIZE(p_DelementArea)-1
                                      
        ! Build the integral
        !   int_Omega p dx
        ! This is approximated by
        !   dpintegral = SUM_Elements P(Element)*Volume(Element)
        ! Because taking the value in the element midpoint approximates the
        ! integral exactly by means of the 1x1 Gauss formula, the implementation
        ! is the same as in the Q0 and P0 case, respectively.
        
        dpintegral=0.0_DP
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
        
      CASE DEFAULT

        CALL output_line('Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_normaliseToL20Sca')
        CALL sys_halt()
        
      END SELECT
      
    ELSE
    
      CALL output_line('Unsupported discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_normaliseToL20Sca')
      CALL sys_halt()
      
    END IF

  END SUBROUTINE

! ***************************************************************************
! Block vector filters
! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposePressureDropBC (rx,dtimeweight,rpdbcStructure)
  
!<description>
  ! Implements discrete pressure drop BC's into a block vector.
!</description>

!<input>
  ! The t_discreteBCpressureDrop that describes the discrete pressure 
  ! drop BC's
  TYPE(t_discreteBCpressureDrop), INTENT(IN), TARGET  :: rpdbcStructure
  
  ! Weighting factor for time-dependent problems.
  ! =1.0 for stationary simulation.
  REAL(DP), INTENT(IN) :: dtimeweight
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
      CALL output_line('Pressure drop BC not configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposePressureDropBC')
      CALL sys_halt()
    END IF
    
    ! Currently, only 2D is supported. 
    ! First handle the X-velocity (j=1), then the Y-velocity (j=2)
    DO j=1,NDIM2D
    
      ! Get the subvector
      icp = rpdbcStructure%Icomponents(j)
      CALL lsyssc_getbase_double (rx%RvectorBlock(icp), p_vec)
    
      IF (.NOT.ASSOCIATED(p_vec)) THEN
        CALL output_line('No vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposePressureDropBC')
        CALL sys_halt()
      END IF

      ! Only handle nDOF DOF's, not the complete array!
      ! Probably, the array is longer (e.g. has the length of the vector), but
      ! contains only some entries...

      ! Is the vector sorted?
      IF (rx%RvectorBlock(j)%isortStrategy .LE. 0) THEN
        ! No. Implement directly.
        DO i=1,rpdbcStructure%nDOF
          p_vec(p_idx(i)) = p_vec(p_idx(i)) + p_val(j,i)*dtimeweight
        END DO
      ELSE
        ! Oops, vector sorted. At first get the permutation how its sorted
        ! - or more precisely, the back-permutation, as we need this one for 
        ! the loop below.
        CALL storage_getbase_int (rx%RvectorBlock(j)%h_IsortPermutation,p_Iperm)
        p_Iperm => p_Iperm(rx%RvectorBlock(j)%NEQ+1:)

        ! And 'filter' each DOF during the boundary value implementation!
        DO i=1,rpdbcStructure%nDOF
          p_vec(p_Iperm(p_idx(i))) = p_vec(p_Iperm(p_idx(i))) + &
                                     p_val(j,i) * dtimeweight
        END DO
      END IF
      
    END DO
  
  END SUBROUTINE
  
! *****************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeNLSlipDefectBC (rx,rslipBCStructure)
  
!<description>
  ! Implements discrete nonlinear slip BC's into a scalar defect vector
  ! as configured in the slip BC structure.
  ! This routine performs a special filtering to the defect vector
  ! of the type $r_m := r_m - (n*r_m)*n$ as described in
  ! [Kuzmin, Turek, Haario: Finite element simulation of turbulent
  ! bubble flows in gas-liquid reactors. Technical Report 298,
  ! September 2005, Chair of mathematics III, University of Dortmund]
!</description>

!<input>
  ! The t_discreteBCSlip that describes the discrete Dirichlet BC's
  TYPE(t_discreteBCSlip), INTENT(IN), TARGET  :: rslipBCStructure
!</input>

!<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT), TARGET :: rx
  
!</inputoutput>
  
!</subroutine>
    
    ! local variables
    INTEGER(PREC_DOFIDX) :: i,idof
    REAL(DP), DIMENSION(:), POINTER :: p_vecX,p_vecY
    INTEGER(I32), DIMENSION(:), POINTER :: p_idx
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dnormals
    REAL(DP) :: d
    
    ! If nDOF=0, there are no DOF's the current boundary condition segment,
    ! so we don't have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    IF (rslipBCStructure%nDOF .EQ. 0) RETURN
    
    ! Only 2D supported at the moment
    IF (rslipBCStructure%ncomponents .NE. NDIM2D) THEN
      CALL output_line('Only 2D supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      CALL sys_halt()
    END IF
    
    ! Only double precision vectors supported.
    IF (rx%cdataType .NE. ST_DOUBLE) THEN 
      CALL output_line('Only double precision supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      CALL sys_halt()
    END IF

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.
    
    IF (rslipBCStructure%h_IslipDOFs .EQ. ST_NOHANDLE) THEN
      CALL output_line('Slip-BC not configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      CALL sys_halt()
    END IF

    CALL storage_getbase_int(rslipBCStructure%h_IslipDOFs,p_idx)
    
    IF (.NOT.ASSOCIATED(p_idx)) THEN
      CALL output_line('Slip-BC not configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      CALL sys_halt()
    END IF
    
    IF (rx%RvectorBlock(rslipBCStructure%Icomponents(1))%isortStrategy .NE.&
        rx%RvectorBlock(rslipBCStructure%Icomponents(2))%isortStrategy) THEN
      CALL output_line('Subectors differently sorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      CALL sys_halt()
    END IF
    
    CALL lsyssc_getbase_double ( &
           rx%RvectorBlock(rslipBCStructure%Icomponents(1)), p_vecX)
    CALL lsyssc_getbase_double ( &
           rx%RvectorBlock(rslipBCStructure%Icomponents(2)), p_vecY)
    
    IF ( (.NOT.ASSOCIATED(p_vecX)) .OR. (.NOT.ASSOCIATED(p_vecX)) )THEN
      CALL output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      CALL sys_halt()
    END IF

    CALL storage_getbase_double2d(rslipBCStructure%h_DnormalVectors,p_Dnormals)

    ! Impose the BC-DOF's directly - more precisely, into the
    ! components of all the subvectors.
    !
    ! Only handle nDOF DOF's, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted? If yes, all vectors are sorted the same way!
    IF (rx%RvectorBlock(1)%isortStrategy .LE. 0) THEN
      ! No. Implement directly.
      DO i=1,rslipBCStructure%nDOF
        ! Get the DOF:
        idof = p_idx(i)

        ! Build n*r
        d = p_Dnormals(1,i)*p_vecX(idof) + p_Dnormals(2,i)*p_vecY(idof)
        
        ! Compute: r := r - (n*r)*n
        p_vecX(idof) = p_vecX(idof) - d*p_Dnormals(1,i)
        p_vecY(idof) = p_vecY(idof) - d*p_Dnormals(2,i)
      END DO
    ELSE
      ! Ups, vector sorted. At first get the permutation how its sorted -
      ! or more precisely, the back-permutation, as we need this one for 
      ! the loop below.
      CALL storage_getbase_int (rx%RvectorBlock(1)%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)
      
      ! And 'filter' each DOF during the boundary value implementation!
      DO i=1,rslipBCStructure%nDOF
        ! Get the DOF:
        idof = p_Iperm(p_idx(i))

        ! Build n*r
        d = p_Dnormals(1,i)*p_vecX(idof) + p_Dnormals(2,i)*p_vecY(idof)
        
        ! Compute: r := r - (n*r)*n
        p_vecX(idof) = p_vecX(idof) - d*p_Dnormals(1,idof)
        p_vecY(idof) = p_vecY(idof) - d*p_Dnormals(2,idof)
      END DO
    END IF
  
  END SUBROUTINE

! *****************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeFeastMirrorBC (rx,rfmbcStructure)
  
!<description>
  ! Implements discrete Feast Mirror BC's into a scalar vector.
!</description>

!<input>
  ! The t_discreteBCFeastMirror that describes the discrete Feast Mirror BC's
  TYPE(t_discreteBCFeastMirror), INTENT(IN), TARGET  :: rfmbcStructure
!</input>

!<inputoutput>
  ! The scalar vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT), TARGET :: rx
!</inputoutput>
  
!</subroutine>
    
  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_ImirrorDOFs
  INTEGER :: i
  REAL(DP), DIMENSION(:), POINTER    :: p_Dvec
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  REAL(DP) :: dmirrorWeight

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF (rx%cdataType .NE. ST_DOUBLE) THEN
    CALL output_line('Matrix must be double precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorBC')
    CALL sys_halt()
  END IF
  
  IF (rfmbcStructure%icomponent .EQ. 0) THEN
    CALL output_line('FMBC not configured!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorBC')
    CALL sys_halt()
  END IF
  
  IF (rfmbcStructure%h_ImirrorDOFs .EQ. ST_NOHANDLE) THEN
    ! No data inside of this structure.
    ! May happen if the region is not large enough to cover at least one DOF.
    RETURN
  END IF
  
  ! Get the weight of the entries.
  ! =2 on finest level, =1.5 on level NLMAX-1,...
  !dmirrorWeight = 1.0_DP+REAL(4**rfmbcStructure%icoarseningLevel,DP)
  dmirrorWeight = 1.0_DP+1.0_DP*REAL(2**rfmbcStructure%icoarseningLevel,DP)
  
  ! Get the vector data
  CALL lsyssc_getbase_double (rx%RvectorBlock(rfmbcStructure%icomponent),p_Dvec)
  
  ! Get pointers to the list of DOF's that belong to that region and have
  ! to be tackled.
  CALL storage_getbase_int(rfmbcStructure%h_ImirrorDOFs,p_ImirrorDOFs)

  IF ((rfmbcStructure%isubtype .EQ. 1) .OR. (rfmbcStructure%isubtype .EQ. 3)) THEN
    ! The DOF's should be treated as Dirichlet-DOF's.
    ! Only the matrix is modified according to the FEAST mirror bondary conditions!
    !
    ! For the implementation, we just set dmirrorWeight to 0.0.
    ! This clears all DOF entries and thus treats the DOF's like Dirichlet DOF's.
    dmirrorWeight = 0.0_DP
  END IF
  
  ! The vector entry corresponds to the DOF. For every DOF decide on
  ! whether it's on the FEAST mirror boundary component or not.
  ! If yes, double the entry entry.
  
  ! Is the vector sorted?
  IF (rx%RvectorBlock(rfmbcStructure%icomponent)%isortStrategy .LE. 0) THEN
    
    ! Loop through the DOF's. Each DOF gives us the number of an entry
    ! which is to be doubled.
    DO i=1,SIZE(p_ImirrorDOFs)
      p_Dvec(p_ImirrorDOFs(i)) = dmirrorWeight * p_Dvec(p_ImirrorDOFs(i))
    END DO
    
  ELSE
  
    ! Ok, vector is sorted, so we have to filter all the DOF's through the
    ! permutation before using them for implementing boundary conditions.
    !
    ! Get the permutation (or more precisely, the inverse permutation)
    ! from the vector to renumber the columns into
    ! the actual DOF numbers.
    CALL storage_getbase_int (&
        rx%RvectorBlock(rfmbcStructure%icomponent)%h_IsortPermutation,p_Iperm)
    p_Iperm => p_Iperm(rx%RvectorBlock(rfmbcStructure%icomponent)%NEQ+1:)
    
    ! Loop through the DOF's. Each DOF gives us the number of an entry
    ! which is to be doubled.
    DO i=1,SIZE(p_ImirrorDOFs)
      p_Dvec(p_Iperm(p_ImirrorDOFs(i))) = dmirrorWeight * p_Dvec(p_Iperm(p_ImirrorDOFs(i)))
    END DO

  END IF
  
  END SUBROUTINE
  
! *****************************************************************************

!<subroutine>

  SUBROUTINE vecfil_imposeFeastMirrorDefBC (rx,rfmbcStructure)
  
!<description>
  ! Implements discrete Feast Mirror BC's into a scalar defect vector.
!</description>

!<input>
  ! The t_discreteBCFeastMirror that describes the discrete Feast Mirror BC's
  TYPE(t_discreteBCFeastMirror), INTENT(IN), TARGET  :: rfmbcStructure
!</input>

!<inputoutput>
  ! The scalar vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT), TARGET :: rx
!</inputoutput>
  
!</subroutine>
    
  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_ImirrorDOFs
  INTEGER :: i
  REAL(DP), DIMENSION(:), POINTER    :: p_Dvec
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.
  
  IF (rx%cdataType .NE. ST_DOUBLE) THEN
    CALL output_line('Matrix must be double precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorDefBC')
    CALL sys_halt()
  END IF
  
  IF (rfmbcStructure%icomponent .EQ. 0) THEN
    CALL output_line('FMBC not configured!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorDefBC')
    CALL sys_halt()
  END IF
  
  IF (rfmbcStructure%h_ImirrorDOFs .EQ. ST_NOHANDLE) THEN
    ! No data inside of this structure.
    ! May happen if the region is not large enough to cover at least one DOF.
    RETURN
  END IF
  
  ! Get the vector data
  CALL lsyssc_getbase_double (rx%RvectorBlock(rfmbcStructure%icomponent),p_Dvec)
  
  ! Get pointers to the list of DOF's that belong to that region and have
  ! to be tackled.
  CALL storage_getbase_int(rfmbcStructure%h_ImirrorDOFs,p_ImirrorDOFs)

  IF ((rfmbcStructure%isubtype .EQ. 1) .OR. (rfmbcStructure%isubtype .EQ. 3)) THEN
    ! The DOF's should be treated as Dirichlet-DOF's.
    !
    ! The vector entry corresponds to the DOF. For every DOF decide on
    ! whether it's on the FEAST mirror boundary component or not.
    ! If yes, put the entry to zero.
    
    ! Is the vector sorted?
    IF (rx%RvectorBlock(rfmbcStructure%icomponent)%isortStrategy .LE. 0) THEN
      
      ! Loop through the DOF's. Each DOF gives us the number of an entry
      ! which is to be doubled.
      DO i=1,SIZE(p_ImirrorDOFs)
        p_Dvec(p_ImirrorDOFs(i)) = 0.0_DP
      END DO
      
    ELSE
    
      ! Ok, vector is sorted, so we have to filter all the DOF's through the
      ! permutation before using them for implementing boundary conditions.
      !
      ! Get the permutation (or more precisely, the inverse permutation)
      ! from the vector to renumber the columns into
      ! the actual DOF numbers.
      CALL storage_getbase_int (&
          rx%RvectorBlock(rfmbcStructure%icomponent)%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%RvectorBlock(rfmbcStructure%icomponent)%NEQ+1:)
      
      ! Loop through the DOF's. Each DOF gives us the number of an entry
      ! which is to be doubled.
      DO i=1,SIZE(p_ImirrorDOFs)
        p_Dvec(p_Iperm(p_ImirrorDOFs(i))) = 0.0_DP
      END DO

    END IF
    
  END IF
  
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
      IF (ASSOCIATED(rx%p_rdiscreteBC)) THEN
        p_RdiscreteBC => rx%p_rdiscreteBC%p_RdiscBCList  
      ELSE
        ! No BC
        NULLIFY(p_RdiscreteBC)
      END IF
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
    !DO i=1,SIZE(p_RdiscreteBC)
    DO i=1, rx%p_rdiscreteBC%inumEntriesUsed
    
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

      CASE (DISCBC_TPSLIP)  
        ! Nothing to do
        
      CASE (DISCBC_TPFEASTMIRROR)
        ! Nothing to do
        
      CASE DEFAULT
        CALL output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCsol')
                
        CALL sys_halt()
        
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
    !DO i=1,SIZE(p_RdiscreteBC)
    DO i=1, rx%p_rdiscreteBC%inumEntriesUsed
    
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
        ! Nothing to do.
        
      CASE (DISCBC_TPSLIP)
        ! Nothing to do.
        
      CASE (DISCBC_TPFEASTMIRROR)
        CALL vecfil_imposeFeastMirrorBC (rx,p_RdiscreteBC(i)%rfeastMirrorBCs)
        
      CASE DEFAULT
        CALL output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCrhs')
                
        CALL sys_halt()
        
      END SELECT
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteNLPDropBCrhs (rx,dtimeWeight,rdiscreteBC)

!<description>
  ! Implements discrete pressure drop BC's into a block RHS vector.
  ! To be called inside of a nonlinear or time-stepping loop.
  ! This routine performs a special filtering to the RHS vector
  ! of the type $r_m := r_m - \sum_j P_j \int_{S_j} \phi n ds$ as described 
  ! on page 257 (235) Turek's book.
  !
  ! The filtering is applied to the boundary components configured
  ! for pressure drop when setting up the BC's.
!</description>
  
!<input>
  ! OPTIONAL: Time-step weight. This weight is multiplied to time-dependent
  ! boundary conditions before these are added to the vector rx.
  ! The parameter can be omitted in stationary simulations.
  REAL(DP), INTENT(IN), OPTIONAL :: dtimeWeight

  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the vector.
  TYPE(t_discreteBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: i
    REAL(DP) :: dtweight
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC

    IF (.NOT. PRESENT(rdiscreteBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      IF (ASSOCIATED(rx%p_rdiscreteBC)) THEN
        p_RdiscreteBC => rx%p_rdiscreteBC%p_RdiscBCList  
      ELSE
        ! No BC
        NULLIFY(p_RdiscreteBC)
      END IF
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
    !DO i=1,SIZE(p_RdiscreteBC)
    DO i=1, rx%p_rdiscreteBC%inumEntriesUsed
    
      ! Only implement discrete pressure drop BC's.
      IF (p_RdiscreteBC(i)%itype .EQ. DISCBC_TPPRESSUREDROP) THEN
        CALL vecfil_imposePressureDropBC (rx,dtweight,&
                                          p_RdiscreteBC(i)%rpressureDropBCs)        
      END IF
      
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block defect vectors
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

    INTEGER :: iblock,i !,icp
    REAL(DP) :: dtweight
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC

    IF (.NOT. PRESENT(rdiscreteBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      IF (ASSOCIATED(rx%p_rdiscreteBC)) THEN
        p_RdiscreteBC => rx%p_rdiscreteBC%p_RdiscBCList  
      ELSE
        ! No BC
        NULLIFY(p_RdiscreteBC)
      END IF
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
    !DO i=1,SIZE(p_RdiscreteBC)
    DO i=1, rx%p_rdiscreteBC%inumEntriesUsed
    
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
        
      CASE (DISCBC_TPSLIP)
        ! Slip boundary conditions in the linear case are implemented
        ! in a nonlinear loop - so there's nothing to do here.
                
      CASE (DISCBC_TPFEASTMIRROR)
        ! Routine is on purpose not commented in! Not used for now!
        ! CALL vecfil_imposeFeastMirrorDefBC (rx,p_RdiscreteBC(i)%rfeastMirrorBCs)
        
      CASE DEFAULT
        CALL output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCdef')
        CALL sys_halt()
        
      END SELECT
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteNLSlipBCdef (rx,dtimeWeight,rdiscreteBC)

!<description>
  ! Implements discrete nonlinear slip BC's into a defect vector.
  ! Nonlinear filter, to be called inside of a nonlinear loop.
  ! This routine performs a special filtering to the defect vector
  ! of the type $r_m := r_m - (n*r_m)*n$ as described in
  ! [Kuzmin, Turek, Haario: Finite element simulation of turbulent
  ! bubble flows in gas-liquid reactors. Technical Report 298,
  ! September 2005, Chair of mathematics III, University of Dortmund]
  !
  ! The filtering is applied to all boundary components configured
  ! as slip in rx or rdiscreteBC (if given), respectively.
!</description>
  
!<input>
  ! OPTIONAL: Time-step weight. This weight is multiplied to time-dependent
  ! boundary conditions before these are added to the vector rx.
  ! The parameter can be omitted in stationary simulations.
  REAL(DP), INTENT(IN), OPTIONAL :: dtimeWeight

  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the vector.
  TYPE(t_discreteBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: i
    REAL(DP) :: dtweight
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscreteBC

    IF (.NOT. PRESENT(rdiscreteBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      IF (ASSOCIATED(rx%p_rdiscreteBC)) THEN
        p_RdiscreteBC => rx%p_rdiscreteBC%p_RdiscBCList  
      ELSE
        ! No BC
        NULLIFY(p_RdiscreteBC)
      END IF
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
    !DO i=1,SIZE(p_RdiscreteBC)
    DO i=1, rx%p_rdiscreteBC%inumEntriesUsed
    
      ! Only implement discrete slip BC's.
      IF (p_RdiscreteBC(i)%itype .EQ. DISCBC_TPSLIP) THEN
        CALL vecfil_imposeNLSlipDefectBC (rx,p_RdiscreteBC(i)%rslipBCs)          
      END IF
      
    END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into 
  ! block solution vectors
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteFBCsol (rx,dtimeWeight,rdiscreteFBC)

!<description>
  ! This routine realises the 'impose discrete fictitious boundary conditions 
  ! to solution' filter. 
  ! This filter imposes the discrete fictitious boundary conditions rdiscreteFBC
  ! (if specified) or (if rdiscreteFBC is not specified) the boundary conditions
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
  ! If not specified, the default fictitious boundary conditions associated 
  ! to the vector rx are imposed to the matrix.
  TYPE(t_discreteFBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,i,j
    REAL(DP) :: dtweight
    TYPE(t_discreteFBCEntry), DIMENSION(:), POINTER :: p_RdiscreteFBC

    IF (.NOT. PRESENT(rdiscreteFBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      IF (ASSOCIATED(rx%p_rdiscreteBCfict)) THEN
        p_RdiscreteFBC => rx%p_rdiscreteBCfict%p_RdiscFBCList  
      ELSE
        ! No BC
        NULLIFY(p_RdiscreteFBC)
      END IF
    ELSE
      p_RdiscreteFBC => rdiscreteFBC%p_RdiscFBCList
    END IF
    
    IF (.NOT. ASSOCIATED(p_RdiscreteFBC)) RETURN
    
    ! If the time-weight is not specified, 1.0 is assumed.
    IF (PRESENT(dtimeWeight)) THEN
      dtweight = dtimeWeight
    ELSE
      dtweight = 1.0_DP
    END IF
    ! Note: Time-step weight not used by any filter up to now!
    ! Perhaps in a later implementation it's needed anywhere...
    
    ! Now loop through all entries in this list:
    DO i=1,SIZE(p_RdiscreteFBC)
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteFBC(i)%itype)
      CASE (DISCFBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCFBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! Loop through all components, these boundary conditions should apply to.
        
        DO j=1,p_RdiscreteFBC(i)%rdirichletFBCs%ncomponents
          iblock = p_RdiscreteFBC(i)%rdirichletFBCs%Icomponents(j)
          ! Implement the Dirichlet boundary conditions into that component
          ! of the vector.
          CALL vecfil_imposeDirichletFBC (rx%RvectorBlock(iblock),j,&
                                          p_RdiscreteFBC(i)%rdirichletFBCs)
        END DO

      CASE (DISCBC_TPSLIP)
        ! Nothing to do.
        
      CASE (DISCBC_TPFEASTMIRROR)
        ! Nothing to do
        
      CASE DEFAULT
        CALL output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteFBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCsol')
        CALL sys_halt()
        
      END SELECT
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into 
  ! block RHS vectors
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteFBCrhs (rx,dtimeWeight,rdiscreteFBC)

!<description>
  ! This routine realises the 'impose discrete fictitious boundary conditions 
  ! to RHS' filter. 
  ! This filter imposes the discrete fictitious boundary conditions rdiscreteFBC
  ! (if specified) or (if rdiscreteFBC is not specified) the boundary conditions
  ! which are  associated to the vector rx (with rx%p_discreteFBC) to this 
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
  TYPE(t_discreteFBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,i,j
    REAL(DP) :: dtweight
    TYPE(t_discreteFBCEntry), DIMENSION(:), POINTER :: p_RdiscreteFBC

    IF (.NOT. PRESENT(rdiscreteFBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      IF (ASSOCIATED(rx%p_rdiscreteBCfict)) THEN
        p_RdiscreteFBC => rx%p_rdiscreteBCfict%p_RdiscFBCList  
      ELSE
        ! No BC
        NULLIFY(p_RdiscreteFBC)
      END IF
    ELSE
      p_RdiscreteFBC => rdiscreteFBC%p_RdiscFBCList
    END IF
    
    IF (.NOT. ASSOCIATED(p_RdiscreteFBC)) RETURN
    
    ! If the time-weight is not specified, 1.0 is assumed.
    IF (PRESENT(dtimeWeight)) THEN
      dtweight = dtimeWeight
    ELSE
      dtweight = 1.0_DP
    END IF
    
    ! Now loop through all entries in this list:
    DO i=1,SIZE(p_RdiscreteFBC)
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteFBC(i)%itype)
      CASE (DISCFBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCFBC_TPDIRICHLET)
        ! Loop through all components, these boundary conditions should apply to.
        
        DO j=1,p_RdiscreteFBC(i)%rdirichletFBCs%ncomponents
          iblock = p_RdiscreteFBC(i)%rdirichletFBCs%Icomponents(j)
          
          ! Implement the Dirichlet boundary conditions into that component
          ! of the vector.
          CALL vecfil_imposeDirichletFBC (rx%RvectorBlock(iblock),j,&
                                          p_RdiscreteFBC(i)%rdirichletFBCs)
        END DO
      
      CASE (DISCBC_TPSLIP)
        ! Nothing to do.
        
      CASE (DISCBC_TPFEASTMIRROR)
        ! Nothing to do
        
      CASE DEFAULT
        CALL output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteFBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteFBCrhs')
        CALL sys_halt()
        
      END SELECT
    END DO
  
  END SUBROUTINE

  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into 
  ! block defect vectors
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vecfil_discreteFBCdef (rx,dtimeWeight,rdiscreteFBC)

!<description>
  ! This routine realises the 'impose discrete fictitious boundary conditions 
  ! to defect' filter. 
  ! This filter imposes the discrete fictitious boundary conditions rdiscretFeBC
  ! (if specified) or (if rdiscreteFBC is not specified) the boundary conditions
  ! which are  associated to the defect vector rx (with rx%p_discreteFBC) to 
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
  TYPE(t_discreteFBC), OPTIONAL, INTENT(IN), TARGET :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  TYPE(t_vectorBlock), INTENT(INOUT),TARGET :: rx
!</inputoutput>

!</subroutine>

    INTEGER :: iblock,i,j
    REAL(DP) :: dtweight
    TYPE(t_discreteFBCEntry), DIMENSION(:), POINTER :: p_RdiscreteFBC

    IF (.NOT. PRESENT(rdiscreteFBC)) THEN
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      IF (ASSOCIATED(rx%p_rdiscreteBCfict)) THEN
        p_RdiscreteFBC => rx%p_rdiscreteBCfict%p_RdiscFBCList  
      ELSE
        ! No BC.
        NULLIFY(p_RdiscreteFBC)
      END IF
    ELSE
      p_RdiscreteFBC => rdiscreteFBC%p_RdiscFBCList
    END IF
    
    IF (.NOT. ASSOCIATED(p_RdiscreteFBC)) RETURN
    
    ! If the time-weight is not specified, 1.0 is assumed.
    IF (PRESENT(dtimeWeight)) THEN
      dtweight = dtimeWeight
    ELSE
      dtweight = 1.0_DP
    END IF
    ! Note: Time-step weight not used by any filter up to now!
    ! Perhaps in a later implementation it's needed anywhere...
    
    ! Now loop through all entries in this list:
    DO i=1,SIZE(p_RdiscreteFBC)
    
      ! What for BC's do we have here?
      SELECT CASE (p_RdiscreteFBC(i)%itype)
      CASE (DISCFBC_TPUNDEFINED)
        ! Do-nothing
        
      CASE (DISCFBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        
        ! Loop over all blocks where to implement these FBC's.
        DO j=1,p_RdiscreteFBC(i)%rdirichletFBCs%ncomponents
          iblock = p_RdiscreteFBC(i)%rdirichletFBCs%Icomponents(j)
        
          ! Implement the Dirichlet boundary conditions into that component
          ! of the vector.
          CALL vecfil_imposeDirichletDefectFBC (rx%RvectorBlock(iblock),&
                                                p_RdiscreteFBC(i)%rdirichletFBCs)
        END DO
        
      CASE (DISCBC_TPSLIP)
        ! Nothing to do.
        
      CASE (DISCBC_TPFEASTMIRROR)
        ! Nothing to do
        
      CASE DEFAULT
        CALL output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteFBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteFBCdef')
        CALL sys_halt()
        
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

    IF ((isubvector .LE. 0) .OR. (isubvector .GT. SIZE(rx%RvectorBlock))) THEN
      CALL output_line('isubvector out of allowed range!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_subvectorToL20')
      CALL sys_halt()
    END IF
      
    ! Don't use one single IF here as this may lead to errors depending
    ! on the compiler (subvector=0 and access to vector(isubvector)).
    IF ((rx%RvectorBlock(isubvector)%h_Ddata .EQ. ST_NOHANDLE) .OR. &
        (rx%RvectorBlock(isubvector)%NEQ .LE. 0)) &
      RETURN
        
    ! Normalise the subvector isubvector
    CALL vecfil_normaliseToL20Sca (rx%RvectorBlock(isubvector))

  END SUBROUTINE

END MODULE
