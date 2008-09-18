!##############################################################################
!# ****************************************************************************
!# <name> timeboundaryconditions </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to assemble and implement boundary conditions
!# into space-time vectors.
!#
!# The central routines in this module are:
!#
!# 1.) tbc_implementBCsolution
!#     -> Implements boundary conditions into a given space-time solution vector.
!#
!# 2.) tbc_implementBCRHS
!#     -> Implements boundary conditions into a given space-time RHS vector.
!#
!# 3.) tbc_implementBCdefect
!#     -> Implements initial and boundary conditions into a fiven space-time
!#        defect vector.
!#
!# 4.) tbc_implementInitCondRHS
!#     -> Implements initial conditions into a given space-time RHS vector.
!#
!# 5.) tbc_implementInitCondDefect
!#     -> Implements initial conditions into a given space-time defect vector.
!#
!# 6.) tbc_pressureToL20
!#     -> Normalise the primal and dual pressure to integral mean value = 0.
!#
!# Auxiliary routines:
!#
!# 1.) tbc_implementInitCondDefSingle
!#     -> Implement the initial condition into a spatial vector
!#
!# 2.) tbc_implementTermCondDefSingle
!#     -> Implement the terminal condition into a spatial vector
!#
!# 3.) tbc_implementSpatialBCtoRHS
!#     -> Assembles and implements the boundary conditions at a given point
!#        in time into a given spatial RHS vector.
!#
!# 4.) tbc_implementSpatialBCdefect
!#     -> Assembles and implements the boundary conditions at a given point
!#        in time into a given spatial defect vector.
!#
!# </purpose>
!##############################################################################

MODULE timeboundaryconditions


  USE fsystem
  USE storage
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
  USE paramlist
  USE timestepping
  USE l2projection
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback

  USE cc2dmediumm2nonlinearcore
  USE cc2dmediumm2nonlinearcoreinit
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  
  USE timediscretisation
  USE spacetimevectors
  USE dofmapping
  
  USE spacetimediscretisation

  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementBCsolution (rproblem,rspaceTimeDiscr,rx,rtempVectorX)

!<description>
  ! Implements the boundary conditions of all timesteps into the solution rx.
!</description>

!<input>
  ! Problem structure of the main problem.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented
  ! to.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rtempVectorX
!</inputoutput>

!</subroutine>

    INTEGER :: isubstep
    REAL(DP) :: dtstep
    TYPE(t_vectorBlock) :: rtempVector
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    IF (PRESENT(rtempVectorX)) THEN
      CALL lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    ELSE
      ! Create a temp vector
      CALL lsysbl_createVecBlockByDiscr (&
          rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.TRUE.)
    END IF
    rtempVector%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC

    CALL lsyssc_getbase_double (rtempVector%RvectorBlock(1),p_Ddata)

    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep

    ! The implementation of the boundary conditions depends on the type
    ! of the time discretisation...
    SELECT CASE (rspaceTimeDiscr%rtimeDiscr%ctype)
    CASE (TDISCR_THETA)

      ! Implement the bondary conditions into all initial solution vectors
      DO isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + &
            isubstep*dtstep
        rproblem%rtimedependence%itimestep = isubstep

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF
        
        ! Implement the boundary conditions into the global solution vector.
        CALL sptivec_getTimestepData(rx, 1+isubstep, rtempVector)
        
        CALL cc_implementBC (rproblem,rvector=rtempVector)
        
        CALL sptivec_setTimestepData(rx, 1+isubstep, rtempVector)
        
      END DO
      
    CASE (TDISCR_DG0)

      ! Implement the bondary conditions into all initial solution vectors
      DO isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals-1
      
        ! Current point in time
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + &
            (REAL(isubstep,DP)+0.5_DP)*dtstep
        rproblem%rtimedependence%itimestep = isubstep

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF
        
        ! Implement the boundary conditions into the global solution vector.
        CALL sptivec_getTimestepData(rx, 1+isubstep, rtempVector)
        
        CALL cc_implementBC (rproblem,rvector=rtempVector)
        
        CALL sptivec_setTimestepData(rx, 1+isubstep, rtempVector)
        
      END DO

    CASE DEFAULT
        
      CALL output_line ('Unsupported time discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tbc_implementBCsolution')
      CALL sys_halt()
    
    END SELECT
    
    ! Release the temp vector
    CALL lsysbl_releaseVector (rtempVector)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementBCRHS (rproblem,rspaceTimeDiscr,rb,rtempVectorX)

!<description>
  ! Implements the boundary conditions of all timesteps into the RHS vector rb.
!</description>

!<input>
  ! Problem structure of the main problem.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rtempVectorX
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented
  ! to.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
!</inputoutput>

!</subroutine>

    INTEGER :: isubstep
    REAL(DP) :: dtstep
    TYPE(t_vectorBlock) :: rtempVector
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    IF (PRESENT(rtempVectorX)) THEN
      CALL lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    ELSE
      ! Create a temp vector
      CALL lsysbl_createVecBlockByDiscr (&
          rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.TRUE.)
    END IF
    rtempVector%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC

    CALL lsyssc_getbase_double (rtempVector%RvectorBlock(1),p_Ddata)

    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep

    ! The implementation of the boundary conditions depends on the type
    ! of the time discretisation...
    SELECT CASE (rspaceTimeDiscr%rtimeDiscr%ctype)
    CASE (TDISCR_THETA)

      ! Implement the bondary conditions into all initial solution vectors
      DO isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + &
            isubstep*dtstep
        rproblem%rtimedependence%itimestep = isubstep

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF
        
        ! Implement the boundary conditions into the global solution vector.
        CALL sptivec_getTimestepData(rb, 1+isubstep, rtempVector)
        
        CALL cc_implementBC (rproblem,rvector=rtempVector)
        
        CALL sptivec_setTimestepData(rb, 1+isubstep, rtempVector)
        
      END DO
      
    CASE (TDISCR_DG0)

      ! Implement the bondary conditions into all initial solution vectors
      DO isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + &
            (REAL(isubstep,DP)+0.5_DP) * dtstep
        rproblem%rtimedependence%itimestep = isubstep

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF
        
        ! Implement the boundary conditions into the global solution vector.
        CALL sptivec_getTimestepData(rb, 1+isubstep, rtempVector)
        
        CALL cc_implementBC (rproblem,rvector=rtempVector)
        
        CALL sptivec_setTimestepData(rb, 1+isubstep, rtempVector)
        
      END DO

    CASE DEFAULT
        
      CALL output_line ('Unsupported time discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tbc_implementBCRHS')
      CALL sys_halt()
    
    END SELECT

    ! Release the temp vector
    CALL lsysbl_releaseVector (rtempVector)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementBCdefect (rproblem,rspaceTimeDiscr,rd,rtempVectorX)

!<description>
  ! Implements the boundary conditions of all timesteps into the defect rd.
!</description>

!<input>
  ! Problem structure of the main problem.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented
  ! to.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rtempVectorX
!</inputoutput>

!</subroutine>

    REAL(DP) :: dtstep
    INTEGER :: isubstep
    TYPE(t_vectorBlock) :: rtempVector
    
    IF (PRESENT(rtempVectorX)) THEN
      CALL lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    ELSE
      ! Create temp vectors
      CALL lsysbl_createVecBlockByDiscr (&
          rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.TRUE.)
    END IF
    rtempVector%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC

    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep

    ! The implementation of the boundary conditions depends on the type
    ! of the time discretisation...
    SELECT CASE (rspaceTimeDiscr%rtimeDiscr%ctype)
    CASE (TDISCR_THETA)

      ! Implement the bondary conditions into all initial solution vectors
      DO isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + &
            isubstep*dtstep
        rproblem%rtimedependence%itimestep = isubstep

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF
        
        ! Implement the boundary conditions into the global solution vector.
        CALL sptivec_getTimestepData(rd, 1+isubstep, rtempVector)
        
        CALL cc_implementBC (rproblem,rdefect=rtempVector)
        
        ! In the very first time step, we have the initial condition for the
        ! solution. The defect is =0 there!
        IF (isubstep .EQ. 0) THEN
          CALL lsyssc_clearVector (rtempVector%RvectorBlock(1))
          CALL lsyssc_clearVector (rtempVector%RvectorBlock(2))
          CALL lsyssc_clearVector (rtempVector%RvectorBlock(3))
        END IF
        
        CALL sptivec_setTimestepData(rd, 1+isubstep, rtempVector)
        
      END DO
    
    CASE (TDISCR_DG0)


      ! Implement the bondary conditions into all initial solution vectors
      DO isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + &
            (REAL(isubstep,DP)+0.5_DP) * dtstep
        rproblem%rtimedependence%itimestep = isubstep

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF
        
        ! Implement the boundary conditions into the global solution vector.
        CALL sptivec_getTimestepData(rd, 1+isubstep, rtempVector)
        
        CALL cc_implementBC (rproblem,rdefect=rtempVector)
        
        ! In the very first time step, we have the initial condition for the
        ! solution. The defect is =0 there!
        IF (isubstep .EQ. 0) THEN
          CALL lsyssc_clearVector (rtempVector%RvectorBlock(1))
          CALL lsyssc_clearVector (rtempVector%RvectorBlock(2))
          CALL lsyssc_clearVector (rtempVector%RvectorBlock(3))
        END IF
        
        CALL sptivec_setTimestepData(rd, 1+isubstep, rtempVector)
        
      END DO

    CASE DEFAULT
        
      CALL output_line ('Unsupported time discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tbc_implementBCdefect')
      CALL sys_halt()
    
    END SELECT

    CALL lsysbl_releaseVector(rtempVector)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_pressureToL20 (rx,rtempVectorX)

!<description>
  ! Normalises the primal and dual pressure in all time steps to have integral
  ! mean value = 0. This routine is typically used to filter an indefinite
  ! solution vector (e.g. in the pure-Dirichlet case).
!</description>

!<inputoutput>
  ! A space-time vector where the pressure vectors whould be normalised.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rtempVectorX
!</inputoutput>

!</subroutine>

    INTEGER :: isubstep
    TYPE(t_vectorBlock) :: rtempVector
    
    IF (PRESENT(rtempVectorX)) THEN
      CALL lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    ELSE
      ! Create temp vectors
      CALL lsysbl_createVecBlockByDiscr (&
          rx%p_rblockDiscretisation,rtempVector,.TRUE.)
    END IF

    ! Normalise the primal and dual pressure to zero.
    DO isubstep = 0,rx%NEQtime-1
    
      CALL sptivec_getTimestepData(rx, 1+isubstep, rtempVector)
      
      CALL vecfil_subvectorToL20 (rtempVectorX,3)
      CALL vecfil_subvectorToL20 (rtempVectorX,6)
      
      CALL sptivec_setTimestepData(rx, 1+isubstep, rtempVector)
      
    END DO
  
    CALL lsysbl_releaseVector(rtempVector)

  END SUBROUTINE

  ! *************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementInitCondRHS (rproblem,rb,rinitCondRHS,rtempVectorD)

!<description>
  ! Implements the initial condition into the RHS vector rb.
  ! Overwrites the rb of the first time step.
  !
  ! Does not implement boundary conditions!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</input>

!<inputoutput>
  ! A space-time vector with the RHS. The initial condition is implemented into
  ! this vector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb

  ! A vector containing the data for the initial condition of the RHS.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rinitCondRHS

  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD
!</inputoutput>

!</subroutine>

    REAL(DP) :: dtheta
    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd
    
    ! Overwrite the primal RHS with the initial primal solution vector.
    ! This realises the inital condition.
    CALL sptivec_getTimestepData(rb, 1+0, rtempVectorD)
    CALL lsyssc_copyVector (rinitCondRHS%RvectorBlock(1),rtempVectorD%RvectorBlock(1))
    CALL lsyssc_copyVector (rinitCondRHS%RvectorBlock(2),rtempVectorD%RvectorBlock(2))
    CALL sptivec_setTimestepData(rb, 1+0, rtempVectorD)
    
!    REAL(DP) :: dtheta
!    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd
!    
!    ! DEBUG!!!    
!    CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
!    CALL lsysbl_getbase_double (rtempVectorD,p_Db)
!
!    ! Theta-scheme identifier.
!    ! =1: impliciz Euler.
!    ! =0.5: Crank Nicolson
!    dtheta = rproblem%rtimedependence%dtimeStepTheta
!
!    ! Overwrite the primal RHS with the initial primal solution vector.
!    ! This realises the inital condition.
!    CALL sptivec_getTimestepData(rx, 1+0, rtempVectorX)
!    
!    CALL sptivec_getTimestepData(rb, 1+0, rtempVectorD)
!    CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(1),rtempVectorD%RvectorBlock(1))
!    CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(2),rtempVectorD%RvectorBlock(2))
!    CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(3),rtempVectorD%RvectorBlock(3))
!
!    ! Save the modified RHS.
!    CALL sptivec_setTimestepData(rb, 1+0, rtempVectorD)
    
    ! DEBUG!!!
    !CALL sptivec_getTimestepData(rx, 2, rtempVectorX)
    !CALL sptivec_getTimestepData(rb, 2, rtempVectorD)
    !CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(4),rtempVectorD%RvectorBlock(4))
    !CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(5),rtempVectorD%RvectorBlock(5))
    !CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(6),rtempVectorD%RvectorBlock(6))
    !CALL sptivec_setTimestepData(rb, 2, rtempVectorD)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementInitCond (rproblem,rx,rinitCondSol,rtempVector)

!<description>
  ! Implements the initial condition into the vector rx.
  ! Overwrites the rx of the first time step.
  !
  ! Does not implement boundary conditions!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</input>

!<inputoutput>
  ! A space-time vector with the RHS. The initial condition is implemented into
  ! this vector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A vector containing the data for the initial condition of the RHS.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rinitCondSol

  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector
!</inputoutput>

!</subroutine>

    ! Overwrite the primal solution with the initial primal solution vector.
    ! This realises the inital condition.
    CALL sptivec_getTimestepData(rx, 1+0, rtempVector)
    CALL lsyssc_copyVector (rinitCondSol%RvectorBlock(1),rtempVector%RvectorBlock(1))
    CALL lsyssc_copyVector (rinitCondSol%RvectorBlock(2),rtempVector%RvectorBlock(2))
    CALL sptivec_setTimestepData(rx, 1+0, rtempVector)
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementInitCondDefect (rspaceTimeDiscr, rd, rtempvectorD)

!<description>
  ! Implements the initial and terminal condition into a defect vector rd.
  ! Overwrites the rd of the first time step.
  !
  ! Does not implement boundary conditions!
!</description>

!<inputoutput>
  ! Discretisation structure that corresponds to rx.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr

  ! A space-time vector containing the defect in the first subvector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd

  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD
!</inputoutput>

!</subroutine>

    REAL(DP), DIMENSION(:),POINTER :: p_Db
    
    ! DEBUG!!!    
    CALL lsysbl_getbase_double (rtempVectorD,p_Db)
!
!    ! Overwrite the primal defect with 0 -- as the solution must not be changed.
!    ! This realises the inital condition.
!    CALL sptivec_getTimestepData(rd, 1+0, rtempVectorD)
!    CALL tbc_implementInitCondDefSingle (rspaceTimeDiscr, rtempVectorD)
!    CALL sptivec_setTimestepData(rd, 1+0, rtempVectorD)

!    REAL(DP), DIMENSION(:),POINTER :: p_Db
!    
!    ! DEBUG!!!    
!    CALL lsysbl_getbase_double (rtempVectorD,p_Db)
!
!    ! Overwrite the primal defect with 0 -- as the solution must not be changed.
!    ! This realises the inital condition.
!    CALL sptivec_getTimestepData(rd, 1+0, rtempVectorD)
!    CALL tbc_implementInitCondDefSingle (rspaceTimeDiscr, rtempVectorD)
!    CALL sptivec_setTimestepData(rd, 1+0, rtempVectorD)
!
!    IF (rspaceTimeDiscr%dgammaC .EQ. 0.0_DP) THEN
!      ! That's a special case, we have the terminal condition "lambda(T)=0".
!      ! This case must be treated like the initial condition, i.e. the
!      ! dual defect in the last timestep must be overwritten by zero.
!      !
!      ! If gamma<>0, the terminal condition is implemented implicitely
!      ! by the equation "lambda(T)=gamma(y(T)-z(T))" which comes into
!      ! play by the mass matrix term in the system matrix of the last timestep,
!      ! so this does not have to be treated explicitly.
!      
!      CALL sptivec_getTimestepData(rd, rd%NEQtime, rtempVectorD)
!      
!      CALL tbc_implementTermCondDefSingle (rspaceTimeDiscr, rtempvectorD)
!      
!      CALL sptivec_setTimestepData(rd, rd%NEQtime, rtempVectorD)
!      
!    END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementInitCondDefSingle (rspaceTimeDiscr, rd)

!<description>
  ! Implements the initial condition into a defect vector rd,
  ! representing the defect in the first timestep.
  !
  ! Does not implement boundary conditions!
!</description>

!<inputoutput>
  ! Discretisation structure that corresponds to rx.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr

  ! A vector containing the defect in the first subvector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

!    REAL(DP), DIMENSION(:),POINTER :: p_Db
!    
!    ! DEBUG!!!    
!    CALL lsysbl_getbase_double (rd,p_Db)
!    
!    CALL lsyssc_clearVector(rd%RvectorBlock(1))
!    CALL lsyssc_clearVector(rd%RvectorBlock(2))
!    CALL lsyssc_clearVector(rd%RvectorBlock(3))

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementTermCondDefSingle (rspaceTimeDiscr, rd)

!<description>
  ! Implements the terminal condition into a defect vector rd,
  ! representing the defect in the last timestep.
  !
  ! Does not implement boundary conditions!
!</description>

!<inputoutput>
  ! Discretisation structure that corresponds to rx.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr

  ! A vector containing the defect in the last subvector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

!    Note: In the current implementation, the terminal condition is
!    imposed weakly, therefore for following lines of code are not used!
!
!    REAL(DP), DIMENSION(:),POINTER :: p_Db
!    
!    ! DEBUG!!!    
!    CALL lsysbl_getbase_double (rd,p_Db)
!
!    IF (rspaceTimeDiscr%dgammaC .EQ. 0.0_DP) THEN
!      ! That's a special case, we have the terminal condition "lambda(T)=0".
!      ! This case must be treated like the initial condition, i.e. the
!      ! dual defect in the last timestep must be overwritten by zero.
!      !
!      ! If gamma<>0, the terminal condition is implemented implicitely
!      ! by the equation "lambda(T)=gamma(y(T)-z(T))" which comes into
!      ! play by the mass matrix term in the system matrix of the last timestep,
!      ! so this does not have to be treated explicitly.
!      !
!      ! These command implement the terminal condition in a strong sense.
!      ! By commenting these lines out, the terminal condition would be
!      ! implemented in a weak sense, as the equation implemented in
!      ! spacetimelinearsystem.f90 reads:
!      !
!      !     -gamma*M*y + (M+dt*nu*L)*lambda = -gamma*z
!      !
!      ! which adds a mass matrix to a 'smoothing' Laplace part.
!      
!      IF (rspaceTimeDiscr%itypeTerminalCondition .EQ. 0) THEN
!        CALL lsyssc_clearVector(rd%RvectorBlock(4))
!        CALL lsyssc_clearVector(rd%RvectorBlock(5))
!        CALL lsyssc_clearVector(rd%RvectorBlock(6))
!      END IF
!      
!    END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE tbc_implementSpatialBCtoRHS (rproblem, isubstep, dtime, rvector)
  
!<description>
  ! Implements the spatial boundary conditions into the spatial RHS vector
  ! rvector.
!</description>
 
!<input>
  ! Time where the BC's should be implemented.
  ! Must not necessarily coincide with the start/end time of the timestep.
  REAL(DP), INTENT(IN) :: dtime
    
  ! Number of the substep where to implement the BC.
  INTEGER, INTENT(IN) :: isubstep
!</input>  

!<inputoutput>  
  ! Problem structure.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! Source and destination RHS vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>
  
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rvector,p_Ddata)

    ! Set the time where we are at the moment
    rproblem%rtimedependence%dtime = dtime
    rproblem%rtimedependence%itimestep = isubstep
        
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Discretise the boundary conditions at the new point in time -- 
    ! if the boundary conditions are nonconstant in time!
    IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
      CALL cc_updateDiscreteBC (rproblem)
    END IF

    ! Implement the boundary conditions into the RHS.
    ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
    ! as Dirichlet values mustn't be multiplied with GAMMA!
    CALL vecfil_discreteBCsol (rvector)
    CALL vecfil_discreteFBCsol (rvector)      
  
    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE tbc_implementSpatialBCdefect (rproblem,isubstep,dtime,rspaceTimeDiscr,rd)

!<description>
  ! Implements the boundary conditions at timestep isubstep into the defect rd.
!</description>

!<input>
  ! Problem structure of the main problem.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr

  ! Time where the BC's should be implemented.
  ! Must not necessarily coincide with the start/end time of the timestep.
  REAL(DP), INTENT(IN) :: dtime
    
  ! Number of the substep where to implement the BC.
  INTEGER, INTENT(IN) :: isubstep
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented to.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    TYPE(t_vectorBlock) :: rtempVector
    
    CALL lsysbl_duplicateVector(rd,rtempVector,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    rtempVector%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC

    ! Current point in time
    rproblem%rtimedependence%dtime = dtime
    rproblem%rtimedependence%itimestep = isubstep

    ! -----
    ! Discretise the boundary conditions at the new point in time -- 
    ! if the boundary conditions are nonconstant in time!
    IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
      CALL cc_updateDiscreteBC (rproblem)
    END IF
    
    CALL cc_implementBC (rproblem,rdefect=rtempVector)
    
    ! In the very first time step, we have the initial condition for the
    ! solution. The defect is =0 there!
    IF (isubstep .EQ. 0) THEN
      CALL lsyssc_clearVector (rd%RvectorBlock(1))
      CALL lsyssc_clearVector (rd%RvectorBlock(2))
      CALL lsyssc_clearVector (rd%RvectorBlock(3))
    END IF
    
    CALL lsysbl_releaseVector(rtempVector)

  END SUBROUTINE

END MODULE
