!##############################################################################
!# ****************************************************************************
!# <name> timerhsevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to evaluate/create the coupled space-time
!# RHS-vector. 
!#
!# The central routines in this module are:
!#
!# 1.) trhsevl_assembleRHS
!#     -> Assemble a space-time RHS vector
!#
!# Auxiliary routines:
!#
!# 1.) trhsevl_assembleThetaRHS
!#     -> Assemble a space-time RHS vector according to a Theta scheme
!#
!# 2.) trhsevl_assembledG0RHS
!#     -> Assemble a space-time RHS vector according to the dG(0)-scheme
!#
!# 3.) trhsevl_assembleSpatialRHS
!#     -> Assembles the spatial RHS at a given point in time.
!#
!# </purpose>
!##############################################################################

MODULE timerhsevaluation


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
  USE cc2dmediumm2matvecassembly
  
  USE timediscretisation
  USE spacetimevectors
  USE dofmapping
  USE timeboundaryconditions
  
  USE spacetimediscretisation

  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)

!<description>
  ! Assembles the space-time RHS vector rb. 
  !
  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
  ! this routine!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time system.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
  
  ! Whether to implement boundary conditions into the RHS or not.
  LOGICAL, INTENT(IN) :: bimplementBC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  ! If this is undefined, a new space-time vector is created.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
!</inputoutput>

!</subroutine>

    IF (rb%NEQtime .EQ. 0) THEN
      ! Create a new vector if rb is undefined.
      CALL sptivec_initVectorDiscr (rb,rspaceTimeDiscr%rtimeDiscr,&
          rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
    END IF

    ! What's the current time discretisation? Depending on that,
    ! we have to call the corresponding RHS calculation routine.
    
    SELECT CASE (rspaceTimeDiscr%rtimeDiscr%ctype)
    CASE (TDISCR_THETA)
      CALL trhsevl_assembleThetaRHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)
    CASE (TDISCR_DG0)
      CALL trhsevl_assembledG0RHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)
    CASE DEFAULT
      CALL output_line ('Unsupported time discretisation', &
                        OU_CLASS_ERROR,OU_MODE_STD,'trhsevl_assembleRHS')
      CALL sys_halt()
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE trhsevl_assembleThetaRHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)

!<description>
  ! Assembles the space-time RHS vector rb for a Theta-Scheme. 
  !
  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
  ! this routine!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time system.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
  
  ! Whether to implement boundary conditions into the RHS or not.
  LOGICAL, INTENT(IN) :: bimplementBC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
 
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubstep,nintervals
    REAL(DP) :: dtheta,dtstep
    
    ! Temporary vectors
    TYPE(t_vectorBlock) :: rtempVector1,rtempVector2,rtempVector3

    ! A temporary vector for the creation of the RHS.
    TYPE(t_vectorBlock) :: rtempVectorRHS
    
    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd, p_Drhs

    ! Create temp vectors for the assembly
    CALL lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,rtempVector1,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,rtempVector2,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,rtempVector3,.TRUE.)

    ! Theta-scheme identifier.
    ! =1: impliciz Euler.
    ! =0.5: Crank Nicolson
    dtheta = rspaceTimeDiscr%rtimeDiscr%dtheta
    
    ! Time step size, number of intervals
    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep
    nintervals = rspaceTimeDiscr%rtimeDiscr%nintervals
    
    ! ----------------------------------------------------------------------
    ! Generate the global RHS vector
    
    CALL lsysbl_getbase_double (rtempVector1,p_Dx)
    CALL lsysbl_getbase_double (rtempVector2,p_Db)
    CALL lsysbl_getbase_double (rtempVector3,p_Dd)

    ! Assemble 1st RHS vector in X temp vector.
    CALL trhsevl_assembleSpatialRHS (rproblem,0,0.0_DP,rtempVector1)
        
    ! Assemble the 2nd RHS vector in the RHS temp vector
    CALL trhsevl_assembleSpatialRHS (rproblem,1,dtstep,rtempVector2)

    ! Assemble the 3rd RHS vector in the defect temp vector
    IF (rspaceTimeDiscr%rtimeDiscr%nintervals .GE. 2) THEN
      CALL trhsevl_assembleSpatialRHS (rproblem,2,2.0_DP*dtstep,rtempVector3)
    ELSE
      CALL lsysbl_copyVector (rtempVector2,rtempVector3)
    END IF
        
    ! Create a copy of the X temp vector (RHS0). That vector will be
    ! our destination vector for assembling the RHS in all timesteps.
    CALL lsysbl_copyVector (rtempVector1,rtempVectorRHS)
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
    
    ! RHS 0,1,2 -> 1-2-3
    
    DO isubstep = 0,nintervals
    
      IF (isubstep .EQ. 0) THEN
      
        ! Primal RHS comes from rtempVector1. The dual from the
        ! isubstep+1'th RHS in rtempVector2.
        !
        ! primal RHS(0) = PRIMALRHS(0)
        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)

        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(1),rtempVectorRHS%RvectorBlock(1))
        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(2),rtempVectorRHS%RvectorBlock(2))
        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))

        CALL lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(4),rtempVector2%RvectorBlock(4),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(4))
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(5),rtempVector2%RvectorBlock(5),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(5))
        ! Pressure is fully implicit
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(6),rtempVector2%RvectorBlock(6),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(6))

        ! In the 0'th timestep, there is no RHS in the dual equation!
        ! That is because of the initial condition, which fixes the primal solution
        ! => dual solution has no influence on the primal one
        ! => setting up a dual RHS in not meaningful as the dual RHS cannot
        !    influence the primal solution
        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(4))
        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(5))
        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(6))
            
      ELSE IF (isubstep .LT. nintervals) THEN
      
        ! We are somewhere 'in the middle'.
        !
        ! Dual RHS comes from rtempVector3. The primal from the
        ! isubstep-1'th RHS.
        !
        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
        
        CALL lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(1),rtempVector2%RvectorBlock(1),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(1))                                        
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(2),rtempVector2%RvectorBlock(2),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(2))                                        
        ! Pressure is fully implicit
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(3),rtempVector2%RvectorBlock(3),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(3))

        CALL lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(4),rtempVector3%RvectorBlock(4),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(4))
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(5),rtempVector3%RvectorBlock(5),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(5))
        ! Pressure is fully implicit
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(6),rtempVector3%RvectorBlock(6),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(6))
        
        IF (isubstep .LT. rspaceTimeDiscr%rtimeDiscr%nintervals-1) THEN
          ! Shift the RHS vectors and generate the RHS for the next time step.
          ! (Yes, I know, this could probably be solved more elegant without copying anything
          ! using a ring buffer ^^)
          CALL lsysbl_copyVector(rtempVector2,rtempVector1)
          CALL lsysbl_copyVector(rtempVector3,rtempVector2)
          CALL trhsevl_assembleSpatialRHS (rproblem,isubstep+2,&
              (isubstep+2)*rspaceTimeDiscr%rtimeDiscr%dtstep,rtempVector3)
        END IF
        
      ELSE
      
        ! We are 'at the end'.
        !
        ! Dual RHS comes from rtempVector3. The primal from the
        ! isubstep-1'th RHS and rtempVector3.
        !
        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
        ! dual RHS(0)   = DUALRHS(0)
      
        CALL lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(1),rtempVector3%RvectorBlock(1),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(1))                                        
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(2),rtempVector3%RvectorBlock(2),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(2))                                        
        ! Pressure is fully implicit
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(3),rtempVector3%RvectorBlock(3),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(3))

        !CALL generateRHS (rproblem,isubstep+1,rspaceTimeDiscr%niterations,&
        !    rtempVector3, .TRUE., .FALSE.)

        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(6),rtempVectorRHS%RvectorBlock(6))

        ! Multiply the last RHS of the dual equation -z by gamma, that's it.
        CALL lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),rspaceTimeDiscr%dgammaC)
        CALL lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),rspaceTimeDiscr%dgammaC)

      END IF

      ! Implement the boundary conditions into the RHS vector        
      IF (bimplementBC) THEN
        CALL tbc_implementSpatialBCtoRHS (rproblem,isubstep,&
            isubstep*dtstep, rtempVectorRHS)
      END IF
      
      ! Save the RHS.
      CALL sptivec_setTimestepData(rb, 1+isubstep, rtempVectorRHS)
      
    END DO
    
    ! Release the temp vectors.
    CALL lsysbl_releaseVector (rtempVectorRHS)
    
    CALL lsysbl_releaseVector (rtempVector3)
    CALL lsysbl_releaseVector (rtempVector2)
    CALL lsysbl_releaseVector (rtempVector1)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE trhsevl_assembledG0RHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)

!<description>
  ! Assembles the space-time RHS vector rb according to the dG(0)-scheme.
  !
  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
  ! this routine!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time system.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
  
  ! Whether to implement boundary conditions into the RHS or not.
  LOGICAL, INTENT(IN) :: bimplementBC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
 
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubstep,nintervals
    REAL(DP) :: dtstep
    
    ! Temporary vector
    TYPE(t_vectorBlock) :: rtempVector

    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd, p_Drhs

    ! Create a temp vector for the assembly
    CALL lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,rtempVector,.TRUE.)

    CALL lsysbl_getbase_double (rtempVector,p_Db)

    ! The assembly for dG0-RHS is rather easy.
    ! The i'th time DOF belongs to the i'th time interval Ti and
    ! has the following form:
    !  ___
    !  f_i  =  1/|Ti| int_Ti f(.,t) dt  ~  f(.,T_i(midpoint))
    !
    ! by the midpoint rule in time! So we just make a loop
    ! over the timesteps and calculate the f()-values in the
    ! midpoints of the time intervals!
    
    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep
    nintervals = rspaceTimeDiscr%rtimeDiscr%nintervals

    DO isubstep = 0,nintervals-1
      ! Assemble at the midpoint of the time interval
      CALL trhsevl_assembleSpatialRHS (rproblem,isubstep,&
        (REAL(isubstep,DP)+0.5_DP)*dtstep,rtempVector)
        
      CALL sptivec_setTimestepData(rb, 1+isubstep, rtempVector)
    END DO
      
    ! Release the temp vector.
    CALL lsysbl_releaseVector (rtempVector)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE trhsevl_assembleSpatialRHS (rproblem, isubstep, dtime, rrhs)
  
!<description>
  ! Generate the RHS vector at the time dtime. isubstep may specify the
  ! timestep where to generate the RHS.
!</description>
 
!<input>
  ! Time where the RHS should be generated.
  ! Must not necessarily coincide with the start/end time of the timestep.
  REAL(DP), INTENT(IN) :: dtime
    
  ! Number of the substep where to generate the RHS vector.
  INTEGER, INTENT(IN) :: isubstep
!</input>  

!<inputoutput>  
  ! Problem structure.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! Destination vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>
  
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rrhs,p_Drhs)

    ! Set the time where we are at the moment
    !rproblem%rtimedependence%dtime = &
    !    rproblem%rtimedependence%dtimeInit + isubstep*dtstep
    !rproblem%rtimedependence%itimestep = isubstep

    rproblem%rtimedependence%dtime = dtime
    rproblem%rtimedependence%itimestep = isubstep

    ! Get a pointer to the RHS on the finest level as well as to the
    ! block discretisation structure:
    p_rdiscretisation => rrhs%p_rblockDiscretisation
    
    ! The vector structure is already prepared, but the entries are missing.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the 
    ! first block.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Note that the vector is unsorted after calling this routine!
    !
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Discretise the X-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(2),rlinform,.TRUE.,&
              rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    CALL lsyssc_clearVector(rrhs%RvectorBlock(3))
    
    ! The RHS terms for the dual equation are calculated similarly using
    ! the desired 'target' flow field.
    !
    ! Discretise the X-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              rrhs%RvectorBlock(4),coeff_TARGET_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(2),rlinform,.TRUE.,&
              rrhs%RvectorBlock(5),coeff_TARGET_y,&
              rproblem%rcollection)
      
    ! Depending on the formulation, to get a reference dual velocity,
    ! it might be necessary to switch the sign of the target velocity field 
    ! because the RHS of the dual equation is '-z'!
    ! Remember that it this case the signs of the mass matrices that couple
    ! primal and dual velocity must be changed, too!
    
    IF (rproblem%roptcontrol%ispaceTimeFormulation .EQ. 0) THEN
      CALL lsyssc_scaleVector (rrhs%RvectorBlock(4),-1.0_DP)
      CALL lsyssc_scaleVector (rrhs%RvectorBlock(5),-1.0_DP)
    END IF

    ! Dual pressure RHS is =0.
    CALL lsyssc_clearVector(rrhs%RvectorBlock(6))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
  
  END SUBROUTINE
    
END MODULE
