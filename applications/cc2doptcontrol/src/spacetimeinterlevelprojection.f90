!##############################################################################
!# ****************************************************************************
!# <name> spacetimeinterevelprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises the interlevel projection (prolongation and 
!# restriction) of space-time coupled solution vectors. The t_sptiProjection
!# structure configures the interlevel projection. For the actual projection,
!# the following routines can be found in this module:
!#
!# 1.) sptipr_initProjection
!#     -> Initialises a space-time interlevel projection structure
!#
!# 2.) sptipr_doneProjection
!#     -> Cleans up a space-time interlevel projection structure
!#
!# 3.) sptipr_performProlongation
!#     -> Prolongation of a space-time coupled solution vector to a higher 
!#        level
!#
!# 4.) sptipr_performRestriction
!#     -> Retriction of a space-time coupled defect vector to a lower level
!#
!# 5.) sptipr_performInterpolation
!#     -> Interpolation of a space-time coupled solution vector to a lower
!#        level
!#
!# </purpose>
!##############################################################################

MODULE spacetimeinterlevelprojection

  USE fsystem
  USE genoutput
  USE storage
  USE linearsystemscalar
  USE linearsystemblock
  
  USE spacetimevectors
  USE multilevelprojection

  USE cc2dmediumm2spacetimediscret
  USE cc2dmediumm2postprocessing

  IMPLICIT NONE

!<types>

!<typeblock>

  ! Space-time projection structure.
  TYPE t_sptiProjection
  
    ! Order of projection in time.
    ! =1: first order (linear), =2: second order (quadratic)
    INTEGER :: itimeOrder = 0
    
    ! Whether restriction in space should be performed.
    ! =FALSE: Only perform prolongation/restriction in time.
    ! =TRUE(default) : Perform prolongation/restriction simultaneously in space/time.
    LOGICAL :: bspaceTimeSimultaneously = .TRUE.
  
    ! Interlevel-projection-structure that describes the spatial prolongation/
    ! restriction.
    TYPE(t_interlevelProjectionBlock) :: rspatialProjection
  
  END TYPE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptipr_initProjection (rprojection,rdiscretisation,bspaceProjection)
  
!<description>
  ! Initialises a space-time interlevel projection structure rprojection based
  ! on the discretisation structure rprojection. By default, the time 
  ! interpolation is configured to be 1st order, and the prolongation/restriction
  ! is done simultaneously in space/time.
  !
  ! The resulting space-time interlevel projection structure can be used to
  ! perform prolongation/restriction in space/time. It defines how to do
  ! prolongation/restriction between level i-1 and level i.
!</description>

!<input>
  ! A block discretisation structure specifying the spatial discretisation.
  TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscretisation
  
  ! Whether prolongation/restriction should be done also in space, not only
  ! on time.
  ! =TRUE: space and time restriction
  ! =FALSE: only time restriction
  LOGICAL, INTENT(IN) :: bspaceProjection
!</input>

!<output>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  TYPE(t_sptiProjection), INTENT(OUT) :: rprojection
!</output>

!</subroutine>

    ! Set itimeOrder=1 -> structure initialised.
    rprojection%itimeOrder = 1

    ! Intialise the spatial interlevel projection structure with the
    ! standard method.
    CALL mlprj_initProjectionDiscr (rprojection%rspatialProjection,rdiscretisation)
    
    rprojection%bspaceTimeSimultaneously = bspaceProjection
    
    ! The other variables are left in their default initialisation as
    ! set by INTENT(OUT).
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptipr_doneProjection (rprojection)
  
!<description>
  ! Releases memory allocated in sptipr_initProjection and cleans up rprojection.
!</description>

!<inputoutput>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time. 
  ! The structure is cleaned up.
  TYPE(t_sptiProjection), INTENT(INOUT) :: rprojection
!</output>

!</inputoutput>

    ! Clean up the spatial projection structure
    CALL mlprj_doneProjection (rprojection%rspatialProjection)
    
    ! Set itimeOrder=0 -> structure not initialised anymore.
    rprojection%itimeOrder = 0

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptipr_performProlongation (rprojection,rcoarseVector, &
                                         rfineVector,rtempVecCoarse,rtempVecFine,&
                                         rdiscrFine,rdiscrCoarse,rproblem)
  
!<description>
  ! Performs a prolongation for a given space/time vector (i.e. a projection
  ! in the primal space where the solution lives). The vector
  ! rcoarseVector on a coarser space/time mesh is projected to the vector
  ! rfineVector on a finer space/time mesh. 
  ! rprojection configures how the transfer is performed.
  ! This projection structure rprojection must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  TYPE(t_sptiProjection), INTENT(IN) :: rprojection

  ! Coarse grid vector
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rcoarseVector
!</input>

  ! A space-time discretisation structure defining the discretisation of
  ! rfineVector.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rdiscrFine,rdiscrCoarse

  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

!<inputoutput>
  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the coarse grid.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecCoarse

  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the fine grid.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecFine
!</inputoutput>

!<output>
  ! Fine grid vector
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rfineVector
!</output>

!</subroutine>

    ! local variables
    INTEGER :: istep
    TYPE(t_vectorBlock) :: rx1,rx3
    TYPE(t_vectorScalar) :: rtempVecFineScalar

    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx1,p_Dx3
    REAL(DP), DIMENSION(:), POINTER :: p_DtempVecCoarse,p_DtempVecFine
    REAL(DP), DIMENSION(:), POINTER :: p_DtempVecFineSca

    ! Linear prolongation:
    !
    ! We need two more temp vectors:
    !
    ! One for the current timestep
    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.FALSE.)
    
    ! And one for the next timestep
    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.FALSE.)
    
    ! We need a scalar representation of the temp vector
    CALL lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rx1,p_Dx1)
    CALL lsysbl_getbase_double (rx3,p_Dx3)
    CALL lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
    CALL lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
    CALL lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)

    ! Load timestep 0 into the temp vector and interpolate to the current level.
    ! Put the result to rx3.
    IF (rprojection%bspaceTimeSimultaneously) THEN
      ! Space + time
      CALL sptivec_getTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
      CALL mlprj_performProlongation (rprojection%rspatialProjection,rtempVecCoarse, &
                                      rx3,rtempVecFineScalar)
    ELSE
      ! Only time
      CALL sptivec_getTimestepData (rcoarseVector, 1+0, rx3)
    END IF
    
    ! Save that vector as initial vector on the fine grid.
    CALL sptivec_setTimestepData (rfineVector, 1+0, rx3)
    
    IF (rcoarseVector%NEQtime .NE. rfineVector%NEQtime) THEN
      
      ! Prolongation in time.
      !
      ! Loop through the time steps.
      DO istep = 1,rcoarseVector%NEQtime-1
    
        ! rx3 was the vector from the last timestep. Shift it to rx1 and load
        ! the vector for the new timestep in rx3
        CALL lsysbl_copyVector (rx3,rx1)

        IF (rprojection%bspaceTimeSimultaneously) THEN
          ! Space + time
          CALL sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
          CALL mlprj_performProlongation (rprojection%rspatialProjection,rtempVecCoarse, &
                                          rx3,rtempVecFineScalar)
        ELSE
          ! Only time
          CALL sptivec_getTimestepData (rcoarseVector, 1+istep, rx3)
        END IF
        
        ! Save that vector as new vector on the fine grid.
        CALL sptivec_setTimestepData (rfineVector, 1+2*istep, rx3)
        
        ! In the temp vector, create the interpolation between rx1 and rx3.
        ! THat's the prolongated vector.
        CALL lsysbl_copyVector (rx1,rtempVecFine)
        CALL lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP,0.5_DP)

        ! Save that vector as new vector on the fine grid.
        CALL sptivec_setTimestepData (rfineVector, 1+2*istep-1, rtempVecFine)
      
      END DO
      
    ELSE
    
      IF (rprojection%bspaceTimeSimultaneously) THEN
      
        ! Prolongation only in space. 
    
        ! Loop through the time steps.
        DO istep = 1,rcoarseVector%NEQtime-1
    
          ! Space prolongation
          CALL sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
          CALL mlprj_performProlongation (rprojection%rspatialProjection,rtempVecCoarse, &
                                          rx3,rtempVecFineScalar)

          ! Save that vector as new vector on the fine grid.
          CALL sptivec_setTimestepData (rfineVector, 1+istep, rx3)

        END DO

      END IF
    
    END IF
    
    ! Release the temp vectors
    CALL lsyssc_releaseVector (rtempVecFineScalar)
    CALL lsysbl_releaseVector (rx3)
    CALL lsysbl_releaseVector (rx1)
    
    ! DEBUG!!!
    !CALL c2d2_postprocSpaceTimeGMV (rproblem,rdiscrFine,rfineVector,'gmv/fine.gmv')
    !CALL c2d2_postprocSpaceTimeGMV (rproblem,rdiscrCoarse,rcoarseVector,'gmv/coarse.gmv')
    
!    ! Constant prolongation
!    !
!    ! We need two more temp vectors:
!    !
!    ! One for the current timestep
!    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.FALSE.)
!    
!    ! And one for the next timestep
!    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.FALSE.)
!    
!    ! We need a scalar representation of the temp vector
!    CALL lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
!    
!    ! DEBUG!!!
!    CALL lsysbl_getbase_double (rx1,p_Dx1)
!    CALL lsysbl_getbase_double (rx3,p_Dx3)
!    CALL lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!    CALL lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!    CALL lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)
!
!    ! Load timestep 0 into the temp vector and interpolate to the current level.
!    ! Put the result to rx3.
!    IF (rprojection%bspaceTimeSimultaneously) THEN
!      ! Space + time
!      CALL sptivec_getTimestepData (rcoarseVector, 0, rtempVecCoarse)
!      CALL mlprj_performProlongation (rprojection%rspatialProjection,rtempVecCoarse, &
!                                      rx3,rtempVecFineScalar)
!    ELSE
!      ! Only time
!      CALL sptivec_getTimestepData (rcoarseVector, 0, rx3)
!    END IF
!    
!    ! Save that vector as initial vector on the fine grid.
!    CALL sptivec_setTimestepData (rfineVector, 0, rx3)
!    
!    IF (rcoarseVector%ntimesteps .NE. rfineVector%ntimesteps) THEN
!      
!      ! Prolongation in time.
!      !
!      ! Loop through the time steps.
!      DO istep = 1,rcoarseVector%ntimesteps
!    
!        ! rx3 was the vector from the last timestep. Shift it to rx1 and load
!        ! the vector for the new timestep in rx3
!        CALL lsysbl_copyVector (rx3,rx1)
!
!        IF (rprojection%bspaceTimeSimultaneously) THEN
!          ! Space + time
!          CALL sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!          CALL mlprj_performProlongation (rprojection%rspatialProjection,rtempVecCoarse, &
!                                          rx3,rtempVecFineScalar)
!        ELSE
!          ! Only time
!          CALL sptivec_getTimestepData (rcoarseVector, istep, rx3)
!        END IF
!        
!        ! Save that vector as new vector on the fine grid.
!        CALL sptivec_setTimestepData (rfineVector, 2*istep, rx3)
!        
!        ! In the temp vector, create the prolongation. Use the constant prolongation.
!        ! Primal vector of rx1 -> fine vector <- Dual vector of rx3.
!        ! That's the prolongated vector.
!        CALL lsysbl_copyVector (rx1,rtempVecFine)
!        CALL lsyssc_copyVector (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4))
!        CALL lsyssc_copyVector (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5))
!        CALL lsyssc_copyVector (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6))
!
!        ! Save that vector as new vector on the fine grid.
!        CALL sptivec_setTimestepData (rfineVector, 2*istep-1, rtempVecFine)
!      
!      END DO
!      
!    ELSE
!    
!      IF (rprojection%bspaceTimeSimultaneously) THEN
!      
!        ! Prolongation only in space. 
!    
!        ! Loop through the time steps.
!        DO istep = 1,rcoarseVector%ntimesteps
!    
!          ! Space prolongation
!          CALL sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!          CALL mlprj_performProlongation (rprojection%rspatialProjection,rtempVecCoarse, &
!                                          rx3,rtempVecFineScalar)
!
!          ! Save that vector as new vector on the fine grid.
!          CALL sptivec_setTimestepData (rfineVector, istep, rx3)
!
!        END DO
!
!      END IF
!    
!    END IF
!    
!    ! Release the temp vectors
!    CALL lsyssc_releaseVector (rtempVecFineScalar)
!    CALL lsysbl_releaseVector (rx3)
!    CALL lsysbl_releaseVector (rx1)
!
  END SUBROUTINE
     
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptipr_performRestriction (rprojection,rcoarseVector, &
                                         rfineVector,rtempVecCoarse,rtempVecFine)
  
!<description>
  ! Performs a restriction for a given space/time vector (i.e. a projection
  ! in the dual space where the RHS vector lives).The vector
  ! rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser space/time mesh. 
  ! rprojection configures how the transfer is performed.
  ! This projection structure rprojection must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  TYPE(t_sptiProjection), INTENT(IN) :: rprojection

  ! Fine grid vector
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rfineVector
!</input>

!<inputoutput>
  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the coarse grid.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecCoarse

  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the fine grid.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecFine
!</inputoutput>

!<output>
  ! Coarse grid vector
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rcoarseVector
!</output>

!</subroutine>

    ! local variables
    INTEGER :: istep
    TYPE(t_vectorBlock) :: rx1,rx3
    TYPE(t_vectorScalar) :: rx1scalar
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx1,p_Dx3
    REAL(DP), DIMENSION(:), POINTER :: p_DtempVecCoarse,p_DtempVecFine

    ! Linear restriction:
    !
    ! Prolongation 'distributes' information from the coarse grid nodes
    ! to the fine grid nodes as follows:
    !
    ! Timestep:    n-1         n          n+1        (fine grid)
    !               x <--1/2-- X --1/2--> x
    !                         ^ \
    !                        /   \
    !                       +--1--+
    !
    ! Restriction is the 'adjoint' of the prolongation, so the corresponding
    ! matrix is the transposed. Written nodewise this means, that restriction
    ! has to 'collect' the distributed information with the same weights:
    !
    ! Timestep:    n-1         n          n+1        (fine grid)
    !               x --1/2--> X <--1/2-- x
    !                         ^ \
    !                        /   \
    !                       +--1--+
    !
    ! But because this is a restriction of a Finite-Difference RHS vector,
    ! the RHS must be divided by h/(h/2)=2 !

    ! We need two more temp vectors:
    !
    ! One for the current timestep
    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.FALSE.)
    
    ! And one for the next timestep
    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.FALSE.)
    
    ! Create a scalar representation of rx1 for auxiliary reasons
    CALL lsysbl_createScalarFromVec (rx1,rx1Scalar)
   
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rx1,p_Dx1)
    CALL lsysbl_getbase_double (rx3,p_Dx3)
    CALL lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
    CALL lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
       
    IF (rcoarseVector%NEQtime .NE. rfineVector%NEQtime) THEN       
    
      ! Load timestep 0,1 (fine grid) into the temp vectors.
      CALL sptivec_getTimestepData (rfineVector, 1+0, rtempVecFine)
      CALL sptivec_getTimestepData (rfineVector, 1+1, rx3)
      
      ! Perform restriction for the first time step.
      !                          X <--1/2-- x
      !                         ^ \
      !                        /   \
      !                       +--1--+
      CALL lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
      
      ! Probably restrict the vector in space to the lower level.
      ! Save the result.
      IF (rprojection%bspaceTimeSimultaneously) THEN
        ! Space + time
        CALL mlprj_performRestriction (rprojection%rspatialProjection,rtempVecCoarse, &
                                      rtempVecFine,rx1Scalar)
        CALL sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
      ELSE
        ! Only time
        CALL sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecFine)
      END IF
      
      ! Loop through the time steps -- on the coarse grid!
      DO istep = 1,rcoarseVector%NEQtime-1-1
      
        ! rx3 was the 'right inner fine grid node x' in the above picture.
        ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
        ! and load the new 'right inner fine grid node x' to rx1
        ! as well as the 'inner fine grid node X' to rtempVecFine.
        CALL lsysbl_copyVector (rx3,rx1)
        
        CALL sptivec_getTimestepData (rfineVector,1+2*istep, rtempVecFine)
        CALL sptivec_getTimestepData (rfineVector,1+2*istep+1, rx3)
        
        ! In rtempVecFine, create the restriction of the fine grid vectors
        ! rx1 and rx3 according to
        !
        ! Timestep:    n-1          n          n+1        (fine grid)
        !              rx3 --1/2--> X <--1/2-- rx3
        !             (old)        ^ \
        !                         /   \
        !                        +--1--+
        
        CALL lsysbl_vectorLinearComb (rx1,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
        CALL lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP)
        
        ! Probably restrict the vector in space to the lower level.
        ! Save the result.
        IF (rprojection%bspaceTimeSimultaneously) THEN
          ! Space + time
          CALL mlprj_performRestriction (rprojection%rspatialProjection,rtempVecCoarse, &
                                        rtempVecFine,rx1Scalar)
          CALL sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
        ELSE
          ! Only time
          CALL sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
        END IF
      
      END DO

      ! Last time step.

      ! rx3 is now the new 'left inner fine grid node x'.
      ! Load the 'inner fine grid node X' to rtempVecFine.
      CALL sptivec_getTimestepData (rfineVector,1+2*(rcoarseVector%NEQtime-1), &
        rtempVecFine)
      
      ! In rtempVecFine, create the restriction of the fine grid vectors
      ! rx1 and rx3 according to
      !
      ! Timestep:    n-1          n       (fine grid)
      !              rx3 --1/2--> X             
      !                          ^ \
      !                         /   \
      !                        +--1--+
      
      CALL lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
      
      ! Probably restrict the vector in space to the lower level.
      ! Save the result.
      IF (rprojection%bspaceTimeSimultaneously) THEN
        ! Space + time
        CALL mlprj_performRestriction (&
            rprojection%rspatialProjection,rtempVecCoarse, rtempVecFine,rx1Scalar)
        CALL sptivec_setTimestepData (rcoarseVector, &
            rcoarseVector%NEQtime, rtempVecCoarse)
      ELSE
        ! Only time
        CALL sptivec_setTimestepData (rcoarseVector, &
            rcoarseVector%NEQtime, rtempVecFine)
      END IF
      
    ELSE
    
      IF (rprojection%bspaceTimeSimultaneously) THEN
      
        ! Restriction only in space.
        !
        ! Loop through the time steps
        DO istep = 0,rcoarseVector%NEQtime-1
        
          ! Load the data
          CALL sptivec_getTimestepData (rfineVector,1+istep, rtempVecFine)

          ! Process the restriction and save the data to the coarse vector.
          CALL mlprj_performRestriction (&
              rprojection%rspatialProjection,rtempVecCoarse, rtempVecFine,rx1Scalar)
              
          CALL sptivec_setTimestepData (rcoarseVector,1+istep, rtempVecCoarse)
        
        END DO
      
      END IF
    
    END IF

    ! Release the temp vectors
    CALL lsyssc_releaseVector (rx1Scalar)
    CALL lsysbl_releaseVector (rx3)
    CALL lsysbl_releaseVector (rx1)

!    ! Constant restriction:
!    !
!    ! Prolongation 'distributes' information from the coarse grid nodes
!    ! to the fine grid nodes as follows:
!    !
!    ! Timestep:    n-1         n          n+1        (fine grid)
!    !               x <--1/2-- X --1/2--> x
!    !                         ^ \
!    !                        /   \
!    !                       +--1--+
!    !
!    ! Restriction is the 'adjoint' of the prolongation, so the corresponding
!    ! matrix is the transposed. Written nodewise this means, that restriction
!    ! has to 'collect' the distributed information with the same weights:
!    !
!    ! Timestep:    n-1         n          n+1        (fine grid)
!    !               x --1/2--> X <--1/2-- x
!    !                         ^ \
!    !                        /   \
!    !                       +--1--+
!    !
!    ! But because this is a restriction of a Finite-Difference RHS vector,
!    ! the RHS must be divided by h/(h/2)=2 !
!
!    ! We need two more temp vectors:
!    !
!    ! One for the current timestep
!    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.FALSE.)
!    
!    ! And one for the next timestep
!    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.FALSE.)
!    
!    ! Create a scalar representation of rx1 for auxiliary reasons
!    CALL lsysbl_createScalarFromVec (rx1,rx1Scalar)
!   
!    ! DEBUG!!!
!    CALL lsysbl_getbase_double (rx1,p_Dx1)
!    CALL lsysbl_getbase_double (rx3,p_Dx3)
!    CALL lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!    CALL lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!       
!    IF (rcoarseVector%ntimesteps .NE. rfineVector%ntimesteps) THEN       
!    
!      ! Load timestep 0,1 (fine grid) into the temp vectors.
!      CALL sptivec_getTimestepData (rfineVector, 0, rtempVecFine)
!      CALL sptivec_getTimestepData (rfineVector, 1, rx3)
!      
!      ! Perform restriction for the first time step.
!      !                          X <--1/2-- x
!      !                         ^ \
!      !                        /   \
!      !                       +--1--+
!      CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4),&
!          0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!      CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5),&
!          0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!      CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6),&
!          0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!      
!      ! Probably restrict the vector in space to the lower level.
!      ! Save the result.
!      IF (rprojection%bspaceTimeSimultaneously) THEN
!        ! Space + time
!        CALL mlprj_performRestriction (rprojection%rspatialProjection,rtempVecCoarse, &
!                                      rtempVecFine,rx1Scalar)
!        CALL sptivec_setTimestepData (rcoarseVector, 0, rtempVecCoarse)
!      ELSE
!        ! Only time
!        CALL sptivec_setTimestepData (rcoarseVector, 0, rtempVecFine)
!      END IF
!      
!      ! Loop through the time steps -- on the coarse grid!
!      DO istep = 1,rcoarseVector%ntimesteps-1
!      
!        ! rx3 was the 'right inner fine grid node x' in the above picture.
!        ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
!        ! and load the new 'right inner fine grid node x' to rx1
!        ! as well as the 'inner fine grid node X' to rtempVecFine.
!        CALL lsysbl_copyVector (rx3,rx1)
!        
!        CALL sptivec_getTimestepData (rfineVector,2*istep, rtempVecFine)
!        CALL sptivec_getTimestepData (rfineVector,2*istep+1, rx3)
!        
!        ! In rtempVecFine, create the restriction of the fine grid vectors
!        ! rx1 and rx3 according to
!        !
!        ! Timestep:    n-1          n          n+1        (fine grid)
!        !              rx3 --1/2--> X <--1/2-- rx3
!        !             (old)        ^ \
!        !             =rx1        /   \
!        !                        +--1--+
!        
!        CALL lsyssc_vectorLinearComb (rx1%RvectorBlock(1),rtempVecFine%RvectorBlock(1),&
!            0.5_DP,1.0_DP/2.0_DP)
!        CALL lsyssc_vectorLinearComb (rx1%RvectorBlock(2),rtempVecFine%RvectorBlock(2),&
!            0.5_DP,1.0_DP/2.0_DP)
!        CALL lsyssc_vectorLinearComb (rx1%RvectorBlock(3),rtempVecFine%RvectorBlock(3),&
!            0.5_DP,1.0_DP/2.0_DP)
!
!        CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4),&
!            0.5_DP,1.0_DP/2.0_DP)
!        CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5),&
!            0.5_DP,1.0_DP/2.0_DP)
!        CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6),&
!            0.5_DP,1.0_DP/2.0_DP)
!
!        ! Probably restrict the vector in space to the lower level.
!        ! Save the result.
!        IF (rprojection%bspaceTimeSimultaneously) THEN
!          ! Space + time
!          CALL mlprj_performRestriction (rprojection%rspatialProjection,rtempVecCoarse, &
!                                        rtempVecFine,rx1Scalar)
!          CALL sptivec_setTimestepData (rcoarseVector, istep, rtempVecCoarse)
!        ELSE
!          ! Only time
!          CALL sptivec_setTimestepData (rcoarseVector, istep, rtempVecFine)
!        END IF
!      
!      END DO
!
!      ! Last time step.
!
!      ! rx3 is now the new 'left inner fine grid node x'.
!      ! Load the 'inner fine grid node X' to rtempVecFine.
!      CALL sptivec_getTimestepData (rfineVector,2*rcoarseVector%ntimesteps, rtempVecFine)
!      
!      ! In rtempVecFine, create the restriction of the fine grid vectors
!      ! rx1 and rx3 according to
!      !
!      ! Timestep:    n-1          n       (fine grid)
!      !              rx3 --1/2--> X             
!      !                          ^ \
!      !                         /   \
!      !                        +--1--+
!      
!      CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(1),rtempVecFine%RvectorBlock(1),&
!          0.5_DP,1.0_DP/2.0_DP)
!      CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(2),rtempVecFine%RvectorBlock(2),&
!          0.5_DP,1.0_DP/2.0_DP)
!      CALL lsyssc_vectorLinearComb (rx3%RvectorBlock(3),rtempVecFine%RvectorBlock(3),&
!          0.5_DP,1.0_DP/2.0_DP)
!      
!      ! Probably restrict the vector in space to the lower level.
!      ! Save the result.
!      IF (rprojection%bspaceTimeSimultaneously) THEN
!        ! Space + time
!        CALL mlprj_performRestriction (&
!            rprojection%rspatialProjection,rtempVecCoarse, rtempVecFine,rx1Scalar)
!        CALL sptivec_setTimestepData (rcoarseVector, &
!            rcoarseVector%ntimesteps, rtempVecCoarse)
!      ELSE
!        ! Only time
!        CALL sptivec_setTimestepData (rcoarseVector, &
!            rcoarseVector%ntimesteps, rtempVecFine)
!      END IF
!      
!    ELSE
!    
!      IF (rprojection%bspaceTimeSimultaneously) THEN
!      
!        ! Restriction only in space.
!        !
!        ! Loop through the time steps
!        DO istep = 0,rcoarseVector%ntimesteps
!        
!          ! Load the data
!          CALL sptivec_getTimestepData (rfineVector,istep, rtempVecFine)
!
!          ! Process the restriction and save the data to the coarse vector.
!          CALL mlprj_performRestriction (&
!              rprojection%rspatialProjection,rtempVecCoarse, rtempVecFine,rx1Scalar)
!              
!          CALL sptivec_setTimestepData (rcoarseVector, istep, rtempVecCoarse)
!        
!        END DO
!      
!      END IF
!    
!    END IF
!
!    ! Release the temp vectors
!    CALL lsyssc_releaseVector (rx1Scalar)
!    CALL lsysbl_releaseVector (rx3)
!    CALL lsysbl_releaseVector (rx1)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptipr_performInterpolation (rprojection,rcoarseVector, &
                                          rfineVector,rtempVecCoarse,rtempVecFine)
  
!<description>
  ! Performs an interpolation for a given space/time vector to a lower level.
  ! The solution vector rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser space/time mesh. 
  ! rprojection configures how the transfer is performed.
  ! This projection structure rprojection must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  TYPE(t_sptiProjection), INTENT(IN) :: rprojection

  ! Fine grid vector
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rfineVector
!</input>

!<inputoutput>
  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the coarse grid.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecCoarse

  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the fine grid.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecFine
!</inputoutput>

!<output>
  ! Coarse grid vector
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rcoarseVector
!</output>

!</subroutine>

    ! local variables
    INTEGER :: istep
    TYPE(t_vectorBlock) :: rx1,rx3
    TYPE(t_vectorScalar) :: rx1scalar
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx1,p_Dx3
    REAL(DP), DIMENSION(:), POINTER :: p_DtempVecCoarse,p_DtempVecFine

    ! Prolongation 'distributes' information from the coarse grid nodes
    ! to the fine grid nodes as follows:
    !
    ! Timestep:    n-1         n          n+1        (fine grid)
    !               x <--1/2-- X --1/2--> x
    !                         ^ \
    !                        /   \
    !                       +--1--+
    !
    ! Interpolation is now taking the mean of adjacent nodes.
    ! Written nodewise this means:
    !
    ! Timestep:    n-1         n          n+1        (fine grid)
    !               x --1/4--> X <--1/4-- x
    !                         ^ \
    !                        /   \
    !                       +-1/2-+

    ! We need two more temp vectors:
    !
    ! One for the current timestep
    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.FALSE.)
    
    ! And one for the next timestep
    CALL lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.FALSE.)
    
    ! Create a scalar representation of rx1 for auxiliary reasons
    CALL lsysbl_createScalarFromVec (rx1,rx1Scalar)
   
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rx1,p_Dx1)
    CALL lsysbl_getbase_double (rx3,p_Dx3)
    CALL lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
    CALL lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
       
    IF (rcoarseVector%NEQtime .NE. rfineVector%NEQtime) THEN       
      ! Load timestep 0,1 (fine grid) into the temp vectors.
      CALL sptivec_getTimestepData (rfineVector, 1+0, rtempVecFine)
      CALL sptivec_getTimestepData (rfineVector, 1+1, rx3)
      
      ! Perform restriction for the first time step.
      !                          X <--1/4-- x
      !                         ^ \
      !                        /   \
      !                       +-3/4-+
      CALL lsysbl_vectorLinearComb (rx3,rtempVecFine,1._DP/4._DP,3._DP/4._DP)
      
      ! Probably restrict the vector in space to the lower level.
      ! Save the result.
      IF (rprojection%bspaceTimeSimultaneously) THEN
        ! Space + time
        CALL mlprj_performInterpolation (rprojection%rspatialProjection,rtempVecCoarse, &
                                        rtempVecFine,rx1Scalar)
        CALL sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
      ELSE
        ! Only time
        CALL sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecFine)
      END IF
      
      ! Loop through the time steps -- on the coarse grid!
      DO istep = 1,rcoarseVector%NEQtime-1-1
      
        ! rx3 was the 'right inner fine grid node x' in the above picture.
        ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
        ! and load the new 'right inner fine grid node x' to rx1
        ! as well as the 'inner fine grid node X' to rtempVecFine.
        CALL lsysbl_copyVector (rx3,rx1)
        
        CALL sptivec_getTimestepData (rfineVector,1+2*istep, rtempVecFine)
        CALL sptivec_getTimestepData (rfineVector,1+2*istep+1, rx3)
        
        ! In rtempVecFine, create the restriction of the fine grid vectors
        ! rx1 and rx3 according to
        !
        ! Timestep:    n-1          n          n+1        (fine grid)
        !              rx3 --1/4--> X <--1/4-- rx3
        !             (old)        ^ \
        !                         /   \
        !                        +-1/2-+
        
        CALL lsysbl_vectorLinearComb (rx1,rtempVecFine,0.25_DP,0.5_DP)
        CALL lsysbl_vectorLinearComb (rx3,rtempVecFine,0.25_DP,1.0_DP)
        
        ! Probably restrict the vector in space to the lower level.
        ! Save the result.
        IF (rprojection%bspaceTimeSimultaneously) THEN
          ! Space + time
          CALL mlprj_performInterpolation (rprojection%rspatialProjection,rtempVecCoarse, &
                                          rtempVecFine,rx1Scalar)
          CALL sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
        ELSE
          ! Only time
          CALL sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
        END IF
      
      END DO

      ! Last time step.

      ! rx3 is now the new 'left inner fine grid node x'.
      ! Load the 'inner fine grid node X' to rtempVecFine.
      CALL sptivec_getTimestepData (rfineVector,1+2*(rcoarseVector%NEQtime-1), &
          rtempVecFine)
      
      ! In rtempVecFine, create the restriction of the fine grid vectors
      ! rx1 and rx3 according to
      !
      ! Timestep:    n-1          n       (fine grid)
      !              rx3 --1/4--> X             
      !                          ^ \
      !                         /   \
      !                        +-3/4-+
      
      CALL lsysbl_vectorLinearComb (rx3,rtempVecFine,1._DP/4._DP,3._DP/4._DP)
      
      ! Probably restrict the vector in space to the lower level.
      ! Save the result.
      IF (rprojection%bspaceTimeSimultaneously) THEN
        ! Space + time
        CALL mlprj_performInterpolation (&
            rprojection%rspatialProjection,rtempVecCoarse, rtempVecFine,rx1Scalar)
        CALL sptivec_setTimestepData (rcoarseVector, &
            rcoarseVector%NEQtime, rtempVecCoarse)
      ELSE
        ! Only time
        CALL sptivec_setTimestepData (rcoarseVector, &
            rcoarseVector%NEQtime, rtempVecFine)
      END IF

    ELSE
    
      IF (rprojection%bspaceTimeSimultaneously) THEN
      
        ! Interpolation only in space.
      
        !
        ! Loop through the time steps
        DO istep = 0,rcoarseVector%NEQtime-1
        
          ! Load the data
          CALL sptivec_getTimestepData (rfineVector,1+istep, rtempVecFine)

          ! Process the restriction and save the data to the coarse vector.
          CALL mlprj_performInterpolation (&
              rprojection%rspatialProjection,rtempVecCoarse, rtempVecFine,rx1Scalar)
              
          CALL sptivec_setTimestepData (rcoarseVector,1+istep, rtempVecCoarse)
        
        END DO

      END IF
    
    END IF

    ! Release the temp vectors
    CALL lsyssc_releaseVector (rx1Scalar)
    CALL lsysbl_releaseVector (rx3)
    CALL lsysbl_releaseVector (rx1)

  END SUBROUTINE

END MODULE
