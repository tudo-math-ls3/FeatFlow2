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

module timerhsevaluation


  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use paramlist
  use timestepping
  
  use collection
  use convection
    
  use basicstructures
  use user_callback

  use spacepreconditioner
  use spacepreconditionerinit
  use timeanalysis
  use spatialbc
  use spacediscretisation
  use spacematvecassembly
  
  use timediscretisation
  use spacetimevectors
  use dofmapping
  use timeboundaryconditions
  
  use spacetimediscretisation

  implicit none

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)

!<description>
  ! Assembles the space-time RHS vector rb.
  !
  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
  ! this routine!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time system.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
  
  ! Whether to implement boundary conditions into the RHS or not.
  logical, intent(IN) :: bimplementBC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  ! If this is undefined, a new space-time vector is created.
  type(t_spacetimeVector), intent(INOUT) :: rb
!</inputoutput>

!</subroutine>

    if (rb%NEQtime .eq. 0) then
      ! Create a new vector if rb is undefined.
      call sptivec_initVectorDiscr (rb,rspaceTimeDiscr%rtimeDiscr,&
          rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
    end if

    ! What's the current time discretisation? Depending on that,
    ! we have to call the corresponding RHS calculation routine.
    
    select case (rspaceTimeDiscr%rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)
      call trhsevl_assembleThetaRHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)
    case (TDISCR_DG0)
      call trhsevl_assembledG0RHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)
    case DEFAULT
      call output_line ('Unsupported time discretisation', &
                        OU_CLASS_ERROR,OU_MODE_STD,'trhsevl_assembleRHS')
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine trhsevl_assembleThetaRHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)

!<description>
  ! Assembles the space-time RHS vector rb for a Theta-Scheme.
  !
  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
  ! this routine!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time system.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
  
  ! Whether to implement boundary conditions into the RHS or not.
  logical, intent(IN) :: bimplementBC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  type(t_spacetimeVector), intent(INOUT) :: rb
 
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: isubstep,nintervals
    real(DP) :: dtheta,dtstep
    
    ! Temporary vectors
    type(t_vectorBlock) :: rtempVector1,rtempVector2,rtempVector3

    ! A temporary vector for the creation of the RHS.
    type(t_vectorBlock) :: rtempVectorRHS
    
    real(DP), dimension(:),pointer :: p_Dx, p_Db, p_Dd, p_Drhs

    ! Create temp vectors for the assembly
    call lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector1,.true.)
    call lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector2,.true.)
    call lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector3,.true.)

    ! Theta-scheme identifier.
    ! =1: impliciz Euler.
    ! =0.5: Crank Nicolson
    dtheta = rspaceTimeDiscr%rtimeDiscr%dtheta
    
    ! Time step size, number of intervals
    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep
    nintervals = rspaceTimeDiscr%rtimeDiscr%nintervals
    
    ! ----------------------------------------------------------------------
    ! Generate the global RHS vector
    
    call lsysbl_getbase_double (rtempVector1,p_Dx)
    call lsysbl_getbase_double (rtempVector2,p_Db)
    call lsysbl_getbase_double (rtempVector3,p_Dd)

    ! Assemble 1st RHS vector in X temp vector.
    call trhsevl_assembleSpatialRHS (rproblem,0,0.0_DP,rtempVector1)
        
    ! Assemble the 2nd RHS vector in the RHS temp vector
    call trhsevl_assembleSpatialRHS (rproblem,1,dtstep,rtempVector2)

    ! Assemble the 3rd RHS vector in the defect temp vector
    if (rspaceTimeDiscr%rtimeDiscr%nintervals .ge. 2) then
      call trhsevl_assembleSpatialRHS (rproblem,2,2.0_DP*dtstep,rtempVector3)
    else
      call lsysbl_copyVector (rtempVector2,rtempVector3)
    end if
        
    ! Create a copy of the X temp vector (RHS0). That vector will be
    ! our destination vector for assembling the RHS in all timesteps.
    call lsysbl_copyVector (rtempVector1,rtempVectorRHS)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
    
    ! RHS 0,1,2 -> 1-2-3
    
    do isubstep = 0,nintervals
    
      if (isubstep .eq. 0) then
      
        ! Primal RHS comes from rtempVector1. The dual from the
        ! mean (depending on the timestep scheme) of 0 and the
        ! isubstep+1'th RHS in rtempVector2.
        !
        ! primal RHS(0) = PRIMALRHS(0)
        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)

        call lsyssc_copyVector (rtempVector1%RvectorBlock(1),rtempVectorRHS%RvectorBlock(1))
        call lsyssc_copyVector (rtempVector1%RvectorBlock(2),rtempVectorRHS%RvectorBlock(2))
        call lsyssc_copyVector (rtempVector1%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))

        call lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(4),rtempVector2%RvectorBlock(4),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(4))
        call lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(5),rtempVector2%RvectorBlock(5),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(5))
        ! Pressure is fully implicit
        call lsyssc_vectorLinearComb (&
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
            
      else if (isubstep .lt. nintervals) then
      
        ! We are somewhere 'in the middle'.
        !
        ! Dual RHS comes from rtempVector3. The primal from the
        ! isubstep-1'th RHS.
        !
        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
        
        call lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(1),rtempVector2%RvectorBlock(1),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(1))
        call lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(2),rtempVector2%RvectorBlock(2),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(2))
        ! Pressure is fully implicit
        call lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(3),rtempVector2%RvectorBlock(3),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(3))

        call lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(4),rtempVector3%RvectorBlock(4),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(4))
        call lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(5),rtempVector3%RvectorBlock(5),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(5))
        ! Pressure is fully implicit
        call lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(6),rtempVector3%RvectorBlock(6),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(6))
        
        if (isubstep .lt. rspaceTimeDiscr%rtimeDiscr%nintervals-1) then
          ! Shift the RHS vectors and generate the RHS for the next time step.
          ! (Yes, I know, this could probably be solved more elegant without copying anything
          ! using a ring buffer ^^)
          call lsysbl_copyVector(rtempVector2,rtempVector1)
          call lsysbl_copyVector(rtempVector3,rtempVector2)
          call trhsevl_assembleSpatialRHS (rproblem,isubstep+2,&
              (isubstep+2)*rspaceTimeDiscr%rtimeDiscr%dtstep,rtempVector3)
        end if
        
      else
      
        ! We are 'at the end'.
        !
        ! Dual RHS comes from rtempVector3. The primal from the
        ! isubstep-1'th RHS and rtempVector3.
        !
        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
        ! dual RHS(0)   = DUALRHS(0)
      
        call lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(1),rtempVector3%RvectorBlock(1),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(1))
        call lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(2),rtempVector3%RvectorBlock(2),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(2))
        ! Pressure is fully implicit
        call lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(3),rtempVector3%RvectorBlock(3),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(3))

        !CALL generateRHS (rproblem,isubstep+1,rspaceTimeDiscr%niterations,&
        !    rtempVector3, .TRUE., .FALSE.)
        call lsyssc_copyVector (rtempVector3%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
        call lsyssc_copyVector (rtempVector3%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
        call lsyssc_copyVector (rtempVector3%RvectorBlock(6),rtempVectorRHS%RvectorBlock(6))

        ! Multiply the last RHS of the dual equation -z by 1+gamma/dtstep, that's it.
        call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),&
            dtheta+rspaceTimeDiscr%dgammaC/dtstep)
        call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),&
            dtheta+rspaceTimeDiscr%dgammaC/dtstep)

      end if

      ! Implement the boundary conditions into the RHS vector
      if (bimplementBC) then
        call tbc_implementSpatialBCtoRHS (rproblem,&
            isubstep*dtstep, rtempVectorRHS)
      end if
      
      ! Save the RHS.
      call sptivec_setTimestepData(rb, 1+isubstep, rtempVectorRHS)
      
    end do
    
    ! Release the temp vectors.
    call lsysbl_releaseVector (rtempVectorRHS)
    
    call lsysbl_releaseVector (rtempVector3)
    call lsysbl_releaseVector (rtempVector2)
    call lsysbl_releaseVector (rtempVector1)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine trhsevl_assembledG0RHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)

!<description>
  ! Assembles the space-time RHS vector rb according to the dG(0)-scheme.
  !
  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
  ! this routine!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time system.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
  
  ! Whether to implement boundary conditions into the RHS or not.
  logical, intent(IN) :: bimplementBC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  type(t_spacetimeVector), intent(INOUT) :: rb
 
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: isubstep,nintervals
    real(DP) :: dtstep
    
    ! Temporary vector
    type(t_vectorBlock) :: rtempVector

    real(DP), dimension(:),pointer :: p_Dx, p_Db, p_Dd, p_Drhs

    ! Create a temp vector for the assembly
    call lsysbl_createVecBlockByDiscr (&
        rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.true.)

    call lsysbl_getbase_double (rtempVector,p_Db)

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

    do isubstep = 0,nintervals-1
      ! Assemble at the midpoint of the time interval
      call trhsevl_assembleSpatialRHS (rproblem,isubstep,&
        (real(isubstep,DP)+0.5_DP)*dtstep,rtempVector)
        
      call sptivec_setTimestepData(rb, 1+isubstep, rtempVector)
    end do
      
    ! Release the temp vector.
    call lsysbl_releaseVector (rtempVector)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine trhsevl_assembleSpatialRHS (rproblem, isubstep, dtime, rrhs)
  
!<description>
  ! Generate the RHS vector at the time dtime. isubstep may specify the
  ! timestep where to generate the RHS.
!</description>
 
!<input>
  ! Time where the RHS should be generated.
  ! Must not necessarily coincide with the start/end time of the timestep.
  real(DP), intent(IN) :: dtime
    
  ! Number of the substep where to generate the RHS vector.
  integer, intent(IN) :: isubstep
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Destination vector
  type(t_vectorBlock), intent(INOUT) :: rrhs
!</inputoutput>
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Drhs
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get a pointer to the RHS on the finest level as well as to the
    ! block discretisation structure:
    p_rdiscretisation => rrhs%p_rblockDiscr
    
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
    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

    ! Discretise the X-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
              rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
              rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhs%RvectorBlock(3))
    
    ! The RHS terms for the dual equation are calculated similarly using
    ! the desired 'target' flow field.
    !
    ! Discretise the X-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(4),rlinform,.true.,&
              rrhs%RvectorBlock(4),coeff_TARGET_x,&
              rproblem%rcollection)
    
    ! And the Y-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(5),rlinform,.true.,&
              rrhs%RvectorBlock(5),coeff_TARGET_y,&
              rproblem%rcollection)
      
    ! Depending on the formulation, to get a reference dual velocity,
    ! it might be necessary to switch the sign of the target velocity field
    ! because the RHS of the dual equation is '-z'!
    ! Remember that it this case the signs of the mass matrices that couple
    ! primal and dual velocity must be changed, too!
    if (rproblem%roptcontrol%ispaceTimeFormulation .eq. 0) then
      call lsyssc_scaleVector (rrhs%RvectorBlock(4),-1.0_DP)
      call lsyssc_scaleVector (rrhs%RvectorBlock(5),-1.0_DP)
    end if

    ! Dual pressure RHS is =0.
    call lsyssc_clearVector(rrhs%RvectorBlock(6))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
  
  end subroutine
    
end module
