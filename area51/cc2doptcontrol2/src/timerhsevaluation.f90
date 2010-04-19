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
  use genoutput
  use paramlist
  use basicgeometry
  use boundary
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  use derivatives
  use scalarpde
  use linearformevaluation
  
  use spatialdiscretisation
  use timediscretisation
  
  use cubature
  use filtersupport
  use collection
  
  use analyticsolution
  use meshhierarchy
  use fespacehierarchy
  use spacetimehierarchy
  
  use constantsoptc
  use structuresoptc
  use structuresspacetimelinsol
  use structuresoptcspacetimenlsol
  use structuresoptflow
  
  use spacediscretisation
  
  use spatialbcdef
  use initmatrices
  use postprocessing
  use user_callback
  
  use timescalehierarchy
  use spacematvecassembly
  use spacetimevectors
  use spacetimelinearsystem
  use forwardbackwardsimulation
  use spacetimeinterlevelprojection
  use nonlinearoneshotspacetimesolver
  use timeboundaryconditions

  implicit none
  
  private
  
  public :: trhsevl_assembleRHS

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine trhsevl_assembleRHS (rglobalData, rrhs, rrhsDiscrete, &
      roptimalControl, roptcBDC)

!<description>
  ! Discretises the RHS according to the space/time discretisation scheme
  ! in rrhsDiscrete. The result is summed up to rrhsDiscrete.
!</description>

!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
  
  ! Analytic solution defining the RHS of the equation.
  type(t_anSolution), intent(inout) :: rrhs
  
  ! Optimal control parameters.
  type(t_settings_optcontrol), intent(inout) :: roptimalControl
  
  ! OPTIONAL: Boundary conditions in the problem. If present, they are implemented
  ! to the RHS.
  type(t_optcBDC), intent(in), optional :: roptcBDC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  ! This vector must already be initialised. The calculated RHS is added to this.
  type(t_spacetimeVector), intent(inout) :: rrhsDiscrete
!</inputoutput>

!</subroutine>

    ! What's the current time discretisation? Depending on that,
    ! we have to call the corresponding RHS calculation routine.
    
    select case (rrhsDiscrete%p_rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)
      call trhsevl_assembleThetaRHS (rglobalData, rrhs, rrhsDiscrete, roptimalControl, roptcBDC)
    !case (TDISCR_DG0)
    !  call trhsevl_assembledG0RHS (rglobalData, rrhs, rrhsDiscrete, roptimalControl, roptcBDC)
    case default
      call output_line ('Unsupported time discretisation', &
                        OU_CLASS_ERROR,OU_MODE_STD,'trhsevl_assembleRHS')
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine trhsevl_assembleThetaRHS (rglobalData, rrhs, rrhsDiscrete, roptimalControl, roptcBDC)

!<description>
  ! Assembles the space-time RHS vector rrhsDiscrete for a Theta-Scheme. 
!</description>

!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
  
  ! Analytic solution defining the RHS of the equation.
  type(t_anSolution), intent(inout) :: rrhs
  
  ! Optimal control parameters.
  type(t_settings_optcontrol), intent(inout) :: roptimalControl

  ! OPTIONAL: Boundary conditions in the problem. If present, they are implemented
  ! to the RHS.
  type(t_optcBDC), intent(in), optional :: roptcBDC
!</input>

!<inputoutput>
  ! A space-time vector that receives the RHS.
  ! If this is undefined, a new space-time vector is created.
  type(t_spacetimeVector), intent(inout) :: rrhsDiscrete
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iiterate,nintervals
    real(DP) :: dtimePrimal,dtstep,dtimePrimal2,dtimeDual,dtimeDual2
    real(DP) :: dweightoldp,dweightnewp,dweightoldd,dweightnewd
    
    ! Temporary vectors
    type(t_vectorBlock) :: rtempVector1,rtempVector2,rtempVector3

    ! A temporary vector for the creation of the RHS.
    type(t_vectorBlock) :: rtempVectorRHS
    
    real(DP), dimension(:),pointer :: p_Dx, p_Db, p_Dd, p_Drhs

    ! Create temp vectors for the assembly
    call lsysbl_createVectorBlock (rrhsDiscrete%p_rspaceDiscr,rtempVector1,.true.)
    call lsysbl_createVectorBlock (rrhsDiscrete%p_rspaceDiscr,rtempVector2,.true.)
    call lsysbl_createVectorBlock (rrhsDiscrete%p_rspaceDiscr,rtempVector3,.true.)

    nintervals = rrhsDiscrete%p_rtimeDiscr%nintervals
    
    ! ----------------------------------------------------------------------
    ! Generate the global RHS vector
    
    call lsysbl_getbase_double (rtempVector1,p_Dx)
    call lsysbl_getbase_double (rtempVector2,p_Db)
    call lsysbl_getbase_double (rtempVector3,p_Dd)

    ! NOTE:
    ! For setting up the RHS, the time in the dual equation matches the time
    ! of the primal equation. For explicit Euler, this is the case anyway.
    ! For Crank-Nicolson / Theta scheme, the RHS must be evaluated in the
    ! time midpoint between two subsequent timesteps; but as the evaluation
    ! points of the dual equation are shifted by dt/2, the time midpoint
    ! in the dual equation coincides with the evaluation point in time of the
    ! primal equation!

    ! Assemble 1st RHS vector in X temp vector.
    call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,0,dtimePrimal,dtstep)
    dtimeDual = dtimePrimal
    call trhsevl_assembleSpatialRHS (rglobalData, rrhs, dtimePrimal, dtimeDual,&
        roptimalControl, rtempVector1)
        
    ! Assemble the 2nd RHS vector in the RHS temp vector
    call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,1,dtimePrimal,dtstep)
    dtimeDual = dtimePrimal
    call trhsevl_assembleSpatialRHS (rglobalData, rrhs, dtimePrimal, dtimeDual,&
        roptimalControl, rtempVector2)

    if (nintervals .gt. 1) then
      ! Assemble the RHS of the 'next' timestep into the 3rd temp vector if possible.
      call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,2,dtimePrimal2,dtstep)
      dtimeDual2 = dtimePrimal2
      call trhsevl_assembleSpatialRHS (rglobalData, rrhs, dtimePrimal2, dtimeDual2,&
          roptimalControl, rtempVector3)
    else
      ! 3rd temp vector is initialised with the rhs in rtempVector2.
      ! Dummy, as the 3rd subvector is usually not used in this case.
      call lsysbl_copyVector (rtempVector2,rtempVector3)
    end if

    ! Create a copy of the X temp vector (RHS0). That vector will be
    ! our destination vector for assembling the RHS in all timesteps.
    call lsysbl_copyVector (rtempVector1,rtempVectorRHS)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
    
    do iiterate = 1,nintervals+1
    
      ! Get the timestep weights for the current interval.
      call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal !- (1.0_DP-rrhsDiscrete%p_rtimeDiscr%dtheta)*dtstep
      call tdiscr_getTimestepWeights(rrhsDiscrete%p_rtimeDiscr,iiterate-1,dweightoldp,dweightnewp)
      call tdiscr_getTimestepWeights(rrhsDiscrete%p_rtimeDiscr,iiterate-1,dweightoldd,dweightnewd)
    
      if (iiterate .eq. 1) then
      
        ! RHS comes from rtempVector1. 
        !
        ! primal RHS(0) = PRIMALRHS(0)
        ! dual RHS(0)   = DUALRHS(0)

        call lsysbl_copyVector (rtempVector1,rtempVectorRHS)
        
        ! Scale the dual RHS according to the timestep scheme.
        call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),1.0_DP-rrhsDiscrete%p_rtimeDiscr%dtheta)
        call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),1.0_DP-rrhsDiscrete%p_rtimeDiscr%dtheta)
        
!        call lsyssc_copyVector (rtempVector1%RvectorBlock(1),rtempVectorRHS%RvectorBlock(1))
!        call lsyssc_copyVector (rtempVector1%RvectorBlock(2),rtempVectorRHS%RvectorBlock(2))
!        call lsyssc_copyVector (rtempVector1%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))
!
!        call lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(4),rtempVector2%RvectorBlock(4),&
!            dweightnewd,dweightoldd,&
!            rtempVectorRHS%RvectorBlock(4))
!        call lsyssc_vectorLinearComb (&                                                   
!            rtempVector1%RvectorBlock(5),rtempVector2%RvectorBlock(5),&
!            dweightnewd,dweightoldd,&
!            rtempVectorRHS%RvectorBlock(5))
!        ! Pressure is fully implicit
!        call lsyssc_vectorLinearComb (&                                                   
!            rtempVector1%RvectorBlock(6),rtempVector2%RvectorBlock(6),&
!            dweightnewd,dweightoldd,&
!            rtempVectorRHS%RvectorBlock(6))

        ! In the 0'th timestep, there is no RHS in the dual equation!
        ! That is because of the initial condition, which fixes the primal solution
        ! => dual solution has no influence on the primal one
        ! => setting up a dual RHS in not meaningful as the dual RHS cannot
        !    influence the primal solution
        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(4))
        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(5))
        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(6))
            
      else if (iiterate .lt. nintervals+1) then
      
        ! We are somewhere 'in the middle'.
        !
        ! Dual RHS comes from rtempVector3. The primal from the
        ! iiterate-1'th RHS.
        !
        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
        ! dual RHS(0)   = DUALRHS(0)
        
        call lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(1),rtempVector2%RvectorBlock(1),&
            dweightoldp,dweightnewp,&
            rtempVectorRHS%RvectorBlock(1))                                        
        call lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(2),rtempVector2%RvectorBlock(2),&
            dweightoldp,dweightnewp,&
            rtempVectorRHS%RvectorBlock(2))                                        
        ! Pressure is fully implicit
        call lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(3),rtempVector2%RvectorBlock(3),&
            dweightoldp,dweightnewp,&
            rtempVectorRHS%RvectorBlock(3))

        call lsyssc_copyVector (rtempVector2%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
        call lsyssc_copyVector (rtempVector2%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
        call lsyssc_copyVector (rtempVector2%RvectorBlock(6),rtempVectorRHS%RvectorBlock(6))
        
!        call lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(4),rtempVector3%RvectorBlock(4),&
!            dweightnewd,dweightoldd,&
!            rtempVectorRHS%RvectorBlock(4))
!        call lsyssc_vectorLinearComb (&                                                   
!            rtempVector2%RvectorBlock(5),rtempVector3%RvectorBlock(5),&
!            dweightnewd,dweightoldd,&
!            rtempVectorRHS%RvectorBlock(5))
!        ! Pressure is fully implicit
!        call lsyssc_vectorLinearComb (&                                                   
!            rtempVector2%RvectorBlock(6),rtempVector3%RvectorBlock(6),&
!            dweightnewd,dweightoldd,&
!            rtempVectorRHS%RvectorBlock(6))
        
        if (iiterate .lt. rrhsDiscrete%p_rtimeDiscr%nintervals) then
          ! Shift the RHS vectors and generate the RHS for the next time step.
          ! (Yes, I know, this could probably be solved more elegant without copying anything
          ! using a ring buffer ^^)
          call lsysbl_copyVector(rtempVector2,rtempVector1)
          call lsysbl_copyVector(rtempVector3,rtempVector2)

          ! Assemble the RHS of the 'next' timestep into the 3rd temp vector if possible.
          call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate+1,dtimePrimal2,dtstep)
          dtimeDual2 = dtimePrimal2
          call trhsevl_assembleSpatialRHS (rglobalData, rrhs, dtimePrimal2, dtimeDual2,&
              roptimalControl, rtempVector3)
        end if
        
      else
      
        ! We are 'at the end'.
        !
        ! Dual RHS comes from rtempVector3. The primal from the
        ! iiterate-1'th RHS and rtempVector3.
        !
        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
        ! dual RHS(0)   = DUALRHS(0)
      
        call lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(1),rtempVector3%RvectorBlock(1),&
            dweightoldp,dweightnewp,&
            rtempVectorRHS%RvectorBlock(1))                                        
        call lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(2),rtempVector3%RvectorBlock(2),&
            dweightoldp,dweightnewp,&
            rtempVectorRHS%RvectorBlock(2))                                        
        ! Pressure is fully implicit
        call lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(3),rtempVector3%RvectorBlock(3),&
            dweightoldp,dweightnewp,&
            rtempVectorRHS%RvectorBlock(3))

!        !CALL generateRHS (rproblem,iiterate+1,niterations,&
!        !    rtempVector3, .TRUE., .FALSE.)
        call lsyssc_copyVector (rtempVector3%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
        call lsyssc_copyVector (rtempVector3%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
        call lsyssc_copyVector (rtempVector3%RvectorBlock(6),rtempVectorRHS%RvectorBlock(6))

        ! Multiply the last RHS of the dual equation -z by theta+gamma/dtstep, that's it.
        call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),&
            dweightnewd + roptimalControl%dgammaC/dtstep)
        call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),&
            dweightnewd + roptimalControl%dgammaC/dtstep)

      end if

      ! Implement the boundary conditions into the RHS vector        
      if (present(roptcBDC)) then
        call tbc_implementSpatialBCtoRHS (roptcBDC, &
            dtimePrimal, dtimeDual, rrhsDiscrete%p_rtimeDiscr, rtempVectorRHS, rglobalData)
      end if
      
      ! Save the RHS.
      call sptivec_setTimestepData(rrhsDiscrete, iiterate, rtempVectorRHS)
      
    end do
    
    ! Release the temp vectors.
    call lsysbl_releaseVector (rtempVectorRHS)
    
    call lsysbl_releaseVector (rtempVector3)
    call lsysbl_releaseVector (rtempVector2)
    call lsysbl_releaseVector (rtempVector1)

  end subroutine
  !2. Rhs ist falsch

!  ! ***************************************************************************
!  
!!<subroutine>
!
!  subroutine trhsevl_assembledG0RHS (rproblem, rspaceTimeDiscr, rb, bimplementBC)
!
!!<description>
!  ! Assembles the space-time RHS vector rb according to the dG(0)-scheme.
!  !
!  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
!  ! this routine!
!!</description>
!
!!<input>
!  ! A problem structure that provides information on all
!  ! levels as well as temporary vectors.
!  type(t_problem), intent(INOUT), target :: rproblem
!  
!  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
!  ! coupled space-time system.
!  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
!  
!  ! Whether to implement boundary conditions into the RHS or not.
!  logical, intent(IN) :: bimplementBC
!!</input>
!
!!<inputoutput>
!  ! A space-time vector that receives the RHS.
!  type(t_spacetimeVector), intent(INOUT) :: rb
! 
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: isubstep,nintervals
!    real(DP) :: dtstep
!    
!    ! Temporary vector
!    type(t_vectorBlock) :: rtempVector
!
!    real(DP), dimension(:),pointer :: p_Dx, p_Db, p_Dd, p_Drhs
!
!    ! Create a temp vector for the assembly
!    call lsysbl_createVecBlockByDiscr (&
!        p_rlevelInfo%rdiscretisation,rtempVector,.true.)
!
!    call lsysbl_getbase_double (rtempVector,p_Db)
!
!    ! The assembly for dG0-RHS is rather easy.
!    ! The i'th time DOF belongs to the i'th time interval Ti and
!    ! has the following form:
!    !  ___
!    !  f_i  =  1/|Ti| int_Ti f(.,t) dt  ~  f(.,T_i(midpoint))
!    !
!    ! by the midpoint rule in time! So we just make a loop
!    ! over the timesteps and calculate the f()-values in the
!    ! midpoints of the time intervals!
!    
!    dtstep = rrhsDiscrete%p_rtimeDiscr%dtstep
!    nintervals = rrhsDiscrete%p_rtimeDiscr%nintervals
!
!    do isubstep = 0,nintervals-1
!      ! Assemble at the midpoint of the time interval
!      call trhsevl_assembleSpatialRHS (rproblem,isubstep,&
!        (real(isubstep,DP)+0.5_DP)*dtstep,rtempVector)
!        
!      call sptivec_setTimestepData(rb, 1+isubstep, rtempVector)
!    end do
!      
!    ! Release the temp vector.
!    call lsysbl_releaseVector (rtempVector)
!
!  end subroutine
!
  ! ***************************************************************************

  subroutine trhsevl_evalFunction (rdiscretisation, rform, &
      nelements, npointsPerElement, Dpoints, &
      IdofsTest, rdomainIntSubset, &
      Dcoefficients, rcollection)
  
  ! Standard evaluation routine. Evaluates component rcollection%Icollection(1)
  ! of the analytical solution identified by the name "RHS" in the collection.
  
  use fsystem
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
  type(t_linearForm), intent(IN) :: rform
  integer, intent(IN) :: nelements
  integer, intent(IN) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
  integer, dimension(:,:), intent(IN) :: IdofsTest
  type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
  type(t_collection), intent(INOUT), optional :: rcollection
  real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
  
  integer :: ierror
  integer, dimension(3) :: Ibounds
  real(dp), dimension(:,:), allocatable :: Dvalues

    ! Prepare a target array.
    Ibounds = ubound(Dcoefficients)
    allocate(Dvalues(Ibounds(2),Ibounds(3)))

    ! Evaluate
    call ansol_evaluate (rcollection,"RHS",rcollection%IquickAccess(1),&
        Dvalues,npointsPerElement,nelements,Dpoints,rdomainIntSubset%p_Ielements,ierror)
        
    ! Check that this was ok. If yes, copy the data to the destination.
    if (ierror .eq. 0) then
      Dcoefficients(1,:,:) = Dvalues(:,:)
    else
      call output_line ('Error evaluating RHS function.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'trhsevl_evalFunction')
      call sys_halt()
    end if
        
    ! Clean up
    deallocate(Dvalues)

  end subroutine
  
!  ! ***************************************************************************
!
!  subroutine trhsevl_evalTargetFlow (rdiscretisation, rform, &
!      nelements, npointsPerElement, Dpoints, &
!      IdofsTest, rdomainIntSubset, &
!      Dcoefficients, rcollection)
!  
!  ! Standard evaluation routine. Evaluates component rcollection%Icollection(1)
!  ! of the analytical solution identified by the name "RHS" in the collection.
!  
!  use fsystem
!  use basicgeometry
!  use triangulation
!  use scalarpde
!  use domainintegration
!  use spatialdiscretisation
!  use collection
!  
!  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
!  type(t_linearForm), intent(IN) :: rform
!  integer, intent(IN) :: nelements
!  integer, intent(IN) :: npointsPerElement
!  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!  integer, dimension(:,:), intent(IN) :: IdofsTest
!  type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
!  type(t_collection), intent(INOUT), optional :: rcollection
!  real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
!  
!  integer :: ierror,i,j
!  integer, dimension(3) :: Ibounds
!  real(dp), dimension(:,:), allocatable :: Dvalues
!
!    ! Prepare a target array.
!    Ibounds = ubound(Dcoefficients)
!    allocate(Dvalues(Ibounds(2),Ibounds(3)))
!
!    ! Evaluate
!    call ansol_evaluate (rcollection,"RHS",rcollection%IquickAccess(1),&
!        Dvalues,npointsPerElement,nelements,Dpoints,rdomainIntSubset%p_Ielements,ierror)
!    
!    do i=1,nelements
!      do j=1,npointsPerElement
!        if ( (Dpoints(1,j,i)-0.5_DP)**2 + (Dpoints(2,j,i)-0.2_DP)**2 .le. 0.05**2) then
!          Dvalues(j,i) = Dvalues(j,i)*max(0.0_DP,1.0_DP-rcollection%DquickAccess(1))
!        end if
!      end do
!    end do
!    
!    ! Check that this was ok. If yes, copy the data to the destination.
!    if (ierror .eq. 0) then
!      Dcoefficients(1,:,:) = Dvalues(:,:)
!    else
!      call output_line ('Error evaluating RHS function.', &
!                        OU_CLASS_ERROR,OU_MODE_STD,'trhsevl_evalFunction')
!      call sys_halt()
!    end if
!        
!    ! Clean up
!    deallocate(Dvalues)
!
!  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine trhsevl_assembleSpatialRHS (rglobalData, rrhs, dtimePrimal, dtimeDual, &
      roptimalControl, rrhsDiscrete)
  
!<description>
  ! Generate the RHS vector at the time dtime. isubstep may specify the
  ! timestep where to generate the RHS.
!</description>
 
!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData

  ! Analytic solution defining the RHS of the equation.
  type(t_anSolution), intent(inout) :: rrhs

  ! Time where the primal/dual RHS should be generated.
  ! Must not necessarily coincide with the start/end time of the timestep.
  real(DP), intent(in) :: dtimePrimal
  real(DP), intent(in) :: dtimeDual

  ! Optimal control parameters.
  type(t_settings_optcontrol), intent(inout) :: roptimalControl
!</input>  

!<inputoutput>  
  ! Destination vector
  type(t_vectorBlock), intent(inout) :: rrhsDiscrete
!</inputoutput>
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A collection structure for the callback routines.
    type(t_collection) :: rcollection
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Drhs
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rrhsDiscrete,p_Drhs)

    ! Get a pointer to the RHS on the finest level as well as to the
    ! block discretisation structure:
    p_rdiscretisation => rrhsDiscrete%p_rblockDiscr
    
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
    call collct_init(rcollection)
    
    if (rrhs%ctype .eq. ANSOL_TP_ANALYTICAL) then
    
      ! Evaluate the RHS using the user defined callback functions
      call user_initCollectForAssembly (rglobalData,dtimePrimal,rcollection)

      ! Discretise the X-velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
                rrhsDiscrete%RvectorBlock(1),user_coeff_RHSprimal_x,rcollection)

      ! And the Y-velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
                rrhsDiscrete%RvectorBlock(2),user_coeff_RHSprimal_y,rcollection)

      call user_doneCollectForAssembly (rglobalData,rcollection)

      ! Dual RHS

      call user_initCollectForAssembly (rglobalData,dtimeDual,rcollection)
                                  
      ! The RHS terms for the dual equation are calculated using
      ! the desired 'target' flow field plus the coefficients of the
      ! dual RHS -- which are normally zero.
      !
      ! Discretise the X-velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(4),rlinform,.true.,&
                rrhsDiscrete%RvectorBlock(4),user_coeff_TARGET_x,rcollection)
      
      ! And the Y-velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(5),rlinform,.true.,&
                rrhsDiscrete%RvectorBlock(5),user_coeff_TARGET_y,rcollection)

      call user_doneCollectForAssembly (rglobalData,rcollection)
                
    else
    
      ! Evaluate the RHS using the analytical definition.
      call ansol_prepareEval (rrhs,rcollection,"RHS",dtimePrimal)
      
      ! Discretise the primal X-velocity part:
      rcollection%IquickAccess(1) = 1
      call linf_buildVectorScalar (&
          p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
          rrhsDiscrete%RvectorBlock(1),trhsevl_evalFunction,rcollection)

      ! And the primal Y-velocity part:
      rcollection%IquickAccess(1) = 2
      call linf_buildVectorScalar (&
          p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
          rrhsDiscrete%RvectorBlock(2),trhsevl_evalFunction,rcollection)
      
      call ansol_doneEval (rcollection,"RHS")
    
      ! Evaluate the target function using the analytical definition.
      call ansol_prepareEval (roptimalControl%rtargetFlow,rcollection,"RHS",dtimeDual)
      
      ! Discretise the X-velocity part:
      rcollection%IquickAccess(1) = 1
      call linf_buildVectorScalar (&
          p_rdiscretisation%RspatialDiscr(4),rlinform,.true.,&
          rrhsDiscrete%RvectorBlock(4),trhsevl_evalFunction,rcollection)

      ! And the Y-velocity part:
      rcollection%IquickAccess(1) = 2
      call linf_buildVectorScalar (&
          p_rdiscretisation%RspatialDiscr(5),rlinform,.true.,&
          rrhsDiscrete%RvectorBlock(5),trhsevl_evalFunction,rcollection)
          
      call ansol_doneEval (rcollection,"RHS")
    
    end if
      
    ! Depending on the formulation, to get a reference dual velocity,
    ! it might be necessary to switch the sign of the target velocity field 
    ! because the RHS of the dual equation is '-z'!
    ! Remember that it this case the signs of the mass matrices that couple
    ! primal and dual velocity must be changed, too!
    if (roptimalControl%ispaceTimeFormulation .ne. 0) then
      call lsyssc_scaleVector (rrhsDiscrete%RvectorBlock(4),-1.0_DP)
      call lsyssc_scaleVector (rrhsDiscrete%RvectorBlock(5),-1.0_DP)
    end if

    ! Now there may exist a 'real' dual RHS (ehich is usually only the case
    ! for analytical test functions). We have to add these to the dual
    ! RHS vectors. This will compute the real dual RHS "f_dual - z".
    if (rrhs%ctype .eq. ANSOL_TP_ANALYTICAL) then
    
      ! Evaluate the RHS using the user defined callback functions
      call user_initCollectForAssembly (rglobalData,dtimeDual,rcollection)

      ! Discretise the X-velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(4),rlinform,.false.,&
                rrhsDiscrete%RvectorBlock(4),user_coeff_RHSdual_x,rcollection)
      
      ! And the Y-velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(5),rlinform,.false.,&
                rrhsDiscrete%RvectorBlock(5),user_coeff_RHSdual_y,rcollection)

      call user_doneCollectForAssembly (rglobalData,rcollection)
                
    else
    
      ! Evaluate the RHS using the analytical definition.
      call ansol_prepareEval (rrhs,rcollection,"RHS",dtimeDual)
      
      ! Discretise the dual X-velocity part:
      rcollection%IquickAccess(1) = 4
      call linf_buildVectorScalar (&
          p_rdiscretisation%RspatialDiscr(4),rlinform,.false.,&
          rrhsDiscrete%RvectorBlock(4),trhsevl_evalFunction,rcollection)

      ! And the dual Y-velocity part:
      rcollection%IquickAccess(1) = 5
      call linf_buildVectorScalar (&
          p_rdiscretisation%RspatialDiscr(5),rlinform,.false.,&
          rrhsDiscrete%RvectorBlock(5),trhsevl_evalFunction,rcollection)

      call ansol_doneEval (rcollection,"RHS")
    
    end if


    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhsDiscrete%RvectorBlock(3))
    
    ! Dual pressure RHS is =0 as well.
    call lsyssc_clearVector(rrhsDiscrete%RvectorBlock(6))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    call collct_done(rcollection)
  
  end subroutine
    
end module
