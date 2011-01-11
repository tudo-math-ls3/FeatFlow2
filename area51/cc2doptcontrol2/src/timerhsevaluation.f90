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
  use fespacehierarchybase
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

  subroutine trhsevl_assembleRHS (rglobalData, rphysics, rrhs, rrhsDiscrete, &
      roptimalControl, roptcBDC)

!<description>
  ! Discretises the RHS according to the space/time discretisation scheme
  ! in rrhsDiscrete. The result is summed up to rrhsDiscrete.
!</description>

!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
  
  ! Physics of the underlying equation.
  type(t_settings_physics), intent(in) :: rphysics

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
      call trhsevl_assembleThetaRHS (rglobalData, rphysics, rrhs, rrhsDiscrete, roptimalControl, roptcBDC)
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

  subroutine trhsevl_assembleThetaRHS (rglobalData, rphysics, rrhs, rrhsDiscrete, roptimalControl, roptcBDC)

!<description>
  ! Assembles the space-time RHS vector rrhsDiscrete for a Theta-Scheme. 
!</description>

!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
  
  ! Physics of the underlying equation.
  type(t_settings_physics), intent(in) :: rphysics

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
    real(DP) :: dtimePrimal,dtstep,dtimeDual
    integer :: ithetaschemetype
    real(DP) :: dtheta
    
    ! Temporary space.
    type(t_spaceTimeVectorAccess) :: raccessPool
    type(t_vectorBlock), pointer :: p_rvecLast,p_rvecCurrent,p_rvecNext
    
    ! A temporary vector for the creation of the RHS.
    type(t_vectorBlock) :: rtempVectorRHS
    
    real(DP), dimension(:),pointer :: p_Dx1, p_Dx2, p_Dx3, p_Drhs

    nintervals = rrhsDiscrete%p_rtimeDiscr%nintervals
    
    ! Type of the theta scheme.
    ! =0: Dual points in time located inbetween primal points in time.
    ! =1: Dual points in time located on primal points in time.
    !     This does not correspond to any known analytical minimisation problem.
    ithetaschemetype = rrhsDiscrete%p_rtimeDiscr%itag
    dtheta = rrhsDiscrete%p_rtimeDiscr%dtheta
    
    ! Temp vector
    call lsysbl_createVectorBlock (rrhsDiscrete%p_rspaceDiscr,rtempVectorRHS,.true.)
    call lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
    
    ! HERE!!!
    ! das ist irgendwie quatsch hier. der vorletzte Vektor wird nicht richtig behandelt!
    
    ! ----------------------------------------------------------------------
    ! Generate the global RHS vector
    
    ! Temporary space, ring buffer. Initialise large enough to hold all
    ! temporary data: Last, current and next vector.
    call sptivec_createAccessPool(rrhsDiscrete%p_rspaceDiscr,raccessPool,3)

    ! NOTE:
    ! For setting up the RHS, the time in the dual equation matches the time
    ! of the primal equation. For explicit Euler, this is the case anyway.
    ! For Crank-Nicolson / Theta scheme, the RHS must be evaluated in the
    ! time midpoint between two subsequent timesteps; but as the evaluation
    ! points of the dual equation are shifted by dt/2, the time midpoint
    ! in the dual equation coincides with the evaluation point in time of the
    ! primal equation!

    do iiterate = 1,nintervals+1
    
      ! Get the previous, current and next rhs, relative to the current iterate.
      ! If necessary, assemble.
      if (iiterate .gt. 1) then
        call sptivec_getVectorFromPool (raccessPool,iiterate-1,p_rvecLast)
        if (.not. associated(p_rvecLast)) then
          call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate-2,dtimePrimal,dtstep)
          dtimeDual = dtimePrimal
          call sptivec_getFreeBufferFromPool (raccessPool,iiterate-1,p_rvecLast)
          call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
              roptimalControl, p_rvecLast)
        end if
        call lsysbl_getbase_double (p_rvecLast,p_Dx1)
      end if

      call sptivec_getVectorFromPool (raccessPool,iiterate,p_rvecCurrent)
      if (.not. associated(p_rvecCurrent)) then
        call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate-1,dtimePrimal,dtstep)
        dtimeDual = dtimePrimal
        call sptivec_getFreeBufferFromPool (raccessPool,iiterate,p_rvecCurrent)
        call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
            roptimalControl, p_rvecCurrent)
      end if
      call lsysbl_getbase_double (p_rvecCurrent,p_Dx2)
      
      if (iiterate .lt. nintervals+1) then
        call sptivec_getVectorFromPool (raccessPool,iiterate+1,p_rvecNext)
        if (.not. associated(p_rvecNext)) then
          call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate,dtimePrimal,dtstep)
          dtimeDual = dtimePrimal
          call sptivec_getFreeBufferFromPool (raccessPool,iiterate+1,p_rvecNext)
          call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
              roptimalControl, p_rvecNext)
        end if
        call lsysbl_getbase_double (p_rvecNext,p_Dx3)
      endif
      
      ! Scale the dual RHS according to the timestep scheme.
      select case (ithetaschemetype)
      case (0)
        select case (rphysics%cequation)
        case (0,1)
          ! Stokes, Navier-Stokes
          if (iiterate .eq. 1) then
          
            ! RHS comes from first vector.
            !
            ! primal RHS(0) = PRIMALRHS(0)
            ! dual RHS(0)   = DUALRHS(0)
        
            call lsysbl_copyVector (p_rvecCurrent,rtempVectorRHS)
            call lsyssc_clearVector (rtempVectorRHS%RvectorBlock(4))
            call lsyssc_clearVector (rtempVectorRHS%RvectorBlock(5))
            
          else if (iiterate .lt. nintervals+1) then
          
            ! Somewhere 'in the middle'
            
            ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
            ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
            
            call lsyssc_vectorLinearComb (&
                p_rvecLast%RvectorBlock(1),p_rvecCurrent%RvectorBlock(1),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(1))                                        
            call lsyssc_vectorLinearComb (&                                                   
                p_rvecLast%RvectorBlock(2),p_rvecCurrent%RvectorBlock(2),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(2))                                        
            ! Divergence equation is fully implicit
            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))

            call lsyssc_vectorLinearComb (&
                p_rvecNext%RvectorBlock(4),p_rvecCurrent%RvectorBlock(4),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(4))                                        
            call lsyssc_vectorLinearComb (&                                                   
                p_rvecNext%RvectorBlock(5),p_rvecCurrent%RvectorBlock(5),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(5))                                        
            ! Divergence equation is fully implicit
            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(6),rtempVectorRHS%RvectorBlock(6))
            
          else
          
            ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
            ! dual RHS(0)   = DUALRHS(0)
          
            call lsyssc_vectorLinearComb (&
                p_rvecLast%RvectorBlock(1),p_rvecCurrent%RvectorBlock(1),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(1))                                        
            call lsyssc_vectorLinearComb (&                                                   
                p_rvecLast%RvectorBlock(2),p_rvecCurrent%RvectorBlock(2),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(2))                                        
            ! Divergence equation is fully implicit
            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))

            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
            call lsyssc_vectorLinearComb (&                                                   
                p_rvecLast%RvectorBlock(6),p_rvecCurrent%RvectorBlock(6),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(6))

            ! Multiply the last RHS of the dual equation -z by theta+gamma/dtstep, that is it.
            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),&
                1.0_DP + roptimalControl%dgammaC/dtstep)
            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),&
                1.0_DP + roptimalControl%dgammaC/dtstep)
          end if
        end select
        
      case (1)
        select case (rphysics%cequation)
        case (0,1)
          ! Stokes, Navier-Stokes
          if (iiterate .eq. 1) then
          
            ! RHS comes from first vector.
            !
            ! primal RHS(0) = PRIMALRHS(0)
            ! dual RHS(0)   = DUALRHS(0)
            call lsysbl_copyVector (p_rvecCurrent,rtempVectorRHS)
            
            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),1.0_DP-dtheta)
            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),1.0_DP-dtheta)
            
          else 
          
            ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
            ! dual RHS(0)   = DUALRHS(0)
          
            call lsyssc_vectorLinearComb (&
                p_rvecLast%RvectorBlock(1),p_rvecCurrent%RvectorBlock(1),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(1))                                        
            call lsyssc_vectorLinearComb (&                                                   
                p_rvecLast%RvectorBlock(2),p_rvecCurrent%RvectorBlock(2),&
                1.0_DP-dtheta,dtheta,&
                rtempVectorRHS%RvectorBlock(2))                                        
            ! Divergence equation is fully implicit
            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))

            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
            call lsyssc_vectorLinearComb (&                                                   
                p_rvecLast%RvectorBlock(6),p_rvecCurrent%RvectorBlock(6),&
                1.0_DP-dtheta,dtheta,rtempVectorRHS%RvectorBlock(6))
          
            ! In the last two steps, we have to take care of the terminal condition!
            if (iiterate .eq. nintervals) then
          
              ! 1st part of the terminal condition
              call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),&
                  (1.0_DP + (1.0_DP-dtheta)*roptimalControl%dgammaC/dtstep))

              call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),&
                  (1.0_DP + (1.0_DP-dtheta)*roptimalControl%dgammaC/dtstep))
          
            else if (iiterate .eq. nintervals+1) then

              ! 2nd part of the terminal condition
              call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),&
                  (dtheta + dtheta*roptimalControl%dgammaC/dtstep))

              call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),&
                  (dtheta + dtheta*roptimalControl%dgammaC/dtstep))
            
            end if          
            
          end if
          
        end select
      end select
      
!    ! Assemble 1st RHS vector in X temp vector.
!    call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,0,dtimePrimal,dtstep)
!    dtimeDual = dtimePrimal
!    call sptivec_getFreeBufferFromPool (raccessPool,1,p_rvecCurrent)
!    call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
!        roptimalControl, p_rvecCurrent)
!        
!    ! Assemble the 2nd RHS vector in the RHS temp vector
!    call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,1,dtimePrimal,dtstep)
!    dtimeDual = dtimePrimal
!    call sptivec_getFreeBufferFromPool (raccessPool,2,p_rvecCurrent)
!    call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
!        roptimalControl, p_rvecCurrent)
!
!    call lsysbl_createVectorBlock (rrhsDiscrete%p_rspaceDiscr,rtempVectorRHS,.true.)
!    dtheta = rrhsDiscrete%p_rtimeDiscr%dtheta
!    
!    ! DEBUG!!!
!    call lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
!    
!    do iiterate = 1,nintervals+1
!    
!      if (iiterate .eq. 1) then
!      
!        ! RHS comes from first vector.
!        !
!        ! primal RHS(0) = PRIMALRHS(0)
!        ! dual RHS(0)   = DUALRHS(0)
!
!        call sptivec_getVectorFromPool (raccessPool,iiterate,p_rvecCurrent)
!        call lsysbl_copyVector (p_rvecCurrent,rtempVectorRHS)
!        
!        ! Scale the dual RHS according to the timestep scheme.
!        select case (ithetaschemetype)
!        case (0)
!          select case (rphysics%cequation)
!          case (0,1)
!            ! Stokes, Navier-Stokes
!            call lsyssc_clearVector (rtempVectorRHS%RvectorBlock(4))
!            call lsyssc_clearVector (rtempVectorRHS%RvectorBlock(5))
!          end select
!        case (1)
!          select case (rphysics%cequation)
!          case (0,1)
!            ! Stokes, Navier-Stokes
!            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),1.0_DP-dtheta)
!            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),1.0_DP-dtheta)
!          end select
!        end select
!        
!      else if (iiterate .lt. nintervals+1) then
!      
!        ! We are somewhere 'in the middle'.
!        !
!        ! Dual RHS comes from rtempVector3. The primal from the
!        ! iiterate-1'th RHS.
!        !
!        ! Get the vectors. Create missing vectors.
!        call sptivec_getVectorFromPool (raccessPool,iiterate-1,p_rvecLast)
!        
!        call sptivec_getVectorFromPool (raccessPool,iiterate,p_rvecCurrent)
!        if (.not. associated(p_rvecCurrent)) then
!          call sptivec_getFreeBufferFromPool (raccessPool,iiterate,p_rvecCurrent)
!          call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate-1,dtimePrimal,dtstep)
!          dtimeDual = dtimePrimal
!          call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
!              roptimalControl, p_rvecCurrent)
!        end if
!
!        call sptivec_getVectorFromPool (raccessPool,iiterate+1,p_rvecNext)
!        if (.not. associated(p_rvecNext)) then
!          call sptivec_getFreeBufferFromPool (raccessPool,iiterate+1,p_rvecNext)
!          call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate,dtimePrimal,dtstep)
!          dtimeDual = dtimePrimal
!          call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
!              roptimalControl, p_rvecNext)
!        end if
!        
!        select case (ithetaschemetype)
!        case (0)
!        
!         
!       case (1)
!        
!          ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!          ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
!        
!          select case (rphysics%cequation)
!          case (0,1)
!            ! Stokes, Navier-Stokes
!        
!          end select
!
!        end select
!        
!      else
!      
!        ! We are 'at the end'.
!        !
!        ! Dual RHS comes from rtempVector3. The primal from the
!        ! iiterate-1'th RHS and rtempVector3.
!        
!        call sptivec_getVectorFromPool (raccessPool,iiterate-1,p_rvecLast)
!        
!        call sptivec_getVectorFromPool (raccessPool,iiterate,p_rvecCurrent)
!        if (.not. associated(p_rvecCurrent)) then
!          call sptivec_getFreeBufferFromPool (raccessPool,iiterate,p_rvecCurrent)
!          call tdiscr_getTimestep(rrhsDiscrete%p_rtimeDiscr,iiterate,dtimePrimal,dtstep)
!          dtimeDual = dtimePrimal
!          call trhsevl_assembleSpatialRHS (rglobalData, rphysics, .true., rrhs, dtimePrimal, dtimeDual,&
!              roptimalControl, p_rvecCurrent)
!        end if
!        
!        select case (ithetaschemetype)
!        case (0)
!        
!          ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!          ! dual RHS(0)   = DUALRHS(0)
!          
!          select case (rphysics%cequation)
!          case (0,1)
!            ! Stokes, Navier-Stokes
!        
!            ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!            ! dual RHS(0)   = DUALRHS(0)
!          
!            call lsyssc_vectorLinearComb (&
!                p_rvecLast%RvectorBlock(1),p_rvecCurrent%RvectorBlock(1),&
!                1.0_DP-dtheta,dtheta,&
!                rtempVectorRHS%RvectorBlock(1))                                        
!            call lsyssc_vectorLinearComb (&                                                   
!                p_rvecLast%RvectorBlock(2),p_rvecCurrent%RvectorBlock(2),&
!                1.0_DP-dtheta,dtheta,&
!                rtempVectorRHS%RvectorBlock(2))                                        
!            ! Divergence equation is fully implicit
!            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))
!
!            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
!            call lsyssc_copyVector (p_rvecCurrent%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
!            call lsyssc_vectorLinearComb (&                                                   
!                p_rvecLast%RvectorBlock(6),p_rvecCurrent%RvectorBlock(6),&
!                1.0_DP-dtheta,dtheta,&
!                rtempVectorRHS%RvectorBlock(6))
!                
!            ! Multiply the last RHS of the dual equation -z by theta+gamma/dtstep, that is it.
!            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),&
!                1.0_DP + roptimalControl%dgammaC/dtstep)
!            call lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),&
!                1.0_DP + roptimalControl%dgammaC/dtstep)
!                
!          end select
!
!        case (1)
!
!          ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!          ! dual RHS(0)   = DUALRHS(0)
!          
!          select case (rphysics%cequation)
!          case (0,1)
!            ! Stokes, Navier-Stokes
!        
!          end select
!
!        end select
!
!      end if

      ! If we work woth a scaled system, scale the RHS by dt.
      if (roptimalControl%csystemScaling .ne. 0) then
        call lsysbl_scaleVector (rtempVectorRHS,dtstep)
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
    call sptivec_releaseAccessPool (raccessPool)

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
  
  subroutine trhsevl_assembleSpatialRHS (rglobalData, rphysics,&
      bclear, rrhs, dtimePrimal, dtimeDual, &
      roptimalControl, rrhsDiscrete)
  
!<description>
  ! Generate the RHS vector at the time dtime. isubstep may specify the
  ! timestep where to generate the RHS.
!</description>
 
!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData

  ! Clear old RHS vector.
  logical, intent(in) :: bclear
  
  ! Physics of the underlying equation.
  type(t_settings_physics), intent(in) :: rphysics

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
    
    select case (rphysics%cequation)
    
    case (0,1)
      ! Stokes and Navier--Stokes, 2D.
      
      ! Component 1: X-velocity.
      ! Component 2: Y-velocity.
      ! Component 3: pressure
      ! Component 4: dual X-velocity.
      ! Component 5: dual Y-velocity.
      ! Component 6: dual pressure.
    
      ! The third subvector must be zero initially - as it represents the RHS of
      ! the equation "div(u) = 0".
      call lsyssc_clearVector(rrhsDiscrete%RvectorBlock(3))
      
      ! Dual pressure RHS is =0 as well.
      call lsyssc_clearVector(rrhsDiscrete%RvectorBlock(6))
      
      if (rrhs%ctype .eq. ANSOL_TP_ANALYTICAL) then
      
        ! Evaluate the RHS using the user defined callback functions
        call user_initCollectForVecAssembly (rglobalData,rrhs%iid,1,dtimePrimal,rcollection)

        ! Discretise the X-velocity part:
        call linf_buildVectorScalar (&
                  p_rdiscretisation%RspatialDiscr(1),rlinform,bclear,&
                  rrhsDiscrete%RvectorBlock(1),user_coeff_RHS,rcollection)

        call user_doneCollectForAssembly (rglobalData,rcollection)
        
        call user_initCollectForVecAssembly (rglobalData,rrhs%iid,2,dtimePrimal,rcollection)

        ! And the Y-velocity part:
        call linf_buildVectorScalar (&
                  p_rdiscretisation%RspatialDiscr(2),rlinform,bclear,&
                  rrhsDiscrete%RvectorBlock(2),user_coeff_RHS,rcollection)

        call user_doneCollectForAssembly (rglobalData,rcollection)

      else
      
        ! Evaluate the RHS using the analytical definition.
        call ansol_prepareEval (rrhs,rcollection,"RHS",dtimePrimal)
        
        ! Discretise the primal X-velocity part:
        rcollection%IquickAccess(1) = 1
        call linf_buildVectorScalar (&
            p_rdiscretisation%RspatialDiscr(1),rlinform,bclear,&
            rrhsDiscrete%RvectorBlock(1),trhsevl_evalFunction,rcollection)

        ! And the primal Y-velocity part:
        rcollection%IquickAccess(1) = 2
        call linf_buildVectorScalar (&
            p_rdiscretisation%RspatialDiscr(2),rlinform,bclear,&
            rrhsDiscrete%RvectorBlock(2),trhsevl_evalFunction,rcollection)

        ! Dual pressure RHS.
        rcollection%IquickAccess(1) = 6
        call linf_buildVectorScalar (&
            p_rdiscretisation%RspatialDiscr(6),rlinform,bclear,&
            rrhsDiscrete%RvectorBlock(6),trhsevl_evalFunction,rcollection)
        
        call ansol_doneEval (rcollection,"RHS")
        
      end if
      
      if (roptimalControl%rtargetFunction%ctype .eq. ANSOL_TP_ANALYTICAL) then
      
        ! Dual RHS, target function

        call user_initCollectForVecAssembly (rglobalData,&
            roptimalControl%rtargetFunction%iid,1,dtimeDual,rcollection)
                                    
        ! The RHS terms for the dual equation are calculated using
        ! the desired 'target' function plus the coefficients of the
        ! dual RHS -- which are normally zero.
        !
        ! Discretise the X-velocity part:
        call linf_buildVectorScalar (&
                  p_rdiscretisation%RspatialDiscr(4),rlinform,bclear,&
                  rrhsDiscrete%RvectorBlock(4),user_coeff_Target,rcollection)

        call user_doneCollectForAssembly (rglobalData,rcollection)

        call user_initCollectForVecAssembly (rglobalData,&
            roptimalControl%rtargetFunction%iid,2,dtimeDual,rcollection)
        
        ! And the Y-velocity part:
        call linf_buildVectorScalar (&
                  p_rdiscretisation%RspatialDiscr(5),rlinform,bclear,&
                  rrhsDiscrete%RvectorBlock(5),user_coeff_Target,rcollection)

        call user_doneCollectForAssembly (rglobalData,rcollection)
        
      else
                        
        ! Evaluate the target function using the analytical definition.
        call ansol_prepareEval (roptimalControl%rtargetFunction,rcollection,"RHS",dtimeDual)
        
        ! Discretise the X-velocity part:
        rcollection%IquickAccess(1) = 1
        call linf_buildVectorScalar (&
            p_rdiscretisation%RspatialDiscr(4),rlinform,bclear,&
            rrhsDiscrete%RvectorBlock(4),trhsevl_evalFunction,rcollection)

        ! And the Y-velocity part:
        rcollection%IquickAccess(1) = 2
        call linf_buildVectorScalar (&
            p_rdiscretisation%RspatialDiscr(5),rlinform,bclear,&
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
      ! for analytical test functions). We have to add this to the dual
      ! RHS vectors. This will compute the real dual RHS "f_dual - z".
      if (rrhs%ctype .eq. ANSOL_TP_ANALYTICAL) then
      
        ! Evaluate the RHS using the user defined callback functions
        call user_initCollectForVecAssembly (rglobalData,rrhs%iid,4,dtimeDual,rcollection)

        ! Discretise the X-velocity part:
        call linf_buildVectorScalar (&
                  p_rdiscretisation%RspatialDiscr(4),rlinform,.false.,&
                  rrhsDiscrete%RvectorBlock(4),user_coeff_RHS,rcollection)
        
        call user_doneCollectForAssembly (rglobalData,rcollection)
        
        call user_initCollectForVecAssembly (rglobalData,rrhs%iid,5,dtimeDual,rcollection)
        
        ! And the Y-velocity part:
        call linf_buildVectorScalar (&
                  p_rdiscretisation%RspatialDiscr(5),rlinform,.false.,&
                  rrhsDiscrete%RvectorBlock(5),user_coeff_RHS,rcollection)

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

        ! Primal pressure:
        rcollection%IquickAccess(1) = 3
        call linf_buildVectorScalar (&
            p_rdiscretisation%RspatialDiscr(3),rlinform,.false.,&
            rrhsDiscrete%RvectorBlock(3),trhsevl_evalFunction,rcollection)

        call ansol_doneEval (rcollection,"RHS")
      
      end if


    end select
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    call collct_done(rcollection)
  
  end subroutine
    
end module
