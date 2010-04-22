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

module timeboundaryconditions


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

  !use spacepreconditioner
  !use spacepreconditionerinit
  use timeanalysis
  use spatialbc
  use spacediscretisation
  
  use timediscretisation
  use spacetimevectors
  use dofmapping
  
  use spacetimediscretisation

  implicit none

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementBCsolution (rproblem,rspaceTimeDiscr,rx,rtempVectorX)

!<description>
  ! Implements the boundary conditions of all timesteps into the solution rx.
!</description>

!<input>
  ! Problem structure of the main problem.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented
  ! to.
  type(t_spacetimeVector), intent(INOUT) :: rx

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  type(t_vectorBlock), intent(INOUT), optional :: rtempVectorX
!</inputoutput>

!</subroutine>

    integer :: isubstep
    real(DP) :: dtstep,dtime
    type(t_vectorBlock) :: rtempVector

    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata

    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    
    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)
    
    if (present(rtempVectorX)) then
      call lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    else
      ! Create a temp vector
      call lsysbl_createVecBlockByDiscr (&
          rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.true.)
    end if
    call lsysbl_assignDiscreteBC(rtempVector,rdiscreteBC)
    call lsysbl_assignDiscreteFBC(rtempVector,rdiscreteFBC)

    call lsyssc_getbase_double (rtempVector%RvectorBlock(1),p_Ddata)

    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep

    ! The implementation of the boundary conditions depends on the type
    ! of the time discretisation...
    select case (rspaceTimeDiscr%rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)

      ! Implement the bondary conditions into all initial solution vectors
      do isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = rproblem%rtimedependence%dtimeInit + isubstep*dtstep

        ! Initialise the collection for the assembly.
        call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)
      
        ! -----
        ! Discretise the boundary conditions at the new point in time.
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call cc_assembleBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
        call cc_assembleFBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
        
        ! Implement the boundary conditions into the global solution vector.
        call sptivec_getTimestepData(rx, 1+isubstep, rtempVector)
        
        call vecfil_discreteBCsol (rtempVector)
        call vecfil_discreteFBCsol (rtempVector)
        
        call sptivec_setTimestepData(rx, 1+isubstep, rtempVector)
        
        ! Clean up the collection (as we are done with the assembly, that's it.
        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

      end do
      
    case (TDISCR_DG0)

      ! Implement the bondary conditions into all initial solution vectors
      do isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals-1
      
        ! Current point in time
        dtime = rproblem%rtimedependence%dtimeInit + (real(isubstep,DP)+0.5_DP)*dtstep

        ! Initialise the collection for the assembly.
        call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

        ! -----
        ! Discretise the boundary conditions at the new point in time.
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call cc_assembleBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
        call cc_assembleFBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
        
        ! Implement the boundary conditions into the global solution vector.
        call sptivec_getTimestepData(rx, 1+isubstep, rtempVector)
        
        call vecfil_discreteBCsol (rtempVector)
        call vecfil_discreteFBCsol (rtempVector)
        
        call sptivec_setTimestepData(rx, 1+isubstep, rtempVector)
        
        ! Clean up the collection (as we are done with the assembly, that's it.
        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
        
      end do

    case DEFAULT
        
      call output_line ('Unsupported time discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tbc_implementBCsolution')
      call sys_halt()
    
    end select
    
    ! Release the temp vector
    call lsysbl_releaseVector (rtempVector)

    ! Release the BC's again.
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
    call bcasm_releaseDiscreteBC(rdiscreteBC)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementBCRHS (rproblem,rspaceTimeDiscr,rb,rtempVectorX)

!<description>
  ! Implements the boundary conditions of all timesteps into the RHS vector rb.
!</description>

!<input>
  ! Problem structure of the main problem.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  type(t_vectorBlock), intent(INOUT), optional :: rtempVectorX
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented
  ! to.
  type(t_spacetimeVector), intent(INOUT) :: rb
!</inputoutput>

!</subroutine>

    integer :: isubstep
    real(DP) :: dtstep,dtime
    type(t_vectorBlock) :: rtempVector
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata

    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    
    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)
    
    if (present(rtempVectorX)) then
      call lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    else
      ! Create a temp vector
      call lsysbl_createVecBlockByDiscr (&
          rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.true.)
    end if
    call lsysbl_assignDiscreteBC(rtempVector,rdiscreteBC)
    call lsysbl_assignDiscreteFBC(rtempVector,rdiscreteFBC)

    call lsyssc_getbase_double (rtempVector%RvectorBlock(1),p_Ddata)

    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep

    ! The implementation of the boundary conditions depends on the type
    ! of the time discretisation...
    select case (rspaceTimeDiscr%rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)

      ! Implement the bondary conditions into all initial solution vectors
      do isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = rproblem%rtimedependence%dtimeInit + isubstep*dtstep

        ! Initialise the collection for the assembly.
        call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

        ! -----
        ! Discretise the boundary conditions at the new point in time.
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call cc_assembleBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
        call cc_assembleFBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
        
        ! Implement the boundary conditions into the global solution vector.
        call sptivec_getTimestepData(rb, 1+isubstep, rtempVector)
        
        call vecfil_discreteBCrhs (rtempVector)
        call vecfil_discreteFBCrhs (rtempVector)
        
        call sptivec_setTimestepData(rb, 1+isubstep, rtempVector)

        ! Clean up the collection (as we are done with the assembly, that's it.
        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
        
      end do
      
    case (TDISCR_DG0)

      ! Implement the bondary conditions into all initial solution vectors
      do isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = rproblem%rtimedependence%dtimeInit + (real(isubstep,DP)+0.5_DP) * dtstep

        ! Initialise the collection for the assembly.
        call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

        ! -----
        ! Discretise the boundary conditions at the new point in time.
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call cc_assembleBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
        call cc_assembleFBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
        
        ! Implement the boundary conditions into the global solution vector.
        call sptivec_getTimestepData(rb, 1+isubstep, rtempVector)
        
        call vecfil_discreteBCrhs (rtempVector)
        call vecfil_discreteFBCrhs (rtempVector)
        
        call sptivec_setTimestepData(rb, 1+isubstep, rtempVector)
        
        ! Clean up the collection (as we are done with the assembly, that's it.
        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

      end do

    case default
        
      call output_line ('Unsupported time discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tbc_implementBCRHS')
      call sys_halt()
    
    end select

    ! Release the temp vector
    call lsysbl_releaseVector (rtempVector)

    ! Release the BC's again.
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
    call bcasm_releaseDiscreteBC(rdiscreteBC)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementBCdefect (rproblem,rspaceTimeDiscr,rd,rtempVectorX)

!<description>
  ! Implements the boundary conditions of all timesteps into the defect rd.
!</description>

!<input>
  ! Problem structure of the main problem.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented
  ! to.
  type(t_spacetimeVector), intent(INOUT) :: rd

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  type(t_vectorBlock), intent(INOUT), optional :: rtempVectorX
!</inputoutput>

!</subroutine>

    real(DP) :: dtstep,dtime
    integer :: isubstep
    type(t_vectorBlock) :: rtempVector
    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    
    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)

    ! Create a temp vector with our BC structures attached.
    if (present(rtempVectorX)) then
      call lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    else
      ! Create temp vectors
      call lsysbl_createVecBlockByDiscr (&
          rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.true.)
    end if
    call lsysbl_assignDiscreteBC(rtempVector,rdiscreteBC)
    call lsysbl_assignDiscreteFBC(rtempVector,rdiscreteFBC)

    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep
    
    ! The implementation of the boundary conditions depends on the type
    ! of the time discretisation...
    select case (rspaceTimeDiscr%rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)

      ! Implement the bondary conditions into all initial solution vectors
      do isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = rproblem%rtimedependence%dtimeInit + isubstep*dtstep

        ! Initialise the collection for the assembly.
        call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

        ! -----
        ! Discretise the boundary conditions at the new point in time.
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call cc_assembleBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
        call cc_assembleFBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
        
        ! Implement the boundary conditions into the global solution vector.
        call sptivec_getTimestepData(rd, 1+isubstep, rtempVector)
        
        call vecfil_discreteBCdef (rtempVector)
        call vecfil_discreteFBCdef (rtempVector)
        
        ! In the very first time step, we have the initial condition for the
        ! solution. The defect is =0 there!
!        IF (isubstep .EQ. 0) THEN
!          CALL lsyssc_clearVector (rtempVector%RvectorBlock(1))
!          CALL lsyssc_clearVector (rtempVector%RvectorBlock(2))
!          CALL lsyssc_clearVector (rtempVector%RvectorBlock(3))
!        END IF
        
        call sptivec_setTimestepData(rd, 1+isubstep, rtempVector)
        
        ! Clean up the collection (as we are done with the assembly, that's it.
        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
        
      end do
    
    case (TDISCR_DG0)


      ! Implement the bondary conditions into all initial solution vectors
      do isubstep = 0,rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = rproblem%rtimedependence%dtimeInit + (real(isubstep,DP)+0.5_DP) * dtstep

        ! Initialise the collection for the assembly.
        call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call cc_assembleBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
        call cc_assembleFBDconditions (rproblem,dtime,&
            rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
        
        ! Implement the boundary conditions into the global solution vector.
        call sptivec_getTimestepData(rd, 1+isubstep, rtempVector)
        
        call vecfil_discreteBCdef (rtempVector)
        call vecfil_discreteFBCdef (rtempVector)
        
        ! In the very first time step, we have the initial condition for the
        ! solution. The defect is =0 there!
!        IF (isubstep .EQ. 0) THEN
!          CALL lsyssc_clearVector (rtempVector%RvectorBlock(1))
!          CALL lsyssc_clearVector (rtempVector%RvectorBlock(2))
!          CALL lsyssc_clearVector (rtempVector%RvectorBlock(3))
!        END IF
        
        call sptivec_setTimestepData(rd, 1+isubstep, rtempVector)
        
        ! Clean up the collection (as we are done with the assembly, that's it.
        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

      end do

    case DEFAULT
        
      call output_line ('Unsupported time discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tbc_implementBCdefect')
      call sys_halt()
    
    end select

    call lsysbl_releaseVector(rtempVector)

    ! Release the BC's again.
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
    call bcasm_releaseDiscreteBC(rdiscreteBC)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_pressureToL20 (rproblem,rx,rtempVectorX)

!<description>
  ! Normalises the primal and dual pressure in all time steps where no Neumann
  ! boundary is present to have integral
  ! mean value = 0. This routine is typically used to filter an indefinite
  ! solution vector (e.g. in the pure-Dirichlet case).
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  type(t_problem), intent(INOUT), target :: rproblem
!</input>

!<inputoutput>
  ! A space-time vector where the pressure vectors whould be normalised.
  type(t_spacetimeVector), intent(INOUT) :: rx

  ! OPTIONAL: A spatial solution vector. If not specified, a vector
  ! is automatically created.
  type(t_vectorBlock), intent(INOUT), optional :: rtempVectorX
!</inputoutput>

!</subroutine>

    integer :: isubstep
    type(t_vectorBlock) :: rtempVector
    type(t_discreteBC) :: rdiscreteBC
    logical :: bhasNeumann
    real(dp) :: dtime,dtstep

    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    if (present(rtempVectorX)) then
      call lsysbl_duplicateVector (&
          rtempVectorX,rtempVector,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    else
      ! Create temp vectors
      call lsysbl_createVecBlockByDiscr (&
          rx%p_rblockDiscretisation,rtempVector,.true.)
    end if

    dtstep = rx%p_rtimeDiscretisation%dtstep

    ! Normalise the primal and dual pressure to zero.
    do isubstep = 0,rx%NEQtime-1

      ! Current point in time
      dtime = rproblem%rtimedependence%dtimeInit + isubstep*dtstep

      ! Initialise the collection for the assembly.
      call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

      ! Assemble the BC's.
      call bcasm_clearDiscreteBC(rdiscreteBC)
      call cc_assembleBDconditions (rproblem,&
          rx%p_rtimeDiscretisation%dtimeInit + isubstep*rx%p_rtimeDiscretisation%dtstep,&
          rx%p_rblockDiscretisation,&
          CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection,bhasNeumann)
    
      if (.not. bhasNeumann) then
        call sptivec_getTimestepData(rx, 1+isubstep, rtempVector)
        
        call vecfil_subvectorToL20 (rtempVectorX,3)
        call vecfil_subvectorToL20 (rtempVectorX,6)
        
        call sptivec_setTimestepData(rx, 1+isubstep, rtempVector)
      end if
      
      ! Clean up the collection (as we are done with the assembly, that's it.
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

    end do
  
    call lsysbl_releaseVector(rtempVector)

    ! Release the boundary conditions
    call bcasm_releaseDiscreteBC(rdiscreteBC)

  end subroutine

  ! *************************************************************************
  
!<subroutine>

  subroutine tbc_implementInitCondRHS (rproblem,rb,rinitCondRHS,rtempVectorD)

!<description>
  ! Implements the initial condition into the RHS vector rb.
  ! Overwrites the rb of the first time step.
  !
  ! Does not implement boundary conditions!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  type(t_problem), intent(INOUT), target :: rproblem
!</input>

!<inputoutput>
  ! A space-time vector with the RHS. The initial condition is implemented into
  ! this vector.
  type(t_spacetimeVector), intent(INOUT) :: rb

  ! A vector containing the data for the initial condition of the RHS.
  type(t_vectorBlock), intent(INOUT) :: rinitCondRHS

  ! A temporary vector in the size of a spatial vector.
  type(t_vectorBlock), intent(INOUT) :: rtempVectorD
!</inputoutput>

!</subroutine>

    ! Overwrite the primal RHS with the initial primal solution vector.
    ! This realises the inital condition.
    call sptivec_getTimestepData(rb, 1+0, rtempVectorD)
    call lsyssc_copyVector (rinitCondRHS%RvectorBlock(1),rtempVectorD%RvectorBlock(1))
    call lsyssc_copyVector (rinitCondRHS%RvectorBlock(2),rtempVectorD%RvectorBlock(2))
    call lsyssc_copyVector (rinitCondRHS%RvectorBlock(3),rtempVectorD%RvectorBlock(3))
    call sptivec_setTimestepData(rb, 1+0, rtempVectorD)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementInitCond (rproblem,rx,rinitCondSol,rtempVector)

!<description>
  ! Implements the initial condition into the vector rx.
  ! Overwrites the rx of the first time step.
  !
  ! Does not implement boundary conditions!
!</description>

!<input>
  ! A problem structure that provides information on all
  ! levels as well as temporary vectors.
  type(t_problem), intent(INOUT), target :: rproblem
!</input>

!<inputoutput>
  ! A space-time vector with the RHS. The initial condition is implemented into
  ! this vector.
  type(t_spacetimeVector), intent(INOUT) :: rx

  ! A vector containing the data for the initial condition of the RHS.
  type(t_vectorBlock), intent(INOUT) :: rinitCondSol

  ! A temporary vector in the size of a spatial vector.
  type(t_vectorBlock), intent(INOUT) :: rtempVector
!</inputoutput>

!</subroutine>

    ! Overwrite the primal solution with the initial primal solution vector.
    ! This realises the inital condition.

!    CALL sptivec_getTimestepData(rx, 1+0, rtempVector)
!    CALL lsyssc_copyVector (rinitCondSol%RvectorBlock(1),rtempVector%RvectorBlock(1))
!    CALL lsyssc_copyVector (rinitCondSol%RvectorBlock(2),rtempVector%RvectorBlock(2))
!    CALL sptivec_setTimestepData(rx, 1+0, rtempVector)
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementInitCondDefect (rspaceTimeDiscr, rd, rtempvectorD)

!<description>
  ! Implements the initial and terminal condition into a defect vector rd.
  ! Overwrites the rd of the first time step.
  !
  ! Does not implement boundary conditions!
!</description>

!<inputoutput>
  ! Discretisation structure that corresponds to rx.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr

  ! A space-time vector containing the defect in the first subvector.
  type(t_spacetimeVector), intent(INOUT) :: rd

  ! A temporary vector in the size of a spatial vector.
  type(t_vectorBlock), intent(INOUT) :: rtempVectorD
!</inputoutput>

!</subroutine>

    real(DP), dimension(:),pointer :: p_Db
    
    ! DEBUG!!!    
    call lsysbl_getbase_double (rtempVectorD,p_Db)

    ! Overwrite the primal defect with 0 -- as the solution must not be changed.
    ! This realises the inital condition.
    call sptivec_getTimestepData(rd, 1+0, rtempVectorD)
    call tbc_implementInitCondDefSingle (rspaceTimeDiscr, rtempVectorD)
    call sptivec_setTimestepData(rd, 1+0, rtempVectorD)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementInitCondDefSingle (rspaceTimeDiscr, rd)

!<description>
  ! Implements the initial condition into a defect vector rd,
  ! representing the defect in the first timestep.
  !
  ! Does not implement boundary conditions!
!</description>

!<inputoutput>
  ! Discretisation structure that corresponds to rx.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr

  ! A vector containing the defect in the first subvector.
  type(t_vectorBlock), intent(INOUT) :: rd
!</inputoutput>

!</subroutine>

    real(DP), dimension(:),pointer :: p_Db
    
    ! DEBUG!!!    
    call lsysbl_getbase_double (rd,p_Db)
    
!    CALL lsyssc_clearVector(rd%RvectorBlock(1))
!    CALL lsyssc_clearVector(rd%RvectorBlock(2))
!    CALL lsyssc_clearVector(rd%RvectorBlock(3))

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementTermCondDefSingle (rspaceTimeDiscr, rd)

!<description>
  ! Implements the terminal condition into a defect vector rd,
  ! representing the defect in the last timestep.
  !
  ! Does not implement boundary conditions!
!</description>

!<inputoutput>
  ! Discretisation structure that corresponds to rx.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr

  ! A vector containing the defect in the last subvector.
  type(t_vectorBlock), intent(INOUT) :: rd
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

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine tbc_implementSpatialBCtoRHS (rproblem, dtime, rvector)
  
!<description>
  ! Implements the spatial boundary conditions into the spatial RHS vector
  ! rvector.
!</description>
 
!<input>
  ! Time where the BC's should be implemented.
  ! Must not necessarily coincide with the start/end time of the timestep.
  real(DP), intent(IN) :: dtime
!</input>  

!<inputoutput>  
  ! Problem structure.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Source and destination RHS vector
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
  
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata
  
    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rvector,p_Ddata)

    ! Initialise the collection for the assembly.
    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)

    ! Assemble the BC's.
    call cc_assembleBDconditions (rproblem,dtime,&
        rvector%p_rblockDiscr,&
        CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
    call cc_assembleFBDconditions (rproblem,dtime,&
        rvector%p_rblockDiscr,&
        CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)

    ! Implement the boundary conditions into the RHS.
    ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
    ! as Dirichlet values mustn't be multiplied with GAMMA!
    call vecfil_discreteBCrhs (rvector,rdiscreteBC)
    call vecfil_discreteFBCrhs (rvector,rdiscreteFBC)

    ! Release the BC's again.
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
    call bcasm_releaseDiscreteBC(rdiscreteBC)
    
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine tbc_implementSpatialBCdefect (rproblem,dtime,rspaceTimeDiscr,rd)

!<description>
  ! Implements the boundary conditions at time dtime into the defect rd.
!</description>

!<input>
  ! Problem structure of the main problem.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Discretisation structure that corresponds to rx.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr

  ! Time where the BC's should be implemented.
  ! Must not necessarily coincide with the start/end time of the timestep.
  real(DP), intent(IN) :: dtime
!</input>

!<inputoutput>
  ! A space-time vector with the solution where the BC's should be implemented to.
  type(t_vectorBlock), intent(INOUT) :: rd
!</inputoutput>

!</subroutine>

  
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata
  
    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rd,p_Ddata)

    ! Initialise the collection for the assembly.
    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)

    ! Assemble the BC's.
    call cc_assembleBDconditions (rproblem,dtime,&
        rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
    call cc_assembleFBDconditions (rproblem,dtime,&
        rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)

    ! Implement the boundary conditions into the RHS.
    ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
    ! as Dirichlet values mustn't be multiplied with GAMMA!
    call vecfil_discreteBCdef (rd,rdiscreteBC)
    call vecfil_discreteFBCdef (rd,rdiscreteFBC)

    ! Release the BC's again.
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
    call bcasm_releaseDiscreteBC(rdiscreteBC)
    
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

end module
