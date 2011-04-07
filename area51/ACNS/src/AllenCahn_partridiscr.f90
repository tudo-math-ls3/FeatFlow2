!##############################################################################
!# ****************************************************************************
!# <name> AllenCahn_partridiscr </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to read the parametrisation, create the
!# triangulation and set up the discretisation for the heat conduction problem.
!# The following routines can be found here:
!#
!# 0.) AC_initSolution
!# 1.) AC_initParamTriang
!#     -> Read .PRM/.TRI files. Generate meshes on all levels.
!#
!# 2.) AC_doneParamTriang
!#     Clean up parametrisation/triangulation, release memory.
!#
!# 3.) AC_initDiscretisation
!#     -> Initialise the spatial discretisation.
!#
!# 4.) AC_doneDiscretisation
!#     -> Clean up the discretisation, release memory.
!# </purpose>
!##############################################################################

MODULE AllenCahn_partridiscr


  
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
  USE sortstrategy
  USE coarsegridcorrection
  USE ucd
  USE timestepping
  USE genoutput
  
  use mprimitives

  USE collection
  USE paramlist
  use convection
  use vectorio    


  USE AllenCahn_callback
  
  USE AllenCahn_basic
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>
  subroutine AC_initSolution(rproblem, rvector)
!<description>
  ! Initialises the initial solution vector into rvector. Depending on the settings
  ! in the DAT file this is either zero or read from a file.
  !
  ! The routine assumes that basic mass matrices have already been constructed.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), intent(INOUT) :: rproblem
!</input>

!<inputoutput>
  ! The solution vector to be initialised. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: istart,ctypeInitialSolution
    type(t_vectorBlock) :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sarray,sfile,sfileString
    integer :: ilev,ierror
    integer(PREC_VECIDX) :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection 
    type(t_linearForm) :: rlinform
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner
    type(t_matrixBlock), dimension(1) :: Rmatrices
!    type(t_matrixBlock) :: Rmatrices
     
    type(t_vectorBlock) :: rsingleRHS,rsingleSol
    type(t_vectorscalar) :: rtemp

!     logical :: mcai
    ! Get the parameter what to do with rvector

    ctypeInitialSolution=3
    
	!istart 
    !sfileString

    ! What is the type of the initial solution?
    select case (ctypeInitialSolution)
    case (0)
      ! Init with zero
      call lsysbl_clearVector (rvector)
      
    case (1,2)
      ! We have to read a file -- formatted or unformatted.
      !
      ! Create a temp vector at level NLMAX-istart+1.
!      if (istart .gt. 0) then
!        ilev = istart
!      else
!        ilev = rproblem%NLMAX-abs(istart)+1
!      end if
      
!      if (ilev .lt. rproblem%NLMIN) then
!        call output_line (&
!            'Level of start vector is < NLMIN! Initialising with zero!', &
!            OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
!        istart = 0
!      end if

!      if (ilev .gt. rproblem%NLMAX) then
!        call output_line (&
!            'Level of start vector is > NLMAX! Initialising with zero!', &
!            OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
!        istart = 0
!      end if
      
      ! Remove possible ''-characters
!      read(sfileString,*) sfile

!      call lsysbl_createVectorBlock (&
!          rproblem%RlevelInfo(ilev)%rdiscretisation,rvector1,.false.)
      
      ! Read in the vector
!      call vecio_readBlockVectorHR (&
!          rvector1, sarray, .true., 0, sfile, istart .gt. 0)
          
      ! If the vector is on level < NLMAX, we have to bring it to level NLMAX
!      do while (ilev .lt. rproblem%NLMAX)
        
        ! Initialise a vector for the higher level and a prolongation structure.
!        call lsysbl_createVectorBlock (&
!            rproblem%RlevelInfo(ilev+1)%rdiscretisation,rvector2,.false.)
        
!        call mlprj_initProjectionVec (rprojection,rvector2)
        
        ! Prolongate to the next higher level.

!        NEQ = mlprj_getTempMemoryVec (rprojection,rvector1,rvector2)
!        if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
!        call mlprj_performProlongation (rprojection,rvector1, &
!                                        rvector2,rvectorTemp)
!        if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)
        
        ! Swap rvector1 and rvector2. Release the coarse grid vector.
!        call lsysbl_swapVectors (rvector1,rvector2)
!        call lsysbl_releaseVector (rvector2)
        
!        call mlprj_doneProjection (rprojection)
        
        ! rvector1 is now on level ilev+1
!        ilev = ilev+1
        
!      end do
      
      ! Copy the resulting vector rvector1 to the output.
!      call lsysbl_copyVector (rvector1,rvector)
      
      ! Release the temp vector
!      call lsysbl_releaseVector (rvector1)

    case (3)
    
      ! We have to create the solution by analytical callback functions.
      ! To do this for an arbitrary finite element, we have to do an
      ! L2 projection for X,Y and P. That means we have to solve:
      !     <u,phi> = <u_analytic,phi>
      ! which means in the FE context:
      !     Mu = f
      ! with f a RHS created by analytical point values.

!      if (rproblem%MSHOW_Initialisation .ge. 1) &
      call output_line('AC problem, Preparing L2 projection of analytical initial solution...')

      ! Get a pointer to the RHS on the finest level as well as to the
      ! block discretisation structure:
      p_rdiscretisation => rvector%p_rblockDiscr
      
      ! The vector structure is already prepared, but the entries are missing.
      !
      ! At first set up the corresponding linear form (u_analytical,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! Create a temp vector as RHS.
      
	  !call lsysbl_createVectorBlock (&
      !    rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,rvector1,.false.)
	  ! if rdiscretisation is not a field of t_ACproblem_lvl, we use
      call lsysbl_createVecBlockIndirect (rvector,rvector1,.TRUE., ST_double)
      
      ! Assemble the RHS.
      !
      ! Initialise the collection for the assembly process with callback routines.
      ! Basically, this stores the simulation time in the collection if the
      ! simulation is nonstationary.
      call AC_initCollectForAssembly (rproblem,rproblem%rcollection)

      ! Discretise the X-velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
                rvector1%RvectorBlock(1),coeff_AnalyticSolution_AC,&
                rproblem%rcollection)
                                  
      ! Clean up the collection (as we are done with the assembly.
      call AC_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
      ! Now prepare a linear solver for solving with the mass matrix.
      ! This is just a simple 1-level Jacobi solver as the mass matrix
      ! is usually well conditioned.
      call linsol_initJacobi (p_rpreconditioner, 0.8_DP)
      call linsol_initDefCorr (p_rsolverNode, p_rpreconditioner)
      
      p_rsolverNode%depsRel = SYS_EPSREAL_DP * 100.0_DP
      p_rsolverNode%nmaxiterations = 1000
      p_rsolverNode%ioutputLevel = 1

      call output_line('AC problem, Solving L2 projection for \phi...')
      
      ! We solve separately for the three components.
      ! Prepare the solver for the velocity.
      ! Prepare the solver for the X-velocity.

! We should make sure that rmatrixMass has been assembled. 


     call lsysbl_createMatFromScalar(&
        rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixMass%RmatrixBlock(1,1),Rmatrices(1))


      call lsyssc_copyvector(rvector%RvectorBlock(1), rtemp)	
      call lsyssc_scalevector(rtemp, 0.0_DP)
!      call lsyssc_synchroniseSortMatVec (Rmatrices(1)%RmatrixBlock(1,1), &
!	                                     rvector%RvectorBlock(1), rtemp)


      call linsol_setMatrices (p_rsolverNode,Rmatrices)
      call linsol_initStructure (p_rsolverNode,ierror)
      call linsol_initData (p_rsolverNode,ierror)

      ! -----
      ! Solve for the X-velocity.
      call lsysbl_createVecFromScalar(rvector1%RvectorBlock(1),rsingleRHS)
      call lsysbl_createVecFromScalar(rvector%RvectorBlock(1),rsingleSol)
      call lsysbl_duplicateVector (rsingleSol,rvector2,&
          LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
     
      call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
!~~~~~We should copy back rsingleSol to rvecto%RvectorBlock(1)
 
      if (p_rsolverNode%iresult .ne. 0) then
        ! Cancel, there is something wrong.
        call output_line('Cannot compute L2 projection, solver broke down. Using zero!')
        call lsysbl_clearVector (rvector)

!      else 
!        call lsysbl_copyVector(rsingleSol, rvector)
      end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! -----
      ! That's it, cleanup.
      call lsysbl_releaseVector (rsingleRHS)
      call lsysbl_releaseVector (rsingleSol)
      call lsysbl_releaseVector (rvector2)
      call lsysbl_releaseMatrix (Rmatrices(1))

	  call lsyssc_releaseVector (rtemp)
      
      call linsol_doneData (p_rsolverNode,ierror)
      call linsol_doneStructure (p_rsolverNode,ierror)
      call linsol_releaseSolver (p_rsolverNode)
    
    end select        

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine AC_initCollectForAssembly (rproblem,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the CC2D framework
  ! and has usually not to be changed by the user.
  !
  ! The subroutine prepares the collection rcollection to be passed to callback
  ! routines for assembling boundary conditions or RHS vectors. It is
  ! called directly prior to the assembly to store problem-specific information
  ! in the quick-access arrays of the collection.
  ! Basically speaking, the routine stores information about whether thw problem
  ! is stationary, nonstationary or about the current simulation time.
!</description>

!<input>
  ! Problem structure with all problem relevant data.
  type(t_ACproblem), intent(IN) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be initialised.
  type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! In a nonstationary simulation, save the simulation time as well as the
    ! minimum and maximum time to the quick-access array of the collection,
    ! so it can be accessed in the callback routines!
!    rcollection%Iquickaccess(1) = rproblem%itimedependence
    rcollection%Dquickaccess(1) = rproblem%rtimedependence%dtime
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>

  subroutine AC_doneCollectForAssembly (rproblem,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the CC2D framework
  ! and has usually not to be changed by the user.
  !
  ! After the assembly process, this subroutine is called to release temporary
  ! information from the collection which was stored there by 
  ! cc_initCollectForAssembly.
!</description>
  
!<input>
  ! Problem structure with all problem relevant data.
  type(t_ACproblem), intent(IN) :: rproblem
!</input>

!<inputoutput>
  ! Collection structure to be cleaned up.
  type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! Currently, this subroutine is empty as all information stored in
    ! the collection in cc_initCollectForAssembly is put to the quick-access
    ! arrays -- which do not have to be cleaned up. 
    ! This might change in future...

  end subroutine
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE AC_initParamTriang (NLMIN,NLMAX,rACproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Minimum refinement level of the mesh; = coarse grid = level 1
  INTEGER, INTENT(IN) :: NLMIN
  
  ! Maximum refinement level
  INTEGER, INTENT(IN) :: NLMAX
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_ACproblem), INTENT(INOUT) :: rACproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  
    ! Initialise the level in the problem structure

    rACproblem%NLMIN = NLMIN
    rACproblem%NLMAX = NLMAX

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL boundary_read_prm(rACproblem%rboundary, './pre/QUAD.prm')
        

    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rACproblem%RlevelInfo(rACproblem%NLMIN)%rtriangulation, &
        './pre/QUAD.tri', rACproblem%rboundary)

    ! Refine the mesh up to the minimum level
    CALL tria_quickRefine2LevelOrdering(rACproblem%NLMIN-1,&
        rACproblem%RlevelInfo(rACproblem%NLMIN)%rtriangulation,rACproblem%rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    CALL tria_initStandardMeshFromRaw (&
        rACproblem%RlevelInfo(rACproblem%NLMIN)%rtriangulation,rACproblem%rboundary)
    
    ! Now, refine to level up to nlmax.
    DO i=rACproblem%NLMIN+1,rACproblem%NLMAX
      CALL tria_refine2LevelOrdering (rACproblem%RlevelInfo(i-1)%rtriangulation,&
          rACproblem%RlevelInfo(i)%rtriangulation, rACproblem%rboundary)
      CALL tria_initStandardMeshFromRaw (rACproblem%RlevelInfo(i)%rtriangulation,&
          rACproblem%rboundary)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE AC_initDiscretisation (rACproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: I
  
    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! An object for the block discretisation on one level
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    DO i=rACproblem%NLMIN,rACproblem%NLMAX
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rACproblem%rboundary
      p_rtriangulation => rACproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      ALLOCATE(p_rdiscretisation)
      CALL spdiscr_initBlockDiscr (p_rdiscretisation,1,&
                                    p_rtriangulation, p_rboundary)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rACproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      ! p_rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:

!       CALL spdiscr_initDiscr_simple ( &
!                   p_rdiscretisation%RspatialDiscr(1), &
!                   EL_E011,CUB_G4X4, &
!                   p_rtriangulation, p_rboundary)

      CALL spdiscr_initDiscr_simple ( &
                  p_rdiscretisation%RspatialDiscr(1), &
                  EL_Q1,CUB_G4X4, &
                  p_rtriangulation, p_rboundary)

!     CALL spdiscr_initDiscr_simple ( &
!                  p_rdiscretisation%RspatialDiscr(1), &
!                  EL_EM30,CUB_G4X4, &
!                  p_rtriangulation, p_rboundary)


    END DO
                                   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE AC_doneDiscretisation (rACproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rACproblem%NLMAX,rACproblem%NLMIN,-1
      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations....
      CALL spdiscr_releaseBlockDiscr(rACproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! and remove the allocated block discretisation structure from the heap.
      DEALLOCATE(rACproblem%RlevelInfo(i)%p_rdiscretisation)
    END DO
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE AC_doneParamTriang (rACproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rACproblem%NLMAX,rACproblem%NLMIN,-1
      ! Release the triangulation
      CALL tria_done (rACproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! Finally release the domain.
    CALL boundary_release (rACproblem%rboundary)
    
  END SUBROUTINE

END MODULE
