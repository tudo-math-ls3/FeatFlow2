!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2timesupersystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2timesupersystem

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
  USE cc2dmediumm2stationary
  USE adaptivetimestep
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  
  USE spacetimevectors
  USE dofmapping
  
  USE matrixio
    
  IMPLICIT NONE

!<types>

!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  TYPE t_ccoptSpaceTimeDiscretisation
  
    ! Spatial refinement level of this matrix.
    INTEGER :: NLMAX
  
    ! Number of time steps
    INTEGER :: niterations         = 0

    ! Absolute start time of the simulation
    REAL(DP) :: dtimeInit          = 0.0_DP     
    
    ! Maximum time of the simulation
    REAL(DP) :: dtimeMax           = 1.0_DP
    
    ! Time-stepping scheme (IFRSTP);
    ! 0=one step scheme
    INTEGER :: ctimeStepScheme     = 0
    
    ! Parameter for one step scheme (THETA) if itimeStepScheme=0;
    ! =0:Forward Euler(instable), =1: Backward Euler, =0.5: Crank-Nicolson
    REAL(DP) :: dtimeStepTheta     = 1.0_DP
    
    ! Time step length of the time discretisation
    REAL(DP) :: dtstep

    ! Regularisation parameter for the control $\alpha$. Must be <> 0.
    ! A value of 0.0 disables the terminal condition.
    REAL(DP) :: dalphaC = 1.0_DP
    
    ! Regularisation parameter for the terminal condition 
    ! $\gamma/2*||y(T)-z(T)||$.
    ! A value of 0.0 disables the terminal condition.
    REAL(DP) :: dgammaC = 0.0_DP

    ! Problem-related structure that provides the templates for
    ! matrices/vectors on the level of the matrix.
    TYPE(t_problem_lvl), POINTER :: p_rlevelInfo
    
  END TYPE

!</types>

CONTAINS

  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_initParamsSupersystem (rproblem,ilevelTime,ilevelSpace,&
      rspaceTimeDiscr, rx, rb, rd)
  
!<description>
  ! Initialises the time stepping scheme according to the parameters in the
  ! DAT file and the problem structure.
!</description>

!<input>
  ! The problem structure describing the whole problem.
  TYPE(t_problem), INTENT(IN), TARGET :: rproblem
  
  ! 'Refinement level in time'. =1: Initialise as described in
  ! rproblem. >1: Refine (irefineInTime-1) times regularly in time.
  ! (-> #timesteps * 2**(irefineInTime-1) )
  INTEGER, INTENT(IN) :: ilevelTime

  ! 'Refinement level in space'. =1: Calculate on the coarse mesh
  ! >0: calculate on level ilevelSpace. Must be <= rproblem%NLMAX,
  ! >= rproblem%NLMIN.
  INTEGER, INTENT(IN) :: ilevelSpace
!</input>

!<inputoutput>
  ! Supersystem-structure to be initialised.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(OUT) :: rspaceTimeDiscr

  ! A space-time vector that is initialised for the current solution.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A space-time vector that is initialised for the current RHS.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb

  ! A space-time vector that is initialised for the current defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP) :: dtstep

    ! Copy most relevant data from the problem structure.
    rspaceTimeDiscr%NLMAX           = ilevelSpace
    rspaceTimeDiscr%p_rlevelInfo    => rproblem%RlevelInfo(ilevelSpace)
    rspaceTimeDiscr%niterations     = &
        rproblem%rtimedependence%niterations * 2**MIN(0,ilevelTime-1)
    rspaceTimeDiscr%dtimeInit       = rproblem%rtimedependence%dtimeInit
    rspaceTimeDiscr%dtimeMax        = rproblem%rtimedependence%dtimeMax
    rspaceTimeDiscr%ctimeStepScheme = rproblem%rtimedependence%ctimeStepScheme
    rspaceTimeDiscr%dtimeStepTheta  = rproblem%rtimedependence%dtimeStepTheta
    
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',rspaceTimeDiscr%dalphaC,1.0_DP)
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dgammaC',rspaceTimeDiscr%dgammaC,0.0_DP)

    ! The complete time interval is divided into niteration iterations
    dtstep = (rspaceTimeDiscr%dtimemax-rspaceTimeDiscr%dtimeInit) &
             / REAL(MAX(rspaceTimeDiscr%niterations,1),DP)
             
    rspaceTimeDiscr%dtstep = dtstep
    
    ! Initialise the global solution- and defect- vector.
    CALL sptivec_initVector (rx,&
        dof_igetNDofGlobBlock(rproblem%RlevelInfo(ilevelSpace)%p_rdiscretisation,.FALSE.),&
        rspaceTimeDiscr%niterations)
    CALL sptivec_initVector (rb,&
        dof_igetNDofGlobBlock(rproblem%RlevelInfo(ilevelSpace)%p_rdiscretisation,.FALSE.),&
        rspaceTimeDiscr%niterations)
    CALL sptivec_initVector (rd,&
        dof_igetNDofGlobBlock(rproblem%RlevelInfo(ilevelSpace)%p_rdiscretisation,.FALSE.),&
        rspaceTimeDiscr%niterations)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_doneParamsSupersystem (rspaceTimeDiscr,rx,rb,rd)
  
!<description>
  ! Cleans up a given supersystem structure.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(OUT) :: rspaceTimeDiscr

  ! A space-time vector that is initialised for the current solution.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A space-time vector that is initialised for the current RHS.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb

  ! A space-time vector that is initialised for the current defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! Release memory.
    CALL sptivec_releaseVector (rb)
    CALL sptivec_releaseVector (rd)
    CALL sptivec_releaseVector (rx)

  END SUBROUTINE

!  ! ***************************************************************************
!  
!!<subroutine>
!
!  SUBROUTINE c2d2_solveSupersysDirect (rproblem, rspaceTimeDiscr, rx, rd, &
!      rtempvectorX, rtempvectorB, rtempvectorD)
!
!!<description>
!  ! This routine assembles and solves the time-space coupled supersystem:
!  ! $Ax=b$. The RHS vector is generated on-the-fly.
!  ! The routine generates the full matrix in memory and solves with UMFPACK, 
!  ! so it should only be used for debugging!
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!
!  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
!  ! coupled space-time matrix.
!  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!!</input>
!
!!<inputoutput>
!  ! A space-time vector defining the current solution.
!  ! Is replaced by the new solution
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx
!
!  ! A temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX
!
!  ! A second temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorB
!
!  ! A third temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD
!
!  ! A space-time vector that receives the defect.
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    INTEGER :: isubstep,ilevel,ierror,i
!    TYPE(t_matrixBlock) :: rblockTemp
!    TYPE(t_matrixScalar) :: rmassLumped
!    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
!    TYPE(t_vectorBlock) :: rxGlobal, rbGlobal, rdGlobal
!    TYPE(t_vectorBlock) :: rxGlobalSolve, rbGlobalSolve, rdGlobalSolve
!    TYPE(t_matrixBlock) :: rglobalA
!    TYPE(t_linsolNode), POINTER :: rsolverNode
!    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
!    INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: Isize
!    INTEGER(PREC_VECIDX), DIMENSION(6) :: Isize2
!    
!    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd
!    
!    ! If the following constant is set from 1.0 to 0.0, the primal system is
!    ! decoupled from the dual system!
!    REAL(DP), PARAMETER :: dprimalDualCoupling = 0.0 !1.0
!
!    ilevel = rspaceTimeDiscr%NLMAX
!    
!    ! Calculate the lumped mass matrix of the FE space -- we need it later!
!    CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass,&
!        rmassLumped, LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
!    CALL lsyssc_lumpMatrixScalar (rmassLumped,LSYSSC_LUMP_STD)
!    
!    ! ----------------------------------------------------------------------
!    ! 1.) Generate the global RHS vector
!    
!    CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
!    CALL lsysbl_getbase_double (rtempVectorB,p_Db)
!    CALL lsysbl_getbase_double (rtempVectorD,p_Dd)
!
!    DO isubstep = 0,rspaceTimeDiscr%niterations
!    
!      ! Current point in time
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%dtstep
!          
!      ! Generate the RHS of that point in time.
!      CALL c2d2_generateBasicRHS (rproblem,rtempVectorB)
!      
!      ! Multiply the RHS of the dual velocity by dtstep according to the time
!      ! discretisation -- except for if we are in the last timestep!
!      !
!      ! In the last timestep, if GAMMA is =0, we have no terminal condition
!      ! and thus have to force the RHS to 0!
!      ! Otherwise, the terminat condition is multiplied with dgammaC.
!      IF (isubstep .NE. rspaceTimeDiscr%niterations) THEN
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(4),rspaceTimeDiscr%dtstep)
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(5),rspaceTimeDiscr%dtstep)
!      ELSE
!        ! Multiply -z by gamma, that's it.
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(4),rspaceTimeDiscr%dgammaC)
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(5),rspaceTimeDiscr%dgammaC)
!      END IF
!      
!      ! Initialise the collection for the assembly process with callback routines.
!      ! Basically, this stores the simulation time in the collection if the
!      ! simulation is nonstationary.
!      CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)
!
!      ! Discretise the boundary conditions at the new point in time -- 
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
!      END IF
!
!      ! Implement the boundary conditions into the RHS.
!      ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
!      ! as Dirichlet values mustn't be multiplied with GAMMA!
!      CALL vecfil_discreteBCsol (rtempVectorB)
!      CALL vecfil_discreteFBCsol (rtempVectorB)      
!      
!      CALL sptivec_setTimestepData(rd, isubstep, rtempVectorB)
!      
!      ! Clean up the collection (as we are done with the assembly, that's it.
!      CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)
!
!    END DO
!
!    ! Release the mass matr, we don't need it anymore
!    CALL lsyssc_releaseMatrix (rmassLumped)
!    
!    ! ----------------------------------------------------------------------
!    ! 2.) Generate the matrix A
!    !
!    ! Create a global matrix:
!    CALL lsysbl_createEmptyMatrix (rglobalA,6*(rspaceTimeDiscr%niterations+1))
!    
!    ! Loop through the substeps
!    
!    DO isubstep = 0,rspaceTimeDiscr%niterations
!    
!      ! Current point in time
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%dtstep
!
!      ! -----
!      ! Discretise the boundary conditions at the new point in time -- 
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
!      END IF
!      
!      ! The first and last substep is a little bit special concerning
!      ! the matrix!
!      IF (isubstep .EQ. 0) THEN
!        
!        ! We are in the first substep
!      
!        ! -----
!      
!        ! Create a matrix that applies "-M" to the dual velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(4,4)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(5,5)%dscaleFactor = -1.0_DP
!        
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!        
!        ! Include "-M" in the global matrix at position (1,2).
!        CALL insertMatrix (rblockTemp,rglobalA,7,1)
!        
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!      
!        ! -----
!        
!        ! Now the hardest -- or longest -- part: The diagonal matrix.
!        !
!        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
!        ! Will be modified by c2d2_assembleLinearisedMatrices later.
!        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!      
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!      
!        CALL c2d2_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
!            rnonlinearIterationTmp,'CC2D-NONLINEAR')
!
!        ! Set up all the weights in the core equation according to the current timestep.
!        rnonlinearIterationTmp%diota1 = 1.0_DP
!        rnonlinearIterationTmp%diota2 = 0.0_DP
!
!        rnonlinearIterationTmp%dkappa1 = 1.0_DP
!        rnonlinearIterationTmp%dkappa2 = 0.0_DP
!        
!        rnonlinearIterationTmp%dalpha1 = 0.0_DP
!        rnonlinearIterationTmp%dalpha2 = 1.0_DP
!        
!        rnonlinearIterationTmp%dtheta1 = 0.0_DP
!        rnonlinearIterationTmp%dtheta2 = rspaceTimeDiscr%dtstep
!        
!        rnonlinearIterationTmp%dgamma1 = 0.0_DP
!        rnonlinearIterationTmp%dgamma2 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        
!        rnonlinearIterationTmp%deta1 = 0.0_DP
!        rnonlinearIterationTmp%deta2 = rspaceTimeDiscr%dtstep
!        
!        rnonlinearIterationTmp%dtau1 = 0.0_DP
!        rnonlinearIterationTmp%dtau2 = 1.0_DP
!        
!        rnonlinearIterationTmp%dmu1 = 0.0_DP
!        rnonlinearIterationTmp%dmu2 = -rspaceTimeDiscr%dtstep
!        
!        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
!        ! Include the boundary conditions into the matrices.
!        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
!            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
!            
!        ! Insert the system matrix for the dual equation to our global matrix.
!        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
!            rglobalA,1,1)
!            
!      ELSE IF (isubstep .EQ. rspaceTimeDiscr%niterations) THEN
!        
!        ! We are in the last substep
!        
!        ! -----
!        
!        ! Create a matrix that applies "-M" to the primal velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(1,1)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(2,2)%dscaleFactor = -1.0_DP
!        
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!
!        ! Include "-M" in the global matrix at position (2,1).
!        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
!        
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!      
!        ! -----
!        
!        ! Now the hardest -- or longest -- part: The diagonal matrix.
!        !
!        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
!        ! Will be modified by c2d2_assembleLinearisedMatrices later.
!        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!      
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!      
!        CALL c2d2_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
!            rnonlinearIterationTmp,'CC2D-NONLINEAR')
!
!        ! Set up all the weights in the core equation according to the current timestep.
!        rnonlinearIterationTmp%diota1 = 0.0_DP
!        rnonlinearIterationTmp%diota2 = 0.0_DP
!
!        rnonlinearIterationTmp%dkappa1 = 0.0_DP
!        rnonlinearIterationTmp%dkappa2 = 1.0_DP
!        
!        rnonlinearIterationTmp%dalpha1 = 1.0_DP
!        rnonlinearIterationTmp%dalpha2 = 1.0_DP
!        
!        rnonlinearIterationTmp%dtheta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%dtheta2 = 0.0_DP
!        
!        rnonlinearIterationTmp%dgamma1 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        rnonlinearIterationTmp%dgamma2 = 0.0_DP
!        
!        rnonlinearIterationTmp%deta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%deta2 = 0.0_DP
!        
!        rnonlinearIterationTmp%dtau1 = 1.0_DP
!        rnonlinearIterationTmp%dtau2 = 0.0_DP
!        
!        rnonlinearIterationTmp%dmu1 = dprimalDualCoupling * &
!            rspaceTimeDiscr%dtstep / rspaceTimeDiscr%dalphaC
!        rnonlinearIterationTmp%dmu2 = -rspaceTimeDiscr%dgammaC
!        
!        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
!        ! Include the boundary conditions into the matrices.
!        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
!            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
!        
!        ! Insert the system matrix for the dual equation to our global matrix.
!        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
!            rglobalA,isubstep*6+1,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!      
!      ELSE
!      
!        ! We are sonewhere in the middle of the matrix. There is a substep
!        ! isubstep+1 and a substep isubstep-1!
!        
!        ! -----
!        
!        ! Create a matrix that applies "-M" to the dual velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(4,4)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(5,5)%dscaleFactor = -1.0_DP
!        
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!
!        ! Include that in the global matrix below the diagonal
!        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+7,isubstep*6+1)
!        
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!      
!        ! -----
!        
!        ! Create a matrix that applies "-M" to the primal velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(1,1)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(2,2)%dscaleFactor = -1.0_DP
!        
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!
!        ! Include that in the global matrix above the diagonal
!        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
!        
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----      
!
!        ! Now the hardest -- or longest -- part: The diagonal matrix.
!        !
!        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
!        ! Will be modified by c2d2_assembleLinearisedMatrices later.
!        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!      
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!      
!        CALL c2d2_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
!            rnonlinearIterationTmp,'CC2D-NONLINEAR')
!
!        ! Set up all the weights in the core equation according to the current timestep.
!        rnonlinearIterationTmp%diota1 = 0.0_DP
!        rnonlinearIterationTmp%diota2 = 0.0_DP
!
!        rnonlinearIterationTmp%dkappa1 = 0.0_DP
!        rnonlinearIterationTmp%dkappa2 = 0.0_DP
!        
!        rnonlinearIterationTmp%dalpha1 = 1.0_DP
!        rnonlinearIterationTmp%dalpha2 = 1.0_DP
!        
!        rnonlinearIterationTmp%dtheta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%dtheta2 = rspaceTimeDiscr%dtstep
!        
!        rnonlinearIterationTmp%dgamma1 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        rnonlinearIterationTmp%dgamma2 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        
!        rnonlinearIterationTmp%deta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%deta2 = rspaceTimeDiscr%dtstep
!        
!        rnonlinearIterationTmp%dtau1 = 1.0_DP
!        rnonlinearIterationTmp%dtau2 = 1.0_DP
!        
!        rnonlinearIterationTmp%dmu1 = dprimalDualCoupling * &
!            rspaceTimeDiscr%dtstep / rspaceTimeDiscr%dalphaC
!        rnonlinearIterationTmp%dmu2 = -rspaceTimeDiscr%dtstep
!      
!        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
!        ! Include the boundary conditions into the matrices.
!        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
!            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
!        
!        ! Insert the system matrix for the dual equation to our global matrix.
!        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
!            rglobalA,isubstep*6+1,isubstep*6+1)
!            
!      END IF
!    
!    END DO
!    
!    ! Update structural information of the global matrix.
!    CALL lsysbl_updateMatStrucInfo (rglobalA)
!
!    ! Write the global matrix to a file.
!    !CALL matio_writeBlockMatrixHR(rglobalA,'MATRIX',.TRUE.,0,'matrix.txt','(E13.2)')
!    
!    ! Get the global solution/rhs/temp vector.
!    CALL sptivec_convertSupervecToVector (rx, rxGlobal)
!    CALL sptivec_convertSupervecToVector (rd, rbGlobal)
!    CALL lsysbl_createVecBlockIndirect (rbGlobal,rdGlobal,.FALSE.)
!    
!    ! Initialise the UMFPACK solver.
!    CALL linsol_initUMFPACK4 (rsolverNode)
!    
!    ! Add the matrices
!    Rmatrices(1) = rglobalA
!    CALL linsol_setMatrices (rsolverNode,Rmatrices)
!    
!    ! Init the solver
!    CALL linsol_initStructure (rsolverNode,ierror)
!    CALL linsol_initData (rsolverNode,ierror)
!    
!    ! Reshape the x,b and d-vector to fit to our matrix.
!    CALL lsysbl_deriveSubvector(rxGlobal,rxGlobalSolve,bshare=.TRUE.)
!    CALL lsysbl_deriveSubvector(rbGlobal,rbGlobalSolve,bshare=.TRUE.)
!    CALL lsysbl_deriveSubvector(rdGlobal,rdGlobalSolve,bshare=.TRUE.)
!    ALLOCATE(Isize(rxGlobal%nblocks*6))
!    Isize2= (/rtempVectorB%RvectorBlock(1)%NEQ,&
!              rtempVectorB%RvectorBlock(2)%NEQ,&
!              rtempVectorB%RvectorBlock(3)%NEQ,&
!              rtempVectorB%RvectorBlock(4)%NEQ,&
!              rtempVectorB%RvectorBlock(5)%NEQ,&
!              rtempVectorB%RvectorBlock(6)%NEQ &
!             /)
!    DO i=0,rxGlobal%nblocks-1
!      Isize(i*6+1:i*6+6) = Isize2(1:6)
!    END DO
!    CALL lsysbl_enforceStructureDirect (Isize,rxGlobalSolve)
!    CALL lsysbl_enforceStructureDirect (Isize,rbGlobalSolve)
!    CALL lsysbl_enforceStructureDirect (Isize,rdGlobalSolve)
!    
!    ! Solve
!    CALL lsysbl_getbase_double (rxGlobalSolve,p_Dx)
!    CALL lsysbl_getbase_double (rbGlobalSolve,p_Db)
!    CALL linsol_solveAdaptively (rsolverNode,rxGlobalSolve,rbGlobalSolve,rdGlobalSolve)
!    
!    ! Release
!    CALL lsysbl_releaseVector (rxGlobalSolve)
!    CALL lsysbl_releaseVector (rbGlobalSolve)
!    CALL lsysbl_releaseVector (rdGlobalSolve)
!
!    CALL linsol_releaseSolver (rsolverNode)
!    CALL lsysbl_releaseVector (rdGlobal)
!    CALL lsysbl_releaseVector (rbGlobal)
!    
!    ! Remember the solution
!    CALL sptivec_convertVectorToSupervec (rxGlobal, rx) 
!    
!    ! Release the global matrix
!    CALL lsysbl_releaseMatrix (rglobalA)
!
!    CALL lsysbl_releaseVector (rxGlobal)
!    
!  CONTAINS
!  
!    SUBROUTINE insertMatrix (rsource,rdest,ileft,itop)
!    
!    ! Includes rsource into rdest at position ileft,itop
!    TYPE(t_matrixBlock), INTENT(IN) :: rsource
!    TYPE(t_matrixBlock), INTENT(INOUT) :: rdest
!    INTEGER, INTENT(IN) :: ileft
!    INTEGER, INTENT(IN) :: itop
!    
!    INTEGER :: i,j
!    
!    DO j=1,rsource%ndiagBlocks
!      DO i=1,rsource%ndiagBlocks
!        IF (lsysbl_isSubmatrixPresent (rsource,i,j)) THEN
!          IF (lsysbl_isSubmatrixPresent (rdest,i+itop-1,j+ileft-1)) THEN
!            CALL lsyssc_releaseMatrix (rdest%RmatrixBlock(j+ileft-1,i+itop-1))
!          END IF
!          CALL lsyssc_duplicateMatrix (rsource%RmatrixBlock(i,j),&
!              rdest%RmatrixBlock(i+itop-1,j+ileft-1),&
!              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        END IF
!      END DO
!    END DO
!        
!    END SUBROUTINE
!    
!  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
      isubstep,irelpos,rnonlinearIteration)

!<description>
  ! This routine sets up the matrix weights in a nonlinear iteration structure
  ! according to the position of the corresponding submatrix.
  ! isubstep defines the timestep which is to be tackled -- which corresponds
  ! to the row in the supermatrix.
  ! irelpos specifies the column in the supermatrix relative to the diagonal --
  ! thus =0 means the diagonal, =1 the submatrix above the diagonal and =-1
  ! the submatrix below the diagonal.
  !
  ! The matrix weights in rnonlinearIteration are initialised based on this
  ! information such that the 'assembleLinearisedMatrix' routine of the core
  ! equation will build the system matrix at that position in the supermatrix.
!</description>

!<input>
  ! Problem structure
  TYPE(t_problem), INTENT(IN) :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr

  ! Theta scheme identifier.
  ! = 0.5: Crank-Nicolson.
  ! = 1.0: Explicit Euler
  REAL(DP), INTENT(IN) :: dtheta
  
  ! Substep in the time-dependent simulation = row in the supermatrix.
  ! Range 0..nsubsteps
  INTEGER, INTENT(IN) :: isubstep
  
  ! Specifies the column in the supermatrix relative to the diagonal.
  ! =0: set matrix weights for the diagonal.
  INTEGER, INTENT(IN) :: irelpos
!</input>

!<inputoutput>
  ! Nonlinear-iteration-structure for the core-equation module. 
  ! Must have been initialised with 'initNonlinearLoop'. The weights
  ! that specify the submatrices of a small 6x6 block matrix system are initialised
  ! depending on the position specified by isubstep and nsubsteps.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! If the following constant is set from 1.0 to 0.0, the primal system is
    ! decoupled from the dual system!
    REAL(DP), PARAMETER :: dprimalDualCoupling = 1.0_DP
    
    ! If the following constant is set from 1.0 to 0.0, the dual system is
    ! decoupled from the primal system!
    REAL(DP), PARAMETER :: ddualPrimalCoupling = 1.0_DP
    
    ! If the following parameter is set from 1.0 to 0.0, the terminal
    ! condition between the primal and dual equation is decoupled, i.e.
    ! the dual equation gets independent from the primal one.
    REAL(DP), PARAMETER :: dterminalCondDecoupled = 1.0_DP

    ! The first and last substep is a little bit special concerning
    ! the matrix!
    IF (isubstep .EQ. 0) THEN
      
      ! We are in the first substep
    
      IF (irelpos .EQ. 0) THEN
      
        ! The diagonal matrix.
        rnonlinearIteration%diota1 = 1.0_DP
        rnonlinearIteration%diota2 = 0.0_DP

        rnonlinearIteration%dkappa1 = 1.0_DP
        rnonlinearIteration%dkappa2 = 0.0_DP
        
        rnonlinearIteration%dalpha1 = 0.0_DP
        rnonlinearIteration%dalpha2 = 1.0_DP
        
        rnonlinearIteration%dtheta1 = 0.0_DP
        rnonlinearIteration%dtheta2 = dtheta * rspaceTimeDiscr%dtstep
        
        rnonlinearIteration%dgamma1 = 0.0_DP
        rnonlinearIteration%dgamma2 = &
            dtheta * rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
        
        rnonlinearIteration%deta1 = 0.0_DP
        rnonlinearIteration%deta2 = rspaceTimeDiscr%dtstep
        
        rnonlinearIteration%dtau1 = 0.0_DP
        rnonlinearIteration%dtau2 = 1.0_DP
        
        rnonlinearIteration%dmu1 = 0.0_DP
        rnonlinearIteration%dmu2 = ddualPrimalCoupling * &
            (-rspaceTimeDiscr%dtstep * dtheta)
        
      ELSE IF (irelpos .EQ. 1) THEN
      
        ! Offdiagonal matrix on the right of the diagonal.

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rnonlinearIteration%diota1 = 0.0_DP
        rnonlinearIteration%diota2 = 0.0_DP

        rnonlinearIteration%dkappa1 = 0.0_DP
        rnonlinearIteration%dkappa2 = 0.0_DP
        
        rnonlinearIteration%dalpha1 = 0.0_DP
        rnonlinearIteration%dalpha2 = -1.0_DP
        
        rnonlinearIteration%dtheta1 = 0.0_DP
        rnonlinearIteration%dtheta2 = rspaceTimeDiscr%dtstep * (1.0_DP-dtheta)
        
        rnonlinearIteration%dgamma1 = 0.0_DP
        rnonlinearIteration%dgamma2 = &
            rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP) * (1.0_DP-dtheta)
        
        rnonlinearIteration%deta1 = 0.0_DP
        rnonlinearIteration%deta2 = 0.0_DP
        
        rnonlinearIteration%dtau1 = 0.0_DP
        rnonlinearIteration%dtau2 = 0.0_DP
        
        rnonlinearIteration%dmu1 = 0.0_DP
        rnonlinearIteration%dmu2 = ddualPrimalCoupling * &
            (-rspaceTimeDiscr%dtstep) * (1.0_DP-dtheta)
            
      END IF
    
    ELSE IF (isubstep .LT. rspaceTimeDiscr%niterations) THEN
      
      ! We are sonewhere in the middle of the matrix. There is a substep
      ! isubstep+1 and a substep isubstep-1!
      
      ! -----
      
      IF (irelpos .EQ. -1) THEN
      
        ! Matrix on the left of the diagonal.
        !
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rnonlinearIteration%diota1 = 0.0_DP
        rnonlinearIteration%diota2 = 0.0_DP

        rnonlinearIteration%dkappa1 = 0.0_DP
        rnonlinearIteration%dkappa2 = 0.0_DP
        
        rnonlinearIteration%dalpha1 = -1.0_DP
        rnonlinearIteration%dalpha2 = 0.0_DP
        
        rnonlinearIteration%dtheta1 = (1.0_DP-dtheta) * rspaceTimeDiscr%dtstep
        rnonlinearIteration%dtheta2 = 0.0_DP
        
        rnonlinearIteration%dgamma1 = &
            (1.0_DP-dtheta) * rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
        rnonlinearIteration%dgamma2 = 0.0_DP
        
        rnonlinearIteration%deta1 = 0.0_DP
        rnonlinearIteration%deta2 = 0.0_DP
        
        rnonlinearIteration%dtau1 = 0.0_DP
        rnonlinearIteration%dtau2 = 0.0_DP
        
        rnonlinearIteration%dmu1 = dprimalDualCoupling * &
            rspaceTimeDiscr%dtstep * (1.0_DP-dtheta) / rspaceTimeDiscr%dalphaC
        rnonlinearIteration%dmu2 = 0.0_DP

      ELSE IF (irelpos .EQ. 0) THEN    

        ! The diagonal matrix.

        rnonlinearIteration%diota1 = 0.0_DP
        rnonlinearIteration%diota2 = 0.0_DP

        rnonlinearIteration%dkappa1 = 0.0_DP
        rnonlinearIteration%dkappa2 = 0.0_DP
        
        rnonlinearIteration%dalpha1 = 1.0_DP
        rnonlinearIteration%dalpha2 = 1.0_DP
        
        rnonlinearIteration%dtheta1 = dtheta * rspaceTimeDiscr%dtstep
        rnonlinearIteration%dtheta2 = dtheta * rspaceTimeDiscr%dtstep
        
        rnonlinearIteration%dgamma1 = &
            dtheta * rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
        rnonlinearIteration%dgamma2 = &
            dtheta * rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
        
        rnonlinearIteration%deta1 = rspaceTimeDiscr%dtstep
        rnonlinearIteration%deta2 = rspaceTimeDiscr%dtstep
        
        rnonlinearIteration%dtau1 = 1.0_DP
        rnonlinearIteration%dtau2 = 1.0_DP
        
        rnonlinearIteration%dmu1 = dprimalDualCoupling * &
            dtheta * rspaceTimeDiscr%dtstep / rspaceTimeDiscr%dalphaC
        rnonlinearIteration%dmu2 = ddualPrimalCoupling * &
            dtheta * (-rspaceTimeDiscr%dtstep)

      ELSE IF (irelpos .EQ. 1) THEN
            
        ! Matrix on the right of the diagonal.
        !
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        rnonlinearIteration%diota1 = 0.0_DP
        rnonlinearIteration%diota2 = 0.0_DP

        rnonlinearIteration%dkappa1 = 0.0_DP
        rnonlinearIteration%dkappa2 = 0.0_DP
        
        rnonlinearIteration%dalpha1 = 0.0_DP
        rnonlinearIteration%dalpha2 = -1.0_DP
        
        rnonlinearIteration%dtheta1 = 0.0_DP
        rnonlinearIteration%dtheta2 = (1.0_DP-dtheta) * rspaceTimeDiscr%dtstep
        
        rnonlinearIteration%dgamma1 = 0.0_DP
        rnonlinearIteration%dgamma2 = &
            (1.0_DP-dtheta) * rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
        
        rnonlinearIteration%deta1 = 0.0_DP
        rnonlinearIteration%deta2 = 0.0_DP
        
        rnonlinearIteration%dtau1 = 0.0_DP
        rnonlinearIteration%dtau2 = 0.0_DP
        
        rnonlinearIteration%dmu1 = 0.0_DP
        rnonlinearIteration%dmu2 = ddualPrimalCoupling * &
            (1.0_DP-dtheta) * (-rspaceTimeDiscr%dtstep)
            
      END IF
    
    ELSE
    
      ! We are in the last substep
      
      IF (irelpos .EQ. -1) THEN
      
        ! Matrix on the left of the diagonal
        !
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        
        rnonlinearIteration%diota1 = 0.0_DP
        rnonlinearIteration%diota2 = 0.0_DP

        rnonlinearIteration%dkappa1 = 0.0_DP
        rnonlinearIteration%dkappa2 = 0.0_DP
        
        rnonlinearIteration%dalpha1 = -1.0_DP
        rnonlinearIteration%dalpha2 = 0.0_DP
        
        rnonlinearIteration%dtheta1 = (1.0_DP-dtheta) * rspaceTimeDiscr%dtstep
        rnonlinearIteration%dtheta2 = 0.0_DP
        
        rnonlinearIteration%dgamma1 = &
            (1.0_DP-dtheta) * rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
        rnonlinearIteration%dgamma2 = 0.0_DP
        
        rnonlinearIteration%deta1 = 0.0_DP
        rnonlinearIteration%deta2 = 0.0_DP
        
        rnonlinearIteration%dtau1 = 0.0_DP
        rnonlinearIteration%dtau2 = 0.0_DP
        
        rnonlinearIteration%dmu1 = dprimalDualCoupling * &
            rspaceTimeDiscr%dtstep * (1.0_DP-dtheta) / rspaceTimeDiscr%dalphaC
        rnonlinearIteration%dmu2 = 0.0_DP
        
      ELSE IF (irelpos .EQ. 0) THEN
      
        ! The diagonal matrix.
        
        rnonlinearIteration%diota1 = 0.0_DP
        rnonlinearIteration%diota2 = (1.0_DP-ddualPrimalCoupling)

        rnonlinearIteration%dkappa1 = 0.0_DP
        rnonlinearIteration%dkappa2 = 1.0_DP
        
        rnonlinearIteration%dalpha1 = 1.0_DP
        rnonlinearIteration%dalpha2 = 1.0_DP
        
        rnonlinearIteration%dtheta1 = dtheta * rspaceTimeDiscr%dtstep
        rnonlinearIteration%dtheta2 = 0.0_DP
        
        rnonlinearIteration%dgamma1 = &
            dtheta * rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
        rnonlinearIteration%dgamma2 = 0.0_DP
        
        rnonlinearIteration%deta1 = rspaceTimeDiscr%dtstep
        rnonlinearIteration%deta2 = 0.0_DP
        
        rnonlinearIteration%dtau1 = 1.0_DP
        rnonlinearIteration%dtau2 = 0.0_DP
        
        rnonlinearIteration%dmu1 = dprimalDualCoupling * &
            dtheta * rspaceTimeDiscr%dtstep / rspaceTimeDiscr%dalphaC
        rnonlinearIteration%dmu2 = ddualPrimalCoupling * &
            (-rspaceTimeDiscr%dgammaC)
        
      END IF

    END IF

  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_solveSupersysDirectCN (rproblem, rspaceTimeDiscr, rx, rd, &
      rtempvectorX, rtempvectorB, rtempvectorD)

!<description>
  ! This routine assembles and solves the time-space coupled supersystem:
  ! $Ax=b$. The RHS vector is generated on-the-fly.
  ! The routine generates the full matrix in memory and solves with UMFPACK, 
  ! so it should only be used for debugging!
  ! The routine assembles the global matrix with the time step scheme
  ! specified in dtimeStepTheta in the problem structure.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector defining the current solution.
  ! Is replaced by the new solution
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX

  ! A second temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorB

  ! A third temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD

  ! A space-time vector that receives the defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubstep,ilevel,ierror,i
    TYPE(t_matrixBlock) :: rblockTemp
    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    TYPE(t_vectorBlock) :: rxGlobal, rbGlobal, rdGlobal,rtempVectorRHS
    TYPE(t_vectorBlock) :: rxGlobalSolve, rbGlobalSolve, rdGlobalSolve
    TYPE(t_matrixBlock) :: rglobalA
    TYPE(t_linsolNode), POINTER :: rsolverNode,p_rpreconditioner
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: Isize
    INTEGER(PREC_VECIDX), DIMENSION(6) :: Isize2
    REAL(DP) :: dtheta
    
    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd

    ! If the following constant is set from 1.0 to 0.0, the primal system is
    ! decoupled from the dual system!
    REAL(DP), PARAMETER :: dprimalDualCoupling = 1.0_DP
    
    ! If the following constant is set from 1.0 to 0.0, the dual system is
    ! decoupled from the primal system!
    REAL(DP), PARAMETER :: ddualPrimalCoupling = 1.0_DP
    
    ! If the following parameter is set from 1.0 to 0.0, the terminal
    ! condition between the primal and dual equation is decoupled, i.e.
    ! the dual equation gets independent from the primal one.
    REAL(DP), PARAMETER :: dterminalCondDecoupled = 1.0_DP

    ilevel = rspaceTimeDiscr%NLMAX
    
    ! Theta-scheme identifier.
    ! =1: impliciz Euler.
    ! =0.5: Crank Nicolson
    dtheta = rproblem%rtimedependence%dtimeStepTheta
    
    ! Assemble the space-time RHS into rd.
    CALL c2d2_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rd, &
      rtempvectorX, rtempvectorB, rtempvectorD, .TRUE.)
      
    ! Implement the initial condition into the RHS.
    CALL c2d2_implementInitCondRHS (rx, rd, rtempvectorX, rtempvectorD)    

    ! ----------------------------------------------------------------------
    ! 2.) Generate the matrix A
    !
    ! Create a global matrix:
    CALL lsysbl_createEmptyMatrix (rglobalA,6*(rspaceTimeDiscr%niterations+1))
    
    ! Loop through the substeps
    
    DO isubstep = 0,rspaceTimeDiscr%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%dtstep

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! The first and last substep is a little bit special concerning
      ! the matrix!
      IF (isubstep .EQ. 0) THEN
        
        ! We are in the first substep
      
        ! -----
        
        ! The diagonal matrix.
        !
        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,0,rnonlinearIterationTmp)
          
        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (&
            rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
            
        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rglobalA,1,1)
            
        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

        ! -----
      
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the primal velocity.

        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,1,rnonlinearIterationTmp)
      
        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        
        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
        ! Include the boundary conditions into that matrix.
        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
        ! main diagonal of the supermatrix.
        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)

        ! We don't need submatrix (1,1) to (2,2).
        rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP

        ! Include that in the global matrix below the diagonal
        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1+6,isubstep*6+1)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)
        
      ELSE IF (isubstep .LT. rspaceTimeDiscr%niterations) THEN
        
        ! We are sonewhere in the middle of the matrix. There is a substep
        ! isubstep+1 and a substep isubstep-1!
        
        ! -----
        
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the primal velocity.

        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,-1,rnonlinearIterationTmp)

        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (&
            rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        
        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
        ! We don't need submatrix (4,4) to (5,5).
        rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP

        ! Include the boundary conditions into that matrix.
        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
        ! main diagonal of the supermatrix.
        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)

        ! Include that in the global matrix below the diagonal
        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)
        
        ! -----      

        ! Now the diagonal matrix.
        !
        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
      
        ! Assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,0,rnonlinearIterationTmp)
            
        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        
        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rglobalA,isubstep*6+1,isubstep*6+1)
            
        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

        ! -----
        
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,1,rnonlinearIterationTmp)
      
        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        
        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
        ! We don't need submatrix (1,1) to (2,2).
        rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP

        ! Include the boundary conditions into that matrix.
        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
        ! main diagonal of the supermatrix.
        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)

        ! Include that in the global matrix above the diagonal
        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1+6,isubstep*6+1)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

      ELSE
      
        ! We are in the last substep
        
        ! -----
        
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,-1,rnonlinearIterationTmp)
      
        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (&
            rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        
        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
        ! We don't need submatrix (4,4) to (5,5).
        rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP

        ! Include the boundary conditions into that matrix.
        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
        ! main diagonal of the supermatrix.
        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)

        ! Include that in the global matrix above the diagonal
        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)
     
        ! -----
        
        ! The diagonal matrix.
        !
        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
      
        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,0,rnonlinearIterationTmp)
        
        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        
        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rglobalA,isubstep*6+1,isubstep*6+1)

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)
      
        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

      END IF
    
    END DO
    
    ! Update structural information of the global matrix.
    CALL lsysbl_updateMatStrucInfo (rglobalA)

    ! Write the global matrix to a file.
    !CALL matio_writeBlockMatrixHR(rglobalA,'MATRIX',.TRUE.,0,'matrixcn.txt','(E13.2)')
    
    ! Get the global solution/rhs/temp vector.
    CALL sptivec_convertSupervecToVector (rx, rxGlobal)
    CALL sptivec_convertSupervecToVector (rd, rbGlobal)
    CALL lsysbl_createVecBlockIndirect (rbGlobal,rdGlobal,.FALSE.)
    
    ! Initialise the UMFPACK solver.
    !CALL linsol_initUMFPACK4 (rsolverNode)
    CALL linsol_initUMFPACK4 (p_rpreconditioner)
    CALL linsol_initDefCorr (rsolverNode,p_rpreconditioner)
    rsolverNode%ioutputLevel = 2
    rsolverNode%depsRel = 1E-10
    rsolverNode%nmaxIterations = 10
    
    ! Add the matrices
    Rmatrices(1) = rglobalA
    CALL linsol_setMatrices (rsolverNode,Rmatrices)
    
    ! Init the solver
    CALL linsol_initStructure (rsolverNode,ierror)
    CALL linsol_initData (rsolverNode,ierror)
    
    ! Reshape the x,b and d-vector to fit to our matrix.
    CALL lsysbl_deriveSubvector(rxGlobal,rxGlobalSolve,bshare=.TRUE.)
    CALL lsysbl_deriveSubvector(rbGlobal,rbGlobalSolve,bshare=.TRUE.)
    CALL lsysbl_deriveSubvector(rdGlobal,rdGlobalSolve,bshare=.TRUE.)
    ALLOCATE(Isize(rxGlobal%nblocks*6))
    Isize2= (/rtempVectorB%RvectorBlock(1)%NEQ,&
              rtempVectorB%RvectorBlock(2)%NEQ,&
              rtempVectorB%RvectorBlock(3)%NEQ,&
              rtempVectorB%RvectorBlock(4)%NEQ,&
              rtempVectorB%RvectorBlock(5)%NEQ,&
              rtempVectorB%RvectorBlock(6)%NEQ &
             /)
    DO i=0,rxGlobal%nblocks-1
      Isize(i*6+1:i*6+6) = Isize2(1:6)
    END DO
    CALL lsysbl_enforceStructureDirect (Isize,rxGlobalSolve)
    CALL lsysbl_enforceStructureDirect (Isize,rbGlobalSolve)
    CALL lsysbl_enforceStructureDirect (Isize,rdGlobalSolve)
    
    ! Solve
    CALL lsysbl_getbase_double (rxGlobalSolve,p_Dx)
    CALL lsysbl_getbase_double (rbGlobalSolve,p_Db)
    CALL linsol_solveAdaptively (rsolverNode,rxGlobalSolve,rbGlobalSolve,rdGlobalSolve)
    
    ! DEBUG!!!
    CALL lsysbl_blockMatVec (rglobalA,rxGlobalSolve,rbGlobalSolve,-1.0_DP,1.0_DP)
    !CALL lsyssc_scalarMatVec (rglobalA%RmatrixBlock(16,13),&
    !                          rxGlobalSolve%RvectorBlock(13),&
    !                          rbGlobalSolve%RvectorBlock(16),1.0_DP,-1.0_DP)
    
    ! Release
    CALL lsysbl_releaseVector (rxGlobalSolve)
    CALL lsysbl_releaseVector (rbGlobalSolve)
    CALL lsysbl_releaseVector (rdGlobalSolve)

    CALL linsol_releaseSolver (rsolverNode)
    CALL lsysbl_releaseVector (rdGlobal)
    CALL lsysbl_releaseVector (rbGlobal)
    
    ! Remember the solution
    CALL sptivec_convertVectorToSupervec (rxGlobal, rx) 
    
    ! Release the global matrix
    CALL lsysbl_releaseMatrix (rglobalA)

    CALL lsysbl_releaseVector (rxGlobal)
    
  CONTAINS
  
    SUBROUTINE insertMatrix (rsource,rdest,ileft,itop)
    
    ! Includes rsource into rdest at position ileft,itop
    TYPE(t_matrixBlock), INTENT(IN) :: rsource
    TYPE(t_matrixBlock), INTENT(INOUT) :: rdest
    INTEGER, INTENT(IN) :: ileft
    INTEGER, INTENT(IN) :: itop
    
    INTEGER :: i,j
    
    DO j=1,rsource%ndiagBlocks
      DO i=1,rsource%ndiagBlocks
        IF (lsysbl_isSubmatrixPresent (rsource,i,j)) THEN
          IF (lsysbl_isSubmatrixPresent (rdest,i+itop-1,j+ileft-1)) THEN
            CALL lsyssc_releaseMatrix (rdest%RmatrixBlock(j+ileft-1,i+itop-1))
          END IF
          CALL lsyssc_duplicateMatrix (rsource%RmatrixBlock(i,j),&
              rdest%RmatrixBlock(i+itop-1,j+ileft-1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        END IF
      END DO
    END DO
        
    END SUBROUTINE
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
      rtempvector1, rtempvector2, rtempvector3, bimplementBC)

!<description>
  ! Assembles the space-time RHS vector rb. Bondary conditions are NOT
  ! implemented!
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
!</input>

!<inputoutput>
  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector1

  ! A second temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector2

  ! A third temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector3

  ! A space-time vector that receives the RHS.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
  
  ! Whether to implement boundary conditions into the RHS or not.
  LOGICAL, INTENT(IN) :: bimplementBC
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubstep
    REAL(DP) :: dtheta,dtstep
    
    ! A temporary vector for the creation of the RHS.
    TYPE(t_vectorBlock) :: rtempVectorRHS
    
    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd, p_Drhs

    ! Theta-scheme identifier.
    ! =1: impliciz Euler.
    ! =0.5: Crank Nicolson
    dtheta = rproblem%rtimedependence%dtimeStepTheta
    dtstep = rspaceTimeDiscr%dtstep
    
    ! ----------------------------------------------------------------------
    ! Generate the global RHS vector
    
    CALL lsysbl_getbase_double (rtempVector1,p_Dx)
    CALL lsysbl_getbase_double (rtempVector2,p_Db)
    CALL lsysbl_getbase_double (rtempVector3,p_Dd)

    ! Assemble 1st RHS vector in X temp vector.
    CALL generateRHS (rproblem,0,rb%ntimesteps,&
        rtempVector1, .TRUE., .FALSE.)
    
    ! Assemble the 2nd RHS vector in the RHS temp vector
    CALL generateRHS (rproblem,1,rb%ntimesteps,&
        rtempVector2, .TRUE., .FALSE.)

    ! Assemble the 3rd RHS vector in the defect temp vector
    CALL generateRHS (rproblem,2,rb%ntimesteps,&
        rtempVector3, .TRUE., .FALSE.)
        
    ! Create a copy of the X temp vector (RHS0). That vector will be
    ! our destination vector for assembling the RHS in all timesteps.
    CALL lsysbl_copyVector (rtempVector1,rtempVectorRHS)
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
    
    ! RHS 0,1,2 -> 1-2-3
    
    DO isubstep = 0,rspaceTimeDiscr%niterations
    
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
            dtstep*dtheta,dtstep*(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(4))
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(5),rtempVector2%RvectorBlock(5),&
            dtstep*dtheta,dtstep*(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(5))
        ! Pressure is fully implicit; no weighting by dtstep!
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(6),rtempVector2%RvectorBlock(6),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(6))
            
      ELSE IF (isubstep .LT. rspaceTimeDiscr%niterations) THEN
      
        ! We are somewhere 'in the middle'.
        !
        ! Dual RHS comes from rtempVector3. The primal from the
        ! isubstep-1'th RHS.
        !
        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
        
        CALL lsyssc_vectorLinearComb (&
            rtempVector1%RvectorBlock(1),rtempVector2%RvectorBlock(1),&
            dtstep*(1.0_DP-dtheta),dtstep*dtheta,&
            rtempVectorRHS%RvectorBlock(1))                                        
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(2),rtempVector2%RvectorBlock(2),&
            dtstep*(1.0_DP-dtheta),dtstep*dtheta,&
            rtempVectorRHS%RvectorBlock(2))                                        
        ! Pressure is fully implicit; no weighting by dtstep!
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector1%RvectorBlock(3),rtempVector2%RvectorBlock(3),&
            (1.0_DP-dtheta),dtheta,&
            rtempVectorRHS%RvectorBlock(3))

        CALL lsyssc_vectorLinearComb (&
            rtempVector2%RvectorBlock(4),rtempVector3%RvectorBlock(4),&
            dtstep*dtheta,dtstep*(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(4))
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(5),rtempVector3%RvectorBlock(5),&
            dtstep*dtheta,dtstep*(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(5))
        ! Pressure is fully implicit; no weighting by dtstep!
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(6),rtempVector3%RvectorBlock(6),&
            dtheta,(1.0_DP-dtheta),&
            rtempVectorRHS%RvectorBlock(6))
        
        IF (isubstep .LT. rspaceTimeDiscr%niterations-1) THEN
          ! Shift the RHS vectors and generate the RHS for the next time step.
          ! (Yes, I know, this could probably be solved more elegant without copying anything
          ! using a ring buffer ^^)
          CALL lsysbl_copyVector(rtempVector2,rtempVector1)
          CALL lsysbl_copyVector(rtempVector3,rtempVector2)
          CALL generateRHS (rproblem,isubstep+2,rspaceTimeDiscr%niterations,&
              rtempVector3, .TRUE., .FALSE.)
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
            dtstep*(1.0_DP-dtheta),dtstep*dtheta,&
            rtempVectorRHS%RvectorBlock(1))                                        
        CALL lsyssc_vectorLinearComb (&                                                   
            rtempVector2%RvectorBlock(2),rtempVector3%RvectorBlock(2),&
            dtstep*(1.0_DP-dtheta),dtstep*dtheta,&
            rtempVectorRHS%RvectorBlock(2))                                        
        ! Pressure is fully implicit; no weighting by dtstep!
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
        CALL generateRHS (rproblem,isubstep,rspaceTimeDiscr%niterations,&
            rtempVectorRHS, .FALSE., .TRUE.)
      END IF
      
      ! Save the RHS.
      CALL sptivec_setTimestepData(rb, isubstep, rtempVectorRHS)
      
    END DO
    
    ! Release the temp vector for generating the RHS.
    CALL lsysbl_releaseVector (rtempVectorRHS)

  CONTAINS
  
    SUBROUTINE generateRHS (rproblem,isubstep,nsubsteps,&
        rvector, bgenerate, bincludeBC)
    
    ! Generate the RHS vector of timestep isubstep and/or include boundary
    ! conditions.
    
    ! Problem structure.
    TYPE(t_problem), INTENT(INOUT) :: rproblem
    
    ! Number of the substep where to generate the RHS vector
    INTEGER, INTENT(IN) :: isubstep
    
    ! Total number of substeps
    INTEGER, INTENT(IN) :: nsubsteps
    
    ! Destination vector
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
    
    ! Whether to generate the RHS vector or not
    LOGICAL, INTENT(IN) :: bgenerate
    
    ! Whether to include boundary conditions
    LOGICAL, INTENT(IN) :: bincludeBC
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
      ! DEBUG!!!
      CALL lsysbl_getbase_double (rvector,p_Ddata)

      ! Set the time where we are at the moment
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*dtstep
          
      ! Assemble the RHS?
      IF (bgenerate) THEN
      
        ! Generate the RHS of that point in time into the vector.
        CALL c2d2_generateBasicRHS (rproblem,rvector)

      END IF
      
      ! Include BC's?
      IF (bincludeBC) THEN
        
        ! Initialise the collection for the assembly process with callback routines.
        ! Basically, this stores the simulation time in the collection if the
        ! simulation is nonstationary.
        CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
        END IF

        ! Implement the boundary conditions into the RHS.
        ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
        ! as Dirichlet values mustn't be multiplied with GAMMA!
        CALL vecfil_discreteBCsol (rvector)
        CALL vecfil_discreteFBCsol (rvector)      
      
        ! Clean up the collection (as we are done with the assembly, that's it.
        CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)
        
      END IF
    
    END SUBROUTINE
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_implementInitCondRHS (rx, rb, rtempvectorX, rtempvectorD)

!<description>
  ! Implements the initial condition into the RHS vector rb.
  ! Overwrites the rb of the first time step.
  !
  ! Does not implement boundary conditions!
!</description>

!<input>
  ! A space-time vector containing the initial condition in the first subvector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
!</input>

!<inputoutput>
  ! A space-time vector with the RHS. The initial condition is implemented into
  ! this vector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX

  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD
!</inputoutput>

!</subroutine>

    ! Overwrite the primal RHS with the initial primal solution vector.
    ! This realises the inital condition.
    CALL sptivec_getTimestepData(rx, 0, rtempVectorX)
    CALL sptivec_getTimestepData(rb, 0, rtempVectorD)
    CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(1),rtempVectorD%RvectorBlock(1))
    CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(2),rtempVectorD%RvectorBlock(2))
    CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(3),rtempVectorD%RvectorBlock(3))
    CALL sptivec_setTimestepData(rb, 0, rtempVectorD)

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

  SUBROUTINE c2d2_assembleSpaceTimeDefect (rproblem, rspaceTimeDiscr, rx, rd, dnorm)

!<description>
  ! This routine assembles the space-time defect d=b-Ax. rd must have been
  ! initialised by the space-time RHS vector b. The routine will then
  ! calculate rd = rd - A rx to get the space-time defect.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector defining the current solution.
  ! Is replaced by the new solution
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A space-time vector with the space-time RHS b. Is overwritten by
  ! d=b-Ax.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!<output>
  ! OPTIONAL: If specified, returns the $l_2$-norm if the defect.
  REAL(DP), INTENT(OUT), OPTIONAL :: dnorm
!<output>

!</subroutine>

    ! local variables
    INTEGER :: isubstep,ilevel
    TYPE(t_vectorBlock) :: rtempVectorD, rtempVector1, rtempVector2, rtempVector3
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscr
    REAL(DP) :: dtemp,dtheta
    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    TYPE(t_matrixBlock) :: rblockTemp
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx1,p_Dx2,p_Dx3,p_Db
    
    ! If the following constant is set from 1.0 to 0.0, the primal system is
    ! decoupled from the dual system!
    REAL(DP), PARAMETER :: dprimalDualCoupling = 1.0_DP
    
    ! If the following constant is set from 1.0 to 0.0, the dual system is
    ! decoupled from the primal system!
    REAL(DP), PARAMETER :: ddualPrimalCoupling = 1.0_DP
    
    ! If the following parameter is set from 1.0 to 0.0, the terminal
    ! condition between the primal and dual equation is decoupled, i.e.
    ! the dual equation gets independent from the primal one.
    REAL(DP), PARAMETER :: dterminalCondDecoupled = 1.0_DP

    ! Level of the discretisation
    ilevel = rspaceTimeDiscr%NLMAX

    ! Theta-scheme identifier.
    ! =1: implicit Euler.
    ! =0.5: Crank Nicolson
    dtheta = rproblem%rtimedependence%dtimeStepTheta
    
    ! Create a temp vector that contains the part of rd which is to be modified.
    p_rdiscr => rspaceTimeDiscr%p_RlevelInfo%p_rdiscretisation
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorD,.FALSE.)
    
    ! The vector will be a defect vector. Assign the boundary conditions so
    ! that we can implement them.
    rtempVectorD%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorD%p_rdiscreteBCfict => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC
    
    ! Create a temp vector for the X-vectors at timestep i-1, i and i+1.
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector1,.FALSE.)
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector2,.FALSE.)
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector3,.FALSE.)
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rtempVector1,p_Dx1)
    CALL lsysbl_getbase_double (rtempVector2,p_Dx2)
    CALL lsysbl_getbase_double (rtempVector3,p_Dx3)
    CALL lsysbl_getbase_double (rtempVectorD,p_Db)
    
    ! Get the parts of the X-vector which are to be modified at first --
    ! subvector 1, 2 and 3.
    CALL sptivec_getTimestepData(rx, 0, rtempVector1)
    IF (rspaceTimeDiscr%niterations .GT. 0) &
      CALL sptivec_getTimestepData(rx, 1, rtempVector2)
    IF (rspaceTimeDiscr%niterations .GT. 1) &
      CALL sptivec_getTimestepData(rx, 2, rtempVector3)
      
    ! If dnorm is specified, clear it.
    IF (PRESENT(dnorm)) THEN
      dnorm = 0.0_DP
    END IF

    ! Loop through the substeps
    
    DO isubstep = 0,rspaceTimeDiscr%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%dtstep

      ! Get the part of rd which is to be modified.
      CALL sptivec_getTimestepData(rd, isubstep, rtempVectorD)

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! The first and last substep is a little bit special concerning
      ! the matrix!
      IF (isubstep .EQ. 0) THEN
        
        ! We are in the first substep. Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( A11 A12   0   0 ... )  ( x1 )  =  ( f1 )
        !  ( ... ... ... ... ... )  ( .. )     ( .. )
        !
        ! So we have to compute:
        !
        !  d1  :=  f1  -  A11 x1  -  A12 x2
        !
        ! -----
        
        ! The diagonal matrix.
        !
        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVector1,rtempVectorD,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,0,rnonlinearIterationTmp)

        ! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        ! Include the boundary conditions into the matrices.
        !CALL c2d2_assembleLinearisedMatrices (&
        !    rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
            
        !! Subtract: rd = rd - A11 x1
        !CALL lsysbl_blockMatVec (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
        !    rtempVector1,rtempVectorD,-1.0_DP,1.0_DP)

        ! Subtract: rd = rd - A11 x1
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
          rtempVector1,rtempVectorD,&
          .FALSE.,.FALSE.,rproblem%rcollection,&
          rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

        ! -----
      
        ! Create the matrix
        !   A12  :=  -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and subtract A12 x2 from rd.
        !
        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVector2,rtempVectorD,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,1,rnonlinearIterationTmp)

        !! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        !! Include the boundary conditions into the matrices.
        !CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        !
        !CALL lsysbl_duplicateMatrix (&
        !    rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
        !    rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        !
        !! Include the boundary conditions into that matrix.
        !! Specify the matrix as 'off-diagonal' matrix because it's not on the
        !! main diagonal of the supermatrix.
        !rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        !CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        !CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
        !
        !! We don't need submatrix (1,1) to (2,2).
        !rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP
        !
        !! Subtract: rd = rd - A12 x2
        !CALL lsysbl_blockMatVec (rblockTemp,rtempVector2,rtempVectorD,-1.0_DP,1.0_DP)
 
        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
        ! We don't need submatrix (1,1) to (2,2).
        rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP

        ! Subtract: rd = rd - A12 x2
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
          rtempVector2,rtempVectorD,&
          .FALSE.,.FALSE.,rproblem%rcollection,rblockTemp)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

      ELSE IF (isubstep .LT. rspaceTimeDiscr%niterations) THEN

        ! We are sonewhere in the middle of the matrix. There is a substep
        ! isubstep+1 and a substep isubstep-1!  Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
        !  ( ... Aii-1 Aii Aii+1 ... )  ( xi )  =  ( fi )
        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
        !
        ! So we have to compute:
        !
        !  dn  :=  fn  -  Aii-1 xi-1  -  Aii xi  -  Aii+1 xi+1
        
        ! -----
        
        ! Create the matrix
        !   Aii-1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the primal velocity.

        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVector1,rtempVectorD,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,-1,rnonlinearIterationTmp)
            
        !! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        !! Include the boundary conditions into the matrices.
        !CALL c2d2_assembleLinearisedMatrices (&
        !    rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        !
        !CALL lsysbl_duplicateMatrix (&
        !    rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
        !    rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        !
        !! We don't need submatrix (4,4) to (5,5).
        !rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP
        !
        !! Include the boundary conditions into that matrix.
        !! Specify the matrix as 'off-diagonal' matrix because it's not on the
        !! main diagonal of the supermatrix.
        !rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        !CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        !CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
        !
        !! Subtract: rd = rd - Aii-1 xi-1
        !CALL lsysbl_blockMatVec (rblockTemp,rtempVector1,rtempVectorD,-1.0_DP,1.0_DP)
        
        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
        ! We don't need submatrix (4,4) to (5,5).
        rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP
        
        ! Subtract: rd = rd - Aii-1 xi-1
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
          rtempVector1,rtempVectorD,&
          .FALSE.,.FALSE.,rproblem%rcollection,rblockTemp)

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

        ! -----      

        ! Now the diagonal matrix.
        !
        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
      
        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVector2,rtempVectorD,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,0,rnonlinearIterationTmp)

        !! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        !! Include the boundary conditions into the matrices.
        !CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        !
        !! Subtract: rd = rd - Aii xi
        !CALL lsysbl_blockMatVec (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
        !    rtempVector2,rtempVectorD,-1.0_DP,1.0_DP)

        ! Subtract: rd = rd - Aii xi
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
          rtempVector2,rtempVectorD,&
          .FALSE.,.FALSE.,rproblem%rcollection,&
          rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)
            
        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

        ! -----
        
        ! Create the matrix
        !   Aii+1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVector3,rtempVectorD,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,1,rnonlinearIterationTmp)

        !! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        !! Include the boundary conditions into the matrices.
        !CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        !
        !CALL lsysbl_duplicateMatrix (&
        !    rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
        !    rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        !
        !! We don't need submatrix (1,1) to (2,2).
        !rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP
        !
        !! Include the boundary conditions into that matrix.
        !! Specify the matrix as 'off-diagonal' matrix because it's not on the
        !! main diagonal of the supermatrix.
        !rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        !CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        !CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
        !
        !! Subtract: rd = rd - Aii+1 xi+1
        !CALL lsysbl_blockMatVec (rblockTemp,rtempVector3,rtempVectorD,-1.0_DP,1.0_DP)

        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
        ! We don't need submatrix (1,1) to (2,2).
        rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP

        ! Subtract: rd = rd - Aii+1 xi+1
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
          rtempVector3,rtempVectorD,&
          .FALSE.,.FALSE.,rproblem%rcollection,rblockTemp)

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)

      ELSE 
        
        ! We are in the last substep. Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( ... ... ... ...   ... )  ( .. )     ( .. )
        !  ( ... ...   0 Ann-1 Ann )  ( xn )  =  ( fn )
        !
        ! So we have to compute:
        !
        !  dn  :=  fn  -  Ann-1 xn-1  -  Ann xn
        
        ! -----
        
        ! Create the matrix
        !   Ann-1 = -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVector2,rtempVectorD,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,-1,rnonlinearIterationTmp)

        !! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        !! Include the boundary conditions into the matrices.
        !CALL c2d2_assembleLinearisedMatrices (&
        !    rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
        !
        !CALL lsysbl_duplicateMatrix (&
        !    rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
        !    rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        !
        !! We don't need submatrix (4,4) to (5,5).
        !rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP
        !
        !! Include the boundary conditions into that matrix.
        !! Specify the matrix as 'off-diagonal' matrix because it's not on the
        !! main diagonal of the supermatrix.
        !rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        !CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        !CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
        !
        !! Subtract: rd = rd - Ann-1 xn-1
        !CALL lsysbl_blockMatVec (rblockTemp,rtempVector2,rtempVectorD,-1.0_DP,1.0_DP)

        CALL lsysbl_duplicateMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
        ! We don't need submatrix (4,4) to (5,5).
        rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP

        ! Subtract: rd = rd - Ann-1 xn-1
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
          rtempVector2,rtempVectorD,&
          .FALSE.,.FALSE.,rproblem%rcollection,rblockTemp)

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)
     
        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)
     
        ! -----
        
        ! Now the diagonal matrix.
        !
        ! Generate the basic system matrix level rspaceTimeDiscr%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
      
        ! Set up a core equation structure and assemble matrices of the
        ! time discretisation scheme.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,rtempVector3,rtempVectorD,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        ! Assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
          isubstep,0,rnonlinearIterationTmp)

        !! Assemble the system matrix on level rspaceTimeDiscr%NLMAX.
        !! Include the boundary conditions into the matrices.
        !CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        ! 
        !! Subtract: rd = rd - Ann xn
        !CALL lsysbl_blockMatVec (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
        !    rtempVector3,rtempVectorD,-1.0_DP,1.0_DP)
        
        ! Subtract: rd = rd - Ann xn
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
          rtempVector3,rtempVectorD,&
          .FALSE.,.FALSE.,rproblem%rcollection,&
          rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)
      
        CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)
      
      END IF
      
      ! Implement the boundary conditions into the defect.
      CALL vecfil_discreteBCdef (rtempVectorD)
      CALL vecfil_discreteFBCdef (rtempVectorD)      
      
      ! Save the defect vector back to rd.
      CALL sptivec_setTimestepData(rd, isubstep, rtempVectorD)
      
      ! If dnorm is specified, calculate the norm of the sub-defect vector and
      ! add it to dnorm.
      IF (PRESENT(dnorm)) THEN
        dnorm = dnorm + lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)**2
      END IF
      
      IF ((isubstep .GT. 0) .AND. &
          (isubstep .LT. rspaceTimeDiscr%niterations-1)) THEN
      
        ! Shift the timestep data: x_n+1 -> x_n -> x_n-1
        IF (rtempVector2%NEQ .NE. 0) &
          CALL lsysbl_copyVector (rtempVector2, rtempVector1)
        
        IF (rtempVector3%NEQ .NE. 0) &
          CALL lsysbl_copyVector (rtempVector3, rtempVector2)
        
        ! Get the new x_n+1 for the next pass through the loop.
        CALL sptivec_getTimestepData(rx, isubstep+2, rtempVector3)
        
      END IF
    
    END DO
    
    ! If dnorm is specified, normalise it.
    ! It was calculated from rspaceTimeDiscr%niterations+1 subvectors.
    IF (PRESENT(dnorm)) THEN
      dnorm = SQRT(dnorm) / REAL(rspaceTimeDiscr%niterations+1,DP)
    END IF
    
    ! Release the nonlinear iteration structure, we don't need it anymore;
    ! the defect assembly is finished.
    CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp)
    
    ! Release the temp vectors.
    CALL lsysbl_releaseVector (rtempVector3)
    CALL lsysbl_releaseVector (rtempVector2)
    CALL lsysbl_releaseVector (rtempVector1)
    CALL lsysbl_releaseVector (rtempVectorD)
    
  END SUBROUTINE 
   
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_precondDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rb, rd, &
      rtempvectorX, rtempvectorB, rtempvectorD, rnonlinearIteration, rnlSolver)

!<description>
  ! This routine performs preconditioning with the nonlinear super-defect
  ! vector rd: $d = C^{-1} d$.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
  
  ! A space-time vector defining the current solution.
  TYPE(t_spacetimeVector), INTENT(IN) :: rx

  ! A space-time vector defining the current RHS.
  TYPE(t_spacetimeVector), INTENT(IN) :: rb

!</input>

!<inputoutput>
  ! Structure for the nonlinear iteration for solving the core equation.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
  
  ! A configuration stucture for the nonlinear solver
  TYPE(t_nlsolNode), INTENT(INOUT) :: rnlSolver
  
  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX

  ! A second temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorB
  
  ! A third temporary vector for the nonlinear iteration
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD
  
  ! A space-time vector that receives the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubstep,ilevel,i
    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    REAL(DP) :: dtheta,dtstep,deta1,deta2,domega
    LOGICAL :: bsuccess
    TYPE(t_matrixBlock) :: rblockTemp
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Db,p_Dd
    
    dtheta = rproblem%rtimedependence%dtimeStepTheta
    dtstep = rspaceTimeDiscr%dtstep

    ! Level of the discretisation
    ilevel = rspaceTimeDiscr%NLMAX
    
    ! ----------------------------------------------------------------------
    ! We use a block-Jacobi scheme for preconditioning...
    !
    ! For this purpose, loop through the substeps.
    
    DO isubstep = 0,rspaceTimeDiscr%niterations
    
      ! Current time step?
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep * dtstep

      CALL output_separator (OU_SEP_MINUS)

      CALL output_line ('Block-Jacobi preconditioning of timestep: '//&
          TRIM(sys_siL(isubstep,10))//&
          ' Time: '//TRIM(sys_sdEL(rproblem%rtimedependence%dtime,10)))
    
      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF

      ! DEBUG!!!      
      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
      CALL lsysbl_getbase_double (rtempVectorB,p_Db)
      CALL lsysbl_getbase_double (rtempVectorD,p_Dd)

      ! Read in the RHS/solution/defect vector of the current timestep.
      CALL sptivec_getTimestepData (rx, isubstep, rtempVectorX)
      CALL sptivec_getTimestepData (rb, isubstep, rtempVectorB)
      CALL sptivec_getTimestepData (rd, isubstep, rtempVectorD)

      ! Initialise the basic system matrices.
      DO i=rproblem%NLMIN,ilevel
        CALL c2d2_generateStaticSystemMatrix (&
            rproblem%RlevelInfo(i),rproblem%RlevelInfo(i)%rmatrix,.FALSE.)
      END DO

      ! Setup the core equation for the nonlinear loop -- in the temporary
      ! rnonlinearIterationTmp structure as copy of rnonlinearIteration.
      rnonlinearIterationTmp = rnonlinearIteration

      ! Finally set up the matrix weights for the diagonal matrix
      CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeDiscr,dtheta,&
        isubstep,0,rnonlinearIterationTmp)
        
      !rnonlinearIterationTmp%deta1 = 1.0_DP
      !rnonlinearIterationTmp%deta2 = 1.0_DP

      ! Call the solver of the core equation to solve it using a nonlinear
      ! iteration. Overwrite rtempVectorX with temporary data from the
      ! nonlinear iteration.
      ! CALL c2d2_solveCoreEquation (rnlSolver,rnonlinearIterationTmp,&
      !     rtempVectorX,rtempVectorB,rproblem%rcollection,rtempVector)             
      
      ! Perform preconditioning of the defect with the method provided by the
      ! core equation module.
      CALL c2d2_precondDefectDirect (rnonlinearIterationTmp,&
        0,rtempVectorX,rtempVectorB,rtempVectorD,domega,bsuccess,rproblem%rcollection)      
    
      ! Save back the preconditioned defect.
      CALL sptivec_setTimestepData (rd, isubstep, rtempVectorD)
      
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_solveSupersystem (rproblem, rspaceTimeDiscr, rx, rb, rd, &
      rtempvectorX, rtempvectorB, rtempVector)
  
!<description>
  ! 
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector defining the initial solution. Is replaced by a new
  ! solution vector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX

  ! A second temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorB

  ! A third temporary vector for the nonlinear iteration
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector

  ! A temporary space-time vector that receives the RHS during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb

  ! A temporary space-time vector that receives the defect during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! The nonlinear solver configuration
    TYPE(t_nlsolNode) :: rnlSol
    TYPE(t_ccNonlinearIteration) :: rnonlinearIteration
    INTEGER :: isubstep,iglobIter
    LOGICAL :: bneumann

    REAL(DP) :: ddefNorm,dinitDefNorm
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx

    ! Some preparations for the nonlinear solver.
    !
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    CALL c2d2_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! or callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    CALL c2d2_initNonlinearLoop (rproblem,rproblem%NLMIN,rspaceTimeDiscr%NLMAX,&
        rtempvectorX, rtempvectorB,&
        rnonlinearIteration,'CC2D-NONLINEAR')

    ! Check the matrices if they are compatible to our
    ! preconditioner. If not, we later have to modify the matrices a little
    ! bit to make it compatible. 
    ! The result of this matrix analysis is saved to the rfinalAssembly structure 
    ! in rnonlinearIteration and allows us later to switch between these two
    ! matrix representations: Compatibility to the discretisation routines
    ! and compatibity to the preconditioner.
    ! The c2d2_checkAssembly routine below uses this information to perform
    ! the actual modification in the matrices.
    CALL c2d2_checkAssembly (rproblem,rnonlinearIteration,&
        rtempVectorB,rnonlinearIteration%rfinalAssembly)
    
    ! Using rfinalAssembly as computed above, make the matrices compatible 
    ! to our preconditioner if they are not.
    CALL c2d2_finaliseMatrices (rnonlinearIteration)
    
    ! Initialise the preconditioner for the nonlinear iteration
    CALL c2d2_preparePreconditioner (rproblem,&
        rnonlinearIteration,rtempvectorX, rtempvectorB)

    ! Implement the bondary conditions into all initial solution vectors
    DO isubstep = 0,rspaceTimeDiscr%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%dtstep

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! Implement the boundary conditions into the global solution vector.
      CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
      
      ! DEBUG!!!
      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
      
      CALL c2d2_implementBC (rproblem,rtempVectorX,rtempVectorB,.FALSE.,&
          .TRUE.,.FALSE.)
      
      CALL sptivec_setTimestepData(rx, isubstep, rtempVectorX)
      
    END DO
    
    ddefNorm = 1.0_DP
    
    ! ---------------------------------------------------------------
    ! Solve the global space-time coupled system.
    !
    ! Get the initial defect: d=b-Ax
    !CALL c2d2_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
    CALL c2d2_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
      rtempvectorX, rtempvectorB, rtempvector, .FALSE.)    

    ! Implement the initial condition into the RHS.
    CALL c2d2_implementInitCondRHS (rx, rb, rtempvectorX, rtempvector)    

    ! Now word with rd, our 'defect' vector
    CALL sptivec_copyVector (rb,rd)

    ! Assemble the defect.
    CALL c2d2_assembleSpaceTimeDefect (rproblem, rspaceTimeDiscr, rx, rd, ddefNorm)
        
    dinitDefNorm = ddefNorm
    
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    CALL output_separator (OU_SEP_EQUAL)        

    iglobIter = 0
    
!    !CALL c2d2_solveSupersysDirect (rproblem, rspaceTimeDiscr, rx, rd, &
!    !  rtempvectorX, rtempvectorB, rtempVector)
!    CALL c2d2_solveSupersysDirectCN (rproblem, rspaceTimeDiscr, rx, rd, &
!      rtempvectorX, rtempvectorB, rtempVector)
    
    DO WHILE ((ddefNorm .GT. 1.0E-5*dinitDefNorm) .AND. (ddefNorm .LT. 1.0E99) .AND. &
              (iglobIter .LT. 10))
    
      iglobIter = iglobIter+1
      
      ! Preconditioning of the defect: d=C^{-1}d
      CALL c2d2_precondDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rb, rd, &
          rtempvectorX, rtempvectorB, rtempVector, rnonlinearIteration, rnlSol)

      ! Add the defect: x = x + omega*d          
      CALL sptivec_vectorLinearComb (rd,rx,1.0_DP,1.0_DP)
          
      ! Assemble the new defect: d=b-Ax
      CALL sptivec_copyVector (rb,rd)
      CALL c2d2_assembleSpaceTimeDefect (rproblem, rspaceTimeDiscr, rx, rd, ddefNorm)
          
      CALL output_separator (OU_SEP_EQUAL)
      CALL output_line ('Iteration: '//sys_si(iglobIter,10)//&
          ' Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
      CALL output_separator (OU_SEP_EQUAL)
      
    END DO
    
    ! ---------------------------------------------------------------
    ! Release the preconditioner of the nonlinear iteration
    CALL c2d2_releasePreconditioner (rnonlinearIteration)
    
    !CALL c2d2_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
    CALL c2d2_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rd, &
      rtempvectorX, rtempvectorB, rtempvector,.FALSE.)

    ! Implement the initial condition into the RHS.
    CALL c2d2_implementInitCondRHS (rx, rd, rtempvectorX, rtempvector)    

    ! Assemble the defect
    CALL c2d2_assembleSpaceTimeDefect (rproblem, rspaceTimeDiscr, rx, rd, ddefNorm)
        
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    CALL output_separator (OU_SEP_EQUAL)        
    
    ! Release parameters of the nonlinear loop, final clean up
    CALL c2d2_doneNonlinearLoop (rnonlinearIteration)

    ! Do we have Neumann boundary?
    bneumann = collct_getvalue_int (rproblem%rcollection, 'INEUMANN') .EQ. YES
    bneumann = .FALSE.
    
    IF (.NOT. bneumann) THEN
      ! Normalise the primal and dual pressure to zero.
      DO isubstep = 0,rspaceTimeDiscr%niterations
      
        CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
        
        CALL vecfil_subvectorToL20 (rtempVectorX,3)
        CALL vecfil_subvectorToL20 (rtempVectorX,6)
        
        CALL sptivec_setTimestepData(rx, isubstep, rtempVectorX)
        
      END DO
      
    END IF
    
  END SUBROUTINE
  
END MODULE
