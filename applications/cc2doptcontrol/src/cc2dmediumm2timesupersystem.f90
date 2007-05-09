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
  TYPE t_ccoptSpaceTimeMatrix
  
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
    
    ! Regularisation parameter for the terminal condition 
    ! $\gamma/2*||y(T)-z(T)||$.
    ! A value of 0.0 disables the terminal condition.
    REAL(DP) :: dgammaC = 0.0_DP

    ! Time stepping scheme template structure that is used to initialise
    ! forward-in-time time stepping.
    TYPE(t_explicitTimeStepping) :: rtimeSteppingForward

    ! Time stepping scheme template structure that is used to initialise
    ! backward-in-time time stepping.
    TYPE(t_explicitTimeStepping) :: rtimeSteppingBackward
    
    ! Problem-related structure that provides the templates for
    ! matrices/vectors on the level of the matrix.
    TYPE(t_problem_lvl), POINTER :: p_rlevelInfo
    
  END TYPE

!</types>

CONTAINS

  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_initParamsSupersystem (rproblem,ilevelTime,ilevelSpace,&
      rsupersystem, rx, rd)
  
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
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(OUT) :: rsupersystem

  ! A space-time vector that is initialised for the current solution.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A space-time vector that is initialised for the current defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP) :: dtstep

    ! Copy most relevant data from the problem structure.
    rsupersystem%NLMAX           = ilevelSpace
    rsupersystem%p_rlevelInfo    => rproblem%RlevelInfo(ilevelSpace)
    rsupersystem%niterations     = &
        rproblem%rtimedependence%niterations * 2**MIN(0,ilevelTime-1)
    rsupersystem%dtimeInit       = rproblem%rtimedependence%dtimeInit
    rsupersystem%dtimeMax        = rproblem%rtimedependence%dtimeMax
    rsupersystem%ctimeStepScheme = rproblem%rtimedependence%ctimeStepScheme
    rsupersystem%dtimeStepTheta  = rproblem%rtimedependence%dtimeStepTheta
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dgammaC',rsupersystem%dgammaC,0.0_DP)

    ! The complete time interval is divided into niteration iterations
    dtstep = (rsupersystem%dtimemax-rsupersystem%dtimeInit) &
             / REAL(MAX(rsupersystem%niterations,1),DP)
             
    rsupersystem%dtstep = dtstep
    
    ! Initialise the time stepping in the problem structure.
    ! We have one scheme for the forward- and one scheme for the backward 
    ! iteration...
    CALL timstp_init (rsupersystem%rtimeSteppingForward, &
                      rsupersystem%ctimeStepScheme, rsupersystem%dtimeInit, &
                     dtstep, rsupersystem%dtimeStepTheta)
    CALL timstp_init (rsupersystem%rtimeSteppingBackward, &
                      rsupersystem%ctimeStepScheme, rsupersystem%dtimemax, &
                      -dtstep, rsupersystem%dtimeStepTheta)
    
    ! Initialise the global solution- and defect- vector.
    CALL sptivec_initVector (rx,&
        dof_igetNDofGlobBlock(rproblem%RlevelInfo(ilevelSpace)%p_rdiscretisation,.FALSE.),&
        rsupersystem%niterations)
    CALL sptivec_initVector (rd,&
        dof_igetNDofGlobBlock(rproblem%RlevelInfo(ilevelSpace)%p_rdiscretisation,.FALSE.),&
        rsupersystem%niterations)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_doneParamsSupersystem (rsupersystem,rx,rd)
  
!<description>
  ! Cleans up a given supersystem structure.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(OUT) :: rsupersystem

  ! A space-time vector that is initialised for the current solution.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A space-time vector that is initialised for the current defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! Release memory.
    CALL sptivec_releaseVector (rd)
    CALL sptivec_releaseVector (rx)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_solveSupersysDirect (rproblem, rsupermatrix, rx, rd, &
      rtempvectorX, rtempvectorB, rtempvectorD)

!<description>
  ! This routine assembles and solves the time-space coupled supersystem:
  ! $Ax=b$. The RHS vector is generated on-the-fly.
  ! The routine generates the full matrix in memory and solves with UMFPACK, 
  ! so it should only be used for debugging!
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A t_ccoptSpaceTimeMatrix structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rsupermatrix
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
    TYPE(t_matrixScalar) :: rmassLumped
    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    TYPE(t_vectorBlock) :: rxGlobal, rbGlobal, rdGlobal
    TYPE(t_vectorBlock) :: rxGlobalSolve, rbGlobalSolve, rdGlobalSolve
    TYPE(t_matrixBlock) :: rglobalA
    TYPE(t_linsolNode), POINTER :: rsolverNode
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: Isize
    INTEGER(PREC_VECIDX), DIMENSION(6) :: Isize2
    TYPE(t_vectorBlock) :: rvectorStructure
    
    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd

    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ilevel = rsupermatrix%NLMAX
    
    ! Calculate the lumped mass matrix of the FE space -- we need it later!
    CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass,&
        rmassLumped, LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
    CALL lsyssc_lumpMatrixScalar (rmassLumped,LSYSSC_LUMP_STD)
    
    ! ----------------------------------------------------------------------
    ! 1.) Generate the global RHS vector
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
    CALL lsysbl_getbase_double (rtempVectorB,p_Db)
    CALL lsysbl_getbase_double (rtempVectorD,p_Dd)

    DO isubstep = 0,rsupermatrix%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rsupermatrix%dtstep
          
      ! Generate the RHS of that point in time.
      CALL c2d2_generateBasicRHS (rproblem,rtempVectorB)
      
      ! Multiply the RHS of the dual velocity by dtstep according to the time
      ! discretisation -- except for if we are in the last timestep!
      !
      ! In the last timestep, if GAMMA is =0, we have no terminal condition
      ! and thus have to force the RHS to 0!
      ! Otherwise, the terminat condition is multiplied with dgammaC.
      IF (isubstep .NE. rsupermatrix%niterations) THEN
        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(4),rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(5),rsupermatrix%dtstep)
      ELSE
        ! Multiply -z by gamma, that's it.
        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(4),rsupermatrix%dgammaC)
        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(5),rsupermatrix%dgammaC)
!      ELSE
!        IF (rsupermatrix%dgammaC .EQ. 0.0_DP) THEN
!          CALL lsyssc_clearVector (rtempVectorB%RvectorBlock(4))
!          CALL lsyssc_clearVector (rtempVectorB%RvectorBlock(5))
!        ELSE
!        
!          ! Load y(T) into rtempVectorX
!          !CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
!
!          ! Calculate GAMMA*(-z(T)). For that purpose, calculate z(T)
!          ! using the L2 projection of z into rtempVectorX(4,5).
!          CALL l2prj_analytL2projectionByMass (rtempVectorB%RvectorBlock(4),&
!              rsupermatrix%p_rlevelInfo%rmatrixMass,rmassLumped,&
!              rtempVectorX%RvectorBlock(4),&
!              rtempVectorD%RvectorBlock(4),coeff_TARGET_x,&
!              rproblem%rcollection)
!          CALL l2prj_analytL2projectionByMass (rtempVectorB%RvectorBlock(5),&
!              rsupermatrix%p_rlevelInfo%rmatrixMass,rmassLumped,&
!              rtempVectorX%RvectorBlock(5),&
!              rtempVectorD%RvectorBlock(5),coeff_TARGET_y,&
!              rproblem%rcollection)
!             
!        END IF
      END IF
      
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
      CALL vecfil_discreteBCsol (rtempVectorB)
      CALL vecfil_discreteFBCsol (rtempVectorB)      
      
      CALL sptivec_setTimestepData(rd, isubstep, rtempVectorB)
      
      ! Clean up the collection (as we are done with the assembly, that's it.
      CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

    END DO

    ! Release the mass matr, we don't need it anymore
    CALL lsyssc_releaseMatrix (rmassLumped)
    
    ! ----------------------------------------------------------------------
    ! 2.) Generate the matrix A
    !
    ! Create a global matrix:
    CALL lsysbl_createEmptyMatrix (rglobalA,6*(rsupermatrix%niterations+1))
    
    ! Loop through the substeps
    
    DO isubstep = 0,rsupermatrix%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rsupermatrix%dtstep

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
      
        ! Create a matrix that applies "-M" to the dual velocity and include it
        ! to the global matrix.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockTemp)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        rblockTemp%RmatrixBlock(4,4)%dscaleFactor = -1.0_DP
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        rblockTemp%RmatrixBlock(5,5)%dscaleFactor = -1.0_DP
        
        ! Include the boundary conditions into that matrix.
        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
        ! main diagonal of the supermatrix.
        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
        
        ! Include "-M" in the global matrix at position (1,2).
        CALL insertMatrix (rblockTemp,rglobalA,7,1)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)
      
        ! -----
        
        ! Now the hardest -- or longest -- part: The diagonal matrix.
        !
        ! Generate the basic system matrix level rsupermatrix%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rsupermatrix%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
          1.0_DP, &
          rsupermatrix%dtstep, &
          rsupermatrix%dtstep * REAL(1-rproblem%iequation,DP))
      
        ! The primal equation consists of identity matrices.
        rnonlinearIterationTmp%bprimalIdentity = .TRUE.

        ! Assemble the system matrix on level rsupermatrix%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
            
        ! Implement boundary conditions into the matrix:
        ! standard boundary conditions        
        CALL matfil_discreteBC (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)
        ! fictitious boundary boundary conditions
        CALL matfil_discreteFBC (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)
            
        ! Replace the primal equation by identity matrices.
        CALL lsyssc_initialiseIdentityMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix%RmatrixBlock(1,1))
        CALL lsyssc_initialiseIdentityMatrix (&
            rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix%RmatrixBlock(2,2))
        
        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rglobalA,1,1)
            
        ! Create an identity matrix for the pressure
        CALL lsyssc_createDiagMatrixStruc (rglobalA%RmatrixBlock(3,3),&
            rtempVectorB%RvectorBlock(3)%NEQ,LSYSSC_MATRIX9)
        CALL lsyssc_initialiseIdentityMatrix (&
            rglobalA%RmatrixBlock(3,3))
            
        ! Scale the B-matrices of the dual pressure.
        rglobalA%RmatrixBlock (4:5,3)%dscaleFactor = rsupermatrix%dtstep
        
        ! Switch off the B/M-matrices of the primal equation
        rglobalA%RmatrixBlock (3,1:2)%dscaleFactor = 0.0_DP
        rglobalA%RmatrixBlock (1:2,3:5)%dscaleFactor = 0.0_DP

      ELSE IF (isubstep .EQ. rsupermatrix%niterations) THEN
        
        ! We are in the last substep
        
        ! -----
        
        ! Create a matrix that applies "-M" to the primal velocity and include it
        ! to the global matrix.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockTemp)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        rblockTemp%RmatrixBlock(1,1)%dscaleFactor = -1.0_DP
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        rblockTemp%RmatrixBlock(2,2)%dscaleFactor = -1.0_DP
        
        ! Include the boundary conditions into that matrix.
        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
        ! main diagonal of the supermatrix.
        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)

        ! Include "-M" in the global matrix at position (2,1).
        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)
      
        ! -----
        
        ! Now the hardest -- or longest -- part: The diagonal matrix.
        !
        ! Generate the basic system matrix level rsupermatrix%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rsupermatrix%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
          1.0_DP, &
          rsupermatrix%dtstep, &
          rsupermatrix%dtstep * REAL(1-rproblem%iequation,DP))
      
        ! The dual equation consists of identity matrices.
        rnonlinearIterationTmp%bdualTerminal = .TRUE.

        ! Assemble the system matrix on level rsupermatrix%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        
        ! Implement boundary conditions into the matrix:
        ! standard boundary conditions        
        CALL matfil_discreteBC (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)
        ! fictitious boundary boundary conditions
        CALL matfil_discreteFBC (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)

        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rglobalA,isubstep*6+1,isubstep*6+1)

        ! Create an identity matrix for the dual pressure
        CALL lsyssc_createDiagMatrixStruc (rglobalA%RmatrixBlock(isubstep*6+6,isubstep*6+6),&
            rtempVectorB%RvectorBlock(3)%NEQ,LSYSSC_MATRIX9)
        CALL lsyssc_initialiseIdentityMatrix (&
            rglobalA%RmatrixBlock(isubstep*6+6,isubstep*6+6))

        ! Scale the B-matrices of the primal pressure.
        rglobalA%RmatrixBlock (isubstep*6+1:isubstep*6+2,isubstep*6+3)%dscaleFactor = &
            rsupermatrix%dtstep

        ! Switch off the B-matrices of the dual equation
        rglobalA%RmatrixBlock (isubstep*6+6,isubstep*6+4:isubstep*6+5)%dscaleFactor = 0.0_DP
        rglobalA%RmatrixBlock (isubstep*6+4:isubstep*6+5,isubstep*6+6)%dscaleFactor = 0.0_DP
        
        ! Switch off / Scale the M-matrices correctly.
        rglobalA%RmatrixBlock (isubstep*6+4,isubstep*6+1)%dscaleFactor = &
            -rsupermatrix%dgammaC
        rglobalA%RmatrixBlock (isubstep*6+5,isubstep*6+2)%dscaleFactor = &
            -rsupermatrix%dgammaC

        ! -----
        ! Create a matrix that applies "M" to the dual velocity and include it
        ! to the global matrix.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockTemp)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
        ! Include the boundary conditions into that matrix.
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
        
        ! Include that in the global matrix below the diagonal
        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1,isubstep*6+1)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)
      
      ELSE
      
        ! We are sonewhere in the middle of the matrix. There is a substep
        ! isubstep+1 and a substep isubstep-1!
        
        ! -----
        
        ! Create a matrix that applies "-M" to the dual velocity and include it
        ! to the global matrix.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockTemp)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        rblockTemp%RmatrixBlock(4,4)%dscaleFactor = -1.0_DP
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        rblockTemp%RmatrixBlock(5,5)%dscaleFactor = -1.0_DP
        
        ! Include the boundary conditions into that matrix.
        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
        ! main diagonal of the supermatrix.
        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)

        ! Include that in the global matrix below the diagonal
        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+7,isubstep*6+1)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)
      
        ! -----
        
        ! Create a matrix that applies "-M" to the primal velocity and include it
        ! to the global matrix.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockTemp)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        rblockTemp%RmatrixBlock(1,1)%dscaleFactor = -1.0_DP
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockTemp%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        rblockTemp%RmatrixBlock(2,2)%dscaleFactor = -1.0_DP
        
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

        ! -----      

        ! Now the hardest -- or longest -- part: The diagonal matrix.
        !
        ! Generate the basic system matrix level rsupermatrix%NLMAX
        ! Will be modified by c2d2_assembleLinearisedMatrices later.
        CALL c2d2_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
            rproblem%RlevelInfo(ilevel)%rmatrix)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rsupermatrix%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
          1.0_DP, &
          rsupermatrix%dtstep, &
          rsupermatrix%dtstep * REAL(1-rproblem%iequation,DP))
      
        ! Assemble the system matrix on level rsupermatrix%NLMAX.
        ! Include the boundary conditions into the matrices.
        CALL c2d2_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        
        ! Implement boundary conditions into the matrix:
        ! standard boundary conditions        
        CALL matfil_discreteBC (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)
        ! fictitious boundary boundary conditions
        CALL matfil_discreteFBC (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix)

        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
            rglobalA,isubstep*6+1,isubstep*6+1)
            
        ! Scale the B-matrices of the primal and dual pressure.
        rglobalA%RmatrixBlock (isubstep*6+1:isubstep*6+2,isubstep*6+3)%dscaleFactor = &
            rsupermatrix%dtstep
        rglobalA%RmatrixBlock (isubstep*6+4:isubstep*6+5,isubstep*6+6)%dscaleFactor = &
            rsupermatrix%dtstep

      END IF
    
    END DO
    
    CALL lsysbl_updateMatStrucInfo (rglobalA)

    ! Write the global matrix to a file.
    !CALL matio_writeBlockMatrixHR(rglobalA,'MATRIX',.TRUE.,0,'matrix.txt','(E10.2)')
    
    ! Get the global solution/rhs/temp vector.
    CALL sptivec_convertSupervecToVector (rx, rxGlobal)
    CALL sptivec_convertSupervecToVector (rd, rbGlobal)
    CALL lsysbl_createVecBlockIndirect (rbGlobal,rdGlobal,.FALSE.)
    
    ! Initialise the UMFPACK solver.
    CALL linsol_initUMFPACK4 (rsolverNode)
    
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

  SUBROUTINE c2d2_assembleDefectSupersystem (rproblem, rsupermatrix, rx, rd, &
      rtempvectorX, rtempvectorB, rtempvectorD, ddefNorm)

!<description>
  ! This routine assembles the defect of the time-space coupled supersystem:
  ! $d=b-Ax$ with d=rd, x=rx, A=rsupermatrix. The RHS vector is generated
  ! on-the-fly.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A t_ccoptSpaceTimeMatrix structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rsupermatrix
  
  ! A space-time vector defining the current solution.
  TYPE(t_spacetimeVector), INTENT(IN) :: rx
!</input>

!<inputoutput>
  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX

  ! A second temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorB

  ! A third temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD

  ! A space-time vector that receives the defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!<output>
  ! The norm of the total defect.
  REAL(DP), INTENT(OUT) :: ddefNorm
!</output>

!</subroutine>

    ! local variables
    INTEGER :: isubstep
    TYPE(t_matrixBlock) :: rblockMass
    TYPE(t_matrixScalar) :: rmassLumped
    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    
    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd

    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ! Calculate the lumped mass matrix of the FE space -- we need it later!
    CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass,&
        rmassLumped, LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
    CALL lsyssc_lumpMatrixScalar (rmassLumped,LSYSSC_LUMP_STD)
    
    ! ----------------------------------------------------------------------
    ! 1.) Generate the global RHS vector
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
    CALL lsysbl_getbase_double (rtempVectorB,p_Db)
    CALL lsysbl_getbase_double (rtempVectorD,p_Dd)

    DO isubstep = 0,rsupermatrix%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rsupermatrix%dtstep

      ! Initialise the collection for the assembly process with callback routines.
      ! Basically, this stores the simulation time in the collection if the
      ! simulation is nonstationary.
      CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

      ! Generate the RHS of that point in time.
      CALL c2d2_generateBasicRHS (rproblem,rtempVectorB)
      
      ! Multiply the RHS of the dual velocity by dtstep according to the time
      ! discretisation -- except for if we are in the last timestep!
      !
      ! In the last timestep, if GAMMA is =0, we have no terminal condition
      ! and thus have to force the RHS to 0!
      IF (isubstep .NE. rsupermatrix%niterations) THEN
        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(4),rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(5),rsupermatrix%dtstep)
      ELSE
        IF (rsupermatrix%dgammaC .EQ. 0.0_DP) THEN
          CALL lsyssc_clearVector (rtempVectorB%RvectorBlock(4))
          CALL lsyssc_clearVector (rtempVectorB%RvectorBlock(5))
        ELSE
        
          ! Load y(T) into rtempVectorX
          CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)

          ! Calculate GAMMA*(y(T)-z(T)). For that purpose, calculate z(T)
          ! using the L2 projection of z into rtempVectorX(4,5).
          CALL l2prj_analytL2projectionByMass (rtempVectorB%RvectorBlock(4),&
              rsupermatrix%p_rlevelInfo%rmatrixMass,rmassLumped,&
              rtempVectorX%RvectorBlock(4),&
              rtempVectorD%RvectorBlock(4),coeff_TARGET_x,&
              rproblem%rcollection)
          CALL l2prj_analytL2projectionByMass (rtempVectorB%RvectorBlock(4),&
              rsupermatrix%p_rlevelInfo%rmatrixMass,rmassLumped,&
              rtempVectorX%RvectorBlock(5),&
              rtempVectorD%RvectorBlock(5),coeff_TARGET_y,&
              rproblem%rcollection)
              
          ! Subtract z from y to get the terminal condition. This is later
          ! transferred to the solution vector without change.
          CALL lsyssc_vectorLinearComb (&
              rtempVectorX%RvectorBlock(1),rtempVectorX%RvectorBlock(4),&
              1.0_DP*rsupermatrix%dgammaC,-1.0_DP*rsupermatrix%dgammaC)
          CALL lsyssc_vectorLinearComb (&
              rtempVectorX%RvectorBlock(2),rtempVectorX%RvectorBlock(5),&
              1.0_DP*rsupermatrix%dgammaC,-1.0_DP*rsupermatrix%dgammaC)
          CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(4),&
              rtempVectorB%RvectorBlock(4))
          CALL lsyssc_copyVector (rtempVectorX%RvectorBlock(5),&
              rtempVectorB%RvectorBlock(5))
        END IF
      END IF
      
      CALL sptivec_setTimestepData(rd, isubstep, rtempVectorB)
      
      ! Clean up the collection (as we are done with the assembly, that's it.
      CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

    END DO
    
    ! ----------------------------------------------------------------------
    ! 2.) Generate rb = rb-Ax
    !
    ! For this purpose, loop through the substeps
    
    ddefNorm = 0.0_DP
    
    DO isubstep = 0,rsupermatrix%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rsupermatrix%dtstep

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! Read in the RHS from the current timestep.
      CALL sptivec_getTimestepData (rd, isubstep, rtempVectorB)
        
      ! The first and last substep is a little bit special concerning
      ! the matrix!
      IF (isubstep .EQ. 0) THEN
        
        ! We are in the first substep
      
        ! -----
      
        ! Read in the solution vector at time $t^{n+1}$
        CALL sptivec_getTimestepData (rx, isubstep+1, rtempVectorX)
        
        ! Multiply the dual velocity (not the primal one!) by "-M" and subtract
        ! from the given RHS. For this purpose, create a block matrix with
        ! "M" on the dual velocity diagonal.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockMass)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
        ! b=b+Mx = b-(-M)x
        CALL lsysbl_blockMatVec (rblockMass,rtempVectorX,rtempVectorB,1.0_DP,1.0_DP)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockMass)
      
        ! -----
        
        ! Now the hardest -- or longest -- part: The nonlinear defect.  
        !
        ! Read in the solution vector at time $t^{n}$
        CALL sptivec_getTimestepData (rx, isubstep, rtempVectorX)

        ! Scale the primal and dual pressure by the step length. this implements
        ! the 'delta*t' in front of these subvectors
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
            rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
            rsupermatrix%dtstep)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rsupermatrix%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
          1.0_DP, &
          rsupermatrix%dtstep, &
          rsupermatrix%dtstep * REAL(1-rproblem%iequation,DP))
      
        ! The primal equation consists of identity matrices.
        ! Copy the RHS to the velocity vectors -- it's the initial condition.
        rnonlinearIterationTmp%bprimalIdentity = .TRUE.
        CALL lsyssc_copyVector (rtempVectorB%RvectorBlock(1),&
            rtempVectorX%RvectorBlock(1))
        CALL lsyssc_copyVector (rtempVectorB%RvectorBlock(2),&
            rtempVectorX%RvectorBlock(2))
        CALL lsyssc_copyVector (rtempVectorB%RvectorBlock(3),&
            rtempVectorX%RvectorBlock(3))

        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
            rtempVectorX,rtempVectorB,&
            .TRUE.,.TRUE.,rproblem%rcollection)
            
        ! Scale the primal and dual pressure back. We only scaled it for setting
        ! up the defect.
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
            1.0_DP/rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
            1.0_DP/rsupermatrix%dtstep)
            
        ! Overwrite the defect of the primal vectors by zero -- since this
        ! is the boundary of the time cylinder!
        CALL lsyssc_clearVector (rtempVectorB%RvectorBlock(1))
        CALL lsyssc_clearVector (rtempVectorB%RvectorBlock(2))

      ELSE IF (isubstep .EQ. rsupermatrix%niterations) THEN
        
        ! We are in the last substep
        
        ! -----
        
        ! Read in the solution vector at time $t^{n-1}$
        CALL sptivec_getTimestepData (rx, isubstep-1, rtempVectorX)
        
        ! Multiply the velocity (not the dual one!) by "-M" and subtract
        ! from the given RHS. For this purpose, create a block matrix with
        ! "M" on the velocity diagonal.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockMass)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
        ! b=b+Mx = b-(-M)x
        CALL lsysbl_blockMatVec (rblockMass,rtempVectorX,rtempVectorB,1.0_DP,1.0_DP)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockMass)

        ! -----
        
        ! Now the hardest -- or longest -- part: The nonlinear defect.  
        !
        ! Read in the solution vector at time $t^{n}$
        CALL sptivec_getTimestepData (rx, isubstep, rtempVectorX)

        ! If GAMMA <> 0, there is a terminal condition. In this case, we have to 
        ! subtract -GAMMA*M*lambda(last timestep) from the RHS.
        IF (rsupermatrix%dgammaC .NE. 0.0_DP) THEN
          CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
              rblockMass)
          CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
              rblockMass%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
              rblockMass%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
          ! b=b+GAMMA*Mx = b-(-GAMMA*M)x
          CALL lsysbl_blockMatVec (rblockMass,rtempVectorX,rtempVectorB,&
              rsupermatrix%dgammaC,1.0_DP)
          
          ! Release the block mass matrix.
          CALL lsysbl_releaseMatrix (rblockMass)
        END IF
        
        ! Scale the primal and dual pressure by the step length. this implements
        ! the 'delta*t' in front of these subvectors
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
            rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
            rsupermatrix%dtstep)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rsupermatrix%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')
            
        CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
          1.0_DP, &
          rsupermatrix%dtstep, &
          rsupermatrix%dtstep * REAL(1-rproblem%iequation,DP))
          
        ! The dual equation consists of identity matrices.
        ! Copy the dual velocity from the RHS to the solution; it's the terminal
        ! condition!
        rnonlinearIterationTmp%bdualTerminal = .TRUE.
        CALL lsyssc_copyVector (rtempVectorB%RvectorBlock(4),&
            rtempVectorX%RvectorBlock(4))
        CALL lsyssc_copyVector (rtempVectorB%RvectorBlock(5),&
            rtempVectorX%RvectorBlock(5))
        CALL lsyssc_copyVector (rtempVectorB%RvectorBlock(6),&
            rtempVectorX%RvectorBlock(6))

        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
            rtempVectorX,rtempVectorB,&
            .TRUE.,.TRUE.,rproblem%rcollection)
            
        ! Scale the primal and dual pressure back. We only scaled it for setting
        ! up the defect.
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
            1.0_DP/rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
            1.0_DP/rsupermatrix%dtstep)

      ELSE
      
        ! We are sonewhere in the middle of the matrix. There is a substep
        ! isubstep+1 and a substep isubstep-1!
        
        ! -----
        
        ! Read in the solution vector at time $t^{n-1}$
        CALL sptivec_getTimestepData (rx, isubstep-1, rtempVectorX)
        
        ! Multiply the velocity (not the dual one!) by "-M" and subtract
        ! from the given RHS. For this purpose, create a block matrix with
        ! "M" on the velocity diagonal.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockMass)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
        ! b=b+Mx = b-(-M)x
        CALL lsysbl_blockMatVec (rblockMass,rtempVectorX,rtempVectorB,1.0_DP,1.0_DP)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockMass)
        
        ! -----
      
        ! Read in the solution vector at time $t^{n+1}$
        CALL sptivec_getTimestepData (rx, isubstep+1, rtempVectorX)
        
        ! Multiply the dual velocity (not the dual one!) by "-M" and subtract
        ! from the given RHS. For this purpose, create a block matrix with
        ! "M" on the dual velocity diagonal.
        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
            rblockMass)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        CALL lsyssc_duplicateMatrix (rsupermatrix%p_rlevelInfo%rmatrixMass, &
            rblockMass%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
        ! b=b+Mx = b-(-M)x
        CALL lsysbl_blockMatVec (rblockMass,rtempVectorX,rtempVectorB,1.0_DP,1.0_DP)
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockMass)
      
        ! -----
        
        ! Now the hardest -- or longest -- part: The nonlinear defect.  
        !
        ! Read in the solution vector at time $t^{n}$
        CALL sptivec_getTimestepData (rx, isubstep, rtempVectorX)

        ! Scale the primal and dual pressure by the step length. this implements
        ! the 'delta*t' in front of these subvectors
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
            rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
            rsupermatrix%dtstep)
      
        ! Set up a core equation structure and assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        CALL c2d2_initNonlinearLoop (&
            rproblem,rproblem%NLMIN,rsupermatrix%NLMAX,rtempVectorX,rtempVectorB,&
            rnonlinearIterationTmp,'CC2D-NONLINEAR')

        CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
          1.0_DP, &
          rsupermatrix%dtstep, &
          rsupermatrix%dtstep * REAL(1-rproblem%iequation,DP))
      
        CALL c2d2_assembleNonlinearDefect (rnonlinearIterationTmp,&
            rtempVectorX,rtempVectorB,&
            .TRUE.,.TRUE.,rproblem%rcollection)
            
        ! Scale the primal and dual pressure back. We only scaled it for setting
        ! up the defect.
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
            1.0_DP/rsupermatrix%dtstep)
        CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
            1.0_DP/rsupermatrix%dtstep)
        
      END IF
    
      ! Save the subdefect of the current timestep.
      CALL sptivec_setTimestepData (rd, isubstep, rtempVectorB)
      
      ! Calculate the norm of the subdefect
      ddefNorm = ddefNorm + lsysbl_vectorNorm (rtempVectorB,LINALG_NORML2)**2
    
    END DO
    
    ! Complete the calculation of the defect norm by dividing by the
    ! number of subvectors and taking the square root.
    ddefNorm = SQRT(ddefNorm/rsupermatrix%niterations)

    ! Release the lumped mass matrix, that's it.
    CALL lsyssc_releaseMatrix (rmassLumped)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_precondDefectSupersystem (rproblem, rsupermatrix, rx, rd, &
      rtempvectorX, rtempvectorB, rtempvector, rnonlinearIteration, rnlSolver)

!<description>
  ! This routine performs preconditioning with the nonlinear super-defect
  ! vector rd: $d = C^{-1} d$.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A t_ccoptSpaceTimeMatrix structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rsupermatrix
  
  ! A space-time vector defining the current solution.
  TYPE(t_spacetimeVector), INTENT(IN) :: rx

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
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector
  
  ! A space-time vector that receives the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubstep
    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    RETURN
    ! ----------------------------------------------------------------------
    ! We use a block-Jacobi scheme for preconditioning...
    !
    ! For this purpose, loop through the substeps.
    
    DO isubstep = 0,rsupermatrix%niterations
    
      ! Current time step?
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*rsupermatrix%dtstep

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
      
      ! Read in the solution- and defect vector of the current timestep.
      CALL sptivec_getTimestepData (rx, isubstep, rtempVectorX)
      CALL sptivec_getTimestepData (rd, isubstep, rtempVectorB)

      ! Scale the primal and dual pressure by the step length. this implements
      ! the 'delta*t' in front of these subvectors
      ! The core equation routines handle the equation
      !   alpha*M*u + theta*nu*Laplace*u + gamma*N(u)u + B*p = ...
      ! but we want to solve
      !   alpha*M*u + theta*nu*Laplace*u + gamma*N(u)u + tstep*B*p = ...
      !
      ! So the trick is to scale p by tstep, solve the core equation
      !   alpha*M*u + theta*nu*Laplace*u + gamma*N(u)u + B*(tstep*p) = ...
      ! and scale it back afterwards.
      CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
          rsupermatrix%dtstep)
      CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
          rsupermatrix%dtstep)
    
      ! Setup the core equation for the nonlinear loop.      
      rnonlinearIterationTmp = rnonlinearIteration

      CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
        1.0_DP, &
        rsupermatrix%dtstep, &
        rsupermatrix%dtstep * REAL(1-rproblem%iequation,DP))

      ! Call the solver of the core equation to solve it using a nonlinear
      ! iteration. Overwrite rtempVectorX with temporary data from the
      ! nonlinear iteration.
      CALL c2d2_solveCoreEquation (rnlSolver,rnonlinearIterationTmp,&
          rtempVectorX,rtempVectorB,rproblem%rcollection,rtempVector)             
    
      ! Scale the pressure back, then we have again the correct solution vector.
      CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(NDIM2D+1),&
          1.0_DP/rsupermatrix%dtstep)
      CALL lsyssc_scaleVector (rtempVectorX%RvectorBlock(2*(NDIM2D+1)),&
          1.0_DP/rsupermatrix%dtstep)
    
      ! Save back the preconditioned defect.
      CALL sptivec_setTimestepData (rd, isubstep, rtempVectorB)
      
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_solveSupersystem (rproblem, rsupermatrix, rx, rd, &
      rtempvectorX, rtempvectorB, rtempVector)
  
!<description>
  ! 
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A t_ccoptSpaceTimeMatrix structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rsupermatrix
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

  ! A temporary space-time vector that receives the defect during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! The nonlinear solver configuration
    TYPE(t_nlsolNode) :: rnlSol
    TYPE(t_ccNonlinearIteration) :: rnonlinearIteration
    INTEGER :: isubstep,iglobIter

    REAL(DP) :: ddefNorm,dinitDefNorm

    ! Some preparations for the nonlinear solver.
    !
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    CALL c2d2_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! or callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    CALL c2d2_initNonlinearLoop (rproblem,rproblem%NLMIN,rsupermatrix%NLMAX,&
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

    DO isubstep = 1,rsupermatrix%niterations
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + (isubstep-1)*rsupermatrix%dtstep

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! Implement the boundary conditions into the global solution vector.
      CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
      
      CALL c2d2_implementBC (rproblem,rtempVectorX,rtempVectorB,.FALSE.,&
          .TRUE.,.FALSE.)
      
      CALL sptivec_setTimestepData(rx, isubstep, rtempVectorX)
      
    END DO
    
    ddefNorm = 1.0_DP
    
    ! ---------------------------------------------------------------
    ! Solve the global space-time coupled system.
    !
    ! Get the initial defect: d=b-Ax
    CALL c2d2_assembleDefectSupersystem (rproblem, rsupermatrix, rx, rd, &
        rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
        
    dinitDefNorm = ddefNorm
    
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    CALL output_separator (OU_SEP_EQUAL)        
    
    iglobIter = 0
    
    CALL c2d2_solveSupersysDirect (rproblem, rsupermatrix, rx, rd, &
      rtempvectorX, rtempvectorB, rtempVector)
    
!    DO WHILE ((ddefNorm .GT. 1.0E-2*dinitDefNorm) .AND. (ddefNorm .LT. 1.0E99) .AND. &
!              (iglobIter .LT. 10))
!    
!      iglobIter = iglobIter+1
!      
!      ! Preconditioning of the defect: d=C^{-1}d
!      CALL c2d2_precondDefectSupersystem (rproblem, rsupermatrix, rx, rd, &
!          rtempvectorX, rtempvectorB, rtempVector, rnonlinearIteration, rnlSol)
!
!      ! Add the defect: x = x + omega*d          
!      CALL sptivec_vectorLinearComb (rd,rx,0.05_DP,1.0_DP)
!          
!      ! Assemble the new defect: d=b-Ax
!      CALL c2d2_assembleDefectSupersystem (rproblem, rsupermatrix, rx, rd, &
!          rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
!          
!      CALL output_separator (OU_SEP_EQUAL)
!      CALL output_line ('Iteration: '//sys_si(iglobIter,10)//&
!          ' Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
!      CALL output_separator (OU_SEP_EQUAL)
!      
!    END DO
    
    ! ---------------------------------------------------------------
    ! Release the preconditioner of the nonlinear iteration
    CALL c2d2_releasePreconditioner (rnonlinearIteration)
    
    ! Release parameters of the nonlinear loop, final clean up
    CALL c2d2_doneNonlinearLoop (rnonlinearIteration)

  END SUBROUTINE
  
END MODULE
