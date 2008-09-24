!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2oneshot </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module aims to solve the fully space-time coupled Stokes and
!# Navier-Stokes optimal control problem with a one-shot approach.
!# As this type of problem is of elliptic nature with a very special
!# structure, a nonlinear defect correction loop is executed for the solution
!# where a space-time coupled multigrid preconditioner does the preconditioning
!# of the nonlinear defect.
!#
!# The following subroutines can be found here:
!#
!# 1.) cc_solveSupersysDirectCN
!#     -> Only for nonstationary Stokes problem. Executes Crank-Nicolson on the
!#        fully space-time coupled system
!#
!# 2.) cc_solveSupersystemMultigrid
!#     -> Solves the Stokes and Navier-Stokes problem with a preconditioned
!#        defect correction approach. Crank Nicolson is used for the time
!#        discretisation.
!#
!# Auxiliary routines:
!#
!# 1.) cc_assembleSpaceTimeDefect
!#     -> Calculates a space-time defect
!#
!# 2.) cc_precondDefectSupersystem
!#     -> Executes preconditioning on a space-time defect vector
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2oneshot

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
  USE statistics
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback

  USE cc2dmediumm2nonlinearcore
  USE cc2dmediumm2nonlinearcoreinit
  USE cc2dmediumm2stationary
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  USE cc2dmediumm2matvecassembly
  USE spacetimediscretisation
  USE cc2dmediumm2spacetimesolver
  
  USE cc2dmediumm2scriptfile
  
  USE spacetimevectors
  USE timerhsevaluation
  USE dofmapping
  
  USE matrixio
    
  IMPLICIT NONE

CONTAINS

!  ! ***************************************************************************
!  
!!<subroutine>
!
!  SUBROUTINE cc_solveSupersysDirect (rproblem, rspaceTimeDiscr, rx, rd, &
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
!    ilevel = rspaceTimeDiscr%ilevel
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
!      CALL cc_generateBasicRHS (rproblem,rtempVectorB)
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
!      CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)
!
!      ! Discretise the boundary conditions at the new point in time -- 
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL cc_updateDiscreteBC (rproblem)
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
!      CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
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
!        CALL cc_updateDiscreteBC (rproblem)
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
!        ! Generate the basic system matrix level rspaceTimeDiscr%ilevel
!        ! Will be modified by cc_assembleLinearisedMatrices later.
!        CALL cc_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!      
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!      
!        CALL cc_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%ilevel,rtempVectorX,rtempVectorB,&
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
!        ! Assemble the system matrix on level rspaceTimeDiscr%ilevel.
!        ! Include the boundary conditions into the matrices.
!        CALL cc_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
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
!        ! Generate the basic system matrix level rspaceTimeDiscr%ilevel
!        ! Will be modified by cc_assembleLinearisedMatrices later.
!        CALL cc_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!      
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!      
!        CALL cc_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%ilevel,rtempVectorX,rtempVectorB,&
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
!        ! Assemble the system matrix on level rspaceTimeDiscr%ilevel.
!        ! Include the boundary conditions into the matrices.
!        CALL cc_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
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
!        ! Generate the basic system matrix level rspaceTimeDiscr%ilevel
!        ! Will be modified by cc_assembleLinearisedMatrices later.
!        CALL cc_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!      
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!      
!        CALL cc_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%ilevel,rtempVectorX,rtempVectorB,&
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
!        ! Assemble the system matrix on level rspaceTimeDiscr%ilevel.
!        ! Include the boundary conditions into the matrices.
!        CALL cc_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
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

  SUBROUTINE cc_applyUmfpackToSupersystem (rproblem, rspaceTimeMatrix, rx, rd, &
      rtempvectorX, rtempvectorB, rtempvectorD)

!<description>
  ! This routine assembles and solves the time-space coupled supersystem:
  ! $Ax=b$. The RHS vector is generated on-the-fly.
  ! The routine generates the full matrix in memory and solves with UMFPACK, 
  ! so it should only be used for debugging!
  ! The routine assembles the global matrix with the time step scheme
  ! specified in dtimeStepTheta in the problem structure.
  !
  ! rspaceTimeMatrix specifies the global space time matrix of the system.
  ! In this routine, there is no processing of any nonlinearity, so only
  ! linear problems can be solved (e.g. Stokes).
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! The definition of the global space time matrix that should be used
  ! for solving the system.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN),TARGET :: rspaceTimeMatrix
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
    TYPE(t_vectorBlock) :: rxGlobal, rbGlobal, rdGlobal
    TYPE(t_vectorBlock) :: rxGlobalSolve, rbGlobalSolve, rdGlobalSolve
    TYPE(t_matrixBlock) :: rglobalA
    TYPE(t_linsolNode), POINTER :: rsolverNode,p_rpreconditioner
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: Isize
    INTEGER(PREC_VECIDX), DIMENSION(6) :: Isize2
    REAL(DP) :: dtheta
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_matrixBlock) :: rmatrix
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    TYPE(t_vectorBlock), DIMENSION(3) :: rtempVectorSol
    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol
    
    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd

    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscretisation
    
    ilevel = p_rspaceTimeDiscr%ilevel
    
    ! Theta-scheme identifier.
    ! =1: impliciz Euler.
    ! =0.5: Crank Nicolson
    dtheta = rproblem%rtimedependence%dtimeStepTheta
    
    ! Generate the RHS for the initial condition.
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rinitialCondRHS,.false.)
    call sptivec_getTimestepData (rx, 1, rtempVectorX)
    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
    call cc_generateInitCondRHS (rproblem,p_rspaceTimeDiscr,&
        rtempVectorX,rinitialCondRHS)

    ! Assemble the space-time RHS into rd.
    CALL trhsevl_assembleRHS (rproblem, p_rspaceTimeDiscr, rd, .TRUE.)
      
    ! Implement the initial condition into the RHS/solution.
    CALL tbc_implementInitCond (rproblem, rd, rinitialCondRHS, rtempvectorD)    
    CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvectorD)    

    ! Release the rhs vector with the init. condition again.
    call lsysbl_releaseVector (rinitialCondRHS)
    call lsysbl_releaseVector (rinitialCondSol)

    ! ----------------------------------------------------------------------
    ! 2.) Generate the matrix A
    !
    ! Create a global matrix:
    CALL lsysbl_createEmptyMatrix (rglobalA,6*(p_rspaceTimeDiscr%rtimeDiscr%nintervals+1))

    ! Basic initialisation of rmatrixComponents with the pointers to the
    ! matrices / discretisation structures on the current level.
    !
    ! The weights in the rmatrixComponents structure are later initialised
    ! according to the actual situation when the matrix is to be used.
    rmatrixComponents%p_rdiscretisation         => &
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation
    rmatrixComponents%p_rmatrixStokes           => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixStokes          
    rmatrixComponents%p_rmatrixB1             => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB1              
    rmatrixComponents%p_rmatrixB2             => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB2              
    rmatrixComponents%p_rmatrixMass           => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixMass            
    rmatrixComponents%p_rmatrixIdentityPressure => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixIdentityPressure
    rmatrixComponents%dnu = collct_getvalue_real (rproblem%rcollection,'NU')
    rmatrixComponents%iupwind1 = collct_getvalue_int (rproblem%rcollection,'IUPWIND1')
    rmatrixComponents%dupsam1 = collct_getvalue_real (rproblem%rcollection,'UPSAM1')
    rmatrixComponents%iupwind2 = collct_getvalue_int (rproblem%rcollection,'IUPWIND2')
    rmatrixComponents%dupsam2 = collct_getvalue_real (rproblem%rcollection,'UPSAM2')

    ! Get a temporary system matrix
    CALL cc_allocSystemMatrix (rproblem,rproblem%RlevelInfo(ilevel),rmatrix)
    
    ! Solution vectors -- for the nonlinearity (if there is one).
    ! For the previous (1), current (2) and next (3) time step.
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(1),.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(2),.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(3),.TRUE.)
    
    ! Load the last solution vectors for handling the nonlinearity.
    ! rtempVectorSol(2) holds the data for the current timestep.
    ! rtempVectorSol(1) holds the data from the previous and
    ! rtempVectorSol(3) that of the next timestep.
    IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
      CALL sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, &
          0, rtempVectorSol(2))
    END IF
    
    ! Loop through the substeps
    
    DO isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*p_rspaceTimeDiscr%rtimeDiscr%dtstep
      rproblem%rtimedependence%itimestep = isubstep

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL cc_updateDiscreteBC (rproblem)
      END IF
      
      ! The first and last substep is a little bit special concerning
      ! the matrix!
      IF (isubstep .EQ. 0) THEN
        
        ! We are in the first substep
      
        ! -----
        
        ! Get the evaluation point for the nonlinearity in the next timestep
        ! into rtempVectorSol(3).
        IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
          CALL sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, &
              1+isubstep+1, rtempVectorSol(3))
        END IF

        ! The diagonal matrix.
      
        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents)
          
        ! Assemble the matrix. No 'previous' solution vector.
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rmatrix,rmatrixComponents,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3)) 
          
        ! Assemble the system matrix on level p_rspaceTimeDiscr%ilevel.
        ! Include the boundary conditions into the matrices.
        !CALL cc_assembleLinearisedMatrices (&
        !    rnonlinearIterationTmp,rproblem%rcollection,&
        !    .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
            
        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rmatrix,rglobalA,1,1)
            
        ! -----
      
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the primal velocity.

        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,1,rmatrixComponents)
      
        ! Assemble the matrix
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rmatrix,rmatrixComponents,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3)) 
        
        CALL lsysbl_duplicateMatrix (&
            rmatrix,&
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

      ELSE IF (isubstep .LT. p_rspaceTimeDiscr%NEQtime-1) THEN
        
        ! We are sonewhere in the middle of the matrix. There is a substep
        ! isubstep+1 and a substep isubstep-1!
        
        ! -----
        
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the primal velocity.

        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,-1,rmatrixComponents)

        ! Assemble the matrix
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rmatrix,rmatrixComponents,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3)) 
        
        CALL lsysbl_duplicateMatrix (&
            rmatrix,rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
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

        ! -----      

        ! Now the diagonal matrix.

        ! Assemble the nonlinear defect.
        ! We use explicit Euler, so the weights are easy.
      
        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents)
            
        ! Assemble the matrix
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rmatrix,rmatrixComponents,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3)) 
        
        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rmatrix,&
            rglobalA,isubstep*6+1,isubstep*6+1)
            
        ! -----
        
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,1,rmatrixComponents)

        ! Assemble the matrix
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rmatrix,rmatrixComponents,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3)) 
        
        CALL lsysbl_duplicateMatrix (&
            rmatrix,rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
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

      ELSE
      
        ! We are in the last substep
        
        ! -----
        
        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,-1,rmatrixComponents)
      
        ! Assemble the matrix
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rmatrix,rmatrixComponents,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3)) 
        
        CALL lsysbl_duplicateMatrix (&
            rmatrix,rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
       
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

        ! -----
        
        ! The diagonal matrix.
      
        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents)
        
        ! Assemble the matrix
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rmatrix,rmatrixComponents,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3)) 
        
        ! Insert the system matrix for the dual equation to our global matrix.
        CALL insertMatrix (rmatrix,rglobalA,isubstep*6+1,isubstep*6+1)

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)
      
      END IF
    
      ! Shift the evaluation vectors: 1 <- 2 <- 3
      IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
        CALL lsysbl_copyVector (rtempVectorSol(2),rtempVectorSol(1))
        CALL lsysbl_copyVector (rtempVectorSol(3),rtempVectorSol(2))
      END IF
    
    END DO
    
    ! Release the temp matrix and vectors
    CALL lsysbl_releaseMatrix (rmatrix)
    
    ! Update structural information of the global matrix.
    CALL lsysbl_updateMatStrucInfo (rglobalA)
    CALL lsysbl_releaseVector (rtempVectorSol(3))
    CALL lsysbl_releaseVector (rtempVectorSol(2))
    CALL lsysbl_releaseVector (rtempVectorSol(1))

    ! Write the global matrix to a file.
    !CALL matio_writeBlockMatrixHR(rglobalA,'MATRIX',.TRUE.,0,'matrixcn.txt','(1X,E20.10)')
    !'(E13.2)')
    
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

!  ! ***************************************************************************
!  
!!<subroutine>
!
!  SUBROUTINE cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, &
!      dnorm, rb)
!
!!<description>
!  ! This routine assembles the space-time defect d=b-Ax. rd must have been
!  ! initialised by the space-time RHS vector b. The routine will then
!  ! calculate rd = rd - A rx to get the space-time defect.
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!
!  ! A structure defining the space time matrix of the corresponding system.
!  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rspaceTimeMatrix
!!</input>
!
!!<inputoutput>
!  ! A space-time vector defining the current solution.
!  ! Is replaced by the new solution
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx
!
!  ! A space-time vector with the space-time RHS b. Is overwritten by
!  ! d=b-Ax.
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!
!  ! OPTIONAL: A space-time vector with the space-time RHS b. 
!  ! If not specified, the routine assumes that rd contains the RHS.
!  TYPE(t_spacetimeVector), INTENT(INOUT), OPTIONAL :: rb
!!</inputoutput>
!
!!<output>
!  ! OPTIONAL: If specified, returns the $l_2$-norm if the defect.
!  REAL(DP), INTENT(OUT), OPTIONAL :: dnorm
!!<output>
!
!!</subroutine>
!
!    ! local variables
!    INTEGER :: isubstep,ilevel,icp,irelnonl
!    TYPE(t_vectorBlock) :: rtempVectorD
!    TYPE(t_vectorBlock), DIMENSION(3) :: rtempVector
!    TYPE(t_vectorBlock), DIMENSION(3) :: rtempVectorEval, rtempVectorEval2, rtempVectorEval3
!    TYPE(t_blockDiscretisation), POINTER :: p_rdiscr
!    REAL(DP) :: dtheta
!    TYPE(t_matrixBlock) :: rblockTemp
!    TYPE(t_ccmatrixComponents) :: rmatrixComponents
!    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
!    
!    ! DEBUG!!!
!    REAL(DP), DIMENSION(:), POINTER :: p_Dx1,p_Dx2,p_Dx3,p_Db
!    REAL(DP), DIMENSION(:), POINTER :: p_Deval1,p_Deval2,p_Deval3
!    
!    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscretisation
!    
!    ! Level of the discretisation
!    ilevel = p_rspaceTimeDiscr%ilevel
!
!    ! Theta-scheme identifier.
!    ! =1: implicit Euler.
!    ! =0.5: Crank Nicolson
!    dtheta = rproblem%rtimedependence%dtimeStepTheta
!    
!    ! Create a temp vector that contains the part of rd which is to be modified.
!    p_rdiscr => p_rspaceTimeDiscr%p_RlevelInfo%p_rdiscretisation
!    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorD,.FALSE.)
!    
!    ! The vector will be a defect vector. Assign the boundary conditions so
!    ! that we can implement them.
!    rtempVectorD%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
!    rtempVectorD%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC
!    
!    ! Create a temp vector for the X-vectors at timestep i-1, i and i+1.
!    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(1),.FALSE.)
!    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(2),.FALSE.)
!    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(3),.FALSE.)
!    
!    ! DEBUG!!!
!    CALL lsysbl_getbase_double (rtempVector(1),p_Dx1)
!    CALL lsysbl_getbase_double (rtempVector(2),p_Dx2)
!    CALL lsysbl_getbase_double (rtempVector(3),p_Dx3)
!    CALL lsysbl_getbase_double (rtempVectorD,p_Db)
!    
!    ! Get the parts of the X-vector which are to be modified at first --
!    ! subvector 1, 2 and 3.
!    CALL sptivec_getTimestepData(rx, 0, rtempVector(1))
!    IF (p_rspaceTimeDiscr%niterations .GT. 0) &
!      CALL sptivec_getTimestepData(rx, 1, rtempVector(2))
!    IF (p_rspaceTimeDiscr%niterations .GT. 1) THEN
!      CALL sptivec_getTimestepData(rx, 2, rtempVector(3))
!    ELSE
!      CALL lsysbl_copyVector (rtempVector(2),rtempVector(3))
!    END IF
!      
!    ! Create temp vectors for the evaluation point of the nonlinearity.
!    ! If ry is not specified, share the vector content with rx, thus
!    ! use ry=rx.
!    IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
!      CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(1),.FALSE.)
!      CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(2),.FALSE.)
!      CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(3),.FALSE.)
!      
!      ! Get the first three evaluation points
!      CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, 0, rtempVectorEval(1))
!      IF (p_rspaceTimeDiscr%niterations .GT. 0) &
!        CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, 1, rtempVectorEval(2))
!      IF (p_rspaceTimeDiscr%niterations .GT. 1) THEN
!        CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, 2, rtempVectorEval(3))
!      ELSE
!        CALL lsysbl_copyVector (rtempVectorEval(2),rtempVectorEval(3))
!      END IF
!      
!      ! DEBUG!!!
!      CALL lsysbl_getbase_double (rtempVectorEval(1),p_Deval1)
!      CALL lsysbl_getbase_double (rtempVectorEval(2),p_Deval2)
!      CALL lsysbl_getbase_double (rtempVectorEval(3),p_Deval3)
!    ELSE
!      CALL lsysbl_duplicateVector (rtempVector(1),rtempVectorEval(1),&
!          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      CALL lsysbl_duplicateVector (rtempVector(2),rtempVectorEval(2),&
!          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      CALL lsysbl_duplicateVector (rtempVector(3),rtempVectorEval(3),&
!          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    END IF
!    
!    ! If dnorm is specified, clear it.
!    IF (PRESENT(dnorm)) THEN
!      dnorm = 0.0_DP
!    END IF
!
!    ! Basic initialisation of rmatrixComponents with the pointers to the
!    ! matrices / discretisation structures on the current level.
!    !
!    ! The weights in the rmatrixComponents structure are later initialised
!    ! according to the actual situation when the matrix is to be used.
!    rmatrixComponents%p_rdiscretisation         => &
!        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation
!    rmatrixComponents%p_rmatrixStokes           => &
!        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixStokes          
!    rmatrixComponents%p_rmatrixB1             => &
!        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB1              
!    rmatrixComponents%p_rmatrixB2             => &
!        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB2              
!    rmatrixComponents%p_rmatrixMass           => &
!        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixMass            
!    rmatrixComponents%p_rmatrixIdentityPressure => &
!        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixIdentityPressure
!    rmatrixComponents%dnu = collct_getvalue_real (rproblem%rcollection,'NU')
!    rmatrixComponents%iupwind1 = collct_getvalue_int (rproblem%rcollection,'IUPWIND1')
!    rmatrixComponents%dupsam1 = collct_getvalue_real (rproblem%rcollection,'UPSAM1')
!    rmatrixComponents%iupwind2 = collct_getvalue_int (rproblem%rcollection,'IUPWIND2')
!    rmatrixComponents%dupsam2 = collct_getvalue_real (rproblem%rcollection,'UPSAM2')
!    
!    ! Loop through the substeps
!    
!    DO isubstep = 0,p_rspaceTimeDiscr%niterations
!    
!      ! Current point in time
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*p_rspaceTimeDiscr%dtstep
!      rproblem%rtimedependence%itimestep = isubstep
!
!      ! Get the part of rd which is to be modified.
!      IF (PRESENT(rb)) THEN
!        ! The RHS is in rb
!        CALL sptivec_getTimestepData(rb, isubstep, rtempVectorD)
!      ELSE
!        ! The RHS is in rd
!        CALL sptivec_getTimestepData(rd, isubstep, rtempVectorD)
!      END IF
!
!      ! -----
!      ! Discretise the boundary conditions at the new point in time -- 
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL cc_updateDiscreteBC (rproblem)
!      END IF
!      
!      ! The first and last substep is a little bit special concerning
!      ! the matrix!
!      IF (isubstep .EQ. 0) THEN
!        
!        ! We are in the first substep. Here, we have to handle the following
!        ! part of the supersystem:
!        !
!        !  ( A11 A12   0   0 ... )  ( x1 )  =  ( f1 )
!        !  ( ... ... ... ... ... )  ( .. )     ( .. )
!        !
!        ! So we have to compute:
!        !
!        !  d1  :=  f1  -  A11 x1  -  A12 x2
!        !
!        ! -----
!        
!        ! The diagonal matrix.
!      
!        ! Set up the matrix weights of that submatrix.
!        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,0,rmatrixComponents)
!          
!        ! Subtract: rd = rd - A11 x1
!        irelnonl = ABS(iidxNonlin)-isubstep + 1
!        CALL cc_assembleDefect (rmatrixComponents,rtempVector(1),rtempVectorD,&
!            1.0_DP,rtempVectorEval(irelnonl))
!
!        ! -----
!      
!        ! Create the matrix
!        !   A12  :=  -M/dt + dtheta*[-nu\Laplace u + u \grad u]
!        ! and subtract A12 x2 from rd.
!        !
!        ! Set up the matrix weights of that submatrix.
!        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,1,rmatrixComponents)
!
!        ! Subtract: rd = rd - A12 x2
!        irelnonl = ABS(iidxNonlin)-isubstep + 1
!        CALL cc_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
!            1.0_DP,rtempVectorEval(irelnonl))
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!        
!        ! The primal defect is =0 because of the initial condition.
!        ! This is important, otherwise the nonlinear loop will not converge
!        ! because of a lack in the defect in the 0th timestep!!!
!        CALL lsyssc_clearVector (rtempVectorD%RvectorBlock(1))
!        CALL lsyssc_clearVector (rtempVectorD%RvectorBlock(2))
!        CALL lsyssc_clearVector (rtempVectorD%RvectorBlock(3))
!
!      ELSE IF (isubstep .LT. p_rspaceTimeDiscr%niterations) THEN
!
!        ! We are sonewhere in the middle of the matrix. There is a substep
!        ! isubstep+1 and a substep isubstep-1!  Here, we have to handle the following
!        ! part of the supersystem:
!        !
!        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
!        !  ( ... Aii-1 Aii Aii+1 ... )  ( xi )  =  ( fi )
!        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
!        !
!        ! So we have to compute:
!        !
!        !  dn  :=  fn  -  Aii-1 xi-1  -  Aii xi  -  Aii+1 xi+1
!        
!        ! -----
!        
!        ! Create the matrix
!        !   Aii-1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
!        ! and include that into the global matrix for the primal velocity.
!
!        ! Set up the matrix weights of that submatrix.
!        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,-1,rmatrixComponents)
!            
!        ! Subtract: rd = rd - Aii-1 xi-1
!        ! Note that at this point, the nonlinearity must be evaluated
!        ! at xi due to the discretisation scheme!!!
!        irelnonl = ABS(iidxNonlin)-isubstep + 2
!        CALL cc_assembleDefect (rmatrixComponents,rtempVector(1),rtempVectorD,&
!            1.0_DP,rtempVectorEval(irelnonl))
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----      
!
!        ! Now the diagonal matrix.
!      
!        ! Assemble the nonlinear defect.
!      
!        ! Set up the matrix weights of that submatrix.
!        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,0,rmatrixComponents)
!
!        ! Subtract: rd = rd - Aii xi
!        irelnonl = ABS(iidxNonlin)-isubstep + 2
!        CALL cc_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
!            1.0_DP,rtempVectorEval(irelnonl))
!            
!        ! -----
!        
!        ! Create the matrix
!        !   Aii+1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
!        ! and include that into the global matrix for the dual velocity.
!
!        ! Set up the matrix weights of that submatrix.
!        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,1,rmatrixComponents)
!          
!        ! Subtract: rd = rd - Aii+1 xi+1
!        irelnonl = ABS(iidxNonlin)-isubstep + 2
!        CALL cc_assembleDefect (rmatrixComponents,rtempVector(3),rtempVectorD,&
!            1.0_DP,rtempVectorEval(irelnonl))
!        
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!      ELSE 
!        
!        ! We are in the last substep. Here, we have to handle the following
!        ! part of the supersystem:
!        !
!        !  ( ... ... ... ...   ... )  ( .. )     ( .. )
!        !  ( ... ...   0 Ann-1 Ann )  ( xn )  =  ( fn )
!        !
!        ! So we have to compute:
!        !
!        !  dn  :=  fn  -  Ann-1 xn-1  -  Ann xn
!        
!        ! -----
!        
!        ! Create the matrix
!        !   Ann-1 = -M/dt + dtheta*[-nu\Laplace u + u \grad u]
!        ! and include that into the global matrix for the dual velocity.
!
!        ! Set up the matrix weights of that submatrix.
!        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,-1,rmatrixComponents)
!          
!        ! Subtract: rd = rd - Ann-1 xn-1
!        ! Note that at this point, the nonlinearity must be evaluated
!        ! at xn due to the discretisation scheme!!!
!        irelnonl = ABS(iidxNonlin)-isubstep + 3
!        CALL cc_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
!            1.0_DP,rtempVectorEval(irelnonl))
!     
!        ! -----
!        
!        ! Now the diagonal matrix.
!      
!        ! Assemble the nonlinear defect.
!      
!        ! Set up the matrix weights of that submatrix.
!        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,0,rmatrixComponents)
!
!        ! Subtract: rd = rd - Ann xn
!        irelnonl = ABS(iidxNonlin)-isubstep + 3
!        CALL cc_assembleDefect (rmatrixComponents,rtempVector(3),rtempVectorD,&
!            1.0_DP,rtempVectorEval(irelnonl))
!      
!      END IF
!      
!      ! Implement the boundary conditions into the defect.
!      CALL vecfil_discreteBCdef (rtempVectorD)
!      CALL vecfil_discreteFBCdef (rtempVectorD)      
!      
!      ! Save the defect vector back to rd.
!      CALL sptivec_setTimestepData(rd, isubstep, rtempVectorD)
!      
!      ! If dnorm is specified, calculate the norm of the sub-defect vector and
!      ! add it to dnorm.
!      IF (PRESENT(dnorm)) THEN
!        IF (rproblem%MT_outputLevel .GE. 2) THEN
!          CALL output_line ('||D_'//TRIM(sys_siL(isubstep,10))//'|| = '//&
!              TRIM(sys_sdEL(&
!                  SQRT(lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)&
!                  ),10)) )
!          DO icp=1,6
!            CALL output_line ('  ||D_'//TRIM(sys_siL(isubstep,10))//'^'//TRIM(sys_siL(icp,2))&
!                //'|| = '//&
!                TRIM(sys_sdEL(&
!                    SQRT(lsyssc_vectorNorm(rtempVectorD%RvectorBlock(icp),LINALG_NORML2)&
!                    ),10)) )
!          END DO
!        END IF
!        dnorm = dnorm + lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)**2
!      END IF
!      
!      IF ((isubstep .GT. 0) .AND. &
!          (isubstep .LT. p_rspaceTimeDiscr%niterations-1)) THEN
!      
!        ! Shift the timestep data: x_n+1 -> x_n -> x_n-1
!        CALL lsysbl_copyVector (rtempVector(2), rtempVector(1))
!        CALL lsysbl_copyVector (rtempVector(3), rtempVector(2))
!        
!        ! Get the new x_n+1 for the next pass through the loop.
!        CALL sptivec_getTimestepData(rx, isubstep+2, rtempVector(3))
!        
!        ! The same for the 'evaluation vector' ry -- if it's specified.
!        ! If not specified, rtempVectorX and rtempVectorEvalX share the
!        ! same memory, so we don't have to do anything.
!        IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
!          CALL lsysbl_copyVector (rtempVectorEval(2), rtempVectorEval(1))
!          CALL lsysbl_copyVector (rtempVectorEval(3), rtempVectorEval(2))
!          CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, isubstep+2, rtempVectorEval(3))
!        END IF
!        
!      ELSE IF ((p_rspaceTimeDiscr%niterations .EQ. 1) .AND. (isubstep .EQ. 0)) THEN
!      
!        ! Only one timestep. Put vector-1 -> vector-2 -> vector-3 so the above
!        ! assembly routine for the last timestep will work.
!        CALL lsysbl_copyVector (rtempVector(2),rtempVector(3))
!        CALL lsysbl_copyVector (rtempVectorEval(2),rtempVectorEval(3))
!
!        CALL lsysbl_copyVector (rtempVector(1),rtempVector(2))
!        CALL lsysbl_copyVector (rtempVectorEval(1),rtempVectorEval(2))
!        
!      END IF
!    
!    END DO
!    
!    ! If dnorm is specified, normalise it.
!    ! It was calculated from p_rspaceTimeDiscr%niterations+1 subvectors.
!    IF (PRESENT(dnorm)) THEN
!      dnorm = SQRT(dnorm) / REAL(p_rspaceTimeDiscr%niterations+1,DP)
!    END IF
!    
!    ! Release the temp vectors.
!    CALL lsysbl_releaseVector (rtempVectorEval(3))
!    CALL lsysbl_releaseVector (rtempVectorEval(2))
!    CALL lsysbl_releaseVector (rtempVectorEval(1))
!    CALL lsysbl_releaseVector (rtempVector(3))
!    CALL lsysbl_releaseVector (rtempVector(2))
!    CALL lsysbl_releaseVector (rtempVector(1))
!    CALL lsysbl_releaseVector (rtempVectorD)
!    
!  END SUBROUTINE 
!   
!  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE cc_solveSupersystemDirect (rproblem, rspaceTimeDiscr, rx, rb, rd)
  
!<description>
  ! This routine assembles and solves the time-space coupled supersystem:
  ! $Ax=b$. The RHS vector is generated on-the-fly.
  ! The routine generates the full matrix in memory and solves with UMFPACK, 
  ! so it should only be used for debugging!
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN),TARGET :: rspaceTimeDiscr
!</input>

!<inputoutput>
  ! A space-time vector defining the initial solution. Is replaced by a new
  ! solution vector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx

  ! A temporary space-time vector that receives the RHS during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb

  ! A temporary space-time vector that receives the defect during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! The nonlinear solver configuration
    INTEGER :: isubstep,iglobIter
    LOGICAL :: bneumann
    CHARACTER(LEN=SYS_STRLEN) :: sstring,slinearSolver
    INTEGER :: nminIterations,nmaxIterations

    REAL(DP) :: ddefNorm,dinitDefNorm,depsRel,depsAbs
    
    TYPE(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx

    ! A temporary vector in the size of a spatial vector.
    TYPE(t_vectorBlock) :: rtempVectorX

    ! A second temporary vector in the size of a spatial vector.
    TYPE(t_vectorBlock) :: rtempVectorB

    ! A third temporary vector for the nonlinear iteration
    TYPE(t_vectorBlock) :: rtempVector

    ! Create temp vectors for X, B and D.
    CALL lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVector,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorB,.TRUE.)
        
    ! Attach the boundary conditions to the temp vectors.
    rtempVector%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVector%p_rdiscreteBCfict => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorX%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorX%p_rdiscreteBCfict => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorB%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorB%p_rdiscreteBCfict => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

!    ! Implement the bondary conditions into all initial solution vectors
!    DO isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals
!    
!      ! Current point in time
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%rtimeDiscr%dtstep
!      rproblem%rtimedependence%itimestep = isubstep
!
!      ! -----
!      ! Discretise the boundary conditions at the new point in time -- 
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL cc_updateDiscreteBC (rproblem)
!      END IF
!      
!      ! Implement the boundary conditions into the global solution vector.
!      CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
!      
!      ! DEBUG!!!
!      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
!      
!      CALL cc_implementBC (rproblem,rvector=rtempVectorX)
!      
!      CALL sptivec_setTimestepData(rx, isubstep, rtempVectorX)
!      
!    END DO

    CALL tbc_implementBCsolution (rproblem,rspaceTimeDiscr,rx)
    
    ! Generate the RHS for the initial condition.
    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rinitialCondRHS,.false.)
    call sptivec_getTimestepData (rx, 1, rtempVectorX)
    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
    call cc_generateInitCondRHS (rproblem,rspaceTimeDiscr,&
        rtempVectorX,rinitialCondRHS)

    ddefNorm = 1.0_DP

    ! ---------------------------------------------------------------
    ! Set up the structure of the global space time matrix.
    ! As we can only handle linear subproblems here, we don't have
    ! to initialise the evaluation point of the nonlinearity.

    rspaceTimeMatrix%p_rspaceTimeDiscretisation => rspaceTimeDiscr
    rspaceTimeMatrix%cmatrixType = 0
    
    ! ---------------------------------------------------------------
    ! Solve the global space-time coupled system.
    !
    ! Get the initial defect: d=b-Ax
    !CALL cc_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
    !  rtempvectorX, rtempvectorB, rtempvector, .FALSE.)    
    CALL trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rb, .FALSE.)

    ! Implement the initial condition into the RHS.
    CALL tbc_implementInitCond (rproblem, rb, rinitialCondRHS, rtempvector)    
    CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)    

    ! Now work with rd, our 'defect' vector
    CALL sptivec_copyVector (rb,rd)

    ! Assemble the defect.
    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
    CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,rproblem%MT_outputLevel .GE. 2)
        
    dinitDefNorm = ddefNorm
    
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    CALL output_separator (OU_SEP_EQUAL)        

    ! Call the routine to generate the global matrix and to solve the system.
    CALL cc_applyUmfpackToSupersystem (rproblem, rspaceTimeMatrix, rx, rd, &
      rtempvectorX, rtempvectorB, rtempVector)

    ! Calculate the final defect    
    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rd, &
    !  rtempvectorX, rtempvectorB, rtempvector,.FALSE.)
    CALL trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rd, .FALSE.)

    ! Implement the initial condition into the RHS.
    CALL tbc_implementInitCond (rproblem, rd, rinitialCondRHS, rtempvector)    
    CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)    

    ! Assemble the defect
    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
    CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,rproblem%MT_outputLevel .GE. 2)
        
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    CALL output_separator (OU_SEP_EQUAL)        
    
    ! Do we have Neumann boundary?
    bneumann = rspaceTimeDiscr%p_rlevelInfo%bhasNeumannBoundary
    
    IF (.NOT. bneumann) THEN
      ! Normalise the primal and dual pressure to integral mean value zero.
      CALL tbc_pressureToL20 (rx,rtempVectorX)
    END IF
    
    CALL lsysbl_releaseVector (rtempVectorB)
    CALL lsysbl_releaseVector (rtempVectorX)
    CALL lsysbl_releaseVector (rtempVector)
    call lsysbl_releaseVector (rinitialCondRHS)
    call lsysbl_releaseVector (rinitialCondSol)
          
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE cc_precondDefectSupersystem (rproblem, rspaceTimeMatrix, rx, rb, rd, &
      rtempvectorX,  rtempvectorD, rpreconditioner)

!<description>
  ! This routine performs preconditioning with the nonlinear super-defect
  ! vector rd: $d = A^{-1} d$.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! space time matrix structure defining the matrix A.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(INout) :: rspaceTimeMatrix
  
  ! A space-time vector defining the current solution.
  TYPE(t_spacetimeVector), INTENT(IN) :: rx

  ! A space-time vector defining the current RHS.
  TYPE(t_spacetimeVector), INTENT(IN) :: rb

!</input>

!<inputoutput>
  ! A spatial preconditioner. This one is applied to each substep in the
  ! global matrix.
  TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner
  
  ! A temporary vector in the size of a spatial vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX

  ! A third temporary vector for the nonlinear iteration
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD
  
  ! A space-time vector that receives the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: isubstep,ilevel,itejac
    REAL(DP) :: dtheta,dtstep,ddefNorm,dinitDef
    LOGICAL :: bsuccess
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    TYPE(t_vectorBlock) :: rtempVectorX1,rtempVectorX3
    TYPE(t_spacetimeVector) :: rcurrentx,rcurrentd
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dd
    
    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscretisation
    
    rspaceTimeMatrix%cmatrixType=1
    
    dtheta = rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep

    ! Level of the discretisation
    ilevel = p_rspaceTimeDiscr%ilevel
    
    ! The weights in the rmatrixComponents structure are later initialised
    ! according to the actual situation when the matrix is to be used.
    rmatrixComponents%p_rdiscretisation         => &
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation
    rmatrixComponents%p_rmatrixStokes           => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixStokes          
    rmatrixComponents%p_rmatrixB1             => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB1              
    rmatrixComponents%p_rmatrixB2             => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB2              
    rmatrixComponents%p_rmatrixMass           => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixMass            
    rmatrixComponents%p_rmatrixIdentityPressure => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixIdentityPressure
    rmatrixComponents%dnu = collct_getvalue_real (rproblem%rcollection,'NU')
    rmatrixComponents%iupwind1 = collct_getvalue_int (rproblem%rcollection,'IUPWIND1')
    rmatrixComponents%dupsam1 = collct_getvalue_real (rproblem%rcollection,'UPSAM1')
    rmatrixComponents%iupwind2 = collct_getvalue_int (rproblem%rcollection,'IUPWIND2')
    rmatrixComponents%dupsam2 = collct_getvalue_real (rproblem%rcollection,'UPSAM2')
    
    ! Create two additional temp vectors
    CALL lsysbl_createVecBlockIndirect (rtempVectorX,rtempVectorX1,.TRUE.)
    CALL lsysbl_createVecBlockIndirect (rtempVectorX,rtempVectorX3,.TRUE.)
    
    ! Create a temp vector for the current iterates of the preconditioner
    call sptivec_initVectorDiscr (rcurrentx,rx%p_rtimeDiscretisation,rx%p_rblockDiscretisation)
    call sptivec_initVectorDiscr (rcurrentd,rd%p_rtimeDiscretisation,rd%p_rblockDiscretisation)

    ! Create the initial defect.
    call sptivec_copyVector (rd,rcurrentd)
    CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rcurrentx, rcurrentd, &
      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, dinitDef,rproblem%MT_outputLevel .GE. 2)
    ddefNorm = dinitDef
    call output_line("Block-JAC: Ite=0, ||res||="//ADJUSTL(sys_sdEP(ddefNorm,20,10)))
    
    ! ----------------------------------------------------------------------
    ! We use a block-Jacobi scheme for preconditioning...
    !
    ! Loop until convergence.
    do itejac = 1,100
    
      ! Sopping criterion: 2 digits.
      if (ddefNorm .lE. dinitDef*1.0E-2_DP) exit
    
      ! For this purpose, loop through the substeps.
      
      DO isubstep = 1,p_rspaceTimeDiscr%NEQtime
      
        ! Current time step?
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + (isubstep-1) * dtstep
        rproblem%rtimedependence%itimestep = isubstep-1

        CALL output_line ('Block-Jacobi preconditioning of timestep: '//&
            TRIM(sys_siL(isubstep,10))//&
            ' Time: '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,10)))
      
        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF

        ! DEBUG!!!      
        CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
        CALL lsysbl_getbase_double (rtempVectorD,p_Dd)

        ! Read in the RHS/solution/defect vector of the current timestep.
        CALL sptivec_getTimestepData (rx, isubstep, rtempVectorX)
        CALL sptivec_getTimestepData (rcurrentd, isubstep, rtempVectorD)

        IF (isubstep .GT. 1) THEN
          CALL sptivec_getTimestepData (rx, isubstep-1, rtempVectorX1)
        END IF

        IF (isubstep .LE. p_rspaceTimeDiscr%NEQtime-1) THEN
          CALL sptivec_getTimestepData (rx, isubstep+1, rtempVectorX3)
        END IF
        
        ! Set up the matrix weights for the diagonal matrix
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents)
          
        ! Perform preconditioning of the defect with the method provided by the
        ! core equation module.
        CALL cc_precondDefect (rpreconditioner,rmatrixComponents,rtempVectorD,&
          rtempVectorX1,rtempVectorX,rtempVectorX3,&
          bsuccess,rproblem%rcollection)      
      
        ! Save back the preconditioned defect.
        CALL sptivec_setTimestepData (rcurrentd, isubstep, rtempVectorD)
        
      END DO
      
      ! Add the correction to the current iterate
      call sptivec_vectorLinearComb(rcurrentd,rcurrentx,1.0_DP,1.0_DP)

      ! Create the new defect
      call sptivec_copyVector (rd,rcurrentd)
      CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rcurrentx, rcurrentd, &
        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,.false.)
      call output_line("Block-JAC: Ite="//trim(sys_siL(itejac,10))//&
          ", ||res||="//ADJUSTL(sys_sdEP(ddefNorm,20,10)))
      
    end do ! itejac
    
    ! Return the preconditioned defect
    call sptivec_copyVector (rcurrentx,rd)
    
    ! Release temp vectors
    CALL lsysbl_releaseVector (rtempVectorX1)
    CALL lsysbl_releaseVector (rtempVectorX3)
    call sptivec_releaseVector (rcurrentx)
    call sptivec_releaseVector (rcurrentd)
    
    rspaceTimeMatrix%cmatrixType=0
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE cc_solveSupersystemDefCorr (rproblem, rspaceTimeDiscr, rx, rb, rd,&
      ctypePreconditioner)
  
!<description>
  ! This is a primitive space-time defect correction solver. It applies the
  ! iteration
  !    $$ rx = rx + C^{-1} ( rb - A rx ) $$
  ! to a space time vector rx and a space time RHS vector rb.
  ! Here, $C^{-1}$ is a spatial preconditioner (linear solver) that is applied 
  ! to each time step during the iteration. The configuration of this defect
  ! correction loop is specified in the '[TIME-DEFCORR]' section in the DAT files.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN),TARGET :: rspaceTimeDiscr
  
  ! Type of preconditioner to use for the space time system.
  ! =1: Standard linear system.
  ! =2: Newton iteration
  INTEGER, INTENT(IN) :: ctypePreconditioner
!</input>

!<inputoutput>
  ! A space-time vector defining the initial solution. Is replaced by a new
  ! solution vector.
  TYPE(t_spacetimeVector), INTENT(INOUT), TARGET :: rx

  ! A temporary space-time vector that receives the RHS during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb

  ! A temporary space-time vector that receives the defect during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! The nonlinear solver configuration
    INTEGER :: isubstep,iglobIter
    LOGICAL :: bneumann
    CHARACTER(LEN=SYS_STRLEN) :: sstring,slinearSolver
    INTEGER :: nminIterations,nmaxIterations
    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol

    REAL(DP) :: ddefNorm,dinitDefNorm,depsRel,depsAbs,dlastDefNorm,depsDiff
    
    TYPE(t_ccspatialPreconditioner) :: rpreconditioner
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx

    ! A temporary vector in the size of a spatial vector.
    TYPE(t_vectorBlock) :: rtempVectorX

    ! A second temporary vector in the size of a spatial vector.
    TYPE(t_vectorBlock) :: rtempVectorB

    ! A third temporary vector for the nonlinear iteration
    TYPE(t_vectorBlock) :: rtempVector
    
    ! A structure for the global space-time matrix on the highest level.
    TYPE(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix,rspaceTimePreconditioner

    ! Create temp vectors for X, B and D.
    CALL lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVector,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorB,.TRUE.)
        
    ! Attach the boundary conditions to the temp vectors.
    rtempVector%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVector%p_rdiscreteBCfict => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorX%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorX%p_rdiscreteBCfict => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorB%p_rdiscreteBC => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorB%p_rdiscreteBCfict => rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    ! Some preparations for the spatial preconditioner.
    !
    ! Get the name of the section containing the linear solver from the DAT file.
    CALL parlst_getvalue_string (rproblem%rparamList, 'TIME-DEFCORR', &
                                 'slinearSolver', sstring, '')
    READ(sstring,*) slinearSolver

    ! Initialise the preconditioner for the preconditioning in every timestep.
    ! Specify slinearSolver as the name of the section that configures the
    ! spatial preconditioner. This is (up to now only) the name of a section
    ! containing configuration data of a linear solver.
    CALL cc_initPreconditioner (rproblem,&
        rproblem%NLMIN,rproblem%NLMAX,rpreconditioner)
    CALL cc_configPreconditioner (rproblem,rpreconditioner,slinearSolver,&
        ctypePreconditioner)

!    ! Implement the bondary conditions into all initial solution vectors
!    DO isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals
!    
!      ! Current point in time
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%rtimeDiscr%dtstep
!      rproblem%rtimedependence%itimestep = isubstep
!
!      ! -----
!      ! Discretise the boundary conditions at the new point in time -- 
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL cc_updateDiscreteBC (rproblem)
!      END IF
!      
!      ! Implement the boundary conditions into the global solution vector.
!      CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
!      
!      ! DEBUG!!!
!      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
!      
!      CALL cc_implementBC (rproblem,rvector=rtempVectorX)
!      
!      CALL sptivec_setTimestepData(rx, isubstep, rtempVectorX)
!      
!    END DO
    CALL tbc_implementBCsolution (rproblem,rspaceTimeDiscr,rx)
    
    ! Generate the RHS for the initial condition.
    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rinitialCondRHS,.false.)
    call sptivec_getTimestepData (rx, 1, rtempVectorX)
    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
    call cc_generateInitCondRHS (rproblem,rspaceTimeDiscr,&
        rtempVectorX,rinitialCondRHS)

    ddefNorm = 1.0_DP
    
    ! ---------------------------------------------------------------
    ! Set up the structure of the global space time matrix.
    ! Set the evaluation point of the matrix to the current solution
    ! vector.

    rspaceTimeMatrix%p_rspaceTimeDiscretisation => rspaceTimeDiscr
    rspaceTimeMatrix%cmatrixType = 0
    rspaceTimeMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
    rspaceTimeMatrix%p_rsolution => rx
    
    ! Set up a structure for the matrix that serves as space-time

    ! preconditioner. This is based on the space time matrix...
    
    rspaceTimePreconditioner = rspaceTimeMatrix
    IF ((rproblem%iequation .EQ. 0) .AND. (ctypePreconditioner .EQ. 1)) THEN
      ! ...but may also be the Newton matrix!
      rspaceTimePreconditioner%cmatrixType = 1
    END IF

    ! ---------------------------------------------------------------
    ! Solve the global space-time coupled system.
    !
    ! Get the initial defect: d=b-Ax
    !CALL cc_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
    !  rtempvectorX, rtempvectorB, rtempvector, .FALSE.)    
    CALL trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rb, .FALSE.)

    ! Implement the initial condition into the RHS.
    CALL tbc_implementInitCond (rproblem, rb, rinitialCondRHS, rtempvector)    
    CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)    

    ! Now work with rd, our 'defect' vector
    CALL sptivec_copyVector (rb,rd)

    ! Assemble the defect.
    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
    CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,rproblem%MT_outputLevel .GE. 2)
        
    dinitDefNorm = ddefNorm
    dlastDefNorm = 0.0_DP
    
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    CALL output_separator (OU_SEP_EQUAL)        

    ! Get some solver parameters for the iteration
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-DEFCORR', &
                              'nminIterations', nminIterations, 0)

    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-DEFCORR', &
                              'nmaxIterations', nmaxIterations, 0)

    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-DEFCORR', &
                                 'depsRel', depsRel, 1.0E-5_DP)

    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-DEFCORR', &
                                 'depsDiff', depsDiff, 1.0E-5_DP)

    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-DEFCORR', &
                                 'depsAbs', depsAbs, 1.0E99_DP)

    iglobIter = 0
    
!    !CALL cc_solveSupersysDirect (rproblem, rspaceTimeDiscr, rx, rd, &
!    !  rtempvectorX, rtempvectorB, rtempVector)
!    CALL cc_solveSupersysDirectCN (rproblem, rspaceTimeDiscr, rx, rd, &
!      rtempvectorX, rtempvectorB, rtempVector)
    
    DO WHILE ((iglobIter .LT. nminIterations) .OR. &
              ((((ddefNorm .GT. depsRel*dinitDefNorm) .OR. (ddefNorm .GE. depsAbs)) .AND.&
                (ABS(ddefNorm-dlastDefNorm) .GE. depsDiff*dlastDefNorm)) .AND. &
               (iglobIter .LT. nmaxIterations)))
    
      iglobIter = iglobIter+1
      
      ! Preconditioning of the defect: d=C^{-1}d
      CALL cc_precondDefectSupersystem (rproblem, rspaceTimePreconditioner, &
          rx, rb, rd, rtempvectorX,  rtempVector, rpreconditioner)

      ! Add the defect: x = x + omega*d          
      CALL sptivec_vectorLinearComb (rd,rx,1.0_DP,1.0_DP)
          
      ! Assemble the new defect: d=b-Ax
      CALL sptivec_copyVector (rb,rd)
      !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
      CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,rproblem%MT_outputLevel .GE. 2)
          
      CALL output_separator (OU_SEP_EQUAL)
      CALL output_line ('Iteration: '//trim(sys_siL(iglobIter,10))//&
          '              Defect of supersystem: '//adjustl(sys_sdEP(ddefNorm,20,10)))
      CALL output_separator (OU_SEP_EQUAL)
      
    END DO
    
    ! ---------------------------------------------------------------
    ! Release the preconditioner of the nonlinear iteration
    CALL cc_donePreconditioner (rpreconditioner)
    
    !CALL cc_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rd, &
    !  rtempvectorX, rtempvectorB, rtempvector,.FALSE.)
    CALL trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rd, .FALSE.)

    ! Implement the initial condition into the RHS.
    CALL tbc_implementInitCond (rproblem, rd, rinitialCondRHS, rtempvector)    
    CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)    

    ! Assemble the defect
    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
    CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,rproblem%MT_outputLevel .GE. 2)
        
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    CALL output_separator (OU_SEP_EQUAL)        
    
    ! Do we have Neumann boundary?
    bneumann = rspaceTimeDiscr%p_rlevelInfo%bhasNeumannBoundary
    
    IF (.NOT. bneumann) THEN
      ! Normalise the primal and dual pressure to integral mean value zero.
      CALL tbc_pressureToL20 (rx,rtempVectorX)
    END IF
    
    CALL lsysbl_releaseVector (rtempVectorB)
    CALL lsysbl_releaseVector (rtempVectorX)
    CALL lsysbl_releaseVector (rtempVector)
    call lsysbl_releaseVector (rinitialCondRHS)
    call lsysbl_releaseVector (rinitialCondSol)
          
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE cc_solveSupersystemMultigrid (rproblem, RspaceTimeDiscr, rx, rb, rd,&
      ctypePreconditioner)
  
!<description>
  ! This subroutine solves the nonstationary space time coupled (Navier-)Stokes
  ! optimal control problem. For this purpose, a nonlinear defect correction
  ! loop is executed. Every defect is preconditioned by a space-time coupled
  ! preconditioner like Block-Jacobi, block SOR or whatever.
  !
  ! The caller must provide problem information in rproblem and a set of
  ! matrix configurations for all time levels which are allowed for the solver
  ! to use. If multiple time-levels are provided, space-time coupled multigrid
  ! is used for preconditioning. matrix configurations must be provided in
  ! RspaceTimeDiscr. The maximum level in this array defines the level where
  ! the system is solved. This level must correspond to the vectors rx, rb and 
  ! rd.
  !
  ! The caller can provide an initial solution in rx. However, rb is
  ! overwritten by a newly generated space-time coupled RHS vector.
  ! rd is used as temporary vector.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! An array of t_ccoptSpaceTimeDiscretisation structure defining all the
  ! levels of the coupled space-time the discretisation.
  ! The solution/rhs rx/rb must correspond to the maximum level in this
  ! array.
  TYPE(t_ccoptSpaceTimeDiscretisation), &
      DIMENSION(:), INTENT(IN), TARGET :: RspaceTimeDiscr

  ! Type of preconditioner to use for the space time system.
  ! =1: Standard linear system.
  ! =2: Newton iteration
  INTEGER, INTENT(IN) :: ctypePreconditioner
!</input>

!<inputoutput>
  ! A space-time vector defining the initial solution. Is replaced by a new
  ! solution vector.
  TYPE(t_spacetimeVector), INTENT(INOUT), TARGET :: rx

  ! A temporary space-time vector that receives the RHS during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT), TARGET :: rb

  ! A temporary space-time vector that receives the defect during the calculation.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!</inputoutput>

!</subroutine>

    ! The nonlinear solver configuration
    INTEGER :: isubstep,iglobIter,ierror,ilev,ispacelev,nsmSteps,cspaceTimeSmoother
    INTEGER :: ilowerSpaceLevel,itemp,ctypeCoarseGridSolver,i,nlmax
    LOGICAL :: bneumann
    REAL(DP), DIMENSION(4) :: Derror
    INTEGER(I32) :: nminIterations,nmaxIterations
    REAL(DP) :: depsRel,depsAbs,domega,domegaPrecond,depsDiff
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    TYPE(t_ccspatialPreconditioner), DIMENSION(SIZE(RspaceTimeDiscr)) :: RspatialPrecond
    TYPE(t_ccspatialPreconditioner), DIMENSION(SIZE(RspaceTimeDiscr)) :: RspatialPrecondPrimal
    TYPE(t_ccspatialPreconditioner), DIMENSION(SIZE(RspaceTimeDiscr)) :: RspatialPrecondDual
    TYPE(t_sptiProjection), DIMENSION(SIZE(RspaceTimeDiscr)) :: RinterlevelProjection
    TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(SIZE(RspaceTimeDiscr)) :: RspaceTimePrecondMatrix
    TYPE(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
    TYPE(t_vectorBlock) :: rtempVecCoarse,rtempVecFine
    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol
    
    ! A solver node that identifies our solver.
    TYPE(t_sptilsNode), POINTER :: p_rprecond,p_rsolverNode
    TYPE(t_sptilsNode), POINTER :: p_rmgSolver,p_rsmoother,p_rcgrSolver
    TYPE(t_sptilsNode), POINTER :: p_rpresmoother,p_rpostsmoother
    
    ! TYPE(t_spaceTimeVector) :: rtemp

    REAL(DP) :: ddefNorm,dinitDefNorm,dlastDefNorm,dtempdef,depsrelLinSol
    REAL(DP) :: dinexactNewtonExponent,dinexactNewtonEpsRel
    
    CHARACTER(LEN=SYS_STRLEN) :: slinearSolver,sstring
    
    INTEGER :: inpos, innextrec
    CHARACTER :: ctest
    
    TYPE(t_timer) :: rtimerMGStep,rtimerNonlinear,rtimerPreconditioner
    INTEGER :: ilinearIterations
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx

    ! A temporary vector in the size of a spatial vector.
    TYPE(t_vectorBlock) :: rtempVectorX

    ! A second temporary vector in the size of a spatial vector.
    TYPE(t_vectorBlock) :: rtempVectorB

    ! A third temporary vector for the nonlinear iteration
    TYPE(t_vectorBlock) :: rtempVector

    ! STATISTICS: Total time needed for smoothing operations
    TYPE(t_timer) :: rtimeSmoothing
        
    ! STATISTICS: Total time needed for the coarse grid solver
    TYPE(t_timer) :: rtimeCoarseGridSolver
    
    ! STATISTICS: Time needed for linear algebra stuff (matrix-vector, 
    ! vector-copy, prolongation/restriction,...)
    TYPE(t_timer) :: rtimeLinearAlgebra
    
    ! STATISTICS: Time needed for prolongation/restriction
    TYPE(t_timer) :: rtimeProlRest
    
    ! STATISTICS: Time for solving problems in space.
    TYPE(t_timer) :: rtimeSpacePrecond

    ! STATISTICS: Time for initialisation / factorisation of the space time system;
    ! global and in one step
    TYPE(t_timer) :: rtimeFactorisation,rtimeFactorisationStep

    ! Get a poiter to the discretisation structure on the maximum
    ! space/time level. That's the level where we solve here.
    p_rspaceTimeDiscr => RspaceTimeDiscr(SIZE(RspaceTimeDiscr))
    
    ! Create temp vectors for X, B and D.
    CALL lsysbl_createVecBlockByDiscr (&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVectorX,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVectorB,.TRUE.)
        
    ! Attach the boundary conditions to the temp vectors.
    rtempVector%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVector%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorX%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorX%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorB%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorB%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    ! Implement the bondary conditions into all initial solution vectors
    CALL tbc_implementBCsolution (rproblem,p_rspaceTimeDiscr,rx,rtempvectorX)

    ! Generate the RHS for the initial condition.
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rinitialCondRHS,.false.)
    call sptivec_getTimestepData (rx, 1, rtempVectorX)
    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
    call cc_generateInitCondRHS (rproblem,p_rspaceTimeDiscr,&
        rtempVectorX,rinitialCondRHS)

    ! We set up a space-time preconditioner, e.g. in the following configuration:
    ! Main Preconditioner: Multigrid
    !     -> Presmoother:  Block Jacobi/GS
    !     -> Postsmoother: Block Jacobi/GS
    !     -> CGr-Solver:   Defect correction (or UMFPACK)
    !        -> Preconditioner: Block Jacobi
    !
    ! So we start creating the main solver: Multigrid.
    !
    ! Create as many time levels as specified by the length of RspatialPrecond.
    CALL sptils_initMultigrid (rproblem,1,SIZE(RspatialPrecond),p_rmgSolver)
    
    ! Loop over the time levels.
    DO ilev=1,SIZE(RspatialPrecond)
    
      ! Type of smoother to use?
      CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
          'cspaceTimeSmoother', cspaceTimeSmoother, 0)
          
      ! Type of coarse grid solver?
      CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                'ctypeCoarseGridSolver', ctypeCoarseGridSolver, 0)

      ! Get the refinement level in space that belongs to this space-time level.    
      ispacelev = RspaceTimeDiscr(ilev)%ilevel
      
      IF (ilev .EQ. 1) THEN
        i = ctypeCoarseGridSolver
        
        ! Get the name of the section containing the linear solver from the DAT file --
        ! in case we use a linear solver as spatial preconditioner.
        CALL parlst_getvalue_string (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                    'slinearSolver', sstring, '')
        READ(sstring,*) slinearSolver

      ELSE
        i = cspaceTimeSmoother
        
        ! Get the name of the section containing the linear solver from the DAT file --
        ! in case we use a linear solver as spatial preconditioner.
        CALL parlst_getvalue_string (rproblem%rparamList, 'TIME-SMOOTHER', &
                                    'slinearSolver', sstring, '')
        READ(sstring,*) slinearSolver

      END IF
      
      SELECT CASE (i)
      CASE (0:)
        ! Initialise the spatial preconditioner for Block Jacobi
        ! (note: this is slightly expensive in terms of memory!
        ! Probably we could use the same preconditioner for all levels,
        ! but this has still to be implemeted and is a little bit harder!)
        CALL cc_initPreconditioner (rproblem,rproblem%nlmin,&
            ispacelev,RspatialPrecond(ilev))
        ! Specify slinearSolver as the name of the section that configures the
        ! spatial preconditioner. This is (up to now only) the name of a section
        ! containing configuration data of a linear solver.
        CALL cc_configPreconditioner (rproblem,RspatialPrecond(ilev),&
            slinearSolver,ctypePreconditioner)
         
      CASE DEFAULT
        PRINT *,'Unknown preconditioner/smoother: ',i
        CALL sys_halt()
      END SELECT

      ! Generate an interlevel projection structure for that level.
      ! Note that space restriction/prolongation must be switched off if
      ! we are on the spatial coarse mesh!
      IF (ilev .GT. 1) THEN
        ilowerSpaceLevel = RspaceTimeDiscr(ilev-1)%ilevel
      ELSE
        ilowerSpaceLevel = ispacelev
      END IF
      CALL sptipr_initProjection (rinterlevelProjection(ilev),&
          RspaceTimeDiscr(ilev)%p_rlevelInfo%rdiscretisation,&
          ilowerSpaceLevel .NE. ispacelev)
         
      IF (ilev .EQ. 1) THEN
        ! ..., on the minimum level, create a coarse grid solver, ...

        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                    'domega', domega, 1.0_DP)
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
                                    'domega', domegaPrecond, 1.0_DP)
        
        SELECT CASE (ctypeCoarseGridSolver)
        CASE (0)
          ! Block Jacobi preconditioner
          CALL sptils_initBlockJacobi (rproblem,p_rprecond,domega,RspatialPrecond(ilev))

          ! Defect correction solver
          CALL sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
          
        CASE (1)
          ! Block SOR preconditioner
          CALL sptils_initBlockSOR (rproblem,p_rprecond,domegaPrecond,RspatialPrecond(ilev))

          ! Defect correction solver
          CALL sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)

        CASE (2)
          ! Forward backward Gauss Seidel
          CALL sptils_initBlockFBGS (rproblem,p_rprecond,&
            domega,domegaPrecond,RspatialPrecond(ilev))

          ! Defect correction solver
          CALL sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
          
        CASE (3)
          ! CG with Block Jacobi as preconditioner
          CALL sptils_initBlockJacobi (rproblem,p_rprecond,domegaPrecond,RspatialPrecond(ilev))

          ! CG solver        
          CALL sptils_initCG (rproblem,p_rcgrSolver,p_rprecond)
          
        CASE (4)
          ! UMFACK Gauss elimination
          CALL sptils_initUMFPACK4 (rproblem,p_rcgrSolver)
          p_rcgrSolver%domega = domega

        CASE (5)
          ! Forward backward Gauss Seidel as preconditioner
          CALL sptils_initBlockFBGS (rproblem,p_rprecond,&
            domegaPrecond,domegaPrecond,RspatialPrecond(ilev))
          !CALL sptils_initBlockJacobi (rproblem,p_rprecond,domegaPrecond,RspatialPrecond(ilev))

          ! BiCGStab solver        
          CALL sptils_initBiCGStab (rproblem,p_rcgrSolver,p_rprecond)

        CASE (6)
          ! UMFPACK
          CALL sptils_initUMFPACK4 (rproblem,p_rprecond)

          ! Defect correction solver
          CALL sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
          p_rcgrSolver%domega = domega
          
        CASE DEFAULT
          PRINT *,'Unknown solver: ',ctypeCoarseGridSolver
          STOP
        END SELECT
        
        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
            'nminIterations', p_rcgrSolver%nminIterations, 1)
        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
            'nmaxIterations', p_rcgrSolver%nmaxIterations, 100)
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                    'depsRel', p_rcgrSolver%depsRel, 1E-5_DP)
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                    'depsAbs', p_rcgrSolver%depsAbs, 1E-5_DP)
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                    'depsDiff', p_rcgrSolver%depsDiff, 0.0_DP)
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                    'ddivRel', p_rcgrSolver%ddivRel, 1.0_DP)
        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
                                'istoppingCriterion', p_rcgrSolver%istoppingCriterion, 0)

        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
            'ioutputLevel', p_rcgrSolver%ioutputLevel, 100)
            
        ! If there's a subsolver, configure it.
        IF (ASSOCIATED(p_rprecond)) THEN
          CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
              'nminIterations', p_rprecond%nminIterations, 1)
          CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
              'nmaxIterations', p_rprecond%nmaxIterations, 100)
          CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
                                      'depsRel', p_rprecond%depsRel, 1E-5_DP)
          CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
                                      'depsAbs', p_rprecond%depsAbs, 1E-5_DP)
          CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
                                      'depsDiff', p_rprecond%depsDiff, 0.0_DP)
          CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
                                      'ddivRel', p_rprecond%ddivRel, 1.0_DP)
          CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
                                  'istoppingCriterion', p_rprecond%istoppingCriterion, 0)

          CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
              'ioutputLevel', p_rprecond%ioutputLevel, 100)
        END IF
        
        ! ...; finally initialise the level with that       
        CALL sptils_setMultigridLevel (p_rmgSolver,ilev,&
                      rinterlevelProjection(ilev),&
                      NULL(),NULL(),p_rcgrSolver)
        
      ELSE
        ! ... on higher levels, create a smoother, ...
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
                                    'domega', domega, 1.0_DP)
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
                                    'domegaPrecond', domegaPrecond, 1.0_DP)
        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
            'nsmSteps', nsmSteps, 1)

        ! DEBUG!!!
!        IF (ilev .LE. SIZE(RspatialPrecond)-1) THEN
!          cspaceTimeSmoother = 7
!          ! nsmSteps = nsmSteps/2
!        END IF
        ! DEBUG!!!

        SELECT CASE (cspaceTimeSmoother)
        CASE (0)
          ! Block Jacobi
          CALL sptils_initBlockJacobi (rproblem,p_rsmoother,&
            domega,RspatialPrecond(ilev))
        CASE (1)
          ! Block SOR
          CALL sptils_initBlockSOR (rproblem,p_rsmoother,&
            domegaPrecond,RspatialPrecond(ilev))
        CASE (2)
          ! Block Forward-Backward Gauss-Seidel
          CALL sptils_initBlockFBGS (rproblem,p_rsmoother,&
            domega,domegaPrecond,RspatialPrecond(ilev))
        CASE (3)
          ! CG with Block Jacobi as preconditioner
          CALL sptils_initBlockJacobi (rproblem,p_rprecond,&
            domegaPrecond,RspatialPrecond(ilev))
          CALL sptils_initCG (rproblem,p_rsmoother,p_rprecond)
        CASE (4)
          ! CG with Block Jacobi as preconditioner
          CALL sptils_initBlockFBGS (rproblem,p_rprecond,domegaPrecond,&
            domegaPrecond,RspatialPrecond(ilev))
          CALL sptils_initCG (rproblem,p_rsmoother,p_rprecond)
        CASE (6)
          ! DefCorr with UMFPACK as preconditioner
          CALL sptils_initUMFPACK4 (rproblem,p_rprecond)
          CALL sptils_initDefCorr (rproblem,p_rsmoother,p_rprecond)
        CASE (7)
          ! BiCGStab with FBGS as preconditioner
          CALL sptils_initBlockFBGS (rproblem,p_rprecond,&
            domega,domegaPrecond,RspatialPrecond(ilev))
            p_rprecond%nminIterations = 1
            p_rprecond%nmaxIterations = 1
          CALL sptils_initBiCGStab (rproblem,p_rsmoother,p_rprecond)
        CASE DEFAULT
          PRINT *,'Unknown smoother: ',cspaceTimeSmoother
          STOP
        END SELECT
        
        CALL sptils_convertToSmoother (p_rsmoother,nsmSteps,domega)
        
        ! Switch off smoothing is set that way in the DAT file
        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
                                  'ioutputLevel', p_rsmoother%ioutputLevel, 1)

        p_rpresmoother => p_rsmoother
        p_rpostsmoother => p_rsmoother
        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
                                  'ipresmoothing', itemp, 1)
        IF (itemp .EQ. 0) NULLIFY(p_rpresmoother)
        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
                                  'ipostsmoothing', itemp, 1)
        IF (itemp .EQ. 0) NULLIFY(p_rpostsmoother)
        
        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
                                    'depsRel', p_rsmoother%depsRel, 0.0_DP)

        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
                                    'depsAbs', p_rsmoother%depsAbs, 0.0_DP)

        CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
                                    'depsDiff', p_rsmoother%depsDiff, 0.0_DP)

        CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
            'istoppingCriterion', p_rsmoother%istoppingCriterion, 1)

        ! DEBUG!!!
!        IF (ilev .EQ. SIZE(RspatialPrecond)-1) THEN
!          p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels(&
!            SIZE(p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels)-1)%depsRelCycle = 1E-10_DP
!          p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels(&
!            SIZE(p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels)-1)%depsAbsCycle = 1E-14_DP
!        END IF
        ! DEBUG!!!

        IF ((.NOT. ASSOCIATED(p_rpresmoother)) .AND. &
            (.NOT. ASSOCIATED(p_rpostsmoother))) THEN
          CALL sptils_releaseSolver(p_rsmoother)
        END IF

        ! ...; finally initialise the level with that       
        CALL sptils_setMultigridLevel (p_rmgSolver,ilev,&
                      rinterlevelProjection(ilev),&
                      p_rpresmoother,p_rpostsmoother,NULL())

      END IF
       
    END DO
    
    ! Our main solver is MG now.
    p_rsolverNode => p_rmgSolver
    
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
                             'nmaxIterations', p_rsolverNode%nmaxIterations, 1)
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
                             'nminIterations', p_rsolverNode%nminIterations, 10)
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
                             'ioutputLevel', p_rsolverNode%ioutputLevel, 0)
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
                             'icycle', p_rmgSolver%p_rsubnodeMultigrid%icycle, 0)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
                                'depsRel', p_rmgSolver%depsRel, 1E-5_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
                                'depsAbs', p_rmgSolver%depsAbs, 1E-5_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
                                'depsDiff', p_rmgSolver%depsDiff, 0.0_DP)
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
                             'istoppingCriterion', p_rmgSolver%istoppingCriterion, 0)

    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
                                'dalphaMin', p_rmgSolver%p_rsubnodeMultigrid%dalphaMin,&
                                1.0_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
                                'dalphaMax', p_rmgSolver%p_rsubnodeMultigrid%dalphaMax,&
                                1.0_DP)
                                
    ! Save the relative stopping criterion of the linear solver; we
    ! probably need it for the adaptive Newton.
    depsrelLinSol = p_rmgSolver%depsRel
                                
    ! For the optimal coarse grid correction, we multiply the residuum of the
    ! 3rd and 6th equation by -1; gives better results (matrix symmetry?!?).
    ! Note that there is currently no INIT/DONE-Routine for the weights, so we
    ! do that manually...
    ALLOCATE(p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(6))
    p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(:) = 1.0_DP
    p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(3) = -1.0_DP
    p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(6) = -1.0_DP
    
    ! Initialise the basic parameters of the system matrices on all levels.
    DO ilev=1,SIZE(RspatialPrecond)
    
      ! Pointer to the corresponding space time discretisation structure
      RspaceTimePrecondMatrix(ilev)%p_rspaceTimeDiscretisation => RspaceTimeDiscr(ilev)
      
      ! Configure the matrix type; standard or Newton matrix
      SELECT CASE (ctypePreconditioner)
      CASE (CCPREC_LINEARSOLVER)
        ! Standard system matrix
        RspaceTimePrecondMatrix(ilev)%cmatrixType = 0
        
      CASE (CCPREC_NEWTON,CCPREC_INEXACTNEWTON)
        ! Newton matrix
        RspaceTimePrecondMatrix(ilev)%cmatrixType = 1
        RspaceTimePrecondMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
      END SELECT
      
    END DO 
    
    ! Allocate space-time vectors on all lower levels that hold the solution vectors
    ! for the evaluation of the nonlinearity (if we have a nonlinearity).
    DO ilev=1,SIZE(RspatialPrecond)-1
      ALLOCATE(RspaceTimePrecondMatrix(ilev)%p_rsolution)
      CALL sptivec_initVector (RspaceTimePrecondMatrix(ilev)%p_rsolution,&
        RspaceTimeDiscr(ilev)%NEQtime,&
        RspaceTimeDiscr(ilev)%p_rlevelInfo%rdiscretisation)
    END DO

    ! On the maximum level attach the solution vector.
    nlmax = UBOUND(RspaceTimeDiscr,1)
    RspaceTimePrecondMatrix(nlmax)%p_rsolution => rx
    
    ! Create a space time matrix for the maximum level which serves as basis for
    ! setting up the defect. This is the actual system matrix, thus it doesn't
    ! contain any Newton parts!
    rspaceTimeMatrix%p_rspaceTimeDiscretisation => RspaceTimeDiscr(ilev)
    rspaceTimeMatrix%p_rsolution => rx
    
    ! That matrix also has to apply projection operators when being applied
    ! to a vector -- in case control constraints are active.
    rspaceTimeMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
    
    ! Attach matrix information to the linear solver
    CALL sptils_setMatrices (p_rsolverNode,RspaceTimePrecondMatrix)
    
    ! Initialise the space-time preconditioner
    CALL stat_clearTimer (rtimeFactorisation)
    CALL stat_startTimer (rtimeFactorisation)
    CALL sptils_initStructure (p_rsolverNode,ierror)
    CALL stat_stopTimer (rtimeFactorisation)
    
    ddefNorm = 1.0_DP
    
    ! ---------------------------------------------------------------
    ! Solve the global space-time coupled system.
    !
    ! Get the initial defect: d=b-Ax
    !CALL cc_assembleSpaceTimeRHS (rproblem, &
    !  RspaceTimePrecondMatrix(nlmax)%p_rspaceTimeDiscretisation, rb, &
    !  rtempvectorX, rtempvectorB, rtempvector, .FALSE.)    
    CALL trhsevl_assembleRHS (rproblem, &
      RspaceTimePrecondMatrix(nlmax)%p_rspaceTimeDiscretisation, rb, .FALSE.)

    ! Implement the initial condition into the RHS.
    CALL tbc_implementInitCond (rproblem, rb, rinitialCondRHS, rtempvector)
    CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)
        
    ! DEBUG!!!
    !CALL sptivec_saveToFileSequence (rb,&
    !    '(''./debugdata/initrhs.txt.'',I5.5)',.TRUE.)

    ! Now work with rd, our 'defect' vector
    CALL sptivec_copyVector (rb,rd)

    ! Assemble the defect.
    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, &
    !    rx, rd, ddefNorm)
    CALL sptivec_copyVector (rb,rd)
    CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,&
      ddefNorm,rproblem%MT_outputLevel .GE. 2)

    ! DEBUG!!!
    !CALL sptivec_saveToFileSequence (rb,&
    !    '(''./debugdata/initdef.txt.'',I5.5)',.TRUE.)
        
    dinitDefNorm = ddefNorm
    dlastDefNorm = 0.0_DP
    if (dinitDefNorm .eq. 0.0_DP) then
      ! Trick to avoid div/0.
      dinitDefNorm = 1.0_DP
    end if
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    ! Value of the functional
    CALL cc_optc_nonstatFunctional (rproblem,&
        RspaceTimePrecondMatrix(SIZE(RspatialPrecond))%p_rsolution,&
        rtempVector,RspaceTimeDiscr(SIZE(RspatialPrecond))%dalphaC,&
        RspaceTimeDiscr(SIZE(RspatialPrecond))%dgammaC,&
        Derror)
    CALL output_line ('||y-z||       = '//TRIM(sys_sdEL(Derror(1),10)))
    CALL output_line ('||u||         = '//TRIM(sys_sdEL(Derror(2),10)))
    CALL output_line ('||y(T)-z(T)|| = '//TRIM(sys_sdEL(Derror(3),10)))
    CALL output_line ('J(y,u)        = '//TRIM(sys_sdEL(Derror(4),10)))
    CALL output_separator (OU_SEP_EQUAL)        

    ! Get configuration parameters from the DAT file
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SOLVER', &
                              'nminIterations', nminIterations, 1)
    CALL parlst_getvalue_int (rproblem%rparamList, 'TIME-SOLVER', &
                             'nmaxIterations', nmaxIterations, 10)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
                                 'depsRel', depsRel, 1E-5_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
                                 'depsAbs', depsAbs, 1E-5_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
                                 'depsDiff', depsDiff, 0.0_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
                                 'domega', domega, 1.0_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
                                 'domega', domega, 1.0_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
                                 'dinexactNewtonEpsRel', dinexactNewtonEpsRel, 1.0E-2_DP)
    CALL parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
                                 'dinexactNewtonExponent', dinexactNewtonExponent, 2.0_DP)

    ! Initalise statistic variables
    CALL stat_clearTimer (rtimeSmoothing)
    CALL stat_clearTimer (rtimeCoarseGridSolver)
    CALL stat_clearTimer (rtimeLinearAlgebra)
    CALL stat_clearTimer (rtimeProlRest)
    CALL stat_clearTimer (rtimeSpacePrecond)
    
    CALL stat_clearTimer (rtimerPreconditioner)
    CALL stat_clearTimer (rtimerNonlinear)
    CALL stat_startTimer (rtimerNonlinear)

    iglobIter = 0
    ilinearIterations = 0
    
    DO WHILE ((iglobIter .LT. nminIterations) .OR. &
              ((((ddefNorm .GT. depsRel*dinitDefNorm) .OR. (ddefNorm .GE. depsAbs)) .AND.&
                (ABS(ddefNorm-dlastDefNorm) .GE. depsDiff*dlastDefNorm)) &
              .AND. (iglobIter .LT. nmaxIterations)))
    
      iglobIter = iglobIter+1
      
      ! Project the solution down to all levels, so the nonlinearity
      ! on the lower levels can be calculated correctly.
      ! Use the memory in rtempvectorX and rtempvectorB as temp memory;
      ! it's large enough.
      CALL lsysbl_duplicateVector (rtempvectorX,rtempVecFine,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsysbl_duplicateVector (rtempvectorB,rtempVecCoarse,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      DO ilev=SIZE(RspatialPrecond)-1,1,-1
        CALL lsysbl_enforceStructureDiscr(&
            RspaceTimeDiscr(ilev+1)%p_rlevelInfo%rdiscretisation,rtempVecFine)
        rtempVecFine%p_rdiscreteBC =>&
            RspaceTimeDiscr(ilev+1)%p_rlevelInfo%p_rdiscreteBC
        rtempVecFine%p_rdiscreteBCfict =>&
            RspaceTimeDiscr(ilev+1)%p_rlevelInfo%p_rdiscreteFBC
            
        CALL lsysbl_enforceStructureDiscr(&
            RspaceTimeDiscr(ilev)%p_rlevelInfo%rdiscretisation,rtempVecCoarse)
        rtempVecCoarse%p_rdiscreteBC =>&
            RspaceTimeDiscr(ilev)%p_rlevelInfo%p_rdiscreteBC
        rtempVecCoarse%p_rdiscreteBCfict =>&
            RspaceTimeDiscr(ilev)%p_rlevelInfo%p_rdiscreteFBC

        ! Interpolate down
        CALL sptipr_performInterpolation (RinterlevelProjection(ilev+1),&
            RspaceTimePrecondMatrix(ilev)%p_rsolution, RspaceTimePrecondMatrix(ilev+1)%p_rsolution,&
            rtempVecCoarse,rtempVecFine)
            
        ! Set boundary conditions
        CALL tbc_implementBCsolution (rproblem,RspaceTimeDiscr(ilev),&
            RspaceTimePrecondMatrix(ilev)%p_rsolution,rtempVecCoarse)
      END DO
      
      CALL lsysbl_releaseVector(rtempVecFine)
      CALL lsysbl_releaseVector(rtempVecCoarse)
          
!      IF (rproblem%MT_outputLevel .GE. 2) THEN
!        CALL output_line ('Writing solution to file '//&
!            './ns/tmpsolution'//TRIM(sys_siL(iglobIter-1,10)))
!        CALL sptivec_saveToFileSequence(&
!            myRspaceTimeDiscr(SIZE(RspatialPrecond))%p_rsolution,&
!            '(''./ns/tmpsolution'//TRIM(sys_siL(iglobIter-1,10))//'.'',I5.5)',&
!            .TRUE.,rtempVectorX)
!          
!        DO ilev=1,SIZE(RspatialPrecond)
!          CALL cc_postprocSpaceTimeGMV(rproblem,myRspaceTimeDiscr(ilev),&
!              myRspaceTimeDiscr(ilev)%p_rsolution,&
!              './gmv/iteration'//TRIM(sys_siL(iglobIter-1,10))//'level'//TRIM(sys_siL(ilev,10))//&
!              '.gmv')
!        END DO
!      END IF
      
      IF (ctypePreconditioner .eq. CCPREC_INEXACTNEWTON) THEN
      
        ! Determine the stopping criterion for the inexact Newton.
        ! This is an adaptive stopping criterion depending on the current
        ! defect. In detail, we calculate:
        !
        !   |b-Ax_{i+1}|         ( |b-Ax_i| ) exp             ( |b-Ax_i| )
        !   ------------ = min { ( -------- )     , depsrel * ( -------- ) }
        !     |b-Ax_0|           ( |b-Ax_0| )                 ( |b-Ax_0| )
        !
        ! see e.g. [Michael Hinze, Habilitation, p. 51]
        !
        ! Switch off the relative stopping criterion in the linear solver:
        
        p_rsolverNode%depsRel = 0.0_DP
        
        ! Calculate the new absolute stopping criterion:
        
        dtempdef = ddefNorm / dinitDefNorm
        
        p_rsolverNode%depsAbs = MIN(dtempDef**dinexactNewtonExponent,&
                                    dinexactNewtonEpsRel*dtempdef) * dinitDefNorm
        
        ! For the coarse grid solver, we choose the same stopping criterion.
        ! But just for safetyness, the coarse grid solver should gain at least
        ! one digit!
        p_rcgrSolver%depsRel = 1.0E-1_DP
        p_rcgrSolver%depsAbs = p_rsolverNode%depsAbs
      
      END IF

      IF (rproblem%MT_outputLevel .GE. 1) THEN
        ! Value of the functional
        CALL cc_optc_nonstatFunctional (rproblem,&
            rspaceTimeMatrix%p_rsolution,&
            rtempVector,RspaceTimeDiscr(SIZE(RspatialPrecond))%dalphaC,&
            RspaceTimeDiscr(SIZE(RspatialPrecond))%dgammaC,&
            Derror)
        CALL output_line ('||y-z||       = '//TRIM(sys_sdEL(Derror(1),10)))
        CALL output_line ('||u||         = '//TRIM(sys_sdEL(Derror(2),10)))
        CALL output_line ('||y(T)-z(T)|| = '//TRIM(sys_sdEL(Derror(3),10)))
        CALL output_line ('J(y,u)        = '//TRIM(sys_sdEL(Derror(4),10)))
        IF (ctypePreconditioner .eq. CCPREC_INEXACTNEWTON) THEN
          CALL output_lbrk ()
          CALL output_line ('Inexact Newton: Stopping criterion = '//&
              TRIM(sys_sdEL(p_rsolverNode%depsAbs,10)))
        END IF
        CALL output_separator (OU_SEP_EQUAL)
      END IF
            
      CALL stat_clearTimer (rtimerMGStep)

      ! DEBUG!!!
      !CALL sptivec_copyVector (rd,rtemp)
      
      ! Preconditioning of the defect: d=C^{-1}d
      IF (ASSOCIATED(p_rsolverNode)) THEN

        CALL stat_clearTimer (rtimeFactorisationStep)
        CALL stat_startTimer (rtimeFactorisationStep)
        CALL sptils_initData (p_rsolverNode,ierror)
        CALL stat_stopTimer (rtimeFactorisationStep)
        
        CALL stat_clearTimer (rtimerMGStep)
        CALL stat_startTimer (rtimerMGStep)
        CALL sptils_precondDefect (p_rsolverNode,rd)
        CALL sptils_doneData (p_rsolverNode)
        
        CALL stat_stopTimer (rtimerMGStep)
        
        ! Sum up time data for statistics.
        CALL stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeSmoothing,&
            rtimeSmoothing)
        CALL stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeCoarseGridSolver,&
            rtimeCoarseGridSolver)
        CALL stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeLinearAlgebra,&
            rtimeLinearAlgebra)
        CALL stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeProlRest,&
            rtimeProlRest)
        CALL stat_addtimers (p_rsolverNode%rtimeSpacePrecond,rtimeSpacePrecond)
        
        CALL output_lbrk ()
        CALL output_line ("Time for smoothing          : "//&
            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeSmoothing%delapsedReal,10))
        CALL output_line ("Time for coarse grid solving: "//&
            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeCoarseGridSolver%delapsedReal,10))
        CALL output_line ("Time for linear algebra     : "//&
            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeLinearAlgebra%delapsedReal,10))
        CALL output_line ("Time for prol/rest          : "//&
            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeProlRest%delapsedReal,10))
        CALL output_lbrk ()
        CALL output_line ("Time for prec. in space     : "//&
            sys_sdL(p_rsolverNode%rtimeSpacePrecond%delapsedReal,10))
        
        ! Count the number of linear iterations and the time for
        ! preconditioning
        ilinearIterations = ilinearIterations + p_rsolverNode%iiterations
        CALL stat_addtimers (rtimerMGStep,rtimerPreconditioner)
        CALL stat_addtimers (rtimeFactorisationStep,rtimeFactorisation)
        
      END IF
      
      !CALL output_line('Linear defect:')
      !CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rd, rtemp, &
      !  -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,.TRUE.)
      !CALL output_separator (OU_SEP_MINUS)
      
      ! Filter the defect for boundary conditions in space and time.
      ! Normally this is done before the preconditioning -- but by doing it
      ! afterwards, the initial conditions can be seen more clearly!
      !CALL tbc_implementInitCondDefect (p_rspaceTimeDiscr,rd,rtempVector)
      !CALL tbc_implementBCdefect (rproblem,p_rspaceTimeDiscr,rd,rtempVector)
      
      ! Add the defect: x = x + omega*d          
      CALL sptivec_vectorLinearComb (rd,rx,domega,1.0_DP)
      
      ! Do we have Neumann boundary?
      bneumann = p_rspaceTimeDiscr%p_rlevelInfo%bhasNeumannBoundary
      
      IF (.NOT. bneumann) THEN
        ! Normalise the primal and dual pressure to integral mean value zero.
        CALL tbc_pressureToL20 (rx,rtempVectorX)
      END IF
      
      ! Implement the initial condition to the new solution vector
      CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvectorX)    
      
      ! Are bounds to the control active? If yes, restrict the control
      ! to the allowed range.
!      if (rproblem%roptcontrol%ccontrolContraints .eq. 1) then
!        call cc_projectControl (rproblem,rx)
!      end if

      ! Call a parser that parses the script file commandfile.txt.
      ! This allows in-program modification of the problem sructure.
      rproblem%rdataOneshot%iglobIter = iglobIter
      rproblem%rdataOneshot%ddefNorm = ddefNorm
      rproblem%rdataOneshot%p_rx => rx
      rproblem%rdataOneshot%p_rb => rb
      
      CALL scr_readScript ('commandfile.txt',0,rproblem)
      
      iglobIter = rproblem%rdataOneshot%iglobIter
      ddefNorm = rproblem%rdataOneshot%ddefNorm
      
      ! Remember the last defect norm for the stopping criterion
      dlastDefNorm = ddefNorm
      
      ! Assemble the new defect: d=b-Ax
      CALL sptivec_copyVector (rb,rd)
      !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, &
      !    rx, rd, ddefNorm)
      IF (rproblem%MT_outputLevel .GE. 2) CALL output_line('Nonlinear defect:')
      
      CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,&
        ddefNorm,rproblem%MT_outputLevel .GE. 2)
          
      ! Filter the defect for boundary conditions in space and time.
      !CALL tbc_implementInitCondDefect (p_rspaceTimeDiscr,rd,rtempVector)
      !CALL tbc_implementBCdefect (rproblem,p_rspaceTimeDiscr,rd,rtempVector)
      
      CALL output_separator (OU_SEP_EQUAL)
      CALL output_line ('Iteration: '//TRIM(sys_siL(iglobIter,10))//&
          '                          Defect of supersystem:  '//ADJUSTL(sys_sdEP(ddefNorm,20,10)))
      CALL output_line ('Time for computation of this iterate: '//&
          TRIM(sys_sdL(rtimerMGStep%delapsedReal+rtimeFactorisationStep%delapsedReal,10)))
      CALL output_separator (OU_SEP_EQUAL)

    END DO
    
    CALL stat_stopTimer (rtimerNonlinear)
    
    ! Decrease iglobIter if the DO-loop was completely processed;
    ! iglobIter = nmaxIterations+1 in that case!
    IF (iglobIter .GT. nmaxIterations) iglobIter = nmaxIterations
    
    !CALL cc_assembleSpaceTimeRHS (rproblem, p_rspaceTimeDiscr, rd, &
    !  rtempvectorX, rtempvectorB, rtempvector,.FALSE.)
    CALL trhsevl_assembleRHS (rproblem, p_rspaceTimeDiscr, rd, .TRUE.)

    ! Implement the initial condition into the RHS/Solution.
    CALL tbc_implementInitCond (rproblem, rd, rinitialCondRHS, rtempvector)    
    CALL tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)

    ! Assemble the defect
    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, &  
    !    rx, rd, ddefNorm)
    CALL cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,rproblem%MT_outputLevel .GE. 2)
        
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
    ! Value of the functional
    CALL cc_optc_nonstatFunctional (rproblem,&
        rspaceTimeMatrix%p_rsolution,&
        rtempVector,RspaceTimeDiscr(SIZE(RspatialPrecond))%dalphaC,&
        RspaceTimeDiscr(SIZE(RspatialPrecond))%dgammaC,&
        Derror)
    CALL output_line ('||y-z||       = '//TRIM(sys_sdEL(Derror(1),10)))
    CALL output_line ('||u||         = '//TRIM(sys_sdEL(Derror(2),10)))
    CALL output_line ('||y(T)-z(T)|| = '//TRIM(sys_sdEL(Derror(3),10)))
    CALL output_line ('J(y,u)        = '//TRIM(sys_sdEL(Derror(4),10)))
    CALL output_separator (OU_SEP_EQUAL)
    CALL output_line ('Total computation time             = '// &
        TRIM(sys_sdL(rtimerNonlinear%delapsedReal,10)))
    CALL output_line ('Total time for factorisation       = '// &
        TRIM(sys_sdL(rtimeFactorisation%delapsedReal,10)))
    CALL output_line ('Total time for preconditioning     = '// &
        TRIM(sys_sdL(rtimerPreconditioner%delapsedReal,10)))
    CALL output_line ('#nonlinear iterations              = '//&
        TRIM(sys_siL(iglobIter,10)))
    CALL output_line ('#iterations preconditioner         = '//&
        TRIM(sys_siL(ilinearIterations,10)))

    CALL output_separator (OU_SEP_EQUAL)

    CALL output_line ("Preconditioner statistics:")
    CALL output_line ("Total time for smoothing           = "//&
        sys_sdL(rtimeSmoothing%delapsedReal,10))
    CALL output_line ("Total time for coarse grid solving = "//&
        sys_sdL(rtimeCoarseGridSolver%delapsedReal,10))
    CALL output_line ("Total time for linear algebra      = "//&
        sys_sdL(rtimeLinearAlgebra%delapsedReal,10))
    CALL output_line ("Total time for prol/rest           = "//&
        sys_sdL(rtimeProlRest%delapsedReal,10))
    CALL output_lbrk ()
    CALL output_line ("Total time for prec. in space      = "//&
        sys_sdL(rtimeSpacePrecond%delapsedReal,10))

    CALL output_separator (OU_SEP_EQUAL)

    ! Do we have Neumann boundary?
    bneumann = p_rspaceTimeDiscr%p_rlevelInfo%bhasNeumannBoundary
    
    IF (.NOT. bneumann) THEN
      ! Normalise the primal and dual pressure to integral mean value zero.
      CALL tbc_pressureToL20 (rx,rtempVectorX)
    END IF
    
    ! Release the multiplication weights for the energy minimisation.
    DEALLOCATE(p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights)
    
    ! Release the space-time and spatial preconditioner. 
    ! We don't need them anymore.
    CALL sptils_releaseSolver (p_rsolverNode)
    
    ! Release temp memory
    ! CALL sptivec_releaseVector (rtemp)
    
    ! Release the spatial preconditioner and temp vector on every level
    DO ilev=1,SIZE(RspatialPrecond)
      CALL cc_donePreconditioner (RspatialPrecondPrimal(ilev))
      CALL cc_donePreconditioner (RspatialPrecondDual(ilev))
      CALL cc_donePreconditioner (RspatialPrecond(ilev))
    END DO

    DO ilev=1,SIZE(RspatialPrecond)-1
      CALL sptivec_releaseVector (RspaceTimePrecondMatrix(ilev)%p_rsolution)
      DEALLOCATE(RspaceTimePrecondMatrix(ilev)%p_rsolution)
    END DO
          
    CALL lsysbl_releaseVector (rtempVectorB)
    CALL lsysbl_releaseVector (rtempVectorX)
    CALL lsysbl_releaseVector (rtempVector)
    call lsysbl_releaseVector (rinitialCondRHS)
    call lsysbl_releaseVector (rinitialCondSol)
    
  END SUBROUTINE
  
END MODULE
