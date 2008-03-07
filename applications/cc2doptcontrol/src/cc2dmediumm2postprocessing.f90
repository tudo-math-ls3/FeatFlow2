!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains postprocessing routines for the cc2dminim2_method2
!# CC2D solver.
!#
!# The following routines can be found here:
!#
!# 1.) cc_postprocessing
!#     -> Evaluate the solution of the stationary solver, write GMV-file.
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2postprocessing

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
  USE pprocerror
  USE pprocgradients
  
  USE collection
  USE convection
  
  USE ucd
  
  USE pprocnavierstokes
  
  USE cc2dmediumm2basic
  USE cc2dmediumm2optcanalysis
  USE cc2dmedium_callback
  USE cc2dmediumm2boundarydef
  
  USE spacetimediscretisation
  
  IMPLICIT NONE
  
!<types>

!<typeblock>

  ! Postprocessing structure. Whenever the postprocessing calculation routines are
  ! called, they fills this structure using the current solution vector (and other
  ! information if necessary). The information in this structure can then be used
  ! for GMV output e.g.
  TYPE t_c2d2postprocessing

    ! A discretisation structure that describes a piecewise constant discretisation
    ! (usually P0 or Q0).
    TYPE(t_spatialDiscretisation) :: rdiscrConstant
    
    ! A discretisation structure that describes a piecewise linear discretisation
    ! (usually P1 or Q1).
    TYPE(t_spatialDiscretisation) :: rdiscrLinear

    ! A discretisation structure that describes a piecewise quadratic discretisation
    ! (usually P2 or Q2).
    TYPE(t_spatialDiscretisation) :: rdiscrQuadratic

    ! A vector that describes the X-velocity field in the vertices
    TYPE(t_vectorScalar) :: rvectorVelX

    ! A vector that describes the Y-velocity field in the vertices
    TYPE(t_vectorScalar) :: rvectorVelY

    ! A vector that describes the pressure field in the vertices
    TYPE(t_vectorScalar) :: rvectorPressure

    ! A vector that describes the pressure field in the cells
    TYPE(t_vectorScalar) :: rvectorPressureCells

    ! A vector that describes the streamfunction
    TYPE(t_vectorScalar) :: rvectorStreamfunction
    
    ! A vector that describes the H1-error of the velocity field in the vertices
    TYPE(t_vectorScalar) :: rvectorH1err

    ! A vector that describes the H1-error of the pressure in the cells
    TYPE(t_vectorScalar) :: rvectorH1errCells
  
  END TYPE

!</typeblock>

!</types>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_postprocessingStationary (rproblem,rvector)
  
!<description>
  ! Postprocessing of solutions of stationary simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
!</input>

!</subroutine>

  ! local variables
    INTEGER :: ieltype
  
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! A vector accepting Q1 data
    TYPE(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    TYPE(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    TYPE(t_discreteBC), TARGET :: rdiscreteBC
    TYPE(t_discreteFBC), TARGET :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    
    ! Forces on the object
    REAL(DP), DIMENSION(NDIM2D) :: Dforces
    REAL(DP) :: df1,df2
    TYPE(t_boundaryRegion) :: rregion
    
    ! Divergence
    TYPE(t_matrixScalar) :: rBmatrix
    TYPE(t_vectorScalar), TARGET :: rtempVector

    ! If we have a uniform discreisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    IF ((rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
         ccomplexity .EQ. SPDISC_UNIFORM) .AND. &
        (boundary_igetNBoundComp(rproblem%p_rboundary) .GE. 2)) THEN

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      CALL boundary_createRegion (rproblem%p_rboundary, &
          2, 0, rregion)
      rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      df1 = 1.0_DP/1000.0_DP
      df2 = 0.1_DP * 0.2_DP**2
      CALL ppns2D_bdforces_uniform (rvector,rregion,Dforces,CUB_G1_1D,df1,df2)

      CALL output_lbrk()
      CALL output_line ('Body forces real bd., bdc/horiz/vert')
      CALL output_line (' 2 / ' &
          //TRIM(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //TRIM(sys_sdEP(Dforces(2),15,6)) )
      
    ENDIF
    
    ! Print out the value of the optimal control functional.
    CALL output_lbrk ()
    CALL cc_printControlFunctionalStat (rproblem,rvector)
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    IF (rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
        ccomplexity .EQ. SPDISC_UNIFORM) THEN
        
      ieltype = rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
                RelementDistribution(1)%itrialElement
                
      IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q1T) THEN
      
        ! Create a temporary vector 
        CALL lsyssc_createVecByDiscr (rvector%RvectorBlock(3)%p_rspatialDiscretisation,&
            rtempVector,.TRUE.)

        ! Calculate divergence = B1^T u1 + B2^T u2
        CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB1,&
            rBmatrix,LSYSSC_TR_VIRTUAL)
        CALL lsyssc_scalarMatVec (&
            rBmatrix, rvector%RvectorBlock(1), &
            rtempVector, 1.0_DP, 0.0_DP)
        CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB2,&
            rBmatrix,LSYSSC_TR_VIRTUAL)
        CALL lsyssc_scalarMatVec (&
            rBmatrix, rvector%RvectorBlock(2), &
            rtempVector, 1.0_DP, 1.0_DP)
        
        CALL output_lbrk()
        CALL output_line ('Divergence = ' &
            //TRIM(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)) )

        CALL lsyssc_releaseVector (rtempVector)
      
      END IF
      
    END IF    
    
    ! The solution vector is probably not in the way, GMV likes it!
    ! GMV for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! GMV understands.
    ! GMV understands only Q1 solutions! So the task is now to create
    ! a Q1 solution from rvector and write that out.
    !
    ! For this purpose, first create a 'derived' simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    CALL spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscretisation,rprjDiscretisation)
    
    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(1))

    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(2))

    ! Also use Q1 for the dual velocity field.
    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(4), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(4))

    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(5), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(5))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    CALL lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.FALSE.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    CALL spdp_projectSolution (rvector,rprjVector)
    
    ! Discretise the boundary conditions for this discretisation
    CALL bcasm_initDiscreteBC(rdiscreteBC)
    CALL cc_assembleBDconditions (rproblem,rprjDiscretisation,&
      rdiscreteBC,rproblem%rcollection)

    CALL bcasm_initDiscreteFBC(rdiscreteFBC)
    CALL cc_assembleFBDconditions (rproblem,rprjDiscretisation,&
      rdiscreteFBC,rproblem%rcollection)
    
    ! Connect the vector with the BC's
    rprjVector%p_rdiscreteBC => rdiscreteBC
    rprjVector%p_rdiscreteBCfict => rdiscreteFBC
    
    ! Filter the solution vector to implement discrete BC's.
    CALL vecfil_discreteBCsol (rprjVector)

    ! Filter the solution vector to implement discrete BC's for fictitious 
    ! boundary components.
    CALL vecfil_discreteFBCsol (rprjVector)
    
    ! Clean up the collection (as we are done with the assembly.
    CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    !
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,'gmv/u.gmv')
    
    ! Write the configuration of the application as comment block
    ! to the output file.
    CALL ucd_addCommentLine (rexport,'Configuration:')
    CALL ucd_addCommentLine (rexport,'---------------')
    CALL ucd_addParameterList (rexport,rproblem%rparamList)
    CALL ucd_addCommentLine (rexport,'---------------')

    ! Write velocity field
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)

    ! Write pressure
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)

    ! Dual velocity field
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
    CALL ucd_addVariableVertexBased (rexport,'X-vel-dual',UCD_VAR_STANDARD, p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'Y-vel-dual',UCD_VAR_STANDARD, p_Ddata2)
    
    ! Dual pressure
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
    CALL ucd_addVariableElementBased (rexport,'pressure-dual',UCD_VAR_STANDARD, p_Ddata)
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    IF (rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
        ccomplexity .EQ. SPDISC_UNIFORM) THEN
        
      ieltype = rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
                RelementDistribution(1)%itrialElement
                
      IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q1T) THEN
          
        CALL ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
        
        CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        CALL ucd_addVariableVertexBased (rexport,'streamfunction',&
            UCD_VAR_STANDARD, p_Ddata)

      END IF
      
    END IF
    
    ! Calculate the derivative of the pressure.
    ! First, project the pressure from Q0 to Q1 (if it's Q0).
    IF (rvector%p_rblockDiscretisation%RspatialDiscretisation(3)% &
        ccomplexity .EQ. SPDISC_UNIFORM) THEN
        
      ieltype = rvector%p_rblockDiscretisation%RspatialDiscretisation(3)% &
                RelementDistribution(1)%itrialElement
                
      IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q0) THEN
        ! rprjVector%RvectorBlock(4) is already prepared to accept Q1 solutions
        CALL spdp_projectSolutionScalar (rvector%RvectorBlock(3),rprjVector%RvectorBlock(4))

        ! Calculate the derivative of the pressure into the first two
        ! subvectors -- overwriting the projected velocity field, which we don't need
        ! anymore.
        CALL ppgrd_calcGradient (rprjVector%RvectorBlock(4),rprjVector)
        
        ! Write pressure derivative field
        CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
        
        CALL ucd_addVariableVertexBased (rexport,'X-deriv-pres',UCD_VAR_STANDARD, p_Ddata)
        CALL ucd_addVariableVertexBased (rexport,'Y-deriv-pres',UCD_VAR_STANDARD, p_Ddata2)
            
      END IF
          
      ieltype = rvector%p_rblockDiscretisation%RspatialDiscretisation(6)% &
                RelementDistribution(1)%itrialElement
                
      IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q0) THEN
        ! rprjVector%RvectorBlock(4) is already prepared to accept Q1 solutions
        CALL spdp_projectSolutionScalar (rvector%RvectorBlock(6),rprjVector%RvectorBlock(4))

        ! Calculate the derivative of the dual pressure into the first two
        ! subvectors -- overwriting the projected velocity field, which we don't need
        ! anymore.
        CALL ppgrd_calcGradient (rprjVector%RvectorBlock(4),rprjVector)
        
        ! Write pressure derivative field
        CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
        
        CALL ucd_addVariableVertexBased (rexport,'X-deriv-dpres',UCD_VAR_STANDARD, p_Ddata)
        CALL ucd_addVariableVertexBased (rexport,'Y-deriv-dpres',UCD_VAR_STANDARD, p_Ddata2)
            
      END IF
          
    END IF
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Release the auxiliary vector
    CALL lsysbl_releaseVector (rprjVector)
    
    ! Release the discretisation strucutre
    CALL spdiscr_releaseBlockDiscr(rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    CALL bcasm_releaseDiscreteBC (rdiscreteBC)
    CALL bcasm_releaseDiscreteFBC (rdiscreteFBC)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_postprocessingNonstat (rproblem,rvector)
  
!<description>
  ! Postprocessing of solutions of stationary simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
!</input>

!</subroutine>

  ! local variables
    INTEGER :: ieltype
  
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! A vector accepting Q1 data
    TYPE(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    TYPE(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    TYPE(t_discreteBC), TARGET :: rdiscreteBC
    TYPE(t_discreteFBC), TARGET :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    
    ! Forces on the object
    REAL(DP), DIMENSION(NDIM2D) :: Dforces
    REAL(DP) :: df1,df2
    TYPE(t_boundaryRegion) :: rregion
    
    ! Divergence
    TYPE(t_matrixScalar) :: rBmatrix
    TYPE(t_vectorScalar), TARGET :: rtempVector
    
    CHARACTER(SYS_STRLEN) :: sgmvName,stemp

    ! If we have a uniform discretisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    IF ((rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
         ccomplexity .EQ. SPDISC_UNIFORM) .AND. &
        (boundary_igetNBoundComp(rproblem%p_rboundary) .GE. 2)) THEN

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      CALL boundary_createRegion (rproblem%p_rboundary, &
          2, 0, rregion)
      rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      df1 = 1.0_DP/1000.0_DP
      df2 = 0.1_DP * 0.2_DP**2
      CALL ppns2D_bdforces_uniform (rvector,rregion,Dforces,CUB_G1_1D,df1,df2)
      
      CALL output_lbrk()
      CALL output_line ('Body forces real bd., bdc/horiz/vert')
      CALL output_line (' 2 / ' &
          //TRIM(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //TRIM(sys_sdEP(Dforces(2),15,6)) )
      
    ENDIF
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
!    IF (rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
!        ccomplexity .EQ. SPDISC_UNIFORM) THEN
!        
!      ieltype = rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
!                RelementDistribution(1)%itrialElement
!                
!      IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q1T) THEN
!      
!        ! Create a temporary vector 
!        CALL lsyssc_createVecByDiscr (rvector%RvectorBlock(3)%p_rspatialDiscretisation,&
!            rtempVector,.TRUE.)
!
!        ! Calculate divergence = B1^T u1 + B2^T u2
!        CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB1,&
!            rBmatrix,LSYSSC_TR_VIRTUAL)
!        CALL lsyssc_scalarMatVec (&
!            rBmatrix, rvector%RvectorBlock(1), &
!            rtempVector, 1.0_DP, 0.0_DP)
!        CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB2,&
!            rBmatrix,LSYSSC_TR_VIRTUAL)
!        CALL lsyssc_scalarMatVec (&
!            rBmatrix, rvector%RvectorBlock(2), &
!            rtempVector, 1.0_DP, 1.0_DP)
!        
!        CALL output_lbrk()
!        CALL output_line ('Divergence = ' &
!            //TRIM(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)) )
!            
!        CALL lsyssc_releaseVector (rtempVector)
!      
!      END IF
!      
!    END IF    
    
    ! The solution vector is probably not in the way, GMV likes it!
    ! GMV for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! GMV understands.
    ! GMV understands only Q1 solutions! So the task is now to create
    ! a Q1 solution from rvector and write that out.
    !
    ! For this purpose, first create a 'derived' simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    CALL spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscretisation,&
        rprjDiscretisation)
    
    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(1))

    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(2))
                 
    ! Also use Q1 for the dual velocity field.
    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(4), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(4))

    CALL spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscretisation(5), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(5))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    CALL lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.FALSE.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    CALL spdp_projectSolution (rvector,rprjVector)
    
    CALL cc_initCollectForAssembly(rproblem,rproblem%rcollection)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0 
    ! discretisation for implementing them into a solution vector.
    ! Discretise the boundary conditions for this discretisation
    CALL bcasm_initDiscreteBC(rdiscreteBC)
    CALL cc_assembleBDconditions (rproblem,rprjDiscretisation,&
      rdiscreteBC,rproblem%rcollection)

    CALL bcasm_initDiscreteFBC(rdiscreteFBC)
    CALL cc_assembleFBDconditions (rproblem,rprjDiscretisation,&
      rdiscreteFBC,rproblem%rcollection)
    
    CALL cc_doneCollectForAssembly(rproblem,rproblem%rcollection)
    
    ! Connect the vector with the BC's
    rprjVector%p_rdiscreteBC => rdiscreteBC
    rprjVector%p_rdiscreteBCfict => rdiscreteFBC
    
    ! Filter the solution vector to implement discrete BC's.
    CALL vecfil_discreteBCsol (rprjVector)
    !DEBUG: CALL vecfil_discreteBCdef (rprjVector)

    ! Filter the solution vector to implement discrete BC's for fictitious 
    ! boundary components.
    CALL vecfil_discreteFBCsol (rprjVector)
    
    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    !
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! Get GMV filename
    CALL parlst_getvalue_string_direct (rproblem%rparamList, 'TIME-POSTPROCESSING', &
                                        'sgmvFileName', stemp, 'gmv/u.gmv')
    READ(stemp,*) sgmvName
    
    IF (sgmvName .NE. '') THEN
    
      ! Start UCD export to GMV file:
      CALL output_lbrk ()
      CALL output_line ('Writing GMV file: ' &
          //TRIM(sgmvName)//'.'//sys_si0(rproblem%rtimedependence%itimeStep,5))
      
      CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
          TRIM(sgmvName)//'.'//sys_si0(rproblem%rtimedependence%itimeStep,5))
      
      ! Write the configuration of the application as comment block
      ! to the output file.
      CALL ucd_addCommentLine (rexport,'Configuration:')
      CALL ucd_addCommentLine (rexport,'---------------')
      CALL ucd_addParameterList (rexport,rproblem%rparamList)
      CALL ucd_addCommentLine (rexport,'---------------')

      ! Write velocity field
      CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      
      CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
      CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
      
      ! Write pressure
      CALL lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
      
      ! Dual velocity field
      CALL lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
      CALL lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
      CALL ucd_addVariableVertexBased (rexport,'X-vel-dual',UCD_VAR_STANDARD, p_Ddata)
      CALL ucd_addVariableVertexBased (rexport,'Y-vel-dual',UCD_VAR_STANDARD, p_Ddata2)
      
      ! Dual pressure
      CALL lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
      CALL ucd_addVariableElementBased (rexport,'pressure-dual',UCD_VAR_STANDARD, p_Ddata)
      
      ! If we have a simple Q1~ discretisation, calculate the streamfunction.
      IF (rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
          ccomplexity .EQ. SPDISC_UNIFORM) THEN
          
        ieltype = rvector%p_rblockDiscretisation%RspatialDiscretisation(1)% &
                  RelementDistribution(1)%itrialElement
                  
        IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q1T) THEN
            
          CALL ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
          
          CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
          CALL ucd_addVariableVertexBased (rexport,'streamfunction',&
              UCD_VAR_STANDARD, p_Ddata)
              
        END IF
        
      END IF
      
      ! Write the file to disc, that's it.
      CALL ucd_write (rexport)
      CALL ucd_release (rexport)
      
    END IF
    
    ! Release the auxiliary vector
    CALL lsysbl_releaseVector (rprjVector)
  
    ! Release the discretisations structure
    CALL spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    CALL bcasm_releaseDiscreteBC (rdiscreteBC)
    CALL bcasm_releaseDiscreteFBC (rdiscreteFBC)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_postprocSpaceTimeGMV (rproblem,rdiscr,rvector,sfilename)
  
!<description>
  ! For every sub-solution in the global space-time vector rvector,
  ! a GMV file is written to disc.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!<input>
  ! A space-time discretisation structure defining the discretisation of
  ! rvector.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rdiscr

  ! A space-time vector. For every timestep, a GMV is written.
  TYPE(t_spaceTimeVector), INTENT(IN) :: rvector
  
  ! A path + basic filename for the GMV-files. A number '.00000','.00001',...
  ! is appended for every timestep.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
!</input>

!</subroutine>

    ! local variables
    TYPE(t_vectorBlock) :: rvectorTmp
    INTEGER :: i,ieltype 
    
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! A vector accepting Q1 data
    TYPE(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    TYPE(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    TYPE(t_discreteBC), TARGET :: rdiscreteBC
    TYPE(t_discreteFBC), TARGET :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    
    ! Create a temp vector
    CALL lsysbl_createVecBlockByDiscr (&
        rdiscr%p_rlevelInfo%rdiscretisation,&
        rvectorTmp,.TRUE.)
    
    ! Attach the boundary conditions to that vector
    rvectorTmp%p_rdiscreteBC => rdiscr%p_rlevelInfo%p_rdiscreteBC
    rvectorTmp%p_rdiscreteBCfict => rdiscr%p_rlevelInfo%p_rdiscreteFBC
      
    ! Postprocessing of all solution vectors.
    DO i = 0,rdiscr%rtimeDiscr%nintervals
    
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + i*rdiscr%rtimeDiscr%dtstep
      rproblem%rtimedependence%itimeStep = i
    
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rvector, i, rvectorTmp)
      CALL tmevl_evaluate(rvector,rproblem%rtimedependence%dtime,rvectorTmp)
    
      ! The solution vector is probably not in the way, GMV likes it!
      ! GMV for example does not understand Q1~ vectors!
      ! Therefore, we first have to convert the vector to a form that
      ! GMV understands.
      ! GMV understands only Q1 solutions! So the task is now to create
      ! a Q1 solution from rvector and write that out.
      !
      ! For this purpose, first create a 'derived' simple discretisation
      ! structure based on Q1 by copying the main guiding block discretisation
      ! structure and modifying the discretisation structures of the
      ! two velocity subvectors:
      
      CALL spdiscr_duplicateBlockDiscr(rvectorTmp%p_rblockDiscretisation,&
          rprjDiscretisation)
      
      CALL spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscretisation%RspatialDiscretisation(1), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscretisation(1))

      CALL spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscretisation%RspatialDiscretisation(2), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscretisation(2))
                   
      ! Also use Q1 for the dual velocity field.
      CALL spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscretisation%RspatialDiscretisation(4), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscretisation(4))

      CALL spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscretisation%RspatialDiscretisation(5), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscretisation(5))
                   
      ! The pressure discretisation substructure stays the old.
      !
      ! Now set up a new solution vector based on this discretisation,
      ! allocate memory.
      CALL lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.FALSE.)
      
      ! Then take our original solution vector and convert it according to the
      ! new discretisation:
      CALL spdp_projectSolution (rvectorTmp,rprjVector)
      
      CALL cc_initCollectForAssembly(rproblem,rproblem%rcollection)
      
      ! Discretise the boundary conditions for this discretisation
      CALL bcasm_initDiscreteBC(rdiscreteBC)
      CALL cc_assembleBDconditions (rproblem,rprjDiscretisation,&
        rdiscreteBC,rproblem%rcollection)

      CALL bcasm_initDiscreteFBC(rdiscreteFBC)
      CALL cc_assembleFBDconditions (rproblem,rprjDiscretisation,&
        rdiscreteFBC,rproblem%rcollection)
      
      ! Connect the vector with the BC's
      rprjVector%p_rdiscreteBC => rdiscreteBC
      rprjVector%p_rdiscreteBCfict => rdiscreteFBC
      
      CALL cc_doneCollectForAssembly(rproblem,rproblem%rcollection)
      
      ! Filter the solution vector to implement discrete BC's.
      CALL vecfil_discreteBCsol (rprjVector)
      !DEBUG: CALL vecfil_discreteBCdef (rprjVector)

      ! Filter the solution vector to implement discrete BC's for fictitious 
      ! boundary components.
      CALL vecfil_discreteFBCsol (rprjVector)
      
      ! Now we have a Q1/Q1/Q0 solution in rprjVector.
      !
      ! From the attached discretisation, get the underlying triangulation
      p_rtriangulation => &
        rvectorTmp%rvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
      
      ! Start UCD export to GMV file:
      CALL output_line ('Writing GMV file: ' &
          //TRIM(sfilename)//'.'//sys_si0(rproblem%rtimedependence%itimeStep,5))
      
      CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
          TRIM(sfilename)//'.'//sys_si0(rproblem%rtimedependence%itimeStep,5))
      
      ! Write the configuration of the application as comment block
      ! to the output file.
      CALL ucd_addCommentLine (rexport,'Configuration:')
      CALL ucd_addCommentLine (rexport,'---------------')
      CALL ucd_addParameterList (rexport,rproblem%rparamList)
      CALL ucd_addCommentLine (rexport,'---------------')

      ! Write velocity field
      CALL lsyssc_getbase_double (rprjVector%rvectorBlock(1),p_Ddata)
      CALL lsyssc_getbase_double (rprjVector%rvectorBlock(2),p_Ddata2)
      
      CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
      CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
      
      ! Write pressure
      CALL lsyssc_getbase_double (rprjVector%rvectorBlock(3),p_Ddata)
      CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
      
      ! Dual velocity field
      CALL lsyssc_getbase_double (rprjVector%rvectorBlock(4),p_Ddata)
      CALL lsyssc_getbase_double (rprjVector%rvectorBlock(5),p_Ddata2)
      CALL ucd_addVariableVertexBased (rexport,'X-vel-dual',UCD_VAR_STANDARD, p_Ddata)
      CALL ucd_addVariableVertexBased (rexport,'Y-vel-dual',UCD_VAR_STANDARD, p_Ddata2)
      
      ! Dual pressure
      CALL lsyssc_getbase_double (rprjVector%rvectorBlock(6),p_Ddata)
      CALL ucd_addVariableElementBased (rexport,'pressure-dual',UCD_VAR_STANDARD, p_Ddata)
      
      ! If we have a simple Q1~ discretisation, calculate the streamfunction.
      IF (rvectorTmp%p_rblockDiscretisation%RspatialDiscretisation(1)% &
          ccomplexity .EQ. SPDISC_UNIFORM) THEN
          
        ieltype = rvectorTmp%p_rblockDiscretisation%RspatialDiscretisation(1)% &
                  RelementDistribution(1)%itrialElement
                  
        IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q1T) THEN
            
          CALL ppns2D_streamfct_uniform (rvectorTmp,rprjVector%rvectorBlock(1))
          
          CALL lsyssc_getbase_double (rprjVector%rvectorBlock(1),p_Ddata)
          CALL ucd_addVariableVertexBased (rexport,'streamfunction',&
              UCD_VAR_STANDARD, p_Ddata)
              
        END IF
        
      END IF
      
      ! Write the file to disc, that's it.
      CALL ucd_write (rexport)
      CALL ucd_release (rexport)
      
      ! Release the auxiliary vector
      CALL lsysbl_releaseVector (rprjVector)
      
      ! Release the discretisation structure
      CALL spdiscr_releaseBlockDiscr(rprjDiscretisation)
      
      ! Throw away the discrete BC's - not used anymore.
      CALL bcasm_releaseDiscreteBC (rdiscreteBC)
      CALL bcasm_releaseDiscreteFBC (rdiscreteFBC)
      
    END DO

    ! Release memory.
    CALL lsysbl_releaseVector (rvectorTmp)
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE cc_initpostprocessing (rproblem,rpostprocessing)

!<description>
  ! Initialises the given postprocessing structure rpostprocessing
  ! according to the main problem rproblem. The structure can then be used
  ! to generate postprocessing data.
!</description>

!<input>
  ! A problem structure that describes the main problem to solve.
  TYPE(t_problem), INTENT(IN),TARGET :: rproblem
!</input>

!<output>  
  TYPE(t_c2d2postprocessing), INTENT(OUT) :: rpostprocessing
!</output>

!</subroutine>

  ! local variables
  TYPE(t_blockDiscretisation), POINTER :: p_rdiscr

    ! For postprocessing, we need discretisation structures in the Q0 and Q1 space,
    ! later perhaps in the Q2 space. For this purpose, derive the corresponding
    ! discretisation structure using the 'main' discretisation structure on the
    ! maximum level.
    !
    ! For simplicity, we use only the discretisation structure of the X-velocity
    ! to derive everything.
    
    p_rdiscr => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation

    ! Piecewise constant space:
    CALL spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscretisation(1), &
                 EL_Q0, CUB_G1X1, &
                 rpostprocessing%rdiscrConstant)

    ! Piecewise linear space:
    CALL spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscretisation(1), &
                 EL_Q1, CUB_G2X2, &
                 rpostprocessing%rdiscrLinear)
  
    ! Piecewise quadratic space:
    CALL spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscretisation(1), &
                 EL_Q2, CUB_G3X3, &
                 rpostprocessing%rdiscrQuadratic)
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>

  SUBROUTINE cc_clearpostprocessing (rpostprocessing)

!<description>
  ! Releases all calculated data in the given postprocessing structure
  ! so that it can be allocated again in the calculation routines.
  ! This routine must be called at the end of the postprocessing routines
  ! to release the temporary memory that was allocated for the vectors
  ! in the postprocessing structure.
!</description>

!<inputoutput>  
  TYPE(t_c2d2postprocessing), INTENT(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors which might be allocated.
    CALL lsyssc_releaseVector (rpostprocessing%rvectorVelX)
    CALL lsyssc_releaseVector (rpostprocessing%rvectorVelY)
    CALL lsyssc_releaseVector (rpostprocessing%rvectorPressure)
    CALL lsyssc_releaseVector (rpostprocessing%rvectorPressureCells)
    CALL lsyssc_releaseVector (rpostprocessing%rvectorStreamfunction)
    CALL lsyssc_releaseVector (rpostprocessing%rvectorH1err)
    CALL lsyssc_releaseVector (rpostprocessing%rvectorH1errCells)

  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>

  SUBROUTINE cc_donepostprocessing (rpostprocessing)

!<description>
  ! Releases a given problem structure. All allocated memory of this structure
  ! is released.
!</description>

!<inputoutput>  
  TYPE(t_c2d2postprocessing), INTENT(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors -- if there are still some allocated
    CALL cc_clearpostprocessing (rpostprocessing)

    ! Release the discretisation structures allocated above.
    CALL spdiscr_releaseDiscr(rpostprocessing%rdiscrQuadratic)
    CALL spdiscr_releaseDiscr(rpostprocessing%rdiscrLinear)
    CALL spdiscr_releaseDiscr(rpostprocessing%rdiscrConstant)
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_printControlFunctionalStat (rproblem,rvector)
  
!<description>
  ! Calculates and prints the value of the optimal control functional J(y,u) 
  ! in the stationary case.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
!</input>

!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(3) :: Derror
    REAL(DP) :: dalphaC

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Analyse the deviation from the target velocity field.
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',dalphaC,0.1_DP)
    CALL cc_optc_stationaryFunctional (rvector,dalphaC,Derror,&
        rproblem%rcollection)

    CALL output_line ('||y-z||_L2: '//TRIM(sys_sdEL(Derror(1),2)))
    CALL output_line ('||u||_L2  : '//TRIM(sys_sdEL(Derror(2),2)))
    CALL output_line ('J(y,u)    : '//TRIM(sys_sdEL(Derror(3),2)))
    
    ! Release assembly-stuff from the collection
    CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
    
  END SUBROUTINE

!  !****************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE ppns2D_fbdforces (rvector,rregion,Dforces,cformulation,CcubU,CcubP)
!
!!<description>
!  ! Calculates the drag-/lift-forces acting on a part of the real
!  ! boundary for a vector rvector with the solution of the 2D
!  ! (Navier-)Stokes equation. It's assumed that
!  !  rvector%rvectorBlock(1) = X-velocity,
!  !  rvector%rvectorBlock(2) = Y-velocity,
!  !  rvector%rvectorBlock(3) = pressure,
!  ! and that X- and Y-velocity is discretised with the same 
!  ! finite element.
!  ! rregion specifies a boundary region where to calculate
!  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
!  ! X- and Y-direction.
!  !
!  ! Double precision version
!!</description>
!
!!<input>
!  ! The FE solution vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!  
!  ! Boundary region where to calculate the boundary forces.
!  ! Can be created e.g. by boundary_createRegion.
!  TYPE(t_boundaryRegion), INTENT(OUT) :: rregion
!  
!  ! Type identifier that specifies the formulation of the Navier Stokes
!  ! equation that should be used for performing the integration.
!  ! ... = gradient formulation,
!  ! ... = deformation formulation
!  INTEGER, INTENT(IN) :: cformulation
!  
!  ! OPTIONAL: Array with cubature formula identifiers when integrating
!  ! the velocity. If specified, there must be a cubature formula identifier
!  ! for each of the element distributions specified in the discretisation
!  ! structure that is associated to rvector.
!  ! If not specified, the standard cubature formula in the discretisation
!  ! structure for integrating linear forms is used.
!  INTEGER, DIMENSION(:), INTENT(IN) :: CcubU
!
!  ! OPTIONAL: Array with cubature formula identifiers when integrating
!  ! the pressure. If specified, there must be a cubature formula identifier
!  ! for each of the element distributions specified in the discretisation
!  ! structure that is associated to rvector.
!  ! If not specified, the standard cubature formula in the discretisation
!  ! structure for integrating linear forms is used.
!  INTEGER, DIMENSION(:), INTENT(IN) :: CcubP
!!</input>
!
!!<output>
!  ! Array receiving the forces acting on the boundary specified by rregion.
!  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
!  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dforces
!!</output>
!
!!</subroutine>
!
!  ! Array to tell the element which derivatives to calculate.
!  LOGICAL, DIMENSION(EL_MAXNDER) :: BderU, BderP
!  
!  ! Cubature point coordinates on the reference element.
!  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: DxiU, DxiP
!
!  ! For every cubature point on the reference element,
!  ! the corresponding cubature weight
!  REAL(DP), DIMENSION(CUB_MAXCUBP) :: DomegaU, DomegaP
!  
!  ! number of cubature points on the reference element
!  INTEGER :: ncubpU, ncubpP
!
!  ! Pointer to the vector entries
!  REAL(DP), DIMENSION(:), POINTER :: p_DdataUX, p_DdataUY, p_DdataP
!
!  ! An allocateable array accepting the DOF's of a set of elements.
!  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialU, IdofsTrialP
!  
!  ! Allocateable arrays for the values of the basis functions.
!  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTrialU,DbasTrialP
!  
!  ! Number of entries in the vector - for quicker access
!  INTEGER(I32) :: neqU, neqP
!
!  ! Arrays for saving Jacobian determinants and matrices
!  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
!
!  ! Pointer to KVERT of the triangulation
!  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
!  
!  ! Pointer to DCORVG of the triangulation
!  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!  
!  ! Current element distribution
!  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
!
!  ! Number of elements in a block. Normally =BILF_NELEMSIM,
!  ! except if there are less elements in the discretisation.
!  INTEGER :: nelementsPerBlock
!  
!
!  BderU = .FALSE.
!  BderP = .FALSE.
!
!  ! Get the vector data
!  neqU = rvector%RvectorBlock(1)%NEQ
!  neqP = rvector%RvectorBlock(3)%NEQ
!  
!  IF (rvector%cdataType .NE. ST_DOUBLE) THEN
!    PRINT *,'ppns2D_fbdforces: Unsupported vector precision.'
!    STOP
!  END IF
!
!  ! We support only uniform and conformal discretisation structures.
!  IF (.NOT. ASSOCIATED(rvector%p_rblockDiscretisation)) THEN
!    PRINT *,'ppns2D_fbdforces: No discretisation structure!'
!    STOP
!  END IF
!
!  IF ((rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .AND. &
!      (rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
!    PRINT *,'ppns2D_fbdforces: Discretisation too complex!'
!    STOP
!  END IF
!  
!  ! Get pointers to the subvectors from the block vector
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)
!
!  ! Now loop over the different element distributions in the discretisation.
!  ! Each element distribution has a different element combination.
!  ! E.g. distribution 1 may describe P1~/P1~/P0 triangular elements, while
!  ! distrbution 2 may describe Q1~/Q1~/Q0.
!
!  !DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
!
!
!  !END DO ! icurrentElementDistr
!  
!  Dforces = 0.0_DP
!  
!  END SUBROUTINE
!
!
!  !****************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE ppns2D_bdforces (rvector,rregion,Dforces,cformulation,CcubU,CcubP)
!
!!<description>
!  ! Calculates the drag-/lift-forces acting on a part of the real
!  ! boundary for a vector rvector with the solution of the 2D
!  ! (Navier-)Stokes equation. It's assumed that
!  !  rvector%rvectorBlock(1) = X-velocity,
!  !  rvector%rvectorBlock(2) = Y-velocity,
!  !  rvector%rvectorBlock(3) = pressure,
!  ! and that X- and Y-velocity is discretised with the same 
!  ! finite element.
!  ! rregion specifies a boundary region where to calculate
!  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
!  ! X- and Y-direction.
!  !
!  ! Double precision version
!!</description>
!
!!<input>
!  ! The FE solution vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!  
!  ! Boundary region where to calculate the boundary forces.
!  ! Can be created e.g. by boundary_createRegion.
!  TYPE(t_boundaryRegion), INTENT(OUT) :: rregion
!  
!  ! Type identifier that specifies the formulation of the Navier Stokes
!  ! equation that should be used for performing the integration.
!  ! ... = gradient formulation,
!  ! ... = deformation formulation
!  INTEGER, INTENT(IN) :: cformulation
!  
!  ! 1D Cubature formula identifier to use for the line integration.
!  ! One of the CUB_xxxx_1D constants in the cubature.f90.
!  INTEGER, DIMENSION(:), INTENT(IN) :: CcubU
!!</input>
!
!!<output>
!  ! Array receiving the forces acting on the boundary specified by rregion.
!  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
!  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dforces
!!</output>
!
!!</subroutine>
!
!  ! Spatial discretisation structure of velocity and pressure
!  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrU, p_rdiscrP
!  
!  ! Element distribution of velocity and pressure
!  TYPE(t_elementDistribution), POINTER :: p_elemDistrU,p_elemDistrP
!
!
!  ! Array to tell the element which derivatives to calculate.
!  LOGICAL, DIMENSION(EL_MAXNDER) :: BderU, BderP
!  
!  ! Cubature point coordinates on the reference element.
!  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D, DXi2D
!
!  ! For every cubature point on the reference element,
!  ! the corresponding cubature weight
!  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
!  
!  ! number of cubature points on the reference element
!  INTEGER :: ncubp
!
!  ! Pointer to the vector entries
!  REAL(DP), DIMENSION(:), POINTER :: p_DdataUX, p_DdataUY, p_DdataP
!
!  ! An allocateable array accepting the DOF's of a set of elements.
!  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialU, IdofsTrialP
!  
!  ! Allocateable arrays for the values of the basis functions.
!  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTrialU,DbasTrialP
!  
!  ! Number of entries in the vector - for quicker access
!  INTEGER(I32) :: neqU, neqP
!
!  ! Arrays for saving Jacobian determinants and matrices
!  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
!
!  ! Pointer to KVERT of the triangulation
!  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
!  
!  ! Pointer to DCORVG of the triangulation
!  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!  
!
!  ! Number of elements in a block. Normally =BILF_NELEMSIM,
!  ! except if there are less elements in the discretisation.
!  INTEGER :: nelementsPerBlock
!  
!  BderU = .FALSE.
!  BderP = .FALSE.
!
!  ! Get the vector data
!  neqU = rvector%RvectorBlock(1)%NEQ
!  neqP = rvector%RvectorBlock(3)%NEQ
!  
!  IF (rvector%cdataType .NE. ST_DOUBLE) THEN
!    PRINT *,'ppns2D_bdforces: Unsupported vector precision.'
!    STOP
!  END IF
!
!  ! We support only uniform and conformal discretisation structures.
!  IF (.NOT. ASSOCIATED(rvector%p_rblockDiscretisation)) THEN
!    PRINT *,'ppns2D_bdforces: No discretisation structure!'
!    STOP
!  END IF
!
!  IF ((rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .AND. &
!      (rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
!    PRINT *,'ppns2D_bdforces: Discretisation too complex!'
!    STOP
!  END IF
!  
!  ! Get pointers to the subvectors from the block vector
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)
!  
!  ! Get pointers to the spatial discretisation structures of the
!  ! veloctiy and pressure
!  p_rdiscrU => rvector%RvectorBlock(1)%p_spatialDiscretisation
!  p_rdiscrP => rvector%RvectorBlock(3)%p_spatialDiscretisation
!
!  ! Now loop over the different element distributions in the discretisation.
!  ! Each element distribution has a different element combination.
!  ! E.g. distribution 1 may describe P1~/P1~/P0 triangular elements, while
!  ! distrbution 2 may describe Q1~/Q1~/Q0.
!  !
!  ! For every element distribution, we have to loop through all elements
!  ! along the boundary component, calculate their contribution to the
!  ! integral and add everything together.
!
!  DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
!
!    ! Get the element distribution structure of that FE-space comination
!    p_relemDistrU => p_rdiscrU%RelementDistribution(RelementDistribution)
!    p_relemDistrP => p_rdiscrP%RelementDistribution(RelementDistribution)
!    
!    ! By that we know
!    ! - what's the element
!    ! - which transformation from the reference element to use
!    
!    
!    
!    
!
!  END DO ! icurrentElementDistr
!  
!  END SUBROUTINE

END MODULE
