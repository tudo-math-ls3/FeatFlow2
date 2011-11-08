!##############################################################################
!# ****************************************************************************
!# <name> postprocessing </name>
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

module postprocessing

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
  use pprocerror
  use pprocgradients
  
  use collection
  use convection
  
  use ucd
  
  use pprocnavierstokes
  
  use basicstructures
  use optcanalysis
  use user_callback
  use spatialbcdef
  use spacematvecassembly
  
  use spacetimediscretisation
  
  implicit none
  
!<types>

!<typeblock>

  ! Postprocessing structure. Whenever the postprocessing calculation routines are
  ! called, they fills this structure using the current solution vector (and other
  ! information if necessary). The information in this structure can then be used
  ! for GMV output e.g.
  type t_c2d2postprocessing

    ! A discretisation structure that describes a piecewise constant discretisation
    ! (usually P0 or Q0).
    type(t_spatialDiscretisation) :: rdiscrConstant
    
    ! A discretisation structure that describes a piecewise linear discretisation
    ! (usually P1 or Q1).
    type(t_spatialDiscretisation) :: rdiscrLinear

    ! A discretisation structure that describes a piecewise quadratic discretisation
    ! (usually P2 or Q2).
    type(t_spatialDiscretisation) :: rdiscrQuadratic

    ! A vector that describes the X-velocity field in the vertices
    type(t_vectorScalar) :: rvectorVelX

    ! A vector that describes the Y-velocity field in the vertices
    type(t_vectorScalar) :: rvectorVelY

    ! A vector that describes the pressure field in the vertices
    type(t_vectorScalar) :: rvectorPressure

    ! A vector that describes the pressure field in the cells
    type(t_vectorScalar) :: rvectorPressureCells

    ! A vector that describes the streamfunction
    type(t_vectorScalar) :: rvectorStreamfunction
    
    ! A vector that describes the H1-error of the velocity field in the vertices
    type(t_vectorScalar) :: rvectorH1err

    ! A vector that describes the H1-error of the pressure in the cells
    type(t_vectorScalar) :: rvectorH1errCells
  
  end type

!</typeblock>

!</types>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_postprocessingStationary (rproblem,rvector)
  
!<description>
  ! Postprocessing of solutions of stationary simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!</subroutine>

  ! local variables
    integer :: ieltype
  
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    type(t_discreteBC), target :: rdiscreteBC
    type(t_discreteFBC), target :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    ! Forces on the object
    real(DP), dimension(NDIM2D) :: Dforces
    real(DP) :: df1,df2
    type(t_boundaryRegion) :: rregion
    
    ! Divergence
    type(t_vectorScalar), target :: rtempVector

    ! If we have a uniform discreisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    if ((rvector%p_rblockDiscr%RspatialDiscr(1)% &
         ccomplexity .eq. SPDISC_UNIFORM) .and. &
        (boundary_igetNBoundComp(rproblem%rboundary) .ge. 2)) then

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      call boundary_createRegion (rproblem%rboundary, &
          2, 0, rregion)
      rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      df1 = 1.0_DP/1000.0_DP
      df2 = 0.1_DP * 0.2_DP**2
      call ppns2D_bdforces_uniform (rvector,rregion,Dforces,CUB_G1_1D,df1,df2)

      call output_lbrk()
      call output_line ('Body forces real bd., bdc/horiz/vert')
      call output_line (' 2 / ' &
          //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //trim(sys_sdEP(Dforces(2),15,6)) )
      
    endif
    
    ! Print out the value of the optimal control functional.
    call output_lbrk ()
    call cc_printControlFunctionalStat (rproblem,0.0_DP,rvector)
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    if (rvector%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(1)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
      
        ! Create a temporary vector
        call lsyssc_createVecByDiscr (rvector%RvectorBlock(3)%p_rspatialDiscr,&
            rtempVector,.true.)

        ! Calculate divergence = B1^T u1 + B2^T u2
        call lsyssc_scalarMatVec (&
            rproblem%RlevelInfo(rproblem%nlmax)%rstaticInfo%rmatrixD1, rvector%RvectorBlock(1), &
            rtempVector, 1.0_DP, 0.0_DP)
        call lsyssc_scalarMatVec (&
            rproblem%RlevelInfo(rproblem%nlmax)%rstaticInfo%rmatrixD2, rvector%RvectorBlock(2), &
            rtempVector, 1.0_DP, 1.0_DP)
        
        call output_lbrk()
        call output_line ('Divergence = ' &
            //trim(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)) )

        call lsyssc_releaseVector (rtempVector)
      
      end if
      
    end if
    
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
    
    call spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscr,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(2))

    ! Also use Q1 for the dual velocity field.
    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(4), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(4))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(5), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(5))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)
    
    ! Discretise the boundary conditions for this discretisation
    call bcasm_initDiscreteBC(rdiscreteBC)
    call cc_assembleBDconditions (rproblem,0.0_DP,rprjDiscretisation,CCDISCBC_PRIMALDUAL,&
      rdiscreteBC,rproblem%rcollection)

    call bcasm_initDiscreteFBC(rdiscreteFBC)
    call cc_assembleFBDconditions (rproblem,0.0_DP,rprjDiscretisation,CCDISCBC_PRIMALDUAL,&
      rdiscreteFBC,rproblem%rcollection)
    
    ! Connect the vector with the BC's
    call lsysbl_assignDiscreteBC (rprjVector,rdiscreteBC)
    call lsysbl_assignDiscreteFBC (rprjVector,rdiscreteFBC)
    
    ! Filter the solution vector to implement discrete BC's.
    call vecfil_discreteBCsol (rprjVector)

    ! Filter the solution vector to implement discrete BC's for fictitious
    ! boundary components.
    call vecfil_discreteFBCsol (rprjVector)
    
    ! Clean up the collection (as we are done with the assembly.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    !
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,'gmv/u.gmv')
    
    ! Write the configuration of the application as comment block
    ! to the output file.
    call ucd_addCommentLine (rexport,'Configuration:')
    call ucd_addCommentLine (rexport,'---------------')
    call ucd_addParameterList (rexport,rproblem%rparamList)
    call ucd_addCommentLine (rexport,'---------------')

    ! Write velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    call ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    call ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)

    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)

    ! Dual velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
    call ucd_addVariableVertexBased (rexport,'X-vel-dual',UCD_VAR_STANDARD, p_Ddata)
    call ucd_addVariableVertexBased (rexport,'Y-vel-dual',UCD_VAR_STANDARD, p_Ddata2)
    
    ! Dual pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure-dual',UCD_VAR_STANDARD, p_Ddata)
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    if (rvector%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(1)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
          
        call ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
        
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'streamfunction',&
            UCD_VAR_STANDARD, p_Ddata)

      end if
      
    end if
    
    ! Calculate the derivative of the pressure.
    ! First, project the pressure from Q0 to Q1 (if it's Q0).
    if (rvector%p_rblockDiscr%RspatialDiscr(3)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(3)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q0) then
        ! rprjVector%RvectorBlock(4) is already prepared to accept Q1 solutions
        call spdp_projectSolutionScalar (rvector%RvectorBlock(3),rprjVector%RvectorBlock(4))

        ! Calculate the derivative of the pressure into the first two
        ! subvectors -- overwriting the projected velocity field, which we don't need
        ! anymore.
        call ppgrd_calcGradient (rprjVector%RvectorBlock(4),rprjVector)
        
        ! Write pressure derivative field
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
        
        call ucd_addVariableVertexBased (rexport,'X-deriv-pres',UCD_VAR_STANDARD, p_Ddata)
        call ucd_addVariableVertexBased (rexport,'Y-deriv-pres',UCD_VAR_STANDARD, p_Ddata2)
            
      end if
          
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(6)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q0) then
        ! rprjVector%RvectorBlock(4) is already prepared to accept Q1 solutions
        call spdp_projectSolutionScalar (rvector%RvectorBlock(6),rprjVector%RvectorBlock(4))

        ! Calculate the derivative of the dual pressure into the first two
        ! subvectors -- overwriting the projected velocity field, which we don't need
        ! anymore.
        call ppgrd_calcGradient (rprjVector%RvectorBlock(4),rprjVector)
        
        ! Write pressure derivative field
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
        
        call ucd_addVariableVertexBased (rexport,'X-deriv-dpres',UCD_VAR_STANDARD, p_Ddata)
        call ucd_addVariableVertexBased (rexport,'Y-deriv-dpres',UCD_VAR_STANDARD, p_Ddata2)
            
      end if
          
    end if
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
    
    ! Release the discretisation strucutre
    call spdiscr_releaseBlockDiscr(rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    call bcasm_releaseDiscreteBC (rdiscreteBC)
    call bcasm_releaseDiscreteFBC (rdiscreteFBC)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_postprocessingNonstat (rproblem,itimestep,dtime,rvector)
  
!<description>
  ! Postprocessing of solutions of stationary simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
  
  ! Number of the current timestep.
  integer, intent(in) :: itimestep
  
  ! Current simulation time.
  real(dp), intent(in) :: dtime
!</input>

!</subroutine>

  ! local variables
    integer :: ieltype
  
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    type(t_discreteBC), target :: rdiscreteBC
    type(t_discreteFBC), target :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    ! Forces on the object
    real(DP), dimension(NDIM2D) :: Dforces
    real(DP) :: df1,df2
    type(t_boundaryRegion) :: rregion
    
    ! Divergence
    type(t_vectorScalar), target :: rtempVector
    
    character(SYS_STRLEN) :: sgmvName,stemp

    ! If we have a uniform discretisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    if ((rvector%p_rblockDiscr%RspatialDiscr(1)% &
         ccomplexity .eq. SPDISC_UNIFORM) .and. &
        (boundary_igetNBoundComp(rproblem%rboundary) .ge. 2)) then

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      call boundary_createRegion (rproblem%rboundary, &
          2, 0, rregion)
      rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      df1 = 1.0_DP/1000.0_DP
      df2 = 0.1_DP * 0.2_DP**2
      call ppns2D_bdforces_uniform (rvector,rregion,Dforces,CUB_G1_1D,df1,df2)
      
      call output_lbrk()
      call output_line ('Body forces real bd., bdc/horiz/vert')
      call output_line (' 2 / ' &
          //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //trim(sys_sdEP(Dforces(2),15,6)) )
      
    endif
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
!    IF (rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
!        ccomplexity .EQ. SPDISC_UNIFORM) THEN
!
!      ieltype = rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
!                RelementDistr(1)%itrialElement
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
    
    call spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscr,&
        rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(2))
                 
    ! Also use Q1 for the dual velocity field.
    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(4), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(4))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(5), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(5))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)
    
    call cc_initCollectForAssembly(rproblem,dtime,rproblem%rcollection)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0
    ! discretisation for implementing them into a solution vector.
    ! Discretise the boundary conditions for this discretisation
    call bcasm_initDiscreteBC(rdiscreteBC)
    call cc_assembleBDconditions (rproblem,dtime,rprjDiscretisation,CCDISCBC_PRIMALDUAL,&
      rdiscreteBC,rproblem%rcollection)

    call bcasm_initDiscreteFBC(rdiscreteFBC)
    call cc_assembleFBDconditions (rproblem,dtime,rprjDiscretisation,CCDISCBC_PRIMALDUAL,&
      rdiscreteFBC,rproblem%rcollection)
    
    call cc_doneCollectForAssembly(rproblem,rproblem%rcollection)
    
    ! Connect the vector with the BC's
    call lsysbl_assignDiscreteBC (rprjVector,rdiscreteBC)
    call lsysbl_assignDiscreteFBC (rprjVector,rdiscreteFBC)
    
    ! Filter the solution vector to implement discrete BC's.
    call vecfil_discreteBCsol (rprjVector)
    !DEBUG: CALL vecfil_discreteBCdef (rprjVector)

    ! Filter the solution vector to implement discrete BC's for fictitious
    ! boundary components.
    call vecfil_discreteFBCsol (rprjVector)
    
    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    !
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! Get GMV filename
    call parlst_getvalue_string (rproblem%rparamList, 'TIME-POSTPROCESSING', &
                                 'sgmvFileName', stemp, 'gmv/u.gmv')
    read(stemp,*) sgmvName
    
    if (sgmvName .ne. '') then
    
      ! Start UCD export to GMV file:
      call output_lbrk ()
      call output_line ('Writing GMV file: ' &
          //trim(sgmvName)//'.'//sys_si0(itimeStep,5))
      
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
          trim(sgmvName)//'.'//sys_si0(itimeStep,5))
      
      ! Write the configuration of the application as comment block
      ! to the output file.
      call ucd_addCommentLine (rexport,'Configuration:')
      call ucd_addCommentLine (rexport,'---------------')
      call ucd_addParameterList (rexport,rproblem%rparamList)
      call ucd_addCommentLine (rexport,'---------------')

      ! Write velocity field
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      
      call ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
      
      ! Write pressure
      call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
      
      ! Dual velocity field
      call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
      call ucd_addVariableVertexBased (rexport,'X-vel-dual',UCD_VAR_STANDARD, p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Y-vel-dual',UCD_VAR_STANDARD, p_Ddata2)
      
      ! Dual pressure
      call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
      call ucd_addVariableElementBased (rexport,'pressure-dual',UCD_VAR_STANDARD, p_Ddata)

      ! Control u = P[min/max](-1/alpha lambda)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
      call lsyssc_scaleVector (rprjVector%RvectorBlock(4),-1.0_DP/rproblem%roptControl%dalphaC)
      if (rproblem%roptControl%ccontrolConstraints .ne. 0) then
        call cc_projectControlTimestep (rprjVector%RvectorBlock(4),&
          rproblem%roptControl%dumin1,rproblem%roptControl%dumax1)
      end if
      call ucd_addVariableVertexBased (rexport,'X-control',UCD_VAR_STANDARD, p_Ddata)

      call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
      call lsyssc_scaleVector (rprjVector%RvectorBlock(5),-1.0_DP/rproblem%roptControl%dalphaC)
      if (rproblem%roptControl%ccontrolConstraints .ne. 0) then
        call cc_projectControlTimestep (rprjVector%RvectorBlock(5),&
          rproblem%roptControl%dumin2,rproblem%roptControl%dumax2)
      end if
      call ucd_addVariableVertexBased (rexport,'Y-control',UCD_VAR_STANDARD, p_Ddata2)
      
      ! If we have a simple Q1~ discretisation, calculate the streamfunction.
      if (rvector%p_rblockDiscr%RspatialDiscr(1)% &
          ccomplexity .eq. SPDISC_UNIFORM) then
          
        ieltype = rvector%p_rblockDiscr%RspatialDiscr(1)% &
                  RelementDistr(1)%celement
                  
        if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
            
          call ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
          
          call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
          call ucd_addVariableVertexBased (rexport,'streamfunction',&
              UCD_VAR_STANDARD, p_Ddata)
              
        end if
        
      end if
      
      ! Write the file to disc, that's it.
      call ucd_write (rexport)
      call ucd_release (rexport)
      
    end if
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
  
    ! Release the discretisations structure
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    call bcasm_releaseDiscreteBC (rdiscreteBC)
    call bcasm_releaseDiscreteFBC (rdiscreteFBC)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_postprocSpaceTimeGMV (rproblem,rdiscr,rvector,sfilename)
  
!<description>
  ! For every sub-solution in the global space-time vector rvector,
  ! a GMV file is written to disc.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!<input>
  ! A space-time discretisation structure defining the discretisation of
  ! rvector.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rdiscr

  ! A space-time vector. For every timestep, a GMV is written.
  type(t_spaceTimeVector), intent(IN) :: rvector
  
  ! A path + basic filename for the GMV-files. A number '.00000','.00001',...
  ! is appended for every timestep.
  character(LEN=*), intent(IN) :: sfilename
!</input>

!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvectorTmp
    integer :: i,ieltype
    real(dp) :: dtime
    
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    type(t_discreteBC), target :: rdiscreteBC
    type(t_discreteFBC), target :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    ! Create a temp vector
    call lsysbl_createVecBlockByDiscr (&
        rdiscr%p_rlevelInfo%rdiscretisation,&
        rvectorTmp,.true.)
    
    ! Postprocessing of all solution vectors.
    do i = 0,rdiscr%rtimeDiscr%nintervals
    
      dtime = rproblem%rtimedependence%dtimeInit + i*rdiscr%rtimeDiscr%dtstep
    
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rvector, i, rvectorTmp)
      call tmevl_evaluate(rvector,dtime,rvectorTmp)
    
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
      
      call spdiscr_duplicateBlockDiscr(rvectorTmp%p_rblockDiscr,&
          rprjDiscretisation)
      
      call spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscr%RspatialDiscr(1), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscr(1))

      call spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscr%RspatialDiscr(2), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscr(2))
                   
      ! Also use Q1 for the dual velocity field.
      call spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscr%RspatialDiscr(4), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscr(4))

      call spdiscr_deriveSimpleDiscrSc (&
                  rvectorTmp%p_rblockDiscr%RspatialDiscr(5), &
                  EL_Q1, CUB_G2X2, &
                  rprjDiscretisation%RspatialDiscr(5))
                   
      ! The pressure discretisation substructure stays the old.
      !
      ! Now set up a new solution vector based on this discretisation,
      ! allocate memory.
      call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
      
      ! Then take our original solution vector and convert it according to the
      ! new discretisation:
      call spdp_projectSolution (rvectorTmp,rprjVector)
      
      call cc_initCollectForAssembly(rproblem,dtime,rproblem%rcollection)
      
      ! Discretise the boundary conditions for this discretisation
      call bcasm_initDiscreteBC(rdiscreteBC)
      call cc_assembleBDconditions (rproblem,dtime,rprjDiscretisation,CCDISCBC_PRIMALDUAL,&
        rdiscreteBC,rproblem%rcollection)

      call bcasm_initDiscreteFBC(rdiscreteFBC)
      call cc_assembleFBDconditions (rproblem,dtime,rprjDiscretisation,CCDISCBC_PRIMALDUAL,&
        rdiscreteFBC,rproblem%rcollection)
      
      ! Connect the vector with the BC's
      call lsysbl_assignDiscreteBC (rprjVector,rdiscreteBC)
      call lsysbl_assignDiscreteFBC (rprjVector,rdiscreteFBC)
      
      call cc_doneCollectForAssembly(rproblem,rproblem%rcollection)
      
      ! Filter the solution vector to implement discrete BC's.
      call vecfil_discreteBCsol (rprjVector)
      !DEBUG: CALL vecfil_discreteBCdef (rprjVector)

      ! Filter the solution vector to implement discrete BC's for fictitious
      ! boundary components.
      call vecfil_discreteFBCsol (rprjVector)
      
      ! Now we have a Q1/Q1/Q0 solution in rprjVector.
      !
      ! From the attached discretisation, get the underlying triangulation
      p_rtriangulation => &
        rvectorTmp%rvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
      
      ! Start UCD export to GMV file:
      call output_line ('Writing GMV file: ' &
          //trim(sfilename)//'.'//sys_si0(i,5))
      
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
          trim(sfilename)//'.'//sys_si0(i,5))
      
      ! Write the configuration of the application as comment block
      ! to the output file.
      call ucd_addCommentLine (rexport,'Configuration:')
      call ucd_addCommentLine (rexport,'---------------')
      call ucd_addParameterList (rexport,rproblem%rparamList)
      call ucd_addCommentLine (rexport,'---------------')

      ! Write velocity field
      call lsyssc_getbase_double (rprjVector%rvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%rvectorBlock(2),p_Ddata2)
      
      call ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
      
      ! Write pressure
      call lsyssc_getbase_double (rprjVector%rvectorBlock(3),p_Ddata)
      call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
      
      ! Dual velocity field
      call lsyssc_getbase_double (rprjVector%rvectorBlock(4),p_Ddata)
      call lsyssc_getbase_double (rprjVector%rvectorBlock(5),p_Ddata2)
      call ucd_addVariableVertexBased (rexport,'X-vel-dual',UCD_VAR_STANDARD, p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Y-vel-dual',UCD_VAR_STANDARD, p_Ddata2)
      
      ! Dual pressure
      call lsyssc_getbase_double (rprjVector%rvectorBlock(6),p_Ddata)
      call ucd_addVariableElementBased (rexport,'pressure-dual',UCD_VAR_STANDARD, p_Ddata)
      
      ! If we have a simple Q1~ discretisation, calculate the streamfunction.
      if (rvectorTmp%p_rblockDiscr%RspatialDiscr(1)% &
          ccomplexity .eq. SPDISC_UNIFORM) then
          
        ieltype = rvectorTmp%p_rblockDiscr%RspatialDiscr(1)% &
                  RelementDistr(1)%celement
                  
        if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
            
          call ppns2D_streamfct_uniform (rvectorTmp,rprjVector%rvectorBlock(1))
          
          call lsyssc_getbase_double (rprjVector%rvectorBlock(1),p_Ddata)
          call ucd_addVariableVertexBased (rexport,'streamfunction',&
              UCD_VAR_STANDARD, p_Ddata)
              
        end if
        
      end if
      
      ! Write the file to disc, that's it.
      call ucd_write (rexport)
      call ucd_release (rexport)
      
      ! Release the auxiliary vector
      call lsysbl_releaseVector (rprjVector)
      
      ! Release the discretisation structure
      call spdiscr_releaseBlockDiscr(rprjDiscretisation)
      
      ! Throw away the discrete BC's - not used anymore.
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      call bcasm_releaseDiscreteFBC (rdiscreteFBC)
      
    end do

    ! Release memory.
    call lsysbl_releaseVector (rvectorTmp)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine cc_initpostprocessing (rproblem,rpostprocessing)

!<description>
  ! Initialises the given postprocessing structure rpostprocessing
  ! according to the main problem rproblem. The structure can then be used
  ! to generate postprocessing data.
!</description>

!<input>
  ! A problem structure that describes the main problem to solve.
  type(t_problem), intent(IN),target :: rproblem
!</input>

!<output>
  type(t_c2d2postprocessing), intent(OUT) :: rpostprocessing
!</output>

!</subroutine>

  ! local variables
  type(t_blockDiscretisation), pointer :: p_rdiscr

    ! For postprocessing, we need discretisation structures in the Q0 and Q1 space,
    ! later perhaps in the Q2 space. For this purpose, derive the corresponding
    ! discretisation structure using the 'main' discretisation structure on the
    ! maximum level.
    !
    ! For simplicity, we use only the discretisation structure of the X-velocity
    ! to derive everything.
    
    p_rdiscr => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation

    ! Piecewise constant space:
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscr(1), &
                 EL_Q0, CUB_G1X1, &
                 rpostprocessing%rdiscrConstant)

    ! Piecewise linear space:
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rpostprocessing%rdiscrLinear)
  
    ! Piecewise quadratic space:
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscr(1), &
                 EL_Q2, CUB_G3X3, &
                 rpostprocessing%rdiscrQuadratic)
  
  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine cc_clearpostprocessing (rpostprocessing)

!<description>
  ! Releases all calculated data in the given postprocessing structure
  ! so that it can be allocated again in the calculation routines.
  ! This routine must be called at the end of the postprocessing routines
  ! to release the temporary memory that was allocated for the vectors
  ! in the postprocessing structure.
!</description>

!<inputoutput>
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors which might be allocated.
    call lsyssc_releaseVector (rpostprocessing%rvectorVelX)
    call lsyssc_releaseVector (rpostprocessing%rvectorVelY)
    call lsyssc_releaseVector (rpostprocessing%rvectorPressure)
    call lsyssc_releaseVector (rpostprocessing%rvectorPressureCells)
    call lsyssc_releaseVector (rpostprocessing%rvectorStreamfunction)
    call lsyssc_releaseVector (rpostprocessing%rvectorH1err)
    call lsyssc_releaseVector (rpostprocessing%rvectorH1errCells)

  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine cc_donepostprocessing (rpostprocessing)

!<description>
  ! Releases a given problem structure. All allocated memory of this structure
  ! is released.
!</description>

!<inputoutput>
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors -- if there are still some allocated
    call cc_clearpostprocessing (rpostprocessing)

    ! Release the discretisation structures allocated above.
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrQuadratic)
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrLinear)
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrConstant)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_printControlFunctionalStat (rproblem,dtime,rvector)
  
!<description>
  ! Calculates and prints the value of the optimal control functional J(y,u)
  ! in the stationary case.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
  
  ! Current simulation time.
  real(dp), intent(in) :: dtime
!</input>

!</subroutine>
    
    ! local variables
    real(DP), dimension(3) :: Derror
    real(DP) :: dalphaC

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

    ! Analyse the deviation from the target velocity field.
    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',dalphaC,0.1_DP)
    call cc_optc_stationaryFunctional (rvector,dalphaC,Derror,&
        rproblem%rcollection)

    call output_line ('||y-z||_L2: '//trim(sys_sdEL(Derror(1),2)))
    call output_line ('||u||_L2  : '//trim(sys_sdEL(Derror(2),2)))
    call output_line ('J(y,u)    : '//trim(sys_sdEL(Derror(3),2)))
    
    ! Release assembly-stuff from the collection
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
    
  end subroutine

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
!  integer, DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialU, IdofsTrialP
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
!  integer, DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialU, IdofsTrialP
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
!  ! velocity and pressure
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
!    p_relemDistrU => p_rdiscrU%RelementDistr(RelementDistr)
!    p_relemDistrP => p_rdiscrP%RelementDistr(RelementDistr)
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

end module
