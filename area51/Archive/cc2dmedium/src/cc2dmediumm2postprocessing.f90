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
!# 1.) c2d2_postprocessing
!#     -> Evaluate the solution of the stationary solver, write GMV-file.
!# </purpose>
!##############################################################################

module cc2dmediumm2postprocessing

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
  
  use collection
  use convection
  
  use ucd
  
  use pprocnavierstokes
  use pprocerror
  
  use cc2dmediumm2basic
  use cc2dmedium_callback
  
  implicit none
  
!<types>

!<typeblock>

  ! Postprocessing structure. Whenever the postprocessing calculation routines are
  ! called, they fill this structure using the current solution vector (and other
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
    
    ! Whether nonstationary postprocessing should be used or not.
    logical              :: bnonstationaryPostprocessing
    
    ! Point in time when the next UCD file is to be written out
    real(DP)             :: dnextTimeUCD = 0.0_DP
    
    ! Next file extension for UCD output file.
    integer              :: inextFileSuffixUCD = 0

    ! Point in time when the next Film file is to be written out
    real(DP)             :: dnextTimeFilm = 0.0_DP
    
    ! Next file extension for Film output file.
    integer              :: inextFileSuffixFilm = 0

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

  subroutine c2d2_postprocessingStationary (rproblem,rvector,rpostprocessing)
  
!<description>
  ! Postprocessing of solutions of stationary simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! Postprocessing structure.
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!</subroutine>

    ! Calculate body forces.
    call c2d2_calculateBodyForces (rvector,rproblem)
    
    ! Calculate the divergence
    call c2d2_calculateDivergence (rvector,rproblem)
    
    ! Error analysis, comparison to reference function.
    call c2d2_errorAnalysis (rvector,rproblem)
    
    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    call c2d2_writeUCD (rpostprocessing, rvector, rproblem)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_postprocessingNonstat (rproblem,rvector,rpostprocessing)
  
!<description>
  ! Postprocessing of solutions of stationary simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! Postprocessing structure. Defines what to do with solution vectors.
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!</subroutine>

    ! Calculate body forces.
    call c2d2_calculateBodyForces (rvector,rproblem)
    
    ! Calculate the divergence
    call c2d2_calculateDivergence (rvector,rproblem)

    ! Error analysis, comparison to reference function.
    call c2d2_errorAnalysis (rvector,rproblem)
    
    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    call c2d2_writeUCD (rpostprocessing, rvector, rproblem, &
        rproblem%rtimedependence%dtime)
    
    ! Write film output (raw data vectors)
    call c2d2_writeFilm (rpostprocessing, rvector, rproblem, &
        rproblem%rtimedependence%dtime)
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine c2d2_errorAnalysis (rsolution,rproblem)

!<description>
  ! Performs error analysis on a given solution rsolution as specified
  ! in the .DAT file.
  ! The result of the error analysis is written to the standard output.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(IN) :: rsolution
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>
    
    ! local variables
    real(DP),dimension(3) :: Derr
    real(DP) :: derrorVel, derrorP
    integer :: icalcL2,icalcH1
    
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                     'IERRORANALYSISL2', icalcL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                     'IERRORANALYSISH1', icalcH1, 0)
    
    if ((icalcL2 .ne. 0) .or. (icalcH1 .ne. 0)) then
      call output_lbrk()
      call output_line ('Error Analysis')
      call output_line ('--------------')
    end if
    
    if (icalcL2 .ne. 0) then
    
      call c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u-z||_{L^2}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      derrorVel = (0.5_DP*(Derr(1)**2+Derr(2)**2))

      call pperr_scalar (rsolution%RvectorBlock(3),PPERR_L2ERROR,Derr(3),&
                         ffunction_TargetP,rproblem%rcollection)

      derrorP = Derr(3)
      
      call output_line ('||u-reference||_L2 = '//trim(sys_sdEP(derrorVel,15,6)) )
      call output_line ('||p-reference||_L2 = '//trim(sys_sdEP(derrorP,15,6)) )
      
      call c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end if

    if (icalcL2 .ne. 0) then
    
      call c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u-z||_{L^2}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_H1ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_H1ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      derrorVel = (0.5_DP*(Derr(1)**2+Derr(2)**2))

      call output_line ('||u-reference||_H1 = '//trim(sys_sdEP(derrorVel,15,6)) )
      
      call c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end if
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine c2d2_calculateBodyForces (rsolution,rproblem)

!<description>
  ! Calculates body forces as configured in the .DAT file.
  ! The result is written to the standard output.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(IN) :: rsolution
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>
    
    ! local variables
    
    ! Forces on the object
    real(DP), dimension(NDIM2D) :: Dforces
    real(DP) :: df1,df2
    type(t_boundaryRegion) :: rregion
    
    ! If we have a uniform discretisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    if ((rsolution%p_rblockDiscretisation%RspatialDiscr(1)% &
         ccomplexity .eq. SPDISC_UNIFORM) .and. &
        (boundary_igetNBoundComp(rproblem%p_rboundary) .ge. 2)) then

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      call boundary_createRegion (rproblem%p_rboundary, &
          2, 0, rregion)
      rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      df1 = 1.0_DP/1000.0_DP
      df2 = 0.1_DP * 0.2_DP**2
      call ppns2D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G1_1D,df1,df2)

      call output_lbrk()
      call output_line ('Body forces')
      call output_line ('-----------')
      call output_line ('Body forces real bd., bdc/horiz/vert')
      call output_line (' 2 / ' &
          //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //trim(sys_sdEP(Dforces(2),15,6)) )
      
    endif
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine c2d2_calculateDivergence (rsolution,rproblem)

!<description>
  ! Calculates the divergence of a solution.
  ! The result is written to the standard output.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(IN) :: rsolution
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ieltype
    type(t_matrixScalar) :: rBmatrix
    type(t_vectorScalar), target :: rtempVector
    
    if (rsolution%p_rblockDiscretisation%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rsolution%p_rblockDiscretisation%RspatialDiscr(1)% &
                RelementDistr(1)%itrialElement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
      
        ! Create a temporary vector
        call lsyssc_createVecByDiscr (rsolution%RvectorBlock(3)%p_rspatialDiscretisation,&
            rtempVector,.true.)

        ! Calculate divergence = B1^T u1 + B2^T u2
        call lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB1,&
            rBmatrix,LSYSSC_TR_VIRTUAL)
        call lsyssc_scalarMatVec (&
            rBmatrix, rsolution%RvectorBlock(1), &
            rtempVector, 1.0_DP, 0.0_DP)
        call lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB2,&
            rBmatrix,LSYSSC_TR_VIRTUAL)
        call lsyssc_scalarMatVec (&
            rBmatrix, rsolution%RvectorBlock(2), &
            rtempVector, 1.0_DP, 1.0_DP)
        
        call output_lbrk()
        call output_line ('Divergence')
        call output_line ('----------')
        call output_line ('Divergence = ' &
            //trim(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)) )
            
        call lsyssc_releaseVector (rtempVector)
      
      end if
      
    end if
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine c2d2_writeUCD (rpostprocessing,rvector,rproblem,dtime)

!<description>
  ! Writes an UCD postprocessing file as configured in the DAT file.
  ! (-> GMV, AVS, Paraview,...)
!</description>
  
!<input>
  ! Solution vector.
  type(t_vectorBlock), intent(IN) :: rvector
  
  ! OPTIONAL: Simulation time.
  ! Must be ommitted in stationary simulations.
  real(DP), intent(IN), optional :: dtime
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! Postprocessing structure. Must have been initialised prior
  ! to calling this routine.
  ! The time stamp of the last written out GMV is updated.
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! local variables

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
    type(t_discreteBC), pointer :: p_rdiscreteBC
    type(t_discreteFBC), pointer :: p_rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    real(DP) :: dminTime, dmaxTime, dtimeDifferenceUCD
    integer :: ioutputUCD,ieltype,ilevelUCD
    
    character(SYS_STRLEN) :: sfile,sfilename
    
    if (present(dtime)) then
      ! In a nonstationary simulation, first check if we are allowed
      ! to write something.

      call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                      'DMINTIMEUCD', dminTime, -1.E100_DP)
      call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                      'DMAXTIMEUCD', dmaxTime, 1.E100_DP)
      call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                      'DTIMEDIFFERENCEUCD', dtimeDifferenceUCD, 0.0_DP)
                                      
      if ((dtime .lt. dminTime) .or. (dtime .gt. dmaxTime)) return
      
      if (dtimeDifferenceUCD .gt. 0.0_DP) then
        if (rpostprocessing%dnextTimeUCD .gt.  dtime) return
      else if (dtimeDifferenceUCD .lt. 0.0_DP) then
        if (rpostprocessing%dnextTimeUCD .lt. dtime) return
      end if
      ! Otherwise: Always write!
      
    end if

    ! Type of output:
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'IOUTPUTUCD', ioutputUCD, 0)
    if (ioutputUCD .eq. 0) return

    ! Level of output:
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'ILEVELUCD', ilevelUCD, 0)
    if (ilevelUCD .le. 0) then
      ilevelUCD = rproblem%NLMAX+ilevelUCD
    end if
    
    ilevelUCD = min(rproblem%NLMAX,max(rproblem%NLMIN,ilevelUCD))
    
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
    
    call spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscretisation,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscretisation%RspatialDiscr(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(2))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Discretise the boundary conditions according to the Q1/Q1/Q0
    ! discretisation for implementing them into a solution vector.
    nullify(p_rdiscreteBC)
    call bcasm_discretiseBC (rprjDiscretisation,p_rdiscreteBC, &
                            .false.,getBoundaryValues,rproblem%rcollection,&
                            BCASM_DISCFORSOL)
                            
    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => p_rdiscreteBC
    
    ! The same way, discretise boundary conditions of fictitious boundary components.
    nullify(p_rdiscreteFBC)
    call bcasm_discretiseFBC (rprjDiscretisation,p_rdiscreteFBC, &
                              .false.,getBoundaryValuesFBC,rproblem%rcollection,&
                              BCASM_DISCFORSOL)
    rprjVector%p_rdiscreteBCfict => p_rdiscreteFBC
    
    ! Filter the solution vector to implement discrete BC's.
    call vecfil_discreteBCsol (rprjVector)

    ! Filter the solution vector to implement discrete BC's for fictitious
    ! boundary components.
    call vecfil_discreteFBCsol (rprjVector)
    
    ! Clean up the collection (as we are done with the assembly, that's it.
    call c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

    ! Basic filename
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                 'SFILENAMEUCD', sfile, '')
                                 
    ! Remove possible ''-characters
    read(sfile,*) sfilename
    
    ! Create the actual filename
    sfile = trim(adjustl(sfilename))//'.'//sys_si0(rpostprocessing%inextFileSuffixUCD,5)
                                 
    ! Now we have a Q1/Q1/Q0 solution in rprjVector -- on the level NLMAX.
    ! The next step is to project it down to level ilevelUCD.
    ! Due to the fact that solutions are usually 2-level-ordered,
    ! this can be shortened by taking only the first NVT vertices
    ! of the solution vector!
    
    ! From the attached discretisation, get the underlying triangulation
    ! of that level
    p_rtriangulation => rproblem%RlevelInfo(ilevelUCD)%rtriangulation
    
    ! Start UCD export to GMV file:
    call output_lbrk ()
    call output_line ('Writing GMV file: '//sfile)
    
    select case (ioutputUCD)
    case (1)
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)

    case (2)
      call ucd_startAVS (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
          
    case (3)
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
          
    case DEFAULT
      call output_line ('Invalid UCD ooutput type.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
      stop
    end select
        
    ! Is there a simulation time?
    if (present(dtime)) &
      call ucd_setSimulationTime (rexport,dtime)
    
    ! Write the configuration of the application as comment block
    ! to the output file.
    call ucd_addCommentLine (rexport,'Configuration:')
    call ucd_addCommentLine (rexport,'---------------')
    call ucd_addParameterList (rexport,rproblem%rparamList)
    call ucd_addCommentLine (rexport,'---------------')

    ! Write velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    ! CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, &
    !     p_Ddata(1:p_rtriangulation%NVT))
    ! CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, &
    !     p_Ddata2(1:p_rtriangulation%NVT))
    call ucd_addVarVertBasedVec (rexport,'velocity',&
        p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
    
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, &
        p_Ddata(1:p_rtriangulation%NEL))
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    if (rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
                RelementDistr(1)%itrialElement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
          
        call ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
        
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'streamfunction',&
            UCD_VAR_STANDARD, p_Ddata(1:p_rtriangulation%NVT))
            
      end if
      
    end if
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
    
    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    call bcasm_releaseDiscreteBC (p_rdiscreteBC)
    call bcasm_releaseDiscreteFBC (p_rdiscreteFBC)
    
    if (present(dtime)) then
      ! Update time stamp of last written out GMV.
      rpostprocessing%dnextTimeUCD = rpostprocessing%dnextTimeUCD+dtimeDifferenceUCD
      rpostprocessing%inextFileSuffixUCD = rpostprocessing%inextFileSuffixUCD + 1
    end if

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine c2d2_writeFilm (rpostprocessing,rvector,rproblem,dtime)

!<description>
  ! Writes Film output (raw data vectors) to a file as configured in the
  ! DAT file.
  !
  ! Note: This file is usually only used in a nonstationary simulation.
  ! In a stationary simulation, Film output makes no sense!
!</description>
  
!<input>
  ! Solution vector.
  type(t_vectorBlock), intent(IN) :: rvector
  
  ! Simulation time.
  real(DP), intent(IN) :: dtime
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! Postprocessing structure. Must have been initialised prior
  ! to calling this routine.
  ! The time stamp of the last written out Film file is updated.
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: dminTime, dmaxTime, dtimeDifferenceFilm
    integer :: ioutputFilm,ilevelFilm
    
    type(t_vectorBlock) :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sfile,sfilename
    integer :: ilev
    integer(PREC_VECIDX) :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection
    logical :: bformatted
    
    ! Type of output:
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'IOUTPUTFILM', ioutputFilm, 0)
    if (ioutputFilm .eq. 0) return

    ! First check if we are allowed to write something.
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                    'DMINTIMEFILM', dminTime, -1.E100_DP)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                    'DMAXTIMEFILM', dmaxTime, 1.E100_DP)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                    'DTIMEDIFFERENCEFILM', dtimeDifferenceFilm, 0.0_DP)
                                    
    if ((dtime .lt. dminTime) .or. (dtime .gt. dmaxTime)) return
    
    if (dtimeDifferenceFilm .gt. 0.0_DP) then
      if (rpostprocessing%dnextTimeFilm .gt.  dtime) return
    else if (dtimeDifferenceFilm .lt. 0.0_DP) then
      if (rpostprocessing%dnextTimeFilm .lt. dtime) return
    end if

    ! Basic filename
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                 'SFILENAMEFILM', sfile, '')
                                 
    ! Remove possible ''-characters
    read(sfile,*) sfilename
    
    ! Create the actual filename
    sfile = trim(adjustl(sfilename))//'.'//sys_si0(rpostprocessing%inextFileSuffixFilm,5)
                                 
    ! Level of output:
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'ILEVELFILM', ilevelFilm, 0)
    if (ilevelFilm .le. 0) then
      ilevelFilm = rproblem%NLMAX+ilevelFilm
    end if
    
    ilevelFilm = min(rproblem%NLMAX,max(rproblem%NLMIN,ilevelFilm))
    
    if (ilevelFilm .lt. rproblem%NLMIN) then
      print *,'Warning: Level for solution vector is < NLMIN! &
              &Writing out at level NLMIN!'
      ilevelFilm = rproblem%NLMIN
    end if
    
    ! Write formatted output?
    bformatted = ioutputFilm .ne. 2

    ! Interpolate the solution down to level istart.
    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    do ilev = rproblem%NLMAX,ilevelFilm+1,-1
      
      ! Initialise a vector for the lower level and a prolongation structure.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev-1)%p_rdiscretisation,rvector2,.false.)
      
      call mlprj_initProjectionVec (rprojection,rvector2)
      
      ! Interpolate to the next higher level.
      ! (Don't 'restrict'! Restriction would be for the dual space = RHS vectors!)

      NEQ = mlprj_getTempMemoryVec (rprojection,rvector2,rvector1)
      if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
      call mlprj_performInterpolation (rprojection,rvector2,rvector1, &
                                       rvectorTemp)
      if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)
      
      ! Swap rvector1 and rvector2. Release the fine grid vector.
      call lsysbl_swapVectors (rvector1,rvector2)
      call lsysbl_releaseVector (rvector2)
      
      call mlprj_doneProjection (rprojection)
      
    end do

    call output_lbrk ()
    call output_line ('Writing Film file: '//sfile)

    ! Write out the solution.
    if (bformatted) then
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,&
         0, sfile, '(E22.15)')
    else
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,0, sfile)
    end if

    ! Release temp memory.
    call lsysbl_releaseVector (rvector1)

    ! Update time stamp of last written out Film file.
    rpostprocessing%dnextTimeFilm = rpostprocessing%dnextTimeFilm+dtimeDifferenceFilm
    rpostprocessing%inextFileSuffixFilm = rpostprocessing%inextFileSuffixFilm + 1

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine c2d2_initpostprocessing (rproblem,rpostprocessing)

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
  ! Postprocessing structure.
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
    
    p_rdiscr => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscretisation

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
  
    ! Initialise the time/file suffix when the first UCD file is to be written out.
    rpostprocessing%bnonstationaryPostprocessing = (rproblem%itimedependence .ne. 0)
    if (rproblem%itimedependence .ne. 0) then
      rpostprocessing%dnextTimeUCD = rproblem%rtimedependence%dtimeInit
      call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
         'ISTARTSUFFIXUCD', rpostprocessing%inextFileSuffixUCD, 1)
      
      rpostprocessing%dnextTimeFilm = rproblem%rtimedependence%dtimeInit
      call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
         'ISTARTSUFFIXFILM', rpostprocessing%inextFileSuffixFilm, 1)
    end if
                                    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine c2d2_copyPostprocessing (rpostprocessingSrc,rpostprocessingDst)

!<description>
  ! Copies the state of the postprocessing from rpostprocessingSrc
  ! to rpostprocessingDst.
  ! This is typically called to back up or to restore the current state
  ! of a running postprocessing structure.
!</description>

!<input>
  ! Source Postprocessing structure.
  type(t_c2d2postprocessing), intent(IN) :: rpostprocessingSrc
!</input>

!<inputoutput>
  ! Destination Postprocessing structure.
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessingDst
!</inputoutput>

!</subroutine>

    ! Initialise the time/file suffix when the first UCD file is to be written out.
    if (rpostprocessingSrc%bnonstationaryPostprocessing) then
      rpostprocessingDst%bnonstationaryPostprocessing = &
          rpostprocessingSrc%bnonstationaryPostprocessing

      rpostprocessingDst%dnextTimeUCD = rpostprocessingSrc%dnextTimeUCD
      rpostprocessingDst%inextFileSuffixUCD = rpostprocessingSrc%inextFileSuffixUCD

      rpostprocessingDst%dnextTimeFilm = rpostprocessingSrc%dnextTimeFilm
      rpostprocessingDst%inextFileSuffixFilm = rpostprocessingSrc%inextFileSuffixFilm
    end if
                                    
  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine c2d2_clearpostprocessing (rpostprocessing)

!<description>
  ! Releases all calculated data in the given postprocessing structure
  ! so that it can be allocated again in the calculation routines.
  ! This routine must be called at the end of the postprocessing routines
  ! to release the temporary memory that was allocated for the vectors
  ! in the postprocessing structure.
!</description>

!<inputoutput>
  ! Postprocessing structure.
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors which might be allocated.
    if (rpostprocessing%rvectorVelX%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorVelX)
    if (rpostprocessing%rvectorVelY%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorVelY)
    if (rpostprocessing%rvectorPressure%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorPressure)
    if (rpostprocessing%rvectorPressureCells%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorPressureCells)
    if (rpostprocessing%rvectorStreamfunction%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorStreamfunction)
    if (rpostprocessing%rvectorH1err%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorH1err)
    if (rpostprocessing%rvectorH1errCells%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorH1errCells)

  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine c2d2_donepostprocessing (rpostprocessing)

!<description>
  ! Releases a given problem structure. All allocated memory of this structure
  ! is released.
!</description>

!<inputoutput>
  type(t_c2d2postprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors -- if there are still some allocated
    call c2d2_clearpostprocessing (rpostprocessing)

    ! Release the discretisation structures allocated above.
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrQuadratic)
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrLinear)
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrConstant)
  
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
