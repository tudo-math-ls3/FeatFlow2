!##############################################################################
!# ****************************************************************************
!# <name> postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains various routines for the postprocessing during a
!# space-time optimisation.
!#
!# The following routines can be found here:
!#
!# 1.) optcpp_initpostprocessing
!#     -> Initialise the postprocessing
!#
!# 2.) optcpp_donepostprocessing
!#     -> Clean up the postprocessing
!#
!# 3.) optcpp_postprocessSingleSol
!#     -> Postprocessing of a solution at a specific timestep
!#
!# 4.) optcpp_postprocessSpaceTimeVec
!#     -> Postprocessing of a space-time vector
!#
!# 5.) optcpp_postprocSpaceVisOutput
!#     -> Write out visualisation files of all solutions in a space-time vector
!# </purpose>
!##############################################################################

module postprocessing

  use fsystem
  use storage
  use genoutput
  use io
  use basicgeometry
  use linearsolver
  use boundary
  use element
  use cubature
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use linearformevaluation
  use matrixfilters
  use vectorfilters
  use discretebc
  use discretefbc
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
  use vectorio
  
  use pprocnavierstokes

  use analyticsolution
  use optcanalysis
  use spatialbcdef
  use spacetimevectors
  use spacematvecassembly
  
  use constantsoptc
  use structuresoptc
  use timediscretisation
  use structuresoptflow
  use user_callback
  
  !use spacetimediscretisation
  
  implicit none
  
  private
  
  public :: t_optcPostprocessing  
  public :: optcpp_initpostprocessing
  public :: optcpp_donepostprocessing
  public :: optcpp_postprocessSingleSol
  public :: optcpp_postprocessSpaceTimeVec
  public :: optcpp_postprocSpaceVisOutput
  
!!<constants>
!
!!<constantblock description="Type identifiers for creating a space-time vector.">
!
!  ! A zero space-time vector.
!  integer, parameter :: CCSTV_ZERO = 0
!
!  ! The space-time vector is created as copy of a single vector.
!  integer, parameter :: CCSTV_STATIONARY = 1
!  
!  ! The space-time vector is created by reading in files.
!  integer, parameter :: CCSTV_READFILES = 2
!  
!  ! The space-time vector is created as forward simulation starting from an
!  ! initial solution.
!  integer, parameter :: CCSTV_FORWARDSIMULATION = 3
!
!!</constantblock>
!
!!</constants>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_initpostprocessing (rpostproc,cspace,rboundaryConditions,&
      rtimeDiscr,rspaceDiscr)
  
!<description>
  ! Initialises the postprocessing
!</description>

!<input>
  ! Space that is available in rsolution. One of the CCSPACE_xxxx constants.
  integer, intent(in) :: cspace
  
  ! Boundary conditions to use.
  type(t_optcBDC), intent(in), target  :: rboundaryConditions

  ! Underlying space discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Underlying time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
!</input>
  
!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Fetch data
    rpostproc%cspace = cspace
    rpostproc%p_rspaceDiscr => rspaceDiscr
    rpostproc%p_rtimeDiscr => rtimeDiscr
    rpostproc%p_rboundaryConditions => rboundaryConditions

    ! Initialise a P1/Q1/P0/Q0 structure for visualisationoutput.
    call spdiscr_duplicateBlockDiscr(rpostproc%p_rspaceDiscr,rpostproc%rspaceDiscrLinear)
    
    call spdiscr_deriveDiscr_triquad (&
                 rpostproc%p_rspaceDiscr%RspatialDiscr(1), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rpostproc%rspaceDiscrLinear%RspatialDiscr(1))

    call spdiscr_deriveDiscr_triquad (&
                 rpostproc%p_rspaceDiscr%RspatialDiscr(2), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rpostproc%rspaceDiscrLinear%RspatialDiscr(2))

    call spdiscr_deriveDiscr_triquad (&
                 rpostproc%p_rspaceDiscr%RspatialDiscr(3), &
                 EL_P0, EL_Q0, CUB_TRZ_T, CUB_G2X2, &
                 rpostproc%rspaceDiscrLinear%RspatialDiscr(3))

    call spdiscr_deriveDiscr_triquad (&
                 rpostproc%p_rspaceDiscr%RspatialDiscr(4), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rpostproc%rspaceDiscrLinear%RspatialDiscr(4))

    call spdiscr_deriveDiscr_triquad (&
                 rpostproc%p_rspaceDiscr%RspatialDiscr(5), &
                 EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                 rpostproc%rspaceDiscrLinear%RspatialDiscr(5))

    call spdiscr_deriveDiscr_triquad (&
                 rpostproc%p_rspaceDiscr%RspatialDiscr(6), &
                 EL_P0, EL_Q0, CUB_TRZ_T, CUB_G2X2, &
                 rpostproc%rspaceDiscrLinear%RspatialDiscr(6))

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_donepostprocessing (rpostproc)
  
!<description>
  ! Clean up the postprocessing structure
!</description>
  
!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Release the discretisations structure
    call spdiscr_releaseBlockDiscr (rpostproc%rspaceDiscrLinear)

    ! Clean up other data.
    rpostproc%cspace = CCSPACE_PRIMAL
    nullify(rpostproc%p_rspaceDiscr)
    nullify(rpostproc%p_rtimeDiscr)
    nullify(rpostproc%p_rboundaryConditions)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocessSingleSol (rpostproc,ifileid,dtime,rvector,&
      roptcontrol,rsettings,bfirstFile)
  
!<description>
  ! Postprocessing of a single solution at a definite time.
  ! Writes the solution into a visualisation file, calculates forces,...
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! The global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
  
  ! Id of the file on the hard disc; added to the filename.
  integer, intent(in) :: ifileid
  
  ! Current simulation time.
  real(dp), intent(in) :: dtime
  
  ! Parameters about the optimal control
  type(t_settings_optcontrol), intent(in) :: roptControl
  
  ! Must be set to TRUE for the first file of a sequence of files or
  ! of teh file does not belong to a sequence.
  logical, intent(in) :: bfirstFile
!</input>

!</subroutine>

  ! local variables
  integer :: ieltype,cflag,iunit
  logical :: bfileexists
  type(t_collection) :: rcollection

  ! We need some more variables for postprocessing - i.e. writing
  ! a GMV file.
  real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

  ! A pointer to the triangulation.
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A vector accepting Q1 data
  type(t_vectorBlock) :: rprjVector
  
  ! Output block for UCD output to GMV file
  type(t_ucdExport) :: rexport
  
  ! Forces on the object
  real(DP), dimension(NDIM2D) :: Dforces
  type(t_boundaryRegion) :: rregion
  
  ! Divergence
  !type(t_vectorScalar), target :: rtempVector
  
  character(SYS_STRLEN) :: sfilename

  ! Discrete boundary conditions
  type(t_discreteBC) :: rdiscreteBC
  
  ! Discrete fictitious boundary conditions
  type(t_discreteFBC) :: rdiscreteFBC

    ! -------------------------------------------------------------------------
    ! Body forces
    ! -------------------------------------------------------------------------

    ! When writing to a file is enabled, delete the file in the first timestep.
    cflag = SYS_APPEND
    if (bfirstFile) cflag = SYS_REPLACE

    ! If we have a uniform discretisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    if ((rpostproc%p_rspaceDiscr%RspatialDiscr(1)% &
         ccomplexity .eq. SPDISC_UNIFORM) .and. &
        (boundary_igetNBoundComp(rvector%p_rblockDiscr%p_rboundary) .ge. &
            rpostproc%ibodyForcesBdComponent) .and.&
        (rpostproc%icalcForces .ne. 0)) then

      call output_lbrk()
      call output_line ('Body forces real bd., bdc/horiz/vert')

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      call boundary_createRegion (rvector%p_rblockDiscr%p_rboundary, &
          2, 0, rregion)
      rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call ppns2D_bdforces_uniform (rvector,rregion,Dforces,CUB_G1_1D,&
          rpostproc%dbdForcesCoeff1,rpostproc%dbdForcesCoeff2)
      
      call output_line (' 2 / ' &
          //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //trim(sys_sdEP(Dforces(2),15,6)) )
      
      if (rpostproc%iwriteBodyForces .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(rpostproc%sfilenameBodyForces, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,'(A)') '# timestep time bdc horiz vert'
        end if
        write (iunit,'(A)') trim(sys_siL(ifileid,10)) // ' ' &
            // trim(sys_sdEL(dtime,10)) // ' ' &
            // trim(sys_siL(rpostproc%ibodyForcesBdComponent,10)) // ' ' &
            // trim(sys_sdEL(Dforces(1),10)) // ' '&
            // trim(sys_sdEL(Dforces(2),10))
        close (iunit)
      end if
      
    endif
    
    ! -------------------------------------------------------------------------
    ! Writing out of the final solution
    ! -------------------------------------------------------------------------
    
    if (rpostproc%sfinalSolutionFileName .ne. "") then
      ! Write the current solution to disc as it is.
      sfilename = trim(rpostproc%sfinalSolutionFileName)//'.'//sys_si0(ifileid,5)
      
      call output_lbrk ()
      call output_line ('Writing solution file: '//trim(sfilename))
      
      call vecio_writeBlockVectorHR (rvector, "vector"//sys_si0(ifileid,5), .false.,&
          0, sfilename,  "(E20.10)")
    end if

    ! -------------------------------------------------------------------------
    ! Visualisation output
    ! -------------------------------------------------------------------------

    if (rpostproc%ioutputUCD .ne. 0) then

      ! If we have a simple Q1~ discretisation, calculate the streamfunction.
!      IF (rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
!          ccomplexity .EQ. SPDISC_UNIFORM) THEN
!          
!        ieltype = rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
!                  RelementDistr(1)%itrialElement
!                  
!        IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q1T) THEN
!        
!          ! Create a temporary vector 
!          CALL lsyssc_createVecByDiscr (rvector%RvectorBlock(3)%p_rspatialDiscretisation,&
!              rtempVector,.TRUE.)
!
!          ! Calculate divergence = B1^T u1 + B2^T u2
!          CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB1,&
!              rBmatrix,LSYSSC_TR_VIRTUAL)
!          CALL lsyssc_scalarMatVec (&
!              rBmatrix, rvector%RvectorBlock(1), &
!              rtempVector, 1.0_DP, 0.0_DP)
!          CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB2,&
!              rBmatrix,LSYSSC_TR_VIRTUAL)
!          CALL lsyssc_scalarMatVec (&
!              rBmatrix, rvector%RvectorBlock(2), &
!              rtempVector, 1.0_DP, 1.0_DP)
!          
!          CALL output_lbrk()
!          CALL output_line ('Divergence = ' &
!              //TRIM(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)) )
!              
!          CALL lsyssc_releaseVector (rtempVector)
!        
!        END IF
!        
!      END IF    

      ! Initialise boundary condition structures    
      call bcasm_initDiscreteBC(rdiscreteBC)
      call bcasm_initDiscreteFBC(rdiscreteFBC)

      ! The pressure discretisation substructure stays the old.
      !
      ! Now set up a new solution vector based on this discretisation,
      ! allocate memory.
      call lsysbl_createVecBlockByDiscr (rpostproc%rspaceDiscrLinear,rprjVector,.false.)
      
      ! Then take our original solution vector and convert it according to the
      ! new discretisation:
      call spdp_projectSolution (rvector,rprjVector)
      
      call collct_init(rcollection)
      call user_initCollectForAssembly(rsettings%rglobalData,dtime,rcollection)
      
      call sbc_assembleBDconditions (rsettings%roptcBDC,dtime,rpostproc%rspaceDiscrLinear,&
          rpostproc%p_rtimeDiscr,CCSPACE_PRIMALDUAL,rdiscreteBC,rsettings%rglobalData)
      call sbc_assembleFBDconditions (dtime,rpostproc%rspaceDiscrLinear,rpostproc%p_rtimeDiscr,&
          CCSPACE_PRIMALDUAL,rdiscreteFBC,rsettings%rglobalData)
      
      call user_doneCollectForAssembly(rsettings%rglobalData,rcollection)
      call collct_done(rcollection)
      
      ! Filter the solution vector to implement discrete BC's.
      call vecfil_discreteBCsol (rprjVector,rdiscreteBC)

      ! Filter the solution vector to implement discrete BC's for fictitious 
      ! boundary components.
      call vecfil_discreteFBCsol (rprjVector,rdiscreteFBC)
      
      ! Release boundary condition structures.
      call bcasm_releaseDiscreteBC(rdiscreteBC)
      call bcasm_releaseDiscreteFBC(rdiscreteFBC)
      
      ! Now we have a Q1/Q1/Q0 solution in rprjVector.
      !
      ! From the attached discretisation, get the underlying triangulation
      p_rtriangulation => rpostproc%p_rspaceDiscr%p_rtriangulation
      
      ! Check if we have a filename where to write GMV output to.
      if (rpostproc%sfilenameUCD .ne. "") then
      
        ! Start UCD export:
        sfilename = trim(rpostproc%sfilenameUCD)//'.'//sys_si0(ifileid,5)
        
        call output_lbrk ()
        call output_line ('Writing visualisation file: '//trim(sfilename))
        
        select case (rpostproc%ioutputUCD)
        case (1)
          call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)

        case (2)
          call ucd_startAVS (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
              
        case (3)
          call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
              
        case default
          call output_line ('Invalid UCD ooutput type.', &
                            OU_CLASS_ERROR,OU_MODE_STD,'fbsim_writeUCD')
          stop
        end select
        
        ! Write the configuration of the application as comment block
        ! to the output file.
        call ucd_addCommentLine (rexport,'Configuration:')
        call ucd_addCommentLine (rexport,'---------------')
        call ucd_addParameterList (rexport,rsettings%p_rparlist)
        call ucd_addCommentLine (rexport,'---------------')

        ! Write velocity field
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
        
        call ucd_addVarVertBasedVec (rexport,'velocity_p',&
            p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
        
        ! Write out cell based or node based pressure.
        ieltype = rprjVector%p_rblockDiscr%RspatialDiscr(3)% &
                  RelementDistr(1)%celement
                  
        if ((elem_getPrimaryElement(ieltype) .eq. EL_Q1) .or. &
            ((elem_getPrimaryElement(ieltype) .eq. EL_P1))) then
          call ucd_addVariableVertexBased (rexport,'pressure_p',UCD_VAR_STANDARD, &
              p_Ddata(1:p_rtriangulation%NVT))
        else
          call ucd_addVariableElementBased (rexport,'pressure_p',UCD_VAR_STANDARD, &
              p_Ddata(1:p_rtriangulation%NEL))
        end if
        
        ! Dual velocity field
        call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
        call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)

        call ucd_addVarVertBasedVec (rexport,'velocity_d',&
            p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
        
        ! Write out cell based or node based dual pressure.
        call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
        ieltype = rprjVector%p_rblockDiscr%RspatialDiscr(6)% &
                  RelementDistr(1)%celement
                  
        if ((elem_getPrimaryElement(ieltype) .eq. EL_Q1) .or. &
            ((elem_getPrimaryElement(ieltype) .eq. EL_P1))) then
          call ucd_addVariableVertexBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
              p_Ddata(1:p_rtriangulation%NVT))
        else
          call ucd_addVariableElementBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
              p_Ddata(1:p_rtriangulation%NEL))
        end if

        ! Control u = P[min/max](-1/alpha lambda)
        call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
        call lsyssc_scaleVector (rprjVector%RvectorBlock(4),-1.0_DP/roptControl%dalphaC)

        call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
        call lsyssc_scaleVector (rprjVector%RvectorBlock(5),-1.0_DP/roptControl%dalphaC)

        if (roptControl%rconstraints%ccontrolConstraints .ne. 0) then
          call smva_projectControlTimestep (rprjVector%RvectorBlock(4),&
              roptControl%rconstraints%dumin1,roptControl%rconstraints%dumax1)

          call smva_projectControlTimestep (rprjVector%RvectorBlock(5),&
              roptControl%rconstraints%dumin2,roptControl%rconstraints%dumax2)
        end if

        call ucd_addVarVertBasedVec (rexport,'control',&
            p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
        
        ! If we have a simple Q1~ discretisation, calculate the streamfunction.
        if (rvector%p_rblockDiscr%RspatialDiscr(1)%ccomplexity .eq. SPDISC_UNIFORM) then
            
          ieltype = rpostproc%p_rspaceDiscr%RspatialDiscr(1)%RelementDistr(1)%celement
                    
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
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocessSpaceTimeVec (rpostproc,rvector,roptcontrol,rsettings)
  
!<description>
  ! Postprocessing of a space-time vector.
  ! Writes the solution into a visualisation file, calculates forces,...
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! The global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_spaceTimeVector), intent(IN) :: rvector
  
  ! Parameters about the optimal control
  type(t_settings_optcontrol), intent(inout) :: roptControl
!</input>

!</subroutine>

  ! local variables
  type(t_vectorBlock) :: rvecTemp
  integer :: istep
  real(dp) :: dtime
  real(DP), dimension(4) :: Derror

    ! -------------------------------------------------------------------------
    ! Visualisation output, force calculation
    ! -------------------------------------------------------------------------

    ! Create a temp vector in space for postprocessing
    call lsysbl_createVectorBlock (rvector%p_rspaceDiscr,rvecTemp)

    ! Write a file for every timestep
    do istep = 1,rvector%NEQtime
      call tdiscr_getTimestep(rpostproc%p_rtimediscr,istep-1,dtime)
      call sptivec_getTimestepData(rvector,istep,rvecTemp)
      call optcpp_postprocessSingleSol (rpostproc,istep-1,dtime,rvecTemp,&
          roptControl,rsettings,istep .eq. 1)
    end do
    
    call lsysbl_releaseVector (rvecTemp)
    
    ! -------------------------------------------------------------------------
    ! Error analysis
    ! -------------------------------------------------------------------------

    ! If error analysis has to be performed, we can calculate
    ! the real error.
    if (rpostproc%icalcError .eq. 1) then
      call output_lbrk()
      call optcana_analyticalError (rsettings%rglobalData,&
          rsettings%rsettingsOptControl%rconstraints,&
          rvector,rpostproc%ranalyticRefFunction,&
          Derror(1),Derror(2),Derror(3),Derror(4),.true.)
      call output_lbrk()
      call output_line ('||y-y0||           = '//trim(sys_sdEL(Derror(1),10)))   
      call output_line ('||p-p0||           = '//trim(sys_sdEL(Derror(2),10)))   
      call output_line ('||lambda-lambda0|| = '//trim(sys_sdEL(Derror(3),10)))   
      call output_line ('||xi-xi0||         = '//trim(sys_sdEL(Derror(4),10)))   
    end if

    ! Should we calculate the functional?
    if (rpostproc%icalcFunctionalValues .ne. 0) then
      call output_lbrk()
      call optcana_nonstatFunctional (rsettings%rglobalData,roptControl%rconstraints,&
          rvector,roptcontrol%rtargetFlow,roptControl%dalphaC,roptControl%dgammaC,Derror)
      call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(1),10)))
      call output_line ('||u||         = '//trim(sys_sdEL(Derror(2),10)))
      call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
      call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(4),10)))
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocSpaceVisOutput (rsettings,rspaceDiscr,rtimeDiscr,rvector,&
      ioutputUCD,sfilename)
  
!<description>
  ! For every sub-solution in the global space-time vector rvector,
  ! a visualisation file is written to disc.
  !
  ! Used for debug purposes, as this routine does not use a
  ! postprocessing structure.
!</description>

!<input>
  ! Global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! A space-time discretisation structure defining the discretisation of
  ! rvector.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
  
  ! Underlying time discretisation
  type(t_timeDiscretisation), intent(in) :: rtimeDiscr

  ! A space-time vector. For every timestep, a GMV is written.
  type(t_spaceTimeVector), intent(in) :: rvector
  
  ! Type of output file to generate from solutions.
  ! 0=disabled
  ! 1=GMV
  ! 2=AVS
  ! 3=Paraview (VTK)
  ! 4=Matlab
  integer, intent(in) :: ioutputUCD

  ! A path + basic filename for the GMV-files. A number '.00000','.00001',...
  ! is appended for every timestep.
  character(LEN=*), intent(in) :: sfilename
!</input>

!</subroutine>

    ! local variables
    type(t_optcPostprocessing) :: rpostproc
    type(t_vectorBlock) :: rvecTemp
    integer :: istep
    real(dp) :: dtime
    
    ! Create a default postprocessing structure
    call optcpp_initpostprocessing (rpostproc,CCSPACE_PRIMALDUAL,&
        rsettings%roptcBDC,rtimeDiscr,rspaceDiscr)
    
    ! Only apply the visualisation output
    rpostproc%ioutputUCD = ioutputUCD
    rpostproc%sfilenameUCD = sfilename
    
    ! Create a temp vector for the output
    call lsysbl_createVectorBlock(rspaceDiscr,rvecTemp)
    
    ! Write a file for every timestep
    do istep = 1,rvector%NEQtime
      call tdiscr_getTimestep(rtimediscr,istep-1,dtime)
      call sptivec_getTimestepData(rvector,istep,rvecTemp)
      call optcpp_postprocessSingleSol (rpostproc,istep,dtime,rvecTemp,&
          rsettings%rsettingsOptControl,rsettings,istep .eq. 1)
    end do
    
    ! Release all created stuff
    call lsysbl_releaseVector (rvecTemp)
    call optcpp_donepostprocessing (rpostproc)
    
  end subroutine

!  !****************************************************************************
!
!!<subroutine>
!
!  subroutine cc_initpostprocessing (rspaceDiscr,rtimeDiscr,rpostprocessing)
!
!!<description>
!  ! Initialises the given postprocessing structure rpostprocessing
!  ! according to the main problem rproblem. The structure can then be used
!  ! to generate postprocessing data.
!!</description>
!
!!<input>
!  ! A space-time discretisation structure defining the discretisation of
!  ! rvector.
!  type(t_spatialDiscretisation), intent(in), target :: rspaceDiscr
!  
!  ! Underlying time discretisation
!  type(t_timeDiscretisation), intent(in),target :: rtimeDiscr
!
!!</input>
!
!!<output>  
!  type(t_optcPostprocessing), intent(out) :: rpostprocessing
!!</output>
!
!!</subroutine>
!
!    ! Rememnber the discretisation structures.
!    rpostprocessing%p_rspaceDiscr => rspaceDiscr
!    rpostprocessing%p_rtimeDiscr => rtimeDiscr
!
!    ! For postprocessing, we need discretisation structures in the Q0 and Q1 space,
!    ! later perhaps in the Q2 space. For this purpose, derive the corresponding
!    ! discretisation structure using the 'main' discretisation structure on the
!    ! maximum level.
!    !
!    ! For simplicity, we use only the discretisation structure of the X-velocity
!    ! to derive everything.
!    
!    p_rdiscr => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation
!
!    ! Piecewise constant space:
!    call spdiscr_deriveSimpleDiscrSc (&
!                 rspaceDiscr%RspatialDiscr(1), &
!                 EL_Q0, CUB_G1X1, &
!                 rpostprocessing%rdiscrConstant)
!
!    ! Piecewise linear space:
!    call spdiscr_deriveSimpleDiscrSc (&
!                 rspaceDiscr%RspatialDiscr(1), &
!                 EL_Q1, CUB_G2X2, &
!                 rpostprocessing%rdiscrLinear)
!  
!    ! Piecewise quadratic space:
!    call spdiscr_deriveSimpleDiscrSc (&
!                 rspaceDiscr%RspatialDiscr(1), &
!                 EL_Q2, CUB_G3X3, &
!                 rpostprocessing%rdiscrQuadratic)
!  
!  end subroutine
!
!  !****************************************************************************
!  
!!<subroutine>
!
!  subroutine cc_clearpostprocessing (rpostprocessing)
!
!!<description>
!  ! Releases all calculated data in the given postprocessing structure
!  ! so that it can be allocated again in the calculation routines.
!  ! This routine must be called at the end of the postprocessing routines
!  ! to release the temporary memory that was allocated for the vectors
!  ! in the postprocessing structure.
!!</description>
!
!!<inputoutput>  
!  type(t_optcPostprocessing), intent(INOUT) :: rpostprocessing
!!</inputoutput>
!
!!</subroutine>
!
!    ! Release all vectors which might be allocated.
!    call lsyssc_releaseVector (rpostprocessing%rvectorVelX)
!    call lsyssc_releaseVector (rpostprocessing%rvectorVelY)
!    call lsyssc_releaseVector (rpostprocessing%rvectorPressure)
!    call lsyssc_releaseVector (rpostprocessing%rvectorPressureCells)
!    call lsyssc_releaseVector (rpostprocessing%rvectorStreamfunction)
!    call lsyssc_releaseVector (rpostprocessing%rvectorH1err)
!    call lsyssc_releaseVector (rpostprocessing%rvectorH1errCells)
!
!  end subroutine
!
!  !****************************************************************************
!  
!!<subroutine>
!
!  subroutine cc_donepostprocessing (rpostprocessing)
!
!!<description>
!  ! Releases a given problem structure. All allocated memory of this structure
!  ! is released.
!!</description>
!
!!<inputoutput>  
!  type(t_optcPostprocessing), intent(INOUT) :: rpostprocessing
!!</inputoutput>
!
!!</subroutine>
!
!    ! Release all vectors -- if there are still some allocated
!    call cc_clearpostprocessing (rpostprocessing)
!
!    ! Release the discretisation structures allocated above.
!    call spdiscr_releaseDiscr(rpostprocessing%rdiscrQuadratic)
!    call spdiscr_releaseDiscr(rpostprocessing%rdiscrLinear)
!    call spdiscr_releaseDiscr(rpostprocessing%rdiscrConstant)
!  
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_printControlFunctionalStat (rproblem,dtime,rvector)
!  
!!<description>
!  ! Calculates and prints the value of the optimal control functional J(y,u) 
!  ! in the stationary case.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!!</inputoutput>
!
!!<input>
!  ! The solution vector which is to be evaluated by the postprocessing routines.
!  type(t_vectorBlock), intent(IN) :: rvector
!  
!  ! Current simulation time.
!  real(dp), intent(in) :: dtime
!!</input>
!
!!</subroutine>
!    
!    ! local variables
!    real(DP), dimension(3) :: Derror
!    real(DP) :: dalphaC
!
!    ! Initialise the collection for the assembly process with callback routines.
!    ! Basically, this stores the simulation time in the collection if the
!    ! simulation is nonstationary.
!    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)
!
!    ! Analyse the deviation from the target velocity field.
!    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
!                                'dalphaC',dalphaC,0.1_DP)
!    call optcana_stationaryFunctional (rvector,dalphaC,Derror,&
!        rproblem%rcollection)
!
!    call output_line ('||y-z||_L2: '//trim(sys_sdEL(Derror(1),2)))
!    call output_line ('||u||_L2  : '//trim(sys_sdEL(Derror(2),2)))
!    call output_line ('J(y,u)    : '//trim(sys_sdEL(Derror(3),2)))
!    
!    ! Release assembly-stuff from the collection
!    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
!    
!  end subroutine
!
!  !****************************************************************************
!  
!!<subroutine>
!
!  subroutine cc_createVector (rvector,cvectorType,rinitialSol)
!
!!<description>
!  ! Creates a space-time vector based on parameters.
!!</description>
!
!!<input>
!  ! Identification number that defines how to create rvector. A CCSTV_xxxx constant.
!  ! =CCSTV_ZERO: A zero space-time vector.
!  ! =CCSTV_STATIONARY: The space-time vector is created as copy of the single vector rinitialSol.
!  ! =CCSTV_READFILES: The space-time vector is created by reading in files.
!  ! =CCSTV_FORWARDSIMULATION: The space-time vector is created as forward simulation starting from an
!  ! initial solution rinitialSol.
!  integer, intent(in) :: cvectorType
!
!  ! OPTIONAL: Initial solution vector. Only primal solution. Must be at the same
!  ! spatial level and discretised with the same discretisation as the timesteps
!  ! in rvector.
!  type(t_vectorBlock), intent(in), optional :: rinitialSol
!!</input>
!
!!<inputoutput>  
!  ! Space-time vector to be created.
!  type(t_spacetimevector), intent(inout) :: rvector
!  
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: istep
!    type(t_vectorBlock) :: rvectorTemp
!
!    select case (cvectorType)
!    case (CCSTV_ZERO)
!      ! Initialise by zero.
!      call sptivec_clearVector (rvector)
!      
!    case (CCSTV_STATIONARY)
!      ! Copy the initial solution to all timesteps.
!      
!    end select
!
!  end subroutine



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
