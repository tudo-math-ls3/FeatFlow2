!##############################################################################
!# ****************************************************************************
!# <name> ccpostprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains postprocessing routines for the CC2D solver.
!#
!# The following routines can be found here:
!#
!# 1.) cc_initpostprocessing
!#     -> Initialise the postprocessing.
!#
!# 2.) cc_donepostprocessing
!#     -> Clean up the postprocessing.
!#
!# 3.) cc_postprocessingStationary
!#     -> Postprocessing of a solution in a stationary solution process.
!#        Evaluate the solution of the stationary solver, write GMV-file.
!#
!# 4.) cc_postprocessingNonstat
!#     -> Postprocessing of a solution in a nonstationary solution process.
!#        Evaluate the solution of the nonstationary solver, write GMV-file.
!#
!# Auxiliary routines:
!#
!# 1.) cc_errorAnalysis
!#     -> Perform error analysis (comparison to analytic solution;
!#        as configured in the DAT files).
!#
!# 2.) cc_calculateBodyForces
!#     -> Calculate the body forces on an object.
!#
!# 3.) cc_calculateDivergence
!#     -> Calculate the divergence.
!#
!# 4.) cc_writeUCD
!#     -> Write UCD (AVS, GMV,...) output.
!#
!# 5.) cc_writeFilm
!#     -> Write film output (=raw solution vectors).
!#
!# </purpose>
!##############################################################################

module ccpostprocessing

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc  
  use bcassembly
  use triangulation
  use linearalgebra  
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use ccboundaryconditionparser
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use linearsolver
  use domainintegration
  use statistics
  use element
  use vectorio
  use transformation
  use geometry
  use collection
  use convection
  use analyticprojection
  use ucd
  use feevaluation
  use pprocnavierstokes
  use pprocerror
  use elementpreprocessing
  use dofmapping
  use ccbasic
  use cccallback
  
  implicit none
  
  integer(I32), parameter :: CMASM_NELEMSIM               = 1000  
  
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
    
    ! here we will store the inner points of the particle
    type(t_vectorScalar) :: rvectorScalarParticle
  
  end type

!</typeblock>

!</types>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_postprocessingStationary (rproblem,rvector,rpostprocessing)
  
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

    ! local variables
    type(t_timer) :: rtimer

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    call cc_velocityParticle(rpostprocessing,rvector,rproblem)   
    ! Calculate body forces.
    call cc_calculateBodyForces (rvector,rproblem)
    
    ! Calculate the divergence
    call cc_calculateDivergence (rvector,rproblem)
    
    ! Error analysis, comparison to reference function.
    call cc_errorAnalysis (rvector,rproblem)
    
    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    call cc_writeUCD (rpostprocessing, rvector, rproblem)
    
    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimePostprocessing = &
      rproblem%rstatistics%dtimePostprocessing + rtimer%delapsedReal
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_postprocessingNonstat (rproblem,rvector,rpostprocessing)
  
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

    ! local variables
    type(t_timer) :: rtimer

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    call cc_velocityParticle(rpostprocessing,rvector,rproblem)
    ! Calculate body forces.
    call cc_calculateBodyForces (rvector,rproblem)
    
    ! Calculate the divergence
    call cc_calculateDivergence (rvector,rproblem)

    ! Error analysis, comparison to reference function.
    call cc_errorAnalysis (rvector,rproblem)
    
    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    call cc_writeUCD (rpostprocessing, rvector, rproblem, &
        rproblem%rtimedependence%dtime)
    
    ! Write film output (raw data vectors)
    call cc_writeFilm (rpostprocessing, rvector, rproblem, &
        rproblem%rtimedependence%dtime)
        
    call cc_updateParticlePosition(rpostprocessing,rproblem,&
         rproblem%rtimedependence%dtimestep)                
    
    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimePostprocessing = &
      rproblem%rstatistics%dtimePostprocessing + rtimer%delapsedReal

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_errorAnalysis (rsolution,rproblem)

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
    real(DP) :: derrorVel, derrorP, denergy
    integer :: icalcL2,icalcH1,icalcEnergy
    
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                     'IERRORANALYSISL2', icalcL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                     'IERRORANALYSISH1', icalcH1, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                     'ICALCKINETICENERGY', icalcEnergy, 1)
    
    if ((icalcL2 .ne. 0) .or. (icalcH1 .ne. 0) .or. (icalcEnergy .ne. 0)) then
      call output_lbrk()
      call output_line ('Error Analysis')
      call output_line ('--------------')
    end if
    
    if (icalcL2 .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u-z||_{L^2}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      derrorVel = sqrt(0.5_DP*(Derr(1)**2+Derr(2)**2))

      call pperr_scalar (rsolution%RvectorBlock(3),PPERR_L2ERROR,Derr(3),&
                         ffunction_TargetP,rproblem%rcollection)

      derrorP = sqrt(Derr(3))
      
      call output_line ('||u-reference||_L2 = '//trim(sys_sdEP(derrorVel,15,6)) )
      call output_line ('||p-reference||_L2 = '//trim(sys_sdEP(derrorP,15,6)) )
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end if

    if (icalcH1 .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add ||u-z||_{H^1}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_H1ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_H1ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      derrorVel = (0.5_DP*(Derr(1)**2+Derr(2)**2))

      call output_line ('||u-reference||_H1 = '//trim(sys_sdEP(derrorVel,15,6)) )
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end if
    
    if (icalcL2 .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u||^2_{L^2}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1))
      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2))
                         
      denergy = 0.5_DP*(Derr(1)**2+Derr(2)**2)

      call output_line ('||u||^2_L2         = '//trim(sys_sdEP(denergy,15,6)) )
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end if

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_calculateBodyForces (rsolution,rproblem)

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
    if ((rsolution%p_rblockDiscr%RspatialDiscr(1)% &
         ccomplexity .eq. SPDISC_UNIFORM) .and. &
        (boundary_igetNBoundComp(rproblem%rboundary) .ge. 2)) then

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      call boundary_createRegion (rproblem%rboundary, &
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

  subroutine cc_calculateDivergence (rsolution,rproblem)

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
    integer(I32) :: ieltype
    type(t_matrixScalar) :: rBmatrix
    type(t_vectorScalar), target :: rtempVector
    
    if (rsolution%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rsolution%p_rblockDiscr%RspatialDiscr(1)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
      
        ! Create a temporary vector 
        call lsyssc_createVecByDiscr (rsolution%RvectorBlock(3)%p_rspatialDiscr,&
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

  subroutine cc_writeUCD (rpostprocessing,rvector,rproblem,dtime)

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
    type(t_discreteBC), target :: rdiscreteBC
    type(t_discreteFBC), target :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    real(DP) :: dminTime, dmaxTime, dtimeDifferenceUCD
    integer :: ioutputUCD,ilevelUCD,ipolyHandle
    integer(I32) :: ieltype,ivt
    integer :: iin
    
    type(t_vectorScalar) :: vecDensity
    
    real(dp), dimension(:),pointer :: pDensity
    real(DP) :: dxcenter,dycenter,dradius,ddist
    real(dp), dimension(:,:),pointer :: p_Dvertices
    
    type(t_geometryObject), pointer :: p_rgeometryObject
    
    real(DP), dimension(:,:), pointer :: p_Dvertices1    
    
    character(SYS_STRLEN) :: sfile,sfilename
    
    ! get a pointer to the geometry object
    p_rgeometryObject => collct_getvalue_geom (rproblem%rcollection, 'mini')        
    
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
    
    call spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscr,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(2), &
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
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Initialise the discrete BC structure
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0 
    ! discretisation for implementing them into a solution vector.
    call cc_assembleBDconditions (rproblem,rprjDiscretisation,&
        rdiscreteBC,rproblem%rcollection,.true.)
                            
    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => rdiscreteBC
    
    ! The same way, discretise boundary conditions of fictitious boundary components.
    call bcasm_initDiscreteFBC(rdiscreteFBC)
    call cc_assembleFBDconditions (rproblem,rprjDiscretisation,&
        rdiscreteFBC,rproblem%rcollection)
    rprjVector%p_rdiscreteBCfict => rdiscreteFBC
    
    ! Filter the solution vector to implement discrete BC's.
    call vecfil_discreteBCsol (rprjVector)

    ! Filter the solution vector to implement discrete BC's for fictitious 
    ! boundary components.
    call vecfil_discreteFBCsol (rprjVector)
    
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

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
    
    call lsyssc_createVector (vecDensity,p_rtriangulation%NVT,.true.,ST_DOUBLE)
    call lsyssc_getbase_double (vecDensity,pDensity)

    call storage_getbase_double2D (p_rtriangulation%h_DvertexCoords,&
        p_DVertices)
    
    
    
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
    
    do ivt=1,p_rtriangulation%NVT
    
      call geom_isInGeometry (p_rgeometryObject, p_DVertices(:,ivt), iin)
      if(iin .eq. 1)then 
        pDensity(ivt) = rproblem%drho2
        p_Ddata(ivt)   = rproblem%du1
        p_Ddata2(ivt)  = rproblem%du2
      else
        pDensity(ivt) = rproblem%drho1
      end if
      
    end do      
    
    
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
        
    call ucd_addVariableVertexBased (rexport,'Density',UCD_VAR_STANDARD, &
        pDensity(1:p_rtriangulation%NVT))
        
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    if (rvector%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(1)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
          
        call ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
        
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'streamfunction',&
            UCD_VAR_STANDARD, p_Ddata(1:p_rtriangulation%NVT))
            
      end if
      
    end if
    
    call geom_polygonise(p_rgeometryObject,ipolyHandle)
    
    ! Get the vertices
    call storage_getbase_double2D(ipolyHandle, p_Dvertices1)
    
    call ucd_addPolygon(rexport,p_Dvertices1)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
    
    call lsyssc_releaseVector (vecDensity)    
    
    call storage_free(ipolyHandle)
    
    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    call bcasm_releaseDiscreteBC (rdiscreteBC)
    call bcasm_releaseDiscreteFBC (rdiscreteFBC)
    
    if (present(dtime)) then
      ! Update time stamp of last written out GMV.
      rpostprocessing%dnextTimeUCD = rpostprocessing%dnextTimeUCD+dtimeDifferenceUCD
      rpostprocessing%inextFileSuffixUCD = rpostprocessing%inextFileSuffixUCD + 1
    end if

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_writeFilm (rpostprocessing,rvector,rproblem,dtime)

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
    integer :: NEQ
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
      call output_line ('Warning: Level for solution vector is < NLMIN! &
          &Writing out at level NLMIN!', &
          OU_CLASS_WARNING,OU_MODE_STD,'cc_releasePreconditioner')
      call sys_halt()
      ilevelFilm = rproblem%NLMIN
    end if
    
    ! Write formatted output?
    bformatted = ioutputFilm .ne. 2

    ! Interpolate the solution down to level istart.
    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    do ilev = rproblem%NLMAX,ilevelFilm+1,-1
      
      ! Initialise a vector for the lower level and a prolongation structure.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev-1)%rdiscretisation,rvector2,.false.)
      
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

  subroutine cc_copyPostprocessing (rpostprocessingSrc,rpostprocessingDst)

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

  subroutine cc_clearpostprocessing (rpostprocessing)

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
    if (rpostprocessing%rvectorScalarParticle%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorScalarParticle)
      

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
  subroutine cc_velocityParticle(rpostprocessing,rvector,rproblem)
  !<description>
  !
  ! In this routine we calculate the velocity of the particle(s)
  !
  !</description>

  ! structure for a geometry object
  !<inputoutput>
  type(t_problem), intent(INOUT) :: rproblem
  type (t_c2d2postprocessing),intent(inout) :: rpostprocessing
  !</inputoutput>  

  !<input>
  type(t_vectorBlock), intent(IN) :: rvector
  !</input>
  
  !</subroutine>

  ! Local variables
  ! pointer to the entries of the alpha vector  
  real(DP), dimension(:), pointer :: p_Dvector,p_Dvector1,p_DvectorU1,p_Dvector2,p_DvectorU2
  real(DP), dimension(:), pointer :: p_DvectorRO1, p_DvectorRO2
  ! pointer to the nodes of the grid
  real(DP), dimension(:,:), pointer :: p_Ddata

  ! pointer to the triangulation structure
  type(t_triangulation), pointer :: p_rtriangulation
  
  type(t_vectorblock) :: rvectorAlpha
  
  type(t_vectorscalar) :: rvectorDupVec1
  
  type(t_vectorscalar) :: rvectorDupVec2    
  
  type(t_vectorscalar) :: rvectorProj1    
  type(t_vectorscalar) :: rvectorProj2    
  
  integer :: i
  
  real(dp) :: dIntU1,dIntU2,dvolume,dmass,dIntU3,dIntU4,dIntR1,dIntR2,angVel,ddmass

  if (rpostprocessing%rvectorScalarParticle%NEQ .ne. 0) &
    call lsyssc_releaseVector (rpostprocessing%rvectorScalarParticle)
  ! create a vector from the velocity discretisation
  call lsyssc_createVecByDiscr(rvector%RvectorBlock(1)%p_rspatialDiscr, &
  rpostprocessing%rvectorScalarParticle,.true.)

  call lsyssc_createVecByDiscr(rvector%RvectorBlock(1)%p_rspatialDiscr, &
  rvectorProj1,.true.)

  call lsyssc_createVecByDiscr(rvector%RvectorBlock(2)%p_rspatialDiscr, &
  rvectorProj2,.true.)

  
  ! make an L2 projection 
  ! insert the cc_particle as the callback function
  call anprj_discrDirect (rpostprocessing%rvectorScalarParticle,cc_Particle,&
                          rproblem%rcollection,iorder=1)  
  
  ! get a pointer to the entries of this vector
  call lsyssc_getbase_double(rpostprocessing%rvectorScalarParticle,p_Dvector)

  ! copy the u1 velocity to rvectorDupVec1
  call lsyssc_duplicateVector (rvector%RvectorBlock(1),rvectorDupVec1,&
                               LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)

  ! copy the u2 velocity to rvectorDupVec2
  call lsyssc_duplicateVector (rvector%RvectorBlock(2),rvectorDupVec2,&
                               LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)

                               
  call lsyssc_getbase_double(rvector%RvectorBlock(1),p_DvectorU1)                                                                
  call lsyssc_getbase_double(rvector%RvectorBlock(2),p_DvectorU2)

  call lsyssc_getbase_double(rvectorProj1,p_DvectorRO1)                                                                
  call lsyssc_getbase_double(rvectorProj2,p_DvectorRO2)

  
  ! get the entries of the dublicated vectors                                                                
  call lsyssc_getbase_double(rvectorDupVec1,p_Dvector1)                                 
  call lsyssc_getbase_double(rvectorDupVec2,p_Dvector2)                                   

  ! loop through all entries and multiply
  ! by the 0-1 vector
!  do i=1,rvectorDupVec1%NEQ
!    p_Dvector1(i) = p_Dvector1(i) * p_Dvector(i)
!    p_Dvector2(i) = p_Dvector2(i) * p_Dvector(i)    
!  end do

  dIntU1 = 0.0_dp
  dIntU2 = 0.0_dp
  dIntU3 = 0.0_dp
  dIntU4 = 0.0_dp

  
  ! Integrate:
  ! int_Particle (\widetilde{u}) dParticle
  call pperr_scalar (rvectorDupVec1,PPERR_L1ERROR,dIntU1)
  call pperr_scalar (rvectorDupVec2,PPERR_L1ERROR,dIntU2)
  
  call cc_velocityUpdate (rvectorDupVec1,PPERR_L1ERROR,dIntU3,&
                         rvector%RvectorBlock(1)%p_rspatialdiscr,&
                         rpostprocessing%rvectorScalarParticle)  
  call cc_velocityUpdate (rvectorDupVec2,PPERR_L1ERROR,dIntU4,&
                         rvector%RvectorBlock(2)%p_rspatialdiscr,&
                         rpostprocessing%rvectorScalarParticle)  
  
  call cc_torque (rvectorDupVec1,rvectorDupVec2,PPERR_L1ERROR, dIntR1,dIntR2,&
                                  rvector%RvectorBlock(1)%p_rspatialdiscr,&
                                  rpostprocessing%rvectorScalarParticle,&
                                  rproblem%rcollection)
                                  
  ! compute the volume and mass of the 
  ! particle
  dvolume  = (rproblem%drad)**2 * SYS_PI
  dmass  = rproblem%drho2 * dvolume 
  ddmass= (rproblem%drho2-rproblem%drho1) * dvolume
  
  ! compute the
  dIntR2=dIntR2*dmass
  
  ! calculate rho2/mass * (int_particle u_tilde)
  ! in words: the mean velocity in the particle area or volume
  dIntU1 = (rproblem%drho2/dmass) * dIntU1
  dIntU2 = (rproblem%drho2/dmass) * dIntU2
  dIntU3 = (rproblem%drho2/dmass) * dIntU3
  dIntU4 = (rproblem%drho2/dmass) * dIntU4

  ! calculate the angular velocity
  angVel = dIntR1/dIntR2
  
  call output_line ('----------------')  
  print *,"Angular Velocity: ",angVel

  rproblem%dAngVel(2) = angVel

  ! put this value into the collection
  !rproblem%rcollection%DQuickaccess(6) = rproblem%dAngVel(1) + angVel
  rproblem%rcollection%DQuickaccess(6) = angVel
  
  
  ! make an L2 projection 
  ! insert the cc_particle as the callback function
  call anprj_discrDirect (rvectorProj1,cc_LinVelX,&
                          rproblem%rcollection,iorder=1)  

  ! make an L2 projection 
  ! insert the cc_particle as the callback function
  call anprj_discrDirect (rvectorProj2,cc_LinVelY,&
                          rproblem%rcollection,iorder=1)  
                   
                   
                          
  !print *, p_DvectorRO1                               
  ! here we write back the corrected velocity
  ! values to the solution vector for the next step
  do i=1,rvectorDupVec1%NEQ
    if(p_Dvector(i) .ne. 0.0_dp)then
      p_DvectorU1(i) = dIntU3 !+ ! + p_DvectorRO1(i)
      p_DvectorU2(i) = dIntU4 !+ rproblem%rtimedependence%dtimestep/dmass * &
      !(ddmass*-9.807_dp) ! + p_DvectorRO2(i)
    end if
  end do
  
  ! save the velocity
  rproblem%duold1 = rproblem%du1
  rproblem%duold2 = rproblem%du2
  
  ! save the current velocity
  rproblem%du1 = dIntU3
  rproblem%du2 = dIntU4  !+ rproblem%rtimedependence%dtimestep/dmass * &
      !(ddmass*-9.807_dp)
  
  ! release the temporary vectors
  call lsyssc_releaseVector(rvectorDupVec1)
  call lsyssc_releaseVector(rvectorDupVec2)
  call lsyssc_releaseVector(rvectorProj1)
  call lsyssc_releaseVector(rvectorProj2)

  call output_lbrk ()
  
  end subroutine  
  
! ***************************************************************************  
  
!<subroutine>
  subroutine cc_updateParticlePosition(rpostprocessing,rproblem,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  
  type (t_c2d2postprocessing),intent(inout) :: rpostprocessing
  type(t_problem), intent(INOUT) :: rproblem
  real(dp), intent(in) :: dtimestep
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,k,icurrentElementDistr, ICUBP, NVE
    integer(I32) :: IEL, IELmax, IELset
    real(DP) :: OM, DN1, DN2, DN3,dpp
    real(DP) :: ah1,ah2,ah3,dvelx,dvely,dvelz,dmasssl,ddmasssl,dvolume
    real(DP) :: dfx,dfy,du3x,du3y,du3z,dalx,daly,dalz,nennerX
    real(DP) :: dCenterX,dCenterY,dCenterXold,dCenterYold,dydt
    type(t_geometryObject), pointer :: p_rgeometryObject        
    
    p_rgeometryObject => collct_getvalue_geom (rproblem%rcollection, 'mini')

    print *,"--------------------------"
    
   
    
    !
    p_rgeometryObject%rcoord2D%Dorigin(1)= p_rgeometryObject%rcoord2D%Dorigin(1) + &
    dtimestep * rproblem%du1
    
    !
    p_rgeometryObject%rcoord2D%Dorigin(2)= p_rgeometryObject%rcoord2D%Dorigin(2) + &
    dtimestep * rproblem%du2 

    print *,"Y-Position: ",p_rgeometryObject%rcoord2D%Dorigin(2)

!    p_rgeometryObject%rcoord2D%drotation = p_rgeometryObject%rcoord2D%drotation + &
!                                         dtimestep *(rproblem%DAngVel(1) + &
!                                         0.5_dp* rproblem%DAngVel(2))
!
!    
!    ! Recalculate SIN and COS values of angle
!    p_rgeometryObject%rcoord2D%dsin_rotation = sin(p_rgeometryObject%rcoord2D%drotation)
!    p_rgeometryObject%rcoord2D%dcos_rotation = cos(p_rgeometryObject%rcoord2D%drotation)
    
    ! update the angular velocity
    !rproblem%DAngVel(1) = rproblem%DAngVel(1) + rproblem%DAngVel(2)
    rproblem%DAngVel(1) = rproblem%DAngVel(2)
    
    print *,"--------------------------"
    
    print *,"Domega: ",rproblem%DAngVel(2)
 
    print *,"--------------------------"
    
    print *,"Angular Velocity: ",rproblem%DAngVel(1)

    print *,"--------------------------"
    
    print *,"U1 Velocity: ",rproblem%du1

    print *,"--------------------------"
 
    print *,"U2 Velocity: ",rproblem%du2 
    
    
  end subroutine

! ***************************************************************************  

  subroutine cc_Particle(cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)                      :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>
  ! local variables
  real :: dCenterX,dCenterY,dRadius
  integer :: i,j,iin

  type(t_geometryObject), pointer :: p_rgeometryObject
  
  ! get a pointer to the geometry object
  p_rgeometryObject => collct_getvalue_geom (rcollection, 'mini')       

  dCenterX = p_rgeometryObject%rcoord2D%Dorigin(1)
  dCenterY = p_rgeometryObject%rcoord2D%Dorigin(2)


  select case (cderivative)
  case (DER_FUNC)
  ! 
  do i=1,nelements
    do j=1,npointsPerElement
      call geom_isInGeometry (p_rgeometryObject, Dpoints(:,j,i), iin)
      if(iin .eq. 1)then 
        Dvalues(j,i) =  1.0_DP 
      else
        Dvalues(j,i) = 0.0_DP
      end if
    end do
  end do
    
  case (DER_DERIV_X)
    ! Not really useful in the case at hand
    Dvalues (:,:) = 0.0_dp
  case (DER_DERIV_Y)
    ! not much better than the above case
    Dvalues (:,:) = 0.0_dp
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  
  end subroutine
  
!****************************************************************************  

  subroutine cc_LinVelX(cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)                      :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>
  ! local variables
  real :: dCenterX,dCenterY,dRadius,angVel,robx,roby,length
  integer :: i,j,iin

  type(t_geometryObject), pointer :: p_rgeometryObject
  
  ! get a pointer to the geometry object
  p_rgeometryObject => collct_getvalue_geom (rcollection, 'mini')
  
  dCenterX = p_rgeometryObject%rcoord2D%Dorigin(1)
  dCenterY = p_rgeometryObject%rcoord2D%Dorigin(2)
  
  angVel=rcollection%DQuickAccess(6)


  select case (cderivative)
  case (DER_FUNC)
  ! 
  do i=1,nelements
    do j=1,npointsPerElement
      call geom_isInGeometry (p_rgeometryObject, Dpoints(:,j,i), iin)
      if(iin .eq. 1)then
        robx = -1.0_dp*(Dpoints(2,j,i)-dCenterY)
        roby = (Dpoints(1,j,i)-dCenterX)
        length = sqrt(robx*robx+roby*roby)
        Dvalues(j,i) = angVel * robx/length
      else
        Dvalues(j,i) = 0.0_DP
      end if
    end do
  end do
    
  case (DER_DERIV_X)
    ! Not really useful in the case at hand
    Dvalues (:,:) = 0.0_dp
  case (DER_DERIV_Y)
    ! not much better than the above case
    Dvalues (:,:) = 0.0_dp
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  
  end subroutine

!****************************************************************************  

  subroutine cc_LinVelY(cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)                      :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>
  ! local variables
  real :: dCenterX,dCenterY,dRadius,angVel,roby,robx,length
  integer :: i,j,iin

  type(t_geometryObject), pointer :: p_rgeometryObject
  
  ! get a pointer to the geometry object
  p_rgeometryObject => collct_getvalue_geom (rcollection, 'mini')
  
  dCenterX = p_rgeometryObject%rcoord2D%Dorigin(1)
  dCenterY = p_rgeometryObject%rcoord2D%Dorigin(2)
  
  angVel=rcollection%DQuickAccess(6)

  select case (cderivative)
  case (DER_FUNC)
  ! 
  do i=1,nelements
    do j=1,npointsPerElement
      call geom_isInGeometry (p_rgeometryObject, Dpoints(:,j,i), iin)
      if(iin .eq. 1)then
        robx = -1.0_dp*(Dpoints(2,j,i)-dCenterY)
        roby = (Dpoints(1,j,i)-dCenterX)
        length = sqrt(robx*robx+roby*roby)
        Dvalues(j,i) = angVel * roby
      else
        Dvalues(j,i) = 0.0_DP
      end if
    end do
  end do
    
  case (DER_DERIV_X)
    ! Not really useful in the case at hand
    Dvalues (:,:) = 0.0_dp
  case (DER_DERIV_Y)
    ! not much better than the above case
    Dvalues (:,:) = 0.0_dp
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  
  end subroutine


!****************************************************************************

!<subroutine>

  subroutine cc_velocityUpdate (rvectorScalar, cerrortype, derror,&
                                  rdiscretisation,&
                                  rvectorParticle,&
                                  rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalar
  type(t_vectorScalar), intent(IN), target :: rvectorParticle
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation

!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional 
  ! information for callback routines.
  type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    integer :: icurrentElementDistr, ICUBP, NVE
    integer :: IEL, IELmax, IELset, IELGlobal
    real(DP) :: OM
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(:), allocatable :: Domega
    
    ! number of cubature points on the reference element
    integer :: ncubp
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial,indofFunc
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    type(t_elementDistribution), pointer :: p_relementDistributionA
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
  
    ! Type of transformation from the reference to the real element 
    integer(I32) :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Pointer to the element-wise error
    real(DP), dimension(:), pointer :: p_Derror

    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(CMASM_NELEMSIM, p_rtriangulation%NEL)
    
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1, rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
      
      p_relementDistributionA =>&
      rvectorParticle%p_rspatialDiscr%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      indofFunc = elem_igetNDofLoc(p_relementDistributionA%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Get the number of cubature points for the cubature formula
      ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeEval)
      
      ! Allocate two arrays for the points and the weights
      allocate(Domega(ncubp))
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
      
      ! Get the cubature formula
      call cub_getCubature(p_relementDistribution%ccubTypeEval, p_DcubPtsRef, Domega)
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial, nelementsPerBlock))

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp, nelementsPerBlock, 4))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, CMASM_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+CMASM_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
                                     
        ! Prepare the call to the evaluation routine of the analytic function.    
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
        rintSubset%celement = p_relementDistribution%celement
    
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        ! L1-error uses only the values of the function.
        Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP

        ! Calculate the values of the FE function in the
        ! cubature points: u_h(x,y).
        ! Save the result to Dcoefficients(:,:,2)
        
        call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                Dcoefficients(:,1:IELmax-IELset+1,2))

        call fevl_evaluate_sim3 (rvectorParticle, revalElementSet,&
                p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                Dcoefficients(:,1:IELmax-IELset+1,3))

        
        ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
        ! the error "u-u_h(cubature pt.)"!
        !        
        ! Loop through elements in the set and for each element,
        ! loop through the DOF's and cubature points to calculate the
        ! integral: int_Omega abs(u-u_h) dx
          
        do IEL=1,IELmax-IELset+1
          
          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp
            
            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            
            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
            
            ! L1-error is:   int_... abs(u-u_h) dx
            
            derror = derror + OM * Dcoefficients(icubp,IEL,2) * Dcoefficients(icubp,IEL,3)
            
          end do ! ICUBP 
          
        end do ! IEL
        
      ! Release the temporary domain integration structure again
      call domint_doneIntegration (rintSubset)
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      deallocate(Domega)

    end do ! icurrentElementDistr

  end subroutine cc_velocityUpdate

!****************************************************************************

!<subroutine>
  subroutine cc_torque (rvectorScalarU1,rvectorScalarU2,cerrortype, int1,int2,&
                                  rdiscretisation,&
                                  rvectorParticle,&
                                  rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalarU1
  type(t_vectorScalar), intent(IN), target :: rvectorScalarU2
  type(t_vectorScalar), intent(IN), target :: rvectorParticle
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation

!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional 
  ! information for callback routines.
  type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(OUT) :: int1
  real(DP), intent(OUT) :: int2
!</output>

!</subroutine>

    ! local variables
    integer :: icurrentElementDistr, ICUBP, NVE
    integer :: IEL, IELmax, IELset, IELGlobal
    real(DP) :: OM
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(:), allocatable :: Domega
    
    ! number of cubature points on the reference element
    integer :: ncubp
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
  
    ! Type of transformation from the reference to the real element 
    integer(I32) :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Pointer to the element-wise error
    real(DP), dimension(:), pointer :: p_Derror
    type(t_geometryObject), pointer :: p_rgeometryObject
    real(dp) :: dcenterx,dcentery,robx,roby,ut1,ut2,ah,rx,ry  
    
    p_rgeometryObject => collct_getvalue_geom(rcollection,'mini')

    dcenterx = p_rgeometryObject%rcoord2D%Dorigin(1)
    dcentery = p_rgeometryObject%rcoord2D%Dorigin(2)

    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(CMASM_NELEMSIM, p_rtriangulation%NEL)
    
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    int1 = 0.0_DP
    int2 = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1, rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Get the number of cubature points for the cubature formula
      ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeEval)
      
      ! Allocate two arrays for the points and the weights
      allocate(Domega(ncubp))
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
      
      ! Get the cubature formula
      call cub_getCubature(p_relementDistribution%ccubTypeEval, p_DcubPtsRef, Domega)
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial, nelementsPerBlock))

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp, nelementsPerBlock, 4))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)
      
      ! get the real points
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, CMASM_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+CMASM_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
                                     
        ! Prepare the call to the evaluation routine of the analytic function.    
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
        rintSubset%celement = p_relementDistribution%celement
    
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        ! L1-error uses only the values of the function.
        Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP

        ! Calculate the values of the FE function in the
        ! cubature points: u_h(x,y).
        ! Save the result to Dcoefficients(:,:,2)
        
        call fevl_evaluate_sim3 (rvectorScalarU1, revalElementSet,&
                p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                Dcoefficients(:,1:IELmax-IELset+1,2))

        call fevl_evaluate_sim3 (rvectorScalarU2, revalElementSet,&
                p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                Dcoefficients(:,1:IELmax-IELset+1,3))

        call fevl_evaluate_sim3 (rvectorParticle, revalElementSet,&
                p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                Dcoefficients(:,1:IELmax-IELset+1,4))

        
        ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
        ! the error "u-u_h(cubature pt.)"!
        !        
        ! Loop through elements in the set and for each element,
        ! loop through the DOF's and cubature points to calculate the
        ! integral: int_Omega abs(u-u_h) dx
          
        do IEL=1,IELmax-IELset+1
          
          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp
            
            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            
            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
            
            ! we want to calculate: int_Particle (r_ob .perpdot. u_tilde) dParticle and
            ! int_Particle |r|^2 dParticle
            
            ! calculate r_ob
            robx = -1.0_dp*(revalElementSet%p_DpointsReal(2,icubp,iel) - dcentery)
            roby = (revalElementSet%p_DpointsReal(1,icubp,iel)  - dcenterx)
            
            ! calculate u_tilde
            ut1 =  Dcoefficients(icubp,IEL,2) * Dcoefficients(icubp,IEL,4)
            ut2 =  Dcoefficients(icubp,IEL,3) * Dcoefficients(icubp,IEL,4)
            
            ! calculate the perpdot
            ah = robx * ut1 + roby * ut2
            
            ! add upp and multiply by the cubature weight
            int1 = int1 + ah * OM

            ! calculate the 2d integral int_Particle r_ob .dot. r_ob
            int2 = int2 + (robx*robx + roby*roby) * OM
            
          end do ! ICUBP 
          
        end do ! IEL
        
      ! Release the temporary domain integration structure again
      call domint_doneIntegration (rintSubset)
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      deallocate(Domega)

    end do ! icurrentElementDistr

  end subroutine cc_torque

  
! ***************************************************************************
  

end module
