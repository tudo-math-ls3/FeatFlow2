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
!# 5.) cc_forces
!#     -> calculate Drag/Lift forces with penalty/FB methods
!#        and volume integration.
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
!# 6.) cc_forcesIntegration
!#     -> performs the actual integration for cc_forces
!#
!# </purpose>
!##############################################################################

module ccpostprocessing

  use fsystem
  use storage
  use linearsolver
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
  use ccboundaryconditionparser
  use analyticprojection
  use collection
  use convection
  use geometry
  use ucd
  
  use pprocnavierstokes
  use pprocerror
  
  use ccbasic
  use cccallback
  
  implicit none
  
!<types>

!<typeblock>

  ! Postprocessing structure. Whenever the postprocessing calculation routines are
  ! called, they fill this structure using the current solution vector (and other
  ! information if necessary). The information in this structure can then be used
  ! for GMV output e.g.
  type t_cc3dpostprocessing

    ! A discretisation structure that describes a piecewise constant discretisation
    ! (usually P0 or Q0).
    type(t_spatialDiscretisation) :: rdiscrConstant
    
    ! A discretisation structure that describes a piecewise linear discretisation
    ! (usually P1 or Q1).
    type(t_spatialDiscretisation) :: rdiscrLinear

    ! A discretisation structure that describes a piecewise quadratic discretisation
    ! (usually P2 or Q2) - currently not implemented in 3D
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

    ! A vector that describes the Z-velocity field in the vertices
    type(t_vectorScalar) :: rvectorVelZ

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
    
    type(t_vectorScalar) :: rvectorScalarFB
    
    type(t_vectorScalar) :: rvectorScalarFBQ1
    
    type(t_vectorScalar) :: rResForceX
    
    type(t_vectorScalar) :: rResForceY
    
    type(t_vectorScalar) :: rResForceZ
    
  
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

  ! Postprocvessing structure. 
  type(t_cc3dpostprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!</subroutine>
    
    ! Drag/Lift Calculation
    call cc_forces(rpostprocessing,rvector,rproblem)

    ! Calculate body forces.
    call cc_calculateBodyForces (rvector,rproblem)
    
    ! Calculate the divergence
    call cc_calculateDivergence (rvector,rproblem)
    
    ! Error analysis, comparison to reference function.
    call cc_errorAnalysis (rvector,rproblem)
    
    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    call cc_writeUCD (rpostprocessing, rvector, rproblem)
    
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

  ! Postprocvessing structure. Defines what to do with solution vectors.
  type(t_cc3dpostprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!</subroutine>

    ! Drag/Lift Calculation
    call cc_forcesNonStat(rpostprocessing,rvector,rproblem)

    ! Calculate body forces.
    call cc_calculateBodyForces (rvector,rproblem)
    
    ! Calculate the divergence
    call cc_calculateDivergence (rvector,rproblem)

    ! Error analysis, comparison to reference function.
    call cc_errorAnalysis (rvector,rproblem)
    
    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    call cc_writeUCD (rpostprocessing, rvector, rproblem, &
        rproblem%rtimedependence%dtime)
        
    call cc_updateParticlePosition(rpostprocessing,rproblem,&
         rproblem%rtimedependence%dtimestep)        
    
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
    real(DP),dimension(4) :: Derr
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
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u-z||_{L^2}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      call pperr_scalar (rsolution%RvectorBlock(3),PPERR_L2ERROR,Derr(3),&
                         ffunction_TargetZ,rproblem%rcollection)

      derrorVel = (0.5_DP*(Derr(1)**2+Derr(2)**2+Derr(3)**2))

      call pperr_scalar (rsolution%RvectorBlock(4),PPERR_L2ERROR,Derr(4),&
                         ffunction_TargetP,rproblem%rcollection)

      derrorP = Derr(4)
      
      call output_line ('||u-reference||_L2 = '//trim(sys_sdEP(derrorVel,15,6)) )
      call output_line ('||p-reference||_L2 = '//trim(sys_sdEP(derrorP,15,6)) )
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end if

    if (icalcH1 .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u-z||_{L^2}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_H1ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_H1ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      call pperr_scalar (rsolution%RvectorBlock(3),PPERR_H1ERROR,Derr(3),&
                         ffunction_TargetY,rproblem%rcollection)

      derrorVel = (0.5_DP*(Derr(1)**2+Derr(2)**2+Derr(3)**2))

      call output_line ('||u-reference||_H1 = '//trim(sys_sdEP(derrorVel,15,6)) )
      
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
    
    ! Forces on the object
    real(DP), dimension(NDIM3D) :: Dforces
    real(DP) :: df1,df2
    type(t_meshRegion) :: rregion
    integer, dimension(1) :: Iregion = (/ 7 /)
    
    ! If we have a uniform discretisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    if ((rsolution%p_rblockDiscr%RspatialDiscr(1)% &
         ccomplexity .eq. SPDISC_UNIFORM)) then
         
      ! Calculate a mesh region for the seventh boundary component
      call ccdc_calcBoundaryMeshRegion(rproblem,rregion,&
        rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation, Iregion)
      
      ! Is the mesh region empty?
      if (rregion%NAT .le. 0) then
        call mshreg_done(rregion)
        return
      end if

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      df1 = 1.0_DP/1000.0_DP
      df2 = 0.041_DP * 0.2_DP**2
      call ppns3D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G1X1,df1,df2)

      call output_lbrk()
      call output_line ('Body forces')
      call output_line ('-----------')
      call output_line ('Body forces real bd., bdc/horiz/vert')
      call output_line (' 2 / ' &
          //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //trim(sys_sdEP(Dforces(2),15,6)) // ' / '&
          //trim(sys_sdEP(Dforces(3),15,6)) )
      
      ! Destroy the mesh region
      call mshreg_done(rregion)
      
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
    integer :: ieltype
    type(t_matrixScalar) :: rBmatrix
    type(t_vectorScalar), target :: rtempVector
    
    if (rsolution%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rsolution%p_rblockDiscr%RspatialDiscr(1)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T_3D) then
      
        ! Create a temporary vector 
        call lsyssc_createVecByDiscr (rsolution%RvectorBlock(4)%p_rspatialDiscr,&
            rtempVector,.true.)

        ! Calculate divergence = B1^T u1 + B2^T u2 + B3^T u3
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
        call lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB3,&
            rBmatrix,LSYSSC_TR_VIRTUAL)
        call lsyssc_scalarMatVec (&
            rBmatrix, rsolution%RvectorBlock(3), &
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
  type(t_cc3dpostprocessing), intent(INOUT) :: rpostprocessing  
!</inputoutput>

!</subroutine>

    ! local variables

    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2,p_Ddata3,p_Ddata4,p_Ddata5

    ! A pointer to the triangulation.
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector
    
    ! A vector accepting Q1 data
    type(t_vectorScalar) :: rprjVectorScalar
    
    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! A discretisation structure for Q1
    type(t_spatialDiscretisation) :: rprjDiscretisation1    
    
    ! Discrete boundary conditions for the output vector
    type(t_discreteBC), target :: rdiscreteBC
    type(t_discreteFBC), target :: rdiscreteFBC
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    real(DP) :: dminTime, dmaxTime, dtimeDifferenceUCD
    integer :: ioutputUCD,ieltype,ilevelUCD,ii,iloop,hpolyHandle
    integer :: ipoly,hpolyHandle1
    
    type(t_geometryObject) :: rgeometryObject 
    
    character(SYS_STRLEN) :: sfile,sfilename
    
    !real(dp), dimension(:,:), pointer :: Dvertices     
    real(dp), dimension(3,3) :: Dvertices     
    
    call geom_init_sphere(rgeometryObject,0.05_dp,(/0.0_dp,0.0_dp,0.0_dp/))
    
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
                 EL_Q1_3D, CUB_G2_3D, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_Q1_3D, CUB_G2_3D, &
                 rprjDiscretisation%RspatialDiscr(2))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(3), &
                 EL_Q1_3D, CUB_G2_3D, &
                 rprjDiscretisation%RspatialDiscr(3))
                 
    call spdiscr_deriveSimpleDiscrSc (&
                 rpostprocessing%rvectorScalarFB%p_rspatialDiscr,&
                 EL_Q1_3D, CUB_G2_3D, &
                 rprjDiscretisation1)
                 
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    call lsyssc_createVecByDiscr(rprjDiscretisation1,rprjVectorScalar,.true.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)
    
    call spdp_projectSolutionScalar(rpostprocessing%rvectorScalarFB,rprjVectorScalar)

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Initialise the discrete BC structure
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0 
    ! discretisation for implementing them into a solution vector.
    call cc_assembleBDconditions (rproblem,rprjDiscretisation,&
        rdiscreteBC,rproblem%rcollection)
                            
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
    
    ! Start UCD export to GMV file:
    call output_lbrk ()
    call output_line ('Writing GMV file: '//sfile)
    
    select case (ioutputUCD)
    case (1)
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
      
      call geom_sphere_polygonise(rgeometryObject, hpolyHandle,hpolyHandle1)

      
      ipoly = rgeometryObject%rsphere%ipoly
      do iloop=1,ipoly
        ! mach das p_DVertices array
        call geom_sphere_getPolygon(rgeometryObject,iloop,hpolyHandle,&
        hpolyHandle1,Dvertices)
        call ucd_addPolygon(rexport,Dvertices,4)
      end do
        call storage_free(hpolyHandle)
        call storage_free(hpolyHandle1)

    case (2)
      call ucd_startAVS (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
          
    case (3)
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
          
    case DEFAULT
      call output_line ('Invalid UCD output type.', &
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
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata3)
    call lsyssc_getbase_double (rpostprocessing%rvectorScalarFBQ1,p_Ddata4)
    call lsyssc_getbase_double (rprjVectorScalar,p_Ddata5)
    
!    do ii=1,p_rtriangulation%NVT
!      if(p_Ddata5(ii) .gt. 0.0_dp)then
!        p_Ddata5(ii)=1.0_dp
!      end if
!    end do
    
    ! CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, &
    !     p_Ddata(1:p_rtriangulation%NVT))
    ! CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, &
    !     p_Ddata2(1:p_rtriangulation%NVT))
    ! CALL ucd_addVariableVertexBased (rexport,'Z-vel',UCD_VAR_ZVELOCITY, &
    !     p_Ddata3(1:p_rtriangulation%NVT))
    call ucd_addVarVertBasedVec (rexport,'velocity',&
        p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT),&
        p_Ddata3(1:p_rtriangulation%NVT))
    
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, &
        p_Ddata(1:p_rtriangulation%NEL))
        
    call ucd_addVariableVertexBased (rexport,'FBM',UCD_VAR_STANDARD, &
         p_Ddata4(1:p_rtriangulation%NVT))        

    call ucd_addVariableVertexBased (rexport,'FBM1',UCD_VAR_STANDARD, &
         p_Ddata5(1:p_rtriangulation%NVT))        
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
!    IF (rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
!        ccomplexity .EQ. SPDISC_UNIFORM) THEN
!        
!      ieltype = rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
!                RelementDistr(1)%itrialElement
!                
!      IF (elem_getPrimaryElement(ieltype) .EQ. EL_Q1T_3D) THEN
!          
!        CALL ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
!        
!        CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
!        CALL ucd_addVariableVertexBased (rexport,'streamfunction',&
!            UCD_VAR_STANDARD, p_Ddata(1:p_rtriangulation%NVT))
!            
!      END IF
!      
!    END IF
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
    call lsyssc_releaseVector (rprjVectorScalar)
    
    
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
  type(t_cc3dpostprocessing), intent(INOUT) :: rpostprocessing  
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
  ! Postprocvessing structure.
  type(t_cc3dpostprocessing), intent(OUT) :: rpostprocessing
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
                 EL_Q0_3D, CUB_G1_3D, &
                 rpostprocessing%rdiscrConstant)

    ! Piecewise linear space:
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscr(1), &
                 EL_Q1_3D, CUB_G2_3D, &
                 rpostprocessing%rdiscrLinear)
  
    ! Piecewise quadratic space:
!    CALL spdiscr_deriveSimpleDiscrSc (&
!                 p_rdiscr%RspatialDiscr(1), &
!                 EL_Q2, CUB_G3X3, &
!                 rpostprocessing%rdiscrQuadratic)
  
    ! init vector
    call lsyssc_createVecByDiscr(rpostprocessing%rdiscrLinear,&
         rpostprocessing%rvectorScalarFBQ1,.true.)
  
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
    
    call lsyssc_createVector(rpostprocessing%rResForceX,6,.true.)
    call lsyssc_createVector(rpostprocessing%rResForceY,6,.true.)
    call lsyssc_createVector(rpostprocessing%rResForceZ,6,.true.)
                                    
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
  type(t_cc3dpostprocessing), intent(IN) :: rpostprocessingSrc
!</input>

!<inputoutput>  
  ! Destination Postprocessing structure.
  type(t_cc3dpostprocessing), intent(INOUT) :: rpostprocessingDst
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
  
! ***************************************************************************  

  !<subroutine>
  subroutine cc_forces(rpostprocessing,rvector,rproblem)
  !<description>
  ! A routine that calculated the forces acting on an object
  ! and moves the object according to these forces
  !</description>

  ! structure for a geometry object
  !<inputoutput>
  type(t_problem), intent(INOUT) :: rproblem
  type (t_cc3dpostprocessing),intent(inout) :: rpostprocessing
  !</inputoutput>  

  !<input>
  type(t_vectorBlock), intent(IN) :: rvector
  !</input>
  
  !</subroutine>

  ! Local variables
  ! pointer to the entries of the alpha vector  
  real(DP), dimension(:), pointer :: p_Dvector  

  ! pointer to the nodes of the grid
  real(DP), dimension(:,:), pointer :: p_Ddata
  
  ! Forces
  real(DP) :: DResForceX,DResForceY,DResForceZ
  real(DP) :: ddist

  ! pointer to the triangulation structure
  type(t_triangulation), pointer :: p_rtriangulation
  
  type(t_vectorblock) :: rvectorAlpha

  integer :: ive,NEL,itest

  real(DP) :: mytime,df1,df2,dxcenter,dycenter,dzcenter,dradius
  
  !mytime = dtime
 
  call lsyssc_createVecByDiscr(rvector%RvectorBlock(1)%p_rspatialDiscr, &
  rpostprocessing%rvectorScalarFB,.true.)                 
  
  itest=rproblem%rcollection%IquickAccess(1)
  
  ! get a pointer to the triangulation
  p_rtriangulation => &
  rpostprocessing%rvectorScalarFB%p_rspatialDiscr%p_rtriangulation
  
  ! make an L2 Projection of the analytic function to a FEM vector
  call anprj_discrDirect (rpostprocessing%rvectorScalarFB,cc_CB,&
                          rproblem%rcollection,iorder=1)  
  
  print *,"Mids inside: ",rproblem%rcollection%IquickAccess(1)
  
  ! get a pointer to the entries of the vector
  call lsyssc_getbase_double(rpostprocessing%rvectorScalarFBQ1,p_Dvector)
  
  ! And get the vertice coordinates
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, p_Ddata)  

  ! Definition of the circle
  dxcenter = 0.5
  dycenter = 0.2
  dzcenter = 0.2
  dradius  = 0.05001

  do ive=1,p_rtriangulation%NVT
    ! Get the distance to the center
    ddist = sqrt( (p_Ddata(1,ive)-dxcenter)**2 + (p_Ddata(2,ive)-dycenter)**2 )
    ! Point inside?
    if (ddist .le. dradius) then
      p_Dvector(ive)=1.0_dp
    else
      p_Dvector(ive)=0.0_dp  
    end if
  end do
  

  call output_lbrk ()
  call output_separator(OU_SEP_EQUAL)
  call output_line ('Q1 Vector recalculated ')
  call output_separator(OU_SEP_EQUAL)   
  
  df1=0.001_dp
  df2=0.041_dp*(0.45_dp*4.0_dp/9.0_dp)**2
  call cc_forcesIntegration(rvector, rpostprocessing%rvectorScalarFB, DResForceX, &
  DResForceY,DResForceZ,rvectorAlpha,df1,df2)
  
  call output_lbrk()
  call output_line ('Drag forces')
  call output_line ('-----------') 
  call output_line (trim(sys_sdL(DResForceX,10))//"/"//trim(sys_sdL(DResForceY,10))&
                    //"/"//trim(sys_sdL(DResForceZ,10)) ) 
  
  end subroutine  

! ***************************************************************************

  !<subroutine>
  subroutine cc_forcesNonStat(rpostprocessing,rvector,rproblem)
  !<description>
  ! A routine that calculated the forces acting on an object
  ! and moves the object according to these forces
  !</description>

  ! structure for a geometry object
  !<inputoutput>
  type(t_problem), intent(INOUT) :: rproblem
  type (t_cc3dpostprocessing),intent(inout) :: rpostprocessing
  !</inputoutput>  

  !<input>
  type(t_vectorBlock), intent(IN) :: rvector
  !</input>
  
  !</subroutine>

  ! Local variables
  ! pointer to the entries of the alpha vector  
  real(DP), dimension(:), pointer :: p_Dvector  
    real(DP), dimension(:), pointer :: p_Dvector1  

  ! pointer to the nodes of the grid
  real(DP), dimension(:,:), pointer :: p_Ddata
  
  ! Forces
  real(DP) :: DResForceX,DResForceY,DResForceZ

  ! pointer to the triangulation structure
  type(t_triangulation), pointer :: p_rtriangulation
  
  type(t_vectorblock) :: rvectorAlpha

  integer :: ive,NEL

  real(DP) :: mytime,df1,df2,dCenterX,dCenterY,dCenterZ,dradius,ddist
  
  !   
  real(dp), dimension(:), pointer :: p_DforceX
  real(dp), dimension(:), pointer :: p_DforceY
  real(dp), dimension(:), pointer :: p_DforceZ
  
  !mytime = dtime

  if (rpostprocessing%rvectorScalarFB%NEQ .ne. 0) &
    call lsyssc_releaseVector (rpostprocessing%rvectorScalarFB)
 
  call lsyssc_createVecByDiscr(rvector%RvectorBlock(1)%p_rspatialDiscr, &
  rpostprocessing%rvectorScalarFB,.true.)                 
  
  ! get a pointer to the triangulation
  p_rtriangulation => &
  rpostprocessing%rvectorScalarFB%p_rspatialDiscr%p_rtriangulation
  
  ! make an L2 Projection of the analytic function to a FEM vector
  call anprj_discrDirect (rpostprocessing%rvectorScalarFB,cc_CBallNonStat,&
                          rproblem%rcollection,iorder=1)  
  
  ! get a pointer to the entries of the vector
  call lsyssc_getbase_double(rpostprocessing%rvectorScalarFB,p_Dvector)
  
  ! get a pointer to the entries of the vector
  call lsyssc_getbase_double(rpostprocessing%rvectorScalarFBQ1,p_Dvector1)
  
  ! And get the vertice coordinates
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, p_Ddata)  

  ! Definition of the circle
  dCenterX = rproblem%rcollection%DQuickaccess(7)
  dCenterY = rproblem%rcollection%DQuickaccess(8)
  dCenterZ = rproblem%rcollection%DQuickaccess(9)
  dradius  = 0.1001

  do ive=1,p_rtriangulation%NVT
    ! Get the distance to the center
    ddist = sqrt( (p_Ddata(1,ive)-dCenterX)**2 + (p_Ddata(2,ive)-dCenterY)**2 &
            + (p_Ddata(3,ive)-dCenterZ)**2 )
    ! Point inside?
    if (ddist .le. dradius) then
      p_Dvector1(ive)=1.0_dp
    else
      p_Dvector1(ive)=0.0_dp  
    end if
  end do

  print *,dCenterX,dCenterY,dCenterZ
  

  call output_lbrk ()
  call output_separator(OU_SEP_EQUAL)
  call output_line ('Q1 Vector recalculated ')
  call output_separator(OU_SEP_EQUAL)   
  
  df1=0.001_dp
  df2=0.041_dp*(0.45_dp*4.0_dp/9.0_dp)**2
  call cc_forcesIntegrationNonStat(rvector, rpostprocessing,rpostprocessing%rvectorScalarFB,DResForceX,&
  DResForceY,DResForceZ,rvectorAlpha,df1,df2)
  
  call output_lbrk()
  call output_line ('Drag Coefficient')
  call output_line ('-----------') 
  call output_line (trim(sys_sdL(DResForceX,10))//"/"//trim(sys_sdL(DResForceY,10))&
                    //"/"//trim(sys_sdL(DResForceZ,10)) ) 
                    
  call lsyssc_getbase_double (rpostprocessing%rResForceX,p_DforceX)
  call lsyssc_getbase_double (rpostprocessing%rResForceY,p_DforceY)
  call lsyssc_getbase_double (rpostprocessing%rResForceZ,p_DforceZ)

  call output_lbrk()
  call output_line ('Pressure Contribution')
  call output_line ('-----------') 
  call output_line (trim(sys_sdL(p_DforceX(2),10))//"/"//trim(sys_sdL(p_DforceY(2),10))&
                    //"/"//trim(sys_sdL(p_DforceZ(2),10)) ) 

  call output_line ('Velocity Contribution')
  call output_line ('-----------') 
  call output_line (trim(sys_sdL(p_DforceX(3),10))//"/"//trim(sys_sdL(p_DforceY(3),10))&
                    //"/"//trim(sys_sdL(p_DforceZ(3),10)) ) 
                    
  call output_line ('Drag Forces')
  call output_line ('-----------') 
  call output_line (trim(sys_sdL(p_DforceX(4),10))//"/"//trim(sys_sdL(p_DforceY(4),10))&
                    //"/"//trim(sys_sdL(p_DforceZ(4),10)) ) 
                    
  
  end subroutine  
  
! ***************************************************************************

    
!<subroutine>
  subroutine cc_CB(cderivative,rdiscretisation, &
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
  integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTest

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
  real :: dCenterX,dCenterY,dDist,eps1
  integer :: i,j,k
  
  k=0
  eps1=-1.0_dp-5
  
  select case (cderivative)
  case (DER_FUNC)
  ! 
  dCenterX = 0.5_DP
  dCenterY = 0.2_DP

  do i=1,nelements
    do j=1,npointsPerElement
        dDist = sqrt( (Dpoints(1,j,i) - dCenterX)**2 + (Dpoints(2,j,i)-dCenterY)**2)
      if(dDist .le. 0.05001_dp)then
        Dvalues(j,i) =  1.0_DP 
        k=k+1
      else
        Dvalues(j,i) = 0.0_DP
      end if
    end do
  end do
  
  rcollection%IquickAccess(1) = &
  rcollection%IquickAccess(1) + k
  
    
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
  
! *************************************************************************** 
    
!<subroutine>
  subroutine cc_CBallNonStat(cderivative,rdiscretisation, &
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
  integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTest

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
  real :: dCenterX,dCenterY,dCenterZ,dRadius,dtime,dtimestep
  integer :: i,j,k
  
  dCenterX = rcollection%DQuickaccess(7)
  dCenterY = rcollection%DQuickaccess(8)
  dCenterZ = rcollection%DQuickaccess(9)
  
  ! Yep, now do some stuff....
  if(present(rcollection)) then
    dtimestep = rcollection%Dquickaccess(4)
    dtime     = rcollection%Dquickaccess(1)
  else
    dtimestep = 0.0_dp
    dtime     = 0.0_dp
  end if

  select case (cderivative)
  case (DER_FUNC)
  ! 
  k=0
  do i=1,nelements
    do j=1,npointsPerElement
        dRadius = sqrt( (Dpoints(1,j,i) - dCenterX)**2 + &
                  (Dpoints(2,j,i)-dCenterY)**2 + (Dpoints(3,j,i)-dCenterZ)**2)
      if(dRadius .le. 0.05_DP)then
        Dvalues(j,i) =  1.0_DP 
        k=k+1
      else
        Dvalues(j,i) = 0.0_DP
      end if
    end do
  end do
  
  rcollection%IquickAccess(1) = &
  rcollection%IquickAccess(1) + k
  
    
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

! ***************************************************************************

!<subroutine>
  subroutine cc_forcesIntegration(rvectorSol,rvectorAlpha,Dfx,Dfy,Dfz,Dalpha,df1,df2)
!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

  ! The body forces are defined as the integrals
  !
  !    Dforces(1) = 2/df2 * int_s [df1 dut/dn n_y - p n_x] ds 
  !    Dforces(2) = 2/df2 * int_s [df1 dut/dn n_x + p n_y] ds 

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorBlock), intent(in), target :: rvectorSol
  
  type(t_vectorBlock), intent(in), target :: Dalpha
  
  type(t_vectorScalar), intent(in), target :: rvectorAlpha

  
  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(IN), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(IN), optional      :: df2
  
!</input>

!<output>
  ! Array receiving the calculated error.
  real(DP),intent(OUT) :: Dfx
  real(DP),intent(OUT) :: Dfy
  real(DP),intent(OUT) :: Dfz
!</output>
  
!</subroutine>

    ! local variables
    real(DP), dimension(:),pointer :: p_Ddata1
    integer :: i,k,icurrentElementDistr, ICUBP, NVE
    integer(I32) :: IEL, IELmax, IELset
    real(DP) :: OM, DN1, DN2, DN3,dpp
    real(DP) :: ah1,ah2,ah3,du1x,du1y,du1z
    real(DP) :: du2x,du2y,du2z,du3x,du3y,du3z,dalx,daly,dalz
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega
    
    ! number of cubature points on the reference element
    integer :: ncubp
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial,indofFunc1,indofFunc2
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer(I32), dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef
    
    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistributionU
    
    type(t_elementDistribution), pointer :: p_relementDistributionP

    type(t_elementDistribution), pointer :: p_relementDistributionA
    
    ! Number of elements in the current element distribution
    integer(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    real(DP) :: negGrad
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_evalElementSet) :: rintSubset
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsTrial
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsFunc1
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsFunc2
    
  
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag
    
    real(dp) :: dpf1,dpf2
    
    ! Prepare the weighting coefficients
    dpf1 = 1.0_dp
    dpf2 = 2.0_dp
    if(present(df1)) dpf1 = df1
    if(present(df2)) dpf2 = df2
    
    !call lsyssc_getbase_double(rvectorAlpha,p_Ddata1)
    !print *,p_ddata1
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rvectorAlpha%p_rspatialDiscr%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)
    
    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.
    do icurrentElementDistr = 1,rvectorSol%p_rblockDiscr%RspatialDiscr(1)%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistributionU => &
      rvectorSol%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(icurrentElementDistr)
      
      p_relementDistributionA =>&
      rvectorAlpha%p_rspatialDiscr%RelementDistr(icurrentElementDistr)
      
      p_relementDistributionP =>&
      rvectorSol%p_rblockDiscr%RspatialDiscr(4)%RelementDistr(icurrentElementDistr)
      
    ! Cancel if this element distribution is empty
    if(p_relementDistributionU%NEL .eq. 0) cycle    
    
    ! get the number of local DOF's for trial functions
    indofTrial = elem_igetNDofLoc(p_relementDistributionU%celement)
    
    indofFunc1 = elem_igetNDofLoc(p_relementDistributionA%celement)
    
    indofFunc2 = elem_igetNDofLoc(p_relementDistributionP%celement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_relementDistributionU%celement)
    
    if (NVE .NE. elem_igetNVE(p_relementDistributionA%celement)) then
      print *,'cc_forcesIntegration: element spaces incompatible!'
      call sys_halt()
    end if
    
    ! get form the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistributionU%celement)
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    ! Now Dxi stores the point coordinates of the cubature points on the reference element    
    call cub_getCubPoints(p_relementDistributionU%ccubTypeEval,ncubp,Dxi,Domega)
    
    ! allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))
    
    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,ncubp
      do k=1,ubound(p_DcubPtsRef,1)
        p_DcubPtsRef(k,i) = Dxi(i,k)
      end do ! end k
    end do ! end i
    
    ! allocate memory for the DOF's of all the elements
    allocate(IdofsTrial(indofTrial,nelementsPerBlock))
    
    allocate(IdofsFunc1(indofFunc1,nelementsPerBlock))
    
    allocate(IdofsFunc2(indofFunc2,nelementsPerBlock))
    
    ! allocate memory for the coefficients
    allocate(Dcoefficients(ncubp, nelementsPerBlock,13))
    
    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = elem_getEvaluationTag(p_relementDistributionU%celement)
    
    cevaluationTag=ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionP%celement))
    
    cevaluationTag=ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionA%celement))
    
    ! Make sure that we have determinants
    cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial functions    
    call storage_getbase_int(p_relementDistributionU%h_IelementList, &
                            p_IelementList)
                            
    ! Get the number of elements there
    NEL = p_relementDistributionU%NEL 
    
    ! Initialize the forces of the element with zero
    Dfx = 0.0_dp
    Dfy = 0.0_dp
    Dfz = 0.0_dp
    
    ! Prepare the call to the evaluation routine of the analytic function.    
    call elprep_init(rintSubset)

    ! Loop over the elements - blockwise.
    do IELset = 1, NEL, PPERR_NELEMSIM
    
      ! We always handle LINF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
      ! Calculate the global DOF's into IdofsTrial.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our LINF_NELEMSIM elements simultaneously.
      IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
      ! Calculate the global DOF's into IdofsTrial
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our LINF_NELEMSIM elements simultaneously.
      !--------------------------------------------------------------------------------            
      call dof_locGlobMapping_mult(rvectorSol%p_rblockDiscr%RspatialDiscr(1),&
                                   p_IelementList(IELset:IELmax),IdofsTrial)
      !-------------------------------------------------------------------------------                                               
      call dof_locGlobMapping_mult(rvectorAlpha%p_rspatialDiscr,&
                                   p_IelementList(IELset:IELmax),IdofsFunc1)
      !--------------------------------------------------------------------------------                                                                                  
      call dof_locGlobMapping_mult(rvectorSol%p_rblockDiscr%RspatialDiscr(4),&
                                   p_IelementList(IELset:IELmax),IdofsFunc2)
      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.                                   
      call elprep_prepareSetForEvaluation(rintSubset,&
           cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax),&
           ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => rintSubset%p_Ddetj
      
      ! In the next loop, we don't have to evaluate the coordinates
      ! on the reference elements anymore.      
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      
      ! At this point, we must select the correct domain integration and coefficient
      ! calculation routine, depending which type of error we should compute!
      
      !----------------------------------------------------------------------------
      !                         EVALUATION PHASE
      !----------------------------------------------------------------------------
       ! We need to build the following system:
      !    /                                                              \
      !   | |-p(x_i)             |   |du1/dx (x_i) du1/dy (x_i) ...|       |   | -dalpha/dx (x_i) |
      !   | |       -p(x_i)      | + |...  you know this works  ...| + u^t | * | -dalpha/dy (x_i) |
      !   | |             -p(x_i)|   |...                       ...|       |   | -dalpha/dz (x_i) |
      !    \                                                              /
      !
      ! 
      ! Get the pressure in the cubature points
      ! Save the result to Dcoefficients(:,:,1)
      
      ! Build the p matrix
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(4), rintSubset,&
           p_relementDistributionP%celement, &
           IdofsFunc2, DER_FUNC, Dcoefficients(:,1:IELmax-IELset+1_I32,1))
      
      ! build the first row of this (u1,u2,u3)
      ! First Row
      ! Save the result to Dcoefficients(:,:,2:4)                        
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(1), rintSubset, &
         p_relementDistributionU%celement, &
         IdofsTrial, DER_DERIV3D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,2))
       
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(1), rintSubset, &
         p_relementDistributionU%celement, &
         IdofsTrial, DER_DERIV3D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,3))
          
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(1), rintSubset, &
         p_relementDistributionU%celement, &
         IdofsTrial, DER_DERIV3D_Z, Dcoefficients(:,1:IELmax-IELset+1_I32,4))
          
      ! Second Row -------------------------------
      ! Save the result to Dcoefficients(:,:,5:7)
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,5))  

      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,6))                
              
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Z, Dcoefficients(:,1:IELmax-IELset+1_I32,7))
              
      ! Third Row -------------------------------
      ! Save the result to Dcoefficients(:,:,8:10)
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(3), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,8))  

      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(3), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,9))                
              
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(3), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Z, Dcoefficients(:,1:IELmax-IELset+1_I32,10))
              
      ! Build the alpha vector
      ! Save the result to Dcoefficients(:,:,6:7)
      call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
              p_relementDistributionA%celement, IdofsFunc1, DER_DERIV3D_X,&
              Dcoefficients(:,1:IELmax-IELset+1_I32,11))
              
      call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
              p_relementDistributionA%celement, IdofsFunc1, DER_DERIV3D_Y,&
              Dcoefficients(:,1:IELmax-IELset+1_I32,12))                                                                                  
              
      call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
              p_relementDistributionA%celement, IdofsFunc1, DER_DERIV3D_Z,&
              Dcoefficients(:,1:IELmax-IELset+1_I32,13))                                                                                                
     
      ! Loop through elements in the set and for each element,
      ! loop through the DOF's and cubature points to calculate the
      ! integral: int_Omega (-p * I + Dj(u)) * (-grad(alpha)) dx
      do IEL=1,IELmax-IELset+1
      
        ! loop over all cubature points on the current element
        do icubp = 1, ncubp
        
          OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
          
          ! get the pressure
          dpp = Dcoefficients(icubp,iel,1)
          
          ! x,y and z derivative of u1
          du1x = Dcoefficients(icubp,iel,2)
          du1y = Dcoefficients(icubp,iel,3)
          du1z = Dcoefficients(icubp,iel,4)
          
          ! x,y and z derivative of u2
          du2x = Dcoefficients(icubp,iel,5)
          du2y = Dcoefficients(icubp,iel,6)
          du2z = Dcoefficients(icubp,iel,7)
          
          ! x,y and z derivative of u3
          du3x = Dcoefficients(icubp,iel,8)
          du3y = Dcoefficients(icubp,iel,9)
          du3z = Dcoefficients(icubp,iel,10)
          
          dalx = Dcoefficients(icubp,iel,11)
          daly = Dcoefficients(icubp,iel,12)                    
          dalz = Dcoefficients(icubp,iel,13)
          
          dn1  = -dalx
          dn2  = -daly          
          dn3  = -dalz                    
          
          ah1 = -dpp*dn1+dpf1*(du1x*dn1+du1y*dn2+du1z*dn3)
          ah2 = -dpp*dn2+dpf1*(du2x*dn1+du2y*dn2+du2z*dn3)          
          ah3 = -dpp*dn3+dpf1*(du3x*dn1+du3y*dn2+du3z*dn3) 
          
!          ah1 = -dpp*dn1+dpf1*((du1x+du1x)*dn1+(du1y+du2x)*dn2 &
!                               +(du1z+du3x)*dn3)
!                               
!          ah2 = -dpp*dn2+dpf1*((du2x+du1y)*dn1+(du2y+du2y)*dn2 &
!                               +(du2z+du3y)*dn3)          
!                               
!          ah3 = -dpp*dn3+dpf1*((du3x+du1z)*dn1+(du3y+du2z)*dn2 &
!                               + (du3z+du3z)*dn3)           
          
          Dfx = Dfx + ah1 * om         
          Dfy = Dfy + ah2 * om                             
          Dfz = Dfz + ah3 * om                                       
        
        end do ! end icubp
      
      end do ! end iel     
      
     end do ! end IELSet
     
     Dfx = Dfx * 2.0_dp/dpf2
     Dfy = Dfy * 2.0_dp/dpf2
     Dfz = Dfz * 2.0_dp/dpf2
      
    ! Release memory
    call elprep_releaseElementSet(rintSubset)

    deallocate(p_DcubPtsRef)
    deallocate(Dcoefficients)
    deallocate(IdofsTrial)
    deallocate(IdofsFunc1)
    deallocate(IdofsFunc2)
     
    
    end do ! end icurrentElementDist
  
  
  end subroutine  
  
! ***************************************************************************  

!<subroutine>
  subroutine cc_forcesIntegrationNonStat(rvectorSol,rpostprocessing,&
             rvectorAlpha,dfx,dfy,dfz,Dalpha,df1,df2)
!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

  ! The body forces are defined as the integrals
  !
  !    Dforces(1) = 2/df2 * int_s [df1 dut/dn n_y - p n_x] ds 
  !    Dforces(2) = 2/df2 * int_s [df1 dut/dn n_x + p n_y] ds 

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorBlock), intent(in), target :: rvectorSol
  
  type(t_vectorBlock), intent(in), target :: Dalpha
  
  type(t_vectorScalar), intent(in), target :: rvectorAlpha
  
  type (t_cc3dpostprocessing),intent(inout) :: rpostprocessing

  
  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(IN), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(IN), optional      :: df2
  
!</input>

!<output>
  ! Array receiving the calculated error.
  real(DP),intent(OUT) :: dfx
  real(DP),intent(OUT) :: dfy
  real(DP),intent(OUT) :: dfz
!</output>
  
!</subroutine>

    ! local variables
    integer :: i,k,icurrentElementDistr, ICUBP, NVE
    integer(I32) :: IEL, IELmax, IELset
    real(DP) :: OM, DN1, DN2, DN3,dpp
    real(DP) :: ah1,ah2,ah3,du1x,du1y,du1z,ah4,ah5,ah6,ah7,ah8,ah9
    real(DP) :: du2x,du2y,du2z,du3x,du3y,du3z,dalx,daly,dalz
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega
    
    ! number of cubature points on the reference element
    integer :: ncubp
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial,indofFunc1,indofFunc2
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer(I32), dimension(:), pointer :: p_IelementList
    
    ! 
    real(dp), dimension(:), pointer :: p_DforceX
    real(dp), dimension(:), pointer :: p_DforceY
    real(dp), dimension(:), pointer :: p_DforceZ
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef
    
    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistributionU
    
    type(t_elementDistribution), pointer :: p_relementDistributionP

    type(t_elementDistribution), pointer :: p_relementDistributionA
    
    ! Number of elements in the current element distribution
    integer(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    real(DP) :: negGrad
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_evalElementSet) :: rintSubset
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsTrial
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsFunc1
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsFunc2
    
  
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag
    
    real(dp) :: dpf1,dpf2
    
    ! Prepare the weighting coefficients
    dpf1 = 1.0_dp
    dpf2 = 2.0_dp
    if(present(df1)) dpf1 = df1
    if(present(df2)) dpf2 = df2
    
    !
    call lsyssc_getbase_double (rpostprocessing%rResForceX,p_DforceX)
    call lsyssc_getbase_double (rpostprocessing%rResForceY,p_DforceY)
    call lsyssc_getbase_double (rpostprocessing%rResForceZ,p_DforceZ)
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rvectorAlpha%p_rspatialDiscr%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)
    
    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.
    do icurrentElementDistr = 1,rvectorSol%p_rblockDiscr%RspatialDiscr(1)%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistributionU => &
      rvectorSol%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(icurrentElementDistr)
      
      p_relementDistributionA =>&
      rvectorAlpha%p_rspatialDiscr%RelementDistr(icurrentElementDistr)
      
      p_relementDistributionP =>&
      rvectorSol%p_rblockDiscr%RspatialDiscr(4)%RelementDistr(icurrentElementDistr)
      
    ! Cancel if this element distribution is empty
    if(p_relementDistributionU%NEL .eq. 0) cycle    
    
    ! get the number of local DOF's for trial functions
    indofTrial = elem_igetNDofLoc(p_relementDistributionU%celement)
    
    indofFunc1 = elem_igetNDofLoc(p_relementDistributionA%celement)
    
    indofFunc2 = elem_igetNDofLoc(p_relementDistributionP%celement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_relementDistributionU%celement)
    
    if (NVE .NE. elem_igetNVE(p_relementDistributionA%celement)) then
      print *,'cc_forcesIntegration: element spaces incompatible!'
      call sys_halt()
    end if
    
    ! get form the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistributionU%celement)
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    ! Now Dxi stores the point coordinates of the cubature points on the reference element    
    call cub_getCubPoints(p_relementDistributionU%ccubTypeEval,ncubp,Dxi,Domega)
    
    ! allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))
    
    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,ncubp
      do k=1,ubound(p_DcubPtsRef,1)
        p_DcubPtsRef(k,i) = Dxi(i,k)
      end do ! end k
    end do ! end i
    
    ! allocate memory for the DOF's of all the elements
    allocate(IdofsTrial(indofTrial,nelementsPerBlock))
    
    allocate(IdofsFunc1(indofFunc1,nelementsPerBlock))
    
    allocate(IdofsFunc2(indofFunc2,nelementsPerBlock))
    
    ! allocate memory for the coefficients
    allocate(Dcoefficients(ncubp, nelementsPerBlock,13))
    
    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = elem_getEvaluationTag(p_relementDistributionU%celement)
    
    cevaluationTag=ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionP%celement))
    
    cevaluationTag=ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionA%celement))
    
    ! Make sure that we have determinants
    cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial functions    
    call storage_getbase_int(p_relementDistributionU%h_IelementList, &
                            p_IelementList)
                            
    ! Get the number of elements there
    NEL = p_relementDistributionU%NEL 
    
    ! Initialize the forces of the element with zero
    Dfx = 0.0_dp
    Dfy = 0.0_dp
    Dfz = 0.0_dp
    p_DforceX(1)= 0.0_dp
    p_DforceY(1)= 0.0_dp
    p_DforceZ(1)= 0.0_dp
    p_DforceX(4)= 0.0_dp
    p_DforceY(4)= 0.0_dp
    p_DforceZ(4)= 0.0_dp

    
    ! Prepare the call to the evaluation routine of the analytic function.    
    call elprep_init(rintSubset)

    ! Loop over the elements - blockwise.
    do IELset = 1, NEL, PPERR_NELEMSIM
    
      ! We always handle LINF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
      ! Calculate the global DOF's into IdofsTrial.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our LINF_NELEMSIM elements simultaneously.
      IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
      ! Calculate the global DOF's into IdofsTrial
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our LINF_NELEMSIM elements simultaneously.
      !--------------------------------------------------------------------------------            
      call dof_locGlobMapping_mult(rvectorSol%p_rblockDiscr%RspatialDiscr(1),&
                                   p_IelementList(IELset:IELmax),IdofsTrial)
      !-------------------------------------------------------------------------------                                               
      call dof_locGlobMapping_mult(rvectorAlpha%p_rspatialDiscr,&
                                   p_IelementList(IELset:IELmax),IdofsFunc1)
      !--------------------------------------------------------------------------------                                                                                  
      call dof_locGlobMapping_mult(rvectorSol%p_rblockDiscr%RspatialDiscr(4),&
                                   p_IelementList(IELset:IELmax),IdofsFunc2)
      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.                                   
      call elprep_prepareSetForEvaluation(rintSubset,&
           cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax),&
           ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => rintSubset%p_Ddetj
      
      ! In the next loop, we don't have to evaluate the coordinates
      ! on the reference elements anymore.      
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      
      ! At this point, we must select the correct domain integration and coefficient
      ! calculation routine, depending which type of error we should compute!
      
      !----------------------------------------------------------------------------
      !                         EVALUATION PHASE
      !----------------------------------------------------------------------------
       ! We need to build the following system:
      !    /                                                              \
      !   | |-p(x_i)             |   |du1/dx (x_i) du1/dy (x_i) ...|       |   | -dalpha/dx (x_i) |
      !   | |       -p(x_i)      | + |...  you know this works  ...| + u^t | * | -dalpha/dy (x_i) |
      !   | |             -p(x_i)|   |...                       ...|       |   | -dalpha/dz (x_i) |
      !    \                                                              /
      !
      ! 
      ! Get the pressure in the cubature points
      ! Save the result to Dcoefficients(:,:,1)
      
      ! Build the p matrix
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(4), rintSubset,&
           p_relementDistributionP%celement, &
           IdofsFunc2, DER_FUNC, Dcoefficients(:,1:IELmax-IELset+1_I32,1))
      
      ! build the first row of this (u1,u2,u3)
      ! First Row
      ! Save the result to Dcoefficients(:,:,2:4)                        
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(1), rintSubset, &
         p_relementDistributionU%celement, &
         IdofsTrial, DER_DERIV3D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,2))
       
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(1), rintSubset, &
         p_relementDistributionU%celement, &
         IdofsTrial, DER_DERIV3D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,3))
          
      call fevl_evaluate_sim3(rvectorSol%RvectorBlock(1), rintSubset, &
         p_relementDistributionU%celement, &
         IdofsTrial, DER_DERIV3D_Z, Dcoefficients(:,1:IELmax-IELset+1_I32,4))
          
      ! Second Row -------------------------------
      ! Save the result to Dcoefficients(:,:,5:7)
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,5))  

      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,6))                
              
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Z, Dcoefficients(:,1:IELmax-IELset+1_I32,7))
              
      ! Third Row -------------------------------
      ! Save the result to Dcoefficients(:,:,8:10)
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(3), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,8))  

      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(3), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,9))                
              
      call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(3), rintSubset, &
              p_relementDistributionU%celement, &
              IdofsTrial, DER_DERIV3D_Z, Dcoefficients(:,1:IELmax-IELset+1_I32,10))
              
      ! Build the alpha vector
      ! Save the result to Dcoefficients(:,:,6:7)
      call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
              p_relementDistributionA%celement, IdofsFunc1, DER_DERIV3D_X,&
              Dcoefficients(:,1:IELmax-IELset+1_I32,11))
              
      call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
              p_relementDistributionA%celement, IdofsFunc1, DER_DERIV3D_Y,&
              Dcoefficients(:,1:IELmax-IELset+1_I32,12))                                                                                  
              
      call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
              p_relementDistributionA%celement, IdofsFunc1, DER_DERIV3D_Z,&
              Dcoefficients(:,1:IELmax-IELset+1_I32,13))                                                                                                
     
      ! Loop through elements in the set and for each element,
      ! loop through the DOF's and cubature points to calculate the
      ! integral: int_Omega (-p * I + Dj(u)) * (-grad(alpha)) dx
      do IEL=1,IELmax-IELset+1
      
        ! loop over all cubature points on the current element
        do icubp = 1, ncubp
        
          OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
          
          ! get the pressure
          dpp = Dcoefficients(icubp,iel,1)
          
          ! x,y and z derivative of u1
          du1x = Dcoefficients(icubp,iel,2)
          du1y = Dcoefficients(icubp,iel,3)
          du1z = Dcoefficients(icubp,iel,4)
          
          ! x,y and z derivative of u2
          du2x = Dcoefficients(icubp,iel,5)
          du2y = Dcoefficients(icubp,iel,6)
          du2z = Dcoefficients(icubp,iel,7)
          
          ! x,y and z derivative of u3
          du3x = Dcoefficients(icubp,iel,8)
          du3y = Dcoefficients(icubp,iel,9)
          du3z = Dcoefficients(icubp,iel,10)
          
          dalx = Dcoefficients(icubp,iel,11)
          daly = Dcoefficients(icubp,iel,12)                    
          dalz = Dcoefficients(icubp,iel,13)
          
          dn1  = -dalx
          dn2  = -daly          
          dn3  = -dalz                    
          
          ah1 = -dpp*dn1+dpf1*(du1x*dn1+du1y*dn2+du1z*dn3)
          ah2 = -dpp*dn2+dpf1*(du2x*dn1+du2y*dn2+du2z*dn3)          
          ah3 = -dpp*dn3+dpf1*(du3x*dn1+du3y*dn2+du3z*dn3) 

          ah4 = -dpp*dn1
          ah5 = -dpp*dn2
          ah6 = -dpp*dn3
          
          ah7 = dpf1*(du1x*dn1+du1y*dn2+du1z*dn3)
          ah8 = dpf1*(du2x*dn1+du2y*dn2+du2z*dn3)
          ah9 = dpf1*(du3x*dn1+du3y*dn2+du3z*dn3)


          p_DforceX(2) = p_DforceX(2) + ah4 * om
          p_DforceY(2) = p_DforceY(2) + ah5 * om
          p_DforceZ(2) = p_DforceZ(2) + ah6 * om          
          
          p_DforceX(3) = p_DforceX(3) + ah7 * om
          p_DforceY(3) = p_DforceY(3) + ah8 * om
          p_DforceZ(3) = p_DforceZ(3) + ah9 * om          
          

          p_DforceX(4) = p_DforceX(4) + (ah7+ah4) * om
          p_DforceY(4) = p_DforceY(4) + (ah8+ah5) * om
          p_DforceZ(4) = p_DforceZ(4) + (ah9+ah6) * om          

          
!          ah1 = -dpp*dn1+dpf1*((du1x+du1x)*dn1+(du1y+du2x)*dn2 &
!                               +(du1z+du3x)*dn3)
!                               
!          ah2 = -dpp*dn2+dpf1*((du2x+du1y)*dn1+(du2y+du2y)*dn2 &
!                               +(du2z+du3y)*dn3)          
!                               
!          ah3 = -dpp*dn3+dpf1*((du3x+du1z)*dn1+(du3y+du2z)*dn2 &
!                               + (du3z+du3z)*dn3)           
          
          Dfx = Dfx + ah1 * om         
          Dfy = Dfy + ah2 * om                             
          Dfz = Dfz + ah3 * om                                       
        
        end do ! end icubp
      
      end do ! end iel     
      
     end do ! end IELSet
     
     p_DforceX(1) = Dfx
     p_DforceY(1) = Dfy
     p_DforceZ(1) = Dfz
     
     Dfx = Dfx * 2.0_dp/dpf2
     Dfy = Dfy * 2.0_dp/dpf2
     Dfz = Dfz * 2.0_dp/dpf2
      
    ! Release memory
    call elprep_releaseElementSet(rintSubset)

    deallocate(p_DcubPtsRef)
    deallocate(Dcoefficients)
    deallocate(IdofsTrial)
    deallocate(IdofsFunc1)
    deallocate(IdofsFunc2)
     
    
    end do ! end icurrentElementDist
  
  
  end subroutine  

! ***************************************************************************

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
  ! Postprocvessing structure.
  type(t_cc3dpostprocessing), intent(INOUT) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors which might be allocated.
    if (rpostprocessing%rvectorVelX%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorVelX)
    if (rpostprocessing%rvectorVelY%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorVelY)
    if (rpostprocessing%rvectorVelZ%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorVelZ)
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
    if (rpostprocessing%rvectorScalarFB%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorScalarFB)
    if (rpostprocessing%rResForceX%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rResForceX)
    if (rpostprocessing%rResForceY%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rResForceY)
    if (rpostprocessing%rResForceZ%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rResForceZ)
    if (rpostprocessing%rvectorScalarFBQ1%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorScalarFBQ1)
      

  end subroutine

  !****************************************************************************
  
  
!<subroutine>
  subroutine cc_updateParticlePosition(rpostprocessing,rproblem,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  
  type (t_cc3dpostprocessing),intent(inout) :: rpostprocessing
  type(t_problem), intent(INOUT) :: rproblem
  real(dp), intent(in) :: dtimestep
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,k,icurrentElementDistr, ICUBP, NVE
    integer(I32) :: IEL, IELmax, IELset
    real(DP) :: OM, DN1, DN2, DN3,dpp
    real(DP) :: ah1,ah2,ah3,dvelx,dvely,dvelz,dmasssl,ddmasssl,dvolume
    real(DP) :: dfx,dfy,dfz,du3x,du3y,du3z,dalx,daly,dalz,nennerX
    real(DP) :: dCenterX,dCenterY,dCenterZ,dCenterXold,dCenterYold,dCenterZold
    
    ! 
    real(dp), dimension(:), pointer :: p_DforceX
    real(dp), dimension(:), pointer :: p_DforceY
    real(dp), dimension(:), pointer :: p_DforceZ
    
    !
    call lsyssc_getbase_double (rpostprocessing%rResForceX,p_DforceX)
    call lsyssc_getbase_double (rpostprocessing%rResForceY,p_DforceY)
    call lsyssc_getbase_double (rpostprocessing%rResForceZ,p_DforceZ)
        
    ! position
    dCenterX=rproblem%rcollection%DQuickaccess(7)
    dCenterY=rproblem%rcollection%DQuickaccess(8)
    dCenterZ=rproblem%rcollection%DQuickaccess(9)
    
    ! old position
    dCenterXold=rproblem%rcollection%DQuickaccess(7)
    dCenterYold=rproblem%rcollection%DQuickaccess(8)
    dCenterZold=rproblem%rcollection%DQuickaccess(9)
    
    ! some physical parameters
    dvolume  = (0.05_dp)**3 * (4.0_dp/3.0_dp) * 3.14159265_dp
    dmasssl  = 5.1_dp * dvolume 
    ddmasssl = 4.1_dp * dvolume 
    
    ! Mean resistance forces for the given time step
    dfx = 0.5_dp * (p_DforceX(5) + p_DforceX(4))
    dfy = 0.5_dp * (p_DforceY(5) + p_DforceY(4))
    dfz = 0.5_dp * (p_DforceZ(5) + p_DforceZ(4))
    
    ! Velocity difference for the given time step
    dvelx   = dtimestep*(1.0_dp*dfx+ddmasssl*0.0_dp)/dmasssl
    dvely   = dtimestep*(1.0_dp*dfy+ddmasssl*0.0_dp)/dmasssl
    dvelz   = dtimestep*(1.0_dp*dfz+ddmasssl*0.0_dp)/dmasssl
    nennerX = dtimestep*(1.0_dp*dfx+ddmasssl*0.0_dp)
    
    print *,"--------------------------"
    
    print *,"nennerX: ",nennerX
 
    print *,"--------------------------"
    
    print *,"dmasssl: ",dmasssl
 
 
    ! save the old forces for the next time step
    p_DforceX(5) = p_DforceX(4)
    p_DforceY(5) = p_DforceY(4)
    p_DforceZ(5) = p_DforceZ(4)
    
    ! Update the position
!    dCenterX=dCenterX+dtime*p_DforceX(4)+0.5_dp*dvelx+0.025_dp !gravity
!    dCenterY=dCenterY+dtime*p_DforceY(4)+0.5_dp*dvely
!    dCenterZ=dCenterZ+dtime*p_DforceZ(4)+0.5_dp*dvelz

    dCenterX=dCenterX+dtimestep*(p_DforceX(6)+0.5_dp*dvelx)
    dCenterY=dCenterY+dtimestep*(p_DforceY(6)+0.5_dp*dvely)
    dCenterZ=dCenterZ+dtimestep*(p_DforceZ(6)+0.5_dp*dvelz)
    
    print *,"New Position X: ",dCenterX
    
    print *,"--------------------------"
    
    print *,"VelocityXfromFlow: ",dvelx

    ! save the current velocity
    p_DforceX(6) = p_DforceX(6) + dvelx
    p_DforceY(6) = p_DforceY(6) + dvely
    p_DforceZ(6) = p_DforceZ(6) + dvelz

!    p_DforceX(6) = dvelx
!    p_DforceY(6) = dvely

    
    ! write back new position
    rproblem%rcollection%DQuickaccess(7)=dCenterX
    rproblem%rcollection%DQuickaccess(8)=dCenterY
    rproblem%rcollection%DQuickaccess(9)=dCenterZ
    
    ! hier ist das problem, diese Velocity ist nicht
    ! sowas wie der VelocityVector...
    rproblem%rcollection%DQuickaccess(10)= p_DforceX(6)
    rproblem%rcollection%DQuickaccess(11)= p_DforceY(6)
    rproblem%rcollection%DQuickaccess(11)= p_DforceZ(6)

!    rproblem%rcollection%DQuickaccess(10)=(dCenterX-dCenterXold)/dtimestep
!    rproblem%rcollection%DQuickaccess(11)=(dCenterY-dCenterYold)/dtimestep
!    rproblem%rcollection%DQuickaccess(12)=(dCenterZ-dCenterZold)/dtimestep

!    rproblem%rcollection%DQuickaccess(10)=(dCenterX-dCenterXold)
!    rproblem%rcollection%DQuickaccess(11)=(dCenterY-dCenterYold)
    
    print *,"--------------------------"
    
    print *,"accumulated velocity: ",p_DforceX(6)
    
    
    print *,"--------------------------"
    
    print *,"VelocityXTotal: ",rproblem%rcollection%DQuickaccess(10)
    
  end subroutine

! ***************************************************************************
  
!<subroutine>

  subroutine cc_donepostprocessing (rpostprocessing)

!<description>
  ! Releases a given problem structure. All allocated memory of this structure
  ! is released.
!</description>

!<inputoutput>  
  type(t_cc3dpostprocessing), intent(INOUT) :: rpostprocessing
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

end module
