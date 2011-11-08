!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2init </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic initialisation routines for cc2dmini_method2:
!# Initialisation of the main structures with parameters:
!#
!# 1.) cc_initParameters
!#     -> Init the problem structure with data from the INI/DAT files
!#
!# 2.) cc_doneParameters
!#     -> Clean up the problem structure
!#
!# 3.) cc_initParamTriang
!#     -> Read parametrisation, read triangulation, refine the mesh
!#
!# 4.) cc_doneParamTriang
!#     -> Remove the meshes of all levels from the heap
!#
!# </purpose>
!##############################################################################

module paramtriainit

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
    
  use basicstructures
  use nonstationaryoptcsolver
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initOutput (rproblem)
  
!<description>
  ! Initialises basic output settings based on the parameters in the DAT file.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initParameters (rproblem)
  
!<description>
  ! Initialises the structure rproblem with data from the initialisation
  ! files.
  !
  ! The parameters in rproblem\%rparameters are evaluated.
  ! Important parameters are written to the problem structure
  ! rproblem and/or the enclosed collection rproblem\%rcollection.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    real(DP) :: dnu,d1
    integer :: ilvmin,ilvmax,i1
    character(LEN=SYS_STRLEN) :: sstring
    character(LEN=SYS_STRLEN), dimension(2) :: SrhsExpressions
    character(LEN=10), dimension(3), parameter :: EXPR_VARIABLES = &
      (/'X    ','Y    ','TIME '/)

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

    ! Get the viscosity parameter, save it to the problem structure
    ! as well as into the collection.
    ! Note that the parameter in the DAT file is 1/nu !
    call parlst_getvalue_double (rproblem%rparamList,'CC-PHYSICSPRIMAL',&
                                 'RE',dnu,1000.0_DP)

    dnu = 1E0_DP/dnu
    rproblem%rphysicsPrimal%dnu = dnu
    
    ! Which type of problem to discretise? (Stokes, Navier-Stokes,...)
    call parlst_getvalue_int (rproblem%rparamList,'CC-PHYSICSPRIMAL',&
                              'iEquation',i1,0)
    rproblem%rphysicsPrimal%cequation = i1

    ! Type of subproblem (gradient tensor, deformation tensor,...)
    call parlst_getvalue_int (rproblem%rparamList,'CC-PHYSICSPRIMAL',&
                              'isubEquation',i1,0)
    rproblem%rphysicsPrimal%isubEquation = i1

    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)
    if (ilvmin .le. 0) ilvmin = ilvmax+ilvmin

    ! Initialise the level in the problem structure
    ilvmin = min(ilvmin,ilvmax)
    ilvmax = max(ilvmin,ilvmax)
    rproblem%NLMIN = ilvmin
    rproblem%NLMAX = ilvmax

    ! Allocate memory for the levels
    allocate(rproblem%RlevelInfo(1:ilvmax))
    
    ! Specification of the RHS
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iRHS',i1,0)
    rproblem%irhs = i1

    ! If the RHS is given as expression, create a parser object for the
    ! expression.
    if (rproblem%irhs .eq. -1) then
      
      ! Create an analytic solution structure based on the
      ! RHS expressions.
      call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
          'srhsExpressionX',SrhsExpressions(1),"", bdequote=.true.)

      call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
          'srhsExpressionY',SrhsExpressions(2),"", bdequote=.true.)
      
      call ansol_init_meshless (rproblem%rrhs)
      call ansol_configAnalytical (rsolution,SrhsExpressions)
      
    else if (rproblem%irhs .eq. 0) then
      
      ! Initialise a zero flow.
      call ansol_init_meshless (rproblem%rrhs)
      call ansol_configAnalytical (rsolution)
    
    end if

    ! Type of boundary conditions
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iBoundary',i1,0)
    call collct_setvalue_int( rproblem%rcollection,'IBOUNDARY',i1,.true.)

    ! Time dependence
    call cc_initParTimeDependence (rproblem,'TIME-DISCRETISATION',&
        rproblem%rparamList)
        
    ! Optimal control
    call cc_initOptControl(rproblem)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneParameters (rproblem)
  
!<description>
  ! Cleans up parameters read from the DAT files. Removes all references to
  ! parameters from the collection rproblem\%rcollection that were
  ! set up in cc_initParameters.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Optimal control
    call cc_doneOptControl(rproblem)

    ! Release memory of the level specification
    deallocate(rproblem%RlevelInfo)

    ! Remove information about boundary conditions
    call collct_deleteValue(rproblem%rcollection,'IBOUNDARY')
    
    ! Remove type of problem to discretise
    call collct_deleteValue(rproblem%rcollection,'ISTOKES')

    ! Release the RHS in case it was initialised.
    call ansol_done(rproblem%rrhs)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initParamTriang (rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the NLMIN/NLMAX parameters
  ! from the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,imeshType,ncellsX
  
    ! Variable for a filename:
    character(LEN=60) :: sPRMFile, sTRIFile

    ! Get the .prm and the .tri file from the parameter list.
    ! note that parlst_getvalue_string returns us exactly what stands
    ! in the parameter file, so we have to apply READ to get rid of
    ! probable ''!
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
        "sParametrisation",sPRMFile,bdequote=.true.)
                              
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
        "sMesh",sTRIFile,bdequote=.true.)
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rproblem%rboundary, sPrmFile)
        
    ! Allocate memory for all the meshes.
    allocate(rproblem%p_rbasicTria(1:rproblem%NLMAX))
        
    ! Now set up the basic triangulation. Which type of triangulation should
    ! be set up?
    call parlst_getvalue_int (rproblem%rparamList,'PARAMTRIANG',&
        'imeshType',imeshType,0)
    select case (imeshType)
    case (0)
      ! Standard mesh, specified by a TRI file.
      call tria_readTriFile2D (rproblem%p_rbasicTria(1), &
          sTRIFile, rproblem%rboundary)
    case (1)
      ! Sliced QUAD mesh with ncellsX cells on the coarse grid.
      call parlst_getvalue_int (rproblem%rparamList,'PARAMTRIANG',&
                                'ncellsX',ncellsX,1)
      call cc_generateSlicedQuadMesh(rproblem%p_rbasicTria(rproblem%NLMIN),&
          ncellsX)
    case default
      call output_line ('Unknown mesh type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
      call sys_halt()
    end select
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (rproblem%p_rbasicTria(1),&
        rproblem%rboundary)

    ! Refine the mesh up to the minimum level.
    ! Skip levels 2..nlmin-1.
    if (rproblem%NLMIN .gt. 1) then
      call tria_duplicate(rproblem%p_rbasicTria(1),&
          rproblem%p_rbasicTria(rproblem%NLMIN),TR_SHARE_ALL)
      call tria_quickRefine2LevelOrdering(rproblem%NLMIN-1,&
          rproblem%p_rbasicTria(rproblem%NLMIN),rproblem%rboundary)
          
      ! Create information about adjacencies and everything one needs from
      ! a triangulation. Afterwards, we have the mesh at level NLMIN.
      call tria_initStandardMeshFromRaw (rproblem%p_rbasicTria(rproblem%NLMIN),&
          rproblem%rboundary)
    end if

    rproblem%RlevelInfo(rproblem%NLMIN)%p_rtriangulation => &
        rproblem%p_rbasicTria(rproblem%NLMIN)
    
    ! Now, refine to level up to nlmax.
    do i=rproblem%NLMIN+1,rproblem%NLMAX
      call tria_refine2LevelOrdering (rproblem%p_rbasicTria(i-1),&
          rproblem%p_rbasicTria(i), rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%p_rbasicTria(i),&
          rproblem%rboundary)
      rproblem%RlevelInfo(i)%p_rtriangulation => rproblem%p_rbasicTria(i)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    ! Release the triangulation on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      call tria_done (rproblem%p_rbasicTria(i))
    end do
    if (rproblem%NLMIN .gt. 1) then
      call tria_done (rproblem%p_rbasicTria(1))
    end if
    deallocate(rproblem%p_rbasicTria)
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)
    
  end subroutine


end module
