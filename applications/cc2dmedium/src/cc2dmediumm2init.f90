!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2init </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic initialisation routines for cc2dmini_method2:
!# Initialisation of the main structures with parameters:
!#
!# 1.) c2d2_initParameters
!#     -> Init the problem structure with data from the INI/DAT files
!#
!# 2.) c2d2_doneParameters
!#     -> Clean up the problem structure
!#
!# 3.) c2d2_initParamTriang
!#     -> Read parametrisation, read triangulation, refine the mesh
!#
!# 4.) c2d2_doneParamTriang
!#     -> Remove the meshes of all levels from the heap
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2init

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
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmediumm2nonstationary
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initOutput (rproblem)
  
!<description>
  ! Initialises basic output settings based on the parameters in the DAT file.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initParameters (rproblem)
  
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
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    REAL(DP) :: dnu,d1
    INTEGER :: ilvmin,ilvmax,i1

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

    ! Get the viscosity parameter, save it to the problem structure
    ! as well as into the collection.
    ! Note that the parameter in the DAT file is 1/nu !
    CALL parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'RE',dnu,1000.0_DP)

    dnu = 1E0_DP/dnu
    rproblem%dnu = dnu
    
    ! By default, X- and Y-velocity matrix are coupled.
    rproblem%bdecoupledXY = .FALSE.
    
    ! Add the (global) viscosity parameter
    CALL collct_setvalue_real(rproblem%rcollection,'NU',dnu,.TRUE.)

    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)

    ! Initialise the level in the problem structure
    rproblem%NLMIN = ilvmin
    rproblem%NLMAX = ilvmax

    CALL collct_setvalue_int (rproblem%rcollection,'NLMIN',ilvmin,.TRUE.)
    CALL collct_setvalue_int (rproblem%rcollection,'NLMAX',ilvmax,.TRUE.)
    
    ! Which type of problem to discretise? (Stokes, Navier-Stokes,...)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iEquation',i1,0)
    rproblem%iequation = i1

    ! Type of subproblem (gradient tensor, deformation tensor,...)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isubEquation',i1,0)
    rproblem%isubEquation = i1

    ! Stabilisation of nonlinearity
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iUpwind',i1,0)
    CALL collct_setvalue_int( rproblem%rcollection,'IUPWIND',i1,.TRUE.)

    CALL parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                'dUpsam',d1,0.0_DP)
    CALL collct_setvalue_real (rproblem%rcollection,'UPSAM',d1,.TRUE.)

    ! Type of boundary conditions
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iBoundary',i1,0)
    CALL collct_setvalue_int( rproblem%rcollection,'IBOUNDARY',i1,.TRUE.)

    ! Time dependence
    CALL c2d2_initParTimeDependence (rproblem,'TIME-DISCRETISATION',&
        rproblem%rparamList)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneParameters (rproblem)
  
!<description>
  ! Cleans up parameters read from the DAT files. Removes all references to
  ! parameters from the collection rproblem\%rcollection that were
  ! set up in c2d2_initParameters.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Remove information about boundary conditions
    CALL collct_deleteValue(rproblem%rcollection,'IBOUNDARY')
    
    ! Remove information about stabilisation
    CALL collct_deleteValue(rproblem%rcollection,'UPSAM')
    CALL collct_deleteValue(rproblem%rcollection,'IUPWIND')

    ! Remove type of problem to discretise
    CALL collct_deleteValue(rproblem%rcollection,'ISTOKES')

    ! Remove min/max level from the collection
    CALL collct_deleteValue(rproblem%rcollection,'NLMAX')
    CALL collct_deleteValue(rproblem%rcollection,'NLMIN')

    ! Remove the viscosity parameter
    CALL collct_deletevalue(rproblem%rcollection,'NU')
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initParamTriang (rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the NLMIN/NLMAX parameters
  ! from the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,ilvmin,ilvmax
  
    ! Variable for a filename:  
    CHARACTER(LEN=SYS_STRLEN) :: sString
    CHARACTER(LEN=SYS_STRLEN) :: sPRMFile, sTRIFile

    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)
    
    ! Get the .prm and the .tri file from the parameter list.
    ! note that parlst_getvalue_string returns us exactly what stands
    ! in the parameter file, so we have to apply READ to get rid of
    ! probable ''!
    CALL parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sParametrisation',sString)
    READ (sString,*) sPRMFile
                              
    CALL parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sMesh',sString)
    READ (sString,*) sTRIFile
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    ! Set p_rboundary to NULL() to create a new structure.
    NULLIFY(rproblem%p_rboundary)
    CALL boundary_read_prm(rproblem%p_rboundary, sPrmFile)
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation, &
        sTRIFile, rproblem%p_rboundary)

    ! Refine the mesh up to the minimum level
    CALL tria_quickRefine2LevelOrdering(rproblem%NLMIN-1,&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%p_rboundary)

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    CALL tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%p_rboundary)
    
    ! Now, refine to level up to nlmax.
    DO i=rproblem%NLMIN+1,rproblem%NLMAX
      CALL tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%p_rboundary)
      CALL tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%p_rboundary)
    END DO
    
    ! Compress the level hierarchy.
    ! Share the vertex coordinates of all levels, so the coarse grid coordinates
    ! are 'contained' in the fine grid coordinates. The effect is:
    ! 1.) Save some memory
    ! 2.) Every change in the fine grid coordinates also affects the coarse
    !     grid coordinates and vice versa.
    DO i=rproblem%NLMAX-1,rproblem%NLMIN
      CALL tria_compress2LevelOrdHierarchy (rproblem%RlevelInfo(i+1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    ! Release the triangulation on all levels
    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      CALL tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
  END SUBROUTINE

END MODULE
