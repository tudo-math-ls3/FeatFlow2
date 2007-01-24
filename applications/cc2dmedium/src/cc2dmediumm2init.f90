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
  
  IMPLICIT NONE
  
CONTAINS

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
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    REAL(DP) :: dnu,d1
    INTEGER :: ilvmin,ilvmax,i1

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
    
    ! Which type of problem to discretise?
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iEquation',i1,0)
    CALL collct_setvalue_int (rproblem%rcollection,'ISTOKES',i1,.TRUE.)

    ! Stabilisation of nonlinearity
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iUpwind',i1,0)
    CALL collct_setvalue_int( rproblem%rcollection,'IUPWIND',i1,.TRUE.)

    CALL parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                'dUpsam',d1,0.0_DP)
    CALL collct_setvalue_real (rproblem%rcollection,'UPSAM',d1,.TRUE.)

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
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

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
  
    INCLUDE 'cout.inc'
    INCLUDE 'cerr.inc'
    INCLUDE 'cmem.inc'
    INCLUDE 'cparametrization.inc'

!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the NLMIN/NLMAX parameters
  ! from the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

  ! local variables
  INTEGER :: i,ilvmin,ilvmax
  
    ! For compatibility to old F77: an array accepting a set of triangulations
    INTEGER, DIMENSION(SZTRIA,NNLEV) :: TRIAS

    ! Variable for a filename:  
    CHARACTER(LEN=256) :: sString
    CHARACTER(LEN=60) :: sPRMFile, sTRIFile

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
        
    ! Remark that this does not read in the parametrisation for FEAT 1.x.
    ! Unfortunately we still need it for creating the initial triangulation!
    ! Therefore, read the file again wihh FEAT 1.x routines.
    IMESH = 1
    CALL GENPAR (.TRUE.,IMESH,sPRMFile)

    ! Now read in the triangulation - in FEAT 1.x syntax.
    ! Refine it to level ilvmin/ilvmax.
    ! This will probably modify ilvmin/ilvmax in case of a level
    ! shift, i.e. if ilvmax > ilvmin+9 !
    ! After this routine, we have to rely on ilvmin/ilvmax in the
    ! problem structure rather than those in the parameters.
    CALL INMTRI (2,TRIAS,rproblem%NLMIN,rproblem%NLMAX,0,0,sTRIfile)
    
    ! ... and create a FEAT 2.0 triangulation for that. Until the point where
    ! we recreate the triangulation routines, this method has to be used
    ! to get a triangulation.
    DO i=rproblem%NLMIN,rproblem%NLMAX
      CALL tria_wrp_tria2Structure(TRIAS(:,i),rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! The TRIAS(,)-array is now part pf the triangulation structure,
    ! we don't need it anymore.
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    ! For compatibility to old F77: an array accepting a set of triangulations
    INTEGER, DIMENSION(SZTRIA,NNLEV) :: TRIAS

    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Release the old FEAT 1.x handles.
      ! Get the old triangulation structure of level ilv from the
      ! FEAT2.0 triangulation:
      TRIAS(:,i) = rproblem%RlevelInfo(i)%rtriangulation%Itria
      CALL DNMTRI (i,i,TRIAS)
      
      ! then the FEAT 2.0 stuff...
      CALL tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
    ! Don't forget to throw away the old FEAT 1.0 boundary definition!
    CALL DISPAR

  END SUBROUTINE

END MODULE
