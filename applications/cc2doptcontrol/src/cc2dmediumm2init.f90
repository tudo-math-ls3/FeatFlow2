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

  SUBROUTINE cc_initOutput (rproblem)
  
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

  SUBROUTINE cc_initParameters (rproblem)
  
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
    IF (ilvmin .LE. 0) ilvmin = ilvmax+ilvmin

    ! Initialise the level in the problem structure
    ilvmin = MIN(ilvmin,ilvmax)
    ilvmax = MAX(ilvmin,ilvmax)
    rproblem%NLMIN = ilvmin
    rproblem%NLMAX = ilvmax

    ! Allocate memory for the levels
    ALLOCATE(rproblem%RlevelInfo(1:ilvmax))
    
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
                              'iUpwind1',i1,0)
    CALL collct_setvalue_int( rproblem%rcollection,'IUPWIND1',i1,.TRUE.)

    CALL parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                'dUpsam1',d1,0.0_DP)
    CALL collct_setvalue_real (rproblem%rcollection,'UPSAM1',d1,.TRUE.)

    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iUpwind2',i1,0)
    CALL collct_setvalue_int( rproblem%rcollection,'IUPWIND2',i1,.TRUE.)

    CALL parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                'dUpsam2',d1,0.0_DP)
    CALL collct_setvalue_real (rproblem%rcollection,'UPSAM2',d1,.TRUE.)

    ! Type of boundary conditions
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iBoundary',i1,0)
    CALL collct_setvalue_int( rproblem%rcollection,'IBOUNDARY',i1,.TRUE.)

    ! Time dependence
    CALL cc_initParTimeDependence (rproblem,'TIME-DISCRETISATION',&
        rproblem%rparamList)
        
    ! Optimal control
    CALL cc_initOptControl(rproblem)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneParameters (rproblem)
  
!<description>
  ! Cleans up parameters read from the DAT files. Removes all references to
  ! parameters from the collection rproblem\%rcollection that were
  ! set up in cc_initParameters.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Optimal control
    CALL cc_doneOptControl(rproblem)

    ! Release memory of the level specification
    DEALLOCATE(rproblem%RlevelInfo)

    ! Remove information about boundary conditions
    CALL collct_deleteValue(rproblem%rcollection,'IBOUNDARY')
    
    ! Remove information about stabilisation
    CALL collct_deleteValue(rproblem%rcollection,'UPSAM1')
    CALL collct_deleteValue(rproblem%rcollection,'IUPWIND1')

    CALL collct_deleteValue(rproblem%rcollection,'UPSAM2')
    CALL collct_deleteValue(rproblem%rcollection,'IUPWIND2')

    ! Remove type of problem to discretise
    CALL collct_deleteValue(rproblem%rcollection,'ISTOKES')

    ! Remove the viscosity parameter
    CALL collct_deletevalue(rproblem%rcollection,'NU')
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initParamTriang (rproblem)
  
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
  INTEGER :: i,imeshType,ncellsX
  
    ! Variable for a filename:  
    CHARACTER(LEN=256) :: sString
    CHARACTER(LEN=60) :: sPRMFile, sTRIFile

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
        
    ! Now set up the basic triangulation. Which type of triangulation should
    ! be set up?
    CALL parlst_getvalue_int (rproblem%rparamList,'PARAMTRIANG',&
                              'imeshType',imeshType,0)
    SELECT CASE (imeshType)
    CASE (0)
      ! Standard mesh, specified by a TRI file.
      CALL tria_readTriFile2D (rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation, &
          sTRIFile, rproblem%p_rboundary)
    CASE (1)
      ! Sliced QUAD mesh with ncellsX cells on the coarse grid.
      CALL parlst_getvalue_int (rproblem%rparamList,'PARAMTRIANG',&
                                'ncellsX',ncellsX,1)
      CALL cc_generateSlicedQuadMesh(rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,&
          ncellsX)
    CASE DEFAULT
      CALL output_line ('Unknown mesh type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
      CALL sys_halt()
    END SELECT

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
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneParamTriang (rproblem)
  
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

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_generateSlicedQuadMesh (rtriangulation,ncellsX)
  
!<description>
  ! This routine generates a standard [0,1]^2 QUAD mesh with ncellsX cells in
  ! X-direction. By nature, these cells show an anisotropy of ncellsX:1.
!</description>
  
!<input>
  ! Number of cells in X-direction.
  INTEGER, INTENT(IN) :: ncellsX
!</input>
    
!<output>
  ! Triangulation structure that receives the triangulation.
  TYPE(t_triangulation), INTENT(OUT) :: rtriangulation
!</output>
    
    ! local variables
    INTEGER :: nbdVertives,nverticexX
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata2D
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_Idata2D
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata,p_IverticesAtBoundary,p_IboundaryCpIdx
    INTEGER(I32) :: ivt, iel
    INTEGER :: idim, ive
    INTEGER(I32), DIMENSION(2) :: Isize
    
    ! Initialise the basic mesh
    rtriangulation%ndim = NDIM2D
    rtriangulation%NEL = ncellsX
    rtriangulation%NVT = (ncellsX+1)*2
    rtriangulation%NMT = 0
    rtriangulation%NNVE = 4
    rtriangulation%NBCT = 1
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(TRIA_NVEQUAD2D) = rtriangulation%NEL
    
    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM2D, NVT)
    Isize = (/NDIM2D,INT(rtriangulation%NVT,I32)/)
    CALL storage_new2D ('tria_read_tri2D', 'DCORVG', Isize, ST_DOUBLE, &
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    CALL storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
    
    ! Initialise the point coordinates.
    ! Odd vertices on the bottom, even vertices on top of the QUAD mesh,
    ! numbered from left to right.
    DO ivt=0,ncellsX
      p_Ddata2D(1,2*ivt+1) = REAL(ivt,DP)/REAL(ncellsX,DP)
      p_Ddata2D(2,2*ivt+1) = 0.0_DP

      p_Ddata2D(1,2*ivt+2) = REAL(ivt,DP)/REAL(ncellsX,DP)
      p_Ddata2D(2,2*ivt+2) = 1.0_DP
    END DO
    
    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriangulation%NNVE,INT(rtriangulation%NEL,I32)/)
    CALL storage_new2D ('tria_read_tri2D', 'KVERT', Isize, ST_INT, &
        rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    CALL storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)

    ! Initialise the connectivity for the cells.
    DO iel=0,ncellsX-1
      p_Idata2D(1,iel+1) = 2*iel+1
      p_Idata2D(2,iel+1) = 2*iel+3
      p_Idata2D(3,iel+1) = 2*iel+4
      p_Idata2D(4,iel+1) = 2*iel+2
    END DO
    
    ! Allocate memory for InodalProperty 
    CALL storage_new ('tria_read_tri2D', 'KNPR', &
        INT(rtriangulation%NVT,I32), ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    CALL storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_Idata)

    ! All vertices are on the boundary
    p_Idata(:) = 1
    
    ! Number of vertices on the boundary -- all of them
    rtriangulation%NVBD = rtriangulation%NVT
    
    ! Allocate memory for IverticesAtBoundary.
    CALL storage_new ('tria_generateBasicBoundary', &
        'KVBD', INT(rtriangulation%NVBD,I32), &
        ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    CALL storage_new ('tria_generateBasicBoundary', &
        'KBCT', INT(rtriangulation%NBCT+1,I32), &
        ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    CALL storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    CALL storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! The first element in p_IboundaryCpIdx is (as the head) always =1.
    p_IboundaryCpIdx(1) = 1
    p_IboundaryCpIdx(2) = 1 + rtriangulation%NVBD

    ! Initialise the numbers of the vertices on the boundary
    DO ivt=0,ncellsX
      p_IverticesAtBoundary (ivt+1) = 2*ivt+1
      p_IverticesAtBoundary (2*(ncellsX+1)-ivt) = 2*ivt+2
    END DO
    
    ! Allocate memory for  and DvertexParameterValue
    CALL storage_new ('tria_generateBasicBoundary', &
        'DVBDP', INT(rtriangulation%NVBD,I32), &
        ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
    
    ! Get the array where to store boundary parameter values.
    CALL storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
    CALL storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
        
    ! Initialise the parameter values of the vertices. For the bottommost
    ! edge, they coincide wit the coordinate. For the topmost edge,
    ! that's 3 - x-coordinate.
    DO ivt=1,ncellsX+1
      p_DvertexParameterValue (ivt) = p_Ddata2D(1,p_IverticesAtBoundary(ivt))
      p_DvertexParameterValue (ncellsX+1+ivt) = &
          3.0_DP - p_Ddata2D(1,p_IverticesAtBoundary(ncellsX+1+ivt))
    END DO

  END SUBROUTINE

END MODULE
