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
    rproblem%rphysicsPrimal%iequation = i1

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

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'srhsExpressionX',sstring,'''''')
    read(sstring,*) rproblem%srhsExpressionX

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'srhsExpressionY',sstring,'''''')
    read(sstring,*) rproblem%srhsExpressionY
    
    ! If the RHS is given as expression, create a parser object for the
    ! expression.
    if (rproblem%irhs .eq. -1) then
      
      ! Create a parser object for the rhs expressions
      call fparser_create (rproblem%rrhsParser,NDIM2D)
      
      ! Compile the two expressions
      call fparser_parseFunction (rproblem%rrhsParser,&
          1, rproblem%srhsExpressionX, EXPR_VARIABLES)
      call fparser_parseFunction (rproblem%rrhsParser,&
          2, rproblem%srhsExpressionY, EXPR_VARIABLES)
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

    ! Probably release the parser object for the RHS
    if (rproblem%irhs .eq. -1) then
      call fparser_release(rproblem%rrhsParser)
    end if
    
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
    character(LEN=256) :: sString
    character(LEN=60) :: sPRMFile, sTRIFile

    ! Get the .prm and the .tri file from the parameter list.
    ! note that parlst_getvalue_string returns us exactly what stands
    ! in the parameter file, so we have to apply READ to get rid of
    ! probable ''!
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sParametrisation',sString)
    read (sString,*) sPRMFile
                              
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sMesh',sString)
    read (sString,*) sTRIFile
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rproblem%rboundary, sPrmFile)
        
    ! Now set up the basic triangulation. Which type of triangulation should
    ! be set up?
    call parlst_getvalue_int (rproblem%rparamList,'PARAMTRIANG',&
                              'imeshType',imeshType,0)
    select case (imeshType)
    case (0)
      ! Standard mesh, specified by a TRI file.
      call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation, &
          sTRIFile, rproblem%rboundary)
    case (1)
      ! Sliced QUAD mesh with ncellsX cells on the coarse grid.
      call parlst_getvalue_int (rproblem%rparamList,'PARAMTRIANG',&
                                'ncellsX',ncellsX,1)
      call cc_generateSlicedQuadMesh(rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,&
          ncellsX)
    case DEFAULT
      call output_line ('Unknown mesh type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
      call sys_halt()
    end select

    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(rproblem%NLMIN-1,&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)
    
    ! Now, refine to level up to nlmax.
    do i=rproblem%NLMIN+1,rproblem%NLMAX
      call tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
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
      call tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateSlicedQuadMesh (rtriangulation,ncellsX)
  
!<description>
  ! This routine generates a standard [0,1]^2 QUAD mesh with ncellsX cells in
  ! X-direction. By nature, these cells show an anisotropy of ncellsX:1.
!</description>
  
!<input>
  ! Number of cells in X-direction.
  integer, intent(IN) :: ncellsX
!</input>
    
!<output>
  ! Triangulation structure that receives the triangulation.
  type(t_triangulation), intent(OUT) :: rtriangulation
!</output>
    
    ! local variables
    integer :: nbdVertives,nverticexX
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer(I32), dimension(:,:), pointer :: p_Idata2D
    integer(I32), dimension(:), pointer :: p_Idata,p_IverticesAtBoundary,p_IboundaryCpIdx
    integer(I32) :: ivt, iel
    integer :: idim, ive
    integer(I32), dimension(2) :: Isize
    
    ! Initialise the basic mesh
    rtriangulation%ndim = NDIM2D
    rtriangulation%NEL = ncellsX
    rtriangulation%NVT = (ncellsX+1)*2
    rtriangulation%NMT = 0
    rtriangulation%NNVE = 4
    rtriangulation%NNEE = 4
    rtriangulation%NBCT = 1
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(TRIA_NVEQUAD2D) = rtriangulation%NEL
    
    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM2D, NVT)
    Isize = (/NDIM2D,int(rtriangulation%NVT,I32)/)
    call storage_new ('tria_read_tri2D', 'DCORVG', Isize, ST_DOUBLE, &
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
    
    ! Initialise the point coordinates.
    ! Odd vertices on the bottom, even vertices on top of the QUAD mesh,
    ! numbered from left to right.
    do ivt=0,ncellsX
      p_Ddata2D(1,2*ivt+1) = real(ivt,DP)/real(ncellsX,DP)
      p_Ddata2D(2,2*ivt+1) = 0.0_DP

      p_Ddata2D(1,2*ivt+2) = real(ivt,DP)/real(ncellsX,DP)
      p_Ddata2D(2,2*ivt+2) = 1.0_DP
    end do
    
    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriangulation%NNVE,int(rtriangulation%NEL,I32)/)
    call storage_new ('tria_read_tri2D', 'KVERT', Isize, ST_INT, &
        rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)

    ! Initialise the connectivity for the cells.
    do iel=0,ncellsX-1
      p_Idata2D(1,iel+1) = 2*iel+1
      p_Idata2D(2,iel+1) = 2*iel+3
      p_Idata2D(3,iel+1) = 2*iel+4
      p_Idata2D(4,iel+1) = 2*iel+2
    end do
    
    ! Allocate memory for InodalProperty
    call storage_new ('tria_read_tri2D', 'KNPR', &
        int(rtriangulation%NVT,I32), ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_Idata)

    ! All vertices are on the boundary
    p_Idata(:) = 1
    
    ! Number of vertices on the boundary -- all of them
    rtriangulation%NVBD = rtriangulation%NVT
    
    ! Allocate memory for IverticesAtBoundary.
    call storage_new ('tria_generateBasicBoundary', &
        'KVBD', int(rtriangulation%NVBD,I32), &
        ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    call storage_new ('tria_generateBasicBoundary', &
        'KBCT', int(rtriangulation%NBCT+1,I32), &
        ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    call storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    call storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! The first element in p_IboundaryCpIdx is (as the head) always =1.
    p_IboundaryCpIdx(1) = 1
    p_IboundaryCpIdx(2) = 1 + rtriangulation%NVBD

    ! Initialise the numbers of the vertices on the boundary
    do ivt=0,ncellsX
      p_IverticesAtBoundary (ivt+1) = 2*ivt+1
      p_IverticesAtBoundary (2*(ncellsX+1)-ivt) = 2*ivt+2
    end do
    
    ! Allocate memory for  and DvertexParameterValue
    call storage_new ('tria_generateBasicBoundary', &
        'DVBDP', int(rtriangulation%NVBD,I32), &
        ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
    
    ! Get the array where to store boundary parameter values.
    call storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
        
    ! Initialise the parameter values of the vertices. For the bottommost
    ! edge, they coincide wit the coordinate. For the topmost edge,
    ! that's 3 - x-coordinate.
    do ivt=1,ncellsX+1
      p_DvertexParameterValue (ivt) = p_Ddata2D(1,p_IverticesAtBoundary(ivt))
      p_DvertexParameterValue (ncellsX+1+ivt) = &
          3.0_DP - p_Ddata2D(1,p_IverticesAtBoundary(ncellsX+1+ivt))
    end do

  end subroutine

end module
