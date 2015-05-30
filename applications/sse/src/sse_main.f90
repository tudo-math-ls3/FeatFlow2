!##############################################################################
!# ****************************************************************************
!# <name> sse_main </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the main routines for solving the elliptic
!# equation for sea surface elevation.
!#
!#  1.) sse_initParamTriang
!#      -> Initialises the triangulation structure
!#
!#  2.) sse_initDiscretisation
!#      -> Initialises the discretisation structure
!#
!#  3.) sse_initMatVec
!#      -> Initialises the matrices and vectors
!#
!#  4.) sse_initDiscreteBC
!#      -> Initialises the discretised boundary conditions
!#
!#  5.) sse_implementBC
!#      -> Implements the discretised boundary conditions into the
!#         matrices and vectors
!#
!#  6.) sse_solve
!#      -> Solves the linear system
!#
!#  7.) sse_postprocessing
!#      -> Post-processes the solution
!#
!#  8.) sse_doneMatVec
!#      -> Releases matrices and vectors
!#
!#  9.) sse_doneDiscreteBC
!#      -> Releases the discretised boundary conditions
!#
!# 10.) sse_doneDiscretisation
!#      -> Releases the discretisation structure
!#
!# 11.) sse_doneParamTriang
!#      -> Releases the triangulation structure
!#
!# 12.) sse_outputTable
!#      -> Output the convergence table
!# </purpose>
!##############################################################################

module sse_main

  use bcassembly
  use bilinearformevaluation
  use boundary
  use collection
  use convergencetable
  use cubature
  use derivatives
  use discretebc
  use dofmapping
  use element
  use filtersupport
  use fparser
  use fsystem
  use genoutput
  use globalsystem
  use linearformevaluation
  use linearsolver
  use linearsystemblock
  use linearsystemscalar
  use matrixfilters
  use meshmodification
  use multileveloperators
  use multilevelprojection
  use paramlist
  use pprocerror
  use pprocgradients
  use scalarpde
  use sortstrategy
  use sortstrategybase
  use spatialdiscretisation
  use triangulation
  use vectorfilters

  use sse_base
  use sse_base_corine
  use sse_base_poisson
  use sse_base_sse
  use sse_callback_corine
  use sse_callback_poisson
  use sse_callback_sse

  implicit none

  public :: sse_initParamTriang
  public :: sse_initDiscretisation
  public :: sse_initMatVec
  public :: sse_initDiscreteBC
  public :: sse_implementBC
  public :: sse_solve
  public :: sse_postprocessing
  public :: sse_doneMatVec
  public :: sse_doneDiscreteBC
  public :: sse_doneDiscretisation
  public :: sse_doneParamTriang
  public :: sse_outputTable

  private

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamTriang(rparlist,ilvmin,ilvmax,rproblem)

!<description>
    ! This routine initialises the parametrisation and triangulation of
    ! the domain. The corresponding .prm/.tri files are read from disc
    ! and the triangulation is refined as described by the parameter
    ! ilv. Moreover, the boundaru conditions are initialised.
!</description>

!<input>
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist

    ! Minimum refinement level of the mesh; = coarse grid = level 1
    integer, intent(in) :: ilvmin

    ! Maximum refinement level
    integer, intent(in) :: ilvmax
!</input>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sconfig,spredir,sprmfile,ssection,strifile,sdata
    character(len=SYS_STRLEN) :: sbdrtype,spartype,sexpression
    real(DP) :: ddisturbMeshFactor,dminParam,dmaxParam
    integer :: i,icomp,ios,isegment,nsubstrings
    logical :: bminParam,bmaxParam

    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax
    allocate(rproblem%RlevelInfo(rproblem%ilvmin:rproblem%ilvmax))

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! Get the section for the configuration
    ssection = sse_getSection(rproblem%cproblemType)
    call parlst_getvalue_string(rparlist, ssection, 'problemconfig', sconfig)

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call parlst_getvalue_string(rparlist, sconfig, 'prmfile', sprmfile, '')
    if (trim(sprmfile) .ne. '') then
      call boundary_read_prm(rproblem%rboundary, trim(spredir)//'/'//trim(sprmfile))
    end if

    ! Now read in the basic triangulation.
    call parlst_getvalue_string(rparlist, sconfig, 'trifile', strifile)
    if (trim(sprmfile) .ne. '') then
      call tria_readTriFile2D(rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
          trim(spredir)//'/'//trim(strifile), rproblem%rboundary)
    else
      call tria_readTriFile2D(rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
          trim(spredir)//'/'//trim(strifile))
    end if

    ! Refine the mesh up to the minimum level
    if (trim(sprmfile) .ne. '') then
      call tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
          rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
    else
      call tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
          rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation)
    end if

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    if (trim(sprmfile) .ne. '') then
      call tria_initStandardMeshFromRaw(&
          rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
    else
      call tria_initStandardMeshFromRaw(&
          rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation)
    end if

    ! Now, refine to level up to nlmax.
    if (trim(sprmfile) .ne. '') then
      do i=rproblem%ilvmin+1,rproblem%ilvmax
        call tria_refine2LevelOrdering(rproblem%RlevelInfo(i-1)%rtriangulation,&
            rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
        call tria_initStandardMeshFromRaw(rproblem%RlevelInfo(i)%rtriangulation,&
            rproblem%rboundary)
      end do
    else
      do i=rproblem%ilvmin+1,rproblem%ilvmax
        call tria_refine2LevelOrdering(rproblem%RlevelInfo(i-1)%rtriangulation,&
            rproblem%RlevelInfo(i)%rtriangulation)
        call tria_initStandardMeshFromRaw(rproblem%RlevelInfo(i)%rtriangulation)
      end do
    end if

    ! Compress the level hierarchy.
    ! Share the vertex coordinates of all levels, so the coarse grid coordinates
    ! are "contained" in the fine grid coordinates. The effect is:
    ! 1.) Save some memory
    ! 2.) Every change in the fine grid coordinates also affects the coarse
    !     grid coordinates and vice versa.
    do i=rproblem%ilvmax-1,rproblem%ilvmin,-1
      call tria_compress2LevelOrdHierarchy(rproblem%RlevelInfo(i+1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation)
    end do

    ! Disturb the mesh if required
    call parlst_getvalue_double(rparlist, '',&
        'disturbmeshfactor', ddisturbMeshFactor, 0.0_DP)
    if (ddisturbMeshFactor .gt. 0.0_DP) then
      call meshmod_disturbMesh(rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
          ddisturbMeshFactor)
    end if

    ! Get number of boundary conditions
    nsubstrings = parlst_querysubstrings(rparlist, sconfig, 'BOUNDARYCONDITION')
    allocate(rproblem%RboundaryCondition(nsubstrings))

    ! Read boundary conditions
    do i=1,nsubstrings
      call parlst_getvalue_string(rparlist, sconfig, 'BOUNDARYCONDITION', sdata, isubstring=i)
      read(sdata,*,iostat=ios) icomp,&
          isegment,&
          spartype,&
          dminParam,&
          bminParam,&
          dmaxParam,&
          bmaxParam,&
          sbdrtype,&
          rproblem%RboundaryCondition(i)%iequation,&
          sexpression
      if (ios .ne. 0) then
        call output_line("An error occured while reading the boundary conditions.", &
            OU_CLASS_ERROR,OU_MODE_STD,"sse_initParamTriang")
        call sys_halt()
      end if

      ! Create boundary region
      if (trim(adjustl(spartype)) .eq. '01') then
        call boundary_createRegion(rproblem%rboundary, icomp, isegment, &
            rproblem%RboundaryCondition(i)%rboundaryRegion, BDR_PAR_01)
      else if (trim(adjustl(spartype)) .eq. 'LEN') then
        call boundary_createRegion(rproblem%rboundary, icomp, isegment, &
            rproblem%RboundaryCondition(i)%rboundaryRegion, BDR_PAR_LENGTH)
      else
        call output_line("Invalid type of boundary parametrisation.", &
            OU_CLASS_ERROR,OU_MODE_STD,"sse_initParamTriang")
        call sys_halt()
      end if

      ! Specify dminParam/dmaxParam
      if ((dminParam .gt. 0.0_DP) .or. (dmaxParam .gt. 0.0_DP)) then
        rproblem%RboundaryCondition(i)%rboundaryRegion%dminParam = dminParam
        rproblem%RboundaryCondition(i)%rboundaryRegion%dmaxParam = dmaxParam
      end if

      ! Specify iproperties
      rproblem%RboundaryCondition(i)%rboundaryRegion%iproperties = 0_I32
      if (bminParam) then
        rproblem%RboundaryCondition(i)%rboundaryRegion%iproperties = ior(&
            rproblem%RboundaryCondition(i)%rboundaryRegion%iproperties,&
            BDR_PROP_WITHSTART)
      end if
      if (bmaxParam) then
        rproblem%RboundaryCondition(i)%rboundaryRegion%iproperties = ior(&
            rproblem%RboundaryCondition(i)%rboundaryRegion%iproperties,&
            BDR_PROP_WITHEND)
      end if

      ! Set type of boundary condition
      if (trim(adjustl(sbdrtype)) .eq. 'DIR') then
        rproblem%RboundaryCondition(i)%cboundaryConditionType = BDRCOND_DIRICHLET
      else if (trim(adjustl(sbdrtype)) .eq. 'NEU') then
        rproblem%RboundaryCondition(i)%cboundaryConditionType = BDRCOND_NEUMANN
      else
        call output_line("Invalid type of boundary condition.", &
            OU_CLASS_ERROR,OU_MODE_STD,"sse_initParamTriang")
        call sys_halt()
      end if

      ! Parse expression and store expression number
      call fparser_parseFunction(rfparser, 'bdrCond'//trim(sys_si0(i,3)),&
          sexpression, (/'x','y'/), icomp=rproblem%RboundaryCondition(i)%iexpression)
    end do

  end subroutine sse_initParamTriang

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initDiscretisation(rparlist,rproblem)

!<description>
    ! This routine initialises the discretisation structure of the underlying
    ! problem and saves it to the problem structure.
!</description>

!<input>
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sconfig,sparameter,ssection
    integer, dimension(:), allocatable :: Ccubaturetypes,Celementtypes
    integer :: i,j,nvar,nsubstrings

    ! Allocate temporal storage
    nvar = sse_getNVAR(rproblem%cproblemType)
    allocate(Ccubaturetypes(nvar),Celementtypes(nvar))

    ! Read config section
    ssection = sse_getSection(rproblem%cproblemType)
    call parlst_getvalue_string(rparlist, ssection, 'problemconfig', sconfig)

    ! Get type of element(s)
    nsubstrings = parlst_querysubstrings(rparlist, sconfig, 'ELEMENTTYPE')
    if (nsubstrings .eq. 0) then
      call parlst_getvalue_string(rparlist, sconfig, 'ELEMENTTYPE', sparameter)
      Celementtypes(1) = elem_igetID(sparameter)
      do j = 2, nvar
        Celementtypes(j) = Celementtypes(1)
      end do
    elseif (nsubstrings .eq. nvar) then
      do j = 1, nvar
        call parlst_getvalue_string(rparlist, sconfig, 'ELEMENTTYPE', sparameter, isubstring=j)
        Celementtypes(j) = elem_igetID(sparameter)
      end do
    else
      call output_line("Mismatch in dimensions.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscretisation")
      call sys_halt()
    end if

    ! Get type of cubature formula(e)
    nsubstrings = parlst_querysubstrings(rparlist, sconfig, 'CUBATURETYPE')
    if (nsubstrings .eq. 0) then
      call parlst_getvalue_string(rparlist, sconfig, 'CUBATURETYPE', sparameter)
      Ccubaturetypes(1) = cub_igetID(sparameter)
      do j = 2, nvar
        Ccubaturetypes(j) = Ccubaturetypes(1)
      end do
    elseif (nsubstrings .eq. nvar) then
      do j = 1, nvar
        call parlst_getvalue_string(rparlist, sconfig, 'CUBATURETYPE', sparameter, isubstring=j)
        Ccubaturetypes(j) = cub_igetID(sparameter)
      end do
    else
      call output_line("Mismatch in dimensions.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscretisation")
      call sys_halt()
    end if

    ! Create spatial discretisations on all levels
    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Now we can start to initialise the discretisation. At first, set
      ! up a block discretisation structure that specifies the blocks in
      ! the solution vector. In this problem, we have three blocks.
      call spdiscr_initBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation,&
          nvar,rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)

      ! rproblem%RlevelInfo(i)%rdiscretisation%Rdiscretisations is a
      ! list of scalar discretisation structures for every component
      ! of the solution vector.  Initialise the first three elements
      ! of the list to specify the elements for the solution
      ! components:
      do j=1,nvar
        call spdiscr_initDiscr_simple(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
            Celementtypes(j),rproblem%RlevelInfo(i)%rtriangulation,&
            rproblem%rboundary)
      end  do

      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      do j=1,nvar
        call spdiscr_createDefCubStructure(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
            rproblem%RlevelInfo(i)%RcubatureInfo(j),&
            Ccubaturetypes(j))
      end do
    end do

  end subroutine sse_initDiscretisation

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initMatVec(rparlist,rproblem)

!<description>
    ! Calculates the system matrix and RHS vector of the linear system
    ! by discretising the problem with the default discretisation structure
    ! in the problem structure.
    ! Sets up a solution vector for the linear system.
!</description>

!<input>
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: ssortstrategy
    integer :: i

    select case(rproblem%cproblemType)
    case (POISSON_SCALAR, POISSON_SYSTEM)
      call sse_initMatVec_Poisson(rproblem)

    case (SSE_SCALAR, SSE_SYSTEM1, SSE_SYSTEM2)
      call sse_initMatVec_SSE(rproblem)

    case (CORINE_1D,CORINE_2D)
      call sse_initMatVec_Corine(rproblem)

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initMatVec")
      call sys_halt()
    end select

    ! Sorting does not work at the moment !!!
    return

    ! Get sort strategy from parameter list
    call parlst_getvalue_string(rparlist, '', 'SORTSTRATEGY', ssortstrategy, '')

    ! Apply sort strategy (if any)
    if (trim(ssortstrategy) .ne. '') then

      do i=rproblem%ilvmin,rproblem%ilvmax

        ! Create a sort strategy structure for our discretisation
        call sstrat_initBlockSorting(rproblem%RlevelInfo(i)%rsortStrategy,&
            rproblem%RlevelInfo(i)%rdiscretisation)

        ! Calculate the resorting strategy.
        if (trim(ssortstrategy) .eq. 'CM') then
          call sstrat_initCuthillMcKee(&
              rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
              rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1))

        elseif (trim(ssortstrategy) .eq. 'RCM') then
          call sstrat_initRevCuthillMcKee(&
              rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
              rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1))

        elseif (trim(ssortstrategy) .eq. 'XYZCOORD') then
          call sstrat_initXYZsorting(&
              rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
              rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1), 0)

        elseif (trim(ssortstrategy) .eq. 'ZYXCOORD') then
          call sstrat_initXYZsorting(&
              rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
              rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1), 1)

        elseif (trim(ssortstrategy) .eq. 'FEAST') then
          call sstrat_initFEASTsorting(&
              rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
              rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(1))

        elseif (trim(ssortstrategy) .eq. 'RANDOM') then
          call sstrat_initRandom(&
              rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1))

        else

          call output_line("Invalid sort strategy.", &
              OU_CLASS_ERROR,OU_MODE_STD,"sse_initMatVec")
          call sys_halt()

        end if

        ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
        call lsysbl_setSortStrategy(rproblem%RlevelInfo(i)%rmatrix,&
            rproblem%RlevelInfo(i)%rsortStrategy,&
            rproblem%RlevelInfo(i)%rsortStrategy)

        ! Resort the matrix.
        call lsysbl_sortMatrix(rproblem%RlevelInfo(i)%rmatrix,.true.)

      end do

      ! Install the resorting strategy in the RHS- and the solution
      ! vector, but do not resort them yet!
      ! We resort the vectors just before solving.

      call lsysbl_setSortStrategy(rproblem%rrhs,&
          rproblem%RlevelInfo(rproblem%ilvmax)%rsortStrategy)
      call lsysbl_setSortStrategy(rproblem%rvector,&
          rproblem%RlevelInfo(rproblem%ilvmax)%rsortStrategy)

    end if

  end subroutine sse_initMatVec

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initDiscreteBC(rproblem)

!<description>
    ! This calculates the discrete version of the boundary conditions and
    ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j,cessential

    ! Prepare quick-access arrays of collection structure
    rproblem%rcollection%IquickAccess(1) = rproblem%cproblemType

    ! Dirichlet boundary conditions are essential in the primary formulation
    select case(rproblem%cproblemType)
    case (POISSON_SCALAR, SSE_SCALAR, CORINE_1D, CORINE_2D)
      cessential = BDRCOND_DIRICHLET
    case (POISSON_SYSTEM, SSE_SYSTEM1, SSE_SYSTEM2)
      cessential = BDRCOND_NEUMANN
    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscreteBC")
      call sys_halt()
    end select

    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Create a t_discreteBC structure where we store all discretised
      ! boundary conditions.
      call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)

      ! Loop over all boundary conditions
      do j=1,size(rproblem%RboundaryCondition)

        if (rproblem%RboundaryCondition(j)%cboundaryConditionType .eq. cessential) then

          ! Set position, equation number and boundary component and assemble linear form
          rproblem%rcollection%IquickAccess(2) = rproblem%RboundaryCondition(j)%iequation
          rproblem%rcollection%IquickAccess(3) = rproblem%RboundaryCondition(j)%iexpression
          rproblem%rcollection%IquickAccess(4) = rproblem%RboundaryCondition(j)%rboundaryRegion%iboundCompIdx

          select case(rproblem%cproblemType)
          case (POISSON_SCALAR, POISSON_SYSTEM)
            call bcasm_newDirichletBConRealBD(&
                rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,&
                rproblem%RboundaryCondition(j)%iequation,&
                rproblem%RboundaryCondition(j)%rboundaryRegion,&
                rproblem%RlevelInfo(i)%rdiscreteBC,&
                getBoundaryValues_Poisson,rproblem%rcollection)

          case (SSE_SCALAR, SSE_SYSTEM1, SSE_SYSTEM2)
            call bcasm_newDirichletBConRealBD(&
                rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,&
                rproblem%RboundaryCondition(j)%iequation,&
                rproblem%RboundaryCondition(j)%rboundaryRegion,&
                rproblem%RlevelInfo(i)%rdiscreteBC,&
                getBoundaryValues_SSE,rproblem%rcollection)

          case (CORINE_1D, CORINE_2D)
            call bcasm_newDirichletBConRealBD(&
                rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,&
                rproblem%RboundaryCondition(j)%iequation,&
                rproblem%RboundaryCondition(j)%rboundaryRegion,&
                rproblem%RlevelInfo(i)%rdiscreteBC,&
                getBoundaryValues_Corine,rproblem%rcollection)

            call output_line("Invalid type of problem.", &
                OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscreteBC")
            call sys_halt()
          end select

        end if
      end do
    end do

  end subroutine sse_initDiscreteBC

  ! ***************************************************************************

!<subroutine>

  subroutine sse_implementBC(rproblem)

!<description>
    ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    call vecfil_discreteBCrhs(rproblem%rrhs,&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscreteBC)
    call vecfil_discreteBCsol(rproblem%rvector,&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscreteBC)

    ! Implement discrete boundary conditions into the matrices on all
    ! levels, too. Call the appropriate matrix filter to modify
    ! all matrices according to the attached discrete boundary conditions.
    do i=rproblem%ilvmin,rproblem%ilvmax
      call matfil_discreteBC(rproblem%RlevelInfo(i)%rmatrix,&
          rproblem%RlevelInfo(i)%rdiscreteBC)
    end do

  end subroutine sse_implementBC

  ! ***************************************************************************

!<subroutine>

  subroutine sse_solve(rparlist,rproblem)

!<description>
    ! Solves the given problem by applying a linear solver.
!</description>

!<input>
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilvmin,ilvmax
    integer :: i,j

    ! Error indicator during initialisation of the solver
    integer :: ierror

    ! A pointer to the system matrix and the RHS vector
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector

    ! Temporal vectors and matrices
    type(t_matrixBlock), target :: rmatTmp
    type(t_vectorBlock), target :: rvecTmp,rvecTmp1,rvecTmp2,rvecTmp3,rrhsTmp2

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNodeA,p_rsmootherA
    type(t_linsolNode), pointer :: p_rcoarseGridSolverA,p_rpreconditionerA

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNodeS,p_rsmootherS
    type(t_linsolNode), pointer :: p_rcoarseGridSolverS,p_rpreconditionerS

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    type(t_matrixBlock), dimension(:), pointer :: RmatricesA
    type(t_matrixBlock), dimension(:), pointer :: RmatricesS

    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfoA
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfoS

    ! Iteration counter
    integer :: ite

    ! Norm of residuals
    real(DP) :: dresNorm0,dresNorm


    type(t_vectorScalar) :: rvectorSc,rvecTmpSc,rrhsSc
    type(t_vectorBlock) :: rvectorBl,rvecTmpBl,rrhsBl

    select case(sse_getSolver(rproblem%cproblemType))
    case (SOLVER_SCALAR)
      !---------------------------------------------------------------------------

      ilvmin = rproblem%ilvmin
      ilvmax = rproblem%ilvmax

      ! Get our right hand side / solution / matrix on the finest
      ! level from the problem structure.
      p_rrhs    => rproblem%rrhs
      p_rvector => rproblem%rvector
      p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix

      ! Set up prolongation matrices
      do i=ilvmin+1,ilvmax
        do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents

          ! Create the matrix structure of the prolongation matrix.
          call mlop_create2LvlMatrixStruct(&
              rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(j),&
              rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
              LSYSSC_MATRIX9, rproblem%RlevelInfo(i)%rmatProl(j))

          ! Assemble the entries of the prolongation matrix.
          call mlop_build2LvlProlMatrix(&
              rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(j),&
              rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
              .true., rproblem%RlevelInfo(i)%rmatProl(j),&
              rcubatureInfoCoarse=rproblem%RlevelInfo(i-1)%RcubatureInfo(1),&
              rcubatureInfoFine=rproblem%RlevelInfo(i)%RcubatureInfo(1))

          ! Assemble the entries of the restriction matrix.
          call lsyssc_transposeMatrix(rproblem%RlevelInfo(i)%rmatProl(j),&
              rproblem%RlevelInfo(i)%rmatRest(j),LSYSSC_TR_ALL)
        end do

        ! Now set up an interlevel projecton structure for this level
        ! based on the system matrix on this level.
        call mlprj_initProjectionMat(rproblem%RlevelInfo(i)%rprojection,&
            rproblem%RlevelInfo(i)%rmatrix)

        ! Initialise the matrix-based projection
        do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents
          call mlprj_initMatrixProjection(&
              rproblem%RlevelInfo(i)%rprojection%RscalarProjection(1,j),&
              rproblem%RlevelInfo(i)%rmatProl(j),&
              rmatrixRest=rproblem%RlevelInfo(i)%rmatRest(j))
        end do
      end do

      ! Set up an interlevel projecton structure for the coarse-most level.
      call mlprj_initProjectionMat(rproblem%RlevelInfo(ilvmin)%rprojection,&
          rproblem%RlevelInfo(ilvmin)%rmatrix)

      ! Create a temporary vector we need that for some preparation.
      call lsysbl_createVector(p_rrhs, rvecTmp, .false.)

      ! Now we have to build up the level information for multigrid.
      !
      ! Create a Multigrid-solver. Attach the above filter chain
      ! to the solver, so that the solver automatically filters
      ! the vector during the solution process.
      call linsol_initMultigrid2(p_rsolverNode,ilvmax-ilvmin+1)

      ! Get parameters of multigrid solver from parameter list
      call parlst_getvalue_int(rparlist, 'MG-SOLVER', 'nminIterations', &
          p_rsolverNode%nminIterations, p_rsolverNode%nminIterations)
      call parlst_getvalue_int(rparlist, 'MG-SOLVER', 'nmaxIterations', &
          p_rsolverNode%nmaxIterations, p_rsolverNode%nmaxIterations)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER', 'depsRel', &
          p_rsolverNode%depsRel, p_rsolverNode%depsRel)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER', 'depsAbs', &
          p_rsolverNode%depsAbs, p_rsolverNode%depsAbs)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER', 'depsDiff', &
          p_rsolverNode%depsDiff, p_rsolverNode%depsDiff)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER', 'ddivRel', &
          p_rsolverNode%ddivRel, p_rsolverNode%ddivRel)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER', 'ddivAbs', &
          p_rsolverNode%ddivAbs, p_rsolverNode%ddivAbs)

      ! Then set up smoothers / coarse grid solver:
      do i=ilvmin,ilvmax

        ! Set up a filter chain for implementing boundary conditions on that level
        call filter_initFilterChain(rproblem%RlevelInfo(i)%RfilterChain,&
            rproblem%RlevelInfo(i)%nfilters)
        call filter_newFilterDiscBCDef(rproblem%RlevelInfo(i)%RfilterChain,&
            rproblem%RlevelInfo(i)%nfilters,rproblem%RlevelInfo(i)%rdiscreteBC)

        ! On the coarsest grid, set up a coarse grid solver and no smoother
        ! On finer grids, set up a smoother but no coarse grid solver.
        nullify(p_rpreconditioner)
        nullify(p_rsmoother)
        nullify(p_rcoarseGridSolver)

        if (i .eq. ilvmin) then
          if ((rproblem%RlevelInfo(ilvmin)%rmatrix%nblocksPerRow .eq. 1) .and.&
              (rproblem%RlevelInfo(ilvmin)%rmatrix%nblocksPerCol .eq. 1)) then
            ! Set up a BiCGStab coarse grid solver (with ILU preconditioning)
            call linsol_initMILUs1x1(p_rpreconditioner,0,0.0_DP)
            call linsol_initBiCGStab(p_rcoarseGridSolver,p_rpreconditioner,&
                rproblem%RlevelInfo(i)%RfilterChain)
          else
            ! Set up a GMRES coarse grid solver
            call linsol_initGMRES(p_rcoarseGridSolver,20,&
                Rfilter=rproblem%RlevelInfo(i)%RfilterChain)
          end if

          ! Get parameters of coarse-grid solver from parameter list
          call parlst_getvalue_int(rparlist, 'COARSE-GRID-SOLVER', 'nminIterations', &
              p_rcoarseGridSolver%nminIterations, p_rcoarseGridSolver%nminIterations)
          call parlst_getvalue_int(rparlist, 'COARSE-GRID-SOLVER', 'nmaxIterations', &
              p_rcoarseGridSolver%nmaxIterations, p_rcoarseGridSolver%nmaxIterations)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER', 'depsRel', &
              p_rcoarseGridSolver%depsRel, p_rcoarseGridSolver%depsRel)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER', 'depsAbs', &
              p_rcoarseGridSolver%depsAbs, p_rcoarseGridSolver%depsAbs)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER', 'depsDiff', &
              p_rcoarseGridSolver%depsDiff, p_rcoarseGridSolver%depsDiff)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER', 'ddivRel', &
              p_rcoarseGridSolver%ddivRel, p_rcoarseGridSolver%ddivRel)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER', 'ddivAbs', &
              p_rcoarseGridSolver%ddivAbs, p_rcoarseGridSolver%ddivAbs)

        else
          if ((rproblem%RlevelInfo(ilvmin)%rmatrix%nblocksPerRow .eq. 1) .and.&
              (rproblem%RlevelInfo(ilvmin)%rmatrix%nblocksPerCol .eq. 1)) then
            ! Set up an ILU smoother for multigrid with damping
            ! parameter 0.7, 4 smoothing steps:
            call linsol_initMILUs1x1(p_rsmoother,0,0.0_DP)
            call linsol_convertToSmoother(p_rsmoother,4,0.7_DP)
          else
            ! Set up GMRES smoother for multigrid with damping parameter
            ! 0.7, 4 smoothing steps:
            call linsol_initGMRES(p_rsmoother,20,&
                Rfilter=rproblem%RlevelInfo(i)%RfilterChain)
            call linsol_convertToSmoother(p_rsmoother,4,0.7_DP)
          end if
        end if

        ! And add this multi-grid level. We will use the same smoother
        ! for pre- and post-smoothing.
        call linsol_getMultigrid2Level(p_rsolverNode,i-ilvmin+1,p_rlevelInfo)
        p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
        p_rlevelInfo%p_rpresmoother      => p_rsmoother
        p_rlevelInfo%p_rpostsmoother     => p_rsmoother

        ! Attach the filter chain which imposes boundary conditions on that level.
        p_rlevelInfo%p_RfilterChain => rproblem%RlevelInfo(i)%RfilterChain

        ! Attach our user-defined projection to the level.
        call linsol_initProjMultigrid2Level(p_rlevelInfo,&
            rproblem%RlevelInfo(i)%rprojection)
      end do

      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 2

      ! Attach the system matrices to the solver.
      !
      ! We copy our matrices to a big matrix array and transfer that
      ! to the setMatrices routines. This intitialises then the matrices
      ! on all levels according to that array. Note that this does not
      ! allocate new memory, we create only "links" to existing matrices
      ! into Rmatrices(:)!
      allocate(Rmatrices(ilvmin:ilvmax))
      do i=ilvmin,ilvmax
        call lsysbl_duplicateMatrix(rproblem%RlevelInfo(i)%rmatrix,&
            Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end do

      call linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))

      ! We can release Rmatrices immediately -- as long as we do not
      ! release rproblem%RlevelInfo(i)%rmatrix!
      do i=ilvmin,ilvmax
        call lsysbl_releaseMatrix(Rmatrices(i))
      end do
      deallocate(Rmatrices)

      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure(p_rsolverNode, ierror)

      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
        call sys_halt()
      end if

      call linsol_initData(p_rsolverNode, ierror)

      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix singular!",OU_CLASS_ERROR)
        call sys_halt()
      end if

      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call linsol_solveAdaptively(p_rsolverNode,p_rvector,p_rrhs,rvecTmp)

      ! Release solver data and structure
      call linsol_doneData(p_rsolverNode)
      call linsol_doneStructure(p_rsolverNode)

      ! Release the solver node and all subnodes attached to it (if at all):
      call linsol_releaseSolver(p_rsolverNode)

      ! Release the prolongation matrices
      do i=ilvmax,ilvmin+1,-1
        ! Release the projection structure itself
        call mlprj_doneProjection(rproblem%RlevelInfo(i)%rprojection)

        ! Release the prolongation matrix
        do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents
          call lsyssc_releaseMatrix(rproblem%RlevelInfo(i)%rmatProl(j))
        end  do

        ! Release the restriction matrix
        do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents
          call lsyssc_releaseMatrix(rproblem%RlevelInfo(i)%rmatRest(j))
        end do
      end do

      ! Release the projection structure on the coarse mesh
      call mlprj_doneProjection(rproblem%RlevelInfo(ilvmin)%rprojection)

      ! Release the filter chain
      do i=ilvmin,ilvmax
        call filter_doneFilterChain(rproblem%RlevelInfo(i)%RfilterChain,&
            rproblem%RlevelInfo(i)%nfilters)
      end do

      ! Release the temporary vector
      call lsysbl_releaseVector(rvecTmp)

    case (SOLVER_SADDLEPOINT)
      !---------------------------------------------------------------------------

      ilvmin = rproblem%ilvmin
      ilvmax = rproblem%ilvmax

      ! Get our right hand side / solution / matrix on the finest
      ! level from the problem structure.
      p_rrhs    => rproblem%rrhs
      p_rvector => rproblem%rvector
      p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix

      ! Set up prolongation matrices
      do i=ilvmin+1,ilvmax
        do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents

          ! Create the matrix structure of the prolongation matrix.
          call mlop_create2LvlMatrixStruct(&
              rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(j),&
              rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
              LSYSSC_MATRIX9, rproblem%RlevelInfo(i)%rmatProl(j))

          ! Assemble the entries of the prolongation matrix.
          call mlop_build2LvlProlMatrix(&
              rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(j),&
              rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
              .true., rproblem%RlevelInfo(i)%rmatProl(j),&
              rcubatureInfoCoarse=rproblem%RlevelInfo(i-1)%RcubatureInfo(j),&
              rcubatureInfoFine=rproblem%RlevelInfo(i)%RcubatureInfo(j))

          ! Assemble the entries of the restriction matrix.
          call lsyssc_transposeMatrix(rproblem%RlevelInfo(i)%rmatProl(j),&
              rproblem%RlevelInfo(i)%rmatRest(j),LSYSSC_TR_ALL)
        end do

        ! Derive submatrix for A-component
        call lsysbl_deriveSubmatrix(rproblem%RlevelInfo(i)%rmatrix,&
            rmatTmp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE,2,3,2,3)

        ! Now set up an interlevel projecton structure for this level
        ! based on the matrix A on this level.
        call mlprj_initProjectionMat(rproblem%RlevelInfo(i)%rprojectionA,&
            rmatTmp)

        ! Release temporal matrix
        call lsysbl_releaseMatrix(rmatTmp)

        ! Initialise the matrix-based projection
        call mlprj_initMatrixProjection(&
            rproblem%RlevelInfo(i)%rprojectionA%RscalarProjection(1,1),&
            rproblem%RlevelInfo(i)%rmatProl(2),&
            rmatrixRest=rproblem%RlevelInfo(i)%rmatRest(2))
        call mlprj_initMatrixProjection(&
            rproblem%RlevelInfo(i)%rprojectionA%RscalarProjection(1,2),&
            rproblem%RlevelInfo(i)%rmatProl(3),&
            rmatrixRest=rproblem%RlevelInfo(i)%rmatRest(3))

        ! Derive submatrix for S-component
        call lsysbl_deriveSubmatrix(rproblem%RlevelInfo(i)%rmatrix,&
            rmatTmp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE,1,1,1,1)

        ! Now set up an interlevel projecton structure for this level
        ! based on the matrix A on this level.
        call mlprj_initProjectionMat(rproblem%RlevelInfo(i)%rprojectionS,&
            rmatTmp)

        ! Initialise the matrix-based projection
        call mlprj_initMatrixProjection(&
            rproblem%RlevelInfo(i)%rprojectionS%RscalarProjection(1,1),&
            rproblem%RlevelInfo(i)%rmatProl(1),&
            rmatrixRest=rproblem%RlevelInfo(i)%rmatRest(1))
      end do

      ! Set up an interlevel projecton structure for the coarse-most level.
      call lsysbl_deriveSubmatrix(rproblem%RlevelInfo(ilvmin)%rmatrix,&
          rmatTmp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE,2,3,2,3)
      call mlprj_initProjectionMat(rproblem%RlevelInfo(ilvmin)%rprojectionA,&
          rmatTmp)
      call lsysbl_releaseMatrix(rmatTmp)

      ! Set up an interlevel projecton structure for the coarse-most level.
      call lsysbl_deriveSubmatrix(rproblem%RlevelInfo(ilvmin)%rmatrix,&
          rmatTmp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE,1,1,1,1)
      call mlprj_initProjectionMat(rproblem%RlevelInfo(ilvmin)%rprojectionS,&
          rmatTmp)
      call lsysbl_releaseMatrix(rmatTmp)

      stop

      ! Create a temporary vector we need that for some preparation.
      call lsysbl_createVector(p_rrhs, rvecTmp, .false.)

      ! Now we have to build up the level information for multigrid.
      !
      ! Create two Multigrid-solver. Attach the above filter chain
      ! to the solver, so that the solver automatically filters
      ! the vector during the solution process.
      call linsol_initMultigrid2(p_rsolverNodeA,ilvmax-ilvmin+1)
      call linsol_initMultigrid2(p_rsolverNodeS,ilvmax-ilvmin+1)

      ! Get parameters of multigrid solver for A from parameter list
      call parlst_getvalue_int(rparlist, 'MG-SOLVER-A', 'nminIterations', &
          p_rsolverNodeA%nminIterations, p_rsolverNodeA%nminIterations)
      call parlst_getvalue_int(rparlist, 'MG-SOLVER-A', 'nmaxIterations', &
          p_rsolverNodeA%nmaxIterations, p_rsolverNodeA%nmaxIterations)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-A', 'depsRel', &
          p_rsolverNodeA%depsRel, p_rsolverNodeA%depsRel)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-A', 'depsAbs', &
          p_rsolverNodeA%depsAbs, p_rsolverNodeA%depsAbs)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-A', 'depsDiff', &
          p_rsolverNodeA%depsDiff, p_rsolverNodeA%depsDiff)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-A', 'ddivRel', &
          p_rsolverNodeA%ddivRel, p_rsolverNodeA%ddivRel)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-A', 'ddivAbs', &
          p_rsolverNodeA%ddivAbs, p_rsolverNodeA%ddivAbs)

      call parlst_getvalue_int(rparlist, 'MG-SOLVER-S', 'nminIterations', &
          p_rsolverNodeS%nminIterations, p_rsolverNodeS%nminIterations)
      call parlst_getvalue_int(rparlist, 'MG-SOLVER-S', 'nmaxIterations', &
          p_rsolverNodeS%nmaxIterations, p_rsolverNodeS%nmaxIterations)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-S', 'depsRel', &
          p_rsolverNodeS%depsRel, p_rsolverNodeS%depsRel)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-S', 'depsAbs', &
          p_rsolverNodeS%depsAbs, p_rsolverNodeS%depsAbs)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-S', 'depsDiff', &
          p_rsolverNodeS%depsDiff, p_rsolverNodeS%depsDiff)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-S', 'ddivRel', &
          p_rsolverNodeS%ddivRel, p_rsolverNodeS%ddivRel)
      call parlst_getvalue_double(rparlist, 'MG-SOLVER-S', 'ddivAbs', &
          p_rsolverNodeS%ddivAbs, p_rsolverNodeS%ddivAbs)

      ! Then set up smoothers / coarse grid solver:
      do i=ilvmin,ilvmax

        ! Set up a filter chain for implementing boundary conditions on that level
        call filter_initFilterChain(rproblem%RlevelInfo(i)%RfilterChain,&
            rproblem%RlevelInfo(i)%nfilters)
        call filter_newFilterDiscBCDef(rproblem%RlevelInfo(i)%RfilterChain,&
            rproblem%RlevelInfo(i)%nfilters,rproblem%RlevelInfo(i)%rdiscreteBC)

        ! On the coarsest grid, set up a coarse grid solver and no smoother
        ! On finer grids, set up a smoother but no coarse grid solver.
        nullify(p_rpreconditionerA,p_rpreconditionerS)
        nullify(p_rsmootherA,p_rsmootherS)
        nullify(p_rcoarseGridSolverA,p_rcoarseGridSolverS)

        stop

        if (i .eq. ilvmin) then
          ! Set up a GMRES coarse grid solver for A-component
          call linsol_initGMRES(p_rcoarseGridSolverA,20,&
              Rfilter=rproblem%RlevelInfo(i)%RfilterChain)

          ! Get parameters of coarse-grid solver from parameter list
          call parlst_getvalue_int(rparlist, 'COARSE-GRID-SOLVER-A', 'nminIterations', &
              p_rcoarseGridSolverA%nminIterations, p_rcoarseGridSolverA%nminIterations)
          call parlst_getvalue_int(rparlist, 'COARSE-GRID-SOLVER-A', 'nmaxIterations', &
              p_rcoarseGridSolverA%nmaxIterations, p_rcoarseGridSolverA%nmaxIterations)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER0S', 'depsRel', &
              p_rcoarseGridSolverA%depsRel, p_rcoarseGridSolverA%depsRel)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-A', 'depsAbs', &
              p_rcoarseGridSolverA%depsAbs, p_rcoarseGridSolverA%depsAbs)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-A', 'depsDiff', &
              p_rcoarseGridSolverA%depsDiff, p_rcoarseGridSolverA%depsDiff)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-A', 'ddivRel', &
              p_rcoarseGridSolverA%ddivRel, p_rcoarseGridSolverA%ddivRel)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-A', 'ddivAbs', &
              p_rcoarseGridSolverA%ddivAbs, p_rcoarseGridSolverA%ddivAbs)

          ! Set up a BiCGStab coarse grid solver (with ILU
          ! preconditioning) for S-component
          call linsol_initMILUs1x1(p_rpreconditionerS,0,0.0_DP)
          call linsol_initBiCGStab(p_rcoarseGridSolverS,p_rpreconditionerS,&
              rproblem%RlevelInfo(i)%RfilterChain)

          ! Get parameters of coarse-grid solver from parameter list
          call parlst_getvalue_int(rparlist, 'COARSE-GRID-SOLVER-S', 'nminIterations', &
              p_rcoarseGridSolverS%nminIterations, p_rcoarseGridSolverS%nminIterations)
          call parlst_getvalue_int(rparlist, 'COARSE-GRID-SOLVER-S', 'nmaxIterations', &
              p_rcoarseGridSolverS%nmaxIterations, p_rcoarseGridSolverS%nmaxIterations)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER0S', 'depsRel', &
              p_rcoarseGridSolverS%depsRel, p_rcoarseGridSolverS%depsRel)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-S', 'depsAbs', &
              p_rcoarseGridSolverS%depsAbs, p_rcoarseGridSolverS%depsAbs)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-S', 'depsDiff', &
              p_rcoarseGridSolverS%depsDiff, p_rcoarseGridSolverS%depsDiff)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-S', 'ddivRel', &
              p_rcoarseGridSolverS%ddivRel, p_rcoarseGridSolverS%ddivRel)
          call parlst_getvalue_double(rparlist, 'COARSE-GRID-SOLVER-S', 'ddivAbs', &
              p_rcoarseGridSolverS%ddivAbs, p_rcoarseGridSolverS%ddivAbs)

        else
          ! Set up GMRES smoother for multigrid with damping parameter
          ! 0.7, 4 smoothing steps for A-component:
          call linsol_initGMRES(p_rsmootherA,20,&
              Rfilter=rproblem%RlevelInfo(i)%RfilterChain)
          call linsol_convertToSmoother(p_rsmootherA,4,0.7_DP)

          ! Set up an ILU smoother for multigrid with damping
          ! parameter 0.7, 4 smoothing steps for S-component:
          call linsol_initMILUs1x1(p_rsmootherS,0,0.0_DP)
          call linsol_convertToSmoother(p_rsmootherS,4,0.7_DP)        
        end if

        ! And add this multi-grid level. We will use the same smoother
        ! for pre- and post-smoothing.
        call linsol_getMultigrid2Level(p_rsolverNodeA,i-ilvmin+1,p_rlevelInfoA)
        p_rlevelInfoA%p_rcoarseGridSolver => p_rcoarseGridSolverA
        p_rlevelInfoA%p_rpresmoother      => p_rsmootherA
        p_rlevelInfoA%p_rpostsmoother     => p_rsmootherA

        call linsol_getMultigrid2Level(p_rsolverNodeS,i-ilvmin+1,p_rlevelInfoS)
        p_rlevelInfoS%p_rcoarseGridSolver => p_rcoarseGridSolverS
        p_rlevelInfoS%p_rpresmoother      => p_rsmootherS
        p_rlevelInfoS%p_rpostsmoother     => p_rsmootherS

        ! Attach the filter chain which imposes boundary conditions on that level.
        p_rlevelInfoA%p_RfilterChain => rproblem%RlevelInfo(i)%RfilterChain
        p_rlevelInfoS%p_RfilterChain => rproblem%RlevelInfo(i)%RfilterChain

        ! Attach our user-defined projection to the level.
        call linsol_initProjMultigrid2Level(p_rlevelInfoS,&
            rproblem%RlevelInfo(i)%rprojection)
      end do

      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 2

      ! Attach the system matrices to the solver.
      !
      ! We copy our matrices to a big matrix array and transfer that
      ! to the setMatrices routines. This intitialises then the matrices
      ! on all levels according to that array. Note that this does not
      ! allocate new memory, we create only "links" to existing matrices
      ! into Rmatrices(:)!
      allocate(Rmatrices(ilvmin:ilvmax))
      do i=ilvmin,ilvmax
!!$      call lsysbl_deriveSubmatrix(rproblem%RlevelInfo(i)%rmatrix,&
!!$          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE,2,3,2,3)
        call glsys_assembleGlobal(rproblem%RlevelInfo(i)%rmatrix,&
            Rmatrices(i), .true., .true.)
      end do

      call linsol_setMatrices(p_rsolverNode,Rmatrices(ilvmin:ilvmax))

      ! We can release Rmatrices immediately -- as long as we do not
      ! release rproblem%RlevelInfo(i)%rmatrix!
!!$    do i=ilvmin,ilvmax
!!$      call lsysbl_releaseMatrix(Rmatrices(i))
!!$    end do
!!$    deallocate(Rmatrices)

      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure(p_rsolverNode, ierror)

      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
        call sys_halt()
      end if

      call linsol_initData(p_rsolverNode, ierror)

      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix singular!",OU_CLASS_ERROR)
        call sys_halt()
      end if

      call lsysbl_createScalarFromVec(p_rvector,rvectorSc,.true.)
      call lsysbl_createScalarFromVec(p_rrhs,rrhsSc,.true.)
      call lsysbl_createScalarFromVec(rvecTmp,rvecTmpSc,.true.)

      call lsysbl_createVecFromScalar(rvectorSc, rvectorBl)
      call lsysbl_createVecFromScalar(rrhsSc, rrhsBl)
      call lsysbl_createVecFromScalar(rvecTmpSc, rvecTmpBl)

      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
!!$    call linsol_solveAdaptively(p_rsolverNode,p_rvector,p_rrhs,rvecTmp)

      call linsol_solveAdaptively(p_rsolverNode,rvectorBl,rrhsBl,rvecTmpBl)

!!$    ! Clear solution vector
!!$    call lsysbl_clearVector(p_rvector)
!!$
!!$    ! Compute norm of initial residual
!!$    dresNorm0 = lsysbl_vectorNorm(p_rrhs, LINALG_NORML2)
!!$
!!$    do ite = 1,5
!!$      ! Compute right-hand side -B[k]*N, k=1,2
!!$      call lsyssc_matVec(rproblem%RlevelInfo(ilvmax)%rmatrix%RmatrixBlock(2,1),&
!!$          rvecTmp1%RvectorBlock(1), rrhsTmp2%RvectorBlock(1), -1.0_DP, 0.0_DP)
!!$      call lsyssc_matVec(rproblem%RlevelInfo(ilvmax)%rmatrix%RmatrixBlock(3,1),&
!!$          rvecTmp1%RvectorBlock(1), rrhsTmp2%RvectorBlock(2), -1.0_DP, 0.0_DP)
!!$
!!$      print *, lsysbl_vectorNorm(rrhsTmp2, LINALG_NORML2)
!!$      pause
!!$
!!$      ! Solve the system A[k]*sigma[k]=rhs[k], k=1,2
!!$      call linsol_solveAdaptively(p_rsolverNode,rvecTmp2,rrhsTmp2,rvecTmp3)
!!$
!!$      print *, lsysbl_vectorNorm(rvecTmp2, LINALG_NORML2)
!!$      pause
!!$
!!$
!!$      ! Update solution
!!$      call lsyssc_vectorLinearComb(p_rrhs%RvectorBlock(1),&
!!$          rvecTmp1%RvectorBlock(1), -domega, 1.0_DP)
!!$
!!$      print *, lsysbl_vectorNorm(p_rrhs, LINALG_NORMMAX)
!!$      stop
!!$
!!$      call lsyssc_matVec(rproblem%RlevelInfo(ilvmax)%rmatrix%RmatrixBlock(1,2),&
!!$          rvecTmp2%RvectorBlock(1), rvecTmp1%RvectorBlock(1), domega, 1.0_DP)
!!$      call lsyssc_matVec(rproblem%RlevelInfo(ilvmax)%rmatrix%RmatrixBlock(1,3),&
!!$          rvecTmp2%RvectorBlock(2), rvecTmp1%RvectorBlock(1), domega, 1.0_DP)
!!$
!!$      ! Compute norm of residual
!!$      call lsysbl_copyVector(p_rrhs, rvecTmp)
!!$      call lsysbl_matVec(rproblem%RlevelInfo(ilvmax)%rmatrix, p_rvector,&
!!$          rvecTmp, -1.0_DP, 1.0_DP)
!!$      dresNorm = lsysbl_vectorNorm(rvecTmp, LINALG_NORML2)
!!$
!!$      print *, "NORM", dresNorm,dresNorm0
!!$      pause
!!$    end do

      ! Release solver data and structure
      call linsol_doneData(p_rsolverNode)
      call linsol_doneStructure(p_rsolverNode)

      ! Release the solver node and all subnodes attached to it (if at all):
      call linsol_releaseSolver(p_rsolverNode)

      ! Release the prolongation matrices
      do i=ilvmax,ilvmin+1,-1
        ! Release the projection structure itself
        call mlprj_doneProjection(rproblem%RlevelInfo(i)%rprojection)

        ! Release the prolongation matrix
        do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents
          call lsyssc_releaseMatrix(rproblem%RlevelInfo(i)%rmatProl(j))
        end  do

        ! Release the restriction matrix
        do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents
          call lsyssc_releaseMatrix(rproblem%RlevelInfo(i)%rmatRest(j))
        end do
      end do

      ! Release the projection structure on the coarse mesh
      call mlprj_doneProjection(rproblem%RlevelInfo(ilvmin)%rprojection)
      call mlprj_doneProjection(rproblem%RlevelInfo(ilvmin)%rprojectionA)
      call mlprj_doneProjection(rproblem%RlevelInfo(ilvmin)%rprojectionS)

      ! Release the filter chain
      do i=ilvmin,ilvmax
        call filter_doneFilterChain(rproblem%RlevelInfo(i)%RfilterChain,&
            rproblem%RlevelInfo(i)%nfilters)
      end do

      ! Release the temporary vector
      call lsysbl_releaseVector(rvecTmp)
!!$    call lsysbl_releaseVector(rvecTmp1)
!!$    call lsysbl_releaseVector(rvecTmp2)
!!$    call lsysbl_releaseVector(rrhsTmp2)

    case(SOLVER_SYSTEM)
      !---------------------------------------------------------------------------

      ilvmin = rproblem%ilvmin
      ilvmax = rproblem%ilvmax

      ! Get our right hand side / solution / matrix on the finest
      ! level from the problem structure.
      p_rrhs    => rproblem%rrhs
      p_rvector => rproblem%rvector
      p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix

      call linsol_initUMFPACK4(p_rsolverNode)

      ! Set up a filter chain for implementing boundary conditions on that level
      call filter_initFilterChain(rproblem%RlevelInfo(ilvmax)%RfilterChain,&
          rproblem%RlevelInfo(ilvmax)%nfilters)
      call filter_newFilterDiscBCDef(rproblem%RlevelInfo(ilvmax)%RfilterChain,&
          rproblem%RlevelInfo(ilvmax)%nfilters,rproblem%RlevelInfo(ilvmax)%rdiscreteBC)

      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 2

      allocate(Rmatrices(1))
      call glsys_assembleGlobal(rproblem%RlevelInfo(ilvmax)%rmatrix,&
          Rmatrices(1), .true., .true.)
      call linsol_setMatrices(p_rsolverNode,Rmatrices(1:1))

      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure(p_rsolverNode, ierror)

      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
        call sys_halt()
      end if

      call linsol_initData(p_rsolverNode, ierror)

      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix singular!",OU_CLASS_ERROR)
        call sys_halt()
      end if

      ! Create temporal vectors
      call lsysbl_createScalarFromVec(p_rvector,rvectorSc,.true.)
      call lsysbl_createScalarFromVec(p_rrhs,rrhsSc,.true.)
      call lsysbl_createScalarFromVec(rvecTmp,rvecTmpSc,.true.)

      call lsysbl_createVecFromScalar(rvectorSc, rvectorBl)
      call lsysbl_createVecFromScalar(rrhsSc, rrhsBl)
      call lsysbl_createVecFromScalar(rvecTmpSc, rvecTmpBl)

      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call linsol_solveAdaptively(p_rsolverNode,rvectorBl,rrhsBl,rvecTmpBl)

      ! Release solver data and structure
      call linsol_doneData(p_rsolverNode)
      call linsol_doneStructure(p_rsolverNode)

      ! Release the solver node and all subnodes attached to it (if at all):
      call linsol_releaseSolver(p_rsolverNode)

      call lsysbl_releaseMatrix(Rmatrices(1))
      deallocate(Rmatrices)

      ! Release temporal vectors
      call lsysbl_releaseVector(rvectorBl)
      call lsysbl_releaseVector(rrhsBl)
      call lsysbl_releaseVector(rvecTmpBl)

      call lsyssc_releaseVector(rvectorSc)
      call lsyssc_releaseVector(rrhsSc)
      call lsyssc_releaseVector(rvecTmpSc)

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_solve")
      call sys_halt()
    end select

  end subroutine sse_solve

  ! ***************************************************************************

!<subroutine>

  subroutine sse_postprocessing(rparlist,rproblem,rtable)

!<description>
    ! Writes the solution into a VTK file.
!</description>

!<input>
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist

    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</input>

!<inputoutput>
    ! Convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    ! Output block for UCD output to VTK file
    character(len=SYS_STRLEN) :: sucddir,sucdfile

    ! Error of FE function to reference function
    real(DP) :: derror

    ! Number of DOFs
    integer :: i,ndof,cucdtype,iucdtype

    ! Method for calculating the gradient
    integer :: cgradType,cgradSubtype

    ! Get method for calculating the gradient
    call parlst_getvalue_int(rparlist, '', 'GRADTYPE', cgradType)
    call parlst_getvalue_int(rparlist, '', 'GRADSUBTYPE', cgradSubtype)

    ! Get parameters for UCD output
    if (.not. sys_getenv_string("UCDDIR", sucddir))&
        call parlst_getvalue_string(rparlist, '', 'SUCDDIR', sucddir, './ucd')
    call parlst_getvalue_string(rparlist, '', 'UCDFILE', sucdfile, '')
    call parlst_getvalue_int(rparlist, '', 'UCDTYPE', cucdtype, 0)

    ! Calculate the total number of DoF's
    ndof = 0
    do i=1,rproblem%rvector%nblocks
      ndof = ndof + dof_igetNDofGlob(rproblem%rvector%RvectorBlock(i)%p_rspatialDiscr)
    end do

    ! Add number of cells and vertices
    call ctab_addValue(rtable, "cells",&
        rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL)
    call ctab_addValue(rtable, "dofs", ndof)

    select case(rproblem%cproblemType)
    case (POISSON_SCALAR,POISSON_SYSTEM)
      call sse_postprocessing_Poisson(rproblem, cgradType, cgradSubtype,&
          sucdDir, sucdFile, cucdType, rtable)

    case (SSE_SCALAR,SSE_SYSTEM1,SSE_SYSTEM2)
      call sse_postprocessing_SSE(rproblem, cgradType, cgradSubtype,&
          sucdDir, sucdFile, cucdType, rtable)

    case (CORINE_1D,CORINE_2D)
      call sse_postprocessing_Corine(rproblem,&
          sucdDir, sucdFile, cucdType, rtable)
      
    end select

  end subroutine sse_postprocessing

  ! ***************************************************************************

!<subroutine>

  subroutine sse_doneMatVec(rproblem)

!<description>
    ! Releases system matrix and vectors.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Delete solution/RHS vector
    call lsysbl_releaseVector(rproblem%rvector)
    call lsysbl_releaseVector(rproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue(rproblem%rcollection,"RHS")
    call collct_deletevalue(rproblem%rcollection,"SOLUTION")

    ! Release matrices and vectors on all levels
    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the matrix
      call lsysbl_releaseMatrix(rproblem%RlevelInfo(i)%rmatrix)

      ! Delete the variables from the collection.
      call collct_deletevalue(rproblem%rcollection,"MATRIX",i)

      ! Release the sorting strategy
      call sstrat_doneBlockSorting(rproblem%RlevelInfo(i)%rsortStrategy)

    end do

  end subroutine sse_doneMatVec

  ! ***************************************************************************

!<subroutine>

  subroutine sse_doneDiscreteBC(rproblem)

!<description>
    ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
    end do

  end subroutine sse_doneDiscreteBC


  ! ***************************************************************************

!<subroutine>

  subroutine sse_doneDiscretisation(rproblem)

!<description>
    ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the cubature info structure.
      do j=1,size(rproblem%RlevelInfo(i)%RcubatureInfo)
        call spdiscr_releaseCubStructure(rproblem%RlevelInfo(i)%RcubatureInfo(j))
      end do

      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation)
    end do

  end subroutine sse_doneDiscretisation

  ! ***************************************************************************

!<subroutine>

  subroutine sse_doneParamTriang(rproblem)

!<description>   
    ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Release the triangulation
    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      call tria_done(rproblem%RlevelInfo(i)%rtriangulation)
    end do

    ! Release the problem level structure
    deallocate(rproblem%RlevelInfo)

    ! Release boundary condition structure
    deallocate(rproblem%RboundaryCondition)

    ! Release the boundary
    call boundary_release(rproblem%rboundary)

  end subroutine sse_doneParamTriang

  ! ***************************************************************************

!<subroutine>

  subroutine sse_outputTable(rproblem,rtable)

!<description>
    ! Exports the convergence table to file.
!</description>

!<input>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(in) :: rproblem
!</input>

!<inputoutput>
    ! A convergence table.
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    select case(rproblem%cproblemType)
    case (POISSON_SCALAR,POISSON_SYSTEM)
      call sse_outputTable_Poisson(rproblem,rtable)

    case (SSE_SCALAR,SSE_SYSTEM1,SSE_SYSTEM2)
      call sse_outputTable_SSE(rproblem,rtable)

    case (CORINE_1D,CORINE_2D)
      call sse_outputTable_Corine(rproblem,rtable)

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_outputTable")
      call sys_halt()
    end select

  end subroutine sse_outputTable

end module sse_main
