!##############################################################################
!# ****************************************************************************
!# <name> sse_main </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the main routines for solving the elliptic
!# equation for sea surface elevation.
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
  use fsystem
  use genoutput
  use globalsystem
  use linearalgebra
  use linearformevaluation
  use linearsolver
  use linearsystemblock
  use linearsystemscalar
  use matrixfilters
  use matrixio
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
  use storage
  use triangulation
  use ucd
  use vectorfilters
  use vectorio

  use sse_base
  use sse_callback

  implicit none

  public :: t_problem
  public :: sse_initParamTriang
  public :: sse_initDiscretisation
  public :: sse_initMatVec
  public :: sse_initDiscreteBC
  public :: sse_implementBC
  public :: sse_solve
  public :: sse_postprocessing
  public :: sse_doneMatVec
  public :: sse_doneBC
  public :: sse_doneDiscretisation
  public :: sse_doneParamTriang
  public :: sse_outputTable

  private

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl

    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo), dimension(3) :: RcubatureInfo

    ! A system matrix for that specific level. The matrix will receive the
    ! discrete system operator.
    type(t_matrixBlock) :: rmatrix

    ! A scalar matrix that will recieve the prolongation matrix for this level.
    type(t_matrixScalar), dimension(6) :: rmatProl,rmatRest

    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC

    ! Sorting strategy for resorting vectors/matrices.
    type(t_blockSortStrategy) :: rsortStrategy

    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1) :: RfilterChain

    ! Number of filters in the filter chain.
    integer :: nfilters

  end type

!</typeblock>

!<typeblock description="Application-specific type block for SSE problem">

  type t_problem

    ! Problem type
    integer :: cproblemtype

    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: ilvmin

    ! Maximum refinement level
    integer :: ilvmax

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation.
    type(t_problem_lvl), dimension(:), pointer :: RlevelInfo

    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamTriang(rparlist,ilvmin,ilvmax,rproblem)

!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
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
  real(DP) :: ddisturbMeshFactor
  integer :: i

  ! Path to the mesh
  character(len=SYS_STRLEN) :: spredir,sconfig,strifile,sprmfile

  ! Initialise the level in the problem structure
  rproblem%ilvmin = ilvmin
  rproblem%ilvmax = ilvmax
  allocate(rproblem%RlevelInfo(rproblem%ilvmin:rproblem%ilvmax))

  ! Get the path $PREDIR from the environment, where to read .prm/.tri files
  ! from. If that does not exist, write to the directory "./pre".
  if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

  ! Get the section for the configuration
  call parlst_getvalue_string(rparlist, '', 'config', sconfig, '')

  ! At first, read in the parametrisation of the boundary and save
  ! it to rboundary.
  call parlst_getvalue_string(rparlist, trim(sconfig), 'prmfile', sprmfile, '')
  if (trim(sprmfile) .ne. '') then
    call boundary_read_prm(rproblem%rboundary, trim(spredir)//'/'//trim(sprmfile))
  end if

  ! Now read in the basic triangulation.
  call parlst_getvalue_string(rparlist, trim(sconfig), 'trifile', strifile)
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

  end subroutine

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
  character(len=SYS_STRLEN) :: sconfig,sparameter
  integer, dimension(3) :: Ccubaturetypes,Celementtypes
  integer :: i,j,ccubaturetype,celementtype

  ! Get the section for the configuration
  call parlst_getvalue_string(rparlist, '', 'config', sconfig, '')

  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR,SSE_SCALAR)

    ! Get type of element
    call parlst_getvalue_string(rparlist, trim(sconfig), 'ELEMENTTYPE', sparameter)
    celementtype = elem_igetID(sparameter)

    ! Get type of cubature formula
    call parlst_getvalue_string(rparlist, trim(sconfig), 'CUBATURETYPE', sparameter)
    ccubaturetype = cub_igetID(sparameter)

    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      call spdiscr_initBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation,&
          merge(1,2,rproblem%cproblemtype .eq. POISSON_SCALAR),&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)

      ! rproblem%RlevelInfo(i)%rdiscretisation%Rdiscretisations is a
      ! list of scalar discretisation structures for every component
      ! of the solution vector. Initialise the first element of the
      ! list to specify the element for this solution component:
      do j=1,merge(1,2,rproblem%cproblemtype .eq. POISSON_SCALAR)
        call spdiscr_initDiscr_simple(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j), celementtype,&
            rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      end do

      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      do j=1,merge(1,2,rproblem%cproblemtype .eq. POISSON_SCALAR)
        call spdiscr_createDefCubStructure(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            ccubaturetype)
      end do
    end do

  case (POISSON_SYSTEM)

    ! Get types of element for each variable
    do j=1,3
      call parlst_getvalue_string(rparlist, trim(sconfig), 'ELEMENTTYPE', sparameter,&
          isubstring=j)
      Celementtypes(j) = elem_igetID(sparameter)
    end do

    ! Get types of cubature formula
    do j=1,3
      call parlst_getvalue_string(rparlist, trim(sconfig), 'CUBATURETYPE', sparameter,&
          isubstring=j)
      Ccubaturetypes(j) = cub_igetID(sparameter)
    end do

    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this problem, we have three blocks.
      call spdiscr_initBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation,&
          3,rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)

      ! rproblem%RlevelInfo(i)%rdiscretisation%Rdiscretisations is a
      ! list of scalar discretisation structures for every component
      ! of the solution vector.  Initialise the first three elements
      ! of the list to specify the elements for the solution
      ! components:
      do j=1,3
        call spdiscr_initDiscr_simple(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
            Celementtypes(j),rproblem%RlevelInfo(i)%rtriangulation,&
            rproblem%rboundary)
      end  do

      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      do j=1,3
        call spdiscr_createDefCubStructure(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
            rproblem%RlevelInfo(i)%RcubatureInfo(j),&
            Ccubaturetypes(j))
      end do
    end do

  case (SSE_SYSTEM1,SSE_SYSTEM2)

    ! Get types of element for each variable
    do j=1,3
      call parlst_getvalue_string(rparlist, trim(sconfig), 'ELEMENTTYPE', sparameter,&
          isubstring=j)
      Celementtypes(j) = elem_igetID(sparameter)
    end do

    ! Get types of cubature formula
    do j=1,3
      call parlst_getvalue_string(rparlist, trim(sconfig), 'CUBATURETYPE', sparameter,&
          isubstring=j)
      Ccubaturetypes(j) = cub_igetID(sparameter)
    end do

    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this problem, we have three blocks.
      call spdiscr_initBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation,&
          6,rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)

      ! rproblem%RlevelInfo(i)%rdiscretisation%Rdiscretisations is a
      ! list of scalar discretisation structures for every component
      ! of the solution vector.  Initialise the first three elements
      ! of the list to specify the elements for the solution
      ! components:
      do j=1,3
        call spdiscr_initDiscr_simple(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(2*j-1),&
            Celementtypes(j),rproblem%RlevelInfo(i)%rtriangulation,&
            rproblem%rboundary)
        call spdiscr_initDiscr_simple(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(2*j),&
            Celementtypes(j),rproblem%RlevelInfo(i)%rtriangulation,&
            rproblem%rboundary)
      end  do

      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      do j=1,3
        call spdiscr_createDefCubStructure(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(2*j-1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(j),&
            Ccubaturetypes(j))
        call spdiscr_createDefCubStructure(&
            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(2*j),&
            rproblem%RlevelInfo(i)%RcubatureInfo(j),&
            Ccubaturetypes(j))
      end do
    end do

  case default
    call output_line("Invalid type of problem.", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscretisation")
    call sys_halt()
  end select

  end subroutine

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

  ! A boundary segment
  type(t_boundaryRegion) :: rboundaryRegion

  ! A bilinear and linear form describing the analytic problem to solve
  type(t_bilinearForm) :: rform
  type(t_linearForm) :: rlinform

  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR)
    !---------------------------------------------------------------------------
    !
    ! Problem formulation in (test,trial)-notation
    !
    ! (grad_x,grad_x)*u + (grad_y,grad_y)*u = (func,f)
    !
    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatrix(rproblem%RlevelInfo(i)%rdiscretisation,&
          rproblem%RlevelInfo(i)%rmatrix)

      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"MATRIX",&
          rproblem%RlevelInfo(i)%rmatrix,.true.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,&
          1, 1, LSYSSC_MATRIX9)

      ! And now to the entries of the matrix. For assembling the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%Dcoefficients(1)  = 1.0
      rform%Dcoefficients(2)  = 1.0

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! For debugging purposes only
      !rform%ballCoeffConstant = .false.
      !
      !call bilf_buildMatrixScalar(rform,.true.,&
      !    rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
      !    rproblem%RlevelInfo(i)%RcubatureInfo(1),&
      !    coeff_Matrix_Poisson,rproblem%rcollection)
    end do

    ! Next step: Create a RHS vector and a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rrhs,.true.)
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rvector,.true.)

    ! Save the solution/RHS vector to the collection.
    call collct_setvalue_vec(rproblem%rcollection,"RHS",rproblem%rrhs,.true.)
    call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",rproblem%rvector,.true.)

    ! The vector structure is ready but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC

    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    call linf_buildVectorScalar(&
        rlinform,.true.,rproblem%rrhs%RvectorBlock(1),&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1),&
        coeff_RHS_Poisson,rproblem%rcollection)

    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)

  case (POISSON_SYSTEM)
    !---------------------------------------------------------------------------
    !
    ! Problem formulation in (test,trial)-notation
    !
    !  /                                         \   /       \   /         \
    ! |      0        (func,grad_x) (func,grad_y) | |    u    | | (func,-f) |
    ! |                                           | |         | |           |
    ! | (grad_x,func)  (func,func)       0        |*| sigma_x |=|     0     |
    ! |                                           | |         | |           |
    ! | (grad_y,func)       0        (func,func)  | | sigma_y | |     0     |
    !  \                                         /   \       /   \         /
    !
    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatrix(rproblem%RlevelInfo(i)%rdiscretisation,&
          rproblem%RlevelInfo(i)%rmatrix)

      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"MATRIX",&
          rproblem%RlevelInfo(i)%rmatrix,.true.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,2,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,3,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,1,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,2,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,1,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,3,LSYSSC_MATRIX9)

      ! (1,2)-block (w,grad_x sigma_x)
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (1,3)-block (w,grad_y sigma_y)
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_DERIV_Y
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,3),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (2,1)-block (grad_x v_x,u)
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Dcoefficients(1)  = 1.0

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(2),&
          rcollection=rproblem%rcollection)

      ! (3,1)-block (grad_y v_y,u)
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y
      rform%Dcoefficients(1)  = 1.0

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(3),&
          rcollection=rproblem%rcollection)

      ! (2,2)- and (3,3)-blocks (v_x,sigma_x) and (v_y,sigma_y)
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(2),&
          rcollection=rproblem%rcollection)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,3),&
          rproblem%RlevelInfo(i)%RcubatureInfo(3),&
          rcollection=rproblem%rcollection)
    end do

    ! Next step: Create a RHS vector and a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rrhs,.true.)
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rvector,.true.)

    ! Save the solution/RHS vector to the collection.
    call collct_setvalue_vec(rproblem%rcollection,"RHS",rproblem%rrhs,.true.)
    call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",rproblem%rvector,.true.)

    ! The vector structure is ready but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC

    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    call linf_buildVectorScalar(&
        rlinform,.true.,rproblem%rrhs%RvectorBlock(1),&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1),&
        coeff_RHS_Poisson,rproblem%rcollection)
    call lsyssc_scaleVector(rproblem%rrhs%RvectorBlock(1),-1.0_DP)

    ! Clear the second and thrid component of the RHS vector
    call lsyssc_clearVector(rproblem%rrhs%RvectorBlock(2))
    call lsyssc_clearVector(rproblem%rrhs%RvectorBlock(3))

    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)

  case (SSE_SCALAR)
    !---------------------------------------------------------------------------
    !
    ! Problem formulation in (test,trial)-notation
    !
    !  /                                                                         \
    ! | (grad,-Re(A)*grad)                   (grad,Im(A)*grad)-\omega*(func,func) |
    ! |                                                                           |
    ! | (grad,-Im(A)*grad)+\omega(func,func) (grad,-Re(A)*grad)                   |
    !  \                                                                         /
    !
    !   /     \   / \
    !  | Re(N) | | 0 |
    ! *|       |=|   |
    !  | Im(N) | | 0 |
    !   \     /   \ /

    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatrix(rproblem%RlevelInfo(i)%rdiscretisation,&
          rproblem%RlevelInfo(i)%rmatrix)

      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"MATRIX",&
          rproblem%RlevelInfo(i)%rmatrix,.true.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,&
          1, 1, LSYSSC_MATRIX9)

      ! Anisotropic diffusion matrix (real part)
      rform%itermCount = 4
      rform%ballCoeffConstant = .false.

      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_X
      rform%Idescriptors(1,3) = DER_DERIV_X
      rform%Idescriptors(2,3) = DER_DERIV_Y
      rform%Idescriptors(1,4) = DER_DERIV_Y
      rform%Idescriptors(2,4) = DER_DERIV_Y

      ! Assemble matrix block(1,1)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          coeff_MatrixA_Real,rproblem%rcollection)

      ! Assemble matrix block(2,2)
      call lsyssc_duplicateMatrix(&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      ! Anisotropic diffusion matrix (imaginary part)
      ! plus consistent mass matrix
      rform%itermCount = 5
      rform%ballCoeffConstant = .false.

      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_X
      rform%Idescriptors(1,3) = DER_DERIV_X
      rform%Idescriptors(2,3) = DER_DERIV_Y
      rform%Idescriptors(1,4) = DER_DERIV_Y
      rform%Idescriptors(2,4) = DER_DERIV_Y
      rform%Idescriptors(1,5) = DER_FUNC
      rform%Idescriptors(2,5) = DER_FUNC

      ! Duplicate matrix structure
      call lsyssc_duplicateMatrix(&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix(&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      ! Assemble matrix block(1,2)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          coeff_MatrixA_Aimag,rproblem%rcollection)
      call lsyssc_scaleMatrix(&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),-1.0_DP)

      ! Assemble matrix block(2,1)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          coeff_MatrixA_Aimag,rproblem%rcollection)
    end do

    ! Next step: Create a RHS vector and a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rrhs,.true.)
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rvector,.true.)

    ! Save the solution/RHS vector to the collection.
    call collct_setvalue_vec(rproblem%rcollection,"RHS",rproblem%rrhs,.true.)
    call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",rproblem%rvector,.true.)

    ! Clear the RHS vector on the finest level.
    call lsysbl_clearVector(rproblem%rrhs)

    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)

  case (SSE_SYSTEM1,SSE_SYSTEM2)
    !---------------------------------------------------------------------------
    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatrix(rproblem%RlevelInfo(i)%rdiscretisation,&
          rproblem%RlevelInfo(i)%rmatrix)

      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"MATRIX",&
          rproblem%RlevelInfo(i)%rmatrix,.true.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,2,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,3,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,5,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,1,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,4,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,6,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,1,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,3,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,4,2,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,4,4,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,5,1,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,5,5,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,6,2,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,6,6,LSYSSC_MATRIX9)

      ! (1,2)- and (2,1)-block -\omega*(w,Re(N))
      !                    and  \omega*(w,Im(N))
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = -dtidalfreq

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      rform%Dcoefficients(1)  = dtidalfreq

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (1,3)- and (2,4)-block (w,grad_x Re(sigma_x))
      !                    and (w,grad_x Im(sigma_x))
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0_DP

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,3),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,4),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (1,5)- and (2,6)-block (w,grad_x Re(sigma_x))
      !                    and (w,grad_x Im(sigma_x))
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_DERIV_Y
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0_DP

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,5),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,6),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (3,1)- and (4,2)-block (grad_x v_x,Re(N))
      !                    and (grad_x v_x,Im(N))
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Dcoefficients(1)  = 1.0_DP

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(4,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (5,1)- and (6,2)-block (grad_y v_y,Re(N))
      !                    and (grad_y v_y,Im(N))
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y
      rform%Dcoefficients(1)  = 1.0_DP

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(6,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (3,3)- and (4,4)-block (v_x,Re(sigma_x))
      !                    and (v_x,Im(sigma_x))
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0_DP

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,3),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(4,4),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)

      ! (5,5)- and (6,6)-block (v_y,Re(sigma_y))
      !                    and (v_y,Im(sigma_y))
      rform%itermCount = 1
      rform%ballCoeffConstant = .true.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0_DP

      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,5),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(6,6),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          rcollection=rproblem%rcollection)
    end do

    ! Next step: Create a RHS vector and a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rrhs,.true.)
    call lsysbl_createVector(&
        rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
        rproblem%rvector,.true.)

    ! Save the solution/RHS vector to the collection.
    call collct_setvalue_vec(rproblem%rcollection,"RHS",rproblem%rrhs,.true.)
    call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",rproblem%rvector,.true.)

    ! Clear the RHS vector on the finest level.
    call lsysbl_clearVector(rproblem%rrhs)

    ! Create boundary region
    call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND

    ! Initialise the linear form along the boundary
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_DERIV_X

    ! Assemble the linear forms
    call linf_buildVectorScalarBdr2d(rlinform, CUB_G3_1D, .false.,&
        rproblem%rrhs%RvectorBlock(3), coeff_RHSBdr_Real, rboundaryRegion)
    call linf_buildVectorScalarBdr2d(rlinform, CUB_G3_1D, .false.,&
        rproblem%rrhs%RvectorBlock(4), coeff_RHSBdr_Aimag, rboundaryRegion)

    ! Initialise the linear form along the boundary
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_DERIV_Y

    ! Assemble the linear forms
    call linf_buildVectorScalarBdr2d(rlinform, CUB_G3_1D, .false.,&
        rproblem%rrhs%RvectorBlock(5), coeff_RHSBdr_Real, rboundaryRegion)
    call linf_buildVectorScalarBdr2d(rlinform, CUB_G3_1D, .false.,&
        rproblem%rrhs%RvectorBlock(6), coeff_RHSBdr_Aimag, rboundaryRegion)

    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)

  case default
    call output_line("Invalid type of problem.", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_initMatVec")
    call sys_halt()
  end select

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

  end subroutine

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
  integer :: i,iboundComp,iboundSeg
  type(t_boundaryRegion) :: rboundaryRegion


  do i=rproblem%ilvmin,rproblem%ilvmax

    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)

    select case(rproblem%cproblemtype)
    case (POISSON_SCALAR)

      do iboundComp=1,boundary_igetNBoundComp(rproblem%rboundary)

        do iboundSeg=1,boundary_igetNsegments(rproblem%rboundary,iboundComp)

          ! We ask the boundary routines to create a "boundary region"
          ! - which is simply a part of the boundary corresponding to
          ! a boundary segment.  A boundary region roughly contains
          ! the type, the min/max parameter value and whether the
          ! endpoints are inside the region or not.
          call boundary_createRegion(rproblem%rboundary,iboundComp,iboundSeg,&
              rboundaryRegion)

          ! We use this boundary region and specify that we want to
          ! have Dirichlet boundary there. The following call does the
          ! following:
          ! - Create Dirichlet boundary conditions on the region
          !   rboundaryRegion.  We specify icomponent="1" to indicate
          !   that we set up the Dirichlet BC`s for the first (here:
          !   one and only) component in the solution vector.
          ! - Discretise the boundary condition so that the BC`s can
          !   be applied to matrices and vectors
          ! - Add the calculated discrete BC`s to rdiscreteBC for
          !   later use.
          call bcasm_newDirichletBConRealBD(&
              rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
              rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
              getBoundaryValues_Poisson,rproblem%rcollection)
        end do
      end do

    case (SSE_SCALAR)

#if defined(CASE_ALEX)

      ! We ask the boundary routines to create a "boundary region"
      ! - which is simply a part of the boundary corresponding to
      ! a boundary segment.  A boundary region roughly contains
      ! the type, the min/max parameter value and whether the
      ! endpoints are inside the region or not.
      call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND

      ! Real part of the solution
      call bcasm_newDirichletBConRealBD(&
          rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues_Real,rproblem%rcollection)

      ! Imaginary part
      call bcasm_newDirichletBConRealBD(&
          rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,2,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues_Aimag,rproblem%rcollection)

#elif defined(CASE_MARCHI)

      do iboundComp=1,boundary_igetNBoundComp(rproblem%rboundary)

        do iboundSeg=1,boundary_igetNsegments(rproblem%rboundary,iboundComp)

          ! We ask the boundary routines to create a "boundary region"
          ! - which is simply a part of the boundary corresponding to
          ! a boundary segment.  A boundary region roughly contains
          ! the type, the min/max parameter value and whether the
          ! endpoints are inside the region or not.
          call boundary_createRegion(rproblem%rboundary,iboundComp,iboundSeg,&
              rboundaryRegion)

          ! Real part of the solution
          call bcasm_newDirichletBConRealBD(&
              rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
              rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
              getBoundaryValues_Real,rproblem%rcollection)
          
          ! Imaginary part
          call bcasm_newDirichletBConRealBD(&
              rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,2,&
              rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
              getBoundaryValues_Aimag,rproblem%rcollection)
        end do
      end do
      
#elif defined(CASE_WALTERS)

      ! We ask the boundary routines to create a "boundary region"
      ! - which is simply a part of the boundary corresponding to
      ! a boundary segment.  A boundary region roughly contains
      ! the type, the min/max parameter value and whether the
      ! endpoints are inside the region or not.
      call boundary_createRegion(rproblem%rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND

      ! Real part of the solution
      call bcasm_newDirichletBConRealBD(&
          rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues_Real,rproblem%rcollection)

      ! Imaginary part
      call bcasm_newDirichletBConRealBD(&
          rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,2,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues_Aimag,rproblem%rcollection)

#elif defined(CASE_WINANT)
      
      ! We ask the boundary routines to create a "boundary region"
      ! - which is simply a part of the boundary corresponding to
      ! a boundary segment.  A boundary region roughly contains
      ! the type, the min/max parameter value and whether the
      ! endpoints are inside the region or not.
      call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND

      ! Real part of the solution
      call bcasm_newDirichletBConRealBD(&
          rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues_Real,rproblem%rcollection)

      ! Imaginary part
      call bcasm_newDirichletBConRealBD(&
          rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,2,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues_Aimag,rproblem%rcollection)

#else
#error 'Test case is undefined.' 
#endif

    case (POISSON_SYSTEM,SSE_SYSTEM1,SSE_SYSTEM2)
      ! No essential boundary conditions

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscreteBC")
      call sys_halt()
    end select
  end do

  end subroutine

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

  end subroutine

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

  ! A pointer to the system matrix and the RHS vector as well as
  ! the discretisation
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_vectorBlock), pointer :: p_rrhs,p_rvector
  type(t_vectorBlock), target :: rvecTmp,rvecTmp1,rvecTmp2,rvecTmp3,rrhsTmp2

  ! A solver node that accepts parameters for the linear solver
  type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
  type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices

  ! One level of multigrid
  type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo

  ! Relaxation parameter
  real(DP), parameter :: domega = 0.1_DP

  ! Iteration counter
  integer :: ite

  ! Norm of residuals
  real(DP) :: dresNorm0,dresNorm


  type(t_vectorScalar) :: rvectorSc,rvecTmpSc,rrhsSc
  type(t_vectorBlock) :: rvectorBl,rvecTmpBl,rrhsBl

  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR,SSE_SCALAR)
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

  case (POISSON_SYSTEM+4711)
    !---------------------------------------------------------------------------

    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax

    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix

!!$    ! Set up prolongation matrices
!!$    do i=ilvmin+1,ilvmax
!!$      do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents
!!$
!!$        ! Create the matrix structure of the prolongation matrix.
!!$        call mlop_create2LvlMatrixStruct(&
!!$            rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(j),&
!!$            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
!!$            LSYSSC_MATRIX9, rproblem%RlevelInfo(i)%rmatProl(j))
!!$
!!$        ! Assemble the entries of the prolongation matrix.
!!$        call mlop_build2LvlProlMatrix(&
!!$            rproblem%RlevelInfo(i-1)%rdiscretisation%RspatialDiscr(j),&
!!$            rproblem%RlevelInfo(i)%rdiscretisation%RspatialDiscr(j),&
!!$            .true., rproblem%RlevelInfo(i)%rmatProl(j),&
!!$            rcubatureInfoCoarse=rproblem%RlevelInfo(i-1)%RcubatureInfo(j),&
!!$            rcubatureInfoFine=rproblem%RlevelInfo(i)%RcubatureInfo(j))
!!$
!!$        ! Assemble the entries of the restriction matrix.
!!$        call lsyssc_transposeMatrix(rproblem%RlevelInfo(i)%rmatProl(j),&
!!$            rproblem%RlevelInfo(i)%rmatRest(j),LSYSSC_TR_ALL)
!!$      end do
!!$
!!$      ! Now set up an interlevel projecton structure for this level
!!$      ! based on the system matrix on this level.
!!$      call mlprj_initProjectionMat(rproblem%RlevelInfo(i)%rprojection,&
!!$          rproblem%RlevelInfo(i)%rmatrix)
!!$
!!$      ! Initialise the matrix-based projection
!!$      do j=1,rproblem%RlevelInfo(i)%rdiscretisation%ncomponents
!!$        call mlprj_initMatrixProjection(&
!!$            rproblem%RlevelInfo(i)%rprojection%RscalarProjection(1,j),&
!!$            rproblem%RlevelInfo(i)%rmatProl(j),&
!!$            rmatrixRest=rproblem%RlevelInfo(i)%rmatRest(j))
!!$      end do
!!$    end do

!!$    ! Set up an interlevel projecton structure for the coarse-most level.
!!$    call mlprj_initProjectionMat(rproblem%RlevelInfo(ilvmin)%rprojection,&
!!$        rproblem%RlevelInfo(ilvmin)%rmatrix)
!!$
!!$    ! Create temporary vectors we need that for some preparation.
!!$    call lsysbl_deriveSubvector(p_rrhs, rrhsTmp2, 2, 3, .false.)
!!$    call lsysbl_deriveSubvector(p_rrhs, rvecTmp3, 2, 3, .false.)
!!$    call lsysbl_deriveSubvector(p_rvector, rvecTmp1, 1, 1, .true.)
!!$    call lsysbl_deriveSubvector(p_rvector, rvecTmp2, 2, 3, .true.)
!!$
!!$    call lsysbl_createVector(p_rrhs, rvecTmp, .false.)

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2(p_rsolverNode,ilvmax-ilvmin+1)

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
        call linsol_initUMFPACK4(p_rcoarseGridSolver)
!!$        ! Setting up a GMRES coarse grid solver
!!$        call linsol_initGMRES(p_rcoarseGridSolver,20,&
!!$              Rfilter=rproblem%RlevelInfo(i)%RfilterChain)
      else
        ! Set up GMRES smoother for multigrid with damping parameter
        ! 0.7, 4 smoothing steps:
        call linsol_initGMRES(p_rsmoother,20,&
            Rfilter=rproblem%RlevelInfo(i)%RfilterChain)
        call linsol_convertToSmoother(p_rsmoother,4,0.7_DP)
      end if

      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level(p_rsolverNode,i-ilvmin+1,p_rlevelInfo)
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
      p_rlevelInfo%p_rpresmoother      => p_rsmoother
      p_rlevelInfo%p_rpostsmoother     => p_rsmoother

      ! Attach the filter chain which imposes boundary conditions on that level.
      p_rlevelInfo%p_RfilterChain => rproblem%RlevelInfo(i)%RfilterChain

!!$      ! Attach our user-defined projection to the level.
!!$      call linsol_initProjMultigrid2Level(p_rlevelInfo,&
!!$          rproblem%RlevelInfo(i)%rprojection)
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

  case(POISSON_SYSTEM,SSE_SYSTEM1,SSE_SYSTEM2)
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

  end subroutine

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

  ! A temporal block discretisation structure
  type(t_blockDiscretisation) :: rblockDiscr

  ! A temporal vector to store the recovered gradient
  type(t_vectorBlock), target :: rvectorBlock,rvectorBlockX,rvectorBlockY
  type(t_vectorBlock), target :: rvectorBlock_Real,rvectorBlock_Aimag
  type(t_vectorBlock), target :: rvectorBlockX_Real,rvectorBlockX_Aimag
  type(t_vectorBlock), target :: rvectorBlockY_Real,rvectorBlockY_Aimag

  ! Pointer to gradient components
  type(t_vectorScalar), pointer :: p_rvectorDerivX,p_rvectorDerivY
  type(t_vectorScalar), pointer :: p_rvectorDerivX_Real,p_rvectorDerivY_Real
  type(t_vectorScalar), pointer :: p_rvectorDerivX_Aimag,p_rvectorDerivY_Aimag

  ! Output block for UCD output to VTK file
  type(t_ucdExport) :: rexport
  character(len=SYS_STRLEN) :: sucddir,sucdfile

  ! Error of FE function to reference function
  real(DP) :: derror

  ! Number of DOFs
  integer :: i,ndof,iucdtype

  ! Method for calculating the gradient
  integer :: cgradType,cgradSubType

  real(DP), dimension(:), pointer :: p_SSE_RE,p_SSE_IM
  real(DP), dimension(:,:), pointer :: p_Dcoords

  complex(DP) :: cC,cvalue
  real(DP) :: dalpha,dh,ds,dAv,dr1,dr2
  integer :: ivt

  ! Get method for calculating the gradient
  call parlst_getvalue_int(rparlist, '', 'GRADTYPE', cgradType)
  call parlst_getvalue_int(rparlist, '', 'GRADSUBTYPE', cgradSubType)

  ! Calculate the total number of DoF's
  ndof = 0
  do i=1,rproblem%rvector%nblocks
    ndof = ndof + dof_igetNDofGlob(rproblem%rvector%RvectorBlock(i)%p_rspatialDiscr)
  end do

  ! Add number of cells and vertices
  call ctab_addValue(rtable, "cells",&
      rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL)
  call ctab_addValue(rtable, "dofs", ndof)

  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR,POISSON_SYSTEM)

    ! --- solution -------------------------------------------------------------

    ! Calculate the error to the reference function.
    call pperr_scalar(PPERR_L2ERROR,derror,rproblem%rvector%RvectorBlock(1),&
        getReferenceFunction_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error u", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rproblem%rvector%RvectorBlock(1),&
        getReferenceFunction_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error u", derror)

    ! --- first derivative -----------------------------------------------------

    select case(rproblem%cproblemtype)
    case (POISSON_SCALAR)
      ! Recover gradient by superconvergent patch recovery
      call spdiscr_initBlockDiscr(rblockDiscr,2,&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
          rproblem%rvector%p_rblockDiscr%p_rboundary)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(1), .true.)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(2), .true.)

      call lsysbl_createVector(rblockDiscr, rvectorBlock, .false.)
      call ppgrd_calcGradient(rproblem%rvector%RvectorBlock(1), rvectorBlock,&
          cgradType, cgradSubType)
      p_rvectorDerivX => rvectorBlock%RvectorBlock(1)
      p_rvectorDerivY => rvectorBlock%RvectorBlock(2)

    case (POISSON_SYSTEM)
      ! Set pointer to scalar solution components
      p_rvectorDerivX => rproblem%rvector%RvectorBlock(2)
      p_rvectorDerivY => rproblem%rvector%RvectorBlock(3)
    end select

    ! Calculate the error to the reference DerivX.
    call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivX,&
        getReferenceDerivX_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error u_x", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivX,&
        getReferenceDerivX_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error u_x", derror)

    ! Calculate the error to the reference DerivY.
    call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivY,&
        getReferenceDerivY_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error u_y", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivY,&
        getReferenceDerivY_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error u_y", derror)

    ! --- second derivative ----------------------------------------------------

    select case(rproblem%cproblemtype)
    case (POISSON_SYSTEM)
      ! Recover gradient by superconvergent patch recovery
      call spdiscr_initBlockDiscr(rblockDiscr,2,&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
          rproblem%rvector%p_rblockDiscr%p_rboundary)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(1), .true.)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(2), .true.)
    end select

    ! Recover second derivative by superconvergent patch recovery
    call lsysbl_createVector(rblockDiscr, rvectorBlockX, .false.)
    call lsysbl_createVector(rblockDiscr, rvectorBlockY, .false.)
    call ppgrd_calcGradient(p_rvectorDerivX, rvectorBlockX,&
        cgradType, cgradSubType)
    call ppgrd_calcGradient(p_rvectorDerivY, rvectorBlockY,&
        cgradType, cgradSubType)

    ! Calculate the error to the reference DerivXX.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX%RvectorBlock(1),&
        getReferenceDerivXX_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error u_xx", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX%RvectorBlock(1),&
        getReferenceDerivXX_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error u_xx", derror)

    ! Calculate the error to the reference DerivXY.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX%RvectorBlock(2),&
        getReferenceDerivXY_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error u_xy", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX%RvectorBlock(2),&
        getReferenceDerivXY_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error u_xy", derror)

    ! Calculate the error to the reference DerivYX.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY%RvectorBlock(1),&
        getReferenceDerivYX_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error u_yx", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY%RvectorBlock(1),&
        getReferenceDerivXY_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error u_yx", derror)

    ! Calculate the error to the reference DerivYY.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY%RvectorBlock(2),&
        getReferenceDerivYY_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error u_yy", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY%RvectorBlock(2),&
        getReferenceDerivYY_Poisson, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error u_yy", derror)
    
    ! Start UCD export to file:
    if (.not. sys_getenv_string("UCDDIR", sucddir))&
        call parlst_getvalue_string(rparlist, '', 'SUCDDIR', sucddir, './ucd')
    call parlst_getvalue_string(rparlist, '', 'UCDFILE', sucdfile, '')
    call parlst_getvalue_int(rparlist, '', 'UCDTYPE', iucdtype, 0)

    if ((trim(adjustl(sucdfile)) .ne. '') .and.&
        (iucdtype .ne. UCD_FORMAT_NONE)) then

      select case(iucdtype)
      case(UCD_FORMAT_GMV)
        call ucd_startGMV(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".gmv")
      case(UCD_FORMAT_BGMV)
        call ucd_startBGMV(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".vtk")
      case(UCD_FORMAT_AVS)
        call ucd_startAVS(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".vtk")
      case(UCD_FORMAT_VTK)
        call ucd_startVTK(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".vtk")
      case default
        call output_line("Invalid type of UCD output file..", &
            OU_CLASS_ERROR,OU_MODE_STD,"sse_postprocessing")
        call sys_halt()
      end select

      ! Add the solution and its (recovered) gradient to the UCD exporter
      call ucd_addVectorByVertex(rexport, "u", UCD_VAR_STANDARD, &
          rproblem%rvector%RvectorBlock(1))
      call ucd_addVectorFieldByVertex(rexport, "grad u", UCD_VAR_STANDARD, &
          (/p_rvectorDerivX,p_rvectorDerivY/))
      call ucd_addVectorFieldByVertex(rexport, "grad u_x", UCD_VAR_STANDARD, &
          (/rvectorBlockX%RvectorBlock(1),rvectorBlockX%RvectorBlock(2)/))
      call ucd_addVectorFieldByVertex(rexport, "grad u_y", UCD_VAR_STANDARD, &
          (/rvectorBlockY%RvectorBlock(1),rvectorBlockY%RvectorBlock(2)/))

      ! Write the file to disc, that is it.
      call ucd_write(rexport)
      call ucd_release(rexport)
    end if

    ! Clean temporal structures
    call lsysbl_releaseVector(rvectorBlock)
    call lsysbl_releaseVector(rvectorBlockX)
    call lsysbl_releaseVector(rvectorBlockY)
    call spdiscr_releaseBlockDiscr(rblockDiscr)

  case (SSE_SCALAR,SSE_SYSTEM1,SSE_SYSTEM2)

    ! --- solution -------------------------------------------------------------

    ! Calculate the error to the reference function.
    call pperr_scalar(PPERR_L2ERROR,derror,rproblem%rvector%RvectorBlock(1),&
        getReferenceFunction_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Re(SSE)", derror)
    call pperr_scalar(PPERR_L2ERROR,derror,rproblem%rvector%RvectorBlock(2),&
        getReferenceFunction_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Im(SSE)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rproblem%rvector%RvectorBlock(1),&
        getReferenceFunction_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Re(SSE)", derror)
    call pperr_scalar(PPERR_H1ERROR,derror,rproblem%rvector%RvectorBlock(2),&
        getReferenceFunction_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Im(SSE)", derror)

    ! --- first derivative -----------------------------------------------------

    select case(rproblem%cproblemtype)
    case (SSE_SCALAR)
      ! Recover gradient by superconvergent patch recovery
      call spdiscr_initBlockDiscr(rblockDiscr,2,&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
          rproblem%rvector%p_rblockDiscr%p_rboundary)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(1), .true.)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(2), .true.)

      call lsysbl_createVector(rblockDiscr, rvectorBlock_Real, .false.)
      call ppgrd_calcGradient(rproblem%rvector%RvectorBlock(1),&
          rvectorBlock_Real, cgradType, cgradSubType)
      p_rvectorDerivX_Real => rvectorBlock_Real%RvectorBlock(1)
      p_rvectorDerivY_Real => rvectorBlock_Real%RvectorBlock(2)

      call lsysbl_createVector(rblockDiscr, rvectorBlock_Aimag, .false.)
      call ppgrd_calcGradient(rproblem%rvector%RvectorBlock(2),&
          rvectorBlock_Aimag, cgradType, cgradSubType)
      p_rvectorDerivX_Aimag => rvectorBlock_Aimag%RvectorBlock(1)
      p_rvectorDerivY_Aimag => rvectorBlock_Aimag%RvectorBlock(2)

    case (SSE_SYSTEM1,SSE_SYSTEM2)
      ! Set pointer to scalar solution components
      p_rvectorDerivX_Real  => rproblem%rvector%RvectorBlock(3)
      p_rvectorDerivX_Aimag => rproblem%rvector%RvectorBlock(4)
      p_rvectorDerivY_Real  => rproblem%rvector%RvectorBlock(5)
      p_rvectorDerivY_Aimag => rproblem%rvector%RvectorBlock(6)
    end select

    ! Calculate the error to the reference DerivX.
    call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivX_Real,&
        getReferenceDerivX_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Re(SSE_x)", derror)

    call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivX_Aimag,&
        getReferenceDerivX_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Im(SSE_x)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivX_Real,&
        getReferenceDerivX_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Re(SSE_x)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivX_Aimag,&
        getReferenceDerivX_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Im(SSE_x)", derror)

    ! Calculate the error to the reference DerivY.
    call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivY_Real,&
        getReferenceDerivY_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Re(SSE_y)", derror)

    call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivY_Aimag,&
        getReferenceDerivY_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Im(SSE_y)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivY_Real,&
        getReferenceDerivY_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Re(SSE_y)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivY_Aimag,&
        getReferenceDerivY_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Im(SSE_y)", derror)

    ! --- second derivative ----------------------------------------------------

    select case(rproblem%cproblemtype)
    case (SSE_SYSTEM1,SSE_SYSTEM2)
      ! Recover gradient by superconvergent patch recovery
      call spdiscr_initBlockDiscr(rblockDiscr,2,&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
          rproblem%rvector%p_rblockDiscr%p_rboundary)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(1), .true.)
      call spdiscr_duplicateDiscrSc(&
          rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(2), .true.)
    end select

    ! Recover second derivative by superconvergent patch recovery
    call lsysbl_createVector(rblockDiscr, rvectorBlockX_Real, .false.)
    call lsysbl_createVector(rblockDiscr, rvectorBlockX_Aimag, .false.)
    call lsysbl_createVector(rblockDiscr, rvectorBlockY_Real, .false.)
    call lsysbl_createVector(rblockDiscr, rvectorBlockY_Aimag, .false.)
    call ppgrd_calcGradient(p_rvectorDerivX_Real, rvectorBlockX_Real,&
        cgradType, cgradSubType)
    call ppgrd_calcGradient(p_rvectorDerivX_Aimag, rvectorBlockX_Aimag,&
        cgradType, cgradSubType)
    call ppgrd_calcGradient(p_rvectorDerivY_Real, rvectorBlockY_Real,&
        cgradType, cgradSubType)
    call ppgrd_calcGradient(p_rvectorDerivY_Aimag, rvectorBlockY_Aimag,&
        cgradType, cgradSubType)

    ! Calculate the error to the reference DerivXX.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX_Real%RvectorBlock(1),&
        getReferenceDerivXX_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Re(SSE_xx)", derror)

    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX_Aimag%RvectorBlock(1),&
        getReferenceDerivXX_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Im(SSE_xx)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX_Real%RvectorBlock(1),&
        getReferenceDerivXX_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Re(SSE_xx)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX_Aimag%RvectorBlock(1),&
        getReferenceDerivXX_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Im(SSE_xx)", derror)

    ! Calculate the error to the reference DerivXY.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX_Real%RvectorBlock(2),&
        getReferenceDerivXY_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Re(SSE_xy)", derror)

    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX_Aimag%RvectorBlock(2),&
        getReferenceDerivXY_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Im(SSE_xy)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX_Real%RvectorBlock(2),&
        getReferenceDerivXY_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Re(SSE_xy)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX_Aimag%RvectorBlock(2),&
        getReferenceDerivXY_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Im(SSE_xy)", derror)

    ! Calculate the error to the reference DerivYX.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY_Real%RvectorBlock(1),&
        getReferenceDerivYX_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Re(SSE_yx)", derror)

    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY_Aimag%RvectorBlock(1),&
        getReferenceDerivYX_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Im(SSE_yx)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY_Real%RvectorBlock(1),&
        getReferenceDerivXY_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Re(SSE_yx)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY_Aimag%RvectorBlock(1),&
        getReferenceDerivXY_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Im(SSE_yx)", derror)

    ! Calculate the error to the reference DerivYY.
    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY_Real%RvectorBlock(2),&
        getReferenceDerivYY_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Re(SSE_yy)", derror)

    call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY_Aimag%RvectorBlock(2),&
        getReferenceDerivYY_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "L2-error Im(SSE_yy)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY_Real%RvectorBlock(2),&
        getReferenceDerivYY_Real, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Re(SSE_yy)", derror)

    call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY_Aimag%RvectorBlock(2),&
        getReferenceDerivYY_Aimag, rcubatureInfo=&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
    call ctab_addValue(rtable, "H1-error Im(SSE_yy)", derror)

    ! Start UCD export to file:
    if (.not. sys_getenv_string("UCDDIR", sucddir))&
        call parlst_getvalue_string(rparlist, '', 'SUCDDIR', sucddir, './ucd')
    call parlst_getvalue_string(rparlist, '', 'UCDFILE', sucdfile, '')
    call parlst_getvalue_int(rparlist, '', 'UCDTYPE', iucdtype, 0)

    if ((trim(adjustl(sucdfile)) .ne. '') .and.&
        (iucdtype .ne. UCD_FORMAT_NONE)) then

      select case(iucdtype)
      case(UCD_FORMAT_GMV)
        call ucd_startGMV(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".gmv")
      case(UCD_FORMAT_BGMV)
        call ucd_startBGMV(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".vtk")
      case(UCD_FORMAT_AVS)
        call ucd_startAVS(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".vtk")
      case(UCD_FORMAT_VTK)
        call ucd_startVTK(rexport,UCD_FLAG_STANDARD,&
            rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
            trim(sucddir)//"/"//trim(sucdfile)//"_"//&
            trim(sys_siL(rproblem%ilvmax,2))//".vtk")
      case default
        call output_line("Invalid type of UCD output file..", &
            OU_CLASS_ERROR,OU_MODE_STD,"sse_postprocessing")
        call sys_halt()
      end select
      
      call lsyssc_getbase_double(rproblem%rvector%RvectorBlock(1), p_SSE_RE)
      call lsyssc_getbase_double(rproblem%rvector%RvectorBlock(2), p_SSE_IM)
      call storage_getbase_double2D(&
          rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation%h_DvertexCoords,&
          p_Dcoords)

      do ivt=1,rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation%NVT
        dh = sse_bottomProfile(p_Dcoords(1,ivt),p_Dcoords(2,ivt))
        ds = sse_bottomStress(p_Dcoords(1,ivt),p_Dcoords(2,ivt))
        dAv = sse_eddyViscosity(p_Dcoords(1,ivt),p_Dcoords(2,ivt))

        ! Compute coefficients
        dalpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*dalpha**3) * &
            ((ds*sin(dalpha*dh))/(dalpha*dAv*sin(dalpha*dh)-ds*cos(dalpha*dh)) + dh*dalpha)
        dr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0/dlengthB**2 - 4.0_DP * cimg*dtidalfreq/cC))
        dr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0/dlengthB**2 - 4.0_DP * cimg*dtidalfreq/cC))

        cvalue = (dr1*exp(dr1*dlength+dr2*p_Dcoords(1,ivt)) -&
                  dr2*exp(dr1*p_Dcoords(1,ivt)+dr2*dlength))/&
                 (dr1*exp(dr1*dlength)-dr2*exp(dr2*dlength))

        p_SSE_RE(ivt) = p_SSE_RE(ivt)-real(cvalue)
        p_SSE_IM(ivt) = p_SSE_IM(ivt)-aimag(cvalue)
        
      end do

      ! Add the solution and its (recovered) gradient to the UCD exporter
      call ucd_addVectorByVertex(rexport, "Re(SSE)", UCD_VAR_STANDARD, &
          rproblem%rvector%RvectorBlock(1))
      call ucd_addVectorByVertex(rexport, "Im(SSE)", UCD_VAR_STANDARD, &
          rproblem%rvector%RvectorBlock(2))
      call ucd_addVectorFieldByVertex(rexport, "Re(grad SSE)", UCD_VAR_STANDARD, &
          (/p_rvectorDerivX_Real,p_rvectorDerivY_Real/))
      call ucd_addVectorFieldByVertex(rexport, "Im(grad SSE)", UCD_VAR_STANDARD, &
          (/p_rvectorDerivX_Aimag,p_rvectorDerivY_Aimag/))
      call ucd_addVectorFieldByVertex(rexport, "Re(grad SSE_x)", UCD_VAR_STANDARD, &
          (/rvectorBlockX_Real%RvectorBlock(1),rvectorBlockX_Real%RvectorBlock(2)/))
      call ucd_addVectorFieldByVertex(rexport, "Im(grad SSE_x)", UCD_VAR_STANDARD, &
          (/rvectorBlockX_Aimag%RvectorBlock(1),rvectorBlockX_Aimag%RvectorBlock(2)/))
      call ucd_addVectorFieldByVertex(rexport, "Re(grad SSE_y)", UCD_VAR_STANDARD, &
          (/rvectorBlockY_Real%RvectorBlock(1),rvectorBlockY_Real%RvectorBlock(2)/))
      call ucd_addVectorFieldByVertex(rexport, "Im(grad SSE_y)", UCD_VAR_STANDARD, &
          (/rvectorBlockY_Aimag%RvectorBlock(1),rvectorBlockY_Aimag%RvectorBlock(2)/))
      
      ! Write the file to disc, that is it.
      call ucd_write(rexport)
      call ucd_release(rexport)
    end if

    ! Clean temporal structures
    call lsysbl_releaseVector(rvectorBlock_Real)
    call lsysbl_releaseVector(rvectorBlock_Aimag)
    call lsysbl_releaseVector(rvectorBlockX_Real)
    call lsysbl_releaseVector(rvectorBlockX_Aimag)
    call lsysbl_releaseVector(rvectorBlockY_Real)
    call lsysbl_releaseVector(rvectorBlockY_Aimag)
    call spdiscr_releaseBlockDiscr(rblockDiscr)

  case default
    call output_line("Invalid type of problem.", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_postprocessing")
    call sys_halt()
  end select

  end subroutine

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

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sse_doneBC(rproblem)

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

  end subroutine


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

  end subroutine

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

  do i=rproblem%ilvmax,rproblem%ilvmin,-1
    ! Release the triangulation
    call tria_done(rproblem%RlevelInfo(i)%rtriangulation)
  end do

  deallocate(rproblem%RlevelInfo)

  ! Finally release the domain.
  call boundary_release(rproblem%rboundary)

  end subroutine

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

  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR,POISSON_SYSTEM)

    ! Compute reduction rates
    call ctab_evalConvergenceRate(rtable,"L2-error u",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error u_x",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error u_y",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error u_xx",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error u_xy",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error u_yx",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error u_yy",CTAB_REDUCTION_RATE)

    call ctab_evalConvergenceRate(rtable,"H1-error u",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error u_x",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error u_y",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error u_xx",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error u_xy",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error u_yx",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error u_yy",CTAB_REDUCTION_RATE)

    ! Adjust format of convergence table
    call ctab_setPrecision(rtable,"L2-error u",3)
    call ctab_setPrecision(rtable,"L2-error u-convrate",3)
    call ctab_setPrecision(rtable,"L2-error u_x",3)
    call ctab_setPrecision(rtable,"L2-error u_y",3)
    call ctab_setPrecision(rtable,"L2-error u_x-convrate",3)
    call ctab_setPrecision(rtable,"L2-error u_y-convrate",3)
    call ctab_setPrecision(rtable,"L2-error u_xx",3)
    call ctab_setPrecision(rtable,"L2-error u_xy",3)
    call ctab_setPrecision(rtable,"L2-error u_yx",3)
    call ctab_setPrecision(rtable,"L2-error u_yy",3)
    call ctab_setPrecision(rtable,"L2-error u_xx-convrate",3)
    call ctab_setPrecision(rtable,"L2-error u_xy-convrate",3)
    call ctab_setPrecision(rtable,"L2-error u_yx-convrate",3)
    call ctab_setPrecision(rtable,"L2-error u_yy-convrate",3)

    call ctab_setPrecision(rtable,"H1-error u",3)
    call ctab_setPrecision(rtable,"H1-error u-convrate",3)
    call ctab_setPrecision(rtable,"H1-error u_x",3)
    call ctab_setPrecision(rtable,"H1-error u_y",3)
    call ctab_setPrecision(rtable,"H1-error u_x-convrate",3)
    call ctab_setPrecision(rtable,"H1-error u_y-convrate",3)
    call ctab_setPrecision(rtable,"H1-error u_xx",3)
    call ctab_setPrecision(rtable,"H1-error u_xy",3)
    call ctab_setPrecision(rtable,"H1-error u_yx",3)
    call ctab_setPrecision(rtable,"H1-error u_yy",3)
    call ctab_setPrecision(rtable,"H1-error u_xx-convrate",3)
    call ctab_setPrecision(rtable,"H1-error u_xy-convrate",3)
    call ctab_setPrecision(rtable,"H1-error u_yx-convrate",3)
    call ctab_setPrecision(rtable,"H1-error u_yy-convrate",3)

    ! Set scientific flag
    call ctab_setScientific(rtable,"L2-error u",.true.)
    call ctab_setScientific(rtable,"L2-error u_x",.true.)
    call ctab_setScientific(rtable,"L2-error u_y",.true.)
    call ctab_setScientific(rtable,"L2-error u_xx",.true.)
    call ctab_setScientific(rtable,"L2-error u_xy",.true.)
    call ctab_setScientific(rtable,"L2-error u_yx",.true.)
    call ctab_setScientific(rtable,"L2-error u_yy",.true.)

    call ctab_setScientific(rtable,"H1-error u",.true.)
    call ctab_setScientific(rtable,"H1-error u_x",.true.)
    call ctab_setScientific(rtable,"H1-error u_y",.true.)
    call ctab_setScientific(rtable,"H1-error u_xx",.true.)
    call ctab_setScientific(rtable,"H1-error u_xy",.true.)
    call ctab_setScientific(rtable,"H1-error u_yx",.true.)
    call ctab_setScientific(rtable,"H1-error u_yy",.true.)

    ! Set Tex captions
    call ctab_setTexCaption(rtable,"cells","\# cells")
    call ctab_setTexCaption(rtable,"dofs","\# dofs")

    call ctab_setTexCaption(rtable,"L2-error u","$L^2(u)$")
    call ctab_setTexCaption(rtable,"H1-error u","$H^1(u)$")
    select case(rproblem%cproblemtype)
    case (POISSON_SCALAR)
      call ctab_setTexCaption(rtable,"L2-error u_x","$L^2(\partial_{x}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"L2-error u_y","$L^2(\partial_{y}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"L2-error u_xx","$L^2(\partial_{xx}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"L2-error u_xy","$L^2(\partial_{xy}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"L2-error u_yx","$L^2(\partial_{yx}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"L2-error u_yy","$L^2(\partial_{yy}^{\rm ZZ}u)$")

      call ctab_setTexCaption(rtable,"H1-error u_x","$H^1(\partial_{x}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"H1-error u_y","$H^1(\partial_{y}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"H1-error u_xx","$H^1(\partial_{xx}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"H1-error u_xy","$H^1(\partial_{xy}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"H1-error u_yx","$H^1(\partial_{yx}^{\rm ZZ}u)$")
      call ctab_setTexCaption(rtable,"H1-error u_yy","$H^1(\partial_{yy}^{\rm ZZ}u)$")

    case (POISSON_SYSTEM)
      call ctab_setTexCaption(rtable,"L2-error u_x","$L^2(\sigma_x)$")
      call ctab_setTexCaption(rtable,"L2-error u_y","$L^2(\sigma_y)$")
      call ctab_setTexCaption(rtable,"L2-error u_xx","$L^2(\partial_{x}^{\rm ZZ}\sigma_x)$")
      call ctab_setTexCaption(rtable,"L2-error u_xy","$L^2(\partial_{y}^{\rm ZZ}\sigma_x)$")
      call ctab_setTexCaption(rtable,"L2-error u_yx","$L^2(\partial_{x}^{\rm ZZ}\sigma_y)$")
      call ctab_setTexCaption(rtable,"L2-error u_yy","$L^2(\partial_{y}^{\rm ZZ}\sigma_y)$")

      call ctab_setTexCaption(rtable,"H1-error u_x","$H^1(\sigma_x)$")
      call ctab_setTexCaption(rtable,"H1-error u_y","$H^1(\sigma_y)$")
      call ctab_setTexCaption(rtable,"H1-error u_xx","$H^1(\partial_{x}^{\rm ZZ}\sigma_x)$")
      call ctab_setTexCaption(rtable,"H1-error u_xy","$H^1(\partial_{y}^{\rm ZZ}\sigma_x)$")
      call ctab_setTexCaption(rtable,"H1-error u_yx","$H^1(\partial_{x}^{\rm ZZ}\sigma_y)$")
      call ctab_setTexCaption(rtable,"H1-error u_yy","$H^1(\partial_{y}^{\rm ZZ}\sigma_y)$")
    end select

    call ctab_setTexCaption(rtable,"L2-error u-convrate","")
    call ctab_setTexCaption(rtable,"L2-error u_x-convrate","")
    call ctab_setTexCaption(rtable,"L2-error u_y-convrate","")
    call ctab_setTexCaption(rtable,"L2-error u_xx-convrate","")
    call ctab_setTexCaption(rtable,"L2-error u_xy-convrate","")
    call ctab_setTexCaption(rtable,"L2-error u_yx-convrate","")
    call ctab_setTexCaption(rtable,"L2-error u_yy-convrate","")

    call ctab_setTexCaption(rtable,"H1-error u-convrate","")
    call ctab_setTexCaption(rtable,"H1-error u_x-convrate","")
    call ctab_setTexCaption(rtable,"H1-error u_y-convrate","")
    call ctab_setTexCaption(rtable,"H1-error u_xx-convrate","")
    call ctab_setTexCaption(rtable,"H1-error u_xy-convrate","")
    call ctab_setTexCaption(rtable,"H1-error u_yx-convrate","")
    call ctab_setTexCaption(rtable,"H1-error u_yy-convrate","")

    ! Set Tex format
    call ctab_setTexFormat(rtable,"cells","r")
    call ctab_setTexFormat(rtable,"dofs","r")

    ! Hide all H1-columns
    call ctab_setHidden(rtable,"H1-error u",.true.)
    call ctab_setHidden(rtable,"H1-error u-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error u_x",.true.)
    call ctab_setHidden(rtable,"H1-error u_y",.true.)
    call ctab_setHidden(rtable,"H1-error u_x-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error u_y-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error u_xx",.true.)
    call ctab_setHidden(rtable,"H1-error u_xy",.true.)
    call ctab_setHidden(rtable,"H1-error u_yx",.true.)
    call ctab_setHidden(rtable,"H1-error u_yy",.true.)
    call ctab_setHidden(rtable,"H1-error u_xx-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error u_xy-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error u_yx-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error u_yy-convrate",.true.)

    select case(rproblem%cproblemtype)
    case (POISSON_SCALAR)
      call ctab_setTexTableCaption(rtable,&
          "$L^2$-Convergence table: Poisson problem solved as scalar equation.")

    case (POISSON_SYSTEM)
      call ctab_setTexTableCaption(rtable,&
          "$L^2$-Convergence table: Poisson problem solved in mixed formulation.")
    end select
    call ctab_setTexTableLabel(rtable,"tab:l2_convergence_rate")

    ! Write convergence table to Tex file
    call ctab_outputTex(rtable,'./table_l2.tex')

    select case(rproblem%cproblemtype)
    case (POISSON_SCALAR)
      call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table: Poisson problem solved as scalar equation.")
    case (POISSON_SYSTEM)
      call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table: Poisson problem solved in mixed formulation.")
    end select
    call ctab_setTexTableLabel(rtable,"tab:h1_convergence_rate")

    ! Unhide all H1-columns
    call ctab_setHidden(rtable,"H1-error u",.false.)
    call ctab_setHidden(rtable,"H1-error u-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error u_x",.false.)
    call ctab_setHidden(rtable,"H1-error u_y",.false.)
    call ctab_setHidden(rtable,"H1-error u_x-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error u_y-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error u_xx",.false.)
    call ctab_setHidden(rtable,"H1-error u_xy",.false.)
    call ctab_setHidden(rtable,"H1-error u_yx",.false.)
    call ctab_setHidden(rtable,"H1-error u_yy",.false.)
    call ctab_setHidden(rtable,"H1-error u_xx-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error u_xy-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error u_yx-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error u_yy-convrate",.false.)

    ! Hide all L2-columns
    call ctab_setHidden(rtable,"L2-error u",.true.)
    call ctab_setHidden(rtable,"L2-error u-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error u_x",.true.)
    call ctab_setHidden(rtable,"L2-error u_y",.true.)
    call ctab_setHidden(rtable,"L2-error u_x-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error u_y-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error u_xx",.true.)
    call ctab_setHidden(rtable,"L2-error u_xy",.true.)
    call ctab_setHidden(rtable,"L2-error u_yx",.true.)
    call ctab_setHidden(rtable,"L2-error u_yy",.true.)
    call ctab_setHidden(rtable,"L2-error u_xx-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error u_xy-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error u_yx-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error u_yy-convrate",.true.)

    ! Write convergence table to Tex file
    call ctab_outputTex(rtable,'./table_h1.tex')

  case (SSE_SCALAR,SSE_SYSTEM1,SSE_SYSTEM2)

    ! Compute reduction rates
    call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_x)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_y)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_xx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_xy)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_yx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_yy)",CTAB_REDUCTION_RATE)

    call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_x)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_y)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_xx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_xy)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_yx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_yy)",CTAB_REDUCTION_RATE)

    call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_x)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_y)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_xx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_xy)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_yx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_yy)",CTAB_REDUCTION_RATE)

    call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_x)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_y)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_xx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_xy)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_yx)",CTAB_REDUCTION_RATE)
    call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_yy)",CTAB_REDUCTION_RATE)

    ! Adjust format of convergence table
    call ctab_setPrecision(rtable,"L2-error Re(SSE)",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_x)",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_y)",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_x)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_y)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_xx)",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_xy)",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_yx)",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_yy)",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_xx)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_xy)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_yx)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Re(SSE_yy)-convrate",3)

    call ctab_setPrecision(rtable,"L2-error Im(SSE)",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_x)",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_y)",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_x)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_y)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_xx)",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_xy)",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_yx)",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_yy)",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_xx)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_xy)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_yx)-convrate",3)
    call ctab_setPrecision(rtable,"L2-error Im(SSE_yy)-convrate",3)

    call ctab_setPrecision(rtable,"H1-error Re(SSE)",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_x)",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_y)",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_x)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_y)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_xx)",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_xy)",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_yx)",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_yy)",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_xx)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_xy)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_yx)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Re(SSE_yy)-convrate",3)

    call ctab_setPrecision(rtable,"H1-error Im(SSE)",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_x)",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_y)",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_x)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_y)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_xx)",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_xy)",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_yx)",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_yy)",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_xx)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_xy)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_yx)-convrate",3)
    call ctab_setPrecision(rtable,"H1-error Im(SSE_yy)-convrate",3)

    ! Set scientific flag
    call ctab_setScientific(rtable,"L2-error Re(SSE)",.true.)
    call ctab_setScientific(rtable,"L2-error Re(SSE_x)",.true.)
    call ctab_setScientific(rtable,"L2-error Re(SSE_y)",.true.)
    call ctab_setScientific(rtable,"L2-error Re(SSE_xx)",.true.)
    call ctab_setScientific(rtable,"L2-error Re(SSE_xy)",.true.)
    call ctab_setScientific(rtable,"L2-error Re(SSE_yx)",.true.)
    call ctab_setScientific(rtable,"L2-error Re(SSE_yy)",.true.)

    call ctab_setScientific(rtable,"L2-error Im(SSE)",.true.)
    call ctab_setScientific(rtable,"L2-error Im(SSE_x)",.true.)
    call ctab_setScientific(rtable,"L2-error Im(SSE_y)",.true.)
    call ctab_setScientific(rtable,"L2-error Im(SSE_xx)",.true.)
    call ctab_setScientific(rtable,"L2-error Im(SSE_xy)",.true.)
    call ctab_setScientific(rtable,"L2-error Im(SSE_yx)",.true.)
    call ctab_setScientific(rtable,"L2-error Im(SSE_yy)",.true.)

    call ctab_setScientific(rtable,"H1-error Re(SSE)",.true.)
    call ctab_setScientific(rtable,"H1-error Re(SSE_x)",.true.)
    call ctab_setScientific(rtable,"H1-error Re(SSE_y)",.true.)
    call ctab_setScientific(rtable,"H1-error Re(SSE_xx)",.true.)
    call ctab_setScientific(rtable,"H1-error Re(SSE_xy)",.true.)
    call ctab_setScientific(rtable,"H1-error Re(SSE_yx)",.true.)
    call ctab_setScientific(rtable,"H1-error Re(SSE_yy)",.true.)

    call ctab_setScientific(rtable,"H1-error Im(SSE)",.true.)
    call ctab_setScientific(rtable,"H1-error Im(SSE_x)",.true.)
    call ctab_setScientific(rtable,"H1-error Im(SSE_y)",.true.)
    call ctab_setScientific(rtable,"H1-error Im(SSE_xx)",.true.)
    call ctab_setScientific(rtable,"H1-error Im(SSE_xy)",.true.)
    call ctab_setScientific(rtable,"H1-error Im(SSE_yx)",.true.)
    call ctab_setScientific(rtable,"H1-error Im(SSE_yy)",.true.)

    ! Set Tex captions
    call ctab_setTexCaption(rtable,"cells","\# cells")
    call ctab_setTexCaption(rtable,"dofs","\# dofs")

    call ctab_setTexCaption(rtable,"L2-error Re(SSE)","$L^2(Re(N))$")
    call ctab_setTexCaption(rtable,"L2-error Im(SSE)","$L^2(Im(N))$")
    select case(rproblem%cproblemtype)
    case (SSE_SCALAR)
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_x)","$L^2(\partial_{x}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_y)","$L^2(\partial_{y}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_xx)","$L^2(\partial_{xx}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_xy)","$L^2(\partial_{xy}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_yx)","$L^2(\partial_{yx}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_yy)","$L^2(\partial_{yy}^{\rm ZZ}Re(N))$")

      call ctab_setTexCaption(rtable,"L2-error Im(SSE_x)","$L^2(\partial_{x}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_y)","$L^2(\partial_{y}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_xx)","$L^2(\partial_{xx}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_xy)","$L^2(\partial_{xy}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_yx)","$L^2(\partial_{yx}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_yy)","$L^2(\partial_{yy}^{\rm ZZ}Im(N))$")

      call ctab_setTexCaption(rtable,"H1-error Re(SSE_x)","$H^1(\partial_{x}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_y)","$H^1(\partial_{y}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_xx)","$H^1(\partial_{xx}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_xy)","$H^1(\partial_{xy}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_yx)","$H^1(\partial_{yx}^{\rm ZZ}Re(N))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_yy)","$H^1(\partial_{yy}^{\rm ZZ}Re(N))$")

      call ctab_setTexCaption(rtable,"H1-error Im(SSE_x)","$H^1(\partial_{x}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_y)","$H^1(\partial_{y}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_xx)","$H^1(\partial_{xx}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_xy)","$H^1(\partial_{xy}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_yx)","$H^1(\partial_{yx}^{\rm ZZ}Im(N))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_yy)","$H^1(\partial_{yy}^{\rm ZZ}Im(N))$")

    case (SSE_SYSTEM1,SSE_SYSTEM2)
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_x)","$L^2(Re(\sigma_x))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_y)","$L^2(Re(\sigma_y))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_xx)","$L^2(\partial_{x}^{\rm ZZ}Re(\sigma_x))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_xy)","$L^2(\partial_{y}^{\rm ZZ}Re(\sigma_x))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_yx)","$L^2(\partial_{x}^{\rm ZZ}Re(\sigma_y))$")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_yy)","$L^2(\partial_{y}^{\rm ZZ}Re(\sigma_y))$")

      call ctab_setTexCaption(rtable,"L2-error Im(SSE_x)","$L^2(Im(\sigma_x))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_y)","$L^2(Im(\sigma_y))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_xx)","$L^2(\partial_{x}^{\rm ZZ}Im(\sigma_x))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_xy)","$L^2(\partial_{y}^{\rm ZZ}Im(\sigma_x))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_yx)","$L^2(\partial_{x}^{\rm ZZ}Im(\sigma_y))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_yy)","$L^2(\partial_{y}^{\rm ZZ}Im(\sigma_y))$")

      call ctab_setTexCaption(rtable,"H1-error Re(SSE_x)","$H^1(Re(\sigma_x))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_y)","$H^1(Re(\sigma_y))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_xx)","$H^1(\partial_{x}^{\rm ZZ}Re(\sigma_x))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_xy)","$H^1(\partial_{y}^{\rm ZZ}Re(\sigma_x))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_yx)","$H^1(\partial_{x}^{\rm ZZ}Re(\sigma_y))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_yy)","$H^1(\partial_{y}^{\rm ZZ}Re(\sigma_y))$")

      call ctab_setTexCaption(rtable,"H1-error Im(SSE_x)","$H^1(Im(\sigma_x))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_y)","$H^1(Im(\sigma_y))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_xx)","$H^1(\partial_{x}^{\rm ZZ}Im(\sigma_x))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_xy)","$H^1(\partial_{y}^{\rm ZZ}Im(\sigma_x))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_yx)","$H^1(\partial_{x}^{\rm ZZ}Im(\sigma_y))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_yy)","$H^1(\partial_{y}^{\rm ZZ}Im(\sigma_y))$")
    end select

    call ctab_setTexCaption(rtable,"L2-error Re(SSE)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Re(SSE_x)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Re(SSE_y)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Re(SSE_xx)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Re(SSE_xy)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Re(SSE_yx)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Re(SSE_yy)-convrate","")

    call ctab_setTexCaption(rtable,"L2-error Im(SSE)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Im(SSE_x)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Im(SSE_y)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Im(SSE_xx)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Im(SSE_xy)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Im(SSE_yx)-convrate","")
    call ctab_setTexCaption(rtable,"L2-error Im(SSE_yy)-convrate","")

    call ctab_setTexCaption(rtable,"H1-error Re(SSE)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Re(SSE_x)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Re(SSE_y)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Re(SSE_xx)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Re(SSE_xy)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Re(SSE_yx)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Re(SSE_yy)-convrate","")

    call ctab_setTexCaption(rtable,"H1-error Im(SSE)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Im(SSE_x)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Im(SSE_y)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Im(SSE_xx)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Im(SSE_xy)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Im(SSE_yx)-convrate","")
    call ctab_setTexCaption(rtable,"H1-error Im(SSE_yy)-convrate","")

    ! Set Tex format
    call ctab_setTexFormat(rtable,"cells","r")
    call ctab_setTexFormat(rtable,"dofs","r")

    ! Hide all H1-columns
    call ctab_setHidden(rtable,"H1-error Re(SSE)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)-convrate",.true.)

    call ctab_setHidden(rtable,"H1-error Im(SSE)",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_x)",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_y)",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_x)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_y)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_xx)",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_xy)",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_yx)",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_yy)",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_xx)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_xy)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_yx)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Im(SSE_yy)-convrate",.true.)

    ! Hide imaginary L2-columns
    call ctab_setHidden(rtable,"L2-error Im(SSE)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_x)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_y)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_x)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_y)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xx)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xy)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yx)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yy)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xy)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yy)-convrate",.true.)

    call ctab_setTexTableCaption(rtable,"$L^2$-Convergence table, real part")
    call ctab_setTexTableLabel(rtable,"tab:l2_convergence_rate_real")

    ! Write convergence table to Tex file
    call ctab_outputTex(rtable,'./table_l2.tex')

    ! Unhide imaginary parts of L2-columns
    call ctab_setHidden(rtable,"L2-error Im(SSE)",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE)-convrate",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_x)",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_y)",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_x)-convrate",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_y)-convrate",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xx)",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xy)",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yx)",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yy)",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xx)-convrate",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xy)-convrate",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yx)-convrate",.false.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yy)-convrate",.false.)

    ! Hide real parts of L2-columns
    call ctab_setHidden(rtable,"L2-error Re(SSE)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_x)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_y)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_x)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_y)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xx)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xy)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yx)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yy)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xy)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yy)-convrate",.true.)

    call ctab_setTexTableCaption(rtable,"$L^2$-Convergence table, imaginary part")
    call ctab_setTexTableLabel(rtable,"tab:l2_convergence_rate_aimag")

    ! Write convergence table to Tex file
    call ctab_outputTex(rtable,'./table_l2.tex',SYS_APPEND)

    ! Hide all L2-columns
    call ctab_setHidden(rtable,"L2-error Re(SSE)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_x)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_y)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_x)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_y)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xx)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xy)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yx)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yy)",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_xy)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Re(SSE_yy)-convrate",.true.)

    call ctab_setHidden(rtable,"L2-error Im(SSE)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_x)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_y)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_x)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_y)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xx)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xy)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yx)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yy)",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_xy)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yx)-convrate",.true.)
    call ctab_setHidden(rtable,"L2-error Im(SSE_yy)-convrate",.true.)

    ! Unhide real parts of H1-columns
    call ctab_setHidden(rtable,"H1-error Re(SSE)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)-convrate",.false.)

    call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table, real part")
    call ctab_setTexTableLabel(rtable,"tab:h1_convergence_rate_real")

    ! Write convergence table to Tex file
    call ctab_outputTex(rtable,'./table_h1.tex')

    ! Hide real parts of H1-columns
    call ctab_setHidden(rtable,"H1-error Re(SSE)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)-convrate",.true.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)-convrate",.true.)

    ! Unhide imaginary parts of H1-columns
    call ctab_setHidden(rtable,"H1-error Re(SSE)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_x)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_y)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xx)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_xy)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yx)-convrate",.false.)
    call ctab_setHidden(rtable,"H1-error Re(SSE_yy)-convrate",.false.)

    call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table, imaginary part")
    call ctab_setTexTableLabel(rtable,"tab:h1_convergence_rate_aimag")

    ! Write convergence table to Tex file
    call ctab_outputTex(rtable,'./table_h1.tex', SYS_APPEND)

  case default
    call output_line("Invalid type of problem.", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_outputTable")
    call sys_halt()
  end select

  end subroutine

end module sse_main
