!##############################################################################
!# ****************************************************************************
!# <name> codire_init </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to initialize
!# the solver for a conservation law for a scalar variable
!#
!# The following routines are available:
!#
!# 1.) codire_initProblem
!#     -> initialize the global problem structure
!#
!# 2.) codire_initParser
!#     -> initialize the global function parser
!#
!# 3.) codire_initConstOperators1d
!#     -> initialize constant coefficient operators for the finite
!#        element discretization in 1D
!#
!# 4.) codire_initConstOperators2d
!#     -> initialize constant coefficient operators for the finite
!#        element discretization in 2D
!#
!# 5.) codire_initConstOperators3d
!#     -> initialize constant coefficient operators for the finite
!#        element discretization in 3D
!#
!# 6.) codire_initSolutionFromParser
!#     -> initialize the solution profile at time t=0
!#
!# 7.) codire_initSolutionFromImage
!#     -> initialize the solution profile at time t=0
!#
!# 8.) codire_initExactSolution
!#     -> initialize the exact solution profile
!#
!# 9.) codire_initRHS
!#      -> initialize the constant right-hand side
!#
!# 10.) codire_initVelocity
!#      -> initialize the velocity profile
!#
!# 11.) codire_initDiffusion1d
!#      -> initialize the diffusion "matrix" in 1D
!#
!# 12.) codire_initDiffusion2d
!#      -> initialize the diffusion matrix in 2D
!#
!# 13.) codire_initDiffusion3d
!#      -> initialize the diffusion matrix in 3D
!# 
!# </purpose>
!##############################################################################

module codire_init

  use afcstabilisation
  use fparser
  use fsystem
  use genoutput
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocsolution
  use spatialdiscretisation
  use stdoperators
  use storage
  use triangulation
  
  use codire_basic  
  use codire_callback
  use boundaryfilter
  use problem
  use solver

  implicit none

  private
  public :: codire_initProblem
  public :: codire_initParser
  public :: codire_initConstOperators1d
  public :: codire_initConstOperators2d
  public :: codire_initConstOperators3d
  public :: codire_initSolutionFromParser
  public :: codire_initSolutionFromImage
  public :: codire_initExactSolution
  public :: codire_initRHS
  public :: codire_initVelocity 
  public :: codire_initDiffusion1d
  public :: codire_initDiffusion2d
  public :: codire_initDiffusion3d

contains
 
   !*****************************************************************************

!<subroutine>

  subroutine codire_initProblem(rproblem, trifile, prmfile, indatfile,&
                             nlmin, nlmax, iconvToTria, ndimension)

!<description>
    ! This subroutine initializes the complete triangulation
    ! structures for a two-dimensional problem, that is,
    !
    ! (1) the boundary description is read from file and generated
    !
    ! (2) the coarse grid is read from file and data structures are generated
    !
    ! (3) all finer grids are assembled from regular refinement
    !
!</description>

!<input>
    ! name of file containing the parametrization
    character(LEN=*), intent(IN) :: prmfile

    ! name of file containing the triangulation
    character(LEN=*), intent(IN) :: trifile
    
    ! name of file containing the initial/boundary data
    character(LEN=*), intent(IN) :: indatfile
    
    ! minimum multigrid level
    integer, intent(IN) :: nlmin

    ! maximum multigrid level
    integer, intent(IN) :: nlmax

    ! convert grid to triangular grid?
    integer, intent(IN) :: iconvToTria

    ! number of spatial dimension
    integer, intent(IN) :: ndimension
!</input>

!<output>
    ! global problem structure
    type(t_problem), intent(OUT) :: rproblem
!</output>
!</subroutine>

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel, p_rproblemLevel
    integer :: ilev
    logical :: berror

    ! Initialize global problem structure
    call problem_createProblem(rproblem)
            
    ! Initialize coarse level
    nullify(rproblemLevel); allocate(rproblemLevel)
    call problem_createLevel(rproblemLevel, nlmin)
    
    select case(ndimension)
    case (NDIM1D)
      ! Read coarse mesh from TRI-file without generating an extended raw mesh
      call tria_readTriFile1D(rproblemLevel%rtriangulation,&
                              trifile, .true.)
      
    case (NDIM2D)
      ! Create new boundary and read from PRM-file
      call boundary_read_prm(rproblem%rboundary, prmfile)
      
      ! Read coarse mesh from TRI-file, convert it to triangular mesh if required
      call tria_readTriFile2D(rproblemLevel%rtriangulation,&
                              trifile, rproblem%rboundary, .true.)
      if (iconvToTria .eq. 1) call tria_rawGridToTri(rproblemLevel%rtriangulation)

    case (NDIM3D)
      call tria_readTriFile3D(rproblemLevel%rtriangulation,&
                              trifile, rproblem%rboundary, .true.)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'codire_initProblem')
      call sys_halt()
    end select
    
    ! Refine coarse mesh to minimum multigrid level and create standard mesh
    call tria_quickRefine2LevelOrdering(nlmin-1, rproblemLevel%rtriangulation,&
                                        rproblem%rboundary)
    call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
                                      rproblem%rboundary)

    ! Append level to global problem
    call problem_appendLevel(rproblem, rproblemLevel)
    p_rproblemLevel => rproblemLevel

    ! Generate fine levels
    do ilev = nlmin+1, nlmax
      
      ! Initialize current level
      nullify(rproblemLevel); allocate(rproblemLevel)
      call problem_createLevel(rproblemLevel, ilev)
      
      ! Generate regularly refined mesh
      call tria_refine2LevelOrdering(p_rproblemLevel%rtriangulation,&
                                     rproblemLevel%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Append current level to global problem
      call problem_appendLevel(rproblem, rproblemLevel)
      p_rproblemLevel => rproblemLevel
    end do
    
    ! Compress triangulation structure
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while (associated(p_rproblemLevel) .and.&
              associated(p_rproblemLevel%p_rproblemLevelCoarse))
      call tria_compress2LevelOrdHierarchy(p_rproblemLevel%rtriangulation,&
                                           p_rproblemLevel%p_rproblemLevelCoarse%rtriangulation)
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine codire_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine codire_initParser(sfilename)

!<description>
    ! This subroutine initializes the global function parser.
    ! The given file is parsed for constants and fixed expressions.
!</description>

!<input>
    ! name of parameter file
    character(LEN=*), intent(IN) :: sfilename
!</input>
!</subroutine>

    ! local variables
    character(SYS_STRLEN) :: keyword,name
    character(LEN=1024) :: sdata,expression
    real(DP) :: dvalue
    integer :: iunit,ipos,jpos,ios,idatalen

    ! Try to open the file
    call io_openFileForReading(sfilename, iunit, .true.)
    
    ! Oops...
    if (iunit .eq. -1) then
      call output_line('Unable to open input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'codire_initParser')
      call sys_halt()
    end if

    ! Read through the complete input file and look for global definitions
    ! of constants and fixed expressions
    ios = 0
    do while(ios .eq. 0)

      ! Read next line in file
      call io_readlinefromfile(iunit, sdata, idatalen, ios)
      
      ! Check for keyword defconst or defexp
      ipos = scan(sdata(1:idatalen), ":")
      if (ipos .eq. 0) cycle
      
      call sys_tolower(sdata(1:max(1,ipos-1)), keyword)
      select case(trim(adjustl(keyword)))
      case ("defconst")
        
        ! Split the line into name and value
        jpos    = scan(sdata(1:idatalen), "=" , .true.)
        name    = trim(adjustl(sdata(ipos+1:jpos-1)))
        keyword = trim(adjustl(sdata(jpos+1:)))
        read(keyword,*) dvalue
        
        ! We found a constant that will be applied to the parser
        call fparser_defineConstant(name, dvalue)
        
      case ("defexpr")
        
        ! Split the line into name and value
        jpos       = scan(sdata(1:idatalen), "=" ,.true.)
        name       = trim(adjustl(sdata(ipos+1:jpos-1)))
        expression = trim(adjustl(sdata(jpos+1:)))
        
        ! We found an expression that will be applied to the parser
        call fparser_defineExpression(name, expression)
        
      end select
    end do

    ! Close file
    close (iunit)
  end subroutine codire_initParser
  
  !*****************************************************************************

!<subroutine>

  subroutine codire_initConstOperators1d(rboundaryCondition, rproblem, ieltype,&
                                      imatrixFormat, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine initializes the discrete operators resulting
    ! from the finite element method in 2D. 
    ! The two optional parameters NLMIN and NLMAX can be used to 
    ! restrict the range of multigrid levels for which the initialization
    ! should be performed. If they are omitted, then all levels from 
    ! coarsest to finest grid are processed. 
!</description>

!<input>
    ! boundary conditions
    type(t_boundaryCondition), intent(IN) :: rboundaryCondition

    ! global problem structure
    type(t_problem), intent(IN) :: rproblem

    ! type of finite elements
    integer, intent(IN) :: ieltype

    ! type of matrix format
    ! LSYSSC_MATRIX7 or LSYSSC_MATRIX9
    integer, intent(IN) :: imatrixFormat

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>
!</subroutine>
    
    ! section name for primal diffusive stabilization
    character(LEN=SYS_STRLEN) :: sprimaldiffName

    ! section name for primal convective stabilization
    character(LEN=SYS_STRLEN) :: sprimalconvName

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel
    type(t_bilinearForm) :: rform
    integer :: nlmin,nlmax,i


    ! Set minimal/maximal levels: If given by the user, 
    ! adopt these values. Otherwise, assemble the constant 
    ! matrices for all levels of the multigrid structure.
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = rproblem%p_rproblemLevelMin%ilev
    end if

    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = rproblem%p_rproblemLevelMax%ilev
    end if
    

    ! Assemble the transport and global operator for all multigrid levels
    rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(rproblemLevel))
      
      ! Do we have to assemble operators for this level?
      if ((rproblemLevel%ilev < nlmin) .or.&
          (rproblemLevel%ilev > nlmax)) then
        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
        cycle
      end if

      ! Check if triangulation is 1D
      if (rproblemLevel%rtriangulation%ndim .ne. NDIM1D) then
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_WARNING, OU_MODE_STD, 'codire_initConstOperators1d')
        call sys_halt()
      end if

      ! Initialize the discretization structure
      call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, 1,&
                                   rproblemLevel%rtriangulation)
      
      select case(ieltype)
      case (-1,1,11)
        ! P1 finite elements
        call spdiscr_initDiscr_simple (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)

      case DEFAULT
        call output_line('Unsupproted element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators1d')
        call sys_halt()
      end select
      

      ! Allocate memory for 8 FE-matrix structures, that is
      !   the template matrix (1x),
      !   the system matrix A (1x),
      !   the consistent/lumped mass matrices (2x),
      !   the spatial derivatives (2x),
      !   the discrete transport operator (1x), and
      !   the Jacobian matrix (1x)
      if (associated(rproblemLevel%Rmatrix)) then
        do i = lbound(rproblemLevel%Rmatrix,1),&
               ubound(rproblemLevel%Rmatrix,1)
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(i))
        end do
        deallocate(rproblemLevel%Rmatrix)
      end if
      allocate(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX))

      
      ! Allocate memory for 2 stabilisation structures
      if (associated(rproblemLevel%Rafcstab)) then
        do i = lbound(rproblemLevel%Rafcstab,1),&
               ubound(rproblemLevel%Rafcstab,1)
          call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
        end do
        deallocate(rproblemLevel%Rafcstab)
      end if
      allocate(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION))


      ! Generate the finite element matrix sparsity structure
      call bilf_createMatrixStructure(&
          rproblemLevel%rdiscretisation%RspatialDiscr(1),&
          imatrixFormat, rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE))

      ! Modify structure to incorporate periodic boundary conditions
      call bdrf_calcMatrixPeriodic(rboundaryCondition, rproblemLevel%rtriangulation,&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE))

      ! Generate structure for scalar matrix A and low-order operator L
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      
      ! Generate the global coefficient matrices for the FE method
      select case(abs(ivelocitytype))
      case (VELOCITY_NONE)
        ! zero velocity, do nothing

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP,&
            VELOCITY_BURGERS1D,&
            VELOCITY_BUCKLEV1D)
        ! non-zero velocity, assemble coefficient matrices CX
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                         DER_DERIV1D_X, DER_FUNC1D)

        ! Initialize AFC stabilisation for convective part
        call parlst_getvalue_string(rparlist, '', "primalconv",  sprimalconvName)
        call afcstab_initFromParameterlist(rparlist, sprimalconvName,&
                                           rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))

      case DEFAULT
        call output_line('Invalid type of velocity!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators1d')
      end select


      ! Generate the diffusion matrix
      select case(idiffusiontype)
      case (DIFF_NONE)
        ! zero diffusion, do nothing


      case (DIFF_ISOTROPIC,&
            DIFF_ANISOTROPIC)
        ! Isotropic diffusion, so generate the standard Laplace matrix and
        ! scale it by the negative value of the diffusion coefficient.
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleLaplaceMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .true.,&
                                         -DdiffusionMatrix1D(1,1))


      case DEFAULT
        call output_line('Invalid type of diffusion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators1d')
        call sys_halt()
      end select


      ! Assemble consistent mass matrix MC
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                      DER_FUNC1D, DER_FUNC1D)
      
      ! Perform row-sum mass lumping ML:=lumped(MC)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
      call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                   LSYSSC_LUMP_DIAG)
      

      ! Switch to next coarser level
      rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine codire_initConstOperators1d

  !*****************************************************************************

!<subroutine>

  subroutine codire_initConstOperators2d(rboundaryCondition, rproblem, ieltype,&
                                      imatrixFormat, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine initializes the discrete operators resulting
    ! from the finite element method in 2D. 
    ! The two optional parameters NLMIN and NLMAX can be used to 
    ! restrict the range of multigrid levels for which the initialization
    ! should be performed. If they are omitted, then all levels from 
    ! coarsest to finest grid are processed. 
!</description>

!<input>
    ! boundary conditions
    type(t_boundaryCondition), intent(IN) :: rboundaryCondition

    ! global problem structure
    type(t_problem), intent(IN) :: rproblem

    ! type of finite elements
    integer, intent(IN) :: ieltype

    ! type of matrix format
    ! LSYSSC_MATRIX7 or LSYSSC_MATRIX9
    integer, intent(IN) :: imatrixFormat

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>
!</subroutine>
    
    ! section name for primal diffusive stabilization
    character(LEN=SYS_STRLEN) :: sprimaldiffName

    ! section name for primal convective stabilization
    character(LEN=SYS_STRLEN) :: sprimalconvName

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel
    type(t_bilinearForm) :: rform
    integer :: nlmin,nlmax,i

    ! Set minimal/maximal levels: If given by the user, 
    ! adopt these values. Otherwise, assemble the constant 
    ! matrices for all levels of the multigrid structure.
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = rproblem%p_rproblemLevelMin%ilev
    end if

    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = rproblem%p_rproblemLevelMax%ilev
    end if
    

    ! Assemble the transport and global operator for all multigrid levels
    rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(rproblemLevel))
      
      ! Do we have to assemble operators for this level?
      if ((rproblemLevel%ilev < nlmin) .or.&
          (rproblemLevel%ilev > nlmax)) then
        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
        cycle
      end if

      ! Check if triangulation is 2D
      if (rproblemLevel%rtriangulation%ndim .ne. NDIM2D) then
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_WARNING, OU_MODE_STD, 'codire_initConstOperators2d')
        call sys_halt()
      end if

      ! Initialize the discretization structure
      call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, 1,&
                                   rproblemLevel%rtriangulation)
      
      select case(ieltype)
      case (1)
        ! P1 finite elements
        call spdiscr_initDiscr_simple (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E001, SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
        
      case (11)
        ! Q1 finite elements
        call spdiscr_initDiscr_simple (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E011, SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
        
      case (-1)
        ! mixed P1/Q1 finite elements
        call spdiscr_initDiscr_triquad (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC,&
                                        SPDISC_CUB_AUTOMATIC,&
                                        rproblemLevel%rtriangulation,&
                                        rproblemLevel%p_rproblem%rboundary)

      case DEFAULT
        call output_line('Unsupproted element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators2d')
        call sys_halt()
      end select
      

      ! Allocate memory for 9 FE-matrix structures, that is
      !   the template matrix (1x),
      !   the system matrix A (1x),
      !   the consistent/lumped mass matrices (2x),
      !   the spatial derivatives (3x),
      !   the discrete transport operator (1x), and
      !   the Jacobian matrix (1x)
      if (associated(rproblemLevel%Rmatrix)) then
        do i = lbound(rproblemLevel%Rmatrix,1),&
               ubound(rproblemLevel%Rmatrix,1)
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(i))
        end do
        deallocate(rproblemLevel%Rmatrix)
      end if
      allocate(rproblemLevel%Rmatrix(CDEQ_MATRIX_CY))

      
      ! Allocate memory for 2 stabilisation structures
      if (associated(rproblemLevel%Rafcstab)) then
        do i = lbound(rproblemLevel%Rafcstab,1),&
               ubound(rproblemLevel%Rafcstab,1)
          call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
        end do
        deallocate(rproblemLevel%Rafcstab)
      end if
      allocate(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION))


      ! Generate the finite element matrix sparsity structure
      call bilf_createMatrixStructure(&
          rproblemLevel%rdiscretisation%RspatialDiscr(1),&
          imatrixFormat, rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE))

      ! Modify structure to incorporate periodic boundary conditions
      call bdrf_calcMatrixPeriodic(rboundaryCondition, rproblemLevel%rtriangulation,&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE))

      ! Generate structure for scalar matrix A and low-order operator L
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      
      ! Generate the global coefficient matrices for the FE method
      select case(abs(ivelocitytype))
      case (VELOCITY_NONE)
        ! zero velocity, do nothing

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP,&
            VELOCITY_BURGERS_SPACETIME,&
            VELOCITY_BUCKLEV_SPACETIME,&
            VELOCITY_BURGERS2D)
        ! non-zero velocity, assemble coefficient matrix CX
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                         DER_DERIV2D_X, DER_FUNC2D)

        ! non-zero velocity, assemble coefficient matrix CY
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                                        DER_DERIV2D_Y, DER_FUNC2D)

        ! Initialize AFC stabilisation for convective part
        call parlst_getvalue_string(rparlist, '', "primalconv",  sprimalconvName)
        call afcstab_initFromParameterlist(rparlist, sprimalconvName,&
                                           rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))

      case DEFAULT
        call output_line('Invalid type of velocity!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators2d')
      end select


      ! Generate the diffusion matrix
      select case(idiffusiontype)
      case (DIFF_NONE)
        ! zero diffusion, do nothing


      case (DIFF_ISOTROPIC)
        ! Isotropic diffusion, so generate the standard Laplace matrix and 
        ! scale it by the negative value of the diffusion coefficient.
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleLaplaceMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .true.,&
                                         -DdiffusionMatrix2D(1,1))

        
      case (DIFF_ANISOTROPIC)
        ! For anisotropic diffusion, things are slightly more complicated.
        ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
        ! scalar system matrix in 2D.
        rform%itermCount = 4
        rform%Idescriptors(1,1) = DER_DERIV2D_X
        rform%Idescriptors(2,1) = DER_DERIV2D_X
        
        rform%Idescriptors(1,2) = DER_DERIV2D_X
        rform%Idescriptors(2,2) = DER_DERIV2D_Y
        
        rform%Idescriptors(1,3) = DER_DERIV2D_Y
        rform%Idescriptors(2,3) = DER_DERIV2D_X
        
        rform%Idescriptors(1,4) = DER_DERIV2D_Y
        rform%Idescriptors(2,4) = DER_DERIV2D_Y
        
        ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = -DdiffusionMatrix2D(1,1)
        rform%Dcoefficients(2)  = -DdiffusionMatrix2D(1,2)
        rform%Dcoefficients(3)  = -DdiffusionMatrix2D(2,1)
        rform%Dcoefficients(4)  = -DdiffusionMatrix2D(2,2)

        ! Now we can build the matrix entries.
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call bilf_buildMatrixScalar(rform, .true.,&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_S))
        
        ! Initialize AFC stabilisation for diffusion
        call parlst_getvalue_string(rparlist, '', "primaldiff",  sprimaldiffName)
        call afcstab_initFromParameterlist(rparlist, sprimaldiffName,&
                                           rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION))

      case DEFAULT
        call output_line('Invalid type of diffusion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators2d')
        call sys_halt()
      end select


      ! Assemble consistent mass matrix MC
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                      DER_FUNC2D, DER_FUNC2D) 
      
      ! Perform row-sum mass lumping ML:=lumped(MC)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
      call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                   LSYSSC_LUMP_DIAG)
      

      ! Switch to next coarser level
      rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine codire_initConstOperators2d

  !*****************************************************************************

!<subroutine>

  subroutine codire_initConstOperators3d(rboundaryCondition, rproblem, ieltype,&
                                      imatrixFormat, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine initializes the discrete operators resulting
    ! from the finite element method in 3D. 
    ! The two optional parameters NLMIN and NLMAX can be used to 
    ! restrict the range of multigrid levels for which the initialization
    ! should be performed. If they are omitted, then all levels from 
    ! coarsest to finest grid are processed. 
!</description>

!<input>
    ! boundary conditions
    type(t_boundaryCondition), intent(IN) :: rboundaryCondition

    ! global problem structure
    type(t_problem), intent(IN) :: rproblem

    ! type of finite elements
    integer, intent(IN) :: ieltype

    ! type of matrix format
    ! LSYSSC_MATRIX7 or LSYSSC_MATRIX9
    integer, intent(IN) :: imatrixFormat

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>
!</subroutine>

    ! section name for primal diffusive stabilization
    character(LEN=SYS_STRLEN) :: sprimaldiffName

    ! section name for primal convective stabilization
    character(LEN=SYS_STRLEN) :: sprimalconvName

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel
    type(t_bilinearForm) :: rform
    integer :: nlmin,nlmax,i

    ! Set minimal/maximal levels: If given by the user, 
    ! adopt these values. Otherwise, assemble the constant 
    ! matrices for all levels of the multigrid structure.
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = rproblem%p_rproblemLevelMin%ilev
    end if

    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = rproblem%p_rproblemLevelMax%ilev
    end if
    

    ! Assemble the transport and global operator for all multigrid levels
    rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(rproblemLevel))
      
      ! Do we have to assemble operators for this level?
      if ((rproblemLevel%ilev < nlmin) .or.& 
          (rproblemLevel%ilev > nlmax)) then
        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
        cycle
      end if

      ! Check if triangulation is 3D
      if (rproblemLevel%rtriangulation%ndim .ne. NDIM3D) then
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_WARNING, OU_MODE_STD, 'codire_initConstOperators3d')
        call sys_halt()
      end if

      ! Initialize the discretization structure
      call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, 1,&
                                   rproblemLevel%rtriangulation)
      
      select case(ieltype)
      case (1)
        ! P1 finite elements
        call spdiscr_initDiscr_simple (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E001_3D, SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
        
      case (11)
        ! Q1 finite elements
        call spdiscr_initDiscr_simple (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E010_3D, SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)
        
      case (-1)
        ! mixed P1/Q1 finite elements
        call spdiscr_initDiscr_triquad (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                        EL_E001_3D, EL_E010_3D, SPDISC_CUB_AUTOMATIC,&
                                        SPDISC_CUB_AUTOMATIC,&
                                        rproblemLevel%rtriangulation,&
                                        rproblemLevel%p_rproblem%rboundary)

      case DEFAULT
        call output_line('Unsupproted element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators3d')
        call sys_halt()
      end select
      

      ! Allocate memory for 10 FE-matrix structures, that is
      !   the template matrix (1x),
      !   the system matrix A (1x),
      !   the consistent/lumped mass matrices (2x),
      !   the spatial derivatives (4x),
      !   the discrete transport operator (1x), and
      !   the Jacobian matrix (1x)
      if (associated(rproblemLevel%Rmatrix)) then
        do i = lbound(rproblemLevel%Rmatrix,1),&
               ubound(rproblemLevel%Rmatrix,1)
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(i))
        end do
        deallocate(rproblemLevel%Rmatrix)
      end if
      allocate(rproblemLevel%Rmatrix(CDEQ_MATRIX_CZ))

      
      ! Allocate memory for 2 stabilisation structures
      if (associated(rproblemLevel%Rafcstab)) then
        do i = lbound(rproblemLevel%Rafcstab,1),&
               ubound(rproblemLevel%Rafcstab,1)
          call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
        end do
        deallocate(rproblemLevel%Rafcstab)
      end if
      allocate(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION))


      ! Generate the finite element matrix sparsity structure
      call bilf_createMatrixStructure(&
          rproblemLevel%rdiscretisation%RspatialDiscr(1),&
          imatrixFormat, rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE))

      ! Modify structure to incorporate periodic boundary conditions
      call bdrf_calcMatrixPeriodic(rboundaryCondition, rproblemLevel%rtriangulation,&
                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE))

      ! Generate structure for scalar matrix A and low-order operator L
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      
      ! Generate the global coefficient matrices for the FE method
      select case(abs(ivelocitytype))
      case (VELOCITY_NONE)
        ! zero velocity, do nothing

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP,&
            VELOCITY_BURGERS_SPACETIME,&
            VELOCITY_BUCKLEV_SPACETIME,&
            VELOCITY_BURGERS2D)
        ! non-zero velocity, assemble coefficient matrices CX
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX),&
                                         DER_DERIV3D_X, DER_FUNC3D)

        ! non-zero velocity, assemble coefficient matrices CY
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CY),&
                                        DER_DERIV3D_Y, DER_FUNC3D)

        ! non-zero velocity, assemble coefficient matrices CZ
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_CZ),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleSimpleMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_CZ),&
                                        DER_DERIV3D_Z, DER_FUNC3D)

        ! Initialize AFC stabilisation for convective part
        call parlst_getvalue_string(rparlist, '', "primalconv",  sprimalconvName)
        call afcstab_initFromParameterlist(rparlist, sprimalconvName,&
                                           rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))

      case DEFAULT
        call output_line('Invalid type of velocity!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators3d')
      end select


      ! Generate the diffusion matrix
      select case(idiffusiontype)
      case (DIFF_NONE)
        ! zero diffusion, do nothing

        
      case (DIFF_ISOTROPIC)
        ! Isotropic diffusion, so generate the standard Laplace matrix and
        ! scale it by  the negative value of the diffusion coefficient.
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call stdop_assembleLaplaceMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .true.,&
                                         -DdiffusionMatrix3D(1,1))

        
      case (DIFF_ANISOTROPIC)
        ! For anisotropic diffusion, things are slightly more complicated.
        ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
        ! scalar system matrix in 3D.
        rform%itermCount = 9
        rform%Idescriptors(1,1) = DER_DERIV3D_X
        rform%Idescriptors(2,1) = DER_DERIV3D_X
        
        rform%Idescriptors(1,2) = DER_DERIV3D_X
        rform%Idescriptors(2,2) = DER_DERIV3D_Y

        rform%Idescriptors(1,3) = DER_DERIV3D_X
        rform%Idescriptors(2,3) = DER_DERIV3D_Z
        
        rform%Idescriptors(1,4) = DER_DERIV3D_Y
        rform%Idescriptors(2,4) = DER_DERIV3D_X
        
        rform%Idescriptors(1,5) = DER_DERIV3D_Y
        rform%Idescriptors(2,5) = DER_DERIV3D_Y

        rform%Idescriptors(1,6) = DER_DERIV3D_Y
        rform%Idescriptors(2,6) = DER_DERIV3D_Z

        rform%Idescriptors(1,7) = DER_DERIV3D_Z
        rform%Idescriptors(2,7) = DER_DERIV3D_X

        rform%Idescriptors(1,8) = DER_DERIV3D_Z
        rform%Idescriptors(2,8) = DER_DERIV3D_Y

        rform%Idescriptors(1,9) = DER_DERIV3D_Z
        rform%Idescriptors(2,9) = DER_DERIV3D_Z

        
        ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = -DdiffusionMatrix3D(1,1)
        rform%Dcoefficients(2)  = -DdiffusionMatrix3D(1,2)
        rform%Dcoefficients(3)  = -DdiffusionMatrix3D(1,3)
        rform%Dcoefficients(4)  = -DdiffusionMatrix3D(2,1)
        rform%Dcoefficients(5)  = -DdiffusionMatrix3D(2,2)
        rform%Dcoefficients(6)  = -DdiffusionMatrix3D(2,3)
        rform%Dcoefficients(7)  = -DdiffusionMatrix3D(3,1)
        rform%Dcoefficients(8)  = -DdiffusionMatrix3D(3,2)
        rform%Dcoefficients(9)  = -DdiffusionMatrix3D(3,3)

        ! Now we can build the matrix entries.
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        call bilf_buildMatrixScalar(rform, .true.,&
                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_S))
        
        ! Initialize AFC stabilisation for diffusion
        call parlst_getvalue_string(rparlist, '', "primaldiff",  sprimaldiffName)
        call afcstab_initFromParameterlist(rparlist, sprimaldiffName,&
                                           rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION))

      case DEFAULT
        call output_line('Invalid type of diffusion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initConstOperators3d')
        call sys_halt()
      end select


      ! Assemble consistent mass matrix MC
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                       DER_FUNC3D, DER_FUNC3D) 
      
      ! Perform row-sum mass lumping ML:=lumped(MC)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
                                  rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
      call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
                                   LSYSSC_LUMP_DIAG)
      

      ! Switch to next coarser level
      rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine codire_initConstOperators3d

  !*****************************************************************************

!<subroutine>

  subroutine codire_initSolutionFromParser(rboundaryCondition, rproblemLevel,&
                                        rsolution, indatfile, time)

!<description>
    ! This subroutine initializes the solution vector from
    ! the parameters specified in the indatfile by calling
    ! a function parser w.r.t. the specified variable names.
!</description>
    
!<input>
    ! boundary conditions
    type(t_boundaryCondition), intent(IN) :: rboundaryCondition

    ! multigrid level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! name of indat file
    character(LEN=*), intent(IN) :: indatfile

    ! OPTIONAL: initial time
    real(DP), intent(IN), optional :: time
!</input>

!<output>
    ! solution vector
    type(t_vectorBlock), intent(OUT) :: rsolution
!</output>
!</subroutine>
    
    ! local constants
    character(LEN=*), dimension(2), parameter ::&
                      SOL_SYMBOLICVARS1D = (/ (/'x'/), (/'t'/) /)
    character(LEN=*), dimension(3), parameter ::&
                      SOL_SYMBOLICVARS2D = (/ (/'x'/), (/'y'/), (/'t'/) /)
    character(LEN=*), dimension(4), parameter ::&
                      SOL_SYMBOLICVARS3D = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)

    ! local variables
    type(t_problem), pointer :: p_rproblem
    type(t_fparser) :: rparser
    integer :: istatus
    
    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem

    ! Allocate storage for solution vector
    call lsysbl_createVectorBlock(rproblemLevel%rdiscretisation, &
                                  rsolution, .true., ST_DOUBLE)

    ! How many spatial dimensions do we have?
    select case(rproblemLevel%rdiscretisation%ndimension)
    case (NDIM1D)
      ! Create solution profile in 1D
      call problem_createProfile(p_rproblem, indatfile, '[initial_solution]',&
                                 SOL_SYMBOLICVARS1D, rparser, istatus, time, rsolution)
      
    case (NDIM2D)
      ! Create solution profile in 2D
      call problem_createProfile(p_rproblem, indatfile, '[initial_solution]',&
                                 SOL_SYMBOLICVARS2D, rparser, istatus, time, rsolution)

    case (NDIM3D)
      ! Create solution profile in 3D
      call problem_createProfile(p_rproblem, indatfile, '[initial_solution]',&
                                 SOL_SYMBOLICVARS3D, rparser, istatus, time, rsolution)
    
    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initSolutionFromParser')
      call sys_halt()
    end select

    ! Check status
    if (istatus .ne. 0) then
      call output_line('Unable to create initial solution profile!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initSolutionFromParser')
      call sys_halt()
    end if
    
    ! Release parser
    call fparser_release(rparser)

    ! Impose boundary conditions explicitely
    call bdrf_filterVectorExplicit(rboundaryCondition,&
                                   rproblemLevel%rtriangulation, rsolution, 0.0_DP)
  end subroutine codire_initSolutionFromParser
  
  !*****************************************************************************

!<subroutine>

  subroutine codire_initSolutionFromImage(rboundaryCondition, rproblemLevel,&
                                       rsolution, pgmfile)

!<description>
    ! This subroutine initializes the solution vector from the
    ! portable graymap image stored in file pgmfile.
!</description>
    
!<input>
    ! boundary conditions
    type(t_boundaryCondition), intent(IN) :: rboundaryCondition

    ! multigrid level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! name of portable graymap file
    character(LEN=*), intent(IN) :: pgmfile
!</input>

!<output>
    ! solution vector
    type(t_vectorBlock), intent(OUT) :: rsolution
!</output>
!</subroutine>

    ! local variables
    type(t_problem), pointer :: p_rproblem
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:),   pointer :: p_Ddata
    type(t_pgm) :: rpgm

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem

    ! Check if problem is in 2D
    if (rproblemLevel%rdiscretisation%ndimension .ne. NDIM2D) then
      call output_line('Initialization from image is only available in 2D!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initSolutionFromImage')
      call sys_halt()
    end if

    ! Allocate storage for solution vector
    call lsysbl_createVectorBlock(rproblemLevel%rdiscretisation, &
                                  rsolution, .true., ST_DOUBLE)

    ! Read PGM image from file
    call ppsol_readPGM(0, pgmfile, rpgm)
    
    call storage_getbase_double2D(rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_double(rsolution%h_Ddata, p_Ddata)
    call ppsol_initArrayPGM_Dble(rpgm, p_DvertexCoords, p_Ddata)
    
    ! Release PGM image
    call ppsol_releasePGM(rpgm)

    ! Impose boundary conditions explicitely
    call bdrf_filterVectorExplicit(rboundaryCondition,&
                                   rproblemLevel%rtriangulation, rsolution, 0.0_DP)
  end subroutine codire_initSolutionFromImage

  !*****************************************************************************

!<subroutine>

  subroutine codire_initExactSolution(rproblemLevel, rsolution,&
                                   indatfile, time, istatus)

!<description>
    ! This subroutine initializes the exact solution vector from the
    ! parameters specified in the indatfile by calling a function
    ! parser w.r.t. the specified variable names
!</description>
    
!<input>
    ! multigrid level structire
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! name of indat file
    character(LEN=*), intent(IN) :: indatfile

    ! simulation time
    real(DP), intent(IN) :: time
!</input>

!<output>
    ! solution vector
    type(t_vectorBlock), intent(OUT) :: rsolution

    ! OPTIONAL: status code
    integer, intent(OUT), optional :: istatus
!</output>
!</subroutine>

    ! local constants
    character(LEN=*), dimension(2), parameter ::&
                      SOL_SYMBOLICVARS1D = (/ (/'x'/), (/'t'/) /)
    character(LEN=*), dimension(3), parameter ::&
                      SOL_SYMBOLICVARS2D = (/ (/'x'/), (/'y'/), (/'t'/) /)
    character(LEN=*), dimension(4), parameter ::&
                      SOL_SYMBOLICVARS3D = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)

    ! local variables
    type(t_problem), pointer :: p_rproblem
    type(t_fparser) :: rparser
    integer :: iistatus

    if (present(istatus)) istatus = 0

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem
    
    ! Allocate storage for solution vector
    call lsysbl_createVectorBlock(rproblemLevel%rdiscretisation, &
                                  rsolution, .true., ST_DOUBLE)

    ! How many spatial dimensions do we have?
    select case(rproblemLevel%rdiscretisation%ndimension)
      
    case (NDIM1D)
      ! Create solution profile in 2D
      call problem_createProfile(p_rproblem, indatfile, '[exact_solution]',&
                                 SOL_SYMBOLICVARS1D, rparser, iistatus, time, rsolution)


    case (NDIM2D)
      ! Create solution profile in 2D
      call problem_createProfile(p_rproblem, indatfile, '[exact_solution]',&
                                 SOL_SYMBOLICVARS2D, rparser, iistatus, time, rsolution)

      
    case (NDIM3D)
      ! Create solution profile in 3D
      call problem_createProfile(p_rproblem, indatfile, '[exact_solution]',&
                                 SOL_SYMBOLICVARS3D, rparser, iistatus, time, rsolution)

      
    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initExactSolution')
      call sys_halt()
    end select
    
    ! No exact solution available; free memory and return
    if (iistatus .ne. 0) then
      call output_line('Unable to create exact solution profile!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initExactSolution')
      call lsysbl_releaseVector(rsolution)
      if (present(istatus)) istatus = -1
      return
    end if
    
    ! Release parser
    call fparser_release(rparser)
  end subroutine codire_initExactSolution

  !*****************************************************************************

!<subroutine>

  subroutine codire_initRHS(rboundaryCondition, rproblemLevel, rrhs, indatfile)

!<description>
    ! This subroutine initializes the right-hand side vector
    ! from the parameters specified in the indatfile by calling
    ! a function parser w.r.t. the specified variable names
!</description>
    
!<input>
    ! boundary conditions
    type(t_boundaryCondition), intent(IN) :: rboundaryCondition

    ! multigrid level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! name of indat file
    character(LEN=*), intent(IN) :: indatfile
!</input>

!<output>
    ! right-hand side vector
    type(t_vectorBlock), intent(OUT) :: rrhs
!</output>
!</subroutine>

    ! local constants
    character(LEN=*), dimension(1), parameter ::&
                      SOL_SYMBOLICVARS1D = (/ (/'x'/) /)
    character(LEN=*), dimension(2), parameter ::&
                      SOL_SYMBOLICVARS2D = (/ (/'x'/), (/'y'/) /)
    character(LEN=*), dimension(3), parameter ::&
                      SOL_SYMBOLICVARS3D = (/ (/'x'/), (/'y'/), (/'z'/) /)

    ! local variables
    type(t_problem), pointer :: p_rproblem
    type(t_collection) :: rcollection
    type(t_linearForm) :: rform
    type(t_vectorScalar) :: rvector
    integer :: istatus

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem

    ! What type of r.h.s. should be applied?
    select case(irhstype)
      
    case (RHS_ZERO)
      ! Ok, do nothing at all

      
    case (RHS_ANALYTIC)
      ! How many spatial dimensions do we have?
      select case(rproblemLevel%rdiscretisation%ndimension)

      case (NDIM1D)
        ! Initialize function parser in 2D
        call problem_createProfile(p_rproblem, indatfile, '[rhs]',&
                                   SOL_SYMBOLICVARS1D, rrhsParser, istatus)
        

      case (NDIM2D)
        ! Initialize function parser in 2D
        call problem_createProfile(p_rproblem, indatfile, '[rhs]',&
                                   SOL_SYMBOLICVARS2D, rrhsParser, istatus)

        
      case (NDIM3D)
        ! Initialize function parser in 3D
        call problem_createProfile(p_rproblem, indatfile, '[rhs]',&
                                   SOL_SYMBOLICVARS2D, rrhsParser, istatus)

        
      case DEFAULT
        call output_line('Invalid number of spatial dimensions',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initRHS')
        call sys_halt()
      end select

      
      if (istatus .ne. 0) then
        call output_line('Unable to create function parser for right-hand side!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_initRHS')
        call sys_halt()
      end if
      
      ! Create a collection; used for passing parameters to the RHS.
      call collct_init (rcollection)
      
      ! Set up the corresponding linear form (f,Phi_j):
      rform%itermCount      = 1
      rform%Idescriptors(1) = DER_FUNC

      ! Build the right-hand side vector
      call linf_buildVectorScalar (rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                   rform, .true., rvector, fcb_coeffRHS, rcollection)

      ! Convert the temporal scalar vector to a 1-block vector
      call lsysbl_convertVecFromScalar(rvector, rrhs, rproblemLevel%rdiscretisation)

      ! Release temporal scalar vector
      call lsyssc_releaseVector(rvector)

      ! Release collection
      call collct_done(rcollection)

      ! Impose boundary conditions
      call bdrf_filterVectorByValue(rboundaryCondition,&
                                    rproblemLevel%rtriangulation, rrhs, 0.0_DP)
      
    case DEFAULT
      call output_line('Invalid type of right-hand side!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initRHS')
      call sys_halt()
    end select
  end subroutine codire_initRHS

  !*****************************************************************************

!<subroutine>

  subroutine codire_initVelocity(rproblem, indatfile, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine initializes the velocity vector from the parameters
    ! specified in the indatfile by calling a function parser.
    ! If the optional parameters NLMINOPT and NLMAXOPT are given, then
    ! the velocity vector is only created for these multigrid levels.
    ! Note that the velocity is only evaluated form the function
    ! parser on the finest level. On all coarser levels it is restricted
    ! from the next finer one.
!</description>

!<input>
    ! name of indat file
    character(LEN=*), intent(IN) :: indatfile

    ! OPTIONAL: minimum multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximum multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>

!<inputoutput>
    ! global problem structure
    type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>
!</subroutine>
    
    ! local constants
    character(LEN=*), dimension(2), parameter ::&
                      VEL_SYMBOLICVARS1D = (/ (/'x'/), (/'t'/) /)
    character(LEN=*), dimension(3), parameter ::&
                      VEL_SYMBOLICVARS2D = (/ (/'x'/), (/'y'/), (/'t'/) /)
    character(LEN=*), dimension(4), parameter ::&
                      VEL_SYMBOLICVARS3D = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel
    integer(I32) :: isize
    integer :: NEQ, istatus, nlmin, nlmax, ndimension


    ! Set minimal/maximal levels: If given by the user, adopt these values.
    ! Otherwise, initialize velocity vector on all levels of the multigrid structure.
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = rproblem%p_rproblemLevelMin%ilev
    end if
    
    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = rproblem%p_rproblemLevelMax%ilev
    end if

    ! Process all multigrid levels
    rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(rproblemLevel))

      ! Do we have to initialize the velocity vector for this level?
      if (rproblemLevel%ilev < nlmin .or. rproblemLevel%ilev > nlmax) then
        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
        cycle
      end if
      
      ! Get total number of vertices
      isize = rproblemLevel%rtriangulation%NVT

      ! Get number of spatial dimensions
      ndimension = rproblem%p_rproblemLevelMax%rdiscretisation%ndimension
      
      ! Check if velocity vector exists 
      if (associated(rproblemLevel%rvectorBlock)) then
        
        ! Check if velocity vector needs to be resized
        NEQ = rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%NEQ
        if (NEQ .eq. 0) then
          ! Create two-dimensional velocity vector
          call lsysbl_createVectorBlock(rproblemLevel%rvectorBlock(CDEQ_VELOCITY),&
                                        isize, ndimension, .true.)
        elseif(NEQ .ne. isize*ndimension) then
          ! Resize two-dimensional velocity vector
          call lsysbl_resizeVectorBlock(rproblemLevel%rvectorBlock(CDEQ_VELOCITY),&
                                        isize, .true.)
        end if
        
      else
        
        ! Allocate block vector for velocity
        allocate(rproblemLevel%rvectorBlock(CDEQ_VELOCITY))
        
        ! Create two-dimensional velocity vector
        call lsysbl_createVectorBlock(rproblemLevel%rvectorBlock(CDEQ_VELOCITY),&
                                      isize, ndimension, .true.)    
      end if
      
      ! Are we on the finest level? 
      if (rproblemLevel%ilev .eq. nlmax) then
        
        ! We are on the finest multigrid level.
        ! Hence, create the velocity profile from INDAT-file

        ! How many spatial dimensions do we have?
        select case(ndimension)

        case (NDIM1D)
          ! Create velocity profile in 1D
          call problem_createProfile(rproblem, indatfile, '[initial_velocity]',&
                                     VEL_SYMBOLICVARS1D, rvelocityParser, istatus,&
                                     rvector=rproblemLevel%rvectorBlock(CDEQ_VELOCITY))


        case (NDIM2D)
          ! Create velocity profile in 2D
          call problem_createProfile(rproblem, indatfile, '[initial_velocity]',&
                                     VEL_SYMBOLICVARS2D, rvelocityParser, istatus,&
                                     rvector=rproblemLevel%rvectorBlock(CDEQ_VELOCITY))


        case (NDIM3D)
          ! Create velocity profile in 3D
          call problem_createProfile(rproblem, indatfile, '[initial_velocity]',&
                                     VEL_SYMBOLICVARS3D, rvelocityParser, istatus,&
                                     rvector=rproblemLevel%rvectorBlock(CDEQ_VELOCITY))


        case DEFAULT
          call output_line('Invalid number of spatial dimensions',&
                           OU_CLASS_ERROR,OU_MODE_STD,'codire_initVelocity')
          call sys_halt()
        end select
        
        if (istatus .ne. 0) then
          call output_line('Unable to create initial velocity!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'codire_initVelocity')
          call sys_halt()
        end if
        
      else
        
        ! We are on a coarse multigrid level.
        ! Hence, create the velocity vector as a restriction of 
        ! the velocity vector from a coarser level
        call solver_restrictionBlock(&
            rproblemLevel%p_rproblemLevelFine%rtriangulation,&
            rproblemLevel%rtriangulation,&
            rproblemLevel%p_rproblemLevelFine%rvectorBlock(CDEQ_VELOCITY),&
            rproblemLevel%rvectorBlock(CDEQ_VELOCITY))
        
      end if
      
      ! Switch to next coarser level
      rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
    end do

    ! Mark velocity for update
    bvelocityUpdate = .true.
  end subroutine codire_initVelocity

  !*****************************************************************************

!<subroutine>
  
  subroutine codire_initDiffusion1d(Ddiffusion)

!<description>
    ! This subroutine initializes the diffusion matrix
    !
    !     A = d11
!</description>

!<input>
    ! Matrix of diffusion coefficients
    real(DP), dimension(1,1), intent(IN) :: Ddiffusion
!</input>
!</subroutine>

    select case(idiffusiontype)
    case (DIFF_NONE)
    
      DdiffusionMatrix1D = 0._DP

      
    case (DIFF_ISOTROPIC,&
          DIFF_ANISOTROPIC)
      
      DdiffusionMatrix1D = Ddiffusion
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initDiffusion1d')
      call sys_halt()
    end select
  end subroutine codire_initDiffusion1d

  !*****************************************************************************

!<subroutine>
  
  subroutine codire_initDiffusion2d(Ddiffusion, drotation)

!<description>
    ! This subroutine initializes the (anisotropic) diffusion matrix
    !
    !     A = (  cos t  sin t )  ( d11 d12 ) ( cos t  -sin t )
    !         ( -sin t  cos t )  ( d21 d22 ) ( sin t   cos t )
!</description>

!<input>
    ! Matrix of diffusion coefficients
    real(DP), dimension(2,2), intent(IN) :: Ddiffusion
    
    ! Rotation angle in DEG
    real(DP), intent(IN) :: drotation
!</input>
!</subroutine>

    ! local variables
    real(DP), dimension(2,2) :: Rot,invRot
    real(DP) :: t

    select case(idiffusiontype)
    case (DIFF_NONE)
    
      DdiffusionMatrix2D      = 0._DP

      
    case (DIFF_ISOTROPIC)
      
      DdiffusionMatrix2D      = 0._DP
      DdiffusionMatrix2D(1,1) = maxval(Ddiffusion)
      
      
    case (DIFF_ANISOTROPIC)

      ! Convert angle to RAD
      t = 2.0_DP*SYS_PI*drotation/360.0_DP
     
      ! Compute auxiliary matrices
      invRot(1,1) =  cos(t)
      invRot(2,1) = -sin(t)
      invRot(1,2) =  sin(t)
      invRot(2,2) =  cos(t)

      Rot(1,1) =  cos(t)
      Rot(2,1) =  sin(t)
      Rot(1,2) = -sin(t)
      Rot(2,2) =  cos(t)

      ! Compute diffusion matrix
      DdiffusionMatrix2D = matmul(invRot, matmul(Ddiffusion, Rot))
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_initDiffusion2d')
      call sys_halt()
    end select
  end subroutine codire_initDiffusion2d

  !*****************************************************************************

!<subroutine>
  
  subroutine codire_initDiffusion3d(Ddiffusion, drotation)

!<description>
    ! This subroutine initializes the (anisotropic) diffusion matrix
    !
    !     A = (  cos t  sin t )  ( d11 d12 ) ( cos t  -sin t )
    !         ( -sin t  cos t )  ( d21 d22 ) ( sin t   cos t )
!</description>

!<input>
    ! Matrix of diffusion coefficients
    real(DP), dimension(3,3), intent(IN) :: Ddiffusion
    
    ! Rotation angle in DEG
    real(DP), intent(IN) :: drotation
!</input>
!</subroutine>

    print *, "Not implemented yet"
    stop
    
  end subroutine codire_initDiffusion3d

end module codire_init
