!##############################################################################
!# ****************************************************************************
!# <name> euler_init </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to initialize
!# the solver for the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!# 1.) euler_initProblem
!#     -> initialize the global problem structure
!#
!# 2.) euler_initParser
!#     -> initialize the global function parser
!#
!# 3.) euler_initConstOperators1d
!#     -> initialize constant coefficient operators for the finite
!#        element discretization in 1D
!#
!# 4.) euler_initConstOperators2d
!#     -> initialize constant coefficient operators for the finite
!#        element discretization in 2D
!#
!# 5.) euler_initConstOperators3d
!#     -> initialize constant coefficient operators for the finite
!#        element discretization in 3D
!#
!# 6.) euler_initSolutionFromParser
!#     -> initialize the solution profile at time t=0
!#
!# 7.) euler_initSolutionFromImage
!#     -> initialize the solution profile at time t=0
!#
!# 8.) euler_initExactSolution
!#     -> initialize the exact solution profile
!#
!# </purpose>
!##############################################################################

module euler_init

  use afcstabilisation
  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocsolution
  use spatialdiscretisation
  use stdoperators
  use storage
  use triangulation

  use boundaryfilter
  use euler_basic
  use euler_callback1d
  use euler_callback2d
  use euler_callback3d
  use problem


  implicit none

  private
  public :: euler_initProblem
  public :: euler_initParser
  public :: euler_initConstOperators1d
  public :: euler_initConstOperators2d
  public :: euler_initConstOperators3d
  public :: euler_initSolutionFromParser
  public :: euler_initSolutionFromImage
  public :: euler_initExactSolution


contains
  
  !*****************************************************************************

!<subroutine>

  subroutine euler_initProblem(rproblem, trifile, prmfile, indatfile,&
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
                       OU_CLASS_WARNING,OU_MODE_STD,'euler_initProblem')
      call sys_halt()
    end select

    ! Refine coarse mesh to minimum multigrid level and create standard mesh
    call tria_quickRefine2LevelOrdering (nlmin-1, rproblemLevel%rtriangulation,&
                                         rproblem%rboundary)
    call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation, rproblem%rboundary)
    
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
    do while ((associated(p_rproblemLevel)) .and.&
              (associated(p_rproblemLevel%p_rproblemLevelCoarse)))
      call tria_compress2LevelOrdHierarchy(p_rproblemLevel%rtriangulation,&
                                           p_rproblemLevel%p_rproblemLevelCoarse%rtriangulation)
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine euler_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine euler_initParser(sfilename)

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
                       OU_CLASS_WARNING,OU_MODE_STD,'euler_initParser')
      call sys_halt()
    end if

    ! Read through the complete input file and look for global
    ! definitions of constants and fixed expressions
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
  end subroutine euler_initParser

  !*****************************************************************************

!<subroutine>

  subroutine euler_initConstOperators1d(rproblem, ieltype, imatrixFormat,&
                                     nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine initializes the discrete operators resulting
    ! from the finite element method in 1D.
    ! The two optional parameters NLMIN and NLMAX can be used to
    ! restrict the range of multigrid levels for which the initialization
    ! should be performed. If they are omitted, then all levels from
    ! coarsest to finest grid are processed.
!</description>

!<input>
    ! global problem structure
    type(t_problem), intent(IN) :: rproblem

    ! type of finite elements
    integer, intent(IN) :: ieltype

    ! type of matrix format
    !   LSYSSC_MATRIX7 or LSYSSC_MATRIX9
    integer, intent(IN) :: imatrixFormat

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>
!</subroutine>

    ! section name for primal inviscid stabilization
    character(LEN=SYS_STRLEN) :: sprimalinviscidName

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel
    integer :: nlmin,nlmax,ivar,jvar,i


    ! Set minimal/maximal levels:
    ! If given by the user, adopt these values.
    ! Otherwise, assemble the constant matrices for 
    ! all levels of the multigrid structure.
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
    
    ! Assemble the divergece and global operator
    ! for all multigrid levels
    rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(rproblemLevel))
      
      ! Do we have to assemble operators for this level?
      if (rproblemLevel%ilev < nlmin .or. rproblemLevel%ilev > nlmax) then
        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
        cycle
      end if
      
      ! Check if triangulation is 1D
      if (rproblemLevel%rtriangulation%ndim .ne. NDIM1D) then
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_WARNING, OU_MODE_STD, 'euler_initConstOperators1d')
        call sys_halt()
      end if

      ! Initialize the discretization structure
      select case(isystemformat)

      case (SYSTEM_INTERLEAVEFORMAT)
        call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, 1,&
                                     rproblemLevel%rtriangulation)
        
      case (SYSTEM_BLOCKFORMAT)
        call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, NVAR1D,&
                                     rproblemLevel%rtriangulation)
        
      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators1d')
        call sys_halt()
      end select
      

      select case(ieltype)
      case (-1,1,11)
        ! P1 finite elements
        call spdiscr_initDiscr_simple (rproblemLevel%rdiscretisation%RspatialDiscr(1), &
                                       EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
                                       rproblemLevel%rtriangulation,&
                                       rproblemLevel%p_rproblem%rboundary)

      case DEFAULT
        call output_line('Unsupported element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators1d')
        call sys_halt()
      end select
      
      ! Duplicate scalar discretisation structure
      ! if system matrix is stored as block matrix.
      if (isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, NVAR1D
          call spdiscr_duplicateDiscrSc(rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                        rproblemLevel%rdiscretisation%RspatialDiscr(ivar), .true.)
        end do
      end if


      ! Allocate memory for 5 FE-matrix structures, that is
      !   the template matrix (1x),
      !   the global system matrix A (1x),
      !   the consistent/lumped mass matrices (2x),
      !   the spatial derivatives (1x),
      if (associated(rproblemLevel%rmatrix)) then
        do i = lbound(rproblemLevel%rmatrix,1),&
               ubound(rproblemLevel%rmatrix,1)
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(i))
        end do
        deallocate(rproblemLevel%rmatrix)
      end if
      allocate(rproblemLevel%Rmatrix(CNSE_MATRIX_CX))
      
      ! Allocate memory for 1 FE-block matrix structure
      if (associated(rproblemLevel%RmatrixBlock)) then
        do i = lbound(rproblemLevel%RmatrixBlock,1),&
               ubound(rproblemLevel%RmatrixBlock,1)
          call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(i))
        end do
        deallocate(rproblemLevel%RmatrixBlock)
      end if
      allocate(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))


      ! Allocate memory for 2 stabilisation structures
      if (associated(rproblemLevel%Rafcstab)) then
        do i = lbound(rproblemLevel%Rafcstab,1),&
               ubound(rproblemLevel%Rafcstab,1)
          call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
        end do
        deallocate(rproblemLevel%Rafcstab)
      end if
      allocate(rproblemLevel%Rafcstab(CNSE_AFCSTAB_VISCOUS))

      ! Initialize AFC stabilisation for inviscid part
      call parlst_getvalue_string(rparlist, '', "primalinviscid",  sprimalinviscidName)
      call afcstab_initFromParameterlist(rparlist, sprimalinviscidName,&
                                         rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID))


      ! Generate the finite element matrix sparsity structure
      call bilf_createMatrixStructure(rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                      imatrixFormat, rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE))

      ! Generate structure for scalar matrices MC and CX
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      
      ! Generate the global coefficient matrices for the FE method
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                       DER_FUNC1D, DER_FUNC1D)
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                       DER_DERIV1D_X, DER_FUNC1D)
      
      ! Perform mass lumping ML:=lumped(MC)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
      call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                   LSYSSC_LUMP_DIAG)
      

      ! Generate structure for system matrix.
      select case(isystemFormat)

      case (SYSTEM_INTERLEAVEFORMAT)
        ! The global operator is stored as an interleave matrix with 
        ! NVAR1D components. However, the row and column structure of
        ! the template matrix can be adopted
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
        rproblemLevel%Rmatrix(CNSE_MATRIX_A)%NVAR = NVAR1D

        
        ! What matrix format should be used
        select case(imatrixFormat)
        case (LSYSSC_MATRIX7)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cmatrixFormat = LSYSSC_MATRIX7INTL
          
        case (LSYSSC_MATRIX9)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cmatrixFormat = LSYSSC_MATRIX9INTL
          
        case DEFAULT
          call output_line('Unsupported matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators1d')
          call sys_halt()
        end select

        
        ! What kind of global operator should be adopted?
        select case(isystemCoupling)
        case (SYSTEM_SEGREGATED)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cinterleavematrixFormat = LSYSSC_MATRIXD
          
        case (SYSTEM_ALLCOUPLED)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cinterleavematrixFormat = LSYSSC_MATRIX1
          
        case DEFAULT
          call output_line('Unsupported interleave matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators1d')
          call sys_halt()
        end select

        
        ! Create global operator physically
        call lsyssc_allocEmptyMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                     LSYSSC_SETM_UNDEFINED)
      
        ! Create pseudo block matrix from global operator
        call lsysbl_createMatFromScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                        rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                        rproblemLevel%rdiscretisation)


      case (SYSTEM_BLOCKFORMAT)
        ! The global operator is stored as a block matrix with NVAR1D x NVAR1D
        ! blocks. Thus, create empty NVAR1D x NVAR1D block matrix directly
        call lsysbl_createEmptyMatrix(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A), NVAR1D)

        ! Define global operator as groupmatrix
        rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX

        ! What kind of global operator should be adopted?
        select case(isystemCoupling)
        case (SYSTEM_SEGREGATED)
          
          ! Create diagonal blocks
          do ivar = 1, NVAR1D
            call lsyssc_duplicateMatrix(&
                rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,ivar),&
                LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

          end do

        case (SYSTEM_ALLCOUPLED)

          ! Create all blocks
          do ivar = 1, NVAR1D
            do jvar = 1, NVAR1D
              call lsyssc_duplicateMatrix(&
                  rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                  rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,jvar),&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
              
            end do
          end do
          
        case DEFAULT
          call output_line('Unsupported block matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators1d')
          call sys_halt()
        end select

        ! Update internal structure of block matrix
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))


      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators1d')
        call sys_halt()
      end select
      
      ! Switch to next coarser level
      rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine euler_initConstOperators1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_initConstOperators2d(rproblem, ieltype, imatrixFormat,&
                                     nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine initializes the discrete operators resulting
    ! from the finite element method in 2D.
    ! The two optional parameters NLMIN and NLMAX can be used to
    ! restrict the range of multigrid levels for which the initialization
    ! should be performed. If they are omitted, then all levels from
    ! coarsest to finest grid are processed.
!</description>

!<input>
    ! global problem structure
    type(t_problem), intent(IN) :: rproblem

    ! type of finite elements
    integer, intent(IN) :: ieltype

    ! type of matrix format
    !   LSYSSC_MATRIX7 or LSYSSC_MATRIX9
    integer, intent(IN) :: imatrixFormat

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>
!</subroutine>

    ! section name for primal inviscid stabilization
    character(LEN=SYS_STRLEN) :: sprimalinviscidName

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel
    integer :: nlmin,nlmax,ivar,jvar,i


    ! Set minimal/maximal levels:
    ! If given by the user, adopt these values.
    ! Otherwise, assemble the constant matrices for 
    ! all levels of the multigrid structure.
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
    
    ! Assemble the divergece and global operator
    ! for all multigrid levels
    rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(rproblemLevel))
      
      ! Do we have to assemble operators for this level?
      if (rproblemLevel%ilev < nlmin .or. rproblemLevel%ilev > nlmax) then
        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
        cycle
      end if
      
      ! Check if triangulation is 1D
      if (rproblemLevel%rtriangulation%ndim .ne. NDIM1D) then
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_WARNING, OU_MODE_STD, 'euler_initConstOperators1d')
        call sys_halt()
      end if

      ! Initialize the discretization structure
      select case(isystemformat)

      case (SYSTEM_INTERLEAVEFORMAT)
        call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, 1,&
                                     rproblemLevel%rtriangulation)
        
      case (SYSTEM_BLOCKFORMAT)
        call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, NVAR2D,&
                                     rproblemLevel%rtriangulation)
        
      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators2d')
        call sys_halt()
      end select
      

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
                                        EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC,&
                                        rproblemLevel%rtriangulation,&
                                        rproblemLevel%p_rproblem%rboundary)

      case DEFAULT
        call output_line('Unsupported element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators2d')
        call sys_halt()
      end select
      
      ! Duplicate scalar discretisation structure
      ! if system matrix is stored as block matrix.
      if (isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, NVAR2D
          call spdiscr_duplicateDiscrSc(rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                        rproblemLevel%rdiscretisation%RspatialDiscr(ivar), .true.)
        end do
      end if


      ! Allocate memory for 6 FE-matrix structures, that is
      !   the template matrix (1x),
      !   the global system matrix A (1x),
      !   the consistent/lumped mass matrices (2x),
      !   the spatial derivatives (2x),
      if (associated(rproblemLevel%rmatrix)) then
        do i = lbound(rproblemLevel%rmatrix,1),&
               ubound(rproblemLevel%rmatrix,1)
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(i))
        end do
        deallocate(rproblemLevel%rmatrix)
      end if
      allocate(rproblemLevel%Rmatrix(CNSE_MATRIX_CY))
      
      ! Allocate memory for 1 FE-block matrix structure
      if (associated(rproblemLevel%RmatrixBlock)) then
        do i = lbound(rproblemLevel%RmatrixBlock,1),&
               ubound(rproblemLevel%RmatrixBlock,1)
          call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(i))
        end do
        deallocate(rproblemLevel%RmatrixBlock)
      end if
      allocate(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))


      ! Allocate memory for 2 stabilisation structures
      if (associated(rproblemLevel%Rafcstab)) then
        do i = lbound(rproblemLevel%Rafcstab,1),&
               ubound(rproblemLevel%Rafcstab,1)
          call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
        end do
        deallocate(rproblemLevel%Rafcstab)
      end if
      allocate(rproblemLevel%Rafcstab(CNSE_AFCSTAB_VISCOUS))

      ! Initialize AFC stabilisation for inviscid part
      call parlst_getvalue_string(rparlist, '', "primalinviscid",  sprimalinviscidName)
      call afcstab_initFromParameterlist(rparlist, sprimalinviscidName,&
                                         rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID))


      ! Generate the finite element matrix sparsity structure
      call bilf_createMatrixStructure(rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                      imatrixFormat, rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE))

      ! Generate structure for scalar matrices MC, CX, and CY
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      
      ! Generate the global coefficient matrices for the FE method
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                       DER_FUNC2D, DER_FUNC2D)
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                       DER_DERIV2D_X, DER_FUNC2D)
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                                       DER_DERIV2D_Y, DER_FUNC2D)
      
      ! Perform mass lumping ML:=lumped(MC)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
      call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                   LSYSSC_LUMP_DIAG)
      

      ! Generate structure for system matrix.
      select case(isystemFormat)

      case (SYSTEM_INTERLEAVEFORMAT)
        ! The global operator is stored as an interleave matrix with 
        ! NVAR2D components. However, the row and column structure of
        ! the template matrix can be adopted
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
        rproblemLevel%Rmatrix(CNSE_MATRIX_A)%NVAR = NVAR2D

        
        ! What matrix format should be used
        select case(imatrixFormat)
        case (LSYSSC_MATRIX7)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cmatrixFormat = LSYSSC_MATRIX7INTL
          
        case (LSYSSC_MATRIX9)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cmatrixFormat = LSYSSC_MATRIX9INTL
          
        case DEFAULT
          call output_line('Unsupported matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators2d')
          call sys_halt()
        end select

        
        ! What kind of global operator should be adopted?
        select case(isystemCoupling)
        case (SYSTEM_SEGREGATED)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cinterleavematrixFormat = LSYSSC_MATRIXD
          
        case (SYSTEM_ALLCOUPLED)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cinterleavematrixFormat = LSYSSC_MATRIX1
          
        case DEFAULT
          call output_line('Unsupported interleave matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators2d')
          call sys_halt()
        end select

        
        ! Create global operator physically
        call lsyssc_allocEmptyMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                     LSYSSC_SETM_UNDEFINED)
      
        ! Create pseudo block matrix from global operator
        call lsysbl_createMatFromScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                        rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                        rproblemLevel%rdiscretisation)


      case (SYSTEM_BLOCKFORMAT)
        ! The global operator is stored as a block matrix with NVAR2D x NVAR2D
        ! blocks. Thus, create empty NVAR2D x NVAR2D block matrix directly
        call lsysbl_createEmptyMatrix(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A), NVAR2D)

        ! Define global operator as groupmatrix
        rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX

        ! What kind of global operator should be adopted?
        select case(isystemCoupling)
        case (SYSTEM_SEGREGATED)
          
          ! Create diagonal blocks
          do ivar = 1, NVAR2D
            call lsyssc_duplicateMatrix(&
                rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,ivar),&
                LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

          end do

        case (SYSTEM_ALLCOUPLED)

          ! Create all blocks
          do ivar = 1, NVAR2D
            do jvar = 1, NVAR2D
              call lsyssc_duplicateMatrix(&
                  rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                  rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,jvar),&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
              
            end do
          end do
          
        case DEFAULT
          call output_line('Unsupported block matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators2d')
          call sys_halt()
        end select

        ! Update internal structure of block matrix
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))


      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators2d')
        call sys_halt()
      end select
      
      ! Switch to next coarser level
      rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine euler_initConstOperators2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_initConstOperators3d(rproblem, ieltype, imatrixFormat,&
                                     nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine initializes the discrete operators resulting
    ! from the finite element method in 3D. 
    ! The two optional parameters NLMIN and NLMAX can be used to 
    ! restrict the range of multigrid levels for which the initialization
    ! should be performed. If they are omitted, then all levels from 
    ! coarsest to finest grid are processed. 
!</description>

!<input>
    ! global problem structure
    type(t_problem), intent(IN) :: rproblem

    ! type of finite elements
    integer, intent(IN) :: ieltype

    ! type of matrix format
    !   LSYSSC_MATRIX7 or LSYSSC_MATRIX9
    integer, intent(IN) :: imatrixFormat

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>
!</subroutine>

    ! section name for primal inviscid stabilization
    character(LEN=SYS_STRLEN) :: sprimalinviscidName

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel
    integer :: nlmin,nlmax,ivar,jvar,i


    ! Set minimal/maximal levels:
    ! If given by the user, adopt these values.
    ! Otherwise, assemble the constant matrices for 
    ! all levels of the multigrid structure.
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
    
    ! Assemble the divergece and global operator
    ! for all multigrid levels
    rproblemLevel => rproblem%p_rproblemLevelMax
    do while(associated(rproblemLevel))
      
      ! Do we have to assemble operators for this level?
      if (rproblemLevel%ilev < nlmin .or. rproblemLevel%ilev > nlmax) then
        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
        cycle
      end if
      
      ! Check if triangulation is 3D
      if (rproblemLevel%rtriangulation%ndim .ne. NDIM3D) then
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_WARNING, OU_MODE_STD, 'euler_initConstOperators3d')
        call sys_halt()
      end if
      
      ! Initialize the discretization structure
      select case(isystemformat)

      case (SYSTEM_INTERLEAVEFORMAT)
        call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, 1,&
                                     rproblemLevel%rtriangulation)
        
      case (SYSTEM_BLOCKFORMAT)
        call spdiscr_initBlockDiscr (rproblemLevel%rdiscretisation, NVAR3D,&
                                     rproblemLevel%rtriangulation)

      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators3d')
        call sys_halt()
      end select
      

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
        call output_line('Unsupported element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators3d')
        call sys_halt()
      end select

      ! Duplicate scalar discretisation structure
      if (isystemFormat .eq. SYSTEM_BLOCKFORMAT) then
        do ivar = 2, NVAR3D
          call spdiscr_duplicateDiscrSc(rproblemLevel%rdiscretisation%RspatialDiscr(1),&
                                        rproblemLevel%rdiscretisation%RspatialDiscr(ivar), .true.)
        end do
      end if


      ! Allocate memory for 7 FE-matrix structures, that is
      !   the template matrix (1x),
      !   the global system matrix A (1x),
      !   the consistent/lumped mass matrices (2x),
      !   the spatial derivatives (3x),
      if (associated(rproblemLevel%rmatrix)) then
        do i = lbound(rproblemLevel%rmatrix,1),&
               ubound(rproblemLevel%rmatrix,1)
          call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(i))
        end do
        deallocate(rproblemLevel%rmatrix)
      end if
      allocate(rproblemLevel%Rmatrix(CNSE_MATRIX_CZ))
      
      ! Allocate memory for 1 FE-block matrix structure
      if (associated(rproblemLevel%RmatrixBlock)) then
        do i = lbound(rproblemLevel%RmatrixBlock,1),&
               ubound(rproblemLevel%RmatrixBlock,1)
          call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(i))
        end do
        deallocate(rproblemLevel%RmatrixBlock)
      end if
      allocate(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))


      ! Allocate memory for 2 stabilisation structures
      if (associated(rproblemLevel%Rafcstab)) then
        do i = lbound(rproblemLevel%Rafcstab,1),&
               ubound(rproblemLevel%Rafcstab,1)
          call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
        end do
        deallocate(rproblemLevel%Rafcstab)
      end if
      allocate(rproblemLevel%Rafcstab(CNSE_AFCSTAB_VISCOUS))

      ! Initialize AFC stabilisation for inviscid part
      call parlst_getvalue_string(rparlist, '', "primalinviscid",  sprimalinviscidName)
      call afcstab_initFromParameterlist(rparlist, sprimalinviscidName,&
                                         rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID))


      ! Generate the finite element matrix sparsity structure
      call bilf_createMatrixStructure(&
          rproblemLevel%rdiscretisation%RspatialDiscr(1),&
          imatrixFormat, rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE))

      ! Generate structure for scalar matrices MC, CX, CY, and CZ
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_CZ),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      
      ! Generate the global coefficient matrices for the FE method
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                       DER_FUNC3D, DER_FUNC3D)
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_CX),&
                                       DER_DERIV3D_X, DER_FUNC3D)
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_CY),&
                                       DER_DERIV3D_Y, DER_FUNC3D)
      call stdop_assembleSimpleMatrix (rproblemLevel%Rmatrix(CNSE_MATRIX_CZ),&
                                       DER_DERIV3D_Z, DER_FUNC3D)
      
      ! Perform mass lumping ML:=lumped(MC)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_MC),&
                                  rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
      call lsyssc_lumpMatrixScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                   LSYSSC_LUMP_DIAG)
      

      ! Generate structure for system matrix.
      select case(isystemFormat)

      case (SYSTEM_INTERLEAVEFORMAT)
        ! The global operator is stored as an interleave matrix with 
        ! NVAR3D components. However, the row and column structure of
        ! the template matrix can be adopted
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                                    rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
        rproblemLevel%Rmatrix(CNSE_MATRIX_A)%NVAR = NVAR3D

        
        ! What matrix format should be used
        select case(imatrixFormat)
        case (LSYSSC_MATRIX7)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cmatrixFormat = LSYSSC_MATRIX7INTL
          
        case (LSYSSC_MATRIX9)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cmatrixFormat = LSYSSC_MATRIX9INTL
          
        case DEFAULT
          call output_line('Unsupported matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators3d')
          call sys_halt()
        end select

        
        ! What kind of global operator should be adopted?
        select case(isystemCoupling)
        case (SYSTEM_SEGREGATED)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cinterleavematrixFormat = LSYSSC_MATRIXD
          
        case (SYSTEM_ALLCOUPLED)
          rproblemLevel%Rmatrix(CNSE_MATRIX_A)%cinterleavematrixFormat = LSYSSC_MATRIX1
          
        case DEFAULT
          call output_line('Unsupported interleave matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators3d')
          call sys_halt()
        end select

        
        ! Create global operator physically
        call lsyssc_allocEmptyMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                     LSYSSC_SETM_UNDEFINED)
      
        ! Create pseudo block matrix from global operator
        call lsysbl_createMatFromScalar(rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                        rproblemLevel%RmatrixBlock(CNSE_MATRIX_A),&
                                        rproblemLevel%rdiscretisation)


      case (SYSTEM_BLOCKFORMAT)
        ! The global operator is stored as a block matrix with NVAR3D x NVAR3D
        ! blocks. Thus, create empty NVAR3D x NVAR3D block matrix directly
        call lsysbl_createEmptyMatrix(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A), NVAR3D)
        rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%imatrixSpec = LSYSBS_MSPEC_GROUPMATRIX

        ! What kind of global operator should be adopted?
        select case(isystemCoupling)
        case (SYSTEM_SEGREGATED)
          
          ! Create diagonal blocks
          do ivar = 1, NVAR3D
            call lsyssc_duplicateMatrix(&
                rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,ivar),&
                LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
            call lsyssc_allocEmptyMatrix(&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,ivar),&
                LSYSSC_SETM_UNDEFINED)
          end do

        case (SYSTEM_ALLCOUPLED)

          ! Create all blocks
          do ivar = 1, NVAR3D
            do jvar = 1, NVAR3D
              call lsyssc_duplicateMatrix(&
                  rproblemLevel%Rmatrix(CNSE_MATRIX_TEMPLATE),&
                  rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,jvar),&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
              call lsyssc_allocEmptyMatrix(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(ivar,jvar),&
                                           LSYSSC_SETM_UNDEFINED)
            end do
          end do
          
        case DEFAULT
          call output_line('Unsupported block matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators3d')
          call sys_halt()
        end select

        ! Update structure of block matrix
        call lsysbl_updateMatStrucInfo(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A))

      case DEFAULT
        call output_line('Unsupported system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_initConstOperators3d')
        call sys_halt()
      end select
      
      ! Switch to next coarser level
      rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
    end do
  end subroutine euler_initConstOperators3d
  
  !*****************************************************************************

!<subroutine>

  subroutine euler_initSolutionFromParser(rboundaryCondition, rproblemLevel,&
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
    integer(I32) :: isize
    integer :: istatus

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem

    ! Allocate storage for solution vector
    call lsysbl_createVectorBlock(rproblemLevel%rdiscretisation,&
                                  rsolution, .false., ST_DOUBLE)

    ! Ok, if the vector is stored in interleave format, 
    ! then some modifications are required
    if (isystemFormat .eq. SYSTEM_INTERLEAVEFORMAT) then

      ! Indicate number of variables in scalar vector:
      ! spatial dimension + 2 
      rsolution%RvectorBlock(1)%NVAR = euler_getNVAR(rproblemLevel)

      ! Resize block vector accordingly
      isize = rsolution%NEQ*euler_getNVAR(rproblemLevel)
      call lsysbl_resizeVectorBlock(rsolution, isize, .false., .false.)
    end if


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
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initSolutionFromParser')
      call sys_halt()
    end select
  
    ! Check status
    if (istatus .ne. 0) then
      call output_line('Unable to create initial solution profile',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initSolutionFromParser')
      call sys_halt()
    end if
      
    ! Release parser
    call fparser_release(rparser)
  
  
    ! Impose boundary conditions explicitely
    select case(rproblemLevel%rdiscretisation%ndimension)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rboundaryCondition, rproblemLevel%rtriangulation,&
                                     rsolution, 0._DP, p_rproblem%rboundary,&
                                     fcb_calcBoundaryvalues1d, istatus)
    case (NDIM2D)
      call bdrf_filterVectorExplicit(rboundaryCondition, rproblemLevel%rtriangulation,&
                                     rsolution, 0._DP, p_rproblem%rboundary,&
                                     fcb_calcBoundaryvalues2d, istatus)
    case (NDIM3D)
      call bdrf_filterVectorExplicit(rboundaryCondition, rproblemLevel%rtriangulation,&
                                     rsolution, 0._DP, p_rproblem%rboundary,&
                                     fcb_calcBoundaryvalues3d, istatus)
    end select
    
    ! Check status
    if (istatus .ne. 0) then
      call output_line('Unable to impose boundary coditions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initSolutionFromParser')
      call sys_halt()
    end if
  end subroutine euler_initSolutionFromParser

  !*****************************************************************************

!<subroutine>

  subroutine euler_initSolutionFromImage(rboundaryCondition, rproblemLevel,&
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
    type(t_pgm) :: rpgm
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    integer(I32) :: isize
    integer :: istatus

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem

    ! Check if problem is in 2D
    if (rproblemLevel%rdiscretisation%ndimension .ne. NDIM2D) then
      call output_line('Initialization from image is only available in 2D!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initSolutionFromImage')
      call sys_halt()
    end if

    ! Allocate storage for solution vector
    call lsysbl_createVectorBlock(rproblemLevel%rdiscretisation,&
                                  rsolution, .true., ST_DOUBLE)

    ! Ok, if the vector is stored in interleave format, 
    ! then some modifications are required
    if (isystemFormat .eq. SYSTEM_INTERLEAVEFORMAT) then
      
      ! Indicate number of variables in scalar vector: NVAR2D
      rsolution%RvectorBlock(1)%NVAR = NVAR2D
      
      ! Resize block vector accordingly
      isize = rsolution%NEQ*NVAR2D
      call lsysbl_resizeVectorBlock(rsolution, isize, .false., .false.)
    end if

    ! Read PGM image from file
    call ppsol_readPGM(0, pgmfile, rpgm)

    call storage_getbase_double2D(rproblemLevel%rtriangulation%h_DvertexCoords,&
                                  p_DvertexCoords)
    call storage_getbase_double(rsolution%h_Ddata, p_Ddata)
    call ppsol_initArrayPGM_Dble(rpgm, p_DvertexCoords, p_Ddata)
    
    ! Release PGM image
    call ppsol_releasePGM(rpgm)

    ! Impose boundary conditions explicitely
    call bdrf_filterVectorExplicit(rboundaryCondition, rproblemLevel%rtriangulation,&
                                   rsolution, 0._DP, p_rproblem%rboundary,&
                                   fcb_calcBoundaryvalues2d, istatus)

    ! Check status
    if (istatus .ne. 0) then
      call output_line('Unable to impose boundary coditions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initSolutionFromImage')
      call sys_halt()
    end if
  end subroutine euler_initSolutionFromImage

  !*****************************************************************************

!<subroutine>

  subroutine euler_initExactSolution(rproblemLevel, rsolution, indatfile, time, istatus)

!<description>
    ! This subroutine initializes the exact solution vector from
    ! the parameters specified in the indatfile by calling
    ! a function parser w.r.t. the specified variable names
!</description>
    
!<input>
    ! multigrid level structure
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
    type(t_fparser):: rparser
    integer(I32) :: isize
    integer:: iistatus

    if (present(istatus)) istatus = 0

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem

    ! Allocate storage for solution vector
    call lsysbl_createVectorBlock(rproblemLevel%rdiscretisation,&
                                  rsolution, .false., ST_DOUBLE)

    ! Ok, if the vector is stored in interleave format, 
    ! then some modifications are required
    if (isystemFormat .eq. SYSTEM_INTERLEAVEFORMAT) then
      
      ! Indicate number of variables in scalar vector:
      ! spatial dimension + 2
      rsolution%RvectorBlock(1)%NVAR = euler_getNVAR(rproblemLevel)

      ! Resize block vector accordingly
      isize = rsolution%NEQ*euler_getNVAR(rproblemLevel)
      call lsysbl_resizeVectorBlock(rsolution, isize, .false., .false.)
    end if

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
      call problem_createProfile(p_rproblem, indatfile, '[exact_solution]',&
                                 SOL_SYMBOLICVARS3D, rparser, iistatus, time, rsolution)
    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initExactSolution')
      call sys_halt()
    end select
    

    ! No exact solution available; free memory and return
    if (iistatus .ne. 0) then
      call output_line('Unable to create exact solution profile',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_initExactSolution')
      call lsysbl_releaseVector(rsolution)
      if (present(istatus)) istatus = -1
      return
    end if

    ! Release parser
    call fparser_release(rparser)
  end subroutine euler_initExactSolution
end module euler_init
