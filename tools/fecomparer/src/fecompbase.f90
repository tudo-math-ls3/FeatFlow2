module fecompbase

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  
  use boundary
  use triangulation
  
  use element
  use linearsystemblock
  
  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  
  use paramlist
  use collection
  
  use ucd
  
  use solutions
  use q1projection

  implicit none
  private

  ! Encapsules an FE space hierarchy with a tag
  type t_feHierarchyMapEntry
  
    ! FE space tag
    ! =   0: 2D mesh, Scalar solution, P0
    ! =   1: 2D mesh, Scalar solution, P1
    ! =  11: 2D mesh, Scalar solution, Q1
    ! =  13: 2D mesh, Scalar solution, Q2
    ! =  14: 2D mesh, Scalar solution, Q3
    ! =1010: 2D mesh, 2D Velocity/pressure, Q0/Q0/Q0
    ! =1011: 2D mesh, 2D Velocity/pressure, Q1/Q1/Q1
    ! =1012: 2D mesh, 2D Velocity/pressure, Q2/Q2/QP1
    ! =1013: 2D mesh, 2D Velocity/pressure, Q3/Q3/QP2 (nonparametric)
    ! =1030: 2D mesh, 2D Velocity/pressure, E030/E030/Q0
    ! =1031: 2D mesh, 2D Velocity/pressure, E031/E031/Q0
    ! =1130: 2D mesh, 2D Velocity/pressure, EM30/EM30/Q0
    ! =1131: 2D mesh, 2D Velocity/pressure, EM31/EM31/Q0
    ! =1211: 2D mesh, 2D Velocity/pressure, Q1/Q1/Q1, all discontinuous
    integer :: cfespace
    
    ! The FE space hierarchy
    type(t_feHierarchy) :: rfeHierarchy

  end type
  
  public :: t_feHierarchyMapEntry
  
  ! Problem structure
  type t_problem
  
    !<!-- Parameters -->
  
    ! Name of the underlying parametrisation
    character(LEN=SYS_STRLEN) :: sparametrisation = ""

    ! Name of the underlying coarse mesh
    character(LEN=SYS_STRLEN) :: striangulation = ""

    ! Type of the underlying mesh.
    ! =1: 1D
    ! =2: 2D, given in sparametrisation and striangulation
    ! =3: 3D, given in sparametrisation
    integer :: cmeshtype = 2

    ! Number of solutions to read.
    ! Each solution has its own section [SOLUTIONx].
    integer :: nsolutions = 0

    ! Algebraic computation to apply to the solutions
    ! =0: result=SOLUTION1
    ! =1: result=SOLUTION2
    ! =2: result=SOLUTION2 - SOLUTION1
    integer :: ccomputation = 0

    ! Type of analysis.
    ! =0: Stationary analysis
    ! =1: Nonstationary analysis
    integer :: canalysis = 0

    ! Type of output.
    ! =0: No output
    ! =1: VTK file
    integer :: cpostprocfile = 0

    ! Output filename
    character(LEN=SYS_STRLEN) :: spostprocfilename = ""

    ! Compute ||result||_L2 at level ilevel
    integer :: cerrorL2 = 0

    ! Compute ||result||_LH1 at level ilevel
    integer :: cerrorH1 = 0

    ! Cubature rule to use for error computation
    character(LEN=SYS_STRLEN) :: scubature = "AUTO_G5"
  
    !<!-- Actual data -->
  
    ! A mesh hierarchy and a FE space hierarchy.
    type(t_meshHierarchy) :: rmeshHierarchy
    
    ! List of FE spaces
    type(t_feHierarchyMapEntry), dimension(32) :: RfeHierarchies
    
    ! Number of FE hierarchies
    integer :: nfeHierarchies = 0
    
    ! List of solutions. Solution1 is the main solution which
    ! is used for export, error analysis, etc.
    ! It is expected to be a Q1 solution.
    type(t_solution), dimension(32) :: Rsolutions
    
    ! Reference solution, alias "SOLUTION0"
    type(t_solution) :: rrefSolution
    
    ! Other variables
    type(t_parlist) :: rparlist
    
    ! Boundary definition
    type(t_boundary) :: rboundary
    
    ! Coarse mesh triangulation
    type(t_triangulation) :: rtriangulation
    
  
  end type
  
  public :: t_problem
  
  public :: base_getFEhierarchy
  public :: base_init
  public :: base_done
  public :: base_createMeshHierarchy
  public :: base_allocateSolutions
  public :: base_readSolutions
  
  public :: base_vectorcopy
  public :: base_vectorLinearComb
  public :: base_writeOutput
  public :: base_computeErrors
  public :: base_vectorDivide
  
contains

  ! ***************************************************************************
  subroutine fgetDiscr(ilevel,rtriangulation,rdiscr,rboundary,rcollection)
  
  ! Creates a discretisation for a triangulation.
  
  integer, intent(in) :: ilevel
  type(t_triangulation), intent(in) :: rtriangulation
  type(t_blockDiscretisation), intent(out) :: rdiscr
  type(t_collection), intent(inout), optional :: rcollection
  type(t_boundary), intent(in), optional :: rboundary
  
    integer :: cfespace
    
    cfespace = rcollection%IquickAccess(1)
    
    select case (cfespace)
    
    ! ----------------------------------------
    ! 2D mesh, Scalar solution, P0
    case (0)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_P0_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, Scalar solution, P1
    case (1)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_P1_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, Scalar solution, P2
    case (2)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_P2_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, Scalar solution, Q0
    case (10)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q0_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, Scalar solution, Q1
    case (11)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q1_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, Scalar solution, Q2
    case (12)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q2_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, Scalar solution, Q3
    case (13)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q3_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, Scalar solution, Q1
    case (111)
      call spdiscr_initBlockDiscr (rdiscr,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_DG_Q1_2D,&
          rtriangulation, rboundary)

    ! ----------------------------------------
    ! 2D mesh, 2D Velocity/pressure, Q1/Q1/Q1
    case (1011)
      call spdiscr_initBlockDiscr (rdiscr,3,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q1_2D,&
          rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (rdiscr%RspatialDiscr(1), rdiscr%RspatialDiscr(2), .true.)
      call spdiscr_duplicateDiscrSc (rdiscr%RspatialDiscr(1), rdiscr%RspatialDiscr(3), .true.)
    
    ! ----------------------------------------
    ! 2D mesh, 2D Velocity/pressure, Q2/Q2/QP1
    case (1012)
      call spdiscr_initBlockDiscr (rdiscr,3,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q2_2D,&
          rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (rdiscr%RspatialDiscr(1), rdiscr%RspatialDiscr(2), .true.)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(3),EL_QP1NP_2D,&
          rtriangulation, rboundary)
    
    ! ----------------------------------------
    ! 2D mesh, 2D Velocity/pressure, Q3/Q3/DCQP2 (nonparametric)
    case (1013)
      call spdiscr_initBlockDiscr (rdiscr,3,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_Q3_2D,&
          rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (rdiscr%RspatialDiscr(1), rdiscr%RspatialDiscr(2), .true.)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(3),EL_DCQP2_2D,&
          rtriangulation, rboundary)
    
    ! ----------------------------------------
    ! 2D mesh, 2D Velocity/pressure, Q1DG/Q1DG/Q1DG 
    case (1211)
      call spdiscr_initBlockDiscr (rdiscr,3,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(1),EL_DG_Q1_2D,&
          rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (rdiscr%RspatialDiscr(1), rdiscr%RspatialDiscr(2), .true.)
      call spdiscr_duplicateDiscrSc (rdiscr%RspatialDiscr(1), rdiscr%RspatialDiscr(3), .true.)
      
    case default
    
      call output_line ("Unknown FE space.")
      call sys_halt()

    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_getFEhierarchy (cfespace,RfeHierarchies,ncount,&
      rmeshHierarchy,iindex,rboundary)

!<description>
  ! Returns a FE space hierarchy according to the mesh hierarchy rmeshHierarchy
  ! and the FE space cfespace. If necessary, a new FE space hierarchy is created.
!</description>

!<input>
  ! Desired FE space
  integer :: cfespace
  
  ! Underlying mesh hierarchy
  type(t_meshHierarchy), intent(in) :: rmeshHierarchy
  
  ! OPTIONAL: Boundary definition
  type(t_boundary), intent(in), optional :: rboundary
!</input>

!<inputoutput>
  ! Array and number of elements of existing FE hierarchies
  type(t_feHierarchyMapEntry), dimension(:), intent(inout) :: RfeHierarchies
  integer, intent(inout) :: ncount
  
  ! Index of the hierarchy
  integer, intent(out) :: iindex
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_collection) :: rcollection
    
    ! Try to find it.
    do i=1,ncount
      if (RfeHierarchies(i)%cfespace .eq. cfespace) then
        iindex = i
        return
      end if
    end do
    
    ! Create a new one based on the mesh hierarchy.
    ncount = ncount + 1
    RfeHierarchies(ncount)%cfespace = cfespace

    rcollection%IquickAccess(1) = cfespace
    call fesph_createHierarchy (RfeHierarchies(ncount)%rfeHierarchy,rmeshHierarchy%nlevels,&
        rmeshHierarchy,fgetDiscr,rcollection,rboundary)
    
    iindex = ncount

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_init (rproblem)

!<description>
  ! Reads parameters, basic initialisation
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(inout) :: rproblem
!</input>

!</subroutine>

    integer :: i

    ! Read the parameter file
    call parlst_init (rproblem%rparlist)
    call parlst_readfromfile (rproblem%rparlist, "data/fecomparer.dat")
    
    call parlst_getvalue_string (rproblem%rparlist, "", &
        "sparametrisation", rproblem%sparametrisation, bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparlist, "", &
        "striangulation", rproblem%striangulation, bdequote=.true.)

    call parlst_getvalue_int (rproblem%rparlist, "", &
        "cmeshtype", rproblem%cmeshtype)

    call parlst_getvalue_int (rproblem%rparlist, "", &
        "nsolutions", rproblem%nsolutions)

    call parlst_getvalue_int (rproblem%rparlist, "", &
        "ccomputation", rproblem%ccomputation)

    call parlst_getvalue_int (rproblem%rparlist, "", &
        "canalysis", rproblem%canalysis)

    call parlst_getvalue_int (rproblem%rparlist, "", &
        "cpostprocfile", rproblem%cpostprocfile)

    call parlst_getvalue_string (rproblem%rparlist, "", &
        "spostprocfilename", rproblem%spostprocfilename, bdequote=.true.)

    call parlst_getvalue_int (rproblem%rparlist, "", &
        "cerrorL2", rproblem%cerrorL2)
        
    call parlst_getvalue_int (rproblem%rparlist, "", &
        "cerrorH1", rproblem%cerrorH1)

    call parlst_getvalue_string (rproblem%rparlist, "", &
        "scubature", rproblem%scubature, bdequote=.true.)
        
    ! Read the solutions
    do i=1,rproblem%nsolutions
      call sol_readparams (rproblem%rparlist,i,rproblem%Rsolutions(i))
    end do
    call sol_readparams (rproblem%rparlist,0,rproblem%rrefSolution)
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_createMeshHierarchy (rproblem)

!<description>
  ! Creates the underlying mesh hierarchy
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(inout) :: rproblem
!</input>

!</subroutine>

    ! Figure out the maximum total level
    integer :: i,nmaxLevels
    
    nmaxLevels = rproblem%rrefSolution%ireflevel
    do i=1,rproblem%nsolutions
      nmaxLevels = max(nmaxLevels,rproblem%Rsolutions(i)%ireflevel)
    end do
    
    select case (rproblem%cmeshtype)
    case (2)
    
      ! Read the coarse mesh.
      call boundary_read_prm(rproblem%rboundary, rproblem%sparametrisation)
      
      ! Read the triangulation
      call tria_readTriFile2D (&
          rproblem%rtriangulation, rproblem%striangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%rtriangulation,rproblem%rboundary)

      ! Initialise the mesh hierarchy
      call mshh_initHierarchy (rproblem%rmeshHierarchy,rproblem%rtriangulation,0,&
          nmaxLevels,rproblem%rboundary)
          
      call mshh_refineHierarchy2lv (rproblem%rmeshHierarchy,nmaxLevels,&
          rboundary=rproblem%rboundary,bprint=.true.)
          
      call output_lbrk()
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_allocateSolutions (rproblem)

!<description>
  ! Creates FE hierarchies for the solutions
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(inout), target :: rproblem
!</input>

!</subroutine>

    integer :: i,iindex
    
    do i=0,rproblem%nsolutions
    
      ! Create a FE hierarchy or get an existing one.
      select case (rproblem%cmeshtype)
      case (2)
        if (i .eq. 0) then
          call base_getFEhierarchy (rproblem%rrefSolution%cfespace,&
              rproblem%RfeHierarchies,rproblem%nfeHierarchies,&
              rproblem%rmeshHierarchy,iindex,rproblem%rboundary)
        else
          call base_getFEhierarchy (rproblem%Rsolutions(i)%cfespace,&
              rproblem%RfeHierarchies,rproblem%nfeHierarchies,&
              rproblem%rmeshHierarchy,iindex,rproblem%rboundary)
        end if

      end select
      
      ! Initialise the solution
      if (i .eq. 0) then
        call sol_initSolution (rproblem%rrefSolution,rproblem%RfeHierarchies(iindex)%rfeHierarchy)
      else
        call sol_initSolution (rproblem%Rsolutions(i),rproblem%RfeHierarchies(iindex)%rfeHierarchy)
      end if
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_readSolutions (rproblem)

!<description>
  ! Reads solutions from hard disc.
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(inout), target :: rproblem
!</input>

!</subroutine>

    integer :: i,iindex
    
    do i=0,rproblem%nsolutions
      ! Read the solution
      if (i .eq. 0) then
        call sol_readSolutionData (rproblem%rrefSolution)
      else
        call sol_readSolutionData (rproblem%Rsolutions(i))
      end if
    
    end do
    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_done (rproblem)

!<description>
  ! Reads parameters, basic initialisation
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(inout) :: rproblem
!</input>

!</subroutine>

    integer :: i

    ! Release solutions
    do i=1,rproblem%nsolutions
      call sol_doneSolution (rproblem%Rsolutions(i))
    end do
    call sol_doneSolution (rproblem%rrefSolution)
    
    ! Release FE hierarchies
    do i=1,rproblem%nfeHierarchies
      call fesph_releaseHierarchy(rproblem%RfeHierarchies(i)%rfeHierarchy)
    end do
    rproblem%nfeHierarchies = 0

    call mshh_releaseHierarchy (rproblem%rmeshHierarchy)

    call parlst_done (rproblem%rparlist)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_vectorcopy (rproblem,isource,idest)

!<description>
  ! Copies a vector to another, applies a type cast if necessary.
!</description>

!<input>
  ! Source vector
  integer, intent(in) :: isource
  
  ! Destination vector
  integer, intent(in) :: idest
!</input>

!<inputoutput>
  ! Problem structure
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_solution), pointer :: p_rsource, p_rdest

    ! Get the source and destination vector.
    p_rsource => rproblem%Rsolutions(isource)
    if (idest .gt. 0) then
      p_rdest => rproblem%Rsolutions(idest)
    else
      p_rdest => rproblem%rrefSolution
    end if
    
    ! Is the destination vector discontinuous Q1?
    if (p_rsource%cfespace .ne. p_rdest%cfespace) then
      call output_line ("Projecting solution.")
    
      select case (p_rdest%cfespace)
      case (1211)
        ! Apply a projection
        call sol_projectToDGQ1 (&
            p_rsource%p_Rvectors(p_rsource%ireflevel),&
            p_rdest%p_Rvectors(p_rdest%ireflevel),&
            rproblem%rmeshHierarchy%p_Rtriangulations(1))
        
      end select
    else
      call lsysbl_copyVector (p_rsource%p_Rvectors(p_rsource%ireflevel),&
          p_rdest%p_Rvectors(p_rdest%ireflevel))
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_vectorLinearComb (rproblem,isource1,isource2,idest,cx,cy)

!<description>
  ! Linear combination of two vectors into a third one.
  ! All vectors must be on the same level.
!</description>

!<input>
  ! Source vectors
  integer, intent(in) :: isource1
  integer, intent(in) :: isource2
  
  ! Multipliers
  real(DP), intent(in) :: cx,cy
  
  ! Destination vector. =0 targets into the reference solution
  integer, intent(in) :: idest
!</input>

!<inputoutput>
  ! Problem structure
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_solution), pointer :: p_rsource1, p_rsource2, p_rdest
    integer :: ilev

    ! Get the source and destination vector.
    p_rsource1 => rproblem%Rsolutions(isource1)
    p_rsource2 => rproblem%Rsolutions(isource2)
    if (idest .gt. 0) then
      p_rdest => rproblem%Rsolutions(idest)
    else
      p_rdest => rproblem%rrefSolution
    end if

    ! Destination level is the crucial one.
    ilev = p_rdest%ireflevel
    
    ! Linear combination
    call lsysbl_vectorLinearComb (p_rsource1%p_Rvectors(ilev),p_rsource2%p_Rvectors(ilev),cx,cy,&
        p_rdest%p_Rvectors(ilev))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_writeOutput (rproblem)

!<description>
  ! Creates an output file for the reference solution
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(in), target :: rproblem
!</input>

!</subroutine>

    ! local variables
    type(t_solution), pointer :: p_rsolution
    type(t_vectorBlock), pointer :: p_rvector
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: i
    character(len=SYS_STRLEN) :: sname

    ! Is the reference solution discontinuous Q1?

    ! Get the source and destination vector.
    p_rsolution => rproblem%rrefSolution
    p_rvector => p_rsolution%p_Rvectors(p_rsolution%ireflevel)
    
    ! Is the destination vector discontinuous Q1?
    select case (p_rsolution%cfespace)
    case (1211)
    
      select case (rproblem%cpostprocfile)
      case (1)
        call output_line ("Writing postprocessing file.")
  
        ! Start an output file.
        call ucd_startVTK (rexport,UCD_FLAG_STANDARD + UCD_FLAG_DISCONTINUOUS,&
            p_rvector%p_rblockDiscr%p_rtriangulation,rproblem%spostprocfilename)
            
        ! Write all components
        do i=1,p_rvector%nblocks
          call lsyssc_getbase_double (p_rvector%RvectorBlock(i),p_Ddata)
          sname = "u"//trim(sys_siL(i,3))
          call ucd_addVariableVertexBased (rexport, sname, &
              UCD_VAR_STANDARD, p_Ddata)
        end do
        
        ! Write the file.
        call ucd_write (rexport)
        call ucd_release (rexport)
      end select
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_computeErrors (rproblem)

!<description>
  ! Creates an output file for the reference solution
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(in), target :: rproblem
!</input>

!</subroutine>

    ! local variables
    real(DP) :: dintvalue
    integer :: icomp
    type(t_fev2Vectors) :: revalVectors
    type(t_scalarCubatureInfo) :: rcubatureInfo
    type(t_vectorBlock), pointer :: p_rsol
    
    ! Get the solution
    p_rsol => rproblem%rrefSolution%p_Rvectors(rproblem%rrefSolution%ireflevel)

    if (rproblem%cerrorL2 .ne. 0) then
      ! Compute the L2 norm
      do icomp = 1,p_rsol%nblocks
        
        ! Cubature structure
        call spdiscr_createDefCubStructure (p_rsol%RvectorBlock(icomp)%p_rspatialDiscr, &
            rcubatureInfo, CUB_GEN_AUTO_G5)
      
        call fev2_addVectorToEvalList(revalVectors,p_rsol%RvectorBlock(icomp),0)
        
        ! Compute the error
        call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
            bma_fcalc_L2norm,revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
        
        ! Print the output
        call output_line ("||solution0_"//trim(sys_siL(icomp,10))//"||_L2 = "//&
            trim(sys_sdE(sqrt(dintvalue),10)))
        
        ! Cleanup
        call fev2_releaseVectorList(revalVectors)
        call spdiscr_releaseCubStructure (rcubatureInfo)
        
      end do
    end if
    
    if (rproblem%cerrorH1 .ne. 0) then
      ! Compute the L2 norm
      do icomp = 1,p_rsol%nblocks
        
        ! Cubature structure
        call spdiscr_createDefCubStructure (p_rsol%RvectorBlock(icomp)%p_rspatialDiscr, &
            rcubatureInfo, CUB_GEN_AUTO_G5)
      
        call fev2_addVectorToEvalList(revalVectors,p_rsol%RvectorBlock(icomp),1)
        
        ! Compute the error
        call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
            bma_fcalc_H1norm,revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
        
        ! Print the output
        call output_line ("||solution0_"//trim(sys_siL(icomp,10))//"||_H1 = "//&
            trim(sys_sdE(sqrt(dintvalue),10)))
        
        ! Cleanup
        call fev2_releaseVectorList(revalVectors)
        call spdiscr_releaseCubStructure (rcubatureInfo)
        
      end do
    end if
    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine base_vectorDivide (rproblem,idest,isource)

!<description>
  ! DIvides vector idest by vector isource, pointwise.
!</description>

!<input>
  ! Source vector
  integer, intent(in) :: isource
  
  ! Destination vector
  integer, intent(in) :: idest
!</input>

!<inputoutput>
  ! Problem structure
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_solution), pointer :: p_rsource, p_rdest
    integer :: ilev,i
    real(DP), dimension(:), pointer :: p_Dsource, p_Ddest

    ! Get the source and destination vector.
    p_rsource => rproblem%Rsolutions(isource)
    if (idest .gt. 0) then
      p_rdest => rproblem%Rsolutions(idest)
    else
      p_rdest => rproblem%rrefSolution
    end if

    ! Destination level is the crucial one.
    ilev = p_rdest%ireflevel
    
    ! Get data
    call lsysbl_getbase_double (p_rsource%p_Rvectors(ilev),p_Dsource)
    call lsysbl_getbase_double (p_rdest%p_Rvectors(ilev),p_Ddest)
    
    ! Divide
    do i=1,size(p_Ddest)
    
      if (abs(p_Dsource(i)) .lt. 1E-8_DP) then
        if (p_Dsource(i) .ge. 0.0_DP) then
          p_Ddest(i) = 100000.0_DP
        else
          p_Ddest(i) = -100000.0_DP
        end if
      else
        p_Ddest(i) = p_Ddest(i) / p_Dsource(i)
      end if
    
    end do

  end subroutine

end module
