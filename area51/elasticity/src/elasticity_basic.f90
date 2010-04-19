!##############################################################################
!# ****************************************************************************
!# <name> elasticity_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!#   This module contains basic routines for solving elasticity problems.
!#   The type t_problem is defined here which contains all necessary information
!#   for the problem to be computed. Finally, some constants needed for the computation
!#   are defined.
!# </purpose>
!##############################################################################

module elasticity_basic

  use fsystem
  use genoutput
  use element
  use cubature
  use linearsystemblock
  use discretebc
  use triangulation
  use spatialdiscretisation
  use paramlist

  implicit none

  integer, parameter :: EQ_POISSON          = 1
  integer, parameter :: EQ_ELASTICITY       = 2

  integer, parameter :: SIMUL_REAL          = 1
  integer, parameter :: SIMUL_ANALYTICAL    = 2

  integer, parameter :: BC_NEUMANN          = 1
  integer, parameter :: BC_DIRICHLET        = 2

  integer, parameter :: SOLVER_DIRECT       = 1
  integer, parameter :: SOLVER_CG           = 2
  integer, parameter :: SOLVER_BICGSTAB     = 3
  integer, parameter :: SOLVER_MG           = 4
  integer, parameter :: SOLVER_CG_MG        = 5
  integer, parameter :: SOLVER_BICGSTAB_MG  = 6
  integer, parameter :: SOLVER_MG_CG        = 7
  integer, parameter :: SOLVER_MG_BICGSTAB  = 8
  integer, parameter :: SOLVER_BICGSTAB_MG_CG       =  9
  integer, parameter :: SOLVER_BICGSTAB_MG_BICGSTAB = 10

  integer, parameter :: SMOOTHER_NO         = 0
  integer, parameter :: SMOOTHER_JACOBI     = 1
  integer, parameter :: SMOOTHER_ILU        = 2


!<types>
!<typeblock description="type for storing all necessary information about the simulation">
  type t_problem
    !   grid file
    character(len=500) :: sgridFileTri
    character(len=500) :: sgridFilePrm

    ! number of boundaries (has to be set manually by the user who has to know
    ! the number of boundaries in the current grid)
    integer :: nboundaries = 1

    ! number of boundary segments (has to be set manually by the user who has to know
    ! the number segments in the current grid)
    integer, dimension(:), pointer :: NboundarySegments

    ! max. number boundary segments over all boundaries
    integer :: nmaxNumBoundSegments
  
    ! kind of equation (possible values: EQ_POISSON, EQ_ELASTICITY)
    integer :: cequation = EQ_ELASTICITY

    ! number of blocks (1 for Poisson equation, 2 for 2D elasticity equation)
    integer :: nblocks
  
    ! material parameters (Poisson ratio nu and shear modulus mu)
    real(DP) :: dnu     = 0.3_DP
    real(DP) :: dmu     = 0.5_DP
    real(DP) :: dlambda = 0.75_DP

    ! definition of boundary conditions (BC_NEUMANN or BC_DIRICHLET)
    ! (dimension  nblocks x max. number segments x nboundaries)
    integer, dimension(:,:,:), pointer :: Cbc

    ! type of simulation (possible values: SIMUL_REAL, SIMUL_ANALYTICAL)
    integer :: csimulation = SIMUL_REAL

    ! given surface forces on Neumann boundary condition segments
    ! (dimension nblocks x max. number segments x nboundaries)
    real(DP), dimension(:,:,:), pointer :: DforceSurface

    ! constant RHS values (only needed in case of SIMUL_REAL)
    real(DP) :: dforceVolumeX   = 0.0_DP
    real(DP) :: dforceVolumeY   = 0.0_DP

    ! function IDs per component (only needed in case of SIMUL_ANALYTICAL)
    integer, dimension(2) :: CfuncID = (/4, 52/)

    ! kind of element used (possible values: EL_Q1, EL_Q2)
    integer :: celement

    ! 1D and 2D cubature formulas (they are automatically chosen according to the
    ! selected finite element)
    integer :: ccubature1D, ccubature2D
  
    ! MAX & MIN level where we want to solve.
    integer :: ilevelMax, ilevelMin

    ! kind of solver (possible values: SOLVER_DIRECT,BICGSTAB_SOLVER,SOLVER_MG,SOLVER_CG)
    integer :: csolver = SOLVER_DIRECT

    ! flag whether MG is involved
    logical :: bmgInvolved = .FALSE.
    
    ! max. number of iterations
    integer :: niterations = 5000

    ! tolerance
    real(DP) :: dtolerance = 1e-8_DP

    ! kind of elementary smoother (possible values: SMOOTHER_JACOBI, SMOOTHER_ILU)
    integer :: celementaryPrec = SMOOTHER_JACOBI

    ! MG cycle (0=F-cycle, 1=V-cycle, 2=W-cycle)
    integer :: ccycle = 0

    ! number of smoothing steps
    integer :: nsmoothingSteps = 2

    ! damping parameter
    real(DP) :: ddamp = 0.7_DP

    ! show deformation in gmv(possible values: YES, NO
    integer :: cshowDeformation = NO

    ! number of points where to evaluate the finite element solution
    integer :: nevalPoints = 0
  
    ! number of provided reference solutions in the first
    ! min(nevalPoints, nrefSols) evaluation points
    integer :: nrefSols = 0

    ! given points where to evaluate the FE solution (dimension ndim x nevalPoints)
    real(DP), dimension(:,:), pointer :: DevalPoints

    ! given reference solutions in the evaluation points (dimension ndim x nevalPoints)
    real(DP), dimension(:,:), pointer :: DrefSols

    ! values and derivatives of the FE function in the evaluation points
    ! dimension nblocks x 3 x nevalPoints (3 since we need FUNC, DERX and DERY)
    real(DP), dimension(:,:,:), pointer :: Dvalues

  end type t_problem
!</typeblock>

!<typeblock description="Type block defining all information about one level">
  type t_level
  
    ! object for saving the triangulation
    type(t_triangulation) :: rtriangulation

    ! object specifying the discretisation
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! system matrix for the specific level
    type(t_matrixBlock) :: rmatrix

    ! variable describing the discrete boundary conditions
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>
!</types>

  ! the one and only instance of the t_problem type which is used in this and the main
  ! module to store all necessary information about the simulation
  type(t_problem), save :: rprob

contains


! ****************************************************************************************


!<subroutine>
  subroutine elast_readParameterFile(rprob, rparams)
  
!<description>
    ! get the following parameters from the parameter file:
    !
    !   - gridFilePRM
    !   - gridFileTRI
    !
    !   - numBoundarySegments(#boundaries)
    !
    !   - equation
    !       'Poisson' or 'elasticity'
    !
    !   - nu
    !   - mu
    !
    !   - bc[i](#segments x #blocks)
    !       'D' or 'N', for each boundary i, each segment and component
    !
    !   - simulation
    !       'real' or 'analytic'
    !
    !   - forceSurface[i](#segments x dim)
    !       in case of real simulation, for each boundary i and each segment the
    !       surface force in x- and y-direction
    !
    !   - forceVolumeX
    !       given vol. force in x-direction in case of real simulation
    !
    !   - forceVolumeY
    !       given vol. force in y-direction in case of real simulation
    !
    !   - funcID_u1
    !       ID of analytical function for u1 in case of analytical simulation
    !
    !   - funcID_u2
    !       ID of analytical function for u2 in case of analytical simulation
    !
    !   - element
    !       'Q1' or 'Q2'
    !
    !   - levelMin
    !   - levelMax
    !
    !   - solver
    !       'DIRECT', 'CG', 'BICGSTAB', 'MG', 'CG_MG', 'MG_CG' or 'MG_BICGSTAB'
    !
    !   - numIter
    !   - tolerance
    !
    !   - smoother
    !       'Jacobi' or 'ILU'
    !   - mgCycle
    !       'V', 'F' or 'W'
    !   - numSmoothingSteps
    !   - damping
    !
    !   - showDeformation
    !       'YES' or 'NO'
    !
    !   - evalPoints(2*numEvalPoints)
    !       x- and y-coordinate of points to evaluate
    !   - refSols(2*numEvalPoints)
    !       ref. solution values for u1 and u2 in eval. points
!</description>
    
!<input>
!</input>
  
!<output>
    ! general problem structure
    type(t_problem), intent(out) :: rprob
    ! parameter list read from the parameter file
    type (t_parlist), intent(out) :: rparams
!</output>
    
!</subroutine>

    character(len=SYS_STRLEN) :: snameDatFile
    integer :: i, j, k
    character(len=SYS_STRLEN) :: sstring

    ! initialise the parameter structure and read the DAT file
    call parlst_init(rparams)
 
    ! get the data file
    call sys_getcommandLineArg(1, snameDatFile, &
                               sdefault='./dat/elast_2d_disp_smallDeform_static.dat')
    call parlst_readfromfile(rparams, snameDatFile)
    call output_line('parsing dat-file '//trim(snameDatFile)//'...')
    
    ! PRM file
    call parlst_getvalue_string(rparams, '', 'gridFilePRM', sstring)
    read(sstring,*) rprob%sgridFilePRM
    call output_line('PRM file: '//trim(rprob%sgridFilePRM))
                                 
    ! TRI file
    call parlst_getvalue_string(rparams, '', 'gridFileTRI', sstring)
    read(sstring,*) rprob%sgridFileTRI
    call output_line('TRI file: '//trim(rprob%sgridFilePRM))
       
    ! get number of boundaries by inquiring the number of items of the
    ! parameter 'numBoundarySegments' 
    rprob%nboundaries = parlst_querysubstrings(rparams, '', 'numBoundarySegments')
    call output_line('number of boundaries: '//trim(sys_siL(rprob%nboundaries,3)))
                       
    ! number of boundary segments per boundary (has to be set manually by the user)
    allocate(rprob%NboundarySegments(rprob%nboundaries))
    do i = 1,rprob%nboundaries
      call parlst_getvalue_int(rparams, '', 'numBoundarySegments', &
                               rprob%NboundarySegments(i), iarrayindex = i)
      call output_line('number of segments in boundary '//trim(sys_siL(i,3))//': ' // &
                       trim(sys_siL(rprob%NboundarySegments(i),4)))
    end do

    ! detect max. number of segments over all boundaries
    rprob%nmaxNumBoundSegments = -1
    do i = 1,rprob%nboundaries
      if (rprob%NboundarySegments(i) .gt. rprob%nmaxNumBoundSegments) then
        rprob%nmaxNumBoundSegments = rprob%NboundarySegments(i)
      end if
    end do
    call output_line('max. number of segments: ' // &
                     trim(sys_siL(rprob%nmaxNumBoundSegments,3)))


!BRAL: Nutze die funktionalitaet der parlst_getvalue routinen, einen Default-Wert zu
!setzen, falls Parameter nicht vorhanden.

    ! kind of equation (possible values: POISSON, ELASTICITY)
    call parlst_getvalue_string(rparams, '', 'equation', sstring)
    if (trim(sstring) .eq. 'Poisson') then
      rprob%cequation = EQ_POISSON
      rprob%nblocks = 1
    else if(trim(sstring) .eq. 'elasticity') then
      rprob%cequation = EQ_ELASTICITY
      rprob%nblocks = 2
    else
      call output_line('invalid equation:' // trim(sstring), OU_CLASS_ERROR, &
                       OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
      call sys_halt()
    end if
    call output_line('equation: '//trim(sstring))
    
    ! material parameters (Poisson ratio nu and shear modulus mu)
    if (rprob%cequation .eq. EQ_ELASTICITY) then
      call parlst_getvalue_double(rparams, '', 'nu', rprob%dnu)
      if (rprob%dnu .le. 0.0_DP .or. rprob%dnu .ge. 0.5) then
        call output_line('invalid value for nu:' // trim(sys_sdL(rprob%dnu,8)), &
                         OU_CLASS_ERROR, OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
        call sys_halt()
      endif
      call parlst_getvalue_double(rparams, '', 'mu', rprob%dmu)
      if (rprob%dmu .le. 0.0_DP) then
        call output_line('invalid value for mu:' // trim(sys_sdL(rprob%dmu,8)), &
                         OU_CLASS_ERROR, OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
        call sys_halt()
      endif
      rprob%dlambda = 2.0_DP * rprob%dmu * rprob%dnu/(1 - 2.0_DP * rprob%dnu)
      call output_line('nu: '//trim(sys_sdL(rprob%dnu,6)))
      call output_line('mu: '//trim(sys_sdEL(rprob%dmu,6)))
      call output_line('lambda: '//trim(sys_sdEL(rprob%dlambda,6)))
    endif
                     
    ! boundary conditions ('D' Dirichlet, 'N' Neumann)
    allocate(rprob%Cbc(rprob%nblocks, rprob%nmaxNumBoundSegments, rprob%nboundaries))
    do i = 1, rprob%nboundaries
      do j = 1,rprob%NboundarySegments(i)
        do k = 1, rprob%nblocks
          call parlst_getvalue_string(rparams, '', 'bc'//trim(sys_siL(i,3)), sstring, &
                                      isubstring = 2*(j-1) + k)
          if (trim(sstring) .eq. "D") then
            rprob%Cbc(k,j,i) = BC_DIRICHLET 
          else if (trim(sstring) .eq. "N") then
            rprob%Cbc(k,j,i) = BC_NEUMANN 
          else
            call output_line('invalid boundary condition:' // trim(sstring) // &
                             ', currently only D (Dirichlet) and N (Neumann) supported!',&
                             OU_CLASS_ERROR, OU_MODE_STD, &
                             'elast_2d_disp_smallDeform_static')
            call sys_halt()
          endif
          call output_line('BC of comp. ' // trim(sys_siL(k,3)) // ' in segment ' // &
                           trim(sys_siL(j,3)) // ' of boundary ' // &
                           trim(sys_siL(i,3))//': '// trim(sstring))
        enddo
      end do
    end do

    ! type of simulation (possible values: REAL, ANALYTICAL)
    call parlst_getvalue_string(rparams, '', 'simulation', sstring)
    if(trim(sstring) .eq. 'analytic') then
      rprob%csimulation = SIMUL_ANALYTICAL
    else if(trim(sstring) .eq. 'real') then
      rprob%csimulation = SIMUL_REAL
    else
      call output_line('invalid simulation:' // trim(sstring), &
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
      call sys_halt()
    end if
    call output_line('simulation: '//trim(sstring))

    ! surface forces (i.e. Neumann BCs) for all segments on all boundaries
    ! (only needed in case of csimulation .eq. SIMUL_REAL)
    if (rprob%csimulation .eq. SIMUL_REAL) then
      allocate(rprob%DforceSurface(rprob%nblocks, rprob%nmaxNumBoundSegments, &
                                   rprob%nboundaries))
      do i = 1, rprob%nboundaries
        do j = 1,rprob%NboundarySegments(i)
          do k = 1,rprob%nblocks
            call parlst_getvalue_double(rparams, '', 'forceSurface'//trim(sys_siL(i,3)), &
                                        rprob%DforceSurface(k,j,i), &
                                        iarrayindex = 2*(j-1)+k)
          enddo
          call output_line('(x,y)-surface force in segment ' // trim(sys_siL(j,3)) // &
                           ' of boundary ' // trim(sys_siL(i,3))//': (' // &
                           trim(sys_sdL(rprob%DforceSurface(1,j,i),4)) // &
                           ', '//trim(sys_sdL(rprob%DforceSurface(2,j,i),4))//')')
        end do
      end do
    endif

    ! constant volume forces (only needed in case of csimulation .eq. SIMUL_REAL)
    if (rprob%csimulation .eq. SIMUL_REAL) then
      call parlst_getvalue_double(rparams, '', 'forceVolumeX', rprob%dforceVolumeX)
      call parlst_getvalue_double(rparams, '', 'forceVolumeY', rprob%dforceVolumeY)
      call output_line('volume force: ('//trim(sys_sdL(rprob%dforceVolumeX,4)) // &
                       ', '//trim(sys_sdL(rprob%dforceVolumeY,4))//')')
    endif
         
    ! function IDs (only needed in case of csimulation .eq. SIMUL_ANALYTICAL)
    if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      call parlst_getvalue_int(rparams, '', 'funcID_u1', rprob%CfuncID(1))
      call parlst_getvalue_int(rparams, '', 'funcID_u2', rprob%CfuncID(2))
      call output_line('function ID for u1: ' // trim(sys_siL(rprob%CfuncID(1),3)))
      call output_line('function ID for u2: ' // trim(sys_siL(rprob%CfuncID(2),3)))
    endif
         
    ! get element type and choose cubature formula accordingly
    call parlst_getvalue_string(rparams, '', 'element', sstring)
    if (trim(sstring) .eq. "Q1") then
      rprob%celement = EL_Q1
      rprob%ccubature1D = CUB_G2_1D      
      rprob%ccubature2D = CUB_G2X2
      call output_line('element Q1, cubature G2 / G2X2')
    else if (trim(sstring) .eq. "Q2") then
      rprob%celement = EL_Q2
      rprob%ccubature1D = CUB_G3_1D      
      rprob%ccubature2D = CUB_G3X3
      call output_line('element Q2, cubature G3 / G3X3')
    else
      call output_line('invalid element:' // trim(sstring) // &
                       ', currently only Q1 and Q2 supported!', &
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
      call sys_halt()
    endif
                              
    ! minimum and maximum level
    call parlst_getvalue_int(rparams, '', 'levelMin', rprob%ilevelMin)   
    call parlst_getvalue_int(rparams, '', 'levelMax', rprob%ilevelMax)
    if (rprob%ilevelMin .le. 0 .or. rprob%ilevelMax .le. 0 .or. &
        rprob%ilevelMin .gt. rprob%ilevelMax) then
      call output_line('invalid combination of min./max. grid level: ' // &
                       trim(sys_siL(rprob%ilevelMin,3)) // &
                       "/" // trim(sys_siL(rprob%ilevelMax,3)))
    endif
        

! BRAL: NEW!
!    ! name of the solver file
!    call parlst_getvalue_string(rparams, '', 'solverFile', sstring)
!    call output_line('solver file: '//trim(sstring))
!    call parlst_readfromfile(rparams, sstring)
!    call output_line('parsing solver file '//trim(sstring)//'...')


    ! type of solver (possible values: DIRECT, CG, BICGSTAB, MG, CG_MG, MG_CG, MG_BICGSTAB)
    call parlst_getvalue_string(rparams, '', 'solver', sstring)

    rprob%bmgInvolved = .FALSE.
    if(trim(sstring) .eq. 'DIRECT') then
      rprob%csolver = SOLVER_DIRECT
    else if(trim(sstring) .eq. 'CG') then
      rprob%csolver = SOLVER_CG
    else if(trim(sstring) .eq. 'BICGSTAB') then
      rprob%csolver = SOLVER_BICGSTAB
    else if(trim(sstring) .eq. 'MG') then
      rprob%csolver = SOLVER_MG
      rprob%bmgInvolved = .TRUE.
    else if(trim(sstring) .eq. 'CG_MG') then
      rprob%csolver = SOLVER_CG_MG
      rprob%bmgInvolved = .TRUE.
    else if(trim(sstring) .eq. 'BICGSTAB_MG') then
      rprob%csolver = SOLVER_BICGSTAB_MG
      rprob%bmgInvolved = .TRUE.
    else if(trim(sstring) .eq. 'MG_CG') then
      rprob%csolver = SOLVER_MG_CG
      rprob%bmgInvolved = .TRUE.
    else if(trim(sstring) .eq. 'MG_BICGSTAB') then
      rprob%csolver = SOLVER_MG_BICGSTAB
      rprob%bmgInvolved = .TRUE.
    else if(trim(sstring) .eq. 'BICGSTAB_MG_CG') then
      rprob%csolver = SOLVER_BICGSTAB_MG_CG
      rprob%bmgInvolved = .TRUE.
    else if(trim(sstring) .eq. 'BICGSTAB_MG_BICGSTAB') then
      rprob%csolver = SOLVER_BICGSTAB_MG_BICGSTAB
      rprob%bmgInvolved = .TRUE.
    else
      call output_line('invalid solver:' // trim(sstring), &
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
      call sys_halt()
    end if
    call output_line('solver: '//trim(sstring))

    ! when no multigrid solver is involved, then set min. level to max. level
    if (.not. rprob%bmgInvolved) then
      rprob%ilevelMin = rprob%ilevelMax
      call output_line('grid level: ' // trim(sys_siL(rprob%ilevelMax,3)))
    else
      call output_line('min./max. grid level: ' // trim(sys_siL(rprob%ilevelMin,3)) // &
                       "/" // trim(sys_siL(rprob%ilevelMax,3)))
    end if

    if (rprob%csolver .ne. SOLVER_DIRECT) then
      ! max number of iterations
      call parlst_getvalue_int(rparams, '', 'numIter', rprob%niterations)
      call output_line('max. number of iterations: '//trim(sys_siL(rprob%niterations,10)))
  
      ! tolerance
      call parlst_getvalue_double(rparams, '', 'tolerance', rprob%dtolerance)
      call output_line('rel. stopping criterion: ' // trim(sys_sdEL(rprob%dtolerance,4)))

      ! type of elementary preconditioner/smoother (possible values: NO, JACOBI, ILU)
      call parlst_getvalue_string(rparams, '', 'elementaryPrec', sstring)
      if(trim(sstring) .eq. 'NO') then
        rprob%celementaryPrec = SMOOTHER_NO
      else if(trim(sstring) .eq. 'JACOBI') then
        rprob%celementaryPrec = SMOOTHER_JACOBI
      else if(trim(sstring) .eq. 'ILU') then
        rprob%celementaryPrec = SMOOTHER_ILU
      else
        call output_line('invalid elementary precond./smoother type:' // trim(sstring), &
                         OU_CLASS_ERROR, OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
        call sys_halt()
      end if
      call output_line('elementary preconditioner/smoother: '//trim(sstring))

      ! additional parameters for multigrid
      if (rprob%bmgInvolved) then
        ! MG cycle (0=F-cycle, 1=V-cycle, 2=W-cycle)
        call parlst_getvalue_string(rparams, '', 'mgCycle', sstring)
        if(trim(sstring) .eq. 'V') then
          rprob%ccycle = 1
        else if(trim(sstring) .eq. 'F') then
          rprob%ccycle = 0
        else if(trim(sstring) .eq. 'W') then
          rprob%ccycle = 2
        else
          call output_line('invalid mgCycle:' // trim(sstring))
          call output_line('Choosing F-cycle!')
          rprob%ccycle = 0
        end if
        ! number of smoothing steps
        call parlst_getvalue_int(rparams, '', 'numSmoothingSteps', rprob%nsmoothingSteps)
        call output_line('MG cycle: '//trim(sstring) // ':' // &
                         trim(sys_siL(rprob%nsmoothingSteps,3)) // ':' // &
                         trim(sys_siL(rprob%nsmoothingSteps,3)))
        ! damping parameter
        call parlst_getvalue_double(rparams, '', 'damping', rprob%ddamp)
        call output_line('damping parameter:' // trim(sys_sdL(rprob%ddamp,2)))
      endif
    endif

    ! show deformation in gmv (possible values: YES, NO)
    call parlst_getvalue_string(rparams, '', 'showDeformation', sstring)
    if(trim(sstring) .eq. 'YES') then
      rprob%cshowDeformation = YES
    else if(trim(sstring) .eq. 'NO') then
      rprob%cshowDeformation = NO
    else
      call output_line('invalid value for showDeformation:' // trim(sstring), &
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_2d_disp_smallDeform_static')
      call sys_halt()
    end if
    call output_line('show deformation: '//trim(sstring))

    ! get number of evaluation points by inquiring the number of items of the
    ! parameter 'evalPoints'
    rprob%nevalPoints = parlst_querysubstrings(rparams, '', 'evalPoints')/2
    call output_line('number of evaluation points: '//trim(sys_siL(rprob%nevalPoints,3)))

    if (rprob%nevalPoints .gt. 0) then
      allocate(rprob%DevalPoints(2,rprob%nevalPoints))
      ! we need to store values for 2 blocks x 3 function value types (FUNC, DERX, DERY)
      ! in each eval point
      allocate(rprob%Dvalues(2, 3, rprob%nevalPoints))

      rprob%DevalPoints = 0.0_DP
      rprob%Dvalues = 0.0_DP
      do i = 1, rprob%nevalPoints
        call parlst_getvalue_double(rparams, '', 'evalPoints', &
                                    rprob%DevalPoints(1,i), iarrayindex = 2*i-1)
        call parlst_getvalue_double(rparams, '', 'evalPoints', &
                                    rprob%DevalPoints(2,i), iarrayindex = 2*i)
        call output_line('eval. point: ('// trim(sys_sdL(rprob%DevalPoints(1,i),4)) &
                         // ', ' // trim(sys_sdL(rprob%DevalPoints(2,i),4)) // ')')
      end do

      ! get number of reference solutions in evaluation points by inquiring the number of
      ! items of the parameter 'refSols'
      rprob%nrefSols = parlst_querysubstrings(rparams, '', 'refSols') / 2
      call output_line('number of reference solutions: '//trim(sys_siL(rprob%nrefSols,3)))
  
      if (rprob%nrefSols .gt. 0) then
        allocate(rprob%DrefSols(2,rprob%nrefSols))
        rprob%DrefSols = 0.0_DP
        do i = 1, rprob%nrefSols
          call parlst_getvalue_double(rparams, '', 'refSols', &
                                      rprob%DrefSols(1,i), iarrayindex = 2*i-1)
          call parlst_getvalue_double(rparams, '', 'refSols', &
                                      rprob%DrefSols(2,i), iarrayindex = 2*i)
          call output_line('ref. sol.: ('// trim(sys_sdL(rprob%DrefSols(1,i),8)) &
                           // ', ' // trim(sys_sdL(rprob%DrefSols(1,i),8)) // ')')
        end do
      end if
    else
      ! when there are no evaluation points, we also need no reference solutions
      rprob%nrefSols = 0
    end if

  end subroutine elast_readParameterFile

end module elasticity_basic

