!#########################################################################################
!# ***************************************************************************************
!# <name> elasticity_basic </name>
!# ***************************************************************************************
!#
!# <purpose>
!#   This module contains basic routines and structures needed for solving elasticity
!#   problems.
!# </purpose>
!#########################################################################################

module elasticity_basic

  use fsystem
  use boundary
  use element
  use cubature
  use linearsystemblock
  use discretebc
  use triangulation
  use spatialdiscretisation

  implicit none

  integer, parameter :: EQ_POISSON          = 1
  integer, parameter :: EQ_ELASTICITY       = 2

  integer, parameter :: FORMULATION_DISPL   = 1
  integer, parameter :: FORMULATION_MIXED   = 2
  integer, parameter :: FORMULATION_STOKES  = 3

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

    ! spatial dimension of the problem (2 for 2D, 3 for 3D)
    integer :: ndim = 2

    ! grid file
    character(len=500) :: sgridFileTri
    character(len=500) :: sgridFilePrm

    ! object for saving the domain
    type(t_boundary) :: rboundary

    ! kind of equation (possible values: EQ_POISSON, EQ_ELASTICITY)
    integer :: cequation = EQ_ELASTICITY

    ! finite element formulation
    ! (possible values: FORMULATION_DISPL, FORMULATION_MIXED, FORMULATION_STOKES)
    integer :: cformulation = FORMULATION_DISPL

    ! number of blocks (1 for Poisson equation, 2 for 2D elasticity equation in pure
    ! displacement formulation, 3 for mixed u/p formulation)
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

    ! given boundray condition values (surface forces for Neumann BC segments, prescribed
    ! displacements for Dirichlet BCs)
    ! (dimension nblocks x max. number segments x nboundaries)
    real(DP), dimension(:,:,:), pointer :: DbcValue

    ! constant RHS values (only needed in case of SIMUL_REAL)
    real(DP) :: dforceVolumeX   = 0.0_DP
    real(DP) :: dforceVolumeY   = 0.0_DP

    ! function IDs per component (only needed in case of SIMUL_ANALYTICAL)
    integer, dimension(3) :: CfuncID = (/0, 0, 0/)

    ! finite element (possible values: EL_Q1, EL_Q2)
    integer :: celement = EL_Q1

    ! finite element for the pressure space (possible values: EL_Q1, EL_Q2,
    ! EL_QP1, EL_QP1NP, EL_QP1NPD) (only necessary in case of the mixed formulation)
    integer :: celementPress = EL_Q1

    ! 1D and 2D cubature formulas (they are automatically chosen according to the
    ! selected finite element)
    integer :: ccubature1D = CUB_G2_1D, ccubature2D = CUB_G2X2

    ! 1D and 2D cubature formulas for the pressure discretisation (they are automatically
    ! chosen according to the selected finite element)
    ! (only necessary in case of the mixed formulation)
    integer :: ccubaturePress1D = CUB_G2_1D, ccubaturePress2D = CUB_G2X2
  
    ! MAX & MIN level where we want to solve.
    integer :: ilevelMax = 3, ilevelMin = 3

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
  
    use fsystem
    use paramlist
    use genoutput

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
    !   - elementPress
    !       'Q1', 'Q2', 'P1', 'P1_NP' or 'P1_NPD'
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
    !
    ! The routine also reads in the boundary. In doing so, the number of boundary
    ! components can be inquired automatically and the parameter file can be parsed
    ! accordingly.
!</description>
    
!<input>
!</input>
  
!<output>
    ! general problem structure
    type(t_problem), intent(out) :: rprob
    ! parameter list read from the parameter file
    type(t_parlist), intent(out) :: rparams
!</output>
    
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: snameDatFile
    integer :: i, j, k, iaux
    character(len=SYS_STRLEN) :: sstring
    character(len=1), dimension(3) :: Sbc
    real(DP), dimension(3) :: Dval
    ! max. number boundary segments over all boundaries
    integer :: nmaxNumBoundSegments
  
    ! initialise the parameter structure and read the DAT file
    call parlst_init(rparams)
 
    ! get the data file
    call sys_getcommandLineArg(1, snameDatFile, &
                               sdefault = './dat/elasticity_2d_smallDef_stat.dat')
    call parlst_readfromfile(rparams, snameDatFile)
    call output_line('parsing dat-file '//trim(snameDatFile)//'...')
    
    ! PRM file
    call parlst_getvalue_string(rparams, '', 'gridFilePRM', sstring)
    read(sstring,*) rprob%sgridFilePRM
    call output_line('PRM file: '//trim(rprob%sgridFilePRM))
                                 
    ! read in the parameterisation of the boundary and save it to rboundary
    ! this is done already here so that the number of boundary components can be
    ! inquired automatically
    call output_line('reading boundary parameterisation from file ' // &
                     trim(rprob%sgridFilePRM) // '...')
    call boundary_read_prm(rprob%rboundary, rprob%sgridFilePRM)

    ! TRI file
    call parlst_getvalue_string(rparams, '', 'gridFileTRI', sstring)
    read(sstring,*) rprob%sgridFileTRI
    call output_line('TRI file: '//trim(rprob%sgridFilePRM))
       
    ! get number of boundaries
    call output_line('number of boundaries: ' // &
                     trim(sys_siL(boundary_igetNBoundComp(rprob%rboundary),3)))

    ! inquire number of boundary segments per boundary and calculate max. number of
    ! segments over all boundaries
    nmaxNumBoundSegments = -1
    do i = 1,boundary_igetNBoundComp(rprob%rboundary)
      iaux = boundary_igetNsegments(rprob%rboundary, i)
      call output_line('number of segments in boundary '//trim(sys_siL(i,3))//': ' // &
                       trim(sys_siL(iaux,4)))
      if (iaux .gt. nmaxNumBoundSegments) then
        nmaxNumBoundSegments = iaux
      end if
    enddo
    call output_line('max. number of segments: ' // trim(sys_siL(nmaxNumBoundSegments,3)))

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
                       OU_MODE_STD, 'elast_readParameterFile')
      call sys_halt()
    end if
    call output_line('equation: '//trim(sstring))
    
    ! type of finite element formulation (pure displacement or mixed)
    if (rprob%cequation .eq. EQ_ELASTICITY) then
      call parlst_getvalue_string(rparams, '', 'formulation', sstring)
      if (trim(sstring) .eq. 'displ') then
        rprob%cformulation = FORMULATION_DISPL
        rprob%nblocks = 2
      else if(trim(sstring) .eq. 'mixed') then
        rprob%cformulation = FORMULATION_MIXED
        rprob%nblocks = 3
      else if(trim(sstring) .eq. 'Stokes') then
        rprob%cformulation = FORMULATION_STOKES
        rprob%nblocks = 3
      else
        call output_line('invalid FE formulation:' // trim(sstring), OU_CLASS_ERROR, &
                         OU_MODE_STD, 'elast_readParameterFile')
        call sys_halt()
      end if
      call output_line('FE formulation: '//trim(sstring))
    endif

    ! material parameters (Poisson ratio nu and shear modulus mu)
    if (rprob%cequation .eq. EQ_ELASTICITY) then
      call parlst_getvalue_double(rparams, '', 'nu', rprob%dnu)
      if (rprob%dnu .le. 0.0_DP .or. rprob%dnu .gt. 0.5) then
        call output_line('invalid value for nu:' // trim(sys_sdL(rprob%dnu,8)), &
                         OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
        call sys_halt()
      endif
      call parlst_getvalue_double(rparams, '', 'mu', rprob%dmu)
      if (rprob%dmu .le. 0.0_DP) then
        call output_line('invalid value for mu:' // trim(sys_sdL(rprob%dmu,8)), &
                         OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
        call sys_halt()
      endif
      call output_line('nu: ' // trim(sys_sdL(rprob%dnu,6)))
      call output_line('mu: ' // trim(sys_sdEL(rprob%dmu,6)))
      ! compute second Lame constant lambda
      if (rprob%dnu .eq. 0.5) then
        if (rprob%cformulation .eq. FORMULATION_DISPL) then
          call output_line('nu = 0.5 not feasible in pure displacement formulation', &
                           OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
          call sys_halt()
        endif
        rprob%dlambda = SYS_INFINITY_DP
        call output_line('lambda: infinity')
      else ! dnu .lt. 0.5
        if (rprob%cformulation .eq. FORMULATION_STOKES) then
          call output_line('nu .ne. 0.5 not feasible in Stokes formulation', &
                           OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
          call sys_halt()
        endif
        rprob%dlambda = 2.0_DP*rprob%dmu * rprob%dnu/(1.0_DP - 2.0_DP*rprob%dnu)
        call output_line('lambda: ' // trim(sys_sdEL(rprob%dlambda,6)))
      endif
    endif
                     
    ! type of simulation (possible values: REAL, ANALYTICAL)
    call parlst_getvalue_string(rparams, '', 'simulation', sstring)
    if(trim(sstring) .eq. 'analytic') then
      rprob%csimulation = SIMUL_ANALYTICAL
    else if(trim(sstring) .eq. 'real') then
      rprob%csimulation = SIMUL_REAL
    else
      call output_line('invalid simulation:' // trim(sstring), &
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
      call sys_halt()
    end if
    call output_line('simulation: '//trim(sstring))


    ! process boundary conditions
    ! Currently, the BCs can only be defined per boundary segment. For each boundary,
    ! the dat file must provide a number of lines corresponding to the segments.
    ! The first two entries determine the type of the boundary condition (('D' Dirichlet,
    ! 'N' Neumann), the last two entries the BC values for the two components. A nonzero
    ! Neumann value means a line force in the corresponding direction, while a zero
    ! Dirichlet value fixes the body in this direction. Zero Neumann means 'do nothing'.
    ! Nonzero Dirichlet (prescribed displacements) are not supported yet!
    ! Example for a boundary with three segments:
    !   'N' 'N' 0.0 0.0             # free in x- and y-direction
    !   'D' 'N' 0.0 -100.0          # fixed in x-direction, negative force in y-direction
    !   'D' 'D' 0.0 0.0             # fixed in x- and y-direction
    allocate(rprob%Cbc(rprob%nblocks, nmaxNumBoundSegments, &
                       boundary_igetNBoundComp(rprob%rboundary)))
    if (rprob%csimulation .eq. SIMUL_REAL) then
      allocate(rprob%DbcValue(rprob%nblocks, nmaxNumBoundSegments, &
                                   boundary_igetNBoundComp(rprob%rboundary)))
    endif
    do i = 1, boundary_igetNBoundComp(rprob%rboundary)
      call output_line("Boundary conditions on boundary " // trim(sys_siL(i,3))//":")
      do j = 1,boundary_igetNsegments(rprob%rboundary, i)
        Dval = 0.0_DP
        call parlst_getvalue_string(rparams, '', 'bc'//trim(sys_siL(i,3)), sstring, &
                                    isubstring = j)
        if (rprob%csimulation .eq. SIMUL_REAL) then
          ! in case of a real simulation read BC types and values
          if (rprob%cformulation .eq. FORMULATION_DISPL) then
            read(sstring,*) Sbc(1), Sbc(2), Dval(1), Dval(2)
          else
            read(sstring,*) Sbc(1), Sbc(2), Sbc(3), Dval(1), Dval(2), Dval(3)
          endif
        else
          ! in case of an analytical simulation read only BC types
          if (rprob%cformulation .eq. FORMULATION_DISPL) then
            read(sstring,*) Sbc(1), Sbc(2)
          else
            read(sstring,*) Sbc(1), Sbc(2), Sbc(3)
          endif
        endif

        ! set type of boundary condition
        do k = 1, rprob%nblocks
          if (trim(Sbc(k)) .eq. "D") then
            rprob%Cbc(k,j,i) = BC_DIRICHLET
          else if (trim(Sbc(k)) .eq. "N") then
            rprob%Cbc(k,j,i) = BC_NEUMANN
          else
            call output_line('invalid boundary condition:' // trim(sstring) // &
                             ', currently only D (Dirichlet) and N (Neumann) supported!',&
                             OU_CLASS_ERROR, OU_MODE_STD, &
                             'elast_readParameterFile')
            call sys_halt()
          endif
          sstring = "  segment " // trim(sys_siL(j,1)) // ", comp " // &
                    trim(sys_siL(k,3)) // ": "// trim(Sbc(k))

          ! set BC values (only needed in case of csimulation .eq. SIMUL_REAL)
          if (rprob%csimulation .eq. SIMUL_REAL) then
            rprob%DbcValue(k,j,i) = Dval(k)
            if (rprob%Cbc(k,j,i) .eq. BC_DIRICHLET) then
              if (k .le. rprob%ndim) then
                sstring = trim(sstring) // "  (displacement: " // &
                          trim(sys_sdL(rprob%DbcValue(k,j,i),4)) // ")"
              else
                sstring = trim(sstring) // "  (    pressure: " // &
                          trim(sys_sdL(rprob%DbcValue(k,j,i),4)) // ")"
              endif
              if (Dval(k) .ne. 0.0_DP) then
                call output_line('invalid boundary condition value:' // &
                                 trim(sys_sdEL(Dval(k),2)) // &
                                 ', currently only zero Dirichlet values supported!',&
                                 OU_CLASS_ERROR, OU_MODE_STD, &
                                 'elast_readParameterFile')
              endif
            else
              if (k .le. rprob%ndim) then
                sstring = trim(sstring) // "  (  line force: " // &
                          trim(sys_sdL(rprob%DbcValue(k,j,i),4)) // ")"
              endif
            endif
          endif
          call output_line(trim(sstring))
        enddo
      end do
    end do

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
      if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
          rprob%cformulation .eq. FORMULATION_STOKES) then
        call parlst_getvalue_int(rparams, '', 'funcID_p', rprob%CfuncID(3))
        call output_line('function ID for p: ' // trim(sys_siL(rprob%CfuncID(3),3)))
      endif
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
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
      call sys_halt()
    endif
    if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
        rprob%cformulation .eq. FORMULATION_STOKES) then
      call parlst_getvalue_string(rparams, '', 'elementPress', sstring)
      if (trim(sstring) .eq. "Q0") then
        rprob%celementPress = EL_Q0
! BRAL: sinnvolle Kubatur-Formel?
        rprob%ccubaturePress1D = CUB_G2_1D
        rprob%ccubaturePress2D = CUB_G2X2
        call output_line('pressure element Q0, cubature G2 / G2X2')
      else if (trim(sstring) .eq. "Q1") then
        rprob%celementPress = EL_Q1
        rprob%ccubaturePress1D = CUB_G2_1D
        rprob%ccubaturePress2D = CUB_G2X2
        call output_line('pressure element Q1, cubature G2 / G2X2')
      else if (trim(sstring) .eq. "Q2") then
        rprob%celementPress = EL_Q2
        rprob%ccubaturePress1D = CUB_G3_1D
        rprob%ccubaturePress2D = CUB_G3X3
        call output_line('pressure element Q2, cubature G3 / G3X3')
      else if (trim(sstring) .eq. "P1") then
        ! QP1-element
        rprob%celementPress =  EL_QP1
        rprob%ccubaturePress1D = CUB_G2_1D
        rprob%ccubaturePress2D = CUB_G2X2
        call output_line('pressure element P1, cubature G2 / G2X2')
      else if (trim(sstring) .eq. "P1_NP") then
        ! QP1-element, nonparametric
        rprob%celementPress =  EL_QP1NP
        rprob%ccubaturePress1D = CUB_G2_1D
        rprob%ccubaturePress2D = CUB_G2X2
        call output_line('pressure element P1 non parametric, cubature G2 / G2X2')
      else if (trim(sstring) .eq. "P1_NPD") then
        ! QP1-element, nonparametric, direct on element
        rprob%celementPress =  EL_QP1NPD
        rprob%ccubaturePress1D = CUB_G2_1D
        rprob%ccubaturePress2D = CUB_G2X2
        call output_line('pressure element P1 non parametric direct, cubature G2 / G2X2')
      else
        call output_line('invalid pressure element:' // trim(sstring) // &
                         ', currently only Q0, Q1, Q2, P1, P1_NP and P1_NPD supported!', &
                         OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
        call sys_halt()
      endif
    endif
    
    ! minimum and maximum level
    call parlst_getvalue_int(rparams, '', 'levelMin', rprob%ilevelMin)
    call parlst_getvalue_int(rparams, '', 'levelMax', rprob%ilevelMax)
    if (rprob%ilevelMax .le. 0) then
      call output_line('invalid max. grid level: ' // trim(sys_siL(rprob%ilevelMax,3)))
      call sys_halt()
    endif
    if (rprob%ilevelMin .le. 0 .or. rprob%ilevelMin .gt. rprob%ilevelMax) then
      ! every "invalid" min. level is replaced by max level (so the user can set
      ! "levelMin = 0" or "levelMin = 666" in his dat file in order to always perform a
      ! one-level computation on level levelMax)
      call output_line('replacing (invalid) min. level ' // &
                       trim(sys_siL(rprob%ilevelMin,3)) // " by max. level " // &
                       trim(sys_siL(rprob%ilevelMax,3)))
    endif
    
    ! name of the solver file
    call parlst_getvalue_string(rparams, '', 'solverFile', sstring)
    call output_line('solver file: '//trim(sstring))
    call parlst_readfromfile(rparams, sstring)
    call output_line('parsing solver file '//trim(sstring)//'...')

    call output_line('min./max. grid level: ' // trim(sys_siL(rprob%ilevelMin,3)) // &
                     "/" // trim(sys_siL(rprob%ilevelMax,3)))

    ! show deformation in gmv (possible values: YES, NO)
    call parlst_getvalue_string(rparams, '', 'showDeformation', sstring)
    if(trim(sstring) .eq. 'YES') then
      rprob%cshowDeformation = YES
    else if(trim(sstring) .eq. 'NO') then
      rprob%cshowDeformation = NO
    else
      call output_line('invalid value for showDeformation:' // trim(sstring), &
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_readParameterFile')
      call sys_halt()
    end if
    call output_line('show deformation: '//trim(sstring))

    ! get number of evaluation points by inquiring the number of items of the
    ! parameter 'evalPoints'
    rprob%nevalPoints = parlst_querysubstrings(rparams, '', 'evalPoints')
    call output_line('number of evaluation points: '//trim(sys_siL(rprob%nevalPoints,3)))

    if (rprob%nevalPoints .gt. 0) then
      allocate(rprob%DevalPoints(2,rprob%nevalPoints))
      ! we need to store values for nblocks blocks x 3 function value types
      ! (FUNC, DERX, DERY) in each eval point
      allocate(rprob%Dvalues(rprob%nblocks, 3, rprob%nevalPoints))

      rprob%DevalPoints = 0.0_DP
      rprob%Dvalues = 0.0_DP
      do i = 1, rprob%nevalPoints
        call parlst_getvalue_string(rparams, '', 'evalPoints', sstring, isubstring = i)
        read(sstring,*) rprob%DevalPoints(1,i), rprob%DevalPoints(2,i)
        call output_line('eval. point: ('// trim(sys_sdL(rprob%DevalPoints(1,i),4)) &
                         // ', ' // trim(sys_sdL(rprob%DevalPoints(2,i),4)) // ')')
      end do

      ! get number of reference solutions in evaluation points by inquiring the number of
      ! items of the parameter 'refSols'
      rprob%nrefSols = parlst_querysubstrings(rparams, '', 'refSols')
      call output_line('number of reference solutions: '//trim(sys_siL(rprob%nrefSols,3)))
  
      if (rprob%nrefSols .gt. 0) then
        allocate(rprob%DrefSols(rprob%nblocks,rprob%nrefSols))
        rprob%DrefSols = 0.0_DP
        do i = 1, rprob%nrefSols
          call parlst_getvalue_string(rparams, '', 'refSols', sstring, isubstring = i)
          if (rprob%cformulation .eq. FORMULATION_DISPL) then
            read(sstring,*) rprob%DrefSols(1,i), rprob%DrefSols(2,i)
            call output_line('ref. sol.: ('// trim(sys_sdL(rprob%DrefSols(1,i),8)) &
                             // ', ' // trim(sys_sdL(rprob%DrefSols(2,i),8)) // ')')
          else
            read(sstring,*) rprob%DrefSols(1,i), rprob%DrefSols(2,i), rprob%DrefSols(3,i)
            call output_line('ref. sol.: ('// trim(sys_sdL(rprob%DrefSols(1,i),8)) // &
                             ', ' // trim(sys_sdL(rprob%DrefSols(2,i),8)) // &
                             ', ' // trim(sys_sdL(rprob%DrefSols(3,i),8)) // ')')
          endif
        end do
      end if
    else
      ! when there are no evaluation points, we also need no reference solutions
      rprob%nrefSols = 0
    end if

  end subroutine elast_readParameterFile


! ****************************************************************************************

  
!<subroutine>
  subroutine elast_calcErrors(rprob, rsol)
  
    use fsystem
    use genoutput
    use derivatives
    use pprocerror
    use feevaluation
    use collection

!<description>
  ! This routine calculates some values and errors. In given evaluation points, it
  ! computes FE solutions, derivatives, absolute error, strains and stresses. In case,
  ! a analytical reference solution is given, it computes the error between the FE and
  ! this reference solution.
!</description>
    
!<input>
    ! general problem structure
    type(t_problem), intent(inout) :: rprob
  
    ! solution vector
    type(t_vectorBlock), intent(in) :: rsol
!</input>
  
!<output>
!</output>
!</subroutine>

    ! error between FE function and reference function
    type(t_errorScVec) :: rerror
    real(DP), dimension(:), pointer :: DerrorL2
    real(DP), dimension(:), pointer :: DerrorH1

    ! which function value types to calculate in evaluation points (required by the
    ! subroutine fevl_evaluate2(...))
    integer, dimension(3) :: CderType = (/DER_FUNC, DER_DERIV_X, DER_DERIV_Y/)

    ! auxiliary variables
    real(DP) ::  ddivu, deps11, deps22, deps12, daux1, daux2
    real(DP) ::  dsigma11, dsigma12, dsigma22, dsigma33, dtrace, dmises
    integer :: i

    ! collection structure to provide additional information
    type(t_collection) :: rcollection

    ! calculate in given evaluation points: FE solutions, derivatives, absolute error,
    ! strains, stresses
    if (rprob%nevalPoints .gt. 0) then

      ! Call the appropriate evaluation routine, that calculates all required function
      ! value types for all FE components in all evaluation points in one sweep. The
      ! FE solution values are stored in rprob%Dvalues(:,:,:) with dimension
      ! nblocks x 3 x nevalPoints (3 since we need FUNC, DERX and DERY)
      if (rprob%cformulation .ne. FORMULATION_DISPL .and. &
          rprob%celement .ne. rprob%celementPress) then
        ! In case of the mixed formulation or Stokes with different discretisations for
        ! u and p, there have to be two calls of the function.
        call fevl_evaluate(CderType, rprob%Dvalues, rsol, rprob%DevalPoints, &
                           iblockMin = 1, iblockMax = rprob%ndim)
        call fevl_evaluate(CderType, rprob%Dvalues, rsol, rprob%DevalPoints, &
                           iblockMin = rprob%nblocks, iblockMax = rprob%nblocks)
      else
        ! otherwise, the discretisation is the same for all components
        call fevl_evaluate(CderType, rprob%Dvalues, rsol, rprob%DevalPoints)
      endif

      call output_line('Values in evaluation points:')
      do i = 1, rprob%nevalPoints
        call output_line('   point: (' // trim(sys_sdL(rprob%DevalPoints(1,i),4)) // &
                         ', ' // trim(sys_sdL(rprob%DevalPoints(2,i),4)) // ')')
        call output_line('     u1h: ' // trim(sys_sdEL(rprob%Dvalues(1,1,i),10)))
        call output_line('     u2h: ' // trim(sys_sdEL(rprob%Dvalues(2,1,i),10)))
        if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
            rprob%cformulation .eq. FORMULATION_STOKES) then
          call output_line('      ph: ' // trim(sys_sdEL(rprob%Dvalues(3,1,i),10)))
        endif
        deps11 = rprob%Dvalues(1,2,i)
        deps22 = rprob%Dvalues(2,3,i)
        deps12 = 0.5_DP*(rprob%Dvalues(2,2,i) + rprob%Dvalues(1,3,i))
        call output_line('   eps11: ' // trim(sys_sdEL(deps11,10)))
        call output_line('   eps22: ' // trim(sys_sdEL(deps22,10)))
        call output_line('   eps12: ' // trim(sys_sdEL(deps12,10)))
        ! divergence of u
        ddivu = (rprob%Dvalues(1,2,i) + rprob%Dvalues(2,3,i))
        call output_line('  div(u): ' // trim(sys_sdEL(ddivu,10)))
        ! von Mises stresses
        dsigma11 = 2.0_DP*rprob%dmu*deps11 + rprob%dlambda*ddivu
        dsigma22 = 2.0_DP*rprob%dmu*deps22 + rprob%dlambda*ddivu
        dsigma12 = rprob%dmu*deps12
        dsigma33 = rprob%dlambda*ddivu
        ! trace of the stress tensor divided by 3
        dtrace = (dsigma11 + dsigma22 + dsigma33)/3.0_DP
        dmises = sqrt(  (dsigma11 - dtrace)**2 + (dsigma22 - dtrace)**2 &
                      + (dsigma33 - dtrace)**2 + 2.0_DP*dsigma12**2)
        call output_line('   sig11: ' // sys_sdEL(dsigma11,10) )
        call output_line('   sig22: ' // sys_sdEL(dsigma22,10) )
        call output_line('   sig12: ' // sys_sdEL(dsigma12,10) )
        call output_line('   sig33: ' // sys_sdEL(dsigma33,10) )
        call output_line('   mises: ' // sys_sdEL(dmises, 10) )
        call output_lbrk()
      enddo

      call output_line('Relative errors between FE and reference solutions:')
      do i = 1, min(rprob%nevalPoints, rprob%nrefSols)
        call output_line('   point: (' // trim(sys_sdL(rprob%DevalPoints(1,i),4)) // &
                         ', ' // trim(sys_sdL(rprob%DevalPoints(2,i),4)) // ')')
        call output_line('     u1h: ' // trim(sys_sdEL(rprob%Dvalues(1,1,i),10)))
        call output_line('     u1*: ' // trim(sys_sdEL(rprob%DrefSols(1,i),10)))
        call output_line('     u2h: ' // trim(sys_sdEL(rprob%Dvalues(2,1,i),10)))
        call output_line('     u2*: ' // trim(sys_sdEL(rprob%DrefSols(2,i),10)))
        if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
            rprob%cformulation .eq. FORMULATION_STOKES) then
          call output_line('      ph: ' // trim(sys_sdEL(rprob%Dvalues(3,1,i),10)))
          call output_line('      p*: ' // trim(sys_sdEL(rprob%DrefSols(3,i),10)))
        endif

        daux1 = rprob%DrefSols(1,i)
        daux2 = rprob%DrefSols(1,i) - rprob%Dvalues(1,1,i)
        if (daux1 .ne. 0.0_DP) then
          daux2 = daux2/daux1
        endif
        call output_line('error u1: ' // trim(sys_sdEL(abs(daux2), 10)))

        daux1 = rprob%DrefSols(2,i)
        daux2 = rprob%DrefSols(2,i) - rprob%Dvalues(2,1,i)
        if (daux1 .ne. 0.0_DP) then
          daux2 = daux2/daux1
        endif
        call output_line('error u2: ' // trim(sys_sdEL(abs(daux2), 10)))

        daux1 = sqrt(rprob%DrefSols(1,i)**2 + rprob%DrefSols(2,i)**2)
        daux2 = sqrt(  (rprob%DrefSols(1,i)-rprob%Dvalues(1,1,i))**2 &
                     + (rprob%DrefSols(2,i)-rprob%Dvalues(2,1,i))**2 )
        if (daux1 .ne. 0.0_DP) then
          daux2 = daux2/daux1
        endif
        call output_line(' error u: ' // trim(sys_sdEL(daux2, 10)))
        if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
            rprob%cformulation .eq. FORMULATION_STOKES) then
          daux1 = rprob%DrefSols(3,i)
          daux2 = rprob%DrefSols(3,i) - rprob%Dvalues(3,1,i)
          if (daux1 .ne. 0.0_DP) then
            daux2 = daux2/daux1
          endif
          call output_line(' error p: ' // trim(sys_sdEL(abs(daux2), 10)))
        endif
        call output_lbrk()
      enddo
    end if

    ! calculate the error between FE and reference solution
    if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      allocate(DerrorL2(rprob%nblocks), DerrorH1(rprob%nblocks))

      ! set pointers
      if (rprob%cformulation .ne. FORMULATION_DISPL .and. &
          rprob%celement .ne. rprob%celementPress) then
        ! In case of the mixed formulation or Stokes with different discretisations for
        ! u and p, there have to be two calls of the error routine pperr_scalarVec():
        ! first one for u, second one for p. (This also requires that one has to pass
        ! some extra information to elast_analFunc() via the collection structure, see
        ! below.)
        rerror%p_RvecCoeff => rsol%RvectorBlock(1:rprob%ndim)
        rerror%p_DerrorL2 => DerrorL2(1:rprob%ndim)
        rerror%p_DerrorH1 => DerrorH1(1:rprob%ndim)
      else
        ! otherwise, the discretisation is the same for all components
        rerror%p_RvecCoeff => rsol%RvectorBlock(1:rprob%nblocks)
        rerror%p_DerrorL2 => DerrorL2(1:rprob%nblocks)
        rerror%p_DerrorH1 => DerrorH1(1:rprob%nblocks)
      endif

      ! call the corresponding error routine
      ! Set rcollection%IquickAccess(4) = 0 in order to tell the routine elast_analFunc()
      ! that the displacement components are to be treated. (This is necessary due to a
      ! deficiency of the routine pperr_scalarVec().)
      rcollection%IquickAccess(4) = 0
      call pperr_scalarVec(rerror, elast_analFunc, rcollection)

      ! print the displacement errors
      call output_line('L2 error for u1: ' // sys_sdEL(DerrorL2(1),10) )
      call output_line('L2 error for u2: ' // sys_sdEL(DerrorL2(2),10) )
      call output_line('L2 error for  u: ' // &
                       sys_sdEL(sqrt(DerrorL2(1)**2 + DerrorL2(2)**2),10) )
      call output_line('H1 error for u1: ' // sys_sdEL(DerrorH1(1),10) )
      call output_line('H1 error for u2: ' // sys_sdEL(DerrorH1(2),10) )
      call output_line('H1 error for  u: ' // &
                       sys_sdEL(sqrt(DerrorH1(1)**2 + DerrorH1(2)**2),10) )

      ! Now treat the pressure errors.
      if (rprob%cformulation .ne. FORMULATION_DISPL) then
        if (rprob%celement .ne. rprob%celementPress) then
          ! If different discretisations are used for u and p, the pressure error has to
          ! be computed now.
          rerror%p_RvecCoeff => rsol%RvectorBlock(rprob%nblocks:rprob%nblocks)
          rerror%p_DerrorL2 => DerrorL2(rprob%nblocks:rprob%nblocks)
          rerror%p_DerrorH1 => DerrorH1(rprob%nblocks:rprob%nblocks)

          ! Set rcollection%IquickAccess(4) = 1 in order to tell the routine
          ! elast_analFunc() that the pressure component is to be treated. (This is
          ! necessary due to a deficiency of the routine pperr_scalarVec().)
          rcollection%IquickAccess(4) = 1
          call pperr_scalarVec(rerror, elast_analFunc, rcollection)
        endif
        ! output the pressure error
        call output_line('L2 error for  p: ' // sys_sdEL(DerrorL2(rprob%nblocks),10) )
        call output_line('H1 error for  p: ' // sys_sdEL(DerrorH1(rprob%nblocks),10) )
      endif

      deallocate(DerrorL2, DerrorH1)
    end if

  end subroutine elast_calcErrors


! ****************************************************************************************


!<subroutine>
  subroutine elast_analFunc(icomp, cderivative, rdiscretisation, nel, nptsPerEl, &
                            Dpoints, rdomainIntSubset, Dvalues, rcollection)
  
    use fsystem
    use spatialdiscretisation
    use derivatives
    use collection
    use domainintegration
    
!<description>
    ! This subroutine computes the (analytical) values of a function in a couple of points
    ! on a couple of elements.
    ! Implements interface 'intf_refFunctionScVec.inc'
!</description>
  
!<input>
    ! component for which the analytical function is to tbe computed
    integer, intent(in) :: icomp
  
    ! DER_xxxx derivative identifier that specifies what to compute: DER_FUNC=function
    ! value, DER_DERIV_X=x-derivative,...
    integer, intent(in) :: cderivative
  
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! number of elements, where the values are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element, where the values are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all the elements where the values are to be computed
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>
  
!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! array for storing the values of the (analytical) function (or derivative) in all
    ! points specified in Dpoints, DIMENSION(nptsPerEl, nel)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! If rcollection%IquickAccess(4) .eq. 0 this means that the displacement components
    ! are to be treated, otherwise the pressure component.
    if (rcollection%IquickAccess(4) .eq. 0) then
      ! compute solution for x- or y-displacement (u1 or u2)

      ! this strange '1:nel' is necessary here, since in some situations the
      ! dimensions of Davlues do *not* coincide with nptsPerEl x nel!
!BRAL: ueberpruefen, ob das nicht zu teuer ist (interne Kopie?)
!BRAL: elast_danalyticFunction() vielleicht besser als subroutine aufziehen...
      Dvalues(:,1:nel) = elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                                 cderivative, rprob%CfuncID(icomp))
    else
      ! compute pressure solution

      if (rprob%dnu .eq. 0.5) then
        ! incompressible case
        Dvalues(:,1:nel) = elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                                   cderivative, rprob%CfuncID(3))
      else
        ! in the compressible case the pressure is related to the displacements
        ! via p = - dlambda * div(u)
        if (cderivative .eq. DER_FUNC) then
          ! p = -dlambda * (u1_x + u2_y)
          Dvalues(:,1:nel) = -rprob%dlambda * &
            (elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                     DER_DERIV_X, rprob%CfuncID(1)) &
           + elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                     DER_DERIV_Y, rprob%CfuncID(2)))
        else if (cderivative .eq. DER_DERIV_X) then
          ! p_x = -dlambda * (u1_xx + u2_xy)
          Dvalues(:,1:nel) = -rprob%dlambda * &
            (elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                     DER_DERIV_XX, rprob%CfuncID(1)) &
           + elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                     DER_DERIV_XY, rprob%CfuncID(2)))
        else if (cderivative .eq. DER_DERIV_Y) then
          ! p_x = -dlambda * (u1_xy + u2_yy)
          Dvalues(:,1:nel) = -rprob%dlambda * &
            (elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                     DER_DERIV_XY, rprob%CfuncID(1)) &
           + elast_danalyticFunction(Dpoints(:,:,1:nel), nel, nptsPerEl, &
                                     DER_DERIV_YY, rprob%CfuncID(2)))
        endif
      endif
    endif

  end subroutine elast_analFunc


! ****************************************************************************************


!<function>
  function elast_danalyticFunction(Dpts, nel, nptsPerEl, cderiv, cselect, dparam) &
                                   result(Dval)

    use fsystem
    use derivatives

!<description>
    ! This function provides some analytic functions, which can be used for
    ! validating your FE code.
!</description>

!<input>
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpts

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nel
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: nptsPerEl
  
    ! derivative of the function to be calculated
    integer(I32), intent(in) :: cderiv
  
    ! selector for the desired function
    integer(I32), intent(in) :: cselect
  
    ! optional parameter to influence the solution
    real(DP), intent(in), optional :: dparam
!</input>

!<!--
    !
    ! cselect   u(x,y)
    !   0        0.0
    !   1        0.1 * x
    !   2        0.1 * x^2
    !   3        0.1 * y
    !   4        0.1 * y^2
    !   5        4 * x * (1 - x)
    !   6        4 * y * (1 - y)
    !   7        x * (1 - x) * y * (1 - y)    (zero on the boundary of the unit square;
    !                                          used in FBENCHMARK)
    !   8        -(x * x + y * y) * x         (used in FBENCHMARK)
    !   9        y^2 * (1 - y)^2 * x
    !  10        y^2 * (1 - y)^2 * (x - 1)
    !  11        0.25 *  (1/sqrt(2) + x - y) * (1/sqrt(2) - x + y)
    !                 * ((1-sqrt(2))/sqrt(2) + x + y) * ((1+sqrt(2))/sqrt(2) - x - y)
    !                                         (zero on the boundary of the unit square
    !                                          rotated by Pi/4)
    !  12        sin(x) * sin(y)
    !  13        0.05 * sin(4 * PI * x) * sin(4 * PI * y)
    !            (zero on the boundary of the unit square)
    !  14        cos(x) * cos(y)
    !  15        cos(PI/2 * (x + y))     (zero divergence)
    !  16        -cos(PI/2 * (x + y))    (zero divergence)
    !  17        2 * cos(x) * sin(y) - 2 * (1 - cos(1)) * sin(1)
    !                        (has zero integral over the boundary of the unit square)
    !  18        8 * (1 - x) (Stokes: pressure for parabolic inflow (6) on unit square)
    !  19        sin(PI/2 * (x - y))  (Stokes: to use as pressure with (15)/(16))
    !  20        -(x * x + y * y) * x * x (slight modification of (8); used in Poisson app.)
    !  21        cos(Pi/2)*x - sin(Pi/2)*y - x ((21) and (22) simulate a rigid body
    !  22        sin(Pi/2)*x + cos(Pi/2)*y - y  rotation of Pi/2 around the origin)
    !  23        1.0
    !  24        -(2x - 1)(2y^3 - 3y^2) / 6    (zero divergence together with 6.0_DP x (7))
    !  25        sin(x) cos(y)
    !  26        -cos(x) sin(y)
    !  27        4 - 8x          (Stokes, unit square, zero mean pressure: with /6/ + /0/)
    !  28        2 cos(x) sin(y) - 2 sin(1) + 2 sin(1) cos(1)
    !  29        xy - 1/4        (Stokes, unit square, zero mean pressure: with /7/ + /24/)
    !
    !  Function triple taken from Bochev/ Gunzburger/ Lehoucq, On stabilized finite
    !  element methods for the Stokes problem in the small time limit (preprint)
    !  30        sin(PI*x - 7/10) * sin(PI*y + 1/5)
    !  31        cos(PI*x - 7/10) * cos(PI*y + 1/5)
    !  32        sin(x) * cos(y) + (cos(1) - 1) * sin(1)
    !
    !  33        -sin(gamma x) * sin(gamma y)
    !  34        -cos(gamma x) * cos(gamma y)
    !            (gamma should be set to a multiple of pi. gamma=pi is just fine.
    !             A setting of gamma=3*Pi gives in combination with stabilised Q1/Q1
    !             approaches serious problems for the pressure solution of coarse grids.)
    !
    !  Function triple taken from Bochev/ Gunzburger/ Lehoucq, On stabilized finite
    !  element methods for the Stokes problem in the small time limit (journal version)
    !  (together with 32 for p). The difference to the preprint version (30,31) is that
    !  the velocities are not only divergence free but also zero on the boundary of the
    !  unitsquare.
    !  35        x^2*(1-x)^2 * 2*PI*sin(PI*y)*cos(PI*y)
    !  36        -(2*x*(1-x)^2 - 2*x^2*(1-x))*sin^2(Pi*y)
    !            (these are the first two components of curl(0,0,psi) where
    !             psi(x,y) = x^2(1-x)^2 * sin^2(PI*y) )
    !
    !  37        16 * x * (1 - x) * y * (1 - y)
    !  38        42*x^2 * (2 - y) + 42*sin(2 - 5*y*x) * cos(x + y + 1)
    !  39        (x^2 - 1)^2 (y^2 - 1)y / 4      (zero divergence together with 40)
    !  40        (y^2 - 1)^2 (1 - x^2)x / 4
    !  41        5 x^3 (y-1) + y^3
    !  42        x*(y-0.01) (for solid beam configuration; use together with 43)
    !  43        -0.5*x^2 (for solid beam configuration; use together with 42)
    !  44        sin(c1*x+c2) * sin(c1*y+c2)  (variant of 12/30)
    !  45        cos(c1*x+c2) * cos(c1*y+c2)  (variant of 14/31)
    !  46        -sin(c1*x+c2) * cos(c1*y+c2)  (variant of 25/32)
    !  47        sin(PI*x+0.4) cos(PI*y-0.3) (variant of 25/32)
    !  48        c*x^2*y*sin(c*(x-0.5*y^2))
    !  49        -2*x*cos(c*(x-0.5*y^2)) + c*x^2*sin(c*(x-0.5*y^2))
    !  50        0.05 * sin(2*PI*x)*sin(2*PI*y)
    !            (same as 13, only different factors)
    !  51        sin(PI/2 (x-1)) sin(PI/2 (y-1))
    !            (Circular harmonic function on x^2 + y^2 - 1 = 0)
    !  52        0.05 * cos(2*PI*x)*cos(2*PI*y)
    !            (same as 50, only sin replaced by cos in order to get nonzero values on
    !             the boundary of the unitsquare)


    ! Stokes pairs (for unit square with zero mean pressure):
    ! / 0,  0,  0/ , /23, 23,  0/ , / 6,  0, 27/ , /15, 16, 19/,
    ! /12, 14, 28/ , / 7, 24, 29/ , /25, 26, 19/, /30, 31, 32/, /33, 34, 32/, /35, 36, 32/
    !
    ! Stokes pairs (for square box [-1,1]x[-1,1], zero mean pressure):
    ! /39, 40, 41/
    !
! -->

!<result>
    ! (The result of the function calculation)
    real(DP), dimension(nptsPerEl, nel) :: Dval

!</result>

!<errors>
    ! none
!</errors>
!</function>

!    real(DP) :: daux, daux1

    real(DP) :: dgamma
#ifdef SOLUTION_CAUSING_SERIOUS_PROBLEMS_FOR_STABILISED_Q1_Q1_STOKES
    ! gamma should be set to a multiple of pi. gamma=pi is just fine.
    ! A setting of gamma=3*Pi gives in combination with stabilised Q1/Q1 approaches
    ! serious problems for the pressure solution of coarse grids.
    dgamma = 3.0_DP*SYS_PI
#else
    dgamma = 1.0_DP*SYS_PI
#endif

    ! avoid misleading warnings about uninitialised variables
    Dval(:,:) = 0.0_DP

    select case (cselect)

    case (0) ! u(x,y) = 0.0
      Dval(:,:) = 0.0_DP

    case (1) ! u(x,y) = 0.1 * x
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.1_DP
      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (2) ! u(x,y) = 0.1 * x**2
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(1,:,:) * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.2_DP * Dpts(1,:,:)
      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_XX); Dval(:,:) = 0.2_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (3) ! u(x,y) = 0.1 * y
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(2,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_Y);  Dval(:,:) = 0.1_DP
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (4) ! u(x,y) = 0.1 * y**2
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(2,:,:) * Dpts(2,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_Y);  Dval(:,:) = 0.2_DP * Dpts(2,:,:)
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.2_DP
      end select

    case (5) ! u(x,y) = 4 * x * (1 - x)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 4.0_DP*Dpts(1,:,:) * (1.0_DP - Dpts(1,:,:))
      case (DER_DERIV_X);  Dval(:,:) = 4.0_DP*(1.0_DP - 2.0_DP*Dpts(1,:,:))
      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_XX); Dval(:,:) = -8.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (6) ! u(x,y) = 4 * y * (1 - y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 4.0_DP*Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_Y);  Dval(:,:) = 4.0_DP*(1.0_DP - 2.0_DP*Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = -8.0_DP
      end select

    case (7) ! u(x,y) = x * (1 - x) * y * (1 - y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        Dpts(1,:,:) * (1.0_DP - Dpts(1,:,:)) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) = &
        (-1.0_DP + 2.0_DP*Dpts(1,:,:)) * Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) = &
        Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:)) * (-1.0_DP + 2.0_DP*Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) = (2.0_DP*Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:)))
      case (DER_DERIV_XY); Dval(:,:) = &
        (-1.0_DP + 2.0_DP*Dpts(1,:,:)) * (-1.0_DP + 2.0_DP*Dpts(2,:,:))
      case (DER_DERIV_YY); Dval(:,:) = 2.0_DP*Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:))
      end select

    case (8) ! u(x,y) = -(x * x + y * y) * x
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        -(Dpts(1,:,:) * Dpts(1,:,:) + Dpts(2,:,:) * Dpts(2,:,:)) * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = &
        -(3.0_DP*Dpts(1,:,:) * Dpts(1,:,:) + Dpts(2,:,:) * Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) = -2.0_DP*Dpts(1,:,:) * Dpts(2,:,:)
      case (DER_DERIV_XX); Dval(:,:) = -6.0_DP*Dpts(1,:,:)
      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP*Dpts(2,:,:)
      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP*Dpts(1,:,:)
      end select

    case (9) ! u(x,y) = y^2 * (1 - y)^2 * x
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2 * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = &
        Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2
      case (DER_DERIV_Y);  Dval(:,:) = &
          2.0_DP*Dpts(1,:,:) * (Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2 - Dpts(2,:,:) &
        * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:)))
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = &
          2.0_DP*(Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2 &
        - Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:)))
      case (DER_DERIV_YY); Dval(:,:) = &
        2.0_DP*Dpts(1,:,:) * ((1.0_DP - Dpts(2,:,:))**2 &
        - 4.0_DP*Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:)) + Dpts(2,:,:) * Dpts(2,:,:))
      end select

    case (10) ! u(x,y) = y^2 * (1 - y)^2 * (x - 1)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) =   Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP &
                                     - Dpts(2,:,:))**2 * (Dpts(1,:,:) - 1.0_DP)
      case (DER_DERIV_X);  Dval(:,:) =   Dpts(2,:,:) * Dpts(2,:,:) &
                                     * (1.0_DP - Dpts(2,:,:))**2
      case (DER_DERIV_Y);  Dval(:,:) =   2.0_DP*Dpts(2,:,:)*(2.0_DP*Dpts(2,:,:) - 1.0_DP)&
                                     * (Dpts(2,:,:) - 1.0_DP) * (Dpts(1,:,:) - 1.0_DP)
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 2.0_DP*Dpts(2,:,:) * (2.0_DP*Dpts(2,:,:) - 1.0_DP)&
                                     * (Dpts(2,:,:) - 1.0_DP)
      case (DER_DERIV_YY); Dval(:,:) = 2.0_DP*(1.0_DP - 6.0_DP*Dpts(2,:,:) &
                                     + 6.0_DP*Dpts(2,:,:) * Dpts(2,:,:)) &
                                     * (Dpts(2,:,:) - 1.0_DP)
      end select

!    case (11) ! u1(x,y) =    0.25 * (1/sqrt(2) + x - y) * (1/sqrt(2) - x + y)
!              !           * ((1-sqrt(2))/sqrt(2) + x + y) * ((1+sqrt(2))/sqrt(2) - x - y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =   0.25_DP*(1.0_DP / sqrt(2.0_DP) + Dpts(1,:,:) &
!                              - Dpts(2,:,:)) * (1.0_DP / sqrt(2.0_DP) - Dpts(1,:,:) + Dpts(2,:,:)) &
!                              * ((1.0_DP - sqrt(2.0_DP)) / sqrt(2.0_DP) + Dpts(1,:,:) + Dpts(2,:,:)) &
!                              * ((1.0_DP + sqrt(2.0_DP)) / sqrt(2.0_DP) - Dpts(1,:,:) - Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  -1.5_DP*Dpts(1,:,:) * Dpts(1,:,:) - 1.0_DP*Dpts(2,:,:) * Dpts(2,:,:) * Dpts(1,:,:) &
!                              + 0.5_DP*Dpts(2,:,:) * Dpts(2,:,:) - 0.5_DP*Dpts(2,:,:) + 0.25_DP &
!                              + 1.0_DP*Dpts(1,:,:) * Dpts(2,:,:) + 1.0_DP*Dpts(1,:,:)**3
!      case (DER_DERIV_Y);  Dval(:,:) =  -1.0_DP*Dpts(1,:,:) * Dpts(1,:,:) * Dpts(2,:,:) + 0.5_DP*Dpts(1,:,:) * Dpts(1,:,:) &
!                              - 0.5_DP*Dpts(1,:,:) - 1.5_DP*Dpts(2,:,:) * Dpts(2,:,:) + 0.25_DP &
!                              + 1.0_DP*Dpts(1,:,:) * Dpts(2,:,:) + 1.0_DP*Dpts(2,:,:)**3
!      case (DER_DERIV_XX); Dval(:,:) = -1.0_DP*Dpts(2,:,:) * Dpts(2,:,:) + 1.0_DP*Dpts(2,:,:) + 3.0_DP*Dpts(1,:,:) * Dpts(1,:,:)&
!                              - 3.0_DP*Dpts(1,:,:)
!      case (DER_DERIV_XY); Dval(:,:) = -0.5_DP + 1.0_DP*Dpts(1,:,:) + 1.0_DP*Dpts(2,:,:) - 2.0_DP*Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_YY); Dval(:,:) =  3.0_DP*Dpts(2,:,:) * Dpts(2,:,:) - 3.0_DP*Dpts(2,:,:) - 1.0_DP*Dpts(1,:,:) * Dpts(1,:,:) &
!                               + 1.0_DP*Dpts(1,:,:)
!      end select
!
    case (12) ! u(x,y) = sin(x) * sin(y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) =  sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) =  cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      case (DER_DERIV_XY); Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      case (DER_DERIV_YY); Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      end select
!
!    case (13) ! u(x,y) = 0.05 * sin(4*PI*x)*sin(4*PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  0.05_DP*sin(4.0_DP*SYS_PI*Dpts(1,:,:)) * sin(4.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =   0.2_DP * SYS_PI &
!                              * cos(4.0_DP*SYS_PI*Dpts(1,:,:)) * sin(4.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =   0.2_DP * SYS_PI &
!                              * sin(4.0_DP*SYS_PI*Dpts(1,:,:)) * cos(4.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  -0.8_DP * SYS_PI*SYS_PI &
!                              * sin(4.0_DP*SYS_PI*Dpts(1,:,:)) * sin(4.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =   0.8_DP * SYS_PI*SYS_PI &
!                              * cos(4.0_DP*SYS_PI*Dpts(1,:,:)) * cos(4.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  -0.8_DP * SYS_PI*SYS_PI &
!                              * sin(4.0_DP*SYS_PI*Dpts(1,:,:)) * sin(4.0_DP*SYS_PI*Dpts(2,:,:))
!      end select
!
    case (14) ! u(x,y) = cos(x) * cos(y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) = -cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      case (DER_DERIV_XY); Dval(:,:) =  sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      case (DER_DERIV_YY); Dval(:,:) = -cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      end select
!
!    case (15) ! u(x,y) = cos(PI/2 * (x + y))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                        cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_X);  Dval(:,:) =         -0.5_DP*SYS_PI*sin(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_Y);  Dval(:,:) =         -0.5_DP*SYS_PI*sin(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XX); Dval(:,:) = -0.25_DP*SYS_PI*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) = -0.25_DP*SYS_PI*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_YY); Dval(:,:) = -0.25_DP*SYS_PI*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      end select
!
!    case (16) ! u(x,y) = -cos(PI/2 * (x + y))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                            -cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_X);  Dval(:,:) =           0.5_DP*SYS_PI*sin(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_Y);  Dval(:,:) =           0.5_DP*SYS_PI*sin(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XX); Dval(:,:) = 0.25_DP*SYS_PI*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) = 0.25_DP*SYS_PI*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_YY); Dval(:,:) = 0.25_DP*SYS_PI*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) + Dpts(2,:,:)))
!      end select
!
!    case (17) ! u(x,y) = 2 * cos(x) * sin(y) - 2 * (1 - cos(1)) * sin(1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  2.0_DP*cos(Dpts(1,:,:)) * sin(Dpts(2,:,:)) &
!                              -2.0_DP*(1.0_DP - cos(1.0_DP)) * sin(1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) = -2.0_DP*sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  2.0_DP*cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -2.0_DP*cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP*sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP*cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      end select
!
!    case (18) ! u(x,y) = 8 * (1 - x)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  8.0_DP*(1.0_DP - Dpts(1,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = -8.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (19) ! u(x,y) = sin(PI/2 * (x - y))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                              sin(0.5_DP*SYS_PI*(Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_X);  Dval(:,:) =            0.5_DP*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_Y);  Dval(:,:) =           -0.5_DP*SYS_PI*cos(0.5_DP*SYS_PI*(Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_XX); Dval(:,:) = -0.25_DP*SYS_PI*SYS_PI*sin(0.5_DP*SYS_PI*(Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) =  0.25_DP*SYS_PI*SYS_PI*sin(0.5_DP*SYS_PI*(Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_YY); Dval(:,:) = -0.25_DP*SYS_PI*SYS_PI*sin(0.5_DP*SYS_PI*(Dpts(1,:,:) - Dpts(2,:,:)))
!      end select
!
!    case (20) ! u(x,y) = -(x * x + y * y) * x * x
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = -(Dpts(1,:,:) * Dpts(1,:,:) + Dpts(2,:,:) * Dpts(2,:,:)) * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = -4.0_DP*Dpts(1,:,:)**3 + 2.0_DP*Dpts(1,:,:) * Dpts(2,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) = -2.0_DP*Dpts(1,:,:) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XX); Dval(:,:) = -12.0_DP*Dpts(1,:,:) * Dpts(1,:,:) - 2.0_DP*Dpts(2,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XY); Dval(:,:) = -4.0_DP*Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP*Dpts(1,:,:) * Dpts(1,:,:)
!      end select
!
!    case (21) ! u(x,y) = cos(Pi/2)*x - sin(Pi/2)*y - x
!      if (present(dparam)) then
!        daux = 0.125_DP*dparam * SYS_PI
!      else
!        daux = 0.5_DP*SYS_PI
!      endif
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = cos(daux) * Dpts(1,:,:) - sin(daux) * Dpts(2,:,:) - Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = cos(daux) - 1.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) = - sin(daux)
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (22) ! u(x,y) = sin(Pi/2)*x + cos(Pi/2)*y - y
!      if (present(dparam)) then
!        daux = 0.125_DP*dparam * SYS_PI
!      else
!        daux = 0.5_DP*SYS_PI
!      endif
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = sin(daux) * Dpts(1,:,:) + cos(daux) * Dpts(2,:,:) - Dpts(2,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = sin(daux)
!      case (DER_DERIV_Y);  Dval(:,:) = cos(daux) - 1.0_DP
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (23) ! u(x,y) = 1.0
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 1.0_DP
!      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (24) ! u(x,y) = -(2x - 1)(2y^3 - 3y^2) / 6
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = - (2.0_DP*Dpts(1,:,:) - 1.0_DP) * (2.0_DP*Dpts(2,:,:)**3 - 3.0_DP*Dpts(2,:,:)**2) / 6.0_DP
!      case (DER_DERIV_X);  Dval(:,:) = - Dpts(2,:,:)**2 * (2.0_DP*Dpts(2,:,:) - 3.0_DP) / 3.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) = -Dpts(2,:,:) * (Dpts(2,:,:) - 1.0_DP) * (2.0_DP*Dpts(1,:,:) - 1.0_DP)
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP*Dpts(2,:,:) * (Dpts(2,:,:) - 1.0_DP)
!      case (DER_DERIV_YY); Dval(:,:) = -(2.0_DP*Dpts(1,:,:) - 1.0_DP) * (2.0_DP*Dpts(2,:,:) - 1.0_DP)
!      end select
!
!    case (25) ! u(x,y) = sin(x) cos(y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      end select
!
!    case (26) ! u(x,y) = -cos(x) sin(y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      end select
!
!    case (27) ! u(x,y) = 4 - 8x
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  4.0_DP - 8.0_DP*Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = -8.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
    case (28) ! u(x,y) = 2 cos(x) sin(y) - 2 sin(1) + 2 sin(1) cos(1)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) =   2.0_DP*cos(Dpts(1,:,:)) * sin(Dpts(2,:,:)) &
                                - 2.0_DP*sin(1.0_DP) + 2.0_DP*sin(1.0_DP) * cos(1.0_DP)
      case (DER_DERIV_X);  Dval(:,:) = -2.0_DP*sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) =  2.0_DP*cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) = -2.0_DP*cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP*sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP*cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
      end select
!
!    case (29) ! u(x,y) = xy - 1/4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = Dpts(1,:,:) * Dpts(2,:,:) - 0.25_DP
!      case (DER_DERIV_X);  Dval(:,:) = Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) = Dpts(1,:,:)
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 1.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (30) ! u(x,y) = sin(PI*x - 7/10) * sin(PI*y + 1/5)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(SYS_PI*Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI*Dpts(2,:,:) + 0.2_DP)
!      case (DER_DERIV_X);  Dval(:,:) =  cos(SYS_PI*Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_Y);  Dval(:,:) =  sin(SYS_PI*Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI*SYS_PI
!      case (DER_DERIV_XY); Dval(:,:) =  cos(SYS_PI*Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI*SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI*SYS_PI
!      end select
!
!    case (31) ! u(x,y) = cos(PI*x - 7/10) cos(PI*y + 1/5)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  cos(SYS_PI*Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI*Dpts(2,:,:) + 0.2_DP)
!      case (DER_DERIV_X);  Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_Y);  Dval(:,:) = -cos(SYS_PI*Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -cos(SYS_PI*Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI*SYS_PI
!      case (DER_DERIV_XY); Dval(:,:) =  sin(SYS_PI*Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI*SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) = -cos(SYS_PI*Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI*Dpts(2,:,:) + 0.2_DP) * SYS_PI*SYS_PI
!      end select
!
!    case (32) ! u(x,y) = sin(x) cos(y) + (cos(1) - 1) sin(1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:)) + (cos(1.0_DP) - 1.0_DP) * sin(1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      end select
!
!    case (33) ! u(x,y) = -sin(gamma x) * sin(gamma y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                   -sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = -dgamma *          cos(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -dgamma *          sin(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  dgamma * dgamma * sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -dgamma * dgamma * cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  dgamma * dgamma * sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      end select
!
!    case (34) ! u(x,y) = -cos(gamma x) * cos(gamma y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                   -cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  dgamma *          sin(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  dgamma *          cos(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  dgamma * dgamma * cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -dgamma * dgamma * sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  dgamma * dgamma * cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      end select
!
!    case (35) ! u(x,y) = x^2*(1-x)^2 * 2*PI*sin(PI*y)*cos(PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))**2*2.0_DP*SYS_PI*sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = 4.0_DP*(Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 - Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) &
!                                * SYS_PI*sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = 2.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))**2 &
!                                * (SYS_PI**2*cos(SYS_PI*Dpts(2,:,:))**2 - SYS_PI**2*sin(SYS_PI*Dpts(2,:,:))**2)
!      case (DER_DERIV_XX); Dval(:,:) = (4.0_DP*(1.0_DP-Dpts(1,:,:))**2 - 16.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:)) + 4.0_DP*Dpts(1,:,:)**2) &
!                                * SYS_PI*sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = 4.0_DP*SYS_PI**2 * ((Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 &
!                                - Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) * cos(SYS_PI*Dpts(2,:,:))**2 &
!                                - (Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 + Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) * sin(SYS_PI*Dpts(2,:,:))**2)
!      case (DER_DERIV_YY); Dval(:,:) =-8.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))**2 &
!                                * SYS_PI**3*cos(SYS_PI*Dpts(2,:,:))*sin(SYS_PI*Dpts(2,:,:))
!      end select
!
!    case (36) ! u(x,y) = -(2*x*(1-x)^2 - 2*x^2*(1-x))*sin^2(PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = -(2.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 &
!                                - 2.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:)))*sin(SYS_PI*Dpts(2,:,:))**2
!      case (DER_DERIV_X);  Dval(:,:) = -(2.0_DP*(1.0_DP-Dpts(1,:,:))**2 - 8.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:)) &
!                                + 2.0_DP*Dpts(1,:,:)**2) * sin(SYS_PI*Dpts(2,:,:))**2
!      case (DER_DERIV_Y);  Dval(:,:) = -2.0_DP*(2.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2-2.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) &
!                                * sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))*SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -(-12.0_DP + 24.0_DP*Dpts(1,:,:))*sin(SYS_PI*Dpts(2,:,:))**2
!      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP*(2.0_DP*(1.0_DP-Dpts(1,:,:))**2 -8.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))&
!                                + 2.0_DP*Dpts(1,:,:)**2) * sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))*SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) =   4.0_DP*(Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 - Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:)))*SYS_PI**2 &
!                                * (-cos(SYS_PI*Dpts(2,:,:))**2 + sin(SYS_PI*Dpts(2,:,:))**2)
!      end select
!
!    case (37) ! u(x,y) = 16 * x * (1 - x) * y * (1 - y)
!      Dval(:,:) = 0_DP
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = Dpts(1,:,:) * (1.0_DP - Dpts(1,:,:)) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = (-1.0_DP + 2.0_DP*Dpts(1,:,:)) * Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:)) * (-1.0_DP + 2.0_DP*Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = (2.0_DP*Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) = (-1.0_DP + 2.0_DP*Dpts(1,:,:)) * (-1.0_DP + 2.0_DP*Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = 2.0_DP*Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:))
!      end select
!      Dval(:,:) = Dval(:,:)*16.0_DP
!
!    case (38) ! u(x,y) = 42*x^2*(2 - y) + 42*sin(2 - 5*y*x)*cos(x + y + 1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 42.0_DP*Dpts(1,:,:)** 2 * (2.0_DP - Dpts(2,:,:)) &
!                                + 42.0_DP*sin(2.0_DP - 5.0_DP*Dpts(2,:,:)* Dpts(1,:,:)) * cos(Dpts(1,:,:) + Dpts(2,:,:) + 1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) =   84.0_DP*Dpts(1,:,:)* (2.0_DP - Dpts(2,:,:)) &
!                                - 210.0_DP*cos(-2.0_DP + 5.0_DP*Dpts(2,:,:)* Dpts(1,:,:)) * Dpts(2,:,:) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                - 42.0_DP*sin(2.0_DP - 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  -42.0_DP*Dpts(1,:,:) ** 2 &
!                                - 210.0_DP*cos(-2.0_DP + 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(1,:,:) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                - 42.0_DP*sin(2.0_DP - 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =   168.0_DP - 84.0_DP*Dpts(2,:,:) &
!                                + 1050.0_DP*sin(-2.0_DP + 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) ** 2 * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 420.0_DP*cos(-2.0_DP + 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                -  42.0_DP*sin(2.0_DP - 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =   -84.0_DP*Dpts(2,:,:) + 1050.0_DP*sin(-2.0_DP &
!                                + 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) * Dpts(1,:,:) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                - 210.0_DP*cos(-2.0_DP + 5.0_DP* Dpts(2,:,:)* Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 210.0_DP*cos(-2.0_DP + 5.0_DP* Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(1,:,:) * sin(1.0_DP +Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 210.0_DP*cos(-2.0_DP + 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                -  42.0_DP*sin(2.0_DP - 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  1050.0_DP*sin(-2.0_DP &
!                                + 5.0_DP*Dpts(2,:,:)* Dpts(1,:,:)) * Dpts(1,:,:) ** 2 * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 420.0_DP*cos(-2.0_DP + 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(1,:,:)* sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                -  42.0_DP*sin(2.0_DP - 5.0_DP*Dpts(2,:,:) * Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      end select
!
!    case (39) ! u(x,y) = (x^2 - 1)^2 (y^2 - 1)y / 4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 0.25_DP*(Dpts(1,:,:)**2 - 1.0_DP)**2 * (Dpts(2,:,:)**2 - 1.0_DP) * Dpts(2,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = (Dpts(1,:,:)**2 - 1.0_DP) * (Dpts(2,:,:)**2 - 1.0_DP) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) = 0.25_DP*(Dpts(1,:,:)**2 - 1.0_DP)**2 * (3.0_DP*Dpts(2,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_XX); Dval(:,:) = Dpts(2,:,:) * (Dpts(2,:,:)**2 - 1.0_DP) * (3.0_DP*Dpts(1,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_XY); Dval(:,:) = Dpts(1,:,:) * (Dpts(1,:,:)**2 - 1.0_DP) * (3.0_DP*Dpts(2,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_YY); Dval(:,:) = 1.50_DP * (Dpts(1,:,:)**2 - 1)**2 * Dpts(2,:,:)
!      end select
!
!    case (40) ! u(x,y) = (y^2 - 1)^2 (1 - x^2)x / 4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 0.25_DP*(Dpts(2,:,:)**2 - 1.0_DP)**2 * (1.0_DP - Dpts(1,:,:)**2) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = -0.25_DP*(Dpts(2,:,:)**2 - 1.0_DP)**2 * (3.0_DP*Dpts(1,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_Y);  Dval(:,:) = (Dpts(2,:,:)**2 - 1.0_DP) * (1.0_DP - Dpts(1,:,:)**2) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XX); Dval(:,:) = -1.50_DP * (Dpts(2,:,:)**2 - 1)**2 * Dpts(1,:,:)
!      case (DER_DERIV_XY); Dval(:,:) = Dpts(2,:,:) * (Dpts(2,:,:)**2 - 1.0_DP) * (1.0_DP - 3.0_DP*Dpts(1,:,:)**2)
!      case (DER_DERIV_YY); Dval(:,:) = Dpts(1,:,:) * (Dpts(1,:,:)**2 - 1.0_DP) * (1.0_DP - 3.0_DP*Dpts(2,:,:)**2)
!      end select
!
!    case (41) ! u(x,y) = 5 x^3 (y-1) + y^3
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 5.0_DP*Dpts(1,:,:)**3 * (Dpts(2,:,:) - 1.0_DP) + Dpts(2,:,:)**3
!      case (DER_DERIV_X);  Dval(:,:) = 15.0_DP*Dpts(1,:,:)**2 * (Dpts(2,:,:) - 1.0_DP)
!      case (DER_DERIV_Y);  Dval(:,:) = 5.0_DP*Dpts(1,:,:)**3 + 3.0_DP*Dpts(2,:,:)**2
!      case (DER_DERIV_XX); Dval(:,:) = 30.0_DP*Dpts(1,:,:) * (Dpts(2,:,:) - 1.0_DP)
!      case (DER_DERIV_XY); Dval(:,:) = 15.0_DP*Dpts(1,:,:)**2
!      case (DER_DERIV_YY); Dval(:,:) = 6.0 * Dpts(2,:,:)
!      end select
!
!    case (42) ! u(x,y) = x*(y-0.01)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  Dpts(1,:,:) * (Dpts(2,:,:)-0.01_DP)
!      case (DER_DERIV_X);  Dval(:,:) =      (Dpts(2,:,:)-0.01_DP)
!      case (DER_DERIV_Y);  Dval(:,:) =  Dpts(1,:,:)
!      case (DER_DERIV_XX); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  1.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (43) ! u(x,y) = -0.5*x^2
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  -0.5_DP*Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) =  -Dpts(1,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  -1.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (44) ! u(x,y) = sin(c1*x + c2) * sin(c1*y + c2)
!      daux = 2.0_DP*SYS_PI
!      daux1 = 0.02_DP*daux
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_X);  Dval(:,:) =  daux *        cos(daux * Dpts(1,:,:)+ daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_Y);  Dval(:,:) =  daux *        sin(daux *Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XX); Dval(:,:) = -daux * daux * sin(daux * Dpts(1,:,:)+ daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XY); Dval(:,:) =  daux * daux * cos(daux * Dpts(1,:,:)+ daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_YY); Dval(:,:) = -daux * daux * sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      end select
!
!    case (45) ! u(x,y) = cos(c1*x + c2) * cos(c1*y + c2)
!      daux = 2.0_DP*SYS_PI
!      daux1 = 0.02_DP*daux
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_X);  Dval(:,:) = -daux *        sin(daux *Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_Y);  Dval(:,:) = -daux *        cos(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XX); Dval(:,:) = -daux * daux * cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:)+ daux1)
!      case (DER_DERIV_XY); Dval(:,:) =  daux * daux * sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_YY); Dval(:,:) = -daux * daux * cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      end select
!
!    case (46) ! u(x,y) = -sin(c1*x + c2) * cos(c1*y + c2)
!      daux = (2.0_DP*SYS_PI)**2
!      daux1 = 0.02_DP*daux
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =               -sin(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_X);  Dval(:,:) = -daux *        cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_Y);  Dval(:,:) =  daux *        sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XX); Dval(:,:) =  daux * daux * sin(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XY); Dval(:,:) =  daux * daux * cos(daux * Dpts(1,:,:)+ daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_YY); Dval(:,:) =  daux * daux * sin(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      end select
!
!    case (47) ! u(x,y) = sin(PI*x+0.4) cos(PI*y-0.3)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP)
!      case (DER_DERIV_X);  Dval(:,:) =  cos(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI
!      case (DER_DERIV_Y);  Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * sin(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI*SYS_PI
!      case (DER_DERIV_XY); Dval(:,:) = -cos(SYS_PI*Dpts(1,:,:) + 0.4_DP) * sin(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI*SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI*SYS_PI
!      end select
!
!
!    case (48) ! u(x,y) = c*x^2*y*sin(c*(x-0.5*y^2))
!      daux = 20.0_DP*SYS_PI
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =   daux*Dpts(1,:,:)**2*Dpts(2,:,:) * sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_X);  Dval(:,:) =   2*daux*Dpts(1,:,:)*Dpts(2,:,:)*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_Y);  Dval(:,:) =  -daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)**2*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + daux*Dpts(1,:,:)**2*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XX); Dval(:,:) =   (2*daux*Dpts(2,:,:) - daux**3*Dpts(1,:,:)**2*Dpts(2,:,:))*sin(daux*(Dpts(1,:,:) &
!                                -0.5_DP*Dpts(2,:,:)**2)) + 4*daux**2*Dpts(1,:,:)*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XY); Dval(:,:) =   (daux**2*Dpts(1,:,:)**2 - 2*daux**2*Dpts(1,:,:)*Dpts(2,:,:)**2)*cos(daux*(Dpts(1,:,:) &
!                                -0.5_DP*Dpts(2,:,:)**2)) + (2*daux*Dpts(1,:,:) + daux**3*Dpts(1,:,:)**2*Dpts(2,:,:)**2)*sin(daux*(Dpts(1,:,:)&
!                                -0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_YY); Dval(:,:) =  -daux**3*Dpts(1,:,:)**2*Dpts(2,:,:)**3*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                - 3*daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      end select
!
!    case (49) ! u(x,y) = -2*x*cos(c*(x-0.5*y^2)) + c*x^2*sin(c*(x-0.5*y^2))
!      daux = 20.0_DP*SYS_PI
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =   -2*Dpts(1,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + daux*Dpts(1,:,:)**2*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_X);  Dval(:,:) =   (daux**2*Dpts(1,:,:)**2 - 2)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + 4*daux*Dpts(1,:,:)*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_Y);  Dval(:,:) =  -2*daux*Dpts(1,:,:)*Dpts(2,:,:)*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                - daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XX); Dval(:,:) =   (6*daux - daux**3*Dpts(1,:,:)**2) * sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + 6*daux**2*Dpts(1,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XY); Dval(:,:) =   (daux**3*Dpts(1,:,:)**2*Dpts(2,:,:) - 2*daux*Dpts(2,:,:))*sin(daux*(Dpts(1,:,:)&
!                                -0.5_DP*Dpts(2,:,:)**2)) - 4*daux**2*Dpts(1,:,:)*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_YY); Dval(:,:) =   (2*daux**2*Dpts(1,:,:)*Dpts(2,:,:)**2  - daux**2*Dpts(1,:,:)**2)*cos(daux*(Dpts(1,:,:)&
!                                -0.5_DP*Dpts(2,:,:)**2)) - (daux**3*Dpts(1,:,:)**2*Dpts(2,:,:)**2 &
!                                + 2*daux*Dpts(1,:,:))*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      end select
!
!    case (50) ! u(x,y) = 0.05 * sin(2*PI*x)*sin(2*PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  0.05_DP*sin(2.0_DP*SYS_PI*Dpts(1,:,:)) * sin(2.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =   0.1_DP * SYS_PI &
!                              * cos(2.0_DP*SYS_PI*Dpts(1,:,:)) * sin(2.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =   0.1_DP * SYS_PI &
!                              * sin(2.0_DP*SYS_PI*Dpts(1,:,:)) * cos(2.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  -0.2_DP * SYS_PI*SYS_PI &
!                              * sin(2.0_DP*SYS_PI*Dpts(1,:,:)) * sin(2.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =   0.2_DP * SYS_PI*SYS_PI &
!                              * cos(2.0_DP*SYS_PI*Dpts(1,:,:)) * cos(2.0_DP*SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  -0.2_DP * SYS_PI*SYS_PI &
!                              * sin(2.0_DP*SYS_PI*Dpts(1,:,:)) * sin(2.0_DP*SYS_PI*Dpts(2,:,:))
!      end select
!
!    case (51) ! u(x,y) = sin(PI/2 (x-1)) sin(PI/2 (y-1))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(0.5_DP*SYS_PI*(Dpts(1,:,:)-1)) * sin(0.5_DP*SYS_PI*(Dpts(2,:,:)-1))
!      case (DER_DERIV_X);  Dval(:,:) =  0.5_DP*SYS_PI &
!                               * cos(0.5_DP*SYS_PI*(Dpts(1,:,:)-1)) * sin(0.5_DP*SYS_PI*(Dpts(2,:,:)-1))
!      case (DER_DERIV_Y);  Dval(:,:) =  0.5_DP*SYS_PI &
!                               * sin(0.5_DP*SYS_PI*(Dpts(1,:,:)-1)) * cos(0.5_DP*SYS_PI*(Dpts(2,:,:)-1))
!      case (DER_DERIV_XX); Dval(:,:) = -0.25_DP*SYS_PI*SYS_PI &
!                               * sin(0.5_DP*SYS_PI*(Dpts(1,:,:)-1)) * sin(0.5_DP*SYS_PI*(Dpts(2,:,:)-1))
!      case (DER_DERIV_XY); Dval(:,:) =  0.25_DP*SYS_PI*SYS_PI &
!                               * cos(0.5_DP*SYS_PI*(Dpts(1,:,:)-1)) * cos(0.5_DP*SYS_PI*(Dpts(2,:,:)-1))
!      case (DER_DERIV_YY); Dval(:,:) = -0.25_DP*SYS_PI*SYS_PI &
!                               * sin(0.5_DP*SYS_PI*(Dpts(1,:,:)-1)) * sin(0.5_DP*SYS_PI*(Dpts(2,:,:)-1))
!      end select
!
    case (52) ! u(x,y) = 0.05 * cos(2*PI*x)*cos(2*PI*y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        0.05_DP*cos(2.0_DP*SYS_PI *Dpts(1,:,:)) * cos(2.0_DP*SYS_PI*Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) = -0.1_DP * SYS_PI &
        * sin(2.0_DP*SYS_PI*Dpts(1,:,:)) * cos(2.0_DP*SYS_PI*Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) =   -0.1_DP * SYS_PI &
        * cos(2.0_DP*SYS_PI*Dpts(1,:,:)) * sin(2.0_DP*SYS_PI*Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) =  -0.2_DP * SYS_PI*SYS_PI &
        * cos(2.0_DP*SYS_PI*Dpts(1,:,:)) * cos(2.0_DP*SYS_PI*Dpts(2,:,:))
      case (DER_DERIV_XY); Dval(:,:) =   0.2_DP * SYS_PI*SYS_PI &
        * sin(2.0_DP*SYS_PI*Dpts(1,:,:)) * sin(2.0_DP*SYS_PI*Dpts(2,:,:))
      case (DER_DERIV_YY); Dval(:,:) =  -0.2_DP * SYS_PI*SYS_PI &
        * cos(2.0_DP*SYS_PI *Dpts(1,:,:)) * cos(2.0_DP*SYS_PI*Dpts(2,:,:))
      end select
!
!    case (53) ! u(x,y) = -x^3*y
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  -Dpts(1,:,:) * Dpts(1,:,:)*Dpts(1,:,:)*Dpts(2,:,:)
!      case (DER_DERIV_X);  Dval(:,:) =  -3.0_DP*Dpts(1,:,:) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) =  -Dpts(1,:,:) * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_XX); Dval(:,:) =  -6.0_DP*Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XY); Dval(:,:) =  -3.0_DP*Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (54) ! u(x,y) = 1/3*x^4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  (1.0_DP/4.0_DP) * Dpts(1,:,:) * Dpts(1,:,:) * &
!                                            Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) =  Dpts(1,:,:) * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  3.0_DP*Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
    end select

  end function elast_danalyticFunction

end module elasticity_basic

