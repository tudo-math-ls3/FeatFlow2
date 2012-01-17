!##############################################################################
!# ****************************************************************************
!# <name> LinearSolverAutoInitialise </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to initialise a linear solver for a given
!# situation. On one hand, it provides initialisation routines for standard
!# solver settings. On the other hand, there is a parser included which can
!# read user-defined text files from hard disc that configure the setting of
!# a solver.
!#
!# The following routines can be found here:
!#
!# 1.) linsolinit_initFromFile
!#     -> Initialise a linear solver node by reading parameters from a
!#        INI/DAT file
!#
!# 2.) linsolinit_initParams
!#     -> Read parameters of a linear solver from a section of a DAT file
!#        and write them into a solver structure.
!#
!# </purpose>
!##############################################################################

module linearsolverautoinitialise

!$use omp_lib
  use fsystem
  use genoutput
  use spatialdiscretisation
  use linearsolver
  use multilevelprojection
  use filtersupport
  use paramlist

  implicit none
  
  private
  
  public :: linsolinit_CC2DMultigrid
  public :: linsolinit_initFromFile
  public :: linsolinit_initParams

contains

  ! ***************************************************************************

!<subroutine>
  
  subroutine linsolinit_CC2DMultigrid(p_rsolverNode, isolverType, nlevels, &
                        niterationsMin, niterationsMax,daccuracyAbs,daccuracyRel,  &
                        nprecSteps, nsmoothingSteps, &
                        nmaxCoarseGridSteps, dcoarseGridAccuracyAbs, &
                        dcoarseGridAccuracyRel)
  
  !<description>
  
  ! This routine builds a solver node for the linear solver for CC2D-like
  ! applications. The following combinations are supported:
  ! 1.) Multigrid solver, VANKA smoother, VANKA coarse grid solver
  ! 2.) BiCGStab-solver, Multigrid preconditioner,
  !     VANKA smoother, VANKA coarse grid solver
  ! The caller must attach the matrices of all levels to the solver manually
  ! by calling linsol_setMatrices. Afterwards the linear solver can be started.
  
  !</description>
  
  !<input>
  
  ! Type of solver structure which should be build up.
  ! 1 = Multigrid solver, VANKA smoother, VANKA coarse grid solver
  ! 2 = BiCGStab-solver, Multigrid preconditioner, VANKA smoother,
  !     VANKA coarse grid solver
  integer, intent(in)               :: isolverType

  ! Number of levels
  integer, intent(in)               :: nlevels

  ! Minimum number of solver iterations
  integer, intent(in)               :: niterationsMin

  ! Maximum number of solver iterations
  integer, intent(in)               :: niterationsMax

  ! If solver conbination is BiCGStab with MG preconditioning, number of
  ! steps multigrid should perform. Otherwise ignored
  integer, intent(in)               :: nprecSteps

  ! Absolute accuracy of the solver
  real(DP), intent(in)                   :: daccuracyAbs

  ! Relativew accuracy of the solver
  real(DP), intent(in)                   :: daccuracyRel

  ! Number of pre- and postsmoothing steps
  integer, intent(in)               :: nsmoothingSteps
  
  ! Absolute accuracy on the coarse grid
  real(DP), intent(in)                   :: dcoarseGridAccuracyAbs

  ! Relative accuracy on the coarse grid
  real(DP), intent(in)                   :: dcoarseGridAccuracyRel
  
  ! Maximum iterations on the coarse grid
  integer, intent(in)               :: nmaxCoarseGridSteps

  ! A list of spatial discretisation structures for the three equations
  ! that are supported by CC2D: x-velocity, y-velocity, pressure
  type(t_spatialDiscretisation), dimension(3) :: RspatialDiscr
  
!</input>
  
!<output>
  
  ! The solver node that identifies the solver. The caller must attach
  ! level-dependent data (matrix information) to it with the standard
  ! routines in the module LinearSolver. Afterwards the structure
  ! can be used to solve the problem.
  
  type(t_linsolNode),pointer :: p_rsolverNode
  
!</output>
  
!</subroutine>

  ! local variables
  type(t_linsolNode),pointer :: p_rmgSolver
  integer               :: ilevel
  type(t_linsolNode),pointer :: p_rpreSmoother
  type(t_linsolNode),pointer :: p_rpostSmoother
  type(t_linsolNode),pointer :: p_rcoarseGridSolver
  type(t_linsolMGLevelInfo), pointer :: p_rlevelInfo
  type(t_interlevelProjectionBlock) :: rprojection
  
  ! Create the solver node - either BiCGStab or Multigrid.
  ! If we create BiCGStab, attach multigrid as preconditioner.
  
  call linsol_initMultigrid (p_rmgSolver)
  
  if (isolverType .eq. 1) then
  
    p_rsolverNode => p_rmgSolver
    
  else
  
    call linsol_initBiCGStab (p_rsolverNode,p_rmgSolver)
    
    ! Configure MG preconditioning as a fixed number of MG steps
    ! without checking the residuum.
    
    p_rmgSolver%nminIterations = nprecSteps
    p_rmgSolver%nmaxIterations = nprecSteps
    p_rmgSolver%iresCheck      = NO
    
  end if
  
  ! Configure the solver
  p_rsolverNode%nminIterations = niterationsMin
  p_rsolverNode%nmaxIterations = niterationsMax
  p_rsolverNode%depsRel = daccuracyRel
  p_rsolverNode%depsAbs = daccuracyAbs
  
  ! Initialise a standard interlevel projection structure for all levels.
  call mlprj_initProjectionDirect (rprojection,RspatialDiscr)
  
  ! Continue to configure MG by accessing p_rmgSolver.
  ! Loop through the levels.
  
  do ilevel = 1,nlevels
  
    ! On the lowest level create a coarse grid solver structure
    if (ilevel .eq. 1) then
      call linsol_initVANKA (p_rcoarseGridSolver)
      p_rcoarseGridSolver%depsRel = dcoarseGridAccuracyRel
      p_rcoarseGridSolver%depsAbs = dcoarseGridAccuracyAbs
      p_rcoarseGridSolver%nmaxIterations = nmaxCoarseGridSteps
    else
      nullify(p_rcoarseGridSolver)
    end if
    
    ! Create pre- and postsmoother structure on the current level
    call linsol_initVANKA (p_rpreSmoother)
    call linsol_initVANKA (p_rpostSmoother)
    
    ! Configure the structures to form a smoother. A smoother is a solver
    ! that iterates a finite number of times without respecting the
    ! residuum.
    
    call linsol_convertToSmoother (p_rpreSmoother,nsmoothingSteps)
    call linsol_convertToSmoother (p_rpostSmoother,nsmoothingSteps)
    
    ! Create the level, attrach it to the solver and proceed to the
    ! next level
    call linsol_addMultigridLevel (p_rlevelInfo,p_rmgSolver, rprojection,&
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver)
    
  end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine linsolinit_initFromFile (p_rsolverNode, rparamList, &
                                                ssolverName, nlevels, RfilterChain, &
                                                rinterlevelProjection)
  
!<description>
  ! This routine creates a new linear solver node p_rsolverNode of the
  ! heap. The configuration on this solver is set up according to the
  ! parameters in the parameter list rparamList (which can be read e.g.
  ! from an INI/DAT file). The string ssolverName specifies the name of a
  ! section in the parameter list rparamList that serves as 'head' of the
  ! solver, i.e. that defines the 'main' linear solver.
  ! The routine automatically creates all sub-solvers (preconditioners,
  ! smoothers,...) by evaluating the parameters in rparamList,
!</description>

!<input>
  ! The parameter list that contains the whole solver configuration.
  type(t_parlist), intent(in) :: rparamList
  
  ! The name of a section in rparamList that contains the configuration of
  ! the main linear solver.
  character(LEN=*), intent(in) :: ssolverName
  
  ! Number of levels in the discretisation.
  integer, intent(in) :: nlevels

  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the
  ! iteration. If not, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  type(t_filterChain), dimension(:), intent(in), target, optional :: RfilterChain

  ! OPTIONAL: An interlevel projection structure that configures the projection
  ! between the solutions on a finer and a coarser grid. The structure
  ! must have been initialised with mlprj_initProjection.
  !
  ! Note that this structure is level-independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway), so the same structure can be used
  ! to initialise all levels!
  !
  ! The structure must be present if any kind of Multigrid solver is used.
  ! If the application knowns that there is no MG solver, the parameter
  ! can be ommitted.
  type(t_interlevelProjectionBlock), optional :: rinterlevelProjection
!</input>

!<output>
  ! A pointer to a new linear solver node on the heap, initialised
  ! by the configuration in the parameter list.
  type(t_linsolNode), pointer :: p_rsolverNode
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    character(LEN=SYS_STRLEN) :: spreconditioner,spresmoother, spostsmoother
    character(LEN=SYS_STRLEN) :: scoarsegridsolver,sString
    integer :: isolverType,isolverSubtype
    type(t_linsolNode), pointer :: p_rpreconditioner,p_rpresmoother,p_rpostsmoother
    type(t_linsolNode), pointer :: p_rcoarsegridsolver
    integer :: i1,ilev,ikrylowDim
    real(DP) :: d1
    type(t_filterChain), dimension(:), pointer :: p_Rfilter
    type(t_linsolMGLevelInfo), pointer     :: p_rlevelInfo
    type(t_linsolMG2LevelInfo), pointer     :: p_rlevelInfo2

    nullify(p_Rfilter)
    if (present(RfilterChain)) then
      p_Rfilter => RfilterChain
    end if

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, ssolverName, p_rsection)
    
    if (.not. associated(p_rsection)) then
      call output_line ('Cannot create linear solver; no section ''' &
          // trim(ssolverName) // '''!', OU_CLASS_ERROR, OU_MODE_STD, &
          'linsolinit_initFromFile')
      call sys_halt()

    end if
    
    ! Ok, we have the section where we can get parameters from.
    ! Get the solver type we should set up.
    ! Let us hope the parameter 'isolverType' exists! That one is
    ! mandatory; if not, the get-routine will stop.
    call parlst_getvalue_int (p_rsection, 'isolverType', isolverType)

    ! Try to get the solver subtype from the parameter list.
    ! This allows switching between different variants of the same
    ! basic algorithm (e.g. VANKA)
    call parlst_getvalue_int (p_rsection, 'isolverSubtype', isolverSubtype,0)

    ! Many solvers support preconditioners - we try to fetch the
    ! name of the preconditioner in-advance to prevent code
    ! duplication
    call parlst_getvalue_string (p_rsection, 'spreconditioner', sString,'')
    spreconditioner = ''
    if (sString .ne. '') read (sString,*) spreconditioner
    
    ! Now comes a biiig select for all the different types of solvers
    ! that are supported by this routine.
    select case (isolverType)
    
    case (LINSOL_ALG_DEFCORR)
      ! Defect correction
      !
      ! Initialise a solver node for the preconditioner - if there is one.
      nullify(p_rpreconditioner)
      if (spreconditioner .ne. '') then
        call linsolinit_initFromFile (p_rpreconditioner,rparamList,&
                                      spreconditioner,nlevels,RfilterChain)
      end if
      ! Init the solver node
      call linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_Rfilter)
      
    case (LINSOL_ALG_JACOBI)
      ! Jacobi solver
      !
      ! Init the solver node
      call linsol_initJacobi (p_rsolverNode)
      
    case (LINSOL_ALG_SOR)
      ! SOR/GS solver
      !
      ! Init the solver node. domega is set to 1.0 as standard
      ! to activate GS. The actual domega is initialised later if specified
      ! in the DAT file.
      call linsol_initSOR (p_rsolverNode,1.0_DP)
      
    case (LINSOL_ALG_SSOR)
      ! SSOR solver
      !
      ! Init the solver node.
      ! Parameters:
      !  iscale = 1   -> Scale preconditioned vector according to actual formula
      !                  in the literature
      !         = 0   -> no scaling, original FEAT implementation
      call parlst_getvalue_int (p_rsection, 'bscale', i1, 0)
      call linsol_initSSOR (p_rsolverNode,bscale = i1 .eq. 1)
      
    case (LINSOL_ALG_BICGSTAB)
      ! BiCGStab solver
      !
      ! Initialise a solver node for the preconditioner - if there is one.
      nullify(p_rpreconditioner)
      if (spreconditioner .ne. '') then
        call linsolinit_initFromFile (p_rpreconditioner,rparamList,&
                                      spreconditioner,nlevels,RfilterChain)
      end if

      ! Init the solver node
      call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_Rfilter)
      
    case (LINSOL_ALG_GMRES)
      ! GMRES solver
      !
      ! Initialise a solver node for the preconditioner - if there is one.
      nullify(p_rpreconditioner)
      if (spreconditioner .ne. '') then
        call linsolinit_initFromFile (p_rpreconditioner,rparamList,&
                                      spreconditioner,nlevels,RfilterChain)
      end if

      ! Try to get the solver subtype from the parameter list.
      ! This allows switching between right- and left preconditioned BiCGStab.
      call parlst_getvalue_int (p_rsection, 'isolverSubtype', isolverSubtype,0)

      ! Krylow space dimension
      call parlst_getvalue_int (p_rsection, 'ikrylovDim', ikrylowDim,40)
      
      ! Apply Gram Schmidt twice
      call parlst_getvalue_int (p_rsection, 'itwiceGS', i1,0)

      ! Init the solver node
      call linsol_initGMRES (p_rsolverNode,ikrylowDim,p_rpreconditioner,p_Rfilter,&
        i1 .eq. 1)
      
    case (LINSOL_ALG_UMFPACK4)
      ! UMFPACK4 solver
      !
      ! Init the solver node
      call linsol_initUMFPACK4 (p_rsolverNode)
      
    case (LINSOL_ALG_ILU0)
      ! ILU0 solver
      !
      ! Init the solver node
      call linsol_initILU0 (p_rsolverNode)

    case (LINSOL_ALG_MILUS1x1)
      ! (M)ILU solver
      !
      ! Parameters:
      !  ifillinLevel = 0 / 1 / 2 / ...   -> Fill-in level for factorisation
      !  drelax       = 0.0               -> Build ILU(s)
      !               > 0.0               -> Build MILU(s)
      !
      ! Get fill-in level
      call parlst_getvalue_int (p_rsection, 'ifill', i1, 0)
      
      ! Get MILU relaxsation parameter
      call parlst_getvalue_double (p_rsection, 'drelax', d1, 0.0_DP)
      
      ! Init the solver node
      call linsol_initMILUs1x1 (p_rsolverNode,i1,d1)
      
    case (LINSOL_ALG_VANKA)
      ! VANKA solver
      !
      ! Init the solver node
      call linsol_initVANKA (p_rsolverNode,1.0_DP,isolverSubtype)
    
    case (LINSOL_ALG_MULTIGRID)
      ! Multigrid solver
      !
      ! Parameters:
      !  icycle         = 0               -> F-cycle
      !                 = 1               -> V-cycle
      !                 = 2               -> W-cycle
      !  dalphaMin      >= 0.0            -> minimum damping parameter; standard = 1.0
      !  dalphaMin      >= 0.0            -> maximum damping parameter; standard = 1.0
      !  spreSmootherName                 -> Name of the presmoother section
      !  spostSmootherName                -> Name of the postsmoother section
      !  scoarseGridSolver                -> Name of the coarse grid solver section
      !
      ! Ok, that is the most complicated one :-)
      ! At first, initialise the solver:
      call linsol_initMultigrid (p_rsolverNode,p_Rfilter)
      
      ! Then, get solver specific data.
      call parlst_getvalue_int (p_rsection, 'icycle', &
                                p_rsolverNode%p_rsubnodeMultigrid%icycle,&
                                p_rsolverNode%p_rsubnodeMultigrid%icycle)

      ! Coarse grid correction parameters
      call parlst_getvalue_int (p_rsection, 'ccorrectionTypeAlpha', &
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType,&
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType)

      call parlst_getvalue_double (p_rsection, 'dalphaMin', &
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin,&
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin)
      
      call parlst_getvalue_double (p_rsection, 'dalphaMax', &
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax,&
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax)
    
      ! Do we have a presmoother?
      call parlst_getvalue_string (p_rsection, 'spreSmootherName', sString,'')
      spresmoother = ''
      if (sString .ne. '') read (sString,*) spresmoother
      
      ! A postsmoother?
      call parlst_getvalue_string (p_rsection, 'spostSmootherName', sString,'')
      spostsmoother = ''
      if (sString .ne. '') read (sString,*) spostsmoother
      
      ! A coarse grid solver?
      call parlst_getvalue_string (p_rsection, 'scoarseGridSolver', sString,'')
      scoarsegridsolver = ''
      if (sString .ne. '') read (sString,*) scoarsegridsolver
      
      ! We need an interlevel projection structure!
      if (.not. present(rinterlevelProjection)) then
        call output_line ('Cannot create linear solver; ' // &
            'no interlevel projection structure!', OU_CLASS_ERROR, OU_MODE_STD, &
            'linsolinit_initFromFile')
        call sys_halt()
      end if
      
      ! Initialise the coarse grid solver - if there is one. There must be one!
      if (scoarsegridsolver .eq. '') then
        call output_line ('Cannot create linear solver; ' // &
            'no coarse grid solver for MG!', OU_CLASS_ERROR, OU_MODE_STD, &
            'linsolinit_initFromFile')
        call sys_halt()
      end if
      call linsolinit_initFromFile (p_rcoarsegridsolver,rparamList,&
                                    scoarsegridsolver,nlevels,RfilterChain)
      
      ! Build all levels. Level 1 receives the coarse grid solver.
      do ilev = 1,nlevels
        
        ! Where we have a coarse grid solver (level 1), we do not need smoothers.
        if (associated(p_rcoarsegridsolver)) then
          nullify(p_rpresmoother)
          nullify(p_rpostsmoother)
        else
          ! Is there a presmoother?
          if (spresmoother .ne. '') then
            call linsolinit_initFromFile (p_rpresmoother,rparamList,&
                                          spresmoother,nlevels,RfilterChain)
          else
            nullify(p_rpresmoother)
          end if
          
          ! Is there a postsmoother?
          if (spostsmoother .ne. '') then
            ! Check if pre- and postsmoother are identical.
            call sys_toupper (spresmoother)
            call sys_toupper (spostsmoother)
            if (spresmoother .ne. spostsmoother) then
              call linsolinit_initFromFile (p_rpostsmoother,rparamList,&
                                            spostsmoother,nlevels,RfilterChain)
            else
              ! Let the pointer point to the presmoother
              p_rpostsmoother => p_rpresmoother
            end if
          
          else
            nullify(p_rpostsmoother)
          end if
        end if
        
        ! Add the level to Multigrid
        call linsol_addMultigridLevel (p_rlevelInfo,p_rsolverNode, &
                    rinterlevelProjection, &
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver)
                    
        ! Reset the coarse grid solver pointer to NULL.
        ! That way, the coarse grid solver is only attached to level 1.
        nullify(p_rcoarsegridsolver)
        
      end do
    
    case (LINSOL_ALG_MULTIGRID2)
      ! Multigrid solver
      !
      ! Parameters:
      !  icycle         = 0               -> F-cycle
      !                 = 1               -> V-cycle
      !                 = 2               -> W-cycle
      !  dalphaMin      >= 0.0            -> minimum damping parameter; standard = 1.0
      !  dalphaMin      >= 0.0            -> maximum damping parameter; standard = 1.0
      !  spreSmootherName                 -> Name of the presmoother section
      !  spostSmootherName                -> Name of the postsmoother section
      !  scoarseGridSolver                -> Name of the coarse grid solver section

      ! At first, initialise the solver:
      call linsol_initMultigrid2(p_rsolverNode, nlevels, p_Rfilter)

      ! Then, get solver specific data.
      call parlst_getvalue_int (p_rsection, 'icycle', &
                                p_rsolverNode%p_rsubnodeMultigrid2%icycle,&
                                p_rsolverNode%p_rsubnodeMultigrid2%icycle)

      ! Coarse grid correction parameters
      call parlst_getvalue_int (p_rsection, 'ccorrectionTypeAlpha', &
           p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType,&
           p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType)

      call parlst_getvalue_double (p_rsection, 'dalphaMin', &
           p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMin,&
           p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMin)
      
      call parlst_getvalue_double (p_rsection, 'dalphaMax', &
           p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMax,&
           p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMax)
    
      ! Do we have a presmoother?
      call parlst_getvalue_string (p_rsection, 'spreSmootherName', sString,'')
      spresmoother = ''
      if (sString .ne. '') read (sString,*) spresmoother
      
      ! A postsmoother?
      call parlst_getvalue_string (p_rsection, 'spostSmootherName', sString,'')
      spostsmoother = ''
      if (sString .ne. '') read (sString,*) spostsmoother
      
      ! A coarse grid solver?
      call parlst_getvalue_string (p_rsection, 'scoarseGridSolver', sString,'')
      scoarsegridsolver = ''
      if (sString .ne. '') read (sString,*) scoarsegridsolver
      
      ! Initialise the coarse grid solver - if there is one. There must be one!
      if (scoarsegridsolver .eq. '') then
        call output_line ('Cannot create linear solver; ' // &
            'no coarse grid solver for MG!', OU_CLASS_ERROR, OU_MODE_STD, &
            'linsolinit_initFromFile')
        call sys_halt()
      end if

      ! set up a coarse grid solver (always gridlevel 1)
      call linsol_getMultigrid2Level(p_rsolverNode, 1, p_rlevelInfo2)

      call linsolinit_initFromFile(p_rlevelInfo2%p_rcoarseGridSolver, rparamList,&
                                   scoarsegridsolver, nlevels, RfilterChain)
      
      ! Build all the  other levels
      do ilev = 2,nlevels

        ! add this multigrid level
        call linsol_getMultigrid2Level(p_rsolverNode, ilev, p_rlevelInfo2)

        ! Is there a presmoother?
        if (spresmoother .ne. '') then
          call linsolinit_initFromFile(p_rlevelInfo2%p_rpresmoother, rparamList, &
                                       spresmoother, nlevels, RfilterChain)
        else
          nullify(p_rlevelInfo2%p_rpresmoother)
        end if
          
        ! Is there a postsmoother?
        if (spostsmoother .ne. '') then
          ! Check if pre- and postsmoother are identical.
          call sys_toupper (spresmoother)
          call sys_toupper (spostsmoother)
          if (spresmoother .ne. spostsmoother) then
            ! they are not identical, so read in the post smoother
            call linsolinit_initFromFile(p_rlevelInfo2%p_rpostsmoother, rparamList, &
                                         spostsmoother, nlevels, RfilterChain)
          else
            ! otherwise, let the pointer point to the presmoother
            p_rlevelInfo2%p_rpostsmoother => p_rlevelInfo2%p_rpresmoother
          end if
        else
          nullify(p_rlevelInfo2%p_rpostsmoother)
        end if
        
      end do
    
    case default
      call output_line ('Cannot create linear solver; ' // &
            'unsupported solver type isolverType = ' // &
            trim(sys_si(isolverType,8)), &
             OU_CLASS_ERROR, OU_MODE_STD, 'linsolinit_initFromFile')
      call sys_halt()
    
    end select ! isolverType

    ! Ok, we should have a solver node p_rsolverNode now, initialised with
    ! standard parameters. The next task is to change all parameters
    ! in the solver node to those which appear in the given configuration.
    ! If a parameter does not exist in the given configuration, the default
    ! value (that which is already stored in the solver configuration)
    ! must be used.
    !
    ! So parse the given parameters now to initialise the solver node:
    
    call linsolinit_initParams (p_rsolverNode,rparamList,ssolverName)
    
    ! Up to now, we initialised for a linear solver. In case this solver is
    ! used as smoother in a Multigrid algorithm, this initialisation
    ! is not comppletely correct - we have to transform the solver into
    ! one that performs a fixed number of iterations.
    !
    ! We identify a smoother by the existence of the parameter nsmoothingSteps:
    
    if (parlst_queryvalue (p_rsection, 'nsmoothingSteps') .ne. 0) then
    
      call parlst_getvalue_int (p_rsection, 'nsmoothingSteps', &
                                i1, p_rsolverNode%nminIterations)
                                
      ! Convert the solver node into a smoother:
      call linsol_convertToSmoother (p_rsolverNode,i1)
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine linsolinit_initParams (rsolverNode,rparamList,ssection,csolverType)
  
!<description>
  ! This section reads the standard parameters for a linear solver from the
  ! section ssection in the parameter list rparamList and changes
  ! specified entries in the solver node rsolverNode. Variables whose
  ! values are not specified in the parameter list are left unchanged
  ! in rsolverNode.
  !
  ! The parameters that can be specified in the parameter list have exactly
  ! the same names as in the linear solver structure rsolverNode.
  ! Note that not all parameters in rsolverNode are changed by this routine.
  ! 'Critical' parameters like solver type etc. are left unchanged.
  !
  ! csolverType allows to specify a special solver type whose parameters should
  ! be initialised. If this is desired, the routine should be called twice --
  ! once with csolverType=LINSOL_ALG_UNDEFINED to initialise the main parameters
  ! and once with csolverType=LINSOL_ALG_xxxx to initialise the special,
  ! solver dependent parameters.
!</description>

!<input>
  ! The parameter list that contains the whole solver configuration.
  type(t_parlist), intent(in) :: rparamList
  
  ! The name of a section in rparamList that contains the configuration of
  ! the linear solver.
  character(LEN=*), intent(in) :: ssection
  
  ! OPTIONAL: Type of solver structure that should be initialised.
  ! If unspecified or if set to LINSOL_ALG_UNDEFINED, the parameters in the
  ! main solver structure are initialised.
  ! Otherwise, csolverType must be a LINSOL_ALG_xxx constant that specifies
  ! a special solver type whose parameters should be initialised (e.g.
  ! the 'ikrylovDim' parameter of the GMRES method).
  integer, intent(in), optional :: csolverType
!</input>

!<output>
  ! A linear solver node whose parameters should be changed according to the
  ! parameters in rparamList.
  type(t_linsolNode), intent(inout) :: rsolverNode
!</output>

!</subroutine>

    ! local variables
    integer :: csolver
    type(t_parlstSection), pointer :: p_rsection
    integer :: i1,ikrylowDim
    character(LEN=SYS_STRLEN) :: sstring

    csolver = LINSOL_ALG_UNDEFINED
    if (present(csolverType)) csolver = csolverType

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, ssection, p_rsection)
    
    if (.not. associated(p_rsection)) then
      call output_line ('Cannot create linear solver; no section ''' &
          // trim(ssection) // '''!', OU_CLASS_ERROR, OU_MODE_STD, &
          'linsolinit_initParams')
      call sys_halt()
    end if
    
    ! Now comes a biiig select for all the different types of solvers
    ! that are supported by this routine.
    select case (csolver)
      
    case (LINSOL_ALG_SSOR)
      ! SSOR solver
      !
      ! Init the solver node.
      ! Parameters:
      !  iscale = 1   -> Scale preconditioned vector according to actual formula
      !                  in the literature
      !         = 0   -> no scaling, original FEAT implementation
      call parlst_getvalue_int (p_rsection, 'iscale', i1, -1)
      if (i1 .ne. -1) then
        rsolverNode%p_rsubnodeSSOR%bscale = i1 .eq. 1
      end if
      
    case (LINSOL_ALG_GMRES)
      ! GMRES solver

      ! Krylow space dimension
      call parlst_getvalue_int (p_rsection, 'ikrylovDim', ikrylowDim,40)
      
      ! Apply Gram Schmidt twice
      call parlst_getvalue_int (p_rsection, 'btwiceGS', i1,-1)

      rsolverNode%p_rsubnodeGMRES%ikrylovDim = ikrylowDim
      if (i1 .ne. -1) then
        rsolverNode%p_rsubnodeGMRES%btwiceGS = i1 .eq. 1
      end if
      
    case (LINSOL_ALG_MILUS1x1)
      ! (M)ILU solver
      !
      ! Parameters:
      !  ifillinLevel = 0 / 1 / 2 / ...   -> Fill-in level for factorisation
      !  drelax       = 0.0               -> Build ILU(s)
      !               > 0.0               -> Build MILU(s)
      !
      ! Get fill-in level
      call parlst_getvalue_int (p_rsection, 'ifill', &
        rsolverNode%p_rsubnodeMILUs1x1%ifill, &
        rsolverNode%p_rsubnodeMILUs1x1%ifill)
      
      ! Get MILU relaxsation parameter
      call parlst_getvalue_double (p_rsection, 'drelax', &
        rsolverNode%p_rsubnodeMILUs1x1%drelax, &
        rsolverNode%p_rsubnodeMILUs1x1%drelax)
      
    case (LINSOL_ALG_MULTIGRID)
      ! Multigrid solver
      !
      ! Parameters:
      !  icycle         = 0               -> F-cycle
      !                 = 1               -> V-cycle
      !                 = 2               -> W-cycle
      !  dalphaMin      >= 0.0            -> minimum damping parameter; standard = 1.0
      !  dalphaMin      >= 0.0            -> maximum damping parameter; standard = 1.0
      
      ! Then, get solver specific data.
      call parlst_getvalue_int (p_rsection, 'icycle', &
                                rsolverNode%p_rsubnodeMultigrid%icycle,&
                                rsolverNode%p_rsubnodeMultigrid%icycle)

      ! Coarse grid correction parameters
      call parlst_getvalue_int (p_rsection, 'ccorrectionTypeAlpha', &
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType,&
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType)

      call parlst_getvalue_double (p_rsection, 'dalphaMin', &
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin,&
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin)
      
      call parlst_getvalue_double (p_rsection, 'dalphaMax', &
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax,&
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax)
    
    case (LINSOL_ALG_MULTIGRID2)
      ! Multigrid solver
      !
      ! Parameters:
      !  icycle         = 0               -> F-cycle
      !                 = 1               -> V-cycle
      !                 = 2               -> W-cycle
      !  dalphaMin      >= 0.0            -> minimum damping parameter; standard = 1.0
      !  dalphaMin      >= 0.0            -> maximum damping parameter; standard = 1.0
      
      ! Then, get solver specific data.
      call parlst_getvalue_int (p_rsection, 'icycle', &
                                rsolverNode%p_rsubnodeMultigrid2%icycle,&
                                rsolverNode%p_rsubnodeMultigrid2%icycle)

      ! Coarse grid correction parameters
      call parlst_getvalue_int (p_rsection, 'ccorrectionTypeAlpha', &
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType,&
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType)

      call parlst_getvalue_double (p_rsection, 'dalphaMin', &
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMin,&
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMin)
      
      call parlst_getvalue_double (p_rsection, 'dalphaMax', &
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMax,&
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMax)
    
    case default
    
      ! Initialise the main solver parameters
    
      call parlst_getvalue_double (p_rsection, 'domega', &
                                  rsolverNode%domega, rsolverNode%domega)

      call parlst_getvalue_int (p_rsection, 'nminIterations', &
                                rsolverNode%nminIterations, rsolverNode%nminIterations)

      call parlst_getvalue_int (p_rsection, 'nmaxIterations', &
                                rsolverNode%nmaxIterations, rsolverNode%nmaxIterations)

      call parlst_getvalue_string (p_rsection, 'iresCheck', sstring, 'YES')
      call sys_toUpper(sstring)
      if (sstring .eq. 'NO') then
        ! if the string is exactly 'NO' then disable checking of residual norm
        rsolverNode%iresCheck = NO
      else
        ! in all other cases enable checking of residual norm
        rsolverNode%iresCheck = YES
      endif
  
      call parlst_getvalue_double (p_rsection, 'depsRel', &
                                  rsolverNode%depsRel, rsolverNode%depsRel)

      call parlst_getvalue_double (p_rsection, 'depsAbs', &
                                  rsolverNode%depsAbs, rsolverNode%depsAbs)

      call parlst_getvalue_double (p_rsection, 'depsDiff', &
                                  rsolverNode%depsDiff, rsolverNode%depsDiff)

      call parlst_getvalue_int (p_rsection, 'iresNorm', &
                                rsolverNode%iresNorm, rsolverNode%iresNorm)

      call parlst_getvalue_int (p_rsection, 'niteAsymptoticCVR', &
                                rsolverNode%niteAsymptoticCVR, &
                                rsolverNode%niteAsymptoticCVR)

      call parlst_getvalue_int (p_rsection, 'istoppingCriterion', &
                                rsolverNode%istoppingCriterion, &
                                rsolverNode%istoppingCriterion)

      call parlst_getvalue_int (p_rsection, 'ioutputLevel', &
                                rsolverNode%ioutputLevel, rsolverNode%ioutputLevel)

      call parlst_getvalue_int (p_rsection, 'niteResOutput', &
                                rsolverNode%niteResOutput, rsolverNode%niteResOutput)
                                
      call parlst_getvalue_int (p_rsection, 'isolverSubgroup', &
                                rsolverNode%isolverSubgroup, &
                                rsolverNode%isolverSubgroup)

    end select ! isolverType

  end subroutine
  
end module
