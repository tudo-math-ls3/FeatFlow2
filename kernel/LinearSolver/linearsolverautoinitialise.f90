!##############################################################################
!# ****************************************************************************
!# <name> LinearSolverAutoInitialise </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to initialise a linear solver for a given
!# situation. On one hand, it provides initialisation routines for standard
!# solver settings. On the other hand, there's a parser included which can
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

MODULE linearsolverautoinitialise

  USE fsystem
  USE linearsolver
  USE paramlist

  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsolinit_CC2DMultigrid(p_rsolverNode, isolverType, nlevels, &
                        niterationsMin, niterationsMax,daccuracyAbs,daccuracyRel,  &
                        nprecSteps, nsmoothingSteps, &
                        nmaxCoarseGridSteps, dcoarseGridAccuracyAbs, &
                        dcoarseGridAccuracyRel)
  
  !<description>
  
  ! This routine builds a solver node for the linear solver for CC2D-like
  ! applications. The following combinations are supported:
  ! 1.) Multigrid solver, VANCA smoother, VANCA coarse grid solver
  ! 2.) BiCGStab-solver, Multigrid preconditioner, 
  !     VANCA smoother, VANCA coarse grid solver
  ! The caller must attach the matrices of all levels to the solver manually
  ! by calling linsol_setMatrices. Afterwards the linear solver can be started.
  
  !</description>
  
  !<input>
  
  ! Type of solver structure which should be build up.
  ! 1 = Multigrid solver, VANCA smoother, VANCA coarse grid solver
  ! 2 = BiCGStab-solver, Multigrid preconditioner, VANCA smoother,
  !     VANCA coarse grid solver
  INTEGER, INTENT(IN)               :: isolverType

  ! Number of levels
  INTEGER, INTENT(IN)               :: nlevels

  ! Minimum number of solver iterations
  INTEGER, INTENT(IN)               :: niterationsMin

  ! Maximum number of solver iterations
  INTEGER, INTENT(IN)               :: niterationsMax

  ! If solver conbination is BiCGStab with MG preconditioning, number of
  ! steps multigrid should perform. Otherwise ignored
  INTEGER, INTENT(IN)               :: nprecSteps

  ! Absolute accuracy of the solver
  REAL(DP), INTENT(IN)                   :: daccuracyAbs

  ! Relativew accuracy of the solver
  REAL(DP), INTENT(IN)                   :: daccuracyRel

  ! Number of pre- and postsmoothing steps
  INTEGER, INTENT(IN)               :: nsmoothingSteps
  
  ! Absolute accuracy on the coarse grid
  REAL(DP), INTENT(IN)                   :: dcoarseGridAccuracyAbs

  ! Relative accuracy on the coarse grid
  REAL(DP), INTENT(IN)                   :: dcoarseGridAccuracyRel
  
  ! Maximum iterations on the coarse grid
  INTEGER, INTENT(IN)               :: nmaxCoarseGridSteps

  ! A list of spatial discretisation structures for the three equations
  ! that are supported by CC2D: x-velocity, y-velocity, pressure
  TYPE(t_spatialDiscretisation), DIMENSION(3) :: RspatialDiscr
  
!</input>
  
!<output>
  
  ! The solver node that identifies the solver. The caller must attach
  ! level-dependent data (matrix information) to it with the standard
  ! routines in the module LinearSolver. Afterwards the structure
  ! can be used to solve the problem.
  
  TYPE(t_linsolNode),POINTER :: p_rsolverNode
  
!</output>
  
!</subroutine>

  ! local variables
  TYPE(t_linsolNode),POINTER :: p_rmgSolver
  INTEGER               :: ilevel
  TYPE(t_linsolNode),POINTER :: p_rpreSmoother
  TYPE(t_linsolNode),POINTER :: p_rpostSmoother
  TYPE(t_linsolNode),POINTER :: p_rcoarseGridSolver
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
  TYPE(t_interlevelProjectionBlock) :: rprojection 
  
  ! Create the solver node - either BiCGStab or Multigrid.
  ! If we create BiCGStab, attach multigrid as preconditioner.
  
  CALL linsol_initMultigrid (p_rmgSolver)
  
  IF (isolverType .EQ. 1) THEN
  
    p_rsolverNode => p_rmgSolver
    
  ELSE
  
    CALL linsol_initBiCGStab (p_rsolverNode,p_rmgSolver)
    
    ! Configure MG preconditioning as a fixed number of MG steps
    ! without checking the residuum.
    
    p_rmgSolver%nminIterations = nprecSteps
    p_rmgSolver%nmaxIterations = nprecSteps
    p_rmgSolver%iresCheck      = NO
    
  END IF
  
  ! Configure the solver
  p_rsolverNode%nminIterations = niterationsMin
  p_rsolverNode%nmaxIterations = niterationsMax
  p_rsolverNode%depsRel = daccuracyRel
  p_rsolverNode%depsAbs = daccuracyAbs
  
  ! Initialise a standard interlevel projection structure for all levels.
  CALL mlprj_initProjectionDirect (rprojection,RspatialDiscr)
  
  ! Continue to configure MG by accessing p_rmgSolver.
  ! Loop through the levels.
  
  DO ilevel = 1,nlevels
  
    ! On the lowest level create a coarse grid solver structure
    IF (ilevel .EQ. 1) THEN
      CALL linsol_initVANCA (p_rcoarseGridSolver)
      p_rcoarseGridSolver%depsRel = dcoarseGridAccuracyRel
      p_rcoarseGridSolver%depsAbs = dcoarseGridAccuracyAbs
      p_rcoarseGridSolver%nmaxIterations = nmaxCoarseGridSteps
    ELSE
      NULLIFY(p_rcoarseGridSolver)
    END IF
    
    ! Create pre- and postsmoother structure on the current level
    CALL linsol_initVANCA (p_rpreSmoother)
    CALL linsol_initVANCA (p_rpostSmoother)
    
    ! Configure the structures to form a smoother. A smoother is a solver
    ! that iterates a finite number of times without respecting the
    ! residuum.
    
    CALL linsol_convertToSmoother (p_rpreSmoother,nsmoothingSteps)
    CALL linsol_convertToSmoother (p_rpostSmoother,nsmoothingSteps)
    
    ! Create the level, attrach it to the solver and proceed to the 
    ! next level
    CALL linsol_addMultigridLevel (p_rlevelInfo,p_rmgSolver, rprojection,&
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver)
    
  END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE linsolinit_initFromFile (p_rsolverNode,&
                                                rparamList,ssolverName,&
                                                nlevels, &
                                                RfilterChain,rinterlevelProjection)
  
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
  TYPE(t_parlist), INTENT(IN) :: rparamList
  
  ! The name of a section in rparamList that contains the configuration of
  ! the main linear solver.
  CHARACTER(LEN=*), INTENT(IN) :: ssolverName
  
  ! Number of levels in the discretisation.
  INTEGER, INTENT(IN) :: nlevels

  ! OPTIONAL: A filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the 
  ! iteration. If not, no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  TYPE(t_filterChain), DIMENSION(:), INTENT(IN), TARGET, OPTIONAL :: RfilterChain

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
  TYPE(t_interlevelProjectionBlock), OPTIONAL :: rinterlevelProjection
!</input>

!<output>
  ! A pointer to a new linear solver node on the heap, initialised
  ! by the configuration in the parameter list.
  TYPE(t_linsolNode), POINTER :: p_rsolverNode
!</output>

!</subroutine>

    ! local variables
    TYPE(t_parlstSection), POINTER :: p_rsection
    CHARACTER(LEN=PARLST_MLDATA) :: spreconditioner,spresmoother, spostsmoother
    CHARACTER(LEN=PARLST_MLDATA) :: scoarsegridsolver,sString
    INTEGER :: isolverType,isolverSubtype
    TYPE(t_linsolNode), POINTER :: p_rpreconditioner,p_rpresmoother,p_rpostsmoother
    TYPE(t_linsolNode), POINTER :: p_rcoarsegridsolver
    INTEGER :: i1,ilev,ikrylowDim
    REAL(DP) :: d1
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_Rfilter
    TYPE(t_linsolMGLevelInfo), POINTER     :: p_rlevelInfo

    NULLIFY(p_Rfilter)
    IF (PRESENT(RfilterChain)) THEN
      p_Rfilter => RfilterChain
    END IF

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    CALL parlst_querysection(rparamList, ssolverName, p_rsection) 
    
    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      PRINT *,'Cannot create linear solver; no section '''&
              //TRIM(ssolverName)//'''!'
      CALL sys_halt()
    END IF
    
    ! Ok, we have the section where we can get parameters from.
    ! Get the solver type we should set up.
    ! Let's hope the parameter 'isolverType' exists! That one is
    ! mandatory; if not, the get-routine will stop.
    CALL parlst_getvalue_int (p_rsection, 'isolverType', isolverType)

    ! Try to get the solver subtype from the parameter list.
    ! This allows switching between different variants of the same
    ! basic algorithm (e.g. VANCA)
    CALL parlst_getvalue_int (p_rsection, 'isolverSubtype', isolverSubtype,0)

    ! Many solvers support preconditioners - we try to fetch the
    ! name of the preconditioner in-advance to prevent code
    ! duplication
    CALL parlst_getvalue_string (p_rsection, 'spreconditioner', sString,'')
    spreconditioner = ''
    IF (sString .NE. '') READ (sString,*) spreconditioner
    
    ! Now comes a biiig select for all the different types of solvers
    ! that are supported by this routine.
    SELECT CASE (isolverType)
    
    CASE (LINSOL_ALG_DEFCORR)
      ! Defect correction
      !
      ! Initialise a solver node for the preconditioner - if there is one.
      NULLIFY(p_rpreconditioner)
      IF (spreconditioner .NE. '') THEN
        CALL linsolinit_initFromFile (p_rpreconditioner,rparamList,&
                                      spreconditioner,nlevels,RfilterChain)
      END IF
      ! Init the solver node
      CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,p_Rfilter)
      
    CASE (LINSOL_ALG_JACOBI)
      ! Jacobi solver
      !
      ! Init the solver node
      CALL linsol_initJacobi (p_rsolverNode)
      
    CASE (LINSOL_ALG_SOR)
      ! SOR/GS solver
      !
      ! Init the solver node. domega is set to 1.0 as standard
      ! to activate GS. The actual domega is initialised later if specified
      ! in the DAT file.
      CALL linsol_initSOR (p_rsolverNode,1.0_DP)
      
    CASE (LINSOL_ALG_SSOR)
      ! SSOR solver
      !
      ! Init the solver node.
      ! Parameters:
      !  iscale = 1   -> Scale preconditioned vector according to actual formula
      !                  in the literature
      !         = 0   -> no scaling, original FEAT implementation
      CALL parlst_getvalue_int (p_rsection, 'bscale', i1, 0)
      CALL linsol_initSSOR (p_rsolverNode,bscale = i1 .EQ. 1)
      
    CASE (LINSOL_ALG_BICGSTAB)
      ! BiCGStab solver
      !
      ! Initialise a solver node for the preconditioner - if there is one.
      NULLIFY(p_rpreconditioner)
      IF (spreconditioner .NE. '') THEN
        CALL linsolinit_initFromFile (p_rpreconditioner,rparamList,&
                                      spreconditioner,nlevels,RfilterChain)
      END IF

      ! Try to get the solver subtype from the parameter list.
      ! This allows switching between right- and left preconditioned BiCGStab.
      CALL parlst_getvalue_int (p_rsection, 'iSolverSubtype', isolverSubtype,0)

      ! Init the solver node
      CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_Rfilter,&
        isolverSubtype)
      
    CASE (LINSOL_ALG_GMRES)
      ! GMRES solver
      !
      ! Initialise a solver node for the preconditioner - if there is one.
      NULLIFY(p_rpreconditioner)
      IF (spreconditioner .NE. '') THEN
        CALL linsolinit_initFromFile (p_rpreconditioner,rparamList,&
                                      spreconditioner,nlevels,RfilterChain)
      END IF

      ! Try to get the solver subtype from the parameter list.
      ! This allows switching between right- and left preconditioned BiCGStab.
      CALL parlst_getvalue_int (p_rsection, 'isolverSubtype', isolverSubtype,0)

      ! Krylow space dimension
      CALL parlst_getvalue_int (p_rsection, 'ikrylovDim', ikrylowDim,40)
      
      ! Apply Gram Schmidt twice
      CALL parlst_getvalue_int (p_rsection, 'itwiceGS', i1,0)

      ! Init the solver node
      CALL linsol_initGMRES (p_rsolverNode,ikrylowDim,p_rpreconditioner,p_Rfilter,&
        i1 .EQ. 1)
      
    CASE (LINSOL_ALG_UMFPACK4)
      ! UMFPACK4 solver
      !
      ! Init the solver node
      CALL linsol_initUMFPACK4 (p_rsolverNode)
      
    CASE (LINSOL_ALG_MILUS1x1)
      ! (M)ILU solver
      !
      ! Parameters:
      !  ifillinLevel = 0 / 1 / 2 / ...   -> Fill-in level for factorisation
      !  drelax       = 0.0               -> Build ILU(s)
      !               > 0.0               -> Build MILU(s)
      !
      ! Get fill-in level
      CALL parlst_getvalue_int (p_rsection, 'ifill', i1, 0)
      
      ! Get MILU relaxsation parameter
      CALL parlst_getvalue_double (p_rsection, 'drelax', d1, 0.0_DP)
      
      ! Init the solver node
      CALL linsol_initMILUs1x1 (p_rsolverNode,i1,d1)
      
    CASE (LINSOL_ALG_VANCA)
      ! VANCA solver
      !
      ! Init the solver node
      CALL linsol_initVANCA (p_rsolverNode,1.0_DP,isolverSubtype)
    
    CASE (LINSOL_ALG_MULTIGRID)
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
      ! Ok, that's the most complicated one :-)
      ! At first, initialise the solver:
      CALL linsol_initMultigrid (p_rsolverNode,p_Rfilter)
      
      ! Then, get solver specific data.
      CALL parlst_getvalue_int (p_rsection, 'icycle', &
                                p_rsolverNode%p_rsubnodeMultigrid%icycle,&
                                p_rsolverNode%p_rsubnodeMultigrid%icycle)

      ! Coarse grid correction parameters
      CALL parlst_getvalue_int (p_rsection, 'ccorrectionTypeAlpha', &
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType,&
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType)

      CALL parlst_getvalue_double (p_rsection, 'dalphaMin', &
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin,&
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin)
      
      CALL parlst_getvalue_double (p_rsection, 'dalphaMax', &
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax,&
           p_rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax)
    
      ! Do we have a presmoother?
      CALL parlst_getvalue_string (p_rsection, 'spreSmootherName', sString,'')
      spresmoother = ''
      IF (sString .NE. '') READ (sString,*) spresmoother
      
      ! A postsmoother?
      CALL parlst_getvalue_string (p_rsection, 'spostSmootherName', sString,'')
      spostsmoother = ''
      IF (sString .NE. '') READ (sString,*) spostsmoother
      
      ! A coarse grid solver?
      CALL parlst_getvalue_string (p_rsection, 'scoarseGridSolver', sString,'')
      scoarsegridsolver = ''
      IF (sString .NE. '') READ (sString,*) scoarsegridsolver
      
      ! We need an interlevel projection structure!
      IF (.NOT. PRESENT(rinterlevelProjection)) THEN
        PRINT *,'Cannot create linear solver; no interlevel projection structure!'
        CALL sys_halt()
      END IF
      
      ! Initialise the coarse grid solver - if there is one. There must be one!
      IF (scoarsegridsolver .EQ. '') THEN
        PRINT *,'Cannot create linear solver; no coarse grid solver for MG!'
        CALL sys_halt()
      END IF
      CALL linsolinit_initFromFile (p_rcoarsegridsolver,rparamList,&
                                    scoarsegridsolver,nlevels,RfilterChain)
      
      ! Build all levels. Level 1 receives the coarse grid solver.
      DO ilev = 1,nlevels
        
        ! Where we have a coarse grid solver (level 1), we don't need smoothers.
        IF (ASSOCIATED(p_rcoarsegridsolver)) THEN
          NULLIFY(p_rpresmoother)
          NULLIFY(p_rpostsmoother)
        ELSE
          ! Is there a presmoother?
          IF (spresmoother .NE. '') THEN
            CALL linsolinit_initFromFile (p_rpresmoother,rparamList,&
                                          spresmoother,nlevels,RfilterChain)
          ELSE
            NULLIFY(p_rpresmoother)
          END IF
          
          ! Is there a postsmoother?
          IF (spostsmoother .NE. '') THEN
            ! Check if pre- and postsmoother are identical.
            CALL sys_toupper (spresmoother) 
            CALL sys_toupper (spostsmoother) 
            IF (spresmoother .NE. spostsmoother) THEN
              CALL linsolinit_initFromFile (p_rpostsmoother,rparamList,&
                                            spostsmoother,nlevels,RfilterChain)
            ELSE 
              ! Let the pointer point to the presmoother
              p_rpostsmoother => p_rpresmoother
            END IF
          
          ELSE
            NULLIFY(p_rpostsmoother)
          END IF
        END IF
        
        ! Add the level to Multigrid
        CALL linsol_addMultigridLevel (p_rlevelInfo,p_rsolverNode, &
                    rinterlevelProjection, &
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver)
                    
        ! Reset the coarse grid solver pointer to NULL.
        ! That way, the coarse grid solver is only attached to level 1.
        NULLIFY(p_rcoarsegridsolver)
        
      END DO
    
    CASE DEFAULT
      PRINT *,'Cannot create linear solver; unsupported solver type isolverType=',&
              isolverType
      CALL sys_halt()
    
    END SELECT ! isolverType

    ! Ok, we should have a solver node p_rsolverNode now, initialised with
    ! standard parameters. The next task is to change all parameters
    ! in the solver node to those which appeat in the given configuration.
    ! If a parameter does not exist in the given configuration, the default
    ! value (that which is already stored in the solver configuration)
    ! must be used.
    !
    ! So parse the given parameters now to initialise the solver node:
    
    CALL linsolinit_initParams (p_rsolverNode,rparamList,ssolverName)
    
    ! Up to now, we initialised for a linear solver. In case this solver is
    ! used as smoother in a Multigrid algorithm, this initialisation
    ! is not comppletely correct - we have to transform the solver into
    ! one that performs a fixed number of iterations.
    !
    ! We identify a smoother by the existence of the parameter nsmoothingSteps:
    
    IF (parlst_queryvalue (p_rsection, 'nsmoothingSteps') .NE. 0) THEN
    
      CALL parlst_getvalue_int (p_rsection, 'nsmoothingSteps', &
                                i1, p_rsolverNode%nminIterations)
                                
      ! Convert the solver node into a smoother:
      CALL linsol_convertToSmoother (p_rsolverNode,i1)
    
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE linsolinit_initParams (rsolverNode,rparamList,ssection,csolverType)
  
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
  TYPE(t_parlist), INTENT(IN) :: rparamList
  
  ! The name of a section in rparamList that contains the configuration of
  ! the linear solver.
  CHARACTER(LEN=*), INTENT(IN) :: ssection
  
  ! OPTIONAL: Type of solver structure that should be initialised.
  ! If unspecified or if set to LINSOL_ALG_UNDEFINED, the parameters in the
  ! main solver structure are initialised.
  ! Otherwise, csolverType must be a LINSOL_ALG_xxx constant that specifies
  ! a special solver type whose parameters should be initialised (e.g.
  ! the 'ikrylovDim' parameter of the GMRES method).
  INTEGER, INTENT(IN), OPTIONAL :: csolverType
!</input>

!<output>
  ! A linear solver node whose parameters should be changed according to the
  ! parameters in rparamList.
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
!</output>

!</subroutine>

    ! local variables
    INTEGER :: csolver
    TYPE(t_parlstSection), POINTER :: p_rsection
    INTEGER :: i1,ikrylowDim
    
    csolver = LINSOL_ALG_UNDEFINED
    IF (PRESENT(csolverType)) csolver = csolverType

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    CALL parlst_querysection(rparamList, ssection, p_rsection) 
    
    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      PRINT *,'Cannot create linear solver; no section '''&
              //TRIM(ssection)//'''!'
      CALL sys_halt()
    END IF
    
    ! Now comes a biiig select for all the different types of solvers
    ! that are supported by this routine.
    SELECT CASE (csolver)
      
    CASE (LINSOL_ALG_SSOR)
      ! SSOR solver
      !
      ! Init the solver node.
      ! Parameters:
      !  iscale = 1   -> Scale preconditioned vector according to actual formula
      !                  in the literature
      !         = 0   -> no scaling, original FEAT implementation
      CALL parlst_getvalue_int (p_rsection, 'iscale', i1, -1)
      IF (i1 .NE. -1) THEN
        rsolverNode%p_rsubnodeSSOR%bscale = i1 .EQ. 1
      END IF
      
    CASE (LINSOL_ALG_BICGSTAB)
      ! BiCGStab solver

      ! Try to get the solver subtype from the parameter list.
      ! This allows switching between right- and left preconditioned BiCGStab.
      CALL parlst_getvalue_int (p_rsection, 'isolverSubtype', &
        rsolverNode%p_rsubnodeBiCGStab%cprecondType,&
        rsolverNode%p_rsubnodeBiCGStab%cprecondType)

    CASE (LINSOL_ALG_GMRES)
      ! GMRES solver

      ! Krylow space dimension
      CALL parlst_getvalue_int (p_rsection, 'ikrylovDim', ikrylowDim,40)
      
      ! Apply Gram Schmidt twice
      CALL parlst_getvalue_int (p_rsection, 'btwiceGS', i1,-1)

      rsolverNode%p_rsubnodeGMRES%ikrylovDim = ikrylowDim
      IF (i1 .NE. -1) THEN
        rsolverNode%p_rsubnodeGMRES%btwiceGS = i1 .EQ. 1
      END IF
      
    CASE (LINSOL_ALG_MILUS1x1)
      ! (M)ILU solver
      !
      ! Parameters:
      !  ifillinLevel = 0 / 1 / 2 / ...   -> Fill-in level for factorisation
      !  drelax       = 0.0               -> Build ILU(s)
      !               > 0.0               -> Build MILU(s)
      !
      ! Get fill-in level
      CALL parlst_getvalue_int (p_rsection, 'ifill', &
        rsolverNode%p_rsubnodeMILUs1x1%ifill, &
        rsolverNode%p_rsubnodeMILUs1x1%ifill)
      
      ! Get MILU relaxsation parameter
      CALL parlst_getvalue_double (p_rsection, 'drelax', &
        rsolverNode%p_rsubnodeMILUs1x1%drelax, &
        rsolverNode%p_rsubnodeMILUs1x1%drelax)
      
    CASE (LINSOL_ALG_MULTIGRID)
      ! Multigrid solver
      !
      ! Parameters:
      !  icycle         = 0               -> F-cycle
      !                 = 1               -> V-cycle
      !                 = 2               -> W-cycle
      !  dalphaMin      >= 0.0            -> minimum damping parameter; standard = 1.0
      !  dalphaMin      >= 0.0            -> maximum damping parameter; standard = 1.0
      
      ! Then, get solver specific data.
      CALL parlst_getvalue_int (p_rsection, 'icycle', &
                                rsolverNode%p_rsubnodeMultigrid%icycle,&
                                rsolverNode%p_rsubnodeMultigrid%icycle)

      ! Coarse grid correction parameters
      CALL parlst_getvalue_int (p_rsection, 'ccorrectionTypeAlpha', &
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType,&
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%ccorrectionType)

      CALL parlst_getvalue_double (p_rsection, 'dalphaMin', &
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin,&
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMin)
      
      CALL parlst_getvalue_double (p_rsection, 'dalphaMax', &
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax,&
           rsolverNode%p_rsubnodeMultigrid%rcoarseGridCorrection%dalphaMax)
    
    CASE (LINSOL_ALG_MULTIGRID2)
      ! Multigrid solver
      !
      ! Parameters:
      !  icycle         = 0               -> F-cycle
      !                 = 1               -> V-cycle
      !                 = 2               -> W-cycle
      !  dalphaMin      >= 0.0            -> minimum damping parameter; standard = 1.0
      !  dalphaMin      >= 0.0            -> maximum damping parameter; standard = 1.0
      
      ! Then, get solver specific data.
      CALL parlst_getvalue_int (p_rsection, 'icycle', &
                                rsolverNode%p_rsubnodeMultigrid2%icycle,&
                                rsolverNode%p_rsubnodeMultigrid2%icycle)

      ! Coarse grid correction parameters
      CALL parlst_getvalue_int (p_rsection, 'ccorrectionTypeAlpha', &
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType,&
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType)

      CALL parlst_getvalue_double (p_rsection, 'dalphaMin', &
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMin,&
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMin)
      
      CALL parlst_getvalue_double (p_rsection, 'dalphaMax', &
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMax,&
           rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%dalphaMax)
    
    CASE DEFAULT
    
      ! Initialise the main solver parameters
    
      CALL parlst_getvalue_double (p_rsection, 'domega', &
                                  rsolverNode%domega, rsolverNode%domega)

      CALL parlst_getvalue_int (p_rsection, 'nminIterations', &
                                rsolverNode%nminIterations, rsolverNode%nminIterations)

      CALL parlst_getvalue_int (p_rsection, 'nmaxIterations', &
                                rsolverNode%nmaxIterations, rsolverNode%nmaxIterations)

      CALL parlst_getvalue_double (p_rsection, 'depsRel', &
                                  rsolverNode%depsRel, rsolverNode%depsRel)

      CALL parlst_getvalue_double (p_rsection, 'depsAbs', &
                                  rsolverNode%depsAbs, rsolverNode%depsAbs)

      CALL parlst_getvalue_int (p_rsection, 'iresNorm', &
                                rsolverNode%iresNorm, rsolverNode%iresNorm)

      CALL parlst_getvalue_int (p_rsection, 'istoppingCriterion', &
                                rsolverNode%istoppingCriterion, &
                                rsolverNode%istoppingCriterion)

      CALL parlst_getvalue_int (p_rsection, 'ioutputLevel', &
                                rsolverNode%ioutputLevel, rsolverNode%ioutputLevel)

      CALL parlst_getvalue_int (p_rsection, 'niteResOutput', &
                                rsolverNode%niteResOutput, rsolverNode%niteResOutput)
                                
      CALL parlst_getvalue_int (p_rsection, 'isolverSubgroup', &
                                rsolverNode%isolverSubgroup, &
                                rsolverNode%isolverSubgroup)

    END SELECT ! isolverType

  END SUBROUTINE
  
END MODULE
