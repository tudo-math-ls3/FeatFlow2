!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2nonlinearcoreinit </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains initialisation routines for the core equation
!# (see also cc2dminim2nonlinearcore).
!# These routines connect the "problem" structure with the "core equation"
!# structure. In detail, we have routines that initialise the preconditioner
!# and all the information structures that are used during the nonlinear
!# iteration.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_allocSystemMatrix
!#     -> Allocates memory for the system matrix representing the
!#        core equation.
!#
!# 2.) c2d2_initNonlinearLoop
!#     -> Initialises a 'nonlinear iteration structure' with parameters from
!#        the DAT file. This is needed for solving the core equation.
!#     -> Extension to c2d2_createNonlinearLoop.
!#
!# 3.) c2d2_doneNonlinearLoop
!#     -> Cleans up a 'nonlinear iteration structure' initialised by
!#        c2d2_initNonlinearLoop.
!#     -> Extension to c2d2_releaseNonlinearLoop.
!#
!# 4.) c2d2_initPreconditioner
!#     -> Prepare preconditioner of nonlinear iteration
!#
!# 5.) c2d2_updatePreconditioner
!#     -> Updates the preconditioner if there was a change in the system matrices
!#
!# 6.) c2d2_releasePreconditioner
!#     -> Clean up preconditioner of nonlinear iteration
!#
!# 7.) c2d2_getProlRest
!#     -> Auxiliary routine: Set up interlevel projection structure
!#        with information from INI/DAT files
!#
!# Auxiliary routines, not to be called from outside:
!#
!# 1.) c2d2_checkAssembly
!#     -> Checks if the system matrices are compatible to the preconditioner.
!#        Set up some situation dependent 'tweak' flags for the assembly
!#        of the matrices.
!#
!# 2.) c2d2_finaliseMatrices
!#     -> Rearranges the structure of the preconditioner matrices if necessary.
!#
!# 3.) c2d2_unfinaliseMatrices
!#     -> Reverts the changes of c2d2_finaliseMatrices and brings preconditioner
!#        matrices into their original form.
!#
!# The module works in tight relationship to cc2dmediumm2nonlinearcore.
!# cc2dmediumm2nonlinearcodeinit provides the routines to initialise
!# preconditioner and important structures using the problem related
!# structure. This module cc2dmediumm2nonlinearcore on the other hand
!# contains only the 'main' worker routines that do the work of the
!# nonlinear iteration -- independent of the problem structure!
!#
!# Note that this module and the "nonlinearcore" module are the only modules
!# that 'know' the actual structure of the system matrix and how to link
!# it to the main problem! For the actual assembly of the matrices and defect
!# vectors, routines of the module cc2dmediumm2matvecassembly are used.
!#
!# </purpose>
!##############################################################################

module cc2dmediumm2nonlinearcoreinit

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
  use linearsolverautoinitialise
  use matrixrestriction
  
  use collection
  use convection
    
  use cc2dmediumm2basic
  use cc2dmediumm2nonlinearcore
  
  implicit none
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_allocSystemMatrix (rproblem,rlevelInfo,rmatrix)
  
!<description>
  ! Allocates memory for the system matrix. rlevelInfo provides information
  ! about the level where the system matrix should be created.
  !
  ! Before this routine is called, the structure of all matrices in
  ! rlevelInfo must have been set up!
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(IN) :: rproblem

  ! A level-info structure specifying the matrices of the problem.
  type(t_problem_lvl), intent(IN), target :: rlevelInfo
!</input>

!<output>
  ! A block matrix that receives the basic system matrix.
  type(t_matrixBlock), intent(OUT) :: rmatrix
!</output>

!</subroutine>

    ! local variables
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A t_ccmatrixComponents used for defining the matrix structure
    type(t_ccmatrixComponents) :: rmatrixAssembly
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%p_rdiscretisation
    
    ! Get a pointer to the template FEM matrix.
    p_rmatrixTemplateFEM => rlevelInfo%rmatrixTemplateFEM

    ! In the global system, there are two gradient matrices B1 and B2.
    ! Get a pointer to the template structure for these.
    p_rmatrixTemplateGradient => rlevelInfo%rmatrixTemplateGradient

    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
      
    ! Let's consider the global system in detail. It has roughly
    ! the following shape:
    !
    !    ( A11       B1  ) = ( A11  A12  A13 )
    !    (      A22  B2  )   ( A21  A22  A23 )
    !    ( B1^T B2^T .   )   ( A31  A32  A33 )
    !
    ! All matices may have multiplication factors in their front.
    !
    ! The structure of the matrices A11 and A22 of the global system matrix
    ! is governed by the template FEM matrix.
    ! Initialise them with the same structure, i.e. A11, A22 share (!) their
    ! structure (not the entries) with that of the template matrix.
    !
    ! We allocate the system matrix using the c2d2_assembleMatrix routine.
    ! For this purpose, we have to initialise a t_ccmatrixComponents structure
    ! which defines the shape of the matrix. We simply set the parameters
    ! of those terms which wshould appear in the matrix to a value <> 0,
    ! that's enough for the memory allocation.
    
    rmatrixAssembly%dtheta = 1.0_DP   ! A velocity block
    rmatrixAssembly%deta = 1.0_DP     ! A gradient block
    rmatrixAssembly%dtau = 1.0_DP     ! A divergence block
    rmatrixAssembly%p_rdiscretisation => rlevelInfo%p_rdiscretisation
    rmatrixAssembly%p_rmatrixTemplateFEM => rlevelInfo%rmatrixTemplateFEM
    rmatrixAssembly%p_rmatrixTemplateGradient => rlevelInfo%rmatrixTemplateGradient
    rmatrixAssembly%p_rmatrixStokes => rlevelInfo%rmatrixStokes
    rmatrixAssembly%p_rmatrixB1 => rlevelInfo%rmatrixB1
    rmatrixAssembly%p_rmatrixB2 => rlevelInfo%rmatrixB2
    rmatrixAssembly%p_rmatrixMass => rlevelInfo%rmatrixMass

    if (.not. rproblem%bdecoupledXY) then
      call c2d2_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_AUTOMATIC,&
        rmatrix,rmatrixAssembly)
    else
      call c2d2_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_DECOUPLED,&
        rmatrix,rmatrixAssembly)
    end if
                                  
    ! That's it, all submatrices are set up.
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_initNonlinearLoop (rproblem,nlmin,nlmax,rvector,rrhs,&
      rnonlinearIteration,sname)
  
!<description>
  ! Initialises the given nonlinear iteration structure rnonlinearIteration.
  ! Creates the structure with c2d2_createNonlinearLoop and saves all
  ! problem dependent parameters and matrices in it.
  ! The routine automatically calls c2d2_createNonlinearLoop to initialise
  ! the structure with internal parameters.
  ! Note: This does not initialise the preconditioner in that structure!
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! Minimum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners.
  integer, intent(IN) :: nlmin
  
  ! Maximum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners. This level must correspond to rvector and rrhs.
  integer, intent(IN) :: nlmax

  ! The solution vector which is modified later during the nonlinear iteration.
  type(t_vectorBlock), intent(IN), target :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN), target :: rrhs

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(IN) :: sname
!</input>

!<output>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Is filled with data.
  type(t_ccnonlinearIteration), intent(OUT) :: rnonlinearIteration
!</output>

!</subroutine>

    ! local variables
    integer :: ilevel
    logical :: bneumann
    type(t_parlstSection), pointer :: p_rsection

    ! Basic initialisation of the nonlinenar iteration structure.
    call c2d2_createNonlinearLoop (rnonlinearIteration,rproblem%NLMIN,rproblem%NLMAX)
    
    rnonlinearIteration%MT_OutputLevel = rproblem%MT_OutputLevel
    
    ! Get the minimum/maximum damping parameter from the parameter list, save
    ! them to the nonlinear iteration structure (which is now initialised).
    call parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMin', rnonlinearIteration%domegaMin, 0.0_DP)
                              
    call parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMax', rnonlinearIteration%domegaMax, 2.0_DP)
    
    ! Save pointers to the RHS and solution vector
    rnonlinearIteration%p_rsolution => rvector
    rnonlinearIteration%p_rrhs => rrhs
    
    ! Set the preconditioner to 'nothing'
    rnonlinearIteration%rpreconditioner%ctypePreconditioning = -1
    
    ! Deactivate any 'tweak' flags in the final-assembly structure
    rnonlinearIteration%rfinalAssembly%iBmatricesTransposed = NO
    rnonlinearIteration%rfinalAssembly%iadaptiveMatrices = 0
    
    ! Assign the matrix pointers in the nonlinear iteration structure to
    ! all our matrices that we want to use.
    do ilevel = nlmin,nlmax
      rnonlinearIteration%RcoreEquation(ilevel)%p_rsystemMatrix => &
        rproblem%RlevelInfo(ilevel)%rpreallocatedSystemMatrix
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixStokes => &
        rproblem%RlevelInfo(ilevel)%rmatrixStokes
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB1 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB1
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB2 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB2

      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixMass => &
        rproblem%RlevelInfo(ilevel)%rmatrixMass

      rnonlinearIteration%RcoreEquation(ilevel)%p_rtempVector => &
        rproblem%RlevelInfo(ilevel)%rtempVector

    end do
      
    ! Clear auxiliary variables for the nonlinear iteration
    rnonlinearIteration%DresidualInit = 0.0_DP
    rnonlinearIteration%DresidualOld  = 0.0_DP
    
    call parlst_querysection(rproblem%rparamList, sname, p_rsection)

    if (.not. associated(p_rsection)) then
      print *,'Nonlinear solver not available; no section '''&
              //trim(sname)//'''!'
      stop
    end if

    ! Get stopping criteria of the nonlinear iteration
    call parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnonlinearIteration%DepsNL(1), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnonlinearIteration%DepsNL(2), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnonlinearIteration%DepsNL(3), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnonlinearIteration%DepsNL(4), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'ddmpD', &
                                 rnonlinearIteration%DepsNL(5), 0.1_DP)
      
    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions. This filter chain is applied to each
    ! defect vector during the linear and nonlinear iteration.
    allocate(rnonlinearIteration%p_RfilterChain(3))
    
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter for defect vectors:
    rnonlinearIteration%p_RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! The second filter filters for boundary conditions of fictitious boundary
    ! components
    rnonlinearIteration%p_RfilterChain(2)%ifilterType = FILTER_DISCBCDEFFICT
    
    ! Do we have Neumann boundary?
    !
    ! The bhasNeumannBoundary flag of the higher level decides about that...
    bneumann = rproblem%RlevelInfo(rproblem%NLMAX)%bhasNeumannBoundary
    rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_DONOTHING
    if (.not. bneumann) then
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_TOL20
      rnonlinearIteration%p_RfilterChain(3)%itoL20component = NDIM2D+1
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_doneNonlinearLoop (rnonlinearIteration)
  
!<description>
  ! Releases memory allocated in c2d2_initNonlinearLoop.
  ! The routine automatically calls c2d2_releaseNonlinearLoop to release
  ! internal parameters.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that should be cleaned up.
  type(t_ccNonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>
    
    ! Release the filter chain for the defect vectors.
    if (associated(rnonlinearIteration%p_RfilterChain)) &
      deallocate(rnonlinearIteration%p_RfilterChain)
      
    call c2d2_releaseNonlinearLoop (rnonlinearIteration)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_initPreconditioner (rproblem,rnonlinearIteration,rvector,rrhs)
  
!<description>
  ! This routine prepares the preconditioner that us used during the
  ! nonlinear iteration. The structure rpreconditioner will be initialised
  ! based on the information in rproblem.
  ! Necessary variables will be added to the collection structure in
  ! rproblem\%rcollection to be available in the callback routines.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  type(t_vectorBlock), intent(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! This is configured according to the preconditioner as specified in
  ! the DAT files.
  type(t_ccnonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
    integer(PREC_VECIDX) :: imaxmem
    character(LEN=PARLST_MLDATA) :: ssolverName,sstring

    ! At first, ask the parameters in the INI/DAT file which type of
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    call parlst_getvalue_int_direct (rproblem%rparamList, 'CC2D-NONLINEAR', &
        'itypePreconditioning', &
        rnonlinearIteration%rpreconditioner%ctypePreconditioning, 1)
    
    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioner
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX
      
      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      call parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      if (sstring .ne. '') read (sstring,*) ssolverName
      if (ssolverName .eq. '') then
        print *,'No linear subsolver!'
        stop
      end if
                                    
      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      allocate(rnonlinearIteration%rpreconditioner%p_rprojection)
      call mlprj_initProjectionVec (&
          rnonlinearIteration%rpreconditioner%p_rprojection,rrhs)
      
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      call c2d2_getProlRest (rnonlinearIteration%rpreconditioner%p_rprojection, &
          rproblem%rparamList,  'CC-PROLREST')
      
      ! Initialise the linear subsolver using the parameters from the INI/DAT
      ! files, the prepared filter chain and the interlevel projection structure.
      ! This gives us the linear solver node rpreconditioner%p_rsolverNode
      ! which identifies the linear solver.
      !
      ! Note that we pass rpreconditioner%p_RfilterChain as filter chain here,
      ! i.e. we use the same filter chain for all levels. This is a lack in the
      ! design of linsolinit_initFromFile as it forces us to use the same
      ! filters on all levels, which is not always advisable! (e.g. if
      ! a higher level 'sees' Neumann boundary while a lower one doesn't,
      ! there has actually another filter chain to be used on the lower level
      ! than on the higher one...)
      ! We probably change this later...
      call linsolinit_initFromFile (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
                                    rproblem%rparamList,ssolverName,&
                                    NLMAX-NLMIN+1,&
                                    rnonlinearIteration%p_RfilterChain,&
                                    rnonlinearIteration%rpreconditioner%p_rprojection)
      
      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...

      imaxmem = 0
      do i=NLMIN+1,NLMAX
        ! Pass the system metrices on the coarse/fine grid to
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = max(imaxmem,mlprj_getTempMemoryDirect (&
            rnonlinearIteration%rpreconditioner%p_rprojection,&
            rproblem%RlevelInfo(i-1)% &
              p_rdiscretisation%RspatialDiscr(1:rrhs%nblocks),&
            rproblem%RlevelInfo(i)% &
              p_rdiscretisation%RspatialDiscr(1:rrhs%nblocks)))
      end do
      
      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      allocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      call lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc,&
                                max(imaxmem,rrhs%NEQ),.false.)
      
      ! Set up a second temporary vector that we need for calculating
      ! the optimal defect correction.
      allocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      call lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2,&
                                rrhs%NEQ,.false.,ST_DOUBLE)
      
      ! Initialise the matrices.
      call c2d2_updatePreconditioner (rproblem,rnonlinearIteration,&
          rvector,rrhs,.true.,.true.)
      
    case DEFAULT
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_updatePreconditioner (rproblem,rnonlinearIteration,&
      rvector,rrhs,binit,bstructuralUpdate)
  
!<description>
  ! This routine has to be called whenever the system matrices change.
  ! It initialises (depending on the system matrices) the matrices of the
  ! preconditioner or performs an update of them.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  type(t_vectorBlock), intent(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs

  ! First initialisation.
  ! Has to be set to TRUE on the first call. Initialises the preconditioner
  ! with the structure of the matrices.
  logical, intent(IN) :: binit

  ! Whether the structure of the system matrices is new.
  ! This variable has to be set to TRUE whenever there was a structure in
  ! the system matrices. This reinitialises the linear solver.
  logical, intent(IN) :: bstructuralUpdate
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Preconditioner data is saved here.
  type(t_ccnonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
    character(LEN=PARLST_MLDATA) :: sstring,snewton

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A pointer to the matrix of the preconditioner
    type(t_matrixBlock), pointer :: p_rmatrixPreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(NNLEV) :: Rmatrices
    
    ! Pointer to the template FEM matrix
    type(t_matrixScalar), pointer :: p_rmatrixTempateFEM
    
    
    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioner
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
    
      if (.not. binit) then
        ! Restore the standard matrix structure in case the matrices had been
        ! modified by c2d2_finaliseMatrices for the preconditioner -- i.e.
        ! temporarily switch to the matrix structure that is compatible to
        ! the discretisation.
        call c2d2_unfinaliseMatrices (rnonlinearIteration, &
            rnonlinearIteration%rfinalAssembly,.false.)
      end if
    
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX
      
      ! Initialise the preconditioner matrices on all levels.
      do i=NLMIN,NLMAX
      
        ! Prepare the preconditioner matrices level i. This is
        ! basically the system matrix...
        allocate(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        p_rmatrixPreconditioner => &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner
      
        ! The preallocated system matrix is the basic building block for
        ! the preconditioning matrix. But probably the preconditioning
        ! matrix contains some more entries... see below!
        call lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix,&
            p_rmatrixPreconditioner,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
            
        ! ----------------------------------------------------
        ! Should the linear solver use the Newton matrix?
        if ((rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CCPREC_NEWTON) .or. &
            (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CCPREC_NEWTONDYNAMIC)) then
          ! That means, our preconditioner matrix must look like
          !
          !  A11  A12  B1
          !  A21  A22  B2
          !  B1^T B2^T 0
          !
          ! With A12, A21, A11, A22 independent of each other!
          ! Do we have that case? If not, we have to allocate memory
          ! for these matrices.
          p_rmatrixTempateFEM => rproblem%RlevelInfo(i)%rmatrixTemplateFEM
          
          ! If we have a Stokes problem, A12 and A21 don't exist.
          if (rproblem%iequation .eq. 0) then

            if (p_rmatrixPreconditioner%RmatrixBlock(1,2)%cmatrixFormat &
                .eq. LSYSSC_MATRIXUNDEFINED) then
                
              call lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                p_rmatrixPreconditioner%RmatrixBlock(1,2), &
                LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
              ! Allocate memory for the entries; don't initialise the memory.
              call lsyssc_allocEmptyMatrix (&
                  p_rmatrixPreconditioner%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
                
            end if

            if (p_rmatrixPreconditioner%RmatrixBlock(2,1)%cmatrixFormat &
                .eq. LSYSSC_MATRIXUNDEFINED) then
                
              ! Create a new matrix A21 in memory. create a new matrix
              ! using the template FEM matrix...
              call lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                p_rmatrixPreconditioner%RmatrixBlock(2,1), &
                LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
              ! Allocate memory for the entries; don't initialise the memory.
              call lsyssc_allocEmptyMatrix (&
                  p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
               
            end if
            
          end if

          ! A22 may share its entries with A11. If that's the case,
          ! allocate additional memory for A22, as the Newton matrix
          ! requires a separate A22!
          if (lsyssc_isMatrixContentShared( &
              p_rmatrixPreconditioner%RmatrixBlock(1,1), &
              p_rmatrixPreconditioner%RmatrixBlock(2,2)) ) then
            ! Release the matrix structure. As the matrix is a copy
            ! of another one, this will clean up the structure but
            ! not release any memory.
            call lsyssc_releaseMatrix ( &
                p_rmatrixPreconditioner%RmatrixBlock(2,2))

            ! Create a new matrix A21 in memory. create a new matrix
            ! using the template FEM matrix...
            call lsyssc_duplicateMatrix ( &
              p_rmatrixPreconditioner%RmatrixBlock(1,1), &
              p_rmatrixPreconditioner%RmatrixBlock(2,2),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
              
            ! ... then allocate memory for the entries;
            ! don't initialise the memory.
            call lsyssc_allocEmptyMatrix (&
                p_rmatrixPreconditioner%RmatrixBlock(2,2),&
                LSYSSC_SETM_UNDEFINED)
          end if
          
          ! ----------------------------------------------------
          ! Should we even use the adaptive Newton?
          if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CCPREC_NEWTONDYNAMIC) then
            
            ! We have even the extended, dynamic Newton as preconditioner.
            ! Put the parameters for the extended Newton from the DAT file
            ! into the Adaptive-Newton configuration block.
            
            call parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                        'spreconditionerAdaptiveNewton', sstring, '')
            snewton = ''
            if (sstring .ne. '') read (sstring,*) snewton
            if (snewton .ne. '') then
              ! Initialise the parameters of the adaptive Newton
              call parlst_getvalue_int (rproblem%rparamList, snewton, &
                  'nminFixPointIterations', rnonlinearIteration%rpreconditioner% &
                  radaptiveNewton%nminFixPointIterations, 0)

              call parlst_getvalue_int (rproblem%rparamList, snewton, &
                  'nmaxFixPointIterations', rnonlinearIteration%rpreconditioner% &
                  radaptiveNewton%nmaxFixPointIterations, 999)

              call parlst_getvalue_double (rproblem%rparamList, snewton, &
                  'depsAbsNewton', rnonlinearIteration%rpreconditioner% &
                  radaptiveNewton%depsAbsNewton, 1E-5_DP)

              call parlst_getvalue_double (rproblem%rparamList, snewton, &
                  'depsRelNewton', rnonlinearIteration%rpreconditioner% &
                  radaptiveNewton%depsRelNewton, 1E99_DP)
            end if
            
          end if
          
        end if
      end do
      
      if (binit) then
        ! Check the matrices if they are compatible to our
        ! preconditioner. If not, we later have to modify the matrices a little
        ! bit to make it compatible.
        ! The result of this matrix analysis is saved to the rfinalAssembly structure
        ! in rnonlinearIteration and allows us later to switch between these two
        ! matrix representations: Compatibility to the discretisation routines
        ! and compatibity to the preconditioner.
        ! The c2d2_checkAssembly routine below uses this information to perform
        ! the actual modification in the matrices.
        call c2d2_checkAssembly (rproblem,rnonlinearIteration,rrhs,&
            rnonlinearIteration%rfinalAssembly)
      end if
      ! Otherwise, checkAssembly was already called and does not have to be
      ! called again.
      
      ! Using rfinalAssembly as computed above, make the matrices compatible
      ! to our preconditioner if they are not.
      call c2d2_finaliseMatrices (rnonlinearIteration)
      
      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
      do i=NLMIN,NLMAX
        call lsysbl_duplicateMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
          Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end do
      
      call linsol_setMatrices(&
          rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(NLMIN:NLMAX))
          
      ! The solver got the matrices; clean up Rmatrices, it was only of temporary
      ! nature...
      do i=NLMIN,NLMAX
        call lsysbl_releaseMatrix (Rmatrices(i))
      end do
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      if (binit) then
        call linsol_initStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
            ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
      else if (bstructuralUpdate) then
        call linsol_updateStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
            ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
      end if
      
    case DEFAULT
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine c2d2_releasePreconditioner (rnonlinearIteration)
  
!<description>
  ! This routine releases the preconditioner for the nonlinear iteration
  ! which was prepared in c2d2_initPreconditioner. Memory is released
  ! from heap.
!</description>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! The preconditioner data is removed from that,
  type(t_ccnonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Which preconditioner do we have?
    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioning
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the preconditioner matrix on every level
      do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        call lsysbl_releaseMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        deallocate(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
      end do
      
      ! Release the temporary vector(s)
      call lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      deallocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      
      call lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      deallocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      
      ! Clean up data about the projection etc.
      call mlprj_doneProjection(rnonlinearIteration%rpreconditioner%p_rprojection)
      deallocate(rnonlinearIteration%rpreconditioner%p_rprojection)

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      call linsol_doneStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode)
      call linsol_releaseSolver (rnonlinearIteration%rpreconditioner%p_rsolverNode)
      
    case DEFAULT
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_getProlRest (rprojection, rparamList, sname)
  
!<description>
  ! Initialises an existing interlevel projection structure rprojection
  ! with parameters from the INI/DAT files. sname is the section in the
  ! parameter list containing parameters about prolongation restriction.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(IN) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the prolongation/restriction.
  character(LEN=*), intent(IN) :: sname
!</input>

!<output>
  ! An interlevel projection block structure containing an initial
  ! configuration of prolongation/restriction. The structure is modified
  ! according to the parameters in the INI/DAT file(s).
  type(t_interlevelProjectionBlock), intent(INOUT) :: rprojection
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: i1
    real(DP) :: d1

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, sname, p_rsection)

    if (.not. associated(p_rsection)) then
      ! We use the default configuration; stop here.
      return
    end if
    
    ! Now take a look which parameters appear in that section.

    ! Prolongation/restriction order for velocity components
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
    end if

    ! Prolongation/restriction order for pressure
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,NDIM2D+1)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%iinterpolationOrder = i1
    end if
    
    ! Prolongation/restriction variant for velocity components
    ! in case of Q1~ discretisation
    call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
    
    if (i1 .ne. -1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1
    end if
    
    ! Aspect-ratio indicator in case of Q1~ discretisation
    ! with extended prolongation/restriction
    call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
    
    if (i1 .ne. 1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
    end if

    ! Aspect-ratio bound for switching to constant prolongation/restriction
    ! in case of Q1~ discretisation with extended prolongation/restriction
    call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
    
    if (d1 .ne. 20.0_DP) then
      rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
    end if

  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_checkAssembly (rproblem,rnonlinearIteration,rrhs,rfinalAssembly)
  
!<description>
  ! This routine checks the matrices against an existing preconditioner.
  ! It may happen that e.g. VANCA does not like our matrices (probably
  ! they have to be saved transposed or whatever). Information about
  ! which things must be changed in the assembly to make 'everything proper'
  ! is saved to rfinalAssembly.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! Nonlinar iteration structure saving data about the actual configuration
  ! of the core equation.
  type(t_ccnonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!<input>
  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs
!</input>

!<inputoutput>
  ! Nonlinear iteration structure.
  ! The t_ccFinalAssemblyInfo substructure that receives information how to
  ! finally assembly the matrices such that everything in the callback routines
  ! will work.
  type(t_ccFinalAssemblyInfo), intent(INOUT) :: rfinalAssembly
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: NLMIN,NLMAX,ccompatible,iprecType
  character(LEN=PARLST_MLDATA) :: ssolverName,sstring
  type(t_interlevelProjectionBlock) :: rprojection

  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  type(t_matrixBlock), dimension(NNLEV) :: Rmatrices
  type(t_linsolNode), pointer :: p_rsolverNode
    
    ! At first, ask the parameters in the INI/DAT file which type of
    ! preconditioner is to be used.
    call parlst_getvalue_int_direct (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                     'itypePreconditioning', &
                                     iprecType, 1)

    select case (iprecType)
    case (CCPREC_NONE)
      ! No preconditioning
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! That preconditioner is a solver for a linear system.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX

      ! Temporarily set up the solver node for the linear subsolver.
      !
      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      call parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      if (sstring .ne. '') read (sstring,*) ssolverName
      if (ssolverName .eq. '') then
        print *,'No linear subsolver!'
        stop
      end if
                                    
      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      call mlprj_initProjectionVec (rprojection,rrhs)
      
      ! Initialise the linear subsolver using the parameters from the INI/DAT
      ! files, the prepared filter chain and the interlevel projection structure.
      ! This gives us the linear solver node p_rsolverNode
      ! which identifies the linear solver.
      call linsolinit_initFromFile (p_rsolverNode,&
                                    rproblem%rparamList,ssolverName,&
                                    NLMAX-NLMIN+1,&
                                    rinterlevelProjection=rprojection)
      
      ! Check the matrices.
      !
      ! We copy our matrices to a big matrix array and transfer that
      ! to the setMatrices routines. This intitialises then the matrices
      ! on all levels according to that array.
      Rmatrices(NLMIN:NLMAX) = rproblem%RlevelInfo(NLMIN:NLMAX)%rpreallocatedSystemMatrix
      call linsol_matricesCompatible(p_rsolverNode, &
          Rmatrices(NLMIN:NLMAX),ccompatible)
      
      select case (ccompatible)
      case (LINSOL_COMP_OK) ! nothing to do
        rfinalAssembly%iBmatricesTransposed = NO
      case (LINSOL_COMP_ERRTRANSPOSED)
        ! The B-matrices must be assembled in a transposed way. Remember that.
        rfinalAssembly%iBmatricesTransposed = YES
      case DEFAULT
        print *,'Preconditioner incompatible to the matrices. Don''t know why!?!'
        stop
      end select
      
      ! Release the solver node again.
      call linsol_releaseSolver (p_rsolverNode)
      
      ! We also don't need the temporary projection structure anymore.
      call mlprj_doneProjection (rprojection)
      
    end select

    ! Add information about adaptive matrix generation from INI/DAT files
    ! to the collection.
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'iAdaptiveMatrix', rfinalAssembly%iadaptiveMatrices, 0)
                              
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'dAdMatThreshold', rfinalAssembly%dAdMatThreshold, 20.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_finaliseMatrices (rnonlinearIteration)
  
!<description>
  ! This routine performs final assembly tasks to the matrices such that they
  ! are compatible to the preconditioner.
  ! It may happen that e.g. VANCA does not like our matrices (probably
  ! they have to be saved transposed or whatever). In that case, we
  ! have to make slight modifications to our matrices in order to
  ! make them compatible.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that contains pointers to the
  ! matrices which are to be used.
  type(t_ccnonlinearIteration), intent(INOUT), target :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev,NLMIN,NLMAX
    type(t_matrixBlock), pointer :: p_rmatrix

    if (rnonlinearIteration%rfinalAssembly%iBmatricesTransposed .eq. YES) then
      ! There is usually a VANCA subsolver in the main linear solver which
      ! cannot deal with our virtually transposed matrices. So we make
      ! a copy of B1/B2 on every level and really transpose them.

      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX

      ! Loop through the levels, transpose the B-matrices
      do ilev=NLMIN,NLMAX
        ! Get the matrix of the preconditioner
        p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
        ! Release the old B1/B2 matrices from the system matrix. This
        ! does not release any memory, as the content of the matrix
        ! is saved elsewhere.
        call lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,1))
        call lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,2))
        
        ! Transpose B1/B2, write result to the system matrix.
        call lsyssc_transposeMatrix (&
            rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1, &
            p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_ALL)

        call lsyssc_transposeMatrix (&
            rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2, &
            p_rmatrix%RmatrixBlock(3,2),LSYSSC_TR_ALL)
                                    
        ! Release the memory that was allocated for the B2 structure by
        ! the matrix-transpose routine.
        ! Replace the structure of B2 by that of B1; more precisely,
        ! share the structure. We can do this as we know that B1 and B2
        ! have exactly the same structure!
        call lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_REMOVE,LSYSSC_DUP_IGNORE)
        call lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
      
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_unfinaliseMatrices (rnonlinearIteration,rfinalAssembly,bdata)
  
!<description>
  ! Reverts the changes that were done by c2d2_finaliseMatrices, brings
  ! matrices back to their 'standard' form. Necessary if multiple simulations
  ! are performed, e.g. in a time-dependent simulation: The assembly
  ! routines expect matrices in the standard form and cannot work with
  ! 'specially tuned' matrices which are used in the actual solving/
  ! preconditioning process.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that contains pointers to the
  ! matrices which are to be used.
  type(t_ccnonlinearIteration), intent(INOUT), target :: rnonlinearIteration
!</inputoutput>

!<input>
  ! The t_ccFinalAssemblyInfo structure that receives information how to finally
  ! assembly the matrices such that everything in the callback routines will
  ! work. Must be set up with c2d2_checkAssembly above.
  type(t_ccFinalAssemblyInfo), intent(IN) :: rfinalAssembly
  
  ! TRUE  = restore the original matrices in their whole -- structure and data
  ! FALSE = Ignore the data, only restore the original matrix structure.
  !         Used to save time if the matrix content is thrown away in the
  !         reassembling process anyway.
  logical, intent(IN) :: bdata
!</input>

!</subroutine>

    ! local variables
    integer :: ilev,NLMIN,NLMAX
    type(t_matrixBlock), pointer :: p_rmatrix

    ! Which levels have we to take care of during the solution process?
    NLMIN = rnonlinearIteration%NLMIN
    NLMAX = rnonlinearIteration%NLMAX

    ! Loop through the levels, transpose the B-matrices
    do ilev=NLMIN,NLMAX
      ! Get the matrix of the preconditioner
      p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
      if ((rfinalAssembly%iBmatricesTransposed .eq. YES) .and. &
          (iand(p_rmatrix%RmatrixBlock(3,1)%imatrixSpec,&
                LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
                
        ! There is usually a VANCA subsolver in the main linear solver which
        ! cannot deal with our virtually transposed matrices.
        ! The B1/B2 matrices are transposed -- so we re-transpose them
        ! to get the original matrices.

        ! Release the old B1/B2 matrices from the system matrix. This
        ! does not release any memory, as the content of the matrix
        ! is saved elsewhere.
        call lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,1))
        call lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,2))
        
        ! Transpose B1/B2, write result to the system matrix.
        if (bdata) then
          ! Structure and data
          call lsyssc_transposeMatrix (&
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1, &
              p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_ALL)

          call lsyssc_transposeMatrix (&
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2, &
              p_rmatrix%RmatrixBlock(3,2),LSYSSC_TR_ALL)
        else
          ! Only the structure; content gets invalid.
          call lsyssc_transposeMatrix (&
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1, &
              p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_STRUCTURE)

          ! No change has to be done to the B2^T block in the global system
          ! matrix since that one will later share the same structure as B1^T!
          ! So the following command is commented out and should
          ! not be commented in!
          ! CALL lsyssc_transposeMatrix (&
          !     rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2, &
          !     p_rmatrix%RmatrixBlock(3,2),LSYSSC_TR_STRUCTURE)
        end if
                                    
        ! Release the memory that was allocated for the B2 structure by
        ! the matrix-transpose routine.
        ! Replace the structure of B2 by that of B1; more precisely,
        ! share the structure. We can do this as we know that B1 and B2
        ! have exactly the same structure!
        call lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_REMOVE,LSYSSC_DUP_IGNORE)
        call lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
      
      end if

    end do
      
  end subroutine

end module
