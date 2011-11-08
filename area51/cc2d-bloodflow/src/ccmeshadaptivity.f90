!##############################################################################
!# ****************************************************************************
!# <name> ccmeshadaptivity </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains postprocessingmesh adaptation routines for the CC2D solver.
!#
!# The following routines can be found here:
!#
!# 1.) cc_initAdaptivity
!#     -> Initialise the mesh adaptation.
!#
!# 2.) cc_doneAdaptivity
!#     -> Clean up the mesh adaptation.
!#
!# 3.) cc_generateMeshIndicator
!#     -> Generate the mesh indicator vector
!#
!# 4.) cc_performMeshAdaptation
!#     -> Perform mesh adaptation
!#
!# Auxiliary routines:
!#
!#
!# </purpose>
!##############################################################################

module ccmeshadaptivity

  use fsystem
  use storage
  use hadaptaux
  use hadaptivity
  use linearsystemscalar
  use triangulation
  use dofmapping

  use ccbasic
  use ccgeneraldiscretisation
  use ccpostprocessing
  use ccboundarycondition
  
  implicit none

contains
   
  ! ***************************************************************************

!<subroutine>

  subroutine cc_initAdaptation (rproblem)
  
!<description>
  ! Initialises the mesh adaptation from the triangulation.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Initialise switch for mesh adaptation
    rproblem%bmeshAdaptation = .false.

    ! Initialise the mesh adaptation on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      
      ! First get parameters from parameter list
      call hadapt_initFromParameterlist(rproblem%RlevelInfo(i)%rhadapt,&
          rproblem%rparamList,'MESHADAPTIVITY')
      
      ! Set up adaptation structure for the given triangulation
      if (rproblem%RlevelInfo(i)%rhadapt%iadaptationStrategy .ne. 0) then
        call hadapt_initFromTriangulation (rproblem%RlevelInfo(i)%rhadapt,&
            rproblem%RlevelInfo(i)%rtriangulation)

        ! Set global switch for mesh adaptation
        rproblem%bmeshAdaptation = .true.
      end if
    end do
    
  end subroutine cc_initAdaptation
 
  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneAdaptation (rproblem)
  
!<description>
  ! Releases the mesh adaptation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    ! Release the mesh adaptation on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      call hadapt_releaseAdaptation (rproblem%RlevelInfo(i)%rhadapt)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateMeshIndicator (rproblem,rindicator)
  
!<description>
  ! Generates the mesh indicator vector for the monitor function
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(in), target :: rproblem
!/<input>

!<inputoutput>
  ! Scalar mesh indicator vector
  type(t_vectorScalar), intent(inout) :: rindicator
!</inputoutput>

!</subroutine>

    ! Pointer to the top-level triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to the top-level mesh adaptation structure
    type(t_hadapt), pointer :: p_rhadapt

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_Ddata

    real(DP) :: dx, dy, dx0, dy0, tol
    real(DP) :: dxRef1, dyRef1, dxRef2, dyRef2, dxRef3, dyRef3
    real(DP) :: dtime, phi1, phi2, phi3, width, height
    integer :: iel, ivt, ive

    ! Set pointers
    p_rtriangulation => rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation
    p_rhadapt => rproblem%RlevelInfo(rproblem%NLMAX)%rhadapt

    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoordinates)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)

    ! Generate fresh indicator vector
    if(rindicator%NEQ .ne. 0) call lsyssc_releaseVector(rindicator)
    call lsyssc_createVector(rindicator, p_rtriangulation%NEL, .true.)
    call lsyssc_getbase_double(rindicator, p_Ddata)

    
    ! Definition of the rotor
    dx0 = 0.0_DP; dy0 = 0.0_DP
    width = 0.1; height = 0.6; tol=0.02
    phi1  = tanh(dtime)*2*SYS_PI*dtime
    phi2  = tanh(dtime)*2*SYS_PI*dtime+2*SYS_PI/3.0
    phi3  = tanh(dtime)*2*SYS_PI*dtime+4*SYS_PI/3.0
    
    ! Loop over all elements
    elem: do iel = 1, p_rtriangulation%NEL

      ! Loop over all vertices of the element
      vert: do ive = 1, tria_getNVE(p_IverticesAtElement, iel)

        ! Get vertex number and its coordinates
        ivt = p_IverticesAtElement(ive, iel)
        dx  = p_DvertexCoordinates(1, ivt)
        dy  = p_DvertexCoordinates(2, ivt)
        
        ! Compute reference coordinates for the rotor
        dxRef1 = cos(phi1)*(dx-dx0)-sin(phi1)*(dy-dy0)
        dyRef1 = sin(phi1)*(dx-dx0)+cos(phi1)*(dy-dy0)
        dxRef2 = cos(phi2)*(dx-dx0)-sin(phi2)*(dy-dy0)
        dyRef2 = sin(phi2)*(dx-dx0)+cos(phi2)*(dy-dy0)
        dxRef3 = cos(phi3)*(dx-dx0)-sin(phi3)*(dy-dy0)
        dyRef3 = sin(phi3)*(dx-dx0)+cos(phi3)*(dy-dy0)
        
        ! Definition of the rotor boundary
        if ((dxRef1 .ge.-width-tol)  .and. (dxRef1 .le. width+tol) .and.&
            (dyRef1 .le. height+tol) .and. (dyRef1 .ge. -tol) .and.&
            ((dxRef1 .le. -width+tol) .or. (dxRef1 .ge. width-tol)&
             .or. (dyRef1 .ge. height-tol)) .or.&
            (dxRef2 .ge.-width-tol)  .and. (dxRef2 .le. width+tol) .and.&
            (dyRef2 .le. height+tol) .and. (dyRef2 .ge. -tol) .or.&
            (dxRef3 .ge.-width+tol)  .and. (dxRef3 .le. width+tol) .and.&
            (dyRef3 .le. height+tol) .and. (dyRef3 .ge. -tol)) then

          ! Mark element for refinement
          p_Ddata(iel) = 1.0
          
          cycle elem
        end if

      end do vert
    end do elem
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_performMeshAdaptation (rproblem,rindicator,rvector,&
      rrhs,rpostprocessing)
  
!<description>
  ! Performs mesh adaptation
!</description>

!<input>
  ! Scalar mesh indicator vector
  type(t_vectorScalar), intent(in) :: rindicator
!/<input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! The initial solution vector.
  type(t_vectorBlock), intent(inout), target :: rvector
  
  ! The initial RHS vector.
  type(t_vectorBlock), intent(inout) :: rrhs

  ! Postprocessing structure. Defines what to do with solution vectors.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Pointer to the top-level triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to the top-level mesh adaptation structure
    type(t_hadapt), pointer :: p_rhadapt

    ! local variables
    type(t_collection) :: rcollection
    integer :: i

    ! Set pointers
    p_rtriangulation => rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation
    p_rhadapt => rproblem%RlevelInfo(rproblem%NLMAX)%rhadapt

    ! Initialise the collection
    call collct_init (rcollection)
    rcollection%p_rvectorQuickAccess1 => rvector

    ! Initialise the callback function
    call performAdaptation(HADAPT_OPR_INITCALLBACK,rcollection)

    ! Perform one step h-adaptivity
    call hadapt_refreshAdaptation(p_rhadapt,p_rtriangulation)
    call hadapt_performAdaptation(p_rhadapt,rindicator,&
        rcollection,performAdaptation)

    ! Release the collection
    call collct_done (rcollection)

    ! Generate raw mesh from adaptivity structure
    call hadapt_generateRawMesh(p_rhadapt,p_rtriangulation)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (p_rtriangulation,rproblem%rboundary)

    ! Print mesh information
    if (rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ('Mesh statistics:')
      call output_lbrk ()
      do i=rproblem%NLMIN,rproblem%NLMAX
        call tria_infoStatistics (rproblem%RlevelInfo(i)%rtriangulation,&
            i .eq. rproblem%NLMIN,i)
      end do
    end if

    ! Re-initialise discretisation on top-level grid
    if (rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising discretisation...')
    end if
    call cc_doneDiscretisation(rproblem,rproblem%NLMAX,rproblem%NLMAX)
    call cc_initDiscretisation(rproblem,rproblem%NLMAX,rproblem%NLMAX)

    if (rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ('Discretisation statistics:')
      do i=rproblem%NLMIN,rproblem%NLMAX
        call output_lbrk ()
        call output_line ('Level '//sys_siL(i,5))
        call dof_infoDiscrBlock (rproblem%RlevelInfo(i)%rdiscretisation,.false.)
      end do
    end if

    ! And all the other stuff...
    if (rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising postprocessing...')
    end if
    i=rpostprocessing%inextFileSuffixUCD
    call cc_donePostprocessing (rpostprocessing)
    call cc_initPostprocessing (rproblem,rpostprocessing)
    rpostprocessing%inextFileSuffixUCD=i

    if (rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising matrices/vectors...')
    end if
    call cc_doneMatVec (rproblem,rvector=rvector,rrhs=rrhs,&
        nlminOpt=rproblem%NLMAX,nlmaxOpt=rproblem%NLMAX)
    call cc_allocMatVec (rproblem,rvector=rvector,rrhs=rrhs,&
        nlminOpt=rproblem%NLMAX,nlmaxOpt=rproblem%NLMAX)

    call cc_initInitialSolution (rproblem,rvector)

    ! Generate the static matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    if (rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Generating basic matrices...')
    end if
    call cc_generateBasicMat (rproblem,rproblem%NLMAX,rproblem%NLMAX)
    
    ! Generate the RHS vector for the first time step.
    if (rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Generating RHS vector...')
    end if
    call cc_generateBasicRHS (rproblem,rrhs)

    ! Initialise the boundary conditions, but
    ! do not implement any boundary conditions as the nonstationary solver
    ! does not like this.
    if (rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Generating discrete boundary conditions...')
    end if
    call cc_initDiscreteBC (rproblem,rvector,rrhs)

  end subroutine
end module
