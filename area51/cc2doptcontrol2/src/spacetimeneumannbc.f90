!##############################################################################
!# ****************************************************************************
!# <name> spacetimeneumannbc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Encapsules the Neumann boundary conditions in space/time.
!# This type of boundary conditions is saved analytically for each time
!# discretisation as set of boundary segments. During the assembly,
!# of matrices in space, the analytical definition is evaluated and the
!# actual matrix entries are created based on this definition. The
!# definition is therefore independent of the space discretisation!
!# </purpose>
!##############################################################################

module spacetimeneumannbc

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use paramlist
  use timestepping
  use dofmapping
  use discretebc
  use discretefbc
  use derivatives
  
  use collection
  use convection
    
  use constantsoptc
  use structuresoptc
  use user_callback

  !use spacepreconditioner
  !use spacepreconditionerinit
  use spatialbcdef
  use spacediscretisation
  use timediscretisation
  use spacetimevectors
  
  use scalarpde
  use assemblytemplates

  !use timeanalysis
  
  !use spacetimediscretisation

  implicit none
  
  private
  
!<types>
  !<typeblock>
  
  ! Contains the definition of the Neumann boundary conditions for all timesteps
  ! in a space-time discretisation.
  type t_sptiNeumannBoundary
  
    ! Type of Neumann boundary conditions.
    ! =0: Position of the Neumann BC's fixed for all timesteps.
    ! =1: Neumann BC change in time.
    ! WARNING: Currently, only cneumannType=0 supported!!!
    integer :: cneumannType = 0
    
    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr => null()

    ! Underlying time-discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()
    
    ! Underlying assembly templates
    type(t_staticSpaceAsmTemplates), pointer :: p_rasmTemplates => null()
    
    ! Number of unknowns in time
    integer :: NEQtime = 0
    
    ! Analytical definition of the Neumann boundary segments
    type(t_neumannBoundary), dimension(:), pointer :: p_RneumannBoundary => null()
    
    ! Neumann boundary integral template matrix. This matrix contains the structure
    ! of an operator that works on the Neumann boundary. To save time, it is
    ! saved as 'differential' matrix to the standard velocity 'template' matrix.
    ! The matrix can be used for all timesteps as it contains a common structure.
    type(t_matrixScalar) :: rneumannBoudaryOperator
    
  end type
  
  !</typeblock>
!</types>

  public :: t_sptiNeumannBoundary
  public :: stnm_createNeumannBoundary
  public :: stnm_releaseNeumannBoundary
  public :: stnm_assembleNeumannBoundary
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine stnm_createNeumannBoundary (rspaceDiscr,rtimeDiscr,rasmTemplates,&
      rsptiNeumannBC)

!<description>
  ! Creates a space-time Neumann boundary definition structure
  ! from a time discretisation.
!</description>

!<input>
  ! Underlying space discretisation.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Underlying time-discretisation.
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
  
  ! Assembly template structure.
  type(t_staticSpaceAsmTemplates), intent(in), target :: rasmTemplates
!</input>

!<output>
  ! Structure to create.
  type(t_sptiNeumannBoundary), intent(out) :: rsptiNeumannBC
!</output>

!</subroutine>

    ! Initialise the structure.
    rsptiNeumannBC%p_rtimeDiscr => rtimeDiscr
    rsptiNeumannBC%p_rspaceDiscr => rspaceDiscr
    rsptiNeumannBC%p_rasmTemplates => rasmTemplates
    rsptiNeumannBC%NEQtime = rtimeDiscr%nintervals+1
    allocate(rsptiNeumannBC%p_RneumannBoundary(rsptiNeumannBC%NEQtime))
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stnm_releaseNeumannBoundary (rsptiNeumannBC)

!<description>
  ! Releases a space-time Neumann boundary definition structure.
!</description>

!<inputoutput>
  ! Structure to release.
  type(t_sptiNeumannBoundary), intent(inout) :: rsptiNeumannBC
!</inputoutput>

!</subroutine>

    ! Initialise the structure.
    deallocate(rsptiNeumannBC%p_RneumannBoundary)
    rsptiNeumannBC%NEQtime = 0
    nullify(rsptiNeumannBC%p_rtimeDiscr)
    nullify(rsptiNeumannBC%p_rspaceDiscr)
    nullify(rsptiNeumannBC%p_rasmTemplates)
    if (rsptiNeumannBC%rneumannBoudaryOperator%NA .ne. 0) then
      call lsyssc_releaseMatrix (rsptiNeumannBC%rneumannBoudaryOperator)
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stnm_assembleNeumannBoundary (roptcBDC,rsptiNeumannBC,rglobalData)

!<description>
  ! Initialises a space-time Neumann boundary definition structure with the
  ! definition of the Neumann boundary based on the analytical
  ! definitions in the callback routines.
!</description>

!<input>
  ! Boundary conditions in the problem.
  type(t_optcBDC), intent(in) :: roptcBDC

  ! Global data, passed to callback routines
  type(t_globalData), intent(inout) :: rglobalData
!</input>

!<inputoutput>
  ! Structure to release.
  type(t_sptiNeumannBoundary), intent(inout) :: rsptiNeumannBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j
    real(DP) :: dtimePrimal, dtimeDual, dtstep
    type(t_bilinearForm) :: rform
    type(t_neumannBdRegion), pointer :: p_rdualNeumannBd
    type(t_matrixScalar) :: rneumannBoudaryOperator
    
    ! Prepare a bilinear form.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff(1) = .true.
    rform%Dcoefficients(1:rform%itermCount) = 1.0_DP
    
    ! Create a temp matrix for boundary integrals. At first: Full structure, zero content.
    call lsyssc_duplicateMatrix (rsptiNeumannBC%p_rasmTemplates%rmatrixTemplateFEM,&
        rneumannBoudaryOperator, LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call lsyssc_clearMatrix (rneumannBoudaryOperator)
    
    ! Loop through all timesteps
    do i=1,rsptiNeumannBC%NEQtime
      ! Current point in time
      call tdiscr_getTimestep(rsptiNeumannBC%p_rtimeDiscr,i-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal - (1.0_DP-rsptiNeumannBC%p_rtimeDiscr%dtheta)*dtstep
    
      ! Calculate the Neumann BC's.
      call sbc_releaseNeumannBoundary(rsptiNeumannBC%p_RneumannBoundary(i))
      call sbc_assembleBDconditions (roptcBDC,dtimePrimal,dtimeDual,&
            CCSPACE_PRIMALDUAL,rglobalData,SBC_NEUMANN,&
            rsptiNeumannBC%p_rtimediscr,rsptiNeumannBC%p_rspacediscr,&
            rneumannBoundary=rsptiNeumannBC%p_RneumannBoundary(i))
            
      ! Assemble the operator on all boundary components in rneumannBoundary.
      p_rdualNeumannBd => rsptiNeumannBC%p_RneumannBoundary(i)%p_rdualNeumannBdHead
      do j = 1,rsptiNeumannBC%p_RneumannBoundary(i)%nregionsDual
        ! Discretise a mass-operator on the boundary to mark the rows which may contain
        ! boundary integrals.
        call bilf_buildMatrixScalarBdr2D (rform, CUB_G2_1D, .false., &
            rneumannBoudaryOperator,&
            rboundaryRegion=p_rdualNeumannBd%rboundaryRegion)

        ! Next segment            
        p_rdualNeumannBd => p_rdualNeumannBd%p_nextNeumannRegion
      end do
    end do
    
    ! Remove the matrix content and compress the matrix.
    ! The result is a matrix with the minimum stencil that can be used for
    ! all timesteps.
    if (rsptiNeumannBC%rneumannBoudaryOperator%NEQ .ne. 0) then
      call lsyssc_releaseMatrix (rsptiNeumannBC%rneumannBoudaryOperator)
    end if
    call lsyssc_createRowCMatrix(rneumannBoudaryOperator,rsptiNeumannBC%rneumannBoudaryOperator)
    !call lsyssc_duplicateMatrix (rneumannBoudaryOperator,rsptiNeumannBC%rneumannBoudaryOperator,&
    !    LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    call lsyssc_releaseMatrixContent(rsptiNeumannBC%rneumannBoudaryOperator)
    
    ! Release the temp matrix
    call lsyssc_releaseMatrix (rneumannBoudaryOperator)
  
  end subroutine

end module
