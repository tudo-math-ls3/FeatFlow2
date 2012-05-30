!##############################################################################
!# ****************************************************************************
!# <name> kktsystemspaces </name>
!# ****************************************************************************
!#
!# <purpose>
!# this module realises the different primal/dual/control spaces that
!# appear in the kkt system.
!# </purpose>
!##############################################################################

module kktsystemspaces

  use fsystem
  use genoutput
  
  use spatialdiscretisation
  use timediscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  
  use spacetimevectors
  use analyticsolution
  
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use assemblytemplates
  
  use spacediscretisation
  
  implicit none
  
  private

!<types>

!<typeblock>
  
  ! Encapsules the space of vectors in the primal space.
  type t_primalSpace
  
    ! Space-time vector specifying the primal space.
    type(t_spacetimeVector), pointer :: p_rvector => null()
    
    ! Vector access structure which buffers access to p_rvector.
    type(t_spaceTimeVectorAccess), pointer :: p_rvectorAccess => null()
  end type
  
!</typeblock>

  public :: t_primalSpace

!<typeblock>
  
  ! Encapsules the space of vectors in the dual space.
  type t_dualSpace
  
    ! Space-time vector specifying the dual space.
    type(t_spacetimeVector), pointer :: p_rvector => null()
    
    ! Vector access structure which buffers access to p_rvector.
    type(t_spaceTimeVectorAccess), pointer :: p_rvectorAccess => null()

  end type
  
!</typeblock>

  public :: t_dualSpace

!<typeblock>
  
  ! Encapsules the space of control vectors.
  type t_controlSpace
  
    ! Space-time vector specifying the fully discrete part of the 
    ! control space. However, depending on the problem, there is no
    ! projection applied to the control here!
    type(t_spacetimeVector), pointer :: p_rvector => null()
    
    ! Vector access structure which buffers access to p_rvector.
    type(t_spaceTimeVectorAccess), pointer :: p_rvectorAccess => null()

  end type
  
!</typeblock>

  public :: t_controlSpace

!</types>

  ! Initialise the space discretisation of the primal space
  public :: kktsp_initPrimalSpaceDiscr

  ! Initialise the space discretisation of the dual space
  public :: kktsp_initDualSpaceDiscr

  ! Initialise the space discretisation of the control space
  public :: kktsp_initControlSpaceDiscr

  ! Initialise a primal vector
  public :: kktsp_initPrimalVector

  ! Initialise a dual vector
  public :: kktsp_initDualVector

  ! Initialise a control vector
  public :: kktsp_initControlVector

  ! Release a primal vector
  public :: kktsp_donePrimalVector

  ! Release a dual vector
  public :: kktsp_doneDualVector

  ! Release a control vector
  public :: kktsp_doneControlVector

  ! Linear combination in the primal spacce
  public :: kktsp_primalLinearComb

  ! Linear combination in the dual spacce
  public :: kktsp_dualLinearComb

  ! Linear combination in the control spacce
  public :: kktsp_controlLinearComb
  
  ! Clear a primal vector
  public :: kktsp_clearPrimal

  ! Clear a dual vector
  public :: kktsp_clearDual

  ! Clear a control vector
  public :: kktsp_clearControl

contains

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_initPrimalSpaceDiscr (rspaceDiscr,rphysics,&
      rsettingsDiscr,rtriangulation,rboundary)
  
!<description>
  ! Initialises the space discretisation of the primal space.
!</description>

!<input>
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsDiscr
  
  ! Triangulation structure encapsuling the mesh.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Boundary structure encapsuling the domain
  type(t_boundary), intent(in) :: rboundary
!</input>

!<output>
  ! Spatial discretisation to be initialied.
  type(t_blockDiscretisation), intent(out) :: rspaceDiscr
!</output>

!</subroutine>

    ! Create the corresponding discretisation.
    call spdsc_getDist1LevelDiscr (rboundary,rtriangulation,&
        rspaceDiscr,rsettingsDiscr%ielementType)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_initDualSpaceDiscr (rspaceDiscr,rspaceDiscrPrimal,rphysics,roptControl)
  
!<description>
  ! Initialises the space discretisation of the dual space.
!</description>

!<input>
  ! Space discretisation of the primal space.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscrPrimal

  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Structure defining the optimal control problem to calculate
  type(t_settings_optcontrol), intent(in) :: roptControl
!</input>

!<output>
  ! Spatial discretisation to be initialied.
  type(t_blockDiscretisation), intent(out) :: rspaceDiscr
!</output>

!</subroutine>

    ! Primal and dual space are identical
    call spdiscr_duplicateBlockDiscr (rspaceDiscrPrimal, rspaceDiscr, .true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_initControlSpaceDiscr (rspaceDiscr,&
      rspaceDiscrPrimal,rsettingsDiscr,rphysics,roptControl)
  
!<description>
  ! Initialises the space discretisation of the dual space.
!</description>

!<input>
  ! Space discretisation of the primal space.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscrPrimal

  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Structure defining the optimal control problem to calculate
  type(t_settings_optcontrol), intent(in) :: roptControl

  ! Discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsDiscr
!</input>

!<output>
  ! Spatial discretisation to be initialied.
  type(t_blockDiscretisation), intent(out) :: rspaceDiscr
!</output>

!</subroutine>

    ! The structure of the control space is a bit confusing, since
    ! it depends on the control applied.
    !
    ! At first, check the equation to be discretised.
    select case (rphysics%cequation)
    
    ! -------------------------------------------------------------
    ! Stokes/Navier Stokes.
    ! -------------------------------------------------------------
    case (0,1)
    
      ! -----------------------------------------------------------
      ! Distributed control
      ! -----------------------------------------------------------
      if (roptControl%dalphaC .ge. 0.0_DP) then
        
        ! This is distributed control in the velocity space.
        ! The discretisation matches the discretisation of the primal velocity.
        call spdiscr_deriveBlockDiscr (rspaceDiscrPrimal, rspaceDiscr, 1, 2)
      
        return
        
      end if
    
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_initPrimalVector (rvector,rspaceDiscr,rtimeDiscr)
  
!<description>
  ! Creates a primal vector.
!</description>

!<input>
  ! Spatial discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
  
  ! Time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
!</input>

!<output>
  ! Vector to be created.
  type(t_primalSpace), intent(out) :: rvector
!</output>

!</subroutine>

    ! Allocate memory for the data
    allocate(rvector%p_rvector)
    call sptivec_initVector(rvector%p_rvector,rtimeDiscr,rspaceDiscr)
    
    ! Allocate temp memory for accessing the data
    call sptivec_createAccessPool(rvector%p_rvector,rvector%p_rvectorAccess,5)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_donePrimalVector (rvector)
  
!<description>
  ! Releases a primal vector.
!</description>

!<inputoutput>
  ! Vector to be released.
  type(t_primalSpace), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! Release temp memory
    call sptivec_releaseAccessPool(rvector%p_rvectorAccess)

    ! Release vetor data
    call sptivec_releaseVector (rvector%p_rvector)
    deallocate(rvector%p_rvector)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_initDualVector (rvector,rspaceDiscr,rtimeDiscr)
  
!<description>
  ! Creates a dual vector.
!</description>

!<input>
  ! Spatial discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
  
  ! Time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
!</input>

!<output>
  ! Vector to be created.
  type(t_primalSpace), intent(out) :: rvector
!</output>

!</subroutine>

    ! Allocate memory for the data
    allocate(rvector%p_rvector)
    call sptivec_initVector(rvector%p_rvector,rtimeDiscr,rspaceDiscr)
    
    ! Allocate temp memory for accessing the data
    call sptivec_createAccessPool(rvector%p_rvector,rvector%p_rvectorAccess,5)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_doneDualVector (rvector)
  
!<description>
  ! Releases a dual vector.
!</description>

!<inputoutput>
  ! Vector to be released.
  type(t_primalSpace), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! Release temp memory
    call sptivec_releaseAccessPool(rvector%p_rvectorAccess)

    ! Release vetor data
    call sptivec_releaseVector (rvector%p_rvector)
    deallocate(rvector%p_rvector)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_initControlVector (rvector,rspaceDiscr,rtimeDiscr)
  
!<description>
  ! Creates a control vector.
!</description>

!<input>
  ! Spatial discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
  
  ! Time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
!</input>

!<output>
  ! Vector to be created.
  type(t_primalSpace), intent(out) :: rvector
!</output>

!</subroutine>

    ! Allocate memory for the data
    allocate(rvector%p_rvector)
    call sptivec_initVector(rvector%p_rvector,rtimeDiscr,rspaceDiscr)
    
    ! Allocate temp memory for accessing the data
    call sptivec_createAccessPool(rvector%p_rvector,rvector%p_rvectorAccess,5)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_doneControlVector (rvector)
  
!<description>
  ! Releases a control vector.
!</description>

!<inputoutput>
  ! Vector to be released.
  type(t_primalSpace), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! Release temp memory
    call sptivec_releaseAccessPool(rvector%p_rvectorAccess)

    ! Release vetor data
    call sptivec_releaseVector (rvector%p_rvector)
    deallocate(rvector%p_rvector)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_primalLinearComb (rx,ry,cx,cy,rdest)
  
!<description>
  ! Performs a linear combination
  !    ry = cx rx + cy ry
  ! in the primal space.
!</description>

!<input>
  ! Source vector in the control space.
  type(t_primalSpace), intent(inout) :: rx
  
  ! Multiplication factor
  real(DP), intent(in) :: cx

  ! Multiplication factor
  real(DP), intent(in) :: cy
!</input>

!<inputoutput>
  ! 2nd source vector.
  type(t_primalSpace), intent(inout) :: ry
  
  ! OPTIONAL: Destination vector; may coincide with ry.
  ! If not specified, ry is used.
  type(t_primalSpace), intent(inout), optional :: rdest
!</inputoutput>

!</subroutine>

    ! Currently, we can do this by space-time vector linear combinations.
    if (associated(rx%p_rvector)) then
      call sptivec_vectorLinearComb (rx%p_rvector,ry%p_rvector,cx,cy,rdest%p_rvector)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_dualLinearComb (rx,ry,cx,cy,rdest)
  
!<description>
  ! Performs a linear combination
  !    ry = cx rx + cy ry
  ! in the dual space.
!</description>

!<input>
  ! Source vector in the control space.
  type(t_dualSpace), intent(inout) :: rx
  
  ! Multiplication factor
  real(DP), intent(in) :: cx

  ! Multiplication factor
  real(DP), intent(in) :: cy
!</input>

!<inputoutput>
  ! 2nd source vector.
  type(t_dualSpace), intent(inout) :: ry
  
  ! OPTIONAL: Destination vector; may coincide with ry.
  ! If not specified, ry is used.
  type(t_dualSpace), intent(inout), optional :: rdest
!</inputoutput>

!</subroutine>

    ! Currently, we can do this by space-time vector linear combinations.
    if (associated(rx%p_rvector)) then
      call sptivec_vectorLinearComb (rx%p_rvector,ry%p_rvector,cx,cy,rdest%p_rvector)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_controlLinearComb (rx,cx,ry,cy,rz,cz,dres)
  
!<description>
  ! Performs a linear combination
  !    ry = cx rx + cy ry
  ! in the control space.
!</description>

!<input>
  ! Source vector in the control space.
  type(t_controlSpace), intent(inout) :: rx
  
  ! Multiplication factor
  real(DP), intent(in) :: cx

  ! Multiplication factor
  real(DP), intent(in) :: cy
  
  ! OPTIONAL: Multiplication factor
  real(DP), intent(in), optional :: cz
!</input>

!<inputoutput>
  ! On input: 2nd source vector.
  ! On output: Receives the result if rz not specified
  type(t_controlSpace), intent(inout), target :: ry

  ! OPTIONAL: On input: 3rd source vector.
  ! On output: Receives the result if specified.
  type(t_controlSpace), intent(inout), target, optional :: rz
  
!</inputoutput>

!<output>
  ! OPTIONAL: If specified, the L2-norm of the vector is returned.
  real(DP), intent(out), optional :: dres
!</output>

!</subroutine>

    ! local variables
    type(t_controlSpace), pointer :: p_rdest

    if (present(dres)) dres = -1.0_DP

    if (associated (rx%p_rvector)) then
    
      ! Linear combination of the discrete controls.
      if (.not. present(rz)) then
        call sptivec_vectorLinearComb (rx%p_rvector,ry%p_rvector,cx,cy)
        p_rdest => ry
      else
        call sptivec_vectorLinearComb (rx%p_rvector,rz%p_rvector,cx,cz)
        call sptivec_vectorLinearComb (ry%p_rvector,rz%p_rvector,cy,1.0_DP)
        
        p_rdest => rz
      end if

      ! Calculate the L2-norm
      if (present(dres)) then
        dres = sptivec_vectorNorm (p_rdest%p_rvector,LINALG_NORML2)
      end if

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_clearPrimal (rx)
  
!<description>
  ! Clears a primal vector
!</description>

!<inputoutput>
  ! Control vector to be cleared.
  type(t_primalSpace), intent(inout) :: rx
!</inputoutput>

!</subroutine>

    if (associated (rx%p_rvector)) then
      call sptivec_clearVector (rx%p_rvector)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_clearDual (rx)
  
!<description>
  ! Clears a dual vector
!</description>

!<inputoutput>
  ! Control vector to be cleared.
  type(t_dualSpace), intent(inout) :: rx
!</inputoutput>

!</subroutine>

    if (associated (rx%p_rvector)) then
      call sptivec_clearVector (rx%p_rvector)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_clearControl (rx)
  
!<description>
  ! Clears a control vector.
!</description>

!<inputoutput>
  ! Control vector to be cleared.
  type(t_controlSpace), intent(inout) :: rx
!</inputoutput>

!</subroutine>

    if (associated (rx%p_rvector)) then
      ! Clear the discrete control
      call sptivec_clearVector (rx%p_rvector)
    end if

  end subroutine

end module
