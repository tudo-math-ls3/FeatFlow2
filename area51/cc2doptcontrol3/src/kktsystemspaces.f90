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
