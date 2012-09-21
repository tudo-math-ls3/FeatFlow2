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
  
  use bcassemblybase
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use subdomainfem
  use collection
  
  use spacetimevectors
  use analyticsolution
  
  use constantsdiscretisation
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
  
  ! Copy a vector in the primal spacce
  public :: kktsp_primalCopy

  ! Copy a vector in the dual spacce
  public :: kktsp_dualCopy

  ! Copy a vector in the control spacce
  public :: kktsp_controlCopy

  ! Clear a primal vector
  public :: kktsp_clearPrimal

  ! Clear a dual vector
  public :: kktsp_clearDual

  ! Clear a control vector
  public :: kktsp_clearControl

  ! Scalar product in the control space
  public :: kktsp_scalarProductControl
  
  ! Calculates the norm of a control vector
  public :: kktsp_getNormControl

  ! Compares the norms of the subvectors on the terminal.
  public :: kktsp_controlCompare
  
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
  type(t_settings_spacediscr), intent(in) :: rsettingsDiscr
  
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
        rspaceDiscr,rphysics,rsettingsDiscr%ielementType)

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
      rspaceDiscrPrimal,rtriaBoundary,rsettingsDiscr,rphysics,roptControl)

!<description>
  ! Initialises the space discretisation of the control space.
!</description>

!<input>
  ! Space discretisation of the primal space.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscrPrimal
  
  ! Triangulation of the boundary.
  type(t_triangulation), intent(in) :: rtriaBoundary

  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Structure defining the optimal control problem to calculate
  type(t_settings_optcontrol), intent(in) :: roptControl

  ! Discretisation settings
  type(t_settings_spacediscr), intent(in) :: rsettingsDiscr
!</input>

!<output>
  ! Spatial discretisation to be initialied.
  type(t_blockDiscretisation), intent(out) :: rspaceDiscr
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Ielements
    integer :: nelements

    ! The structure of the control space is a bit confusing, since
    ! it depends on the control applied.
    !
    ! At first, check the equation to be discretised.
    select case (rphysics%cequation)
    
    ! -------------------------------------------------------------
    ! Stokes/Navier Stokes.
    ! -------------------------------------------------------------
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
    
      ! -----------------------------------------------------------
      ! Distributed control
      ! -----------------------------------------------------------
      if (roptControl%dalphaDistC .ge. 0.0_DP) then
        
        ! This is distributed control in the velocity space.
        ! The discretisation matches the discretisation of the primal velocity.
        call spdiscr_deriveBlockDiscr (rspaceDiscrPrimal, rspaceDiscr, 1, 2)
      
        return
        
      end if
      
      if (roptControl%dalphaL2BdC .ge. 0.0_DP) then
      
        ! L2 Dirichlet boundary control. This is a bit harder.
        !
        ! Create a discretisation of the boundary based on the current
        ! discretisation. It should be compatible to the current discretisation.
        
        call spdiscr_initBlockDiscr (rspaceDiscr,2,rtriaBoundary,&
            rspaceDiscrPrimal%p_rboundary)
        
        call spdiscr_createCompDiscrManif1D(rspaceDiscr%RspatialDiscr(1),&
            rspaceDiscrPrimal%RspatialDiscr(1),rtriaBoundary)

        call spdiscr_duplicateDiscrSc (rspaceDiscr%RspatialDiscr(1), &
            rspaceDiscr%RspatialDiscr(2), .true.)
        
        return
      
      end if
    
    ! -------------------------------------------------------------
    ! Heat equation
    ! -------------------------------------------------------------
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
    
      ! -----------------------------------------------------------
      ! Distributed control
      ! -----------------------------------------------------------
      if (roptControl%dalphaDistC .ge. 0.0_DP) then
        
        ! This is distributed control in the velocity space.
        ! The discretisation matches the discretisation of the primal velocity.
        call spdiscr_deriveBlockDiscr (rspaceDiscrPrimal, rspaceDiscr, 1, 1)
      
        return
        
      end if

      if (roptControl%dalphaL2BdC .ge. 0.0_DP) then
      
        ! L2 Dirichlet boundary control. This is a bit harder.
        ! 
        !
        ! Create a discretisation of the boundary based on the current
        ! discretisation. It should be compatible to the current discretisation.
        
        call spdiscr_initBlockDiscr (rspaceDiscr,2,rtriaBoundary,&
            rspaceDiscrPrimal%p_rboundary)
        
        call spdiscr_createCompDiscrManif1D(rspaceDiscr%RspatialDiscr(1),&
            rspaceDiscrPrimal%RspatialDiscr(1),rtriaBoundary)
            
        return
      
      end if
    
    end select

  end subroutine

  ! ***************************************************************************
  
  subroutine kktsp_getElementsOnBd (rboundary,rtriangulation,p_Ielements,nelements)
  
  ! Auxiliary routine. Determines all elements on the boundary of the domain.
  
  ! Underlying boundary
  type(t_boundary), intent(in) :: rboundary
  
  ! underlying trinagulation
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! OUT: Pointer to a list of elements on the boundary
  integer, dimension(:), pointer :: p_Ielements
  
  ! OUT: number of elements on the boundary
  integer, intent(out) :: nelements
  
    ! local variables
    integer :: nel,i,j,ibct
    type(t_boundaryRegion) :: rregion
    
    ! Determine maximum total number of elements on the boundary
    nelements = 0
    do i=1,boundary_igetNBoundComp(rboundary)
      call boundary_createRegion (rboundary, i, 0, rregion)
      call bcasm_getElementsInBdRegion (rtriangulation,rregion,nel)
      nelements = nelements + nel
    end do

    ! Calcel here if there are no elements
    nullify(p_Ielements)
    if (nelements .eq. 0) return

    ! Allocate a list where to save the element numbers.
    allocate(p_Ielements(nelements))
    
    ! Get all the elements
    nelements = 0
    do ibct=1,boundary_igetNBoundComp(rboundary)
      call boundary_createRegion (rboundary, ibct, 0, rregion)
      call bcasm_getElementsInBdRegion (rtriangulation,rregion,nel,IelList=p_Ielements(nelements+1:))
      
      if (nel .gt. 1) then
        ! Cancel out duplicates as far as possible.
        i = 1
        j = 1
        do i = 2,nel
          if (p_Ielements(nelements+j) .ne. p_Ielements(nelements+i)) then
            j = j + 1
            p_Ielements(nelements+j) = p_Ielements(nelements+i)
          end if
        end do

        ! Actual number of elements
        nel = j
        
        ! Reduce by one if first and last element match.
        if (nel .gt. 1) then
          if (p_Ielements(nelements+1) .eq. p_Ielements(nelements+nel)) then
            nel = nel - 1
          end if
        end if

      end if      
      
      nelements = nelements + nel
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_initPrimalVector (&
      rvector,rspaceDiscr,rtimeDiscr,istartidx,iendidx)
  
!<description>
  ! Creates a primal vector.
!</description>

!<input>
  ! Spatial discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
  
  ! Time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr

  ! OPTIONAL: First subvector. If not specified, this defaults to 1.
  ! If specified, the vectors 1..istartidx are not created.
  integer, intent(in), optional :: istartidx

  ! OPTIONAL: Last subvector. If not specified, this defaults to #NEQ in time.
  ! If specified, the vectors iendidx+1..#NEQ in time are not created.
  integer, intent(in), optional :: iendidx
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
    allocate(rvector%p_rvectorAccess)
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
    deallocate(rvector%p_rvectorAccess)
    
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
  type(t_dualSpace), intent(out) :: rvector
!</output>

!</subroutine>

    ! Allocate memory for the data
    allocate(rvector%p_rvector)
    call sptivec_initVector(rvector%p_rvector,rtimeDiscr,rspaceDiscr)
    
    ! Allocate temp memory for accessing the data
    allocate(rvector%p_rvectorAccess)
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
  type(t_dualSpace), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! Release temp memory
    call sptivec_releaseAccessPool(rvector%p_rvectorAccess)
    deallocate(rvector%p_rvectorAccess)

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
  type(t_controlSpace), intent(out) :: rvector
!</output>

!</subroutine>

    ! Allocate memory for the data
    allocate(rvector%p_rvector)
    call sptivec_initVector(rvector%p_rvector,rtimeDiscr,rspaceDiscr)
    
    ! Allocate temp memory for accessing the data
    allocate(rvector%p_rvectorAccess)
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
  type(t_controlSpace), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! Release temp memory
    call sptivec_releaseAccessPool(rvector%p_rvectorAccess)
    deallocate(rvector%p_rvectorAccess)

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
      
      ! Invalidate the buffer
      call sptivec_invalidateVecInPool (rdest%p_rvectorAccess)
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
      if (present(rdest)) then
        call sptivec_vectorLinearComb (rx%p_rvector,ry%p_rvector,cx,cy,rdest%p_rvector)

        ! Invalidate the buffer
        call sptivec_invalidateVecInPool (rdest%p_rvectorAccess)
      else
        call sptivec_vectorLinearComb (rx%p_rvector,ry%p_rvector,cx,cy)

        ! Invalidate the buffer
        call sptivec_invalidateVecInPool (ry%p_rvectorAccess)
      end if

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_controlLinearComb (rx,cx,ry,cy,rz,cz)
  
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

!</subroutine>

    ! local variables
    type(t_controlSpace), pointer :: p_rdest

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

      ! Invalidate the buffer
      call sptivec_invalidateVecInPool (p_rdest%p_rvectorAccess)

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  real(DP) function kktsp_getNormControl (rx,iresnorm)
  
!<description>
  ! Calculates the norm of a control vector
!</description>

!<input>
  ! Source vector in the control space.
  type(t_controlSpace), intent(inout) :: rx

  ! type of norm. A LINALG_NORMxxxx constant.
  ! Must be specified if dres is specified.
  integer, intent(in) :: iresnorm
!</input>

!<result>
  ! Norm of the vector.
!</result>

!</subroutine>
    
    real(DP) :: dtemp

    dtemp = rx%p_rvector%p_Dscale(1)
    rx%p_rvector%p_Dscale(1) = 0.0_DP
    kktsp_getNormControl = sptivec_vectorNorm (rx%p_rvector,iresnorm)
    rx%p_rvector%p_Dscale(1) = dtemp

  end function

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
      call sptivec_invalidateVecInPool (rx%p_rvectorAccess)
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
      call sptivec_invalidateVecInPool (rx%p_rvectorAccess)
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
      call sptivec_invalidateVecInPool (rx%p_rvectorAccess)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_primalCopy (rx,ry)
  
!<description>
  ! Copies a vector, ry = rx,
  ! in the primal space.
!</description>

!<input>
  ! Source vector in the control space.
  type(t_primalSpace), intent(inout) :: rx
!</input>

!<inputoutput>
  ! Destination vector.
  type(t_primalSpace), intent(inout) :: ry
!</inputoutput>

!</subroutine>

    call kktsp_primalLinearComb (rx,ry,1.0_DP,0.0_DP)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_dualCopy (rx,ry)
  
!<description>
  ! Copies a vector, ry = rx,
  ! in the dual space.
!</description>

!<input>
  ! Source vector.
  type(t_dualSpace), intent(inout) :: rx
!</input>

!<inputoutput>
  ! Destination vector.
  type(t_dualSpace), intent(inout) :: ry
!</inputoutput>

!</subroutine>

    call kktsp_dualLinearComb (rx,ry,1.0_DP,0.0_DP)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_controlCopy (rx,ry)
  
!<description>
  ! Copies a vector, ry = rx,
  ! in the dual space.
!</description>

!<input>
  ! Source vector.
  type(t_controlSpace), intent(inout) :: rx
!</input>

!<inputoutput>
  ! Destination vector.
  type(t_controlSpace), intent(inout) :: ry
!</inputoutput>

!</subroutine>

    call kktsp_controlLinearComb (rx,1.0_DP,ry,0.0_DP)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  real(DP) function kktsp_scalarProductControl (rx,ry)
  
!<description>
  ! Calculates the scalar product in the control space: dres = rx^T rx.
!</description>

!<input>
  ! Source vector.
  type(t_controlSpace), intent(inout) :: rx

  ! 2nd source vector.
  type(t_controlSpace), intent(inout) :: ry
!</input>

!<result>
  ! Scalar product.
!</result>

!</subroutine>

    real(DP) :: dtemp

    kktsp_scalarProductControl = 0.0_DP
    if (associated (rx%p_rvector)) then
      ! Calculate the scalar product.

      dtemp = rx%p_rvector%p_Dscale(1)
      rx%p_rvector%p_Dscale(1) = 0.0_DP
      kktsp_scalarProductControl = sptivec_scalarProduct (rx%p_rvector, ry%p_rvector)
      rx%p_rvector%p_Dscale(1) = dtemp

    end if
    
  end function

  ! ***************************************************************************

!<subroutine>

  subroutine kktsp_controlCompare (rx,ry)
  
!<description>
  ! Prints the norms of the subvectors in rx and ry to the terminal.
!</description>

!<input>
  ! Source vector in the control space.
  type(t_controlSpace), intent(inout) :: rx
  
  ! Multiplication factor
  type(t_controlSpace), intent(inout) :: ry
!</input>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rx,p_ry
    integer :: i
    real(DP) :: dnorm1, dnorm2

    if (associated (rx%p_rvector)) then
    
      do i=1,rx%p_rvector%NEQtime
        call sptivec_getVectorFromPool (rx%p_rvectorAccess,i,p_rx)
        call sptivec_getVectorFromPool (ry%p_rvectorAccess,i,p_ry)
        
        dnorm1 = lsysbl_vectorNorm (p_rx,LINALG_NORML2)
        dnorm2 = lsysbl_vectorNorm (p_ry,LINALG_NORML2)
        call output_line (&
            "Vector comparison: ||rx("//trim(sys_siL(i,10))//")|| = "//trim(sys_sdEL(dnorm1,10))//&
            ", ||ry("//trim(sys_siL(i,10))//")|| = "//trim(sys_sdEL(dnorm2,10)))
      end do

    end if

  end subroutine

end module
